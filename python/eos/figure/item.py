# Copyright (c) 2023-2026 Danny van Dyk
# Copyright (c) 2023      Philip Lueghausen
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from abc import abstractmethod
from dataclasses import dataclass, field
from eos.analysis_file_description import AnalysisFileContext
from eos.deserializable import Deserializable

import inspect
import eos
import matplotlib as _matplotlib
import matplotlib.collections
import matplotlib.container
import matplotlib.lines
import matplotlib.patches
import matplotlib.transforms
import numpy as _np
import os
import scipy as _scipy
import yaml as _yaml

class ItemColorCycler:
    """Cycles through a fixed, colorblind-friendly palette of colors.

    Items that are not assigned an explicit color draw their next color from this shared cycle, so that
    consecutive items in a plot are automatically given distinct colors. The cycle is reset, via
    :meth:`reset`, whenever a new plot is created.
    """

    _colors = [
        # 10-colors colorblind-friendly cycle from 2107.02270
        (0.247, 0.565, 0.855), (1.000, 0.663, 0.055), (0.741, 0.122, 0.004), (0.580, 0.643, 0.635), (0.514, 0.176, 0.714), (0.663, 0.420, 0.349), (0.906, 0.388, 0.000), (0.725, 0.675, 0.439), (0.443, 0.459, 0.506), (0.573, 0.855, 0.867)
    ]
    _color_idx = 0

    @classmethod
    def reset(cls):
        """Reset the color cycle so that the next call to :meth:`next_color` returns the first color again."""
        cls._color_idx = 0

    @classmethod
    def next_color(cls):
        """Return the next color in the cycle.

        Successive calls cycle through a fixed, colorblind-friendly palette of ten colors,
        wrapping around to the first color once the palette is exhausted.

        :returns: An RGB color as a tuple of three floats in the range :math:`[0, 1]`.
        :rtype: tuple[float, float, float]
        """
        result = cls._colors[cls._color_idx]
        cls._color_idx = (cls._color_idx + 1) % len(cls._colors)
        return result

@dataclass(kw_only=True)
class Item(Deserializable):
    r"""Base class for items to be drawn in a figure."""

    alpha:float=field(default=0.5)
    color:str=field(default=None)
    label:str=field(default=None)
    linestyle:str=field(default='solid')
    linewidth:float=field(default=1.0)

    def __post_init__(self):
        if self.color is None:
            self.color = ItemColorCycler.next_color()

    @abstractmethod
    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the item for drawing.

        Subclasses override this method to load any required data files and to evaluate and
        cache the quantities that :meth:`draw` will render. It is called once before :meth:`draw`.

        :param context: The analysis file context used to resolve relative paths to data files.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        raise NotImplementedError

    @abstractmethod
    def draw(self, ax):
        """Draw the item on the provided axes.

        Subclasses override this method to render the quantities cached by :meth:`prepare`.

        :param ax: The matplotlib axes onto which the item is drawn.
        :type ax: matplotlib.axes.Axes
        """
        raise NotImplementedError

    def _legend_line(self, alpha=None):
        "A line swatch matching a plotted curve; pass the item's alpha if it draws with one."
        if not self.label:
            return []
        handle = _matplotlib.lines.Line2D((0, 1), (0, 0), color=self.color, lw=self.linewidth, ls=self.linestyle, alpha=alpha)
        return [(handle, self.label)]

    def _legend_patch(self):
        "A translucent filled swatch matching a band or histogram."
        if not self.label:
            return []
        handle = _matplotlib.patches.Rectangle((0, 0), 1, 1, color=self.color, alpha=self.alpha)
        return [(handle, self.label)]

    def _legend_marker(self, marker, alpha=None):
        "A point-marker swatch matching error bars; pass the item's alpha if it draws with one."
        if not self.label:
            return []
        handle = _matplotlib.lines.Line2D((0,), (0,), color=self.color, marker=marker, linestyle='none', alpha=alpha)
        return [(handle, self.label)]

    def _legend_errorbar(self, marker='_', has_xerr=False, has_yerr=True, alpha=None):
        """An error-bar swatch (capped bar) matching :class:`ErrorBarsItem`.

        Returns an :class:`matplotlib.container.ErrorbarContainer`, which the default legend
        renders via ``HandlerErrorbar`` as a central marker with caps and a connecting bar.
        The central element is a marker (not a full-width line) and the caps are sized to
        ``2 * errorbar.capsize`` -- matplotlib's own convention -- so the swatch reproduces
        the error bars as the items actually draw them (no connecting line, capped bar). The
        template artists below carry the visual properties (colour, width, marker size, alpha);
        the handler copies them onto the swatch geometry it builds.
        """
        if not self.label:
            return []
        capsize  = 2.0 * _matplotlib.rcParams['errorbar.capsize']
        line     = _matplotlib.lines.Line2D((0,), (0,), color=self.color, marker=marker, linestyle='none', alpha=alpha)
        capline  = _matplotlib.lines.Line2D((0,), (0,), color=self.color, marker='_', linestyle='none', markersize=capsize, alpha=alpha)
        barlines = _matplotlib.collections.LineCollection([[(0, 0), (0, 1)]], colors=self.color, linewidths=self.linewidth, alpha=alpha)
        handle   = _matplotlib.container.ErrorbarContainer((line, [capline], [barlines]), has_xerr=has_xerr, has_yerr=has_yerr)
        return [(handle, self.label)]

    def legend(self):
        """Return the item's legend entry as a list of handle/label pairs.

        The default implementation returns an empty list, i.e. the item contributes no
        dedicated legend entry. Subclasses that require a custom legend handle override this

        :returns: A list of ``(handle, label)`` pairs to be added to the plot legend.
        :rtype: list[tuple[matplotlib.artist.Artist, str]]
        """
        return []


@dataclass(kw_only=True)
class ObservableItem(Item):
    """Plots an observables from the EOS library of builtin observables.

    This item type is used to plot observables as a function of one kinematic variable or one parameter.

    :param observable: The name of the observable to be plotted.
    :type observable: :class:`eos.QualifiedName`
    :param fixed_kinematics: A dictionary of additional kinematic variables and the values to which they shall be fixed.
    :type fixed_kinematics: dict[str, float] | None
    :param fixed_parameters: A dictionary of EOS parameters and the values to which they shall be fixed.
    :type fixed_parameters: dict[eos.QualifiedName, float] | None
    :param fixed_parameters_from_file: Path to a file containing EOS parameters to be used for the observable.
    :type fixed_parameters_from_file: str | None
    :param label: The label for the observable, which will be used in the plot legend.
    :type label: str | None
    :param options: A dictionary of options to be passed to the observable.
    :type options: dict[str, str] | None
    :param range: A tuple of two float values (min, max) representing the range of the kinematic variable to be plotted on the x-axis.
    :type range: tuple[float, float]
    :param resolution: The number of points to be used for the evaluation of the observable. Defaults to 100.
    :type resolution: int
    :param variable: The name of the kinematic or parameter variable to which the x-axis will be mapped; see
        `the complete list of observables <../reference/observables.html>`_ and
        `the complete list of parameters <../reference/parameters.html>`_.
    :type variable: str | :class:`eos.QualifiedName`

    Example:

    .. code-block::

        figure_args = '''
        plot:
          legend: { position: 'upper center' }
          xaxis: { label: r'$q^2$',                range: [0.0, 11.6]   }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$', range: [0.0, 5.0e-3] }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'mu' },  label: r'$\\ell = \\mu$',
                variable: 'q2', range: [0.02, 11.6 ], resolution: 5800
              }
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'tau' }, label: r'$\\ell = \\tau$',
                variable: 'q2', range: [3.17, 11.6 ], resolution: 421
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    observable:eos.QualifiedName
    fixed_kinematics:dict=None
    fixed_parameters:dict=None
    fixed_parameters_from_file:str=None
    label:str=None
    options:dict=None
    range:tuple[float,float]
    resolution:int=field(default=100)
    variable:str

    _api_doc = inspect.cleandoc("""\
    Plotting Observables
    --------------------

    Plot items of type ``observable`` are used to display one of the built-in `observables <../reference/observables.html>`_.

    The following keys are mandatory:

        * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable that will be plotted.
        Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_.
        * ``range`` (*list* or *tuple* of two *float*) --The tuple of [minimal, maximal] values of the specified kinematic variable
        for which the observable will be evaluated.
        * ``variable`` (*str*) -- The name of the kinematic or parameter variable to which the x-axis will be mapped; see
        `the complete list of parameters <../reference/parameters.html>`_.

    The following keys are optional:

        * ``fixed_kinematics`` (*dict[str, float]*) -- A dictionary of additional kinematic variables and the values to which they shall be fixed.
        * ``fixed_parameters`` (*dict[eos.QualifiedName, float]*) -- A dictionary of EOS parameters and the values to which they shall be fixed.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          legend: { position: 'upper center' }
          xaxis: { label: r'$q^2$',                range: [0.0, 11.6]   }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$', range: [0.0, 5.0e-3] }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'mu' },  label: r'$\\ell = \\mu$',
                variable: 'q2', range: { min: 0.02, max: 11.6, num: 5800 }
              }
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'tau' }, label: r'$\\ell = \\tau$',
                variable: 'q2', range: { min: 3.17, max: 11.6, num: 421 }
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """)

    def __post_init__(self):
        super().__post_init__()
        eos.info(f'Handling item to plot {self.observable}')
        self._observable_entry = eos.Observables()[self.observable]
        valid_kinematic_variables = {kv for kv in self._observable_entry.kinematic_variables()}

        # Create kinematics
        self._kinematics = eos.Kinematics()
        if self.fixed_kinematics:
            for k, v in self.fixed_kinematics.items():
                if k not in valid_kinematic_variables:
                    raise ValueError(f"Kinematic variable '{k}' does not match any of the declared kinematic variables '{self.observable}': {valid_kinematic_variables.__repr__()}")
                self._kinematics.declare(k, v)

        # Create parameters
        self._parameters = eos.Parameters.Defaults()
        if self.fixed_parameters_from_file:
            eos.warn('Overriding parameters from file')
            self._parameters.override_from_file(self.fixed_parameters_from_file)

        if self.fixed_parameters:
            if self.fixed_parameters_from_file:
                eos.warn('Overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            for key, value in self.fixed_parameters.items():
                self._parameters.set(key, value)

        # Declare variable that is either parameter or kinematic
        self._variable = None

        # Does the variable correspond to one of the kinematic variables?
        if self.variable in valid_kinematic_variables:
            self._variable = self._kinematics.declare(self.variable, _np.nan)
        else:
            # Is the variable name a QualifiedName?
            try:
                qn = eos.QualifiedName(self.variable)
                # Continues only if no failure occurs
                self._variable = self.parameters[qn]
                self._variable.set(_np.nan)
            except RuntimeError:
                raise ValueError(f"Value of 'variable' for observable '{self.observable}' is neither a valid kinematic variable nor a valid parameter: '{self.variable}")

        # Create options
        self._options = eos.Options()
        if self.options:
            for key, value in self.options.items():
                self._options.declare(key, value)

        self._observable = eos.Observable.make(self.observable, self._parameters, self._kinematics, self._options)

        self._xvalues = _np.linspace(self.range[0], self.range[1], self.resolution)


    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the drawing by evaluating the observable at the sample points.

        Evaluates the observable on a grid of ``resolution`` points spanning ``range`` of the
        chosen kinematic or parameter ``variable`` and caches the results for :meth:`draw`.

        :param context: The analysis file context used to resolve relative paths. If ``None``,
            a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context
        self._yvalues = _np.empty((len(self._xvalues),))
        for i, x in enumerate(self._xvalues):
            self._variable.set(x)
            self._yvalues[i] = self._observable.evaluate()


    def draw(self, ax, **kwargs):
        """Draw a curve of the observable as a function of the chosen variable.

        :param ax: The matplotlib axes onto which the curve is drawn.
        :type ax: matplotlib.axes.Axes
        :param kwargs: Additional keyword arguments forwarded to :meth:`matplotlib.axes.Axes.plot`.
        """
        ax.plot(self._xvalues, self._yvalues, label=self.label, color=self.color, lw=self.linewidth, ls=self.linestyle, **kwargs)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_line()


# Namespace exposing NumPy's functions and constants (e.g. ``sin``, ``exp``, ``pi``) to the
# expressions evaluated by :class:`ExpressionItem`, so that they can be used without a ``np.``
# prefix. ``np`` itself is also made available for any name not exported at NumPy's top level.
_EXPRESSION_NAMESPACE = {name: getattr(_np, name) for name in dir(_np) if not name.startswith('_')}
_EXPRESSION_NAMESPACE['np'] = _np


@dataclass(kw_only=True)
class ExpressionItem(Item):
    r"""Plots an arbitrary mathematical expression as a function of the x-axis variable.

    This item type evaluates a user-supplied Python ``expression`` of the free variable ``x`` on a
    grid spanning ``range`` and draws the result as a curve. The expression is evaluated with
    Python's :func:`eval`, with NumPy's functions and constants (e.g. ``sin``, ``exp``, ``sqrt``,
    ``pi``) available by name as well as the ``np`` module itself; for example,
    ``'exp(-x**2) * sin(2 * pi * x)'``. Because the expression is evaluated as code, it should only
    be populated from trusted figure descriptions.

    :param expression: The expression to be plotted, given as a Python expression in the free variable ``x``.
    :type expression: str
    :param range: A tuple of two float values (min, max) representing the range of the x-axis variable over which the expression is evaluated.
    :type range: tuple[float, float]
    :param resolution: The number of points to be used for the evaluation of the expression. Defaults to 100.
    :type resolution: int

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$x$', range: [0.0, 6.28] }
          yaxis: { label: r'$f(x)$' }
          items:
            - { type: 'expression', expression: 'sin(x)',           label: r'$\sin x$'       }
            - { type: 'expression', expression: 'exp(-x) * cos(x)', label: r'$e^{-x}\cos x$' }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    expression:str
    range:tuple[float,float]
    resolution:int=field(default=100)

    def __post_init__(self):
        super().__post_init__()

        if not isinstance(self.expression, str) or not self.expression.strip():
            raise ValueError("'expression' must be a non-empty string")

        if self.range is None or len(self.range) != 2:
            raise ValueError(f"Invalid range '{self.range}'. It must be a tuple of two float values (min, max).")

        if self.range[0] >= self.range[1]:
            raise ValueError(f"Invalid range '{self.range}'. The first value must be less than the second value.")

        if self.resolution <= 0:
            raise ValueError(f"Invalid resolution '{self.resolution}'. It must be a positive integer.")

        # compile the expression eagerly so that syntax errors surface at construction time
        try:
            self._code = compile(self.expression, '<eos.figure.ExpressionItem>', 'eval')
        except SyntaxError as e:
            raise ValueError(f"Could not parse expression '{self.expression}': {e}")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the drawing by evaluating the expression at the sample points.

        Evaluates ``expression`` on a grid of ``resolution`` points spanning ``range`` and caches the
        results for :meth:`draw`. The free variable ``x`` is passed as a NumPy array, so vectorized
        expressions are evaluated in a single pass; constant expressions are broadcast to the grid.

        :param context: The analysis file context. Accepted for interface consistency and unused.
        :type context: AnalysisFileContext | None
        """
        self._xvalues = _np.linspace(self.range[0], self.range[1], self.resolution)
        # wrap evaluation as well as the conversion/broadcast so that shape/dtype issues (e.g. a
        # complex or non-broadcastable result) are reported with the offending expression too
        try:
            yvalues = eval(self._code, {'__builtins__': {}}, {**_EXPRESSION_NAMESPACE, 'x': self._xvalues})
            self._yvalues = _np.broadcast_to(_np.asarray(yvalues, dtype=float), self._xvalues.shape)
        except Exception as e:
            raise ValueError(f"Could not evaluate expression '{self.expression}': {e}")

    def draw(self, ax, **kwargs):
        """Draw a curve of the expression as a function of the x-axis variable.

        :param ax: The matplotlib axes onto which the curve is drawn.
        :type ax: matplotlib.axes.Axes
        :param kwargs: Additional keyword arguments forwarded to :meth:`matplotlib.axes.Axes.plot`.
        """
        ax.plot(self._xvalues, self._yvalues, label=self.label, alpha=self.alpha, color=self.color,
                lw=self.linewidth, ls=self.linestyle, **kwargs)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_line(alpha=self.alpha)


@dataclass(kw_only=True)
class UncertaintyBandItem(Item):
    """Plots the 68% uncertainty band as a function of one kinematic variable or one parameter.

    This routine requires EOS data file to have been produced, typically by running ``predict-observables`` task.

    :param band: A set of strings that determines which parts of the uncertainty band to draw. Can be any combination of ``'area'``, ``'outer'``, or ``'median'``.
        Defaults to ``{'area', 'outer', 'median'}``.
        ``'area'`` fills the area between the lower and upper bounds of the uncertainty band,
        ``'outer'`` draws the outer lines of the uncertainty band, and ``'median'`` draws the median line of the uncertainty band.
    :type band: set[str] | list[str] | str
    :param datafile: The path to an existing data file of type :class:`eos.data.Prediction` that contains the uncertainty estimates.
    :type datafile: str
    :param interpolation: The type of interpolation to be used for the band. Can be either ``linear`` (default) or ``cubic``.
    :type interpolation: str
    :param range: A tuple of two float values (min, max) representing the range of the kinematic variable to be plotted on the x-axis.
    :type range: tuple[float, float] | None
    :param resolution: The number of points to be used for the interpolation of the band. Defaults to 100.
    :type resolution: int
    :param variable: The name of the kinematic variable that is plotted on the x-axis. Defaults to the first kinematic variable recorded in the data file.
    :type variable: str | None
    :param observable: The name of the observable that is plotted. This name can include options. Defaults to the first observable recorded in the data file.
    :type observable: str | None

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                  range: [0.0,  5e-3] }
          legend: { position: 'upper center' }
          items:
            - { type: 'uncertainty', label: r'$\\ell=\\mu$',
                variable: 'q2', range: [0.02, 11.63],
                datafile: './predictions-data/FF-LQCD-SSE/pred-B-to-D-mu-nu'
              }
            - { type: 'uncertainty', label: r'$\\ell=\\tau$',
                variable: 'q2', range: [3.17, 11.63], resolution: 1000,
                datafile: './predictions-data/FF-LQCD-SSE/pred-B-to-D-tau-nu'
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    band:set[str]|list[str]|str=field(default_factory=lambda : {'area', 'outer', 'median'})
    datafile:str
    interpolation:str=field(default='linear')
    range:tuple[float, float]|None=field(default=None)
    resolution:int=field(default=100)
    variable:str|None=field(default=None)
    observable:str|None=field(default=None)

    _api_doc = inspect.cleandoc("""\
    Plotting Uncertainty Bands
    --------------------------

    Plot items of type ``uncertainty`` are used to display uncertainty bands at 68% probability
    for one of the built-in `observables <../reference/observables.html>`_. This item type requires that
    a predictive distribution of the observables has been previously produced.

    The following key is mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.Prediction*) -- The path to
        a data file that was generated with the ``predict-observables`` task.

    The following keys are optional:

        * ``band`` (a *list* or *set* containing any combination of ``'area'``, ``'outer'``, or ``'median'``) -- The setting for the illustration of the uncertainty band.
        If ``'area'`` is provided, the band's areas are filled.
        If ``'outer'`` is provided, the band's outer lines are drawn.
        If ``'median'`` is provided, the band's median line is drawn.
        Defaults to ``{'area', 'outer', 'median'}``.
        * ``interpolation`` (*str*) -- The type of interpolation to be used for the band. Can be either ``linear`` (default) or ``cubic``.
        * ``range`` (*tuple* of two *float* values) -- The range of the kinematic variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.
        * ``resolution`` (*int*) -- The number of points to be used for the interpolation of the band. Defaults to 100.
        * ``variable`` (*str*) -- The name of the kinematic variable that is plotted on the x-axis. Defaults to the first kinematic variable in the data file.
        * ``observable`` (*str*) -- The name of the observable that is plotted. This name can include options. Defaults to the first observable in the data file.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
          legend: { position: 'upper center' }
          items:
            - { type: 'uncertainty', label: r'$\\ell=\\mu$',
                variable: 'q2', range: [0.02, 11.63],
                datafile: './predictions-data/FF-LQCD-SSE/pred-B-to-D-mu-nu'
              }
            - { type: 'uncertainty', label: r'$\\ell=\\tau$',
                variable: 'q2', range: [3.17, 11.63], resolution: 1000,
                datafile: './predictions-data/FF-LQCD-SSE/pred-B-to-D-tau-nu'
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """)

    def __post_init__(self):
        super().__post_init__()

        if isinstance(self.band, str):
            self.band = {self.band}
        elif isinstance(self.band, list):
            self.band = set(self.band)
        elif not isinstance(self.band, set):
            raise TypeError(f"Parameter 'band' must be a string, list of string, or a set of strings, not {type(self.band).__name__}")

        for band_type in self.band:
            if band_type not in ['area', 'outer', 'median']:
                raise ValueError(f"Unrecognized band type '{band_type}'; must be one of 'area', 'outer', or 'median'")

        if self.interpolation not in ['linear', 'cubic']:
            raise ValueError(f"Unrecognized interpolation type '{self.interpolation}'; must be either 'linear' or 'cubic'")

        if self.range is not None and len(self.range) != 2:
            raise ValueError(f"Range must be a tuple of two values (min, max), not {self.range}")

        if self.range is not None:
            self.range = tuple(float(x) for x in self.range)
            if self.range[0] >= self.range[1]:
                raise ValueError(f"Range must be a tuple of two values (min, max) with max > min, not {self.range}")


    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the uncertainty band for drawing.

        Loads the predictive samples from the data file, selects the requested observable and
        kinematic variable, computes the central 68% interval at each sample point, and
        interpolates the lower, central, and upper bounds onto a regular grid for :meth:`draw`.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context

        self._datafile = eos.data.Prediction(context.data_path(self.datafile))
        if len(self._datafile.varied_parameters) < 1:
            raise ValueError(f"Data file '{self.datafile}' does not contain any predictions")

        self._variable = self.variable
        if self._variable is None:
            # assume that the first prediction has only one kinematic variable,
            # which we adopt as the variable to plot
            kinematics = self._datafile.varied_parameters[0]['kinematics']
            if len(kinematics) > 1:
                raise ValueError(f"Data file '{self.datafile}' contains more than one kinematic variable; specify 'variable' to determine which is supposed to be used")
            self._variable = list(kinematics.keys())[0]

        self._observable = self.observable
        if self._observable is None:
            # assume that only one observable was predicted
            observables = {obs.split(';')[0] for obs in list(self._datafile.lookup_table.keys())}
            if len(observables) > 1:
                raise ValueError(f"Data file '{self.datafile}' contains more than one predicted observable; specify 'observable' to determine which is supposed to be used")
            if len(observables) == 0:
                raise ValueError(f"Data file '{self.datafile}' contains no observable, check predictions")
            self._observable = observables.pop()
        self._observable = eos.QualifiedName(self._observable)

        # Filter the x values and the predictions based on _variable and _observable
        _xvalues = []
        _observable_mask = []
        for idx, vp in enumerate(self._datafile.varied_parameters):
            kinematics = vp['kinematics']
            if self._variable not in kinematics:
                raise ValueError(f"Prediction #{idx} in data file '{self.datafile}' does not depend on the chosen kinematic variable '{self._variable}'")
            _xvalues.append(kinematics[self._variable])

            valid_observable = (self._observable == eos.QualifiedName(vp["name"])) # Check that the observable names match
            for k, v in self._observable.options_part(): # Check that thre specified options match
                if k in vp["options"] and not vp["options"][k] == v:
                    valid_observable = False

            _observable_mask.append(valid_observable)

        if not _np.any(_observable_mask):
            raise ValueError(f"No observable of the data file matches the specified observable '{self.observable}'")

        _xvalues = _np.array(_xvalues)[_observable_mask]
        _samples = self._datafile.samples[:, _observable_mask]
        _weights = self._datafile.weights

        _ovalues_lower   = []
        _ovalues_central = []
        _ovalues_higher  = []
        INTERVAL = [0.15865, 0.5, 0.84135]  # central 68% interval
        for i in range(len(_samples[0])):
            lower, central, higher = _np.quantile(_samples[:, i], q = INTERVAL, weights = _weights, method='inverted_cdf', axis=0)
            _ovalues_lower.append(lower)
            _ovalues_central.append(central)
            _ovalues_higher.append(higher)

        self._xvalues = _np.linspace(_np.min(_xvalues), _np.max(_xvalues), self.resolution)

        if self.interpolation == "linear":
            interpolate = lambda x, y, xv: _np.interp(xv, x, y)
        elif self.interpolation == "cubic":
            from scipy.interpolate import CubicSpline
            interpolate = lambda x, y, xv: CubicSpline(x, y)(xv)

        self._ovalues_lower   = interpolate(_xvalues, _ovalues_lower,   self._xvalues)
        self._ovalues_central = interpolate(_xvalues, _ovalues_central, self._xvalues)
        self._ovalues_higher  = interpolate(_xvalues, _ovalues_higher,  self._xvalues)

        if self.range is not None:
            self._xvalues         = _np.ma.masked_outside(self._xvalues,         self.range[0], self.range[1])

    def draw(self, ax):
        """Draw the uncertainty band on the provided axes.

        Depending on the ``band`` setting, fills the area between the lower and upper bounds,
        draws the outer lines, and/or draws the median line.

        :param ax: The matplotlib axes onto which the band is drawn.
        :type ax: matplotlib.axes.Axes
        """
        label = self.label
        if 'area' in self.band:
            ax.fill_between(self._xvalues, self._ovalues_lower, self._ovalues_higher, alpha=self.alpha, color=self.color, label=label, lw=0)
            label = None # do not label anything else if we fill the band area
        if 'outer' in self.band:
            ax.plot(self._xvalues, self._ovalues_lower,                               alpha=self.alpha, color=self.color, label=label, lw=self.linewidth, ls=self.linestyle)
            ax.plot(self._xvalues, self._ovalues_higher,                              alpha=self.alpha, color=self.color,              lw=self.linewidth, ls=self.linestyle)
            label = None # do not label anything else if we plot the outer lines
        if 'median' in self.band:
            ax.plot(self._xvalues, self._ovalues_central,                             alpha=self.alpha, color=self.color, label=label, lw=self.linewidth, ls=self.linestyle)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_patch()


@dataclass(kw_only=True)
class BinnedUncertaintyItem(Item):
    """Plots one or more uncertainty band integrated over one kinematic variable.

    This routine expects the uncertainty propagation to have produced a data file.

    :param datafile: The path to an existing data file of type :class:`eos.data.Prediction` that contains the uncertainty estimates.
    :type datafile: str
    :param range: A tuple of two float values (min, max) representing the range of the kinematic variable to be plotted on the x-axis.
    :type range: tuple[float, float] | None
    :param rescale_by_width: If set to ``True``, the uncertainty band will be rescaled by the inverse of the width of the bin. Defaults to ``False``.
    :type rescale_by_width: bool
    :param variable: The name of the kinematic variable that is plotted on the x-axis.
    :type variable: str

    Example:

    .. code-block::

        figure_args = '''
        plot:
        xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
        yaxis: { label: r'$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
        legend: { position: 'upper center' }
        items:
            - { type: 'uncertainty-binned', label: r'$\\ell=\\mu$',
                variable: 'q2', range: [0.00, 11.63],
                datafile: './predictions-data/FF-LQCD-SSE/pred-B-to-D-mu-nu-binned'
            }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    datafile:str
    range:tuple[float, float]|None=field(default=None)
    rescale_by_width:bool=field(default=False)
    variable:str

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the binned uncertainty band for drawing.

        Loads the predictive samples from the data file and reads, for each bin, the lower and
        upper edges of the integrated kinematic ``variable``. The samples themselves are reduced
        to quantiles in :meth:`draw`.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context

        self._datafile = eos.data.Prediction(context.data_path(self.datafile))

        try:
            self._xvalues = _np.array([(p['kinematics'][self.variable + '_min'], p['kinematics'][self.variable + '_max']) for p in self._datafile.varied_parameters])
        except KeyError as e:
            raise RuntimeError(f'both \'{self.variable}_min\' and \'{self.variable}_max\' must be present in the kinematics of each prediction in the data file') from e

        if self.range is not None:
            xmin, xmax   = self.range
            self._xvalues = _np.ma.masked_outside(self._xvalues, xmin, xmax)

    def draw(self, ax):
        """Draw the binned uncertainty band on the provided axes.

        For each bin, the central 68% interval and the median are drawn as a filled box with
        outer and median lines, optionally rescaled by the inverse of the bin width.

        :param ax: The matplotlib axes onto which the band is drawn.
        :type ax: matplotlib.axes.Axes
        """

        ovalues_lower   = []
        ovalues_central = []
        ovalues_higher  = []
        for i in range(len(self._xvalues)):
            lower, central, higher = _np.quantile(self._datafile.samples[:, i], [0.15865, 0.5, 0.84135], weights=self._datafile.weights, method='inverted_cdf')
            ovalues_lower.append(lower)
            ovalues_central.append(central)
            ovalues_higher.append(higher)

        label = self.label

        for [xmin, xmax], olo, ocentral, ohi in zip(self._xvalues, ovalues_lower, ovalues_central, ovalues_higher):
            width = (xmax - xmin) if self.rescale_by_width else 1
            olo      /= width
            ocentral /= width
            ohi      /= width
            eos.debug(f"{xmin} ... {xmax} -> {ocentral} with interval {olo} .. {ohi}")
            ax.fill_between([xmin, xmax], [olo, olo], [ohi, ohi], lw=0, color=self.color, alpha=self.alpha, label=label)
            ax.plot([xmin, xmax], [olo,      olo],      color=self.color, alpha=self.alpha, lw=self.linewidth, ls=self.linestyle)
            ax.plot([xmin, xmax], [ocentral, ocentral], color=self.color, alpha=self.alpha, lw=self.linewidth, ls=self.linestyle)
            ax.plot([xmin, xmax], [ohi,      ohi],      color=self.color, alpha=self.alpha, lw=self.linewidth, ls=self.linestyle)
            label = None

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_patch()


@dataclass(kw_only=True)
class OneDimensionalHistogramItem(Item):
    """Plots a one-dimensional histogram.

    :param bins: The number of histogram bins. Defaults to 50.
    :type bins: int
    :param datafile: The path to an existing data file of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction` that contains the samples.
    :type datafile: str
    :param variable: The name of the variable that is plotted on the x-axis. Defaults to the first variable in the data file.
    :type variable: str
    :param range: The range of the variable to be plotted on the x-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type range: tuple[float, float] | None

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$|V_{cb}|$', range: [38.e-3, 47.e-3] }
          legend: { position: 'upper left' }
          items:
            - { type: 'histogram1D', variable: 'CKM::abs(V_cb)', datafile: './inference-data/CKM/samples', color: 'C0' }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    bins:int=field(default=50)
    datafile:str
    variable:str
    range:tuple[float, float]|None=field(default=None)

    _api_doc = inspect.cleandoc("""\
    Plotting One-Dimensional Histograms
    -----------------------------------

    Plot items of type ``histogram1D`` are used to display samples of a probability density, be it a prior, a posterior, or a signal PDF.
    The following key is mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.ImportanceSamples* or *eos.data.Prediction*) -- The path to
        a data file that was generated with one of the ``sample-nested`` or ``predict-observables`` tasks.

    The following keys are optional:

        * ``bins`` (*int*) -- The number of histogram bins. Defaults to 50.
        * ``variable`` (*str*) -- The name of the variable that is plotted on the x-axis. Defaults to the first variable in the data file.
        * ``range`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$|V_{cb}|$', range: [38.e-3, 47.e-3] }
          legend: { position: 'upper left' }
          items:
            - { type: 'histogram1D', variable: 'CKM::abs(V_cb)', datafile: './inference-data/CKM/samples', color: 'C0' }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()

    """)
    def __post_init__(self):
        super().__post_init__()
        if self.bins < 2:
            raise ValueError(f"Number of bins '{self.bins}' is smaller than 2")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the histogram for drawing.

        Loads the data file (of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction`),
        extracts the samples of the chosen ``variable``, and determines the histogram range if it was
        not provided explicitly.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 1D histogram")
        name = os.path.split(datafile)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(datafile)
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(datafile)
        else:
            eos.error(f"Data file '{datafile}' has an unsupported format")
            raise NotImplementedError

        if self.variable not in self._datafile.lookup_table:
            raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variable}'")

        idx = self._datafile.lookup_table[self.variable]
        self.samples = self._datafile.samples[:, idx]

        if self.range is None:
            self.range = (self.samples.min(), self.samples.max())

        pass

    def draw(self, ax):
        """Draw the one-dimensional histogram on the provided axes.

        :param ax: The matplotlib axes onto which the histogram is drawn.
        :type ax: matplotlib.axes.Axes
        """
        ax.hist(self.samples, weights=self._datafile.weights, range=self.range,
                alpha=self.alpha, bins=self.bins, color=self.color, density=True, label=self.label)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_patch()

@dataclass
class TwoDimensionalHistogramItem(Item):
    """Plots a two-dimensional histogram.

    :param datafile: The path to an existing data file of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction` that contains the histogram data.
    :type datafile: str
    :param variables: The tuple of names of the two variables that shall be histogrammed on the x- and y-axis, respectively.
    :type variables: tuple[str, str]
    :param bins: Number of histogram bins. Defaults to 50.
    :type bins: int
    :param xrange: The range of the variable to be plotted on the x-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type xrange: tuple[float, float] | None
    :param yrange: The range of the variable to be plotted on the y-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type yrange: tuple[float, float] | None
    """

    datafile:str
    variables:tuple[str, str]
    bins:int=field(default=50)
    xrange:tuple[float,float]|None=field(default=None)
    yrange:tuple[float,float]|None=field(default=None)

    _api_doc = inspect.cleandoc("""
    Plotting Two-Dimensional Histograms
    ------------------------------------

    Plot items of type ``histogram2D`` are used to display samples of a two-dimensional probability density, be it a prior, a posterior, or a prediction.

    The following keys are mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.ImportanceSamples* or *eos.data.Prediction*) -- The path to
          a data file that was generated with one of the ``sample-nested`` or ``predict-observables`` tasks.
        * ``variables`` (*tuple* of two *str*) -- The names of the two variables that are plotted on the x- and y-axis, respectively.

    The following keys are optional:
        * ``bins`` (*int*) -- The number of histogram bins. Defaults to 50.
        * ``xrange`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.
        * ``yrange`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the y-axis. Defaults to the full range of the variable in the data file.
    """)

    def __post_init__(self):
        super().__post_init__()
        if self.bins < 2:
            raise ValueError(f"Number of bins '{self.bins}' is smaller than 2")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the two-dimensional histogram for drawing.

        Loads the data file (of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction`),
        extracts the samples of the two chosen ``variables``, and determines the x- and y-ranges if they
        were not provided explicitly.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 2D histogram")
        name = os.path.split(datafile)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(datafile)
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(datafile)
        else:
            eos.error(f"Data file '{datafile}' has an unsupported format")
            raise NotImplementedError

        for v in self.variables:
            if v not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{v}'")

        xidx = self._datafile.lookup_table[self.variables[0]]
        yidx = self._datafile.lookup_table[self.variables[1]]
        self.samples = self._datafile.samples[:, (xidx, yidx)]
        self.weights = self._datafile.weights[:]
        if self.xrange is None:
            self.xrange = (self.samples[:,0].min(), self.samples[:,0].max())
        if self.yrange is None:
            self.yrange = (self.samples[:,1].min(), self.samples[:,1].max())

    def draw(self, ax):
        """Draw the two-dimensional histogram on the provided axes.

        :param ax: The matplotlib axes onto which the histogram is drawn.
        :type ax: matplotlib.axes.Axes
        """
        ax.hist2d(self.samples[:, 0], self.samples[:, 1], weights=self.weights, range=[self.xrange, self.yrange],
                  bins=self.bins, label=self.label, rasterized=True, cmap='Greys')

    def legend(self):
        # a 2D density (colormap) has no faithful single-swatch representation,
        # so it deliberately contributes no legend entry (use a colorbar instead).
        return []


@dataclass(kw_only=True)
class OneDimensionalKernelDensityEstimateItem(Item):
    """Plotting a one-dimensional kernel density estimate (KDE).

    This item type is used to display a one-dimensional kernel density estimate (KDE) of a variable, be it a prior, a posterior, or a prediction.

    :param bandwidth: Relative bandwidth factor that multiplies the automatically determined bandwidth of the KDE. Defaults to None, which means that the default bandwidth is used.
    :type bandwidth: float | None
    :param datafile: Path to the file that contains the KDE data, which must be of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction`.
    :type datafile: str
    :param level: Credibility level that shall be visualized in percent (optional). Defaults to None, which means that no credibility level is visualized.
    :type level: float | None
    :param range: The range of the variable to be plotted on the x-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type range: tuple[float, float] | None
    :param variable: Name of the datafile's variable that shall be displayed.
    :type variable: str
    :param xsamples: The number of samples to be used for the x-axis. Defaults to 100.
    :type xsamples: int

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$|V_{cb}|$', range: [38.e-3, 47.e-3] }
          legend: { position: 'upper left' }
          items:
            - { type: 'kde1D',       variable: 'CKM::abs(V_cb)', datafile: './inference-data/CKM/samples', color: 'C0',
                label: 'posterior'
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    bandwidth:float|None=field(default=None)
    datafile:str
    level:float|None=field(default=None)
    range:tuple[float, float]|None=field(default=None)
    variable:str
    xsamples:int=field(default=100)

    _api_doc = inspect.cleandoc("""
    Plotting One-Dimensional Kernel Density Estimates
    -------------------------------------------------

    Plot items of type ``kde1D`` are used to display a one-dimensional kernel density estimate (KDE) of a variable, be it a prior, a posterior, or a prediction.

    The following keys are mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.ImportanceSamples* or *eos.data.Prediction*) -- The path to
          a data file that was generated with one of the ``sample-nested`` or ``predict-observables`` tasks.
        * ``variable`` (*str*) -- The name of the variable that is plotted on the x-axis.

    The following keys are optional:

        * ``bandwidth`` (*float*) -- The relative bandwidth factor that multiplies the automatically determined bandwidth of the KDE (optional).
        * ``level`` (*float*) -- The credibility level that shall be visualized in percent (optional). Defaults to 68.3%
        * ``range`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.
        * ``xsamples`` (*int*) -- The number of samples to be used for the x-axis. Defaults to 100.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$|V_{cb}|$', range: [38.e-3, 47.e-3] }
          legend: { position: 'upper left' }
          items:
            - { type: 'kde1D',       variable: 'CKM::abs(V_cb)', datafile: './inference-data/CKM/samples', color: 'C0',
                label: 'posterior'
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()

    """)

    def __post_init__(self):
        super().__post_init__()
        if self.level is None:
            self.level = 68.3
        if self.level is not None and (self.level <= 0 or self.level >= 100):
            raise ValueError(f"Credibility level '{self.level}' is not in the interval (0, 100)")

        if self.bandwidth is not None and self.bandwidth <= 0.0:
            raise ValueError(f"Bandwidth factor '{self.bandwidth}' is not positive")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the kernel density estimate for drawing.

        Loads the data file, extracts the samples of the chosen ``variable``, fits a Gaussian KDE
        (optionally rescaling the automatically determined bandwidth), and evaluates the resulting
        probability density on a grid of ``xsamples`` points for :meth:`draw`.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 1D KDE")
        name = os.path.split(datafile)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(datafile)

            if self.variable not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variable}'")

            self.idx = self._datafile.lookup_table[self.variable]
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(datafile)
            stripped_lookup_table = { k.split(';')[0]: v for k, v in self._datafile.lookup_table.items() }
            if self.variable in stripped_lookup_table:
                if len(stripped_lookup_table.keys()) != len(self._datafile.lookup_table.keys()):
                    # variable name matches when stripping potential kinematic info from prediction variable names
                    raise ValueError(f"Data file '{datafile}' contains multiple predictions for variable '{self.variable}'; specify the full variable name including options and kinematics")
                self.idx = stripped_lookup_table[self.variable]
            else:
                if self.variable not in self._datafile.lookup_table:
                    raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variable}'")
                self.idx = self._datafile.lookup_table[self.variable]
        else:
            eos.error(f"Data file '{datafile}' has an unsupported format")
            raise NotImplementedError

        samples = self._datafile.samples[:, self.idx]
        weights = self._datafile.weights

        if self.range is None:
            self.range = (samples.min(), samples.max())

        eos.inprogress(f"Computing KDE for samples of variable '{self.variable}'")
        self.kde = _scipy.stats.gaussian_kde(samples, weights=weights)
        self.kde.set_bandwidth(bw_method='silverman')
        if self.bandwidth is not None:
            self.kde.set_bandwidth(bw_method=self.kde.factor * self.bandwidth)

        self.xvalues = _np.linspace(self.range[0], self.range[1], self.xsamples)
        self.pdf = self.kde(self.xvalues)
        self.pdf /= self.pdf.sum()

    def draw(self, ax):
        """Draw the kernel density estimate on the provided axes.

        Plots the estimated probability density and, if ``level`` is set, shades the smallest
        region containing the requested credibility level.

        :param ax: The matplotlib axes onto which the KDE is drawn.
        :type ax: matplotlib.axes.Axes
        """
        # find the PDF value corresponding to a given cumulative probability
        if self.level is not None:
            plevelf = lambda x, pdf, P: self.pdf[self.pdf > x].sum() - P
            plevel = _scipy.optimize.brentq(plevelf, 0., 1., args=(self.pdf, self.level / 100.0 ))
            ax.fill_between(_np.ma.masked_array(self.xvalues, mask=self.pdf < plevel),
                                            _np.ma.masked_array(self.pdf, mask=self.pdf < plevel, fill_value=_np.nan),
                                            facecolor=self.color, alpha=self.alpha)

        ax.plot(self.xvalues, self.pdf, color=self.color, lw=self.linewidth, ls=self.linestyle, label=self.label)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_line()

@dataclass(kw_only=True)
class TwoDimensionalKernelDensityEstimateItem(Item):
    """Plots a two-dimensional kernel density estimate (KDE).

    This item type is used to display a two-dimensional kernel density estimate (KDE) of two variables, be it a prior, a posterior, or a prediction.

    :param bandwidth: Relative bandwidth factor that multiplies the automatically determined bandwidth of the KDE. Defaults to None, which means that the default bandwidth is used.
    :type bandwidth: float | None
    :param contours: The types of contours to be drawn. Can be any combination of ``'lines'``, ``'areas'``, or ``'labels'``. Defaults to ``{'lines'}``.
    :type contours: set[str]
    :param datafile: Path to the file that contains the KDE data, which must be of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction`.
    :type datafile: str
    :param levels: The credibility levels that shall be visualized in percent. Defaults to ``[0, 68, 95]``.
    :type levels: list[float] | None
    :param variables: The tuple of names of the two variables that shall be displayed on the x- and y-axis, respectively.
    :type variables: tuple[str, str]
    :param xrange: The range of the variable to be plotted on the x-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type xrange: tuple[float, float] | None
    :param yrange: The range of the variable to be plotted on the y-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type yrange: tuple[float, float] | None

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$|V_{cb}|$', range: [38e-3, 47e-3] }
          yaxis: { label: '$f_+(0)$',   range: [0.6, 0.75]    }
          items:
            - { type: 'kde2D', label: 'posterior', color: 'C1',
                levels: [68, 95], contours: ['lines', 'areas'], bandwidth: 3.0,
                datafile: './inference-data/CKM/samples',
                variables: ['CKM::abs(V_cb)', 'B->D::alpha^f+_0@BSZ2015']
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    bandwidth:float|None=field(default=None)
    contours:set[str]=field(default_factory=lambda : {'lines'})
    datafile:str
    levels:list[float]|None=field(default=None)
    variables:tuple[str,str]
    xrange:tuple[float,float]|None=field(default=None)
    yrange:tuple[float,float]|None=field(default=None)

    _api_doc = inspect.cleandoc("""
    Plotting Two-Dimensional Kernel Density Estimates
    -------------------------------------------------

    Plot items of type ``kde2D`` are used to display a two-dimensional kernel density estimate (KDE) of two variables, be it a prior, a posterior, or a prediction.

    The following keys are mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.ImportanceSamples* or *eos.data.Prediction*) -- The path to
          a data file that was generated with one of the ``sample-nested`` or ``predict-observables`` tasks.
        * ``variables`` (*tuple* of two *str*) -- The names of the two variables that are plotted on the x- and y-axis, respectively.

    The following keys are optional:

        * ``bandwidth`` (*float*) -- The relative bandwidth factor that multiplies the automatically determined bandwidth of the KDE (optional).
        * ``contours`` (*set* of *str*) -- The types of contours to be drawn. Can be any combination of ``'lines'``, ``'areas'``, or ``'labels'``. Defaults to ``{'lines'}``.
        * ``levels`` (*list* of *float*) -- The credibility levels that shall be visualized in percent (optional). Defaults to ``[0, 68, 95]``.
        * ``xrange`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.
        * ``yrange`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the y-axis. Defaults to the full range of the variable in the data file.
    """)

    def __post_init__(self):
        super().__post_init__()

        if self.levels is None:
            self.levels = [0, 68, 95]
        elif 0 not in self.levels:
            self.levels = [0] + self.levels

        for level in self.levels:
            if level < 0 or level >= 100:
                raise ValueError(f"Credibility level '{level}' is not in the interval (0, 100)")

        if self.bandwidth is not None and self.bandwidth <= 0.0:
            raise ValueError(f"Bandwidth factor '{self.bandwidth}' is not positive")

        for contour_type in self.contours:
            if contour_type not in ['lines', 'areas', 'labels']:
                raise ValueError(f"Contour type '{contour_type}' is not supported")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the two-dimensional kernel density estimate for drawing.

        Loads the data file, extracts the samples of the two chosen ``variables``, fits a Gaussian KDE
        (optionally rescaling the automatically determined bandwidth), and evaluates the resulting
        probability density on a regular grid spanning the x- and y-ranges for :meth:`draw`.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 2D KDE")
        name = os.path.split(datafile)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(datafile)

            if self.variables[0] not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variables[0]}'")
            if self.variables[1] not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variables[1]}'")

            self._xidx = self._datafile.lookup_table[self.variables[0]]
            self._yidx = self._datafile.lookup_table[self.variables[1]]
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(datafile)

            stripped_lookup_table = { k.split(';')[0]: v for k, v in self._datafile.lookup_table.items() }

            if self.variables[0] in stripped_lookup_table:
                if len(stripped_lookup_table.keys()) != len(self._datafile.lookup_table.keys()):
                    # variable name matches when stripping potential kinematic info from prediction variable names
                    raise ValueError(f"Data file '{datafile}' contains multiple predictions for variable '{self.variables[0]}'; specify the full variable name including options and kinematics")
                self._xidx = stripped_lookup_table[self.variables[0]]
            else:
                if self.variables[0] not in self._datafile.lookup_table:
                    raise ValueError(f"Data file '{datafile}' does not contain predictions for variable '{self.variables[0]}'")
                self._xidx = self._datafile.lookup_table[self.variables[0]]

            if self.variables[1] in stripped_lookup_table:
                if len(stripped_lookup_table.keys()) != len(self._datafile.lookup_table.keys()):
                    # variable name matches when stripping potential kinematic info from prediction variable names
                    raise ValueError(f"Data file '{datafile}' contains multiple predictions for variable '{self.variables[1]}'; specify the full variable name including options and kinematics")
                self._yidx = stripped_lookup_table[self.variables[1]]
            else:
                if self.variables[1] not in self._datafile.lookup_table:
                    raise ValueError(f"Data file '{datafile}' does not contain predictions for variable '{self.variables[1]}'")
                self._yidx = self._datafile.lookup_table[self.variables[1]]
        else:
            eos.error(f"Data file '{datafile}' has an unsupported format")
            raise NotImplementedError

        samples = self._datafile.samples[:, (self._xidx, self._yidx)]
        weights = self._datafile.weights

        eos.inprogress(f"Computing KDE for samples of variables '{self.variables[0]}' and '{self.variables[1]}'")
        self._kde = _scipy.stats.gaussian_kde(samples.T, weights=weights)
        self._kde.set_bandwidth(bw_method='silverman')
        if self.bandwidth is not None:
            self._kde.set_bandwidth(bw_method=self._kde.factor * self.bandwidth)

        # determine the extent of the plot
        if self.xrange is None:
            self.xrange = (samples[:, 0].min(), samples[:, 0].max())
        if self.yrange is None:
            self.yrange = (samples[:, 1].min(), samples[:, 1].max())

        # compute the PDF on a grid
        xx,yy = _np.mgrid[self.xrange[0]:self.xrange[1]:100j, self.yrange[0]:self.yrange[1]:100j]
        self._positions = _np.vstack([xx.ravel(), yy.ravel()])
        self._pdf = _np.reshape(self._kde(self._positions).T, xx.shape)
        self._pdf /= self._pdf.sum()

    def _plevels(self):
        """Return the PDF threshold values corresponding to the requested credibility ``levels``.

        Each threshold is the density value above which the requested fraction of the total
        probability lies. The 0% (peak) level is handled explicitly as the maximum density, since
        solving for it numerically would return the upper bracket (~1.0) and produce an out-of-range
        contour for typical PDFs whose maximum is much smaller than one.

        :returns: The threshold values in the same order as ``levels``.
        :rtype: list[float]
        """
        # find the PDF value corresponding to a given cummulative probability
        plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
        plevels = []
        for level in self.levels:
            if level == 0:
                plevels.append(self._pdf.max())
            else:
                plevels.append(_scipy.optimize.brentq(plevel, 0., 1., args=(self._pdf, level / 100.0)))
        return plevels

    def draw(self, ax):
        """Draw the two-dimensional kernel density estimate on the provided axes.

        Draws contour lines at the requested credibility ``levels`` and, depending on the
        ``contours`` setting, optionally fills the contour areas and/or labels the contour lines.

        :param ax: The matplotlib axes onto which the KDE is drawn.
        :type ax: matplotlib.axes.Axes
        """
        plevels = self._plevels()
        labels = [f'{level}%' for level in self.levels]

        if 'areas' in self.contours:
            colors = [_matplotlib.colors.to_rgba(self.color, alpha) for alpha in _np.linspace(0.50, 1.00, len(self.levels))]
            ax.contourf(self._pdf.transpose(),
                        colors=colors,
                        extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                        levels=plevels[::-1])

        CS = ax.contour(self._pdf.transpose(),
                        colors=self.color,
                        extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                        levels=plevels[::-1],
                        linestyles=self.linestyle)

        if 'labels' in self.contours:
            fmt = {}
            for level, label in zip(CS.levels, labels[::-1]):
                fmt[level] = label

            ax.clabel(CS, inline=1, fmt=fmt, fontsize=10)

    def legend(self):
        """Return the item's legend entry as a list of handle/label pairs.

        Provides a filled rectangle handle when contour areas are drawn, and a line handle otherwise.

        :returns: A list containing a single ``(handle, label)`` pair if a label is set, otherwise empty.
        :rtype: list[tuple[matplotlib.artist.Artist, str]]
        """
        entries = []

        if self.label:
            handle = None
            if 'areas' in self.contours:
                handle = _matplotlib.pyplot.Rectangle((0,0),1,1, color=self.color)
            else:
                handle = _matplotlib.pyplot.Line2D((0,1),(0.5,0.), color=self.color, linestyle=self.linestyle)

            entries.append((handle, self.label))

        return entries

@dataclass(kw_only=True)
class TwoDimensionalContoursItem(Item):
    """Plots two-dimensional probability contours from a histogram of pre-existing samples.

    This item type is used to display the credibility contours of two variables, be it a prior, a
    posterior, or a prediction. In contrast to :class:`TwoDimensionalKernelDensityEstimateItem`, the
    probability density is estimated from a two-dimensional histogram of the samples rather than from
    a kernel density estimate, so that no smoothing is applied.

    :param bins: The number of histogram bins along each axis used to estimate the density. Defaults to 100.
    :type bins: int
    :param contours: The types of contours to be drawn. Can be any combination of ``'lines'``, ``'areas'``, or ``'labels'``. Defaults to ``{'lines'}``.
    :type contours: set[str]
    :param datafile: Path to the file that contains the samples, which must be of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction`.
    :type datafile: str
    :param levels: The credibility levels that shall be visualized in percent. Defaults to ``[68, 95, 99]``.
    :type levels: list[float] | None
    :param variables: The tuple of names of the two variables that shall be displayed on the x- and y-axis, respectively.
    :type variables: tuple[str, str]
    :param xrange: The range of the variable to be plotted on the x-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type xrange: tuple[float, float] | None
    :param yrange: The range of the variable to be plotted on the y-axis, given as a tuple of two float values (min, max). Defaults to the full range of the variable in the data file.
    :type yrange: tuple[float, float] | None

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$|V_{cb}|$', range: [38e-3, 47e-3] }
          yaxis: { label: '$f_+(0)$',   range: [0.6, 0.75]    }
          items:
            - { type: 'contours2D', label: 'posterior', color: 'C1',
                levels: [68, 95, 99], contours: ['lines', 'labels'],
                datafile: './inference-data/CKM/samples',
                variables: ['CKM::abs(V_cb)', 'B->D::alpha^f+_0@BSZ2015']
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    bins:int=field(default=100)
    contours:set[str]=field(default_factory=lambda : {'lines'})
    datafile:str
    levels:list[float]|None=field(default=None)
    variables:tuple[str,str]
    xrange:tuple[float,float]|None=field(default=None)
    yrange:tuple[float,float]|None=field(default=None)

    _api_doc = inspect.cleandoc("""
    Plotting Two-Dimensional Contours
    ---------------------------------

    Plot items of type ``contours2D`` are used to display the credibility contours of a two-dimensional
    probability density, be it a prior, a posterior, or a prediction. The density is estimated from a
    two-dimensional histogram of the samples, i.e. without the smoothing applied by a kernel density estimate.

    The following keys are mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.ImportanceSamples* or *eos.data.Prediction*) -- The path to
          a data file that was generated with one of the ``sample-nested`` or ``predict-observables`` tasks.
        * ``variables`` (*tuple* of two *str*) -- The names of the two variables that are plotted on the x- and y-axis, respectively.

    The following keys are optional:

        * ``bins`` (*int*) -- The number of histogram bins along each axis used to estimate the density. Defaults to 100.
        * ``contours`` (*set* of *str*) -- The types of contours to be drawn. Can be any combination of ``'lines'``, ``'areas'``, or ``'labels'``. Defaults to ``{'lines'}``.
        * ``levels`` (*list* of *float*) -- The credibility levels that shall be visualized in percent (optional). Defaults to ``[68, 95, 99]``.
        * ``xrange`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.
        * ``yrange`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the y-axis. Defaults to the full range of the variable in the data file.
    """)

    def __post_init__(self):
        super().__post_init__()

        if self.bins < 2:
            raise ValueError(f"Number of bins '{self.bins}' is smaller than 2")

        if self.levels is None:
            self.levels = [0, 68, 95, 99]
        elif 0 not in self.levels:
            # the 0% level (the peak of the density) is required as the innermost boundary so that
            # the credibility areas can be filled correctly
            self.levels = [0] + self.levels

        for level in self.levels:
            if level < 0 or level >= 100:
                raise ValueError(f"Credibility level '{level}' is not in the interval (0, 100)")

        for contour_type in self.contours:
            if contour_type not in ['lines', 'areas', 'labels']:
                raise ValueError(f"Contour type '{contour_type}' is not supported")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the two-dimensional contours for drawing.

        Loads the data file (of type :class:`eos.data.ImportanceSamples` or :class:`eos.data.Prediction`),
        extracts the samples of the two chosen ``variables``, and estimates the probability density on a
        regular grid from a weighted two-dimensional histogram spanning the x- and y-ranges for :meth:`draw`.

        :param context: The analysis file context used to resolve the relative path to ``datafile``.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 2D contours")
        name = os.path.split(datafile)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(datafile)

            if self.variables[0] not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variables[0]}'")
            if self.variables[1] not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variables[1]}'")

            self._xidx = self._datafile.lookup_table[self.variables[0]]
            self._yidx = self._datafile.lookup_table[self.variables[1]]
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(datafile)

            stripped_lookup_table = { k.split(';')[0]: v for k, v in self._datafile.lookup_table.items() }

            if self.variables[0] in stripped_lookup_table:
                if len(stripped_lookup_table.keys()) != len(self._datafile.lookup_table.keys()):
                    # variable name matches when stripping potential kinematic info from prediction variable names
                    raise ValueError(f"Data file '{datafile}' contains multiple predictions for variable '{self.variables[0]}'; specify the full variable name including options and kinematics")
                self._xidx = stripped_lookup_table[self.variables[0]]
            else:
                if self.variables[0] not in self._datafile.lookup_table:
                    raise ValueError(f"Data file '{datafile}' does not contain predictions for variable '{self.variables[0]}'")
                self._xidx = self._datafile.lookup_table[self.variables[0]]

            if self.variables[1] in stripped_lookup_table:
                if len(stripped_lookup_table.keys()) != len(self._datafile.lookup_table.keys()):
                    # variable name matches when stripping potential kinematic info from prediction variable names
                    raise ValueError(f"Data file '{datafile}' contains multiple predictions for variable '{self.variables[1]}'; specify the full variable name including options and kinematics")
                self._yidx = stripped_lookup_table[self.variables[1]]
            else:
                if self.variables[1] not in self._datafile.lookup_table:
                    raise ValueError(f"Data file '{datafile}' does not contain predictions for variable '{self.variables[1]}'")
                self._yidx = self._datafile.lookup_table[self.variables[1]]
        else:
            eos.error(f"Data file '{datafile}' has an unsupported format")
            raise NotImplementedError

        samples = self._datafile.samples[:, (self._xidx, self._yidx)]
        weights = self._datafile.weights

        # determine the extent of the plot
        if self.xrange is None:
            self.xrange = (samples[:, 0].min(), samples[:, 0].max())
        if self.yrange is None:
            self.yrange = (samples[:, 1].min(), samples[:, 1].max())

        # estimate the PDF from a weighted 2D histogram; ``self._pdf`` holds the probability mass
        # per bin and sums to one, matching the convention used by the KDE-based contour item
        H, _, _ = _np.histogram2d(samples[:, 0], samples[:, 1], bins=self.bins,
                                  range=[self.xrange, self.yrange], weights=weights)
        self._pdf = H / H.sum()

    def _plevels(self):
        """Return the PDF threshold values corresponding to the requested credibility ``levels``.

        Each threshold is the density value above which the requested fraction of the total
        probability lies. The 0% (peak) level is handled explicitly as the maximum density, since
        solving for it numerically would return the upper bracket (~1.0) and produce an out-of-range
        contour for typical PDFs whose maximum is much smaller than one.

        :returns: The threshold values in the same order as ``levels``.
        :rtype: list[float]
        """
        # find the PDF value corresponding to a given cummulative probability
        plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
        plevels = []
        for level in self.levels:
            if level == 0:
                plevels.append(self._pdf.max())
            else:
                plevels.append(_scipy.optimize.brentq(plevel, 0., 1., args=(self._pdf, level / 100.0)))
        return plevels

    def draw(self, ax):
        """Draw the two-dimensional contours on the provided axes.

        Draws contour lines at the requested credibility ``levels`` and, depending on the
        ``contours`` setting, optionally fills the contour areas and/or labels the contour lines.

        :param ax: The matplotlib axes onto which the contours are drawn.
        :type ax: matplotlib.axes.Axes
        """
        plevels = self._plevels()
        labels = [f'{level}%' for level in self.levels]

        if 'areas' in self.contours:
            colors = [_matplotlib.colors.to_rgba(self.color, alpha) for alpha in _np.linspace(0.50, 1.00, len(self.levels))]
            ax.contourf(self._pdf.transpose(),
                        colors=colors,
                        extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                        levels=plevels[::-1])

        CS = ax.contour(self._pdf.transpose(),
                        colors=self.color,
                        extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                        levels=plevels[::-1],
                        linestyles=self.linestyle)

        if 'labels' in self.contours:
            fmt = {}
            for level, label in zip(CS.levels, labels[::-1]):
                fmt[level] = label

            ax.clabel(CS, inline=1, fmt=fmt, fontsize=10)

    def legend(self):
        """Return the item's legend entry as a list of handle/label pairs.

        Provides a filled rectangle handle when contour areas are drawn, and a line handle otherwise.

        :returns: A list containing a single ``(handle, label)`` pair if a label is set, otherwise empty.
        :rtype: list[tuple[matplotlib.artist.Artist, str]]
        """
        entries = []

        if self.label:
            handle = None
            if 'areas' in self.contours:
                handle = _matplotlib.pyplot.Rectangle((0,0),1,1, color=self.color)
            else:
                handle = _matplotlib.pyplot.Line2D((0,1),(0.5,0.), color=self.color, linestyle=self.linestyle)

            entries.append((handle, self.label))

        return entries

@dataclass(kw_only=True)
class ConstraintItem(Item):
    """Plots statistical constraints from the EOS library of experimental and theoretical likelihoods.

    :param constraints: The name or the list of names of the constraints that will be plotted. Must identify at least one of the constraints known to EOS; see `the complete list of constraints <../reference/constraints.html>`_.
    :type constraints: eos.QualifiedName | list[eos.QualifiedName]
    :param variable: The name of the kinematic variable to which the x axis will be mapped.
    :type variable: str
    :param observable: The name of the observable whose constraints will be plotted. Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_. This is only mandatory in multivariate constraints, since these can constrain more than one observable simultaneously.
    :type observable: eos.QualifiedName | None
    :param range: The interval in which the observable is plotted in the case of a multivariate constraint.
    :type range: tuple[int, int] | None
    :param rescale_by_width: Rescales binned constraints by the inverse of the bin width. This is often required to compare theory (integrated) predictions and experimental (averaged) measurements. Defaults to false.
    :type rescale_by_width: bool

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                  range: [0.0,  5e-3] }
          legend: { position: 'lower left' }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2;l=e,q=d', label: r'$\\ell=e$',
                variable: 'q2', range: [ 0.02, 11.63 ], color: 'black'
              }
            - { type: 'constraint', 'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',  observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=e,\\, q=d$',
                variable: 'q2', rescale_by_width: true
              }
            - { type: 'constraint', 'constraints': 'B^0->D^+mu^-nu::BRs@Belle:2015A', observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=\\mu,\\, q=d$',
                variable: 'q2', rescale_by_width: true
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()

    """

    constraints:eos.QualifiedName|list[eos.QualifiedName]
    variable:str
    observable:eos.QualifiedName|None=None
    range:tuple[int, int]|None=field(default=None)
    rescale_by_width:bool=False

    _api_doc = inspect.cleandoc("""\
    Plotting Constraints
    --------------------

    Plot items of type ``constraints`` are used to display one of the built-in `experimental or theoretical constraints <../reference/constraints.html>`_.
    The following keys are mandatory:

        * ``constraints`` (:class:`QualifiedName <eos.QualifiedName>` or iterable thereof) -- The name or the list of names of the constraints
        that will be plotted. Must identify at least one of the constraints known to EOS; see `the complete list of constraints <../reference/constraints.html>`_.
        * ``variable`` (*str*) -- The name of the kinematic variable to which the x axis will be mapped.

    When plotting multivariate constraints, the following key is also mandatory:

        * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable whose constraints will be plotted.
        Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_.
        This is only mandatory in multivariate constraints, since these can constrain more than one observable simultaneously.

    The following keys are optional:

        * ``range`` (list of int) -- The interval in which the observable is plotted in the case of a multivariate constraint.
        * ``rescale_by_width`` (*bool*) -- Rescales binned constraints by the inverse of the bin width. This is often required
        to compare theory (integrated) predictions and experimental (averaged) measurements. Defaults to false.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
          legend: { position: 'lower left' }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2;l=e,q=d', label: r'$\\ell=e$',
                variable: 'q2', range: { min: 0.02, max: 11.63, num: 100}, color: 'black'
              }
            - { type: 'constraint', 'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',  observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=e,\\, q=d$',
                variable: 'q2', rescale_by_width: true
              }
            - { type: 'constraint', 'constraints': 'B^0->D^+mu^-nu::BRs@Belle:2015A', observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=\\mu,\\, q=d$',
                variable: 'q2', rescale_by_width: true
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()

    """)

    def __post_init__(self):
        super().__post_init__()

        if type(self.constraints) == str:
            self.constraints = [eos.QualifiedName(self.constraints)]
        elif type(self.constraints) == eos.QualifiedName:
            self.constraints = [self.constraints]
        elif type(self.constraints) is not list:
            raise TypeError(f'constraints must be a QualifiedName or a list of QualifiedNames, not {type(self.constraints)}')

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the constraint for drawing.

        Looks up and deserializes each named constraint, then extracts the central values and
        (a)symmetric uncertainties of the chosen ``observable`` as a function of the kinematic
        ``variable``, optionally rescaling binned constraints by the inverse of the bin width.
        The supported constraint types are ``Gaussian``, ``MultivariateGaussian(Covariance)``,
        and ``MultivariateGaussian``.

        :param context: The analysis file context. Accepted for interface consistency; this item
            does not read any data files. If ``None``, a default context is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context

        import yaml

        constraints = eos.Constraints()

        self._xvalues = []
        self._xerrors = []
        self._yvalues = []
        self._yerrors = []

        for constraint_name in self.constraints:
            entry = constraints[constraint_name]
            if not entry:
                raise ValueError(f'unknown constraint {constraint_name}')

            constraint = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)

            xvalues = None
            xerrors = None
            yvalues = None
            yerrors = None

            if constraint['type'] == 'Gaussian':
                kinematics = constraint['kinematics']
                width = 1
                if self.variable in kinematics:
                    xvalues = [kinematics[self.variable]]
                    xerrors = [None]

                elif (self.variable + '_min' in kinematics) and (self.variable + '_max' in kinematics):
                    var_max = float(kinematics[self.variable + '_max'])
                    var_min = float(kinematics[self.variable + '_min'])
                    xvalues = [(var_max + var_min) / 2]
                    xerrors = [(var_max - var_min) / 2]
                    if self.rescale_by_width:
                        width = (var_max - var_min)

                yvalues = [float(constraint['mean']) / width]
                sigma_hi = _np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2) / width
                sigma_lo = _np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2) / width
                yerrors = [(sigma_hi, sigma_lo)]
            elif constraint['type'] == 'MultivariateGaussian(Covariance)':
                if not self.observable:
                    raise KeyError('observable needs to be specified for MultivariateGaussian(Covariance) constraints')

                covariance = _np.array(constraint['covariance'])
                observables = constraint['observables']
                means = constraint['means']
                options = constraint['options']
                dim = len(means)
                kinematics = constraint['kinematics']

                xvalues = []
                xerrors = []
                yvalues = []
                yerrors = []
                for i in range(0, dim):
                    width = 1

                    # Check that the observable match and that the provided options match with those of the constraint
                    if not (observables[i] == eos.QualifiedName(self.observable).full().split(';')[0] and
                               _np.all([eos.Options(options[i])[k] == v for k, v in eos.QualifiedName(self.observable).options_part()])):
                        eos.warn(f'    skipping the observable {observables[i]}, with options {options[i]}, because of name or option mismatch')
                        continue
                    _kinematics = kinematics[i]
                    if self.variable in _kinematics:
                        xvalues.append(_kinematics[self.variable])
                        xerrors.append(None)

                    elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                        var_max = float(_kinematics[self.variable + '_max'])
                        var_min = float(_kinematics[self.variable + '_min'])
                        xvalues.append((var_max + var_min) / 2)
                        xerrors.append((var_max - var_min) / 2)
                        if self.rescale_by_width:
                            width = (var_max - var_min)

                    yvalues.append(_np.double(means[i]) / width)
                    yerrors.append(_np.sqrt(_np.double(covariance[i, i])) / width)
            elif constraint['type'] == 'MultivariateGaussian':
                if not self.observable:
                    raise KeyError('observable needs to be specified for MultivariateGaussian constraints')
                sigma_stat_hi = _np.array(constraint['sigma-stat-hi'])
                sigma_stat_lo = _np.array(constraint['sigma-stat-lo'])
                sigma_sys = _np.array(constraint['sigma-sys'])
                sigma = _np.sqrt(_np.power(sigma_sys, 2) + 0.25 * _np.power(sigma_stat_hi + sigma_stat_lo, 2))
                observables = constraint['observables']
                means = constraint['means']
                dim = len(means)
                kinematics = constraint['kinematics']

                xvalues = []
                xerrors = []
                yvalues = []
                yerrors = []
                for i in range(0, dim):
                    width = 1

                    if not (observables[i] == eos.QualifiedName(self.observable).full().split(';')[0] and
                               _np.all([eos.Options(options[i])[k] == v for k, v in eos.QualifiedName(self.observable).options_part()])):
                        eos.warn('    skipping the observable {observables[i]}, with options {options[i]}, because of name or option mismatch')
                        continue
                    _kinematics = kinematics[i]
                    if self.variable in _kinematics:
                        xvalues.append(_kinematics[self.variable])
                        xerrors.append(None)

                    elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                        var_max = float(_kinematics[self.variable + '_max'])
                        var_min = float(_kinematics[self.variable + '_min'])
                        xvalues.append((var_max + var_min) / 2)
                        xerrors.append((var_max - var_min) / 2)
                        if self.rescale_by_width:
                            width = (var_max - var_min)

                    yvalues.append(means[i] / width)
                    yerrors.append(sigma[i] / width)
            else:
                raise ValueError(f'constraint type {constraint["type"]} presently not supported')

        self._xvalues = _np.array(xvalues)
        self._xerrors = _np.array(xerrors)
        self._yvalues = _np.array(yvalues)
        self._yerrors = _np.array(yerrors)

        if len(xvalues) == 0:
            eos.info(f'   skipping plot for constraint {constraint_name} since it does not contain the requested observable')
            return

        if self.range:
            self._mask = _np.logical_and(self._xvalues > min(self.range), self._xvalues < max(self.range))
        else:
            self._mask = _np.array([True] * len(self._xvalues))

    def draw(self, ax):
        """Draw the constraint on the provided axes.

        Renders each data point as an error bar, using horizontal error bars to indicate the bin
        width of binned constraints.

        :param ax: The matplotlib axes onto which the constraint is drawn.
        :type ax: matplotlib.axes.Axes
        """

        if len(self._xvalues) == 0:
            return

        xvalues = self._xvalues[self._mask]
        yvalues = self._yvalues[self._mask]
        xerrors = self._xerrors[self._mask]
        yerrors = self._yerrors[self._mask]

        label = self.label
        for xv, xerr, yv, yerr in zip(xvalues, xerrors, yvalues, yerrors):
            if xerr is not None:
                ax.errorbar(x=xv, y=yv, xerr=xerr, yerr=yerr, color=self.color,
                                         elinewidth=1.0, fmt='_', linestyle='none', label=label)
            else:
                ax.errorbar(x=xv, y=yv, yerr=yerr, color=self.color,
                                         elinewidth=1.0, fmt='_', linestyle='none', label=label)

            # disable the label for subsequent plots
            label = None

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        # the constraint is drawn as capped error bars (always a y error, plus an x error
        # for binned constraints), so the swatch must match. The x error is only known
        # after prepare(); fall back to a y-only swatch if the data is not yet available.
        xerrors = getattr(self, '_xerrors', None)
        has_xerr = xerrors is not None and any(xe is not None for xe in xerrors)
        return self._legend_errorbar(has_xerr=has_xerr, has_yerr=True)


@dataclass(kw_only=True)
class TwoDimensionalConstraintItem(Item):
    r"""Plots the 2D contours of two correlated observables from a single constraint.

    This item type displays the uncertainty region of a single constraint from the EOS library in
    the plane of two of its observables. For a bivariate (``MultivariateGaussian`` or
    ``MultivariateGaussian(Covariance)``) constraint, the chosen pair of observables is drawn as one
    or more covariance ellipses, one per requested confidence level. For a univariate ``Gaussian``
    constraint, the single observable is drawn as a band spanning the orthogonal axis; in that case
    exactly one of ``x`` or ``y`` is given.

    Each of ``x`` and ``y`` is a dictionary with a single mandatory key ``observable`` that names the
    observable mapped to that axis. For the time being, the observable is matched against the
    constraint by its qualified name only.

    :param constraint: The name of the constraint to be plotted. Must identify one of the constraints known to EOS; see `the complete list of constraints <../reference/constraints.html>`_.
    :type constraint: eos.QualifiedName
    :param x: The specification of the observable mapped to the x axis, given as a dictionary with the mandatory key ``observable``.
    :type x: dict | None
    :param y: The specification of the observable mapped to the y axis, given as a dictionary with the mandatory key ``observable``.
    :type y: dict | None
    :param sigmas: The list of confidence levels (in multiples of the standard deviation) for which a contour is drawn. Defaults to ``[1.0]``.
    :type sigmas: list[float]

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$|V_{ub}|$' }
          yaxis: { label: r'$|V_{cb}|$' }
          items:
            - { type: 'constraint2D', constraint: 'B->pilnu::|V_ub|@HFLAV:2024A',
                x: { observable: 'CKM::abs(V_ub)' }, y: { observable: 'CKM::abs(V_cb)' },
                sigmas: [1.0, 2.0], color: 'C0', label: 'HFLAV 2024'
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    constraint:eos.QualifiedName
    x:dict|None=field(default=None)
    y:dict|None=field(default=None)
    sigmas:list[float]=field(default_factory=lambda: [1.0])

    def __post_init__(self):
        super().__post_init__()

        if type(self.constraint) is str:
            self.constraint = eos.QualifiedName(self.constraint)
        elif type(self.constraint) is not eos.QualifiedName:
            raise TypeError(f'constraint must be a QualifiedName, not {type(self.constraint)}')

        if self.x is None and self.y is None:
            raise ValueError("at least one of 'x' or 'y' must be specified for a 2D constraint")

        for axis_name, spec in (('x', self.x), ('y', self.y)):
            if spec is not None and (not isinstance(spec, dict) or 'observable' not in spec):
                raise ValueError(f"'{axis_name}' must be a dictionary containing an 'observable' key")

        # draw the largest contour first (lowest opacity) and the smallest last (highest opacity),
        # so that nested contours remain visible
        self.sigmas = sorted(self.sigmas, reverse=True)
        self._alphas = _np.linspace(0.0, self.alpha, len(self.sigmas) + 1)[1:]

    @staticmethod
    def _name(observable):
        "Return the bare observable name (without options) used to match an axis observable."
        return eos.QualifiedName(observable).full().split(';')[0]

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the 2D constraint for drawing.

        Looks up and deserializes the named constraint and caches the geometry of its uncertainty
        region: the means, standard deviations, and correlation of the chosen pair of observables for
        a bivariate constraint, or the mean and (a)symmetric uncertainties of the single observable
        for a univariate ``Gaussian`` constraint. The supported constraint types are ``Gaussian``,
        ``MultivariateGaussian(Covariance)``, and ``MultivariateGaussian``.

        :param context: The analysis file context. Accepted for interface consistency; this item
            does not read any data files. If ``None``, a default context is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context

        import yaml

        entry = eos.Constraints()[self.constraint]
        if not entry:
            raise ValueError(f'unknown constraint {self.constraint}')

        constraint = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)
        ctype = constraint['type']
        if ctype not in ('Gaussian', 'MultivariateGaussian', 'MultivariateGaussian(Covariance)'):
            raise ValueError(f"constraint type '{ctype}' presently not supported")

        if ctype in ('MultivariateGaussian(Covariance)', 'MultivariateGaussian'):
            if self.x is None or self.y is None:
                raise ValueError("both 'x' and 'y' must be specified for a multivariate constraint")

            # for the time being, observables are matched by their qualified name only,
            # ignoring any kinematics or options given in the axis specifications
            observables = [self._name(o) for o in constraint['observables']]
            xname = self._name(self.x['observable'])
            yname = self._name(self.y['observable'])
            if xname not in observables:
                raise ValueError(f"x-axis observable {self.x['observable']} not contained in constraint {self.constraint}")
            if yname not in observables:
                raise ValueError(f"y-axis observable {self.y['observable']} not contained in constraint {self.constraint}")
            xidx = observables.index(xname)
            yidx = observables.index(yname)

            means = _np.array(constraint['means'])
            self._xmean = means[xidx]
            self._ymean = means[yidx]

            if ctype == 'MultivariateGaussian(Covariance)':
                cov = _np.array(constraint['covariance'])
                self._xsigma = _np.sqrt(cov[xidx, xidx])
                self._ysigma = _np.sqrt(cov[yidx, yidx])
                self._rho    = cov[xidx, yidx] / (self._xsigma * self._ysigma)
            else:
                sigmas = _np.sqrt(
                      (_np.abs(_np.array(constraint['sigma-stat-hi'])) + _np.abs(_np.array(constraint['sigma-stat-lo'])))**2 / 4.0
                    +  _np.array(constraint['sigma-sys'])**2
                )
                self._xsigma = sigmas[xidx]
                self._ysigma = sigmas[yidx]
                self._rho    = _np.array(constraint['correlations'])[xidx, yidx]

            self._shape = 'ellipse'
        else: # Gaussian
            if self.x is not None and self.y is not None:
                raise ValueError("only one of 'x' or 'y' may be specified for a univariate (Gaussian) constraint")

            axis = 'x' if self.x is not None else 'y'
            spec = self.x if self.x is not None else self.y

            # for the time being, the observable is matched by its qualified name only
            if self._name(spec['observable']) != self._name(constraint['observable']):
                raise ValueError(f"{axis}-axis observable {spec['observable']} not contained in constraint {self.constraint}")

            self._mean     = float(constraint['mean'])
            self._sigma_hi = _np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2)
            self._sigma_lo = _np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2)
            self._shape    = 'rect-x' if axis == 'x' else 'rect-y'

    def draw(self, ax):
        """Draw the 2D constraint on the provided axes.

        Renders a covariance ellipse for each requested confidence level for a bivariate constraint,
        or a band spanning the orthogonal axis for a univariate ``Gaussian`` constraint.

        :param ax: The matplotlib axes onto which the 2D constraint is drawn.
        :type ax: matplotlib.axes.Axes
        """
        if self._shape == 'ellipse':
            # the ellipse is built in a frame in which the principal axes are at 45 degrees, then
            # scaled by the standard deviations and rotated/translated into the data frame
            xwidth = 2.0 * _np.sqrt(1.0 + self._rho)
            ywidth = 2.0 * _np.sqrt(1.0 - self._rho)
            for sigma, alpha in zip(self.sigmas, self._alphas):
                ellipse = _matplotlib.patches.Ellipse((0.0, 0.0), width=xwidth, height=ywidth,
                                                       alpha=alpha, color=self.color, linewidth=self.linewidth, fill=True)
                transf = _matplotlib.transforms.Affine2D() \
                    .rotate_deg(45) \
                    .scale(self._xsigma * sigma, self._ysigma * sigma) \
                    .translate(self._xmean, self._ymean)
                ellipse.set_transform(transf + ax.transData)
                ax.add_patch(ellipse)
        elif self._shape == 'rect-x':
            # span the full orthogonal (y) axis: x is in data coordinates, y in axes coordinates,
            # so the band is independent of the y limits at the time this item is drawn
            transform = ax.get_xaxis_transform()
            for sigma, alpha in zip(self.sigmas, self._alphas):
                xlo = self._mean - sigma * self._sigma_lo
                xhi = self._mean + sigma * self._sigma_hi
                ax.add_patch(_matplotlib.patches.Rectangle((xlo, 0), width=xhi - xlo, height=1,
                                                           transform=transform, alpha=alpha, color=self.color))
        elif self._shape == 'rect-y':
            # span the full orthogonal (x) axis: y is in data coordinates, x in axes coordinates,
            # so the band is independent of the x limits at the time this item is drawn
            transform = ax.get_yaxis_transform()
            for sigma, alpha in zip(self.sigmas, self._alphas):
                ylo = self._mean - sigma * self._sigma_lo
                yhi = self._mean + sigma * self._sigma_hi
                ax.add_patch(_matplotlib.patches.Rectangle((0, ylo), width=1, height=yhi - ylo,
                                                           transform=transform, alpha=alpha, color=self.color))

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_patch()


@dataclass(kw_only=True)
class ConstraintResidueItem(Item):
    """Plots the residues of a statistical constraints from the EOS library of experimental and theoretical likelihoods.

    :param constraints: The name or the list of names of the constraints that will be plotted. Must identify at least one of the constraints known to EOS; see `the complete list of constraints <../reference/constraints.html>`_.
    :type constraints: eos.QualifiedName | list[eos.QualifiedName]
    :param variable: The name of the kinematic variable to which the x axis will be mapped.
    :type variable: str
    :param observable: The name of the observable whose residues will be plotted. Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_. This is only mandatory in multivariate constraints, since these can constrain more than one observable simultaneously.
    :type observable: eos.QualifiedName | None
    :param range: The interval in which the observable is plotted in the case of a multivariate constraint.
    :type range: tuple[int, int] | None
    :param parameters: The set of parameters used to evaluate the constraint's observables, given as a dictionary mapping parameter names to their values.
    :type parameters: dict[eos.QualifiedName, float] | None
    :param rescale_by_width: Rescales binned constraints by the inverse of the bin width. This is often required to compare theory (integrated) predictions and experimental (averaged) measurements. Defaults to false.
    :type rescale_by_width: bool

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                  range: [0.0,  5e-3] }
          legend: { position: 'lower left' }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2;l=e,q=d', label: r'$\\ell=e$',
                variable: 'q2', range: [ 0.02, 11.63 ], color: 'black'
              }
            - { type: 'constraint-residue', 'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',  observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=e,\\, q=d$',
                variable: 'q2', 'parameters': {"mass::e": 1.0}, rescale_by_width: true
              }
            - { type: 'constraint-residue', 'constraints': 'B^0->D^+mu^-nu::BRs@Belle:2015A', observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=\\mu,\\, q=d$',
                variable: 'q2', 'parameters': {"mass::e": 1.0}, rescale_by_width: true
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()

    """

    constraints:eos.QualifiedName|list[eos.QualifiedName]
    variable:str
    observable:eos.QualifiedName|None=None
    range:tuple[int, int]|None=field(default=None)
    parameters:dict[eos.QualifiedName,float]|None=field(default=None)
    rescale_by_width:bool=False

    _api_doc = inspect.cleandoc("""\
    Plotting Constraints
    --------------------

    Plot items of type ``constraints`` are used to display one of the built-in `experimental or theoretical constraints <../reference/constraints.html>`_.
    The following keys are mandatory:

        * ``constraints`` (:class:`QualifiedName <eos.QualifiedName>` or iterable thereof) -- The name or the list of names of the constraints
        that will be plotted. Must identify at least one of the constraints known to EOS; see `the complete list of constraints <../reference/constraints.html>`_.
        * ``variable`` (*str*) -- The name of the kinematic variable to which the x axis will be mapped.

    When plotting multivariate constraints, the following key is also mandatory:

        * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable whose constraints will be plotted.
        Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_.
        This is only mandatory in multivariate constraints, since these can constrain more than one observable simultaneously.

    The following keys are optional:

        * ``range`` (list of int) -- The interval in which the observable is plotted in the case of a multivariate constraint.
        * ``parameters`` (dict of [eos.QualifiedName, float]) -- The set of parameters used to evaluate the constraint's observables,
        given as a dictionary mapping parameter names to their values. If None, EOS default parameters will be used.
        * ``rescale_by_width`` (*bool*) -- Rescales binned constraints by the inverse of the bin width. This is often required
        to compare theory (integrated) predictions and experimental (averaged) measurements. Defaults to false.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$ residues',         range: [-2e-3,  2e-3] }
          legend: { position: 'lower left' }
          items:
            - { type: 'constraint-residue', 'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',  observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=e,\\, q=d$',
                variable: 'q2', 'parameters': {"mass::e": 1.0}, rescale_by_width: true
              }
            - { type: 'constraint-residue', 'constraints': 'B^0->D^+mu^-nu::BRs@Belle:2015A', observable: 'B->Dlnu::BR', label: r'Belle 2015 $\\ell=\\mu,\\, q=d$',
                variable: 'q2', 'parameters': {"mass::e": 1.0}, rescale_by_width: true
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()

    """)

    def __post_init__(self):
        super().__post_init__()

        if type(self.constraints) == str:
            self.constraints = [eos.QualifiedName(self.constraints)]
        elif type(self.constraints) == eos.QualifiedName:
            self.constraints = [self.constraints]
        elif type(self.constraints) is not list:
            raise TypeError(f'constraints must be a QualifiedName or a list of QualifiedNames, not {type(self.constraints)}')

        self._parameters = eos.Parameters.Defaults()
        if self.parameters is not None and type(self.parameters) is dict:
            for key, value in self.parameters.items():
                self._parameters.set(key, value)

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the constraint residues for drawing.

        Looks up and deserializes each named constraint, evaluates the corresponding ``observable``
        at the given ``parameters``, and computes the residues (measured value minus prediction)
        together with their uncertainties as a function of the kinematic ``variable``, optionally
        rescaling binned constraints by the inverse of the bin width.

        :param context: The analysis file context. Accepted for interface consistency; this item
            does not read any data files. If ``None``, a default context is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context

        import yaml

        constraints = eos.Constraints()
        obs_kinematics = eos.Kinematics()

        self._xvalues = []
        self._xerrors = []
        self._yvalues = []
        self._yerrors = []

        for constraint_name in self.constraints:
            entry = constraints[constraint_name]
            if not entry:
                raise ValueError(f'unknown constraint {constraint_name}')

            constraint = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)

            xvalues = None
            xerrors = None
            yvalues = None
            yerrors = None

            if constraint['type'] == 'Gaussian':
                kinematics = constraint['kinematics']
                options = constraint['options']
                width = 1
                if self.variable in kinematics:
                    xvalues = [kinematics[self.variable]]
                    xerrors = [None]
                    obs_kinematics.declare(self.variable, kinematics[self.variable])

                elif (self.variable + '_min' in kinematics) and (self.variable + '_max' in kinematics):
                    var_max = float(kinematics[self.variable + '_max'])
                    var_min = float(kinematics[self.variable + '_min'])
                    xvalues = [(var_max + var_min) / 2]
                    xerrors = [(var_max - var_min) / 2]
                    if self.rescale_by_width:
                        width = (var_max - var_min)
                    obs_kinematics.declare(self.variable + '_max', var_max)
                    obs_kinematics.declare(self.variable + '_min', var_min)

                observable_value = eos.Observable.make(self.observable, self._parameters, obs_kinematics, eos.Options(options)).evaluate()
                if _np.isnan(observable_value):
                    eos.warn(f'    observable {self.observable} evaluated to NaN')
                yvalues = [(float(constraint['mean']) - observable_value) / width]
                sigma_hi = _np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2) / width
                sigma_lo = _np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2) / width
                yerrors = [(sigma_hi, sigma_lo)]
            elif constraint['type'] == 'MultivariateGaussian(Covariance)':
                if not self.observable:
                    raise KeyError('observable needs to be specified for MultivariateGaussian(Covariance) constraints')

                covariance = _np.array(constraint['covariance'])
                observables = constraint['observables']
                means = constraint['means']
                options = constraint['options']
                dim = len(means)
                kinematics = constraint['kinematics']

                xvalues = []
                xerrors = []
                yvalues = []
                yerrors = []
                for i in range(0, dim):
                    width = 1

                    # Check that the observable match and that the provided options match with those of the constraint
                    if not (observables[i] == eos.QualifiedName(self.observable).full().split(';')[0] and
                               _np.all([eos.Options(options[i])[k] == v for k, v in eos.QualifiedName(self.observable).options_part()])):
                        eos.warn('    skipping the observable {observables[i]}, with options {options[i]}, because of name or option mismatch')
                        continue
                    _kinematics = kinematics[i]
                    if self.variable in _kinematics:
                        xvalues.append(_kinematics[self.variable])
                        xerrors.append(None)
                        obs_kinematics.declare(self.variable, _kinematics[self.variable])
                    elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                        var_max = float(_kinematics[self.variable + '_max'])
                        var_min = float(_kinematics[self.variable + '_min'])
                        xvalues.append((var_max + var_min) / 2)
                        xerrors.append((var_max - var_min) / 2)
                        if self.rescale_by_width:
                            width = (var_max - var_min)
                        obs_kinematics.declare(self.variable + '_max', var_max)
                        obs_kinematics.declare(self.variable + '_min', var_min)

                    observable_value = eos.Observable.make(self.observable, self._parameters, obs_kinematics, eos.Options(options[i])).evaluate()
                    if _np.isnan(observable_value):
                        eos.warn(f'    observable {self.observable} evaluated to NaN')
                    yvalues.append((_np.double(means[i]) - observable_value) / width)
                    yerrors.append(_np.sqrt(_np.double(covariance[i, i])) / width)
            elif constraint['type'] == 'MultivariateGaussian':
                if not self.observable:
                    raise KeyError('observable needs to be specified for MultivariateGaussian constraints')
                sigma_stat_hi = _np.array(constraint['sigma-stat-hi'])
                sigma_stat_lo = _np.array(constraint['sigma-stat-lo'])
                sigma_sys = _np.array(constraint['sigma-sys'])
                sigma = _np.sqrt(_np.power(sigma_sys, 2) + 0.25 * _np.power(sigma_stat_hi + sigma_stat_lo, 2))
                observables = constraint['observables']
                means = constraint['means']
                options = constraint['options']
                dim = len(means)
                kinematics = constraint['kinematics']

                xvalues = []
                xerrors = []
                yvalues = []
                yerrors = []
                for i in range(0, dim):
                    width = 1

                    if not (observables[i] == eos.QualifiedName(self.observable).full().split(';')[0] and
                               _np.all([eos.Options(options[i])[k] == v for k, v in eos.QualifiedName(self.observable).options_part()])):
                        eos.warn('    skipping the observable {observables[i]}, with options {options[i]}, because of name or option mismatch')
                        continue
                    _kinematics = kinematics[i]
                    if self.variable in _kinematics:
                        xvalues.append(_kinematics[self.variable])
                        xerrors.append(None)
                        obs_kinematics.declare(self.variable, _kinematics[self.variable])
                    elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                        var_max = float(_kinematics[self.variable + '_max'])
                        var_min = float(_kinematics[self.variable + '_min'])
                        xvalues.append((var_max + var_min) / 2)
                        xerrors.append((var_max - var_min) / 2)
                        if self.rescale_by_width:
                            width = (var_max - var_min)
                        obs_kinematics.declare(self.variable + '_max', var_max)
                        obs_kinematics.declare(self.variable + '_min', var_min)

                    observable_value = eos.Observable.make(self.observable, self._parameters, obs_kinematics, eos.Options(options[i])).evaluate()
                    if _np.isnan(observable_value):
                        eos.warn(f'    observable {self.observable} evaluated to NaN')
                    yvalues.append((_np.double(means[i]) - observable_value) / width)
                    yerrors.append(_np.double(sigma[i]) / width)
            else:
                raise ValueError(f'constraint type {constraint["type"]} presently not supported')

            if len(xvalues) == 0:
                eos.info(f'   skipping plot for constraint {constraint_name} since it does not contain the requested observable')
                return

        self._xvalues = _np.array(xvalues)
        self._xerrors = _np.array(xerrors)
        self._yvalues = _np.array(yvalues)
        self._yerrors = _np.array(yerrors)

        if self.range:
            self._mask = _np.logical_and(self._xvalues > min(self.range), self._xvalues < max(self.range))
        else:
            self._mask = _np.array([True] * len(self._xvalues))

    def draw(self, ax):
        """Draw the constraint residues on the provided axes.

        Renders each residue as an error bar, using horizontal error bars to indicate the bin
        width of binned constraints.

        :param ax: The matplotlib axes onto which the residues are drawn.
        :type ax: matplotlib.axes.Axes
        """

        if len(self._xvalues) == 0:
            return

        xvalues = self._xvalues[self._mask]
        yvalues = self._yvalues[self._mask]
        xerrors = self._xerrors[self._mask]
        yerrors = self._yerrors[self._mask]

        label = self.label
        for xv, xerr, yv, yerr in zip(xvalues, xerrors, yvalues, yerrors):
            if xerr is not None:
                ax.errorbar(x=xv, y=yv, xerr=xerr, yerr=yerr, color=self.color,
                                         elinewidth=1.0, fmt='_', linestyle='none', label=label)
            else:
                ax.errorbar(x=xv, y=yv, yerr=yerr, color=self.color,
                                         elinewidth=1.0, fmt='_', linestyle='none', label=label)

            # disable the label for subsequent plots
            label = None

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        # the residues are drawn as capped error bars (always a y error, plus an x error
        # for binned constraints), so the swatch must match. The x error is only known
        # after prepare(); fall back to a y-only swatch if the data is not yet available.
        xerrors = getattr(self, '_xerrors', None)
        has_xerr = xerrors is not None and any(xe is not None for xe in xerrors)
        return self._legend_errorbar(has_xerr=has_xerr, has_yerr=True)


@dataclass(kw_only=True)
class BandItem(Item):
    """Plots a shaded band.

    This item type is used to display a shaded band, which can be used to represent uncertainties or confidence intervals.

    :param x: The size of the band in the x direction, given as a tuple of two float values (min, max).
    :type x: tuple[float, float] | None
    :param y: The size of the band in the y direction, given as a tuple of two float values (min, max).
    :type y: tuple[float, float] | None

    Either of ``x`` or ``y`` must be specified to define the band.

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$x$', range: [0, 10] }
          yaxis: { label: '$y$', range: [0, 10] }
          items:
            - { type: 'band', x: [2, 8], color: 'black', alpha: 0.5 }
            - { type: 'band', y: [2, 8], color: 'blue',  alpha: 0.5 }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    x:tuple[float,float]|None=field(default=None)
    y:tuple[float,float]|None=field(default=None)

    def __post_init__(self):
        super().__post_init__()

        if self.x is None and self.y is None:
            raise ValueError("Either 'x' or 'y' must be specified to define the band")

        if self.x is not None:
            if len(self.x) != 2:
                raise ValueError("'x' must be a tuple of two float values (min, max)")

            if self.x[0] >= self.x[1]:
                raise ValueError(f"Invalid 'x' range: {self.x}. The first value must be less than the second value.")

        if self.y is not None:
            if len(self.y) != 2:
                raise ValueError("'y' must be a tuple of two float values (min, max)")

            if self.y[0] >= self.y[1]:
                raise ValueError(f"Invalid 'y' range: {self.y}. The first value must be less than the second value.")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the band for drawing.

        This item requires no preparation, since its extent is fully determined by the ``x`` and
        ``y`` attributes (or the axes limits) at draw time.

        :param context: The analysis file context. Accepted for interface consistency and unused.
        :type context: AnalysisFileContext | None
        """
        pass

    def draw(self, ax):
        """Draw the shaded band on the provided axes.

        Draws a rectangle spanning ``x`` horizontally and ``y`` vertically. Whichever of the two is
        left unspecified spans the full corresponding axis range, yielding a vertical or horizontal band.

        :param ax: The matplotlib axes onto which the band is drawn.
        :type ax: matplotlib.axes.Axes
        """

        if self.x is not None:
            _xmin, _xmax = tuple(self.x)
        else:
            _xmin, _xmax = ax.get_xlim()

        if self.y is not None:
            _ymin, _ymax = tuple(self.y)
        else:
            _ymin, _ymax = ax.get_ylim()

        rect = _matplotlib.pyplot.Rectangle((_xmin, _ymin), _xmax - _xmin, _ymax - _ymin,
                                            alpha=self.alpha, color=self.color,
                                            linewidth=0, fill=True)
        ax.add_patch(rect)

    def legend(self):
        """Return the item's legend entry as a list of handle/label pairs.

        Provides a filled rectangle handle matching the band's color and transparency.

        :returns: A list containing a single ``(handle, label)`` pair if a label is set, otherwise empty.
        :rtype: list[tuple[matplotlib.artist.Artist, str]]
        """
        return self._legend_patch()

@dataclass(kw_only=True)
class VerticalLineItem(Item):
    r"""Plots a vertical line at a fixed position on the x axis.

    This item type is used to mark a particular value on the x axis, e.g. a
    kinematic threshold. The line's appearance is controlled by the inherited
    ``color``, ``linestyle``, ``linewidth``, and ``alpha`` fields.

    :param x: The x-axis position of the vertical line.
    :type x: float

    Example:

    .. code-block::

        plot:
          xaxis: { label: '$\sqrt{s}$ [GeV]', range: [3.73, 4.085] }
          yaxis: { label: '$\sigma$ [nb]' }
          items:
            - { type: 'vertical', x: 3.875, color: 'gray', linestyle: 'dashed', label: r'$D\bar{D}^*$ threshold' }
    """

    x:float|None=field(default=None)

    def __post_init__(self):
        super().__post_init__()

        if self.x is None:
            raise ValueError("'x' must be specified to define a vertical line")

    def prepare(self, context:AnalysisFileContext=None):
        "Prepare the vertical line for drawing."
        pass

    def draw(self, ax):
        "Draw the vertical line on the axes."
        ax.axvline(self.x, alpha=self.alpha, color=self.color,
                   linestyle=self.linestyle, linewidth=self.linewidth)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_line(alpha=self.alpha)

@dataclass(kw_only=True)
class SignalPDFItem(Item):
    """Plots a single EOS signal PDF w/o uncertainties as a function of one kinematic variable

    :param kinematics: The set of optional kinematic variables, given as a dictionary mapping variable names to their values.
    :type kinematics: dict[str,float] | None
    :param options: The set of optional options for the signal PDF, given as a dictionary mapping option names to their values.
    :type options: dict[str,str] | None
    :param pdf: The name of the signal PDF to be plotted.
    :type pdf: eos.QualifiedName
    :param parameters: The set of optional parameters for the signal PDF, given as a dictionary mapping parameter names to their values.
    :type parameters: dict[eos.QualifiedName,float] | None
    :param parameters_from_file: The path to a parameter file containing the optional parameters of the signal PDF.
    :type parameters_from_file: str | None
    :param range: The range of the variable to be plotted, given as a tuple of two float values (min, max).
    :type range: tuple[float,float]
    :param resolution: The number of points to be used for plotting the signal PDF. Defaults to 100.
    :type resolution: int
    :param variable: The name of the kinematic variable to which the x axis will be mapped.
    :type variable: str

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', unit: r'$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
          legend: { position: 'lower left' }
          items:
            - { type: 'signal-pdf', pdf: 'B->Dlnu::dBR/dq2;l=e,q=d', label: 'Signal PDF',
                variable: 'q2', range: [0.02, 11.63], color: 'black', resolution: 100,
                kinematics: { q2_min: 0.02, q2_max: 11.63 },
                parameters: { B->Dlnu::dBR/dq2;l=e,q=d@BSZ2015:alpha^f+_0@BSZ2015:0.75 }
              }
    """

    kinematics:dict[str,float]|None=field(default=None)
    options:dict[str,str]|None=field(default=None)
    pdf:eos.QualifiedName
    parameters:dict[eos.QualifiedName,float]|None=field(default=None)
    parameters_from_file:str|None=field(default=None)
    range:tuple[float,float]
    resolution:int=field(default=100)
    variable:str

    def __post_init__(self):
        super().__post_init__()

        if self.range is None or len(self.range) != 2:
            raise ValueError(f"Invalid range '{self.range}'. It must be a tuple of two float values (min, max).")

        if self.range[0] >= self.range[1]:
            raise ValueError(f"Invalid range '{self.range}'. The first value must be less than the second value.")

        if self.resolution <= 0:
            raise ValueError(f"Invalid resolution '{self.resolution}'. It must be a positive integer.")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the signal PDF for plotting.

        Builds the parameters (optionally overridden from file), kinematics, and options, instantiates
        the signal PDF, and caches its normalization for :meth:`draw`.

        :param context: The analysis file context used to resolve the relative path to
            ``parameters_from_file``. If ``None``, a default context rooted at the current working
            directory is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context

        self._parameters = eos.Parameters.Defaults()
        # create parameters
        if type(self.parameters_from_file) is str:
            eos.warn('    overriding parameters from file')
            self._parameters.override_from_file(context.data_path(self.parameters_from_file))

        if self.parameters is not None and self.parameters_from_file is not None:
            eos.warn('    overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

        if self.parameters is not None and type(self.parameters) is dict:
            for key, value in self.parameters.items():
                self._parameters.set(key, value)

        # create kinematics
        self._kinematics = eos.Kinematics()
        if self.variable is None:
            raise ValueError('no kinematic variable specified; do not know how to map x to a variable')
        else:
            self._variable = self._kinematics.declare(self.variable, _np.nan)

        if type(self.kinematics) is dict:
            for k, v in self.kinematics.items():
                self._kinematics.declare(k, v)

        # create options
        self._options = eos.Options()
        if self.options is not None and type(self.options) is dict:
            for key, value in self.options.items():
                self._options.declare(key, value)

        # create PDF
        self._pdf = eos.SignalPDF.make(self.pdf, self._parameters, self._kinematics, self._options)
        if self._pdf is None:
            raise ValueError(f"Signal PDF '{self.pdf}' could not be created. Please check the PDF name, the parameters, the kinematics and the options.")

        self._norm = self._pdf.normalization()

    def draw(self, ax):
        """Draw the signal PDF on the provided axes.

        Evaluates the normalized signal PDF on a grid of ``resolution`` points spanning ``range`` of
        the chosen kinematic ``variable`` and plots the resulting curve.

        :param ax: The matplotlib axes onto which the signal PDF is drawn.
        :type ax: matplotlib.axes.Axes
        """

        xvalues = _np.linspace(self.range[0], self.range[1], self.resolution)
        pvalues = _np.full(xvalues.shape, -self._norm)
        for i, xvalue in enumerate(xvalues):
            self._variable.set(xvalue)
            pvalues[i] += self._pdf.evaluate()

        ax.plot(xvalues, _np.exp(pvalues), alpha=self.alpha, color=self.color, label=self.label, lw=self.linewidth, ls=self.linestyle)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_line(alpha=self.alpha)

@dataclass(kw_only=True)
class ComplexPlaneItem(Item):
    """Plots a single observable as a function on the complex plan

    :param fixed_kinematics: A dictionary of additional kinematic variables and the values to which they shall be fixed.
    :type fixed_kinematics: dict[str, float] | None
    :param fixed_parameters: A dictionary of EOS parameters and the values to which they shall be fixed.
    :type fixed_parameters: dict[eos.QualifiedName, float] | None
    :param fixed_parameters_from_file: Path to a file containing EOS parameters to be used for the observable.
    :type fixed_parameters_from_file: str | None
    :param observable: The name of the observable to be plotted.
    :type observable: :class:`eos.QualifiedName`
    :param options: A dictionary of options to be passed to the observable.
    :type options: dict[str, str] | None
    :param ranges: A tuple of two float values (min, max) representing the range of the kinematic variable to be plotted on the x- and y-axes.
    :type ranges: tuple[tuple[float, float], tuple[float, float]]
    :param resolution: The number of points to use for the plot.
    :type resolution: int
    :param variables: The names of the variables to be used as x and y coordinates.
    :type variables: tuple[str, str]
    """
    fixed_kinematics:dict[str,float]|None=None
    fixed_parameters:dict[eos.QualifiedName,float]|None=None
    fixed_parameters_from_file:str|None=None
    observable:eos.QualifiedName
    options:dict[str, str]|None=None
    variables:tuple[str,str]
    ranges:tuple[tuple[float,float],tuple[float,float]]=field(default=((-1.0, +1.0), (-1.0, +1.0)))
    resolution:int=100

    def __post_init__(self):
        eos.info(f'Handling item to plot {self.observable} in the complex plane')
        self._observable_entry = eos.Observables()[self.observable]
        valid_kinematic_variables = {kv for kv in self._observable_entry.kinematic_variables()}

        # Create kinematics
        self._kinematics = eos.Kinematics()
        if self.fixed_kinematics:
            for k, v in self.fixed_kinematics.items():
                if k not in valid_kinematic_variables:
                    raise ValueError(f"Kinematic variable '{k}' does not match any of the declared kinematic variables '{self.observable}': {valid_kinematic_variables.__repr__()}")
                self._kinematics.declare(k, v)

        # Create parameters
        self._parameters = eos.Parameters.Defaults()
        if self.fixed_parameters_from_file:
            eos.warn('Overriding parameters from file')
            self._parameters.override_from_file(self.fixed_parameters_from_file)

        if self.fixed_parameters:
            if self.fixed_parameters_from_file:
                eos.warn('Overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            for key, value in self.fixed_parameters.items():
                self._parameters.set(key, value)

        # Declare variables
        if self.variables[0] not in valid_kinematic_variables:
            raise ValueError(f"Value of the first element of 'variables' ('{self.variables[0]}') is not the name of a kinematic variable for observable '{self.observable}'")

        if self.variables[1] not in valid_kinematic_variables:
            raise ValueError(f"Value of the second element of 'variables' ('{self.variables[1]}') is not the name of a kinematic variable for observable '{self.observable}'")

        self._kvx = self._kinematics.declare(self.variables[0], _np.nan)
        self._kvy = self._kinematics.declare(self.variables[1], _np.nan)

        # Create options
        self._options = eos.Options()
        if self.options:
            for key, value in self.options.items():
                self._options.declare(key, value)

        self._observable = eos.Observable.make(self.observable, self._parameters, self._kinematics, self._options)

        self._xrange = self.ranges[0]
        if self._xrange[0] >= self._xrange[1]:
            raise ValueError(f"Invalid x range '{self._xrange}'. The first value must be less than the second value.")

        self._yrange = self.ranges[1]
        if self._yrange[0] >= self._yrange[1]:
            raise ValueError(f"Invalid y range '{self._yrange}'. The first value must be less than the second value.")

        if self.resolution <= 0:
            raise ValueError(f"Invalid resolution '{self.resolution}'. It must be a positive integer.")

        self._xvalues = _np.linspace(self._xrange[0], self._xrange[1], self.resolution)
        self._yvalues = _np.linspace(self._yrange[0], self._yrange[1], self.resolution)
        self._cvalues = _np.reshape([[x + 1.0j * y for y in self._yvalues] for x in self._xvalues], (self.resolution**2))

    def evaluate(self, c:complex):
        """Evaluate the observable at a complex-valued kinematic variable.

        Sets the two kinematic variables given by ``variables`` to the real and imaginary parts of
        ``c``, respectively, and evaluates the observable.

        :param c: The complex value at which to evaluate the observable.
        :type c: complex
        :returns: The value of the observable at ``c``.
        :rtype: float
        """

        self._kvx.set(_np.real(c))
        self._kvy.set(_np.imag(c))

        return self._observable.evaluate()

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the observable for plotting in the complex plane.

        Evaluates the observable on a regular ``resolution`` x ``resolution`` grid of complex values
        spanning the two ranges given by ``ranges`` and caches the results for :meth:`draw`.

        :param context: The analysis file context. Accepted for interface consistency and unused.
        :type context: AnalysisFileContext | None
        """
        self._ovalues = _np.reshape(list(map(self.evaluate, self._cvalues)), (self.resolution, self.resolution)).T

    def draw(self, ax):
        """Draw the observable on the provided axes as a pseudocolor plot over the complex plane.

        :param ax: The matplotlib axes onto which the observable is drawn.
        :type ax: matplotlib.axes.Axes
        """
        ax.pcolor(self._xvalues, self._yvalues, self._ovalues, cmap='viridis', rasterized=True)

    def legend(self):
        # a pseudocolor plot has no faithful single-swatch representation,
        # so it deliberately contributes no legend entry (use a colorbar instead).
        return []


@dataclass(kw_only=True)
class ErrorBarsItem(Item):
    """Plots one or more error bars at specified position(s).

    :param positions: The list of (x, y) positions where the error bars will be plotted.
    :type positions: list[tuple[float, float]]
    :param xerrors: The list of errors to be used in the plotting. Tuples of errors are interpreted as asymmetric errors (i.e., (error_minus, error_plus)).
    :type xerrors: list[float, tuple[float, float]]
    :param yerrors: The list of errors to be used in the plotting. Tuples of errors are interpreted as asymmetric errors (i.e., (error_minus, error_plus)).
    :type yerrors: list[float, tuple[float, float]]

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$x$', range: [0, 4] }
          yaxis: { label: '$y$', range: [1, 6] }
          items:
            - type: 'errorbars'
              positions: [[1, 2], [2, 3], [3, 5]]
              xerrors: [0.5, 0.5, 0.5]
              yerrors: [0.2, [0.2, 0.3], 0.5]
              color: 'black'
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    marker:str='o'
    positions:list[tuple[float,float]]
    xerrors:list[float|tuple[float, float]]|None=None
    yerrors:list[float|tuple[float, float]]|None=None


    def __post_init__(self):
        if len(self.positions) == 0:
            raise ValueError('At least one position must be specified for error bars.')

        if self.xerrors is None and self.yerrors is None:
            raise ValueError('At least one of xerrors or yerrors must be specified for error bars.')

        if self.xerrors is not None and len(self.xerrors) != len(self.positions):
            raise ValueError(f'The number of x errors ({len(self.xerrors)}) must match the number of positions ({len(self.positions)}).')

        if self.yerrors is not None and len(self.yerrors) != len(self.positions):
            raise ValueError(f'The number of y errors ({len(self.yerrors)}) must match the number of positions ({len(self.positions)}).')

        return super().__post_init__()


    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the error bars for drawing.

        Converts the ``positions`` into x- and y-coordinate arrays and normalizes the ``xerrors`` and
        ``yerrors`` into the asymmetric (lower, upper) form expected by :meth:`draw`.

        :param context: The analysis file context. Accepted for interface consistency and unused.
        :type context: AnalysisFileContext | None
        """
        self._x = _np.array([pos[0] for pos in self.positions])
        self._y = _np.array([pos[1] for pos in self.positions])

        if self.xerrors is None:
            self._xerr = None
        else:
            self._xerr = _np.zeros((self._x.size, 2))
            for i, xerr in enumerate(self.xerrors):
                if type(xerr) is tuple or type(xerr) is list:
                    if len(xerr) != 2:
                        raise ValueError(f'Invalid x error specification {xerr} at index {i}. Must be a float or a tuple/list of two float values (error_minus, error_plus).')
                    self._xerr[i, 0] = xerr[0]
                    self._xerr[i, 1] = xerr[1]
                elif type(xerr) is float or type(xerr) is int:
                    self._xerr[i, 0] = xerr
                    self._xerr[i, 1] = xerr
                else:
                    raise ValueError(f'Invalid x error specification {xerr} at index {i} of type {type(xerr)}. Must be a float or a tuple/list of two float values (error_minus, error_plus).')
            self._xerr = self._xerr.T

        if self.yerrors is None:
            self._yerr = None
        else:
            self._yerr = _np.zeros((self._y.size, 2))
            for i, yerr in enumerate(self.yerrors):
                if type(yerr) is tuple or type(yerr) is list:
                    if len(yerr) != 2:
                        raise ValueError(f'Invalid y error specification {yerr} at index {i}. Must be a float or a tuple/list of two float values (error_minus, error_plus).')
                    self._yerr[i, 0] = yerr[0]
                    self._yerr[i, 1] = yerr[1]
                elif type(yerr) is float or type(yerr) is int:
                    self._yerr[i, 0] = yerr
                    self._yerr[i, 1] = yerr
                else:
                    raise ValueError(f'Invalid y error specification {yerr} at index {i} of type {type(yerr)}. Must be a float or a tuple/list of two float values (error_minus, error_plus).')
            self._yerr = self._yerr.T

    def draw(self, ax):
        """Draw the error bars on the provided axes.

        :param ax: The matplotlib axes onto which the error bars are drawn.
        :type ax: matplotlib.axes.Axes
        """
        ax.errorbar(self._x, self._y, xerr=self._xerr, yerr=self._yerr, fmt='none', color=self.color,
                    alpha=self.alpha, elinewidth=self.linewidth, linestyle=self.linestyle, label=self.label)
        ax.plot(self._x, self._y, marker=self.marker, linestyle='none', color=self.color, alpha=self.alpha)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return self._legend_errorbar(marker=self.marker, has_xerr=self.xerrors is not None, has_yerr=self.yerrors is not None, alpha=self.alpha)


@dataclass(kw_only=True)
class PointItem(Item):
    r"""Plots a single point at a fixed position.

    This item type is used to mark a single data point manually, e.g. a value taken from the
    literature. The point is drawn as an open (unfilled) marker whose edge uses the inherited
    ``color`` and ``alpha`` fields. To draw several points at once, or points with error bars, use
    an :class:`ErrorBarsItem` instead.

    :param x: The x-axis position of the point.
    :type x: float
    :param y: The y-axis position of the point.
    :type y: float
    :param marker: The marker style, given as a valid matplotlib marker. Defaults to ``'x'``.
    :type marker: str
    :param markersize: The marker size in points. Defaults to None, i.e. matplotlib's default size.
    :type markersize: float | None

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: r'$q^2$', range: [0.0, 1.0] }
          yaxis: { label: r'$f_+(q^2)$' }
          items:
            - { type: 'point', x: 0.0, y: 0.261, marker: 'o', markersize: 12, color: 'C0',
                label: r'LCSR (Bharucha 2012)'
              }
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    x:float
    y:float
    marker:str=field(default='x')
    markersize:float|None=field(default=None)

    def __post_init__(self):
        super().__post_init__()

        if self.x is None:
            raise ValueError("'x' must be specified to define a point")

        if self.y is None:
            raise ValueError("'y' must be specified to define a point")

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the point for drawing.

        :param context: The analysis file context. Accepted for interface consistency and unused.
        :type context: AnalysisFileContext | None
        """
        pass

    def draw(self, ax):
        """Draw the point on the provided axes.

        :param ax: The matplotlib axes onto which the point is drawn.
        :type ax: matplotlib.axes.Axes
        """
        ax.plot(self.x, self.y, alpha=self.alpha, color=self.color, label=self.label, linestyle='none',
                marker=self.marker, markeredgecolor=self.color, markerfacecolor='none', markersize=self.markersize)

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        # the point is drawn as an open marker, so its swatch must be open (unfilled) too
        if not self.label:
            return []
        handle = _matplotlib.lines.Line2D((0,), (0,), color=self.color, marker=self.marker, linestyle='none',
                                          markeredgecolor=self.color, markerfacecolor='none', markersize=self.markersize,
                                          alpha=self.alpha)
        return [(handle, self.label)]


class ItemFactory:
    """Factory that creates :class:`Item` instances from their YAML or dictionary description.

    The concrete item class is selected from the mandatory ``type`` key using the :attr:`registry`,
    which maps each supported type string (e.g. ``'observable'``, ``'kde1D'``, ``'constraint'``) to
    its corresponding :class:`Item` subclass.
    """

    registry = {
        'observable': ObservableItem,
        'expression': ExpressionItem,
        'uncertainty': UncertaintyBandItem,
        'uncertainty-binned': BinnedUncertaintyItem,
        'constraint': ConstraintItem,
        'constraint-residue': ConstraintResidueItem,
        'constraint2D': TwoDimensionalConstraintItem,
        'histogram1D': OneDimensionalHistogramItem,
        'histogram2D': TwoDimensionalHistogramItem,
        'kde1D': OneDimensionalKernelDensityEstimateItem,
        'kde2D': TwoDimensionalKernelDensityEstimateItem,
        'contours2D': TwoDimensionalContoursItem,
        'band': BandItem,
        'vertical': VerticalLineItem,
        'signal-pdf': SignalPDFItem,
        'complex-plane': ComplexPlaneItem,
        'errorbars': ErrorBarsItem,
        'point': PointItem,
    }

    @staticmethod
    def from_yaml(yaml_data:str):
        """Create an item from a YAML description.

        :param yaml_data: A YAML string describing a single item, including its ``type`` key.
        :type yaml_data: str
        :returns: The instantiated item.
        :rtype: Item
        """
        kwargs = _yaml.safe_load(yaml_data)
        return ItemFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        """Create an item from its keyword description.

        The mandatory ``type`` key selects the concrete :class:`Item` subclass from :attr:`registry`;
        the remaining keyword arguments are forwarded to that subclass.

        :param kwargs: The item description. Must contain a ``type`` key identifying a registered item type.
        :returns: The instantiated item.
        :rtype: Item
        :raises ValueError: If the ``type`` key is missing or names an unknown item type.
        """
        if 'type' not in kwargs:
            raise ValueError(f"Content item { kwargs } does not contain element 'type'")

        if kwargs['type'] not in ItemFactory.registry:
            raise ValueError(f'Unknown content item type: {kwargs["type"]}')

        item_type = kwargs.pop('type')

        return ItemFactory.registry[item_type].from_dict(**kwargs)
