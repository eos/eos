# Copyright (c) 2023-2025 Danny van Dyk
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

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from eos.analysis_file_description import AnalysisFileContext
from eos.deserializable import Deserializable

import copy as _copy
import inspect
import eos
import matplotlib as _matplotlib
import numpy as _np
import os
import scipy as _scipy
import yaml as _yaml

class ItemColorCycler:
    _colors = [
        # 10-colors colorblind-friendly cycle from 2107.02270
        (0.247, 0.565, 0.855), (1.000, 0.663, 0.055), (0.741, 0.122, 0.004), (0.580, 0.643, 0.635), (0.514, 0.176, 0.714), (0.663, 0.420, 0.349), (0.906, 0.388, 0.000), (0.725, 0.675, 0.439), (0.443, 0.459, 0.506), (0.573, 0.855, 0.867)
    ]
    _color_idx = 0

    @classmethod
    def reset(cls):
        cls._color_idx = 0

    @classmethod
    def next_color(cls):
        """Returns the next available color"""
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
        "Prepare the item for drawing"
        raise NotImplementedError

    @abstractmethod
    def draw(self, ax):
        "Draw the item on the axes"
        raise NotImplementedError

    def legend(self):
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        return ()


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
        "Prepare the drawing by evaluating the observable at the sample points."
        context = AnalysisFileContext() if context is None else context
        self._yvalues = _np.empty((len(self._xvalues),))
        for i, x in enumerate(self._xvalues):
            self._variable.set(x)
            self._yvalues[i] = self._observable.evaluate()


    def draw(self, ax, **kwargs):
        "Draw a curve of the observable."
        ax.plot(self._xvalues, self._yvalues, label=self.label, color=self.color, **kwargs)


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
    """

    band:set[str]|list[str]|str=field(default_factory=lambda : {'area', 'outer', 'median'})
    datafile:str
    interpolation:str=field(default='linear')
    range:tuple[float, float]|None=field(default=None)
    resolution:int=field(default=100)
    variable:str|None=field(default=None)

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
        context = AnalysisFileContext() if context is None else context

        self._datafile = eos.data.Prediction(context.data_path(self.datafile))
        if len(self._datafile.varied_parameters) < 1:
            raise ValueError(f"Data file '{self.datafile}' does not contain any predictions")

        self._variable = self.variable
        if self._variable is not None:
            # assume that the first prediction has only one kinematic variable,
            # which we adopt as the variable to plot
            kinematics = self._datafile.varied_parameters[0]['kinematics']
            if len(kinematics) > 1:
                raise ValueError(f"Data file '{self.datafile}' contains more than one kinematic variable; specify 'variable' to determine which is supposed to be used")
            self._variable = list(kinematics.keys())[0]

        _xvalues = []
        for idx, vp in enumerate(self._datafile.varied_parameters):
            kinematics = vp['kinematics']
            if self._variable not in kinematics:
                raise ValueError(f"Prediction #{idx} in data file '{self.datafile}' does not depend on the chosen kinematic variable '{self._variable}'")
            _xvalues.append(kinematics[self._variable])

        _xvalues = _np.array(_xvalues)
        _samples = self._datafile.samples
        _weights = self._datafile.weights

        _ovalues_lower   = []
        _ovalues_central = []
        _ovalues_higher  = []
        INTERVAL = [0.15865, 0.5, 0.84135]  # central 68% interval
        for i in range(len(_samples[0])):
            lower, central, higher = eos.plot.Plotter._weighted_quantiles(_samples[:, i], INTERVAL, _weights)
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
        """Prepare the item for drawing."""
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
        """Draw the uncertainty band."""

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
            ax.plot([xmin, xmax], [olo,      olo],      color=self.color, alpha=self.alpha)
            ax.plot([xmin, xmax], [ocentral, ocentral], color=self.color, alpha=self.alpha)
            ax.plot([xmin, xmax], [ohi,      ohi],      color=self.color, alpha=self.alpha)
            label = None


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
        "Prepare the histogram for drawing"

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
        "Draw the histogram"
        ax.hist(self.samples, weights=self._datafile.weights, range=self.range,
                alpha=self.alpha, bins=self.bins, color=self.color, density=True, label=self.label)

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
        "Prepare the histogram for drawing."

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
        "Draw the two-dimensional histogram."
        ax.hist2d(self.samples[:, 0], self.samples[:, 1], weights=self.weights, range=[self.xrange, self.yrange],
                  bins=self.bins, label=self.label, rasterized=True, cmap='Greys')


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
        "Prepare the KDE for drawing."

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 1D KDE")
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

        self.idx = self._datafile.lookup_table[self.variable]

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
        "Draw the KDE."
        # find the PDF value corresponding to a given cumulative probability
        if self.level is not None:
            plevelf = lambda x, pdf, P: self.pdf[self.pdf > x].sum() - P
            plevel = _scipy.optimize.brentq(plevelf, 0., 1., args=(self.pdf, self.level / 100.0 ))
            ax.fill_between(_np.ma.masked_array(self.xvalues, mask=self.pdf < plevel),
                                            _np.ma.masked_array(self.pdf, mask=self.pdf < plevel, fill_value=_np.nan),
                                            facecolor=self.color, alpha=self.alpha)

        ax.plot(self.xvalues, self.pdf, color=self.color, linestyle=self.linestyle, label=self.label)

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
        "Prepare the KDE for drawing."

        context = AnalysisFileContext() if context is None else context

        # These checks are necessary to ensure that the data file is in the correct format,
        # but they are not possible in the __post_init__ method, because the data file might not yet exist.
        datafile = context.data_path(self.datafile)
        os.path.exists(datafile) or eos.error(f"Data file '{datafile}' does not exist when preparing 2D KDE")
        name = os.path.split(datafile)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(datafile)
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(datafile)
        else:
            eos.error(f"Data file '{datafile}' has an unsupported format")
            raise NotImplementedError

        if self.variables[0] not in self._datafile.lookup_table:
            raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variables[0]}'")
        if self.variables[1] not in self._datafile.lookup_table:
            raise ValueError(f"Data file '{datafile}' does not contain samples of variable '{self.variables[1]}'")

        self._xidx = self._datafile.lookup_table[self.variables[0]]
        self._yidx = self._datafile.lookup_table[self.variables[1]]

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

    def draw(self, ax):
        "Draw the KDE."
        # find the PDF value corresponding to a given cummulative probability
        plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
        plevels = []
        labels = []
        for level in self.levels:
            plevels.append(_scipy.optimize.brentq(plevel, 0., 1., args=(self._pdf, level / 100.0)))
            labels.append(f'{level}%')

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
        """Return the item's legend entry in form of its handle(s) and label(s)."""
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
          yaxis: { label: r'$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
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
        """Prepare the constraint for drawing."""
        context = AnalysisFileContext() if context is None else context

        import yaml

        constraints = eos.Constraints()
        obs_kinematics = eos.Kinematics()

        self._xvalues = []
        self._xerrors = []
        self._yvalues = []
        self._yerrors = []
        self._skip_observable = False

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
                    xvalues = [(kinematics[self.variable + '_max'] + kinematics[self.variable + '_min']) / 2]
                    xerrors = [(kinematics[self.variable + '_max'] - kinematics[self.variable + '_min']) / 2]
                    if self.rescale_by_width:
                        width = (kinematics[self.variable + '_max'] - kinematics[self.variable + '_min'])

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
                        self._skip_observable = True
                        continue
                    _kinematics = kinematics[i]
                    if self.variable in _kinematics:
                        xvalues.append(_kinematics[self.variable])
                        xerrors.append(None)
                        if self.plot_residues:
                            obs_kinematics.declare(self.variable, _kinematics[self.variable])
                            obs_values.append(eos.Observable.make(self.observable, self.parameters, obs_kinematics, self.options).evaluate())

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

                    if not observables[i] == self.observable:
                        continue
                    _kinematics = kinematics[i]
                    if self.variable in _kinematics:
                        xvalues.append(_kinematics[self.variable])
                        xerrors.append(None)
                        if self.plot_residues:
                            obs_values = []
                            for x in xvalues:
                                obs_kinematics.declare(self.variable, x)
                                obs_values.append(eos.Observable.make(self.observable, self.parameters, obs_kinematics, self.options).evaluate())

                    elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                        xvalues.append((_kinematics[self.variable + '_max'] + _kinematics[self.variable + '_min']) / 2)
                        xerrors.append((_kinematics[self.variable + '_max'] - _kinematics[self.variable + '_min']) / 2)
                        if self.rescale_by_width:
                            width = (_kinematics[self.variable + '_max'] - _kinematics[self.variable + '_min'])
                        if self.plot_residues:
                            obs_values = []
                            for xmax, xmin in zip(kinematics[self.variable + '_max'], kinematics[self.variable + '_min']):
                                obs_kinematics.declare(self.variable + '_max', xmax)
                                obs_kinematics.declare(self.variable + '_min', xmin)
                                obs_values.append(eos.Observable.make(self.observable, self.parameters, obs_kinematics, self.options).evaluate() / width)

                    yvalues.append(means[i] / width)
                    yerrors.append(sigma[i] / width)
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
        "Draw the constraint on the axes."

        if self._skip_observable:
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
        """Prepare the band for drawing."""
        pass

    def draw(self, ax):
        "Draw the band on the axes."

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
        """Return the item's legend entry in form of its handle(s) and label(s)."""
        entries = []

        if self.label:
            handle = _matplotlib.pyplot.Rectangle((0,0),1,1, alpha=self.alpha, color=self.color)
            entries.append((handle, self.label))

        return entries

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
        """Prepare the signal PDF for plotting."""
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
        """Draw the signal PDF on the axes."""

        xvalues = _np.linspace(self.range[0], self.range[1], self.resolution)
        pvalues = _np.full(xvalues.shape, -self._norm)
        for i, xvalue in enumerate(xvalues):
            self._variable.set(xvalue)
            pvalues[i] += self._pdf.evaluate()

        ax.plot(xvalues, _np.exp(pvalues), alpha=self.alpha, color=self.color, label=self.label, ls=self.linestyle)

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
        """Evaluate the observable at a complex-valued kinematic variable."""

        self._kvx.set(_np.real(c))
        self._kvy.set(_np.imag(c))

        return self._observable.evaluate()

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the observable for plotting in the complex plane."""
        self._ovalues = _np.reshape(list(map(self.evaluate, self._cvalues)), (self.resolution, self.resolution)).T

    def draw(self, ax):
        """Draw the observable on the axes."""
        ax.pcolor(self._xvalues, self._yvalues, self._ovalues, cmap='viridis', rasterized=True)


@dataclass(kw_only=True)
class ErrorBars(Item):
    """Plots one or more error bars at specified position(s).

    :param positions: The list of (x, y) positions where the error bars will be plotted.
    :type positions: list[tuple[float, float]]
    :param xerrors: The list of errors to be used in the plotting. Tuples of errors are interpreted as asymmetric errors (i.e., (error_minus, error_plus)).
    :type xerrors: list[float] | list[tuple[float, float]]
    :param yerrors: The list of errors to be used in the plotting. Tuples of errors are interpreted as asymmetric errors (i.e., (error_minus, error_plus)).
    :type yerrors: list[float] | list[tuple[float, float]]

    Example:

    .. code-block::

        figure_args = '''
        plot:
          xaxis: { label: '$x$', range: [0, 4] }
          yaxis: { label: '$y$', range: [1, 6] }
          items:
            - type: 'errorbar'
              positions: [(1, 2), (2, 3), (3, 5)]
              xerrors: [0.5, 0.5, 0.5]
              yerrors: [0.2, (0.2, 0.3), 0.5]
              color: 'black'
        '''
        figure = eos.figure.FigureFactory.from_yaml(figure_args)
        figure.draw()
    """

    marker:str='o'
    positions:list[tuple[float,float]]
    xerrors:list[float]|list[tuple[float,float]]|None=None
    yerrors:list[float]|list[tuple[float,float]]|None=None


    def __post_init__(self):
        if len(self.positions) == 0:
            raise ValueError('At least one position must be specified for error bars.')

        if self.xerrors is None and self.yerrors is None:
            raise ValueError('At least one of xerrors or yerrors must be specified for error bars.')

        if self.xerrors is not None and len(self.xerrors) != len(self.positions):
            raise ValueError('The number of x errors must match the number of positions.')

        if self.yerrors is not None and len(self.yerrors) != len(self.positions):
            raise ValueError('The number of y errors must match the number of positions.')

        return super().__post_init__()


    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the error bars for drawing."""
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
        "Draw the error bars on the axes."
        ax.errorbar(self._x, self._y, xerr=self._xerr, yerr=self._yerr, fmt='none', color=self.color,
                    alpha=self.alpha, elinewidth=self.linewidth, linestyle=self.linestyle, label=self.label)
        ax.plot(self._x, self._y, marker=self.marker, linestyle='none', color=self.color, alpha=self.alpha)

class ItemFactory:
    registry = {
        'observable': ObservableItem,
        'uncertainty': UncertaintyBandItem,
        'uncertainty-binned': BinnedUncertaintyItem,
        'constraint': ConstraintItem,
        'histogram1D': OneDimensionalHistogramItem,
        'histogram2D': TwoDimensionalHistogramItem,
        'kde1D': OneDimensionalKernelDensityEstimateItem,
        'kde2D': TwoDimensionalKernelDensityEstimateItem,
        'band': BandItem,
        'signal-pdf': SignalPDFItem,
        'complex-plane': ComplexPlaneItem,
        'errorbars': ErrorBars,
    }

    @staticmethod
    def from_yaml(yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        return ItemFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        "Factory method to create an item from a dictionary"
        if 'type' not in kwargs:
            raise ValueError(f"Content item { kwargs } does not contain element 'type'")

        if kwargs['type'] not in ItemFactory.registry:
            raise ValueError(f'Unknown content item type: {kwargs["type"]}')

        item_type = kwargs.pop('type')

        return ItemFactory.registry[item_type].from_dict(**kwargs)
