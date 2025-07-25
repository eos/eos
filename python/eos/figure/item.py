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


@dataclass(kw_only=True)
class ObservableItem(Item):
    """Plots an observables from the EOS library of builtin observables"""

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
          xaxis: { label: '$q^2$',                range: [0.0, 11.6]   }
          yaxis: { label: '$d\\mathcal{B}/dq^2$', range: [0.0, 5.0e-3] }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'mu' },  label: '$\\ell = \\mu$',
                variable: 'q2', range: { min: 0.02, max: 11.6, num: 5800 }
              }
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'tau' }, label: '$\\ell = \\tau$',
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
                # Continues only if no failure occures
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
        "Evaluate the observable at the sample points"
        context = AnalysisFileContext() if context is None else context
        self._yvalues = _np.empty((len(self._xvalues),))
        for i, x in enumerate(self._xvalues):
            self._variable.set(x)
            self._yvalues[i] = self._observable.evaluate()


    def draw(self, ax, **kwargs):
        "Draw a line plot of the observable"
        ax.plot(self._xvalues, self._yvalues, label=self.label, color=self.color, **kwargs)


@dataclass(kw_only=True)
class UncertaintyBandItem(Item):
    """Plots an uncertainty band as a function of one kinematic variable or one parameter.

    This routine expects the uncertainty propagation to have produced an EOS data file"""

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
          xaxis: { label: '$q^2$', unit: '$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: '$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
          legend: { position: 'upper center' }
          items:
            - { type: 'uncertainty', label: '$\\ell=\\mu$',
                variable: 'q2', range: [0.02, 11.63],
                datafile: './predictions-data/FF-LQCD-SSE/pred-B-to-D-mu-nu'
              }
            - { type: 'uncertainty', label: '$\\ell=\\tau$',
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
class OneDimensionalHistogramItem(Item):
    """Plots a one-dimensional histogram."""

    bins:int=field(default=100)
    datafile:str
    variable:str

    _api_doc = inspect.cleandoc("""\
    Plotting One-Dimensional Histograms
    -----------------------------------

    Plot items of type ``histogram1D`` are used to display samples of a probability density, be it a prior, a posterior, or a signal PDF.
    The following key is mandatory:

        * ``datafile`` (*str*, path to an existing data file of type *eos.data.ImportanceSamples* or *eos.data.Prediction*) -- The path to
        a data file that was generated with one of the ``sample-nested`` or ``predict-observables`` tasks.

    The following keys are optional:

        * ``bins`` (*int*) -- The number of histogram bins. Defaults to 100.
        * ``variable`` (*str*) -- The name of the variable that is plotted on the x-axis. Defaults to the first variable in the data file.

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

        pass

    def draw(self, ax):
        "Draw the histogram"
        idx = self._datafile.lookup_table[self.variable]
        ax.hist(self._datafile.samples[:, idx], weights=self._datafile.weights,
                alpha=self.alpha, bins=self.bins, color=self.color, density=True, label=self.label)

@dataclass
class TwoDimensionalHistogramItem(Item):
    """Plots a two-dimensional histogram."""

#    :param datafile: Path to the file that contains the histogram data.
#    :type datafile: str
#    :param variables: Names of two of the datafile's variables that shall be histogrammed.
#    :type variables: tuple[str, str]
#    :param bins: Number of histogram bins.
#    :type bins: int
#    """

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
    """)

    datafile:str
    variables:tuple[str, str]
    bins:int=field(default=50)

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

    def draw(self, ax):
        xidx = self._datafile.lookup_table[self.variables[0]]
        yidx = self._datafile.lookup_table[self.variables[1]]
        ax.hist2d(self._datafile.samples[:, xidx], self._datafile.samples[:, yidx],
                  bins=self.bins, cmin=1, label=self.label)


@dataclass(kw_only=True)
class OneDimensionalKernelDensityEstimateItem(Item):
    """Plotting a one-dimensional kernel density estimate (KDE)."""

#    :param datafile: Path to the file that contains the KDE data.
#    :type datafile: str
#    :param variable: Name of the datafile's variable that shall be visualized.
#    :type variable: str
#    :param level: Credibility level that shall be visualized in percent (optional).
#    :type level: float
#    :param bandwidth: Relative bandwidth factor that multiplies the automatically determined bandwidth of the KDE (optional).
#    :type bandwidth: float

    bandwidth:float=field(default=None)
    datafile:str
    level:float=field(default=None)
    range:tuple[float, float]=field(default=None)
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
        * ``level`` (*float*) -- The credibility level that shall be visualized in percent (optional).
        * ``range`` (*tuple* of two *float* values) -- The range of the variable to be plotted on the x-axis. Defaults to the full range of the variable in the data file.
        * ``xsamples`` (*int*) -- The number of samples to be used for the x-axis. Defaults to 100.

    """)

    def __post_init__(self):
        super().__post_init__()
        if self.level is not None and (self.level <= 0 or self.level >= 100):
            raise ValueError(f"Credibility level '{self.level}' is not in the interval (0, 100)")

        if self.bandwidth is not None and self.bandwidth <= 0.0:
            raise ValueError(f"Bandwidth factor '{self.bandwidth}' is not positive")

    def prepare(self, context:AnalysisFileContext=None):
        "Prepare the KDE for drawing"

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

        if self.range is None:
            self.range = (self._datafile.samples[:, self.idx].min(), self._datafile.samples[:, self.idx].max())

        self.kde = _scipy.stats.gaussian_kde(self._datafile.samples[:, self.idx], weights=self._datafile.weights)
        self.kde.set_bandwidth(bw_method='silverman')
        if self.bandwidth is not None:
            self.kde.set_bandwidth(bw_method=self.kde.factor * self.bandwidth)

        self.xvalues = _np.linspace(self.range[0], self.range[1], self.xsamples)
        self.pdf = self.kde(self.xvalues)
        self.pdf_norm = self.pdf.sum()

    def draw(self, ax):
        "Draw the histogram"
        # find the PDF value corresponding to a given cummulative probability
        if self.level is not None:
            plevelf = lambda x, pdf, P: self.pdf[self.pdf > x * self.pdf_norm].sum() - P * self.pdf_norm
            plevel = _scipy.optimize.brentq(plevelf, 0., 1., args=(self.pdf, self.level / 100.0 ))
            ax.fill_between(_np.ma.masked_array(self.xvalues, mask=self.pdf < plevel * self.pdf_norm),
                                            _np.ma.masked_array(self.pdf, mask=self.pdf < plevel * self.pdf_norm, fill_value=_np.nan),
                                            facecolor=self.color, alpha=self.alpha)

        ax.plot(self.xvalues, self.pdf, color=self.color, linestyle=self.linestyle, label=self.label)

@dataclass(kw_only=True)
class TwoDimensionalKernelDensityEstimateItem(Item):
    """Plots a two-dimensional kernel density estimate (KDE)."""

#    :param datafile: Path to the file that contains the KDE data.
#    :type datafile: str
#    :param variables: Name of the datafile's variable that shall be visualized.
#    :type variables: str
#    :param level: Credibility level that shall be visualized in percent (optional).
#    :type level: float
#    :param bandwidth: Relative bandwidth factor that multiplies the automatically determined bandwidth of the KDE (optional).
#    :type bandwidth: float

    bandwidth:float=field(default=None)
    contours:set[str]=field(default_factory=lambda : {'lines'})
    datafile:str
    levels:list[float]=field(default=None)
    variables:tuple[str,str]
    xrange:tuple[float,float]=field(default=None)
    yrange:tuple[float,float]=field(default=None)

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
        "Prepare the KDE for drawing"

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

        self._kde = _scipy.stats.gaussian_kde(samples.T, weights=weights)
        self._kde.set_bandwidth(bw_method='silverman')
        if self.bandwidth is not None:
            self._kde.set_bandwidth(bw_method=self._kde.factor * self.bandwidth)

    def draw(self, ax):
        "Draw the histogram"

        # determine the extent of the plot
        xrange = ax.get_xlim() if self.xrange is None else self.xrange
        yrange = ax.get_ylim() if self.yrange is None else self.yrange

        # compute the PDF on a grid
        xx,yy = _np.mgrid[xrange[0]:xrange[1]:100j, yrange[0]:yrange[1]:100j]
        self._positions = _np.vstack([xx.ravel(), yy.ravel()])
        self._pdf = _np.reshape(self._kde(self._positions).T, xx.shape)
        self._pdf /= self._pdf.sum()

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
                        extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                        levels=plevels[::-1])

        CS = ax.contour(self._pdf.transpose(),
                        colors=self.color,
                        extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
                        levels=plevels[::-1],
                        linestyles=self.linestyle)

        if 'labels' in self.contours:
            fmt = {}
            for level, label in zip(CS.levels, labels[::-1]):
                fmt[level] = label

            ax.clabel(CS, inline=1, fmt=fmt, fontsize=10)

@dataclass(kw_only=True)
class ConstraintItem(Item):
    """Plots statistical constraints from the EOS library of experimental and theoretical likelihoods"""

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
          xaxis: { label: '$q^2$', unit: '$\\textnormal{GeV}^2$', range: [0.0, 11.63] }
          yaxis: { label: '$d\\mathcal{B}/dq^2$',                 range: [0.0,  5e-3] }
          legend: { position: 'lower left' }
          items:
            - { type: 'observable', observable: 'B->Dlnu::dBR/dq2;l=e,q=d', label: '$\\ell=e$',
                variable: 'q2', range: { min: 0.02, max: 11.63, num: 100}, color: 'black'
              }
            - { type: 'constraint', 'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',  observable: 'B->Dlnu::BR', label: 'Belle 2015 $\\ell=e,\\, q=d$',
                variable: 'q2', rescale_by_width: true
              }
            - { type: 'constraint', 'constraints': 'B^0->D^+mu^-nu::BRs@Belle:2015A', observable: 'B->Dlnu::BR', label: 'Belle 2015 $\\ell=\\mu,\\, q=d$',
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
        """Prepare the constraint for drawing"""
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
                sigma = _np.sqrt(np.power(sigma_sys, 2) + 0.25 * np.power(sigma_stat_hi + sigma_stat_lo, 2))
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
                eos.info(f'   skipping plot for constraint {name} since it does not contain the requested observable')
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
        "Draw the constraint on the axes"

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

class ItemFactory:
    registry = {
        'observable': ObservableItem,
        'uncertainty': UncertaintyBandItem,
        'constraint': ConstraintItem,
        'histogram1D': OneDimensionalHistogramItem,
        'histogram2D': TwoDimensionalHistogramItem,
        'kde1D': OneDimensionalKernelDensityEstimateItem,
        'kde2D': TwoDimensionalKernelDensityEstimateItem,
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
