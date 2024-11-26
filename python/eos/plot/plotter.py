# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2017, 2018, 2021 Danny van Dyk
# Copyright (c) 2021 Philip Lüghausen, Méril Reboud
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

import eos
import inspect
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import gaussian_kde
import sys
import os

class Plotter:
    """Produces publication-quality plots

    Plots can contain EOS observables, EOS constraints, and :class:`Analysis <eos.Analysis>` results.
    See `Plot description format`_ for documentation of how to create a plot.

    :param description: Description of the plot and its contents, see `Plot description format`_.
    :type description: dict
    :param output: Name of the output file. The file format is automatically determined based on the file's extension.
    :type output: string, optional
    """
    def __init__(self, description, output=None):
        self.instructions = description
        self.output = output
        self.fig = None
        self.ax = None
        self.xrange = None
        self.yrange = None
        self.__next_z_order = 0
        self.__next_color = 0
        self.colors = [
            'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'
        ]

    @property
    def next_z_order(self):
        """Returns the next available z-order value, incremented for each plot w/o pre-defined z-order value"""
        result = self.__next_z_order
        self.__next_z_order += 1
        return(result)

    @property
    def next_color(self):
        """Returns the next available color"""
        result = self.colors[self.__next_color]
        self.__next_color = (self.__next_color + 1) % len(self.colors)
        return(result)

    @staticmethod
    def _weighted_quantiles(values, quantiles, weights=None):
        """ Compute the quantiles of a weighted sample, ignoring NaN values.
        :param values: (array-like of *float*) -- Sample from which the quantiles are computed.
        :param quantiles: (array-like of *float*) -- Quantiles to compute, the values must be in [0, 1]
        :param weights: (array-like of *float*) -- Array of weights, which should be of the same length as values
        :return: numpy.array with computed quantiles.
        """
        if weights is None:
            weights = np.ones(len(values))
        valid_indices = ~(np.isnan(values) + np.isnan(weights))
        values = np.array(values)[valid_indices]
        weights = np.array(weights)[valid_indices]
        quantiles = np.array(quantiles)
        if (np.any(quantiles < 0) or np.any(quantiles > 1)):
            eos.error('Quantiles should be in [0, 1]')
        if (np.any(weights < 0)):
            eos.error('The sample weights cannot be negative')
        if (np.sum(weights) == 0):
            eos.error('The sum of the sample weights evaluated to zero')

        sorter = np.argsort(values)
        values = values[sorter]
        weights = weights[sorter]

        # Each sample's weight is half assigned to the the left and half assigned to the right of the point
        weighted_quantiles = np.cumsum(weights) - 0.5 * weights
        weighted_quantiles /= np.sum(weights)
        return np.interp(quantiles, weighted_quantiles, values)


    def setup_plot(self):
        """Setting up the plot based on the provided instruction"""
        if not 'plot' in self.instructions:
            raise KeyError('no plot metadata specified')

        if 'axis' in self.instructions:
            self.ax = self.instructions['axis']
        else:
            self.fig, self.ax = plt.subplots()

        myplot = self.instructions['plot']

        mytitle = ''
        myylabel = ''
        myxlabel = ''
        myyscale = 'linear'
        myxscale = 'linear'

        if 'title' in myplot:
            mytitle = myplot['title']

        if 'size' in myplot:
            xwidth, ywidth = myplot['size']
            # size is specified in cm, matplotlib expects inches
            # convert from cm to inches
            xwidth /= 2.54 # cm / inch
            ywidth /= 2.54 # cm / inch
            plt.gcf().set_size_inches((xwidth, ywidth))

        if 'axes' in myplot:
            self.ax.axis(myplot['axes'])

        if 'x' in myplot:
            myx = myplot['x']

            if 'label' in myx:
                myxlabel = myx['label']

            if 'unit' in myx:
                myxlabel += r'\,[' + myx['unit'] + r']'

            if 'range' in myx:
                self.xrange = myx['range']
                self.ax.set_xlim(tuple(self.xrange))

            self.ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            self.ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            self.ax.xaxis.set_ticks_position('both')

            if 'scale' in myx:
                myxscale = myx['scale']
                self.ax.set_xscale(myxscale)

            if 'scaling_factor' in myx:

                if 'format' in myx:
                    self.xformat = myx['format']
                else:
                    self.xformat = '${x}$'
                    eos.warn("Argument plot:x:format might be required when using plot:x:scale to avoid side effects")

                self.x_scaling_factor = float(myx['scaling_factor'])
                self.xticks = matplotlib.ticker.FuncFormatter(lambda x, pos, xscale=self.x_scaling_factor: self.xformat.format(x=x / xscale))
                self.ax.xaxis.set_major_formatter(self.xticks)

            else:
                if 'format' in myx:
                    eos.warn("Argument plot:x:format is only used when plot:x:scale is used")

        if 'y' in myplot:
            myy = myplot['y']

            if 'label' in myy:
                myylabel = myy['label']

            if 'unit' in myy:
                myylabel += r'\,[' + myy['unit'] + r']'

            if 'range' in myy:
                self.yrange = myy['range']
                self.ax.set_ylim(tuple(self.yrange))

            self.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            self.ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            self.ax.yaxis.set_ticks_position('both')

            if 'scale' in myy:
                myyscale = myy['scale']
                self.ax.set_yscale(myyscale)

            if 'scaling_factor' in myy:

                if 'format' in myy:
                    self.yformat = myy['format']
                else:
                    self.yformat = '${x}$'
                    eos.warn("Argument plot:y:format might be required when using plot:y:scale to avoid side effects")

                self.y_scaling_factor = float(myy['scaling_factor'])
                self.yticks = matplotlib.ticker.FuncFormatter(lambda y, pos, yscale=self.y_scaling_factor: self.yformat.format(x=y / yscale))
                self.ax.yaxis.set_major_formatter(self.yticks)

            else:
                if 'format' in myy:
                    eos.warn("Argument plot:y:format is only used when plot:y:scale is used")

        if 'grid' in myplot:
            self.ax.grid(b=True, which=myplot['grid'])

        self.ax.set(xlabel=myxlabel, ylabel=myylabel, title=mytitle)

    class BasePlot:
        """Base class for any of the plots supported by Plotter"""

        def __init__(self, plotter, item):
            self.plotter = plotter
            self.item    = item
            self.z_order = int(item['z-order']) if 'z-order' in item else plotter.next_z_order

            # obtain the defaults if otherwise unavailable
            self.alpha         = item['opacity'] if 'opacity' in item else 0.5
            self.color         = item['color']   if 'color'   in item else self.plotter.next_color
            self.label         = item['label']   if 'label'   in item else None
            self.xsamples      = item['samples'] if 'samples' in item else 100
            self.style         = item['style']   if 'style'   in item else 'solid'
            self.lw            = item['lw']      if 'lw'      in item else None
            self.xlo, self.xhi = item['range']   if 'range'   in item else self.plotter.xrange

        def __lt__(self, other):
            return(self.z_order < other.z_order)

        def handles_labels(self):
            return ([], [])


    class Point(BasePlot):
        """Plots a single point"""

        _api_doc = inspect.cleandoc("""\
        Plotting Points
        ---------------

        Content items of type ``point`` are used to display a single data point manually.
        The following keys are mandatory:

         * ``x`` (*number*) -- The point's x coordinate.
         * ``y`` (*number*) -- The point's y coordinate.

        Beside the common set of optional keys, this item type recognizes the following optional
        keys:

         * ``marker`` (*str*, a valid Matplotlib marker style) -- The point's marker style.
         * ``markersize`` (*number*, a valid Matplotlib marker size) -- The point's marker size in pts.

        Example:

        .. code-block::

           plot_args = {
               'plot': { ... },
               'contents': [
                   {
                       'label': r'LCSR (Bharucha 2012)',
                       'type': 'point',
                       'color': 'C0',
                       'x': 0,
                       'y': 0.261,
                       'marker': 'o',
                       'markersize': 12
                   }
               ]
           }
        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'x' not in item:
                raise KeyError('x coordinate not provided')
            self.x = item['x']

            if 'y' not in item:
                raise KeyError('y coordinate not provided')
            self.y = item['y']

            self.marker     = item['marker']     if 'marker'     in item else 'x'
            self.markersize = item['markersize'] if 'markersize' in item else None

        def plot(self):
            self.plotter.ax.plot(self.x, self.y,
                alpha=self.alpha, color=self.color,
                label=self.label, linestyle='None',
                marker=self.marker, markeredgecolor=self.color, markerfacecolor='None',
                markersize=self.markersize)


    class ErrorBar(BasePlot):
        """Plots a single errors bar"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'x' not in item:
                raise KeyError('x coordinate not provided')
            self.x    = item['x']
            self.xerr = item['xerr'] if 'xerr' in item else None

            if 'y' not in item:
                raise KeyError('y coordinate not provided')
            self.y    = item['y']
            self.yerr = item['yerr'] if 'yerr' in item else None

            self.elinewidth = item['elinewidth'] if 'elinewidth' in item else 1.0

        def plot(self):
            self.plotter.ax.errorbar(x=self.x, y=self.y, xerr=self.xerr, yerr=self.yerr,
                color=self.color, elinewidth=self.elinewidth, fmt='_', linestyle='none', label=self.label)


    class Band(BasePlot):
        """Plots a shaded band"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'x' not in item and 'y' not in item:
                raise KeyError('neither x nor y coordinates provided')

            self.x = tuple(item['x']) if 'x' in item else self.plotter.xrange
            self.y = tuple(item['y']) if 'y' in item else self.plotter.yrange
            self.xmin, self.xmax = self.x
            self.ymin, self.ymax = self.y

        def plot(self):
            rect = plt.Rectangle((self.xmin, self.ymin), self.xmax - self.xmin, self.ymax - self.ymin,
                                 alpha=self.alpha, color=self.color,
                                 linewidth=0, fill=True)
            self.plotter.ax.add_patch(rect)

        def handles_labels(self):
            if self.label:
                handle = plt.Rectangle((0,0),1,1, color=self.color)
                label  = self.label
                return ([handle], [label])
            else:
                return ([], [])


    class Graph(BasePlot):
        """Plots the graph of a function"""

        _api_doc = inspect.cleandoc("""\
        Plotting a Function Graph
        -------------------------

        The graph of a function can be easily plotted by providing the coordinates (x, f(x)) of the function.
        The coordinates are provided using a ``data`` object containing:

         * ``xvalues`` (array-like of *float*) -- The values on the x axis at which the function has been evaluated.
         * ``yvalues`` (array-like of *float*) -- The values oof the function at the points provided by ``xvalues``.

        Example:

        .. code-block::

           xvalues = np.linspace(0, 5, 20)
           plot_args = {
               'plot': {
                   'x': { 'label': r'$q^2$' },
                   'y': { 'label': r'$q^4$' }
               },
               'contents': [
                   {
                       'label': r'$\\ell=\\mu$',
                       'type': 'graph',
                       'data': { 'yvalues': xvalues**2, 'xvalues': xvalues },
                       'range': [0, 5],
                   },
               ]
           }
        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item:
                raise KeyError('No data was provided')

            self.xvalues = np.array(item['data']['xvalues'])
            self.yvalues = np.array(item['data']['yvalues'])
            self.xrange = item['range'] if 'range' in item else None

        def plot(self):
            if self.xrange:
                xvalues = np.ma.masked_outside(self.xvalues, float(self.xrange[0]), float(self.xrange[1]))
                yvalues = np.ma.masked_array(self.yvalues, mask=xvalues.mask)
            else:
                xvalues = self.xvalues
                yvalues = self.yvalues

            self.plotter.ax.plot(xvalues, yvalues, alpha=self.alpha, color=self.color, ls=self.style, label=self.label, lw=self.lw)


    class Observable(BasePlot):
        """Plots a single EOS observable w/o uncertainties as a function of one kinematic variable or one parameter"""

        _api_doc = inspect.cleandoc("""\
        Plotting Observables
        --------------------

        Contents items of type ``observable`` are used to display one of the built-in `observables <../reference/observables.html>`_.
        The following keys are mandatory:

         * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable that will be plotted.
           Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_.
         * ``range`` (*list* or *tuple* of two *float*) --The tuple of [minimal, maximal] values of the specified kinematic variable
           for which the observable will be evaluated.
         * ``variable`` (*str*) -- The name of the kinematic or parameter variable to which the x axis will be mapped; see
           `the complete list of parameters <../reference/parameters.html>`_.

        Example:

        .. code-block::

           plot_args = {
               'plot': { ... },
               'contents': [
                   {
                       'label': r'$\\ell=\\mu$',
                       'type': 'observable',
                       'observable': 'B->Dlnu::dBR/dq2;l=mu',
                       'variable': 'q2',
                       'range': [0.02, 11.60],
                   },
               ]
           }
        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            oname = item['observable']

            obs_entry = eos.Observables()[oname]
            valid_kin_vars = [kv for kv in obs_entry.kinematic_variables()]
            eos.info(f'   plotting EOS observable "{oname}"')

            # create kinematics
            kinematics = eos.Kinematics()
            if 'kinematics' in item:
                for k, v in item['kinematics'].items():
                    if k not in valid_kin_vars:
                        raise ValueError("Kinematic quantity '" + k + "' does not " +
                        "match known ones for observable '" + oname + "': " + valid_kin_vars.__repr__())
                    kinematics.declare(k, v)

            # create parameters
            parameters = eos.Parameters.Defaults()
            if 'parameters' in item and 'parameters-from-file' in item:
                eos.warn('    overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            if 'parameters' in item and 'parameters-from-mode' in item:
                eos.warn('    overriding values read from \'parameters-from-mode\' with explicit values in \'parameters\'')

            if 'parameters-from-file' in item and type(item['parameters-from-file']) is str:
                eos.warn('    overriding parameters from file')
                parameters.override_from_file(item['parameters-from-file'])

            if 'parameters-from-mode' in item and type(item['parameters-from-mode']) is str:
                eos.warn('    overriding parameters from mode')
                mode = eos.Mode(item['parameters-from-mode'])
                for p, v in zip(mode.varied_parameters, mode.mode):
                    parameters.set(p['name'], v)

            if 'parameters' in item and type(item['parameters']) is dict:
                for key, value in item['parameters'].items():
                    parameters.set(key, value)

            # create options
            options = eos.Options()
            if 'options' in item and type(item['options']) is dict:
                for key, value in item['options'].items():
                    options.declare(key, value)


            #
            # handle plot variable, create kinematic or parameter
            #

            if not 'variable' in item:
                raise ValueError("Missing key for plot of observable '" + oname + "': 'variable'")

            # variable is passed *again* as kinematic or parameter?
            if 'parameters' in item and item['variable'] in item['parameters']:
                val = item['parameters'].get(item['variable'])
                raise ValueError("Variable '" + item['variable'] + "' of observable '" + oname + "' is " +
                        "also specified as a fix parameter with value " + str(val))
            if 'kinematics' in item and item['variable'] in item['kinematics']:
                val = item['kinematics'].get(item['variable'])
                raise ValueError("Variable '" + item['variable'] + "' of observable '" + oname + "' is " +
                        "also specified as a fix kinematic with value " + str(val))

            # Declare variable that is either parameter or kinematic
            var = None

            # does the variable correspond to one of the kinematic variables?
            if item['variable'] in valid_kin_vars:
                var = kinematics.declare(item['variable'], np.nan)
            else: # if not,
                # is the variable name a QualifiedName?
                try:
                    qn = eos.QualifiedName(item['variable'])
                    if not parameters.has(qn):
                        var = parameters.declare(qn, '', eos.Unit.Undefined(), np.nan, -np.max, np.max)
                    else:
                        var = parameters[qn]
                except RuntimeError:
                    raise ValueError("Value of 'variable' for observable '" + oname +
                        "' is neither a valid kinematic variable nor parameter: '" + item['variable'] + "'")

            # create observable
            observable = eos.Observable.make(oname, parameters, kinematics, options)

            xvalues = np.linspace(self.xlo, self.xhi, self.xsamples + 1)
            ovalues = np.array([])
            for xvalue in xvalues:
                var.set(xvalue)
                ovalues = np.append(ovalues, observable.evaluate())

            self.plotter.ax.plot(xvalues, ovalues, alpha=self.alpha, color=self.color, label=self.label, ls=self.style, lw=self.lw)

    class Uncertainty(BasePlot):
        """Plots an uncertainty band as a function of one kinematic variable

        This routine expects the uncertainty propagation to have produced an EOS data file"""

        _api_doc = inspect.cleandoc("""\
        Plotting Uncertainty Bands
        --------------------------

        Contents items of type ``uncertainty`` are used to display uncertainty bands at 68% probability
        for one of the built-in `observables <../reference/observables.html>`_. This item type requires that either
        a predictive distribution of the observables has been previously produced.

        Exactly one of the following keys is mandatory:

         * ``data`` (*dict*, see below) -- The data on predictive distribution of the observable whose uncertainty band will be plotted.

         * ``data-file`` (*str*, path to an existing data file of type *eos.data.Prediction*) -- The path to
           a data file that was generated with the ``eos-analysis`` command-line client.

        For ``data`` object, the following keys are mandatory:

         * ``xvalues`` (array-like of *float*) -- The values on the x axis at which the observable has been evaluated.
         * ``samples`` (*list* of tuples of *float*) -- The list of samples of the predictive distribution. Each tuple of samples
           corresponds to an evaluation of the observables at the kinematic configuration corresponding to the ``xvalues``
           entry with the same index.
         * ``weights`` (array-like of *float*, optional) -- The weights of the samples, on a linear scale. Defaults to uniform weights.

        The following keys are optional:

         * ``band`` (a *list* containing ``'lines'``, ``'areas'``, or both) -- The setting for the illustration of the band.
           If ``'outer'`` is provided, the band's outer lines are drawn.
           If ``'median'`` is provided, the band's median line is drawn.
           If ``'area'`` is provided, the band's areas are filled.
           Defaults to ``['area', 'outer', 'median']``.
         * ``interpolation-type`` (*str*) -- The type of interpolation to be used for the band. Can be either ``linear`` (default) or ``cubic``.
         * ``plot-data`` (*bool*) -- If set to ``True``, the data points are plotted on top of the band.

        Example:

        .. code-block::

           analysis = ... # eos.Analysis object as discussed in the example notebook `predictions.ipynb`
           mu_q2values = numpy.unique(numpy.concatenate((numpy.linspace(0.02,1.00,20), numpy.linspace(1.00,11.60,20))))
           mu_obs      = [eos.Observable.make('B->Dlnu::dBR/dq2', prior.parameters,
                                              eos.Kinematics(q2=q2), eos.Options({'form-factors':'BSZ2015','l':'mu'}))
                          for q2 in mu_q2values]
           _, _, mu_samples = analysis.sample(N=5000, pre_N=1000, observables=mu_obs)
           plot_args = {
               'plot': {
                   'x': { 'label': r'$q^2$' },
                   'y': { 'label': r'$d\\mathcal{B}/dq^2$' }
               },
               'contents': [
                   {
                       'label': r'$\\ell=\\mu$',
                       'type': 'uncertainty',
                       'data': { 'samples': mu_samples, 'xvalues': mu_q2values },
                       'range': [0.02, 11.60],
                   },
               ]
           }
        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'data-file' not in item:
                raise KeyError('neither data nor data-file specified')

            if 'data' in item and 'data-file' in item:
                eos.warn('   both data and data-file specified; assuming interactive use is intended')

            self.xvalues = None
            self.samples = None
            self.weights = None
            if 'data' in item:
                self.xvalues = np.array(item['data']['xvalues'])
                self.samples = item['data']['samples']
                if 'weights' in item['data']:
                    self.weights = item['data']['weights']
            else:
                dfname = item['data-file']
                eos.info(f'   plotting uncertainty propagation from "{dfname}"')
                df = eos.data.Prediction(dfname)
                _xvalues = []
                for vp in df.varied_parameters:
                    if len(vp['kinematics']) > 1:
                        raise ValueError('more than one kinematic variable specified')
                    _xvalues.append(list(vp['kinematics'].values())[0])

                self.xvalues = np.array(_xvalues)
                self.samples = df.samples
                self.weights = df.weights

            self.band   = item['band']  if 'band'  in item else ['area', 'outer', 'median']
            self.xrange = item['range'] if 'range' in item else None

            self.interpolation_type = item['interpolation-type'] if 'interpolation-type' in item else 'linear'
            self.plot_data = item['plot-data'] if 'plot-data' in item else False

        def plot(self):
            _ovalues_lower   = []
            _ovalues_central = []
            _ovalues_higher  = []
            for i in range(len(self.samples[0])):
                lower, central, higher = self.plotter._weighted_quantiles(self.samples[:, i],
                                                                         [0.15865, 0.5, 0.84135],
                                                                         self.weights)
                _ovalues_lower.append(lower)
                _ovalues_central.append(central)
                _ovalues_higher.append(higher)

            xvalues = np.linspace(np.min(self.xvalues), np.max(self.xvalues), self.xsamples)
            if self.xrange:
                xvalues = np.ma.masked_outside(xvalues, float(self.xrange[0]), float(self.xrange[1]))

            if self.interpolation_type == "linear":
                interpolate = lambda x, y, xv: np.interp(xv, x, y)
            elif self.interpolation_type == "cubic":
                # work around CubicSpline missing in SciPy version < 0.18
                if scipy.__version__ >= '0.18':
                    from scipy.interpolate import CubicSpline
                    interpolate = lambda x, y, xv: CubicSpline(x, y)(xv)
                else:
                    from scipy.interpolate import spline
                    interpolate = spline

            ovalues_lower   = interpolate(self.xvalues, _ovalues_lower,   xvalues)
            ovalues_central = interpolate(self.xvalues, _ovalues_central, xvalues)
            ovalues_higher  = interpolate(self.xvalues, _ovalues_higher,  xvalues)

            if 'area' in self.band:
                self.plotter.ax.fill_between(xvalues, ovalues_lower, ovalues_higher, alpha=self.alpha, color=self.color, label=self.label, lw=0)
                if 'outer' in self.band:
                    self.plotter.ax.plot(xvalues, ovalues_lower,                     alpha=self.alpha, color=self.color, ls=self.style, lw=self.lw)
                    self.plotter.ax.plot(xvalues, ovalues_higher,                    alpha=self.alpha, color=self.color, ls=self.style, lw=self.lw)
                if 'median' in self.band:
                    self.plotter.ax.plot(xvalues, ovalues_central,                   alpha=self.alpha, color=self.color, ls=self.style, lw=self.lw)
            elif 'outer' in self.band:
                self.plotter.ax.plot(xvalues, ovalues_lower,                         alpha=self.alpha, color=self.color, ls=self.style, label=self.label, lw=self.lw)
                self.plotter.ax.plot(xvalues, ovalues_higher,                        alpha=self.alpha, color=self.color, ls=self.style, lw=self.lw)
                if 'median' in self.band:
                    self.plotter.ax.plot(xvalues, ovalues_central,                   alpha=self.alpha, color=self.color, ls=self.style, lw=self.lw)
            elif 'median' in self.band:
                self.plotter.ax.plot(xvalues, ovalues_central,                       alpha=self.alpha, color=self.color, ls=self.style, label=self.label, lw=self.lw)

            if self.plot_data:
                self.plotter.ax.scatter(self.xvalues, _ovalues_central,              alpha=self.alpha, color=self.color, lw=self.lw)


    class UncertaintyBinned(BasePlot):
        """Plots one or more uncertainty band integrated over one kinematic variable

        This routine expects the uncertainty propagation to have produced an data file"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data-file' not in item and 'data' not in item:
                raise KeyError('neither data nor data-file specified')


            self.weights = None
            if 'data' in item:
                self.xvalues = item['data']['xvalues']
                self.samples = item['data']['samples']
                if 'weights' in item['data']:
                    self.weights = item['data']['weights']
            else:
                if 'variable' not in item:
                    raise KeyError('\'variable\' is mandatory when using \'data-file\' specified')
                variable = item['variable']
                data = eos.data.Prediction(item['data-file'])
                try:
                    self.xvalues = np.array([(p['kinematics'][variable + '_min'], p['kinematics'][variable + '_max']) for p in data.varied_parameters])
                except KeyError as e:
                    raise RuntimeError(f'both \'{variable}_min\' and \'{variable}_max\' must be present in the kinematics of each prediction in the data file') from e
                self.samples = data.samples
                self.weights = data.weights

            self.xvalues = np.array(self.xvalues)

            if 'range' in item:
                xmin, xmax   = item['range']
                xmin, xmax   = (float(xmin), float(xmax))
                self.xvalues = np.ma.masked_outside(self.xvalues, xmin, xmax)

            self.rescale_by_width = item['rescale-by-width'] if 'rescale-by-width' in item else False


        def plot(self):
            item = self.item

            ovalues_lower   = []
            ovalues_central = []
            ovalues_higher  = []
            for i in range(len(self.xvalues)):
                lower, central, higher = self.plotter._weighted_quantiles(self.samples[:, i],
                                                                         [0.15865, 0.5, 0.84135],
                                                                         self.weights)
                ovalues_lower.append(lower)
                ovalues_central.append(central)
                ovalues_higher.append(higher)

            for [xmin, xmax], olo, ocentral, ohi in zip(self.xvalues, ovalues_lower, ovalues_central, ovalues_higher):
                width = (xmax - xmin) if self.rescale_by_width else 1
                olo      /= width
                ocentral /= width
                ohi      /= width
                print(f"{xmin} ... {xmax} -> {ocentral} with interval {olo} .. {ohi}")
                self.plotter.ax.fill_between([xmin, xmax], [olo, olo], [ohi, ohi], lw=0, color=self.color, alpha=self.alpha, label=self.label)
                label = None
                self.plotter.ax.plot([xmin, xmax], [olo,      olo],      color=self.color, alpha=self.alpha)
                self.plotter.ax.plot([xmin, xmax], [ocentral, ocentral], color=self.color, alpha=self.alpha)
                self.plotter.ax.plot([xmin, xmax], [ohi,      ohi],      color=self.color, alpha=self.alpha)
                self.label = None


    class UncertaintyOverview(BasePlot):
        """Plots an overview of uncertainty estimates

        This routine expects the uncertainty propagation to have produced an HDF5 file"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'hdf5-file' not in item:
                raise KeyError('hdf5-file not specified')

            if 'observables' not in item:
                raise KeyError('observables not specified')

            h5fname = item['hdf5-file']
            eos.info(f'   plotting uncertainty propagation from file "{h5fname}"')
            uncfile = eos.data.UncertaintyDataFile(h5fname)

            self.observables    = item['observables']
            self.observable_map = {}
            for idx, p in enumerate(uncfile.parameters):
                name = p[0].decode('ascii')
                self.observable_map[name] = idx

            for observable in self.observables:
                if observable not in self.observable_map:
                    raise ValueError(f'observable \'{observable}\' not contained in HDF5 file')

            self.samples = uncfile.data()

        def plot(self):
            for xvalue, observable in enumerate(self.observables):
                data_idx = self.observable_map[observable]
                lower    = np.percentile(self.samples[:, data_idx], q=15.865)
                upper    = np.percentile(self.samples[:, data_idx], q=84.135)

                self.plotter.ax.fill_between([xvalue - 0.5, xvalue + 0.5], [lower, lower], [upper, upper], alpha=self.alpha, color=self.color, label=self.label, lw=0)
                self.label = None


    class Constraint(BasePlot):
        """Plots constraints from the EOS library of experimental and theoretical likelihoods"""

        _api_doc = inspect.cleandoc("""\
        Plotting Constraints
        --------------------

        Contents items of type ``constraints`` are used to display one of the built-in `experimental or theoretical constraints <../reference/constraints.html>`_.
        The following keys are mandatory:

         * ``constraints`` (:class:`QualifiedName <eos.QualifiedName>` or iterable thereof) -- The name or the list of names of the constraints
           that will be plotted. Must identify at least one of the constraints known to EOS; see `the complete list of constraints <../reference/constraints.html>`_.
         * ``variable`` (*str*) -- The name of the kinematic variable to which the x axis will be mapped.

        When plotting multivariate constraints, the following key is also mandatory:

         * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable whose constraints will be plotted.
           Must identify one of the observables known to EOS; see `the complete list of observables <../reference/observables.html>`_.
           This is only mandatory in multivariate constraints, since these can constrain more than one observable simultaneously.

        The following keys are optional:

         * ``xrange`` (list of int) -- The interval in which the observable is plotted in the case of a multivariate constraint.
         * ``rescale-by-width`` (*bool*) -- Rescales binned constraints by the inverse of the bin width. This is often required
           to compare theory (integrated) predictions and experimental (averaged) measurements. Defaults to false.

        Example:

        .. code-block::

           plot_args = {
               'plot': { ... },
               'contents': [
                   {
                       'label': r'Belle 2015 $\\ell=e,\\, q=d$',
                       'type': 'constraint',
                       'color': 'C0',
                       'constraints': 'B^0->D^+e^-nu::BRs@Belle:2015A',
                       'observable': 'B->Dlnu::BR',
                       'variable': 'q2',
                       'rescale-by-width': False
                   }
               ]
           }

        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'constraints' not in item:
                raise KeyError('no constraints specified')

            if 'variable' not in item:
                raise KeyError('no variable specified')

            # extract information
            self.names            = item['constraints']
            self.observable       = item['observable']       if 'observable'       in item else None
            self.rescale_by_width = item['rescale-by-width'] if 'rescale-by-width' in item else False
            self.variable         = item['variable']
            self.xrange           = item['xrange']           if 'xrange'           in item else None
            self.plot_residues    = item['plot_residues']    if 'plot_residues'    in item else False

            if type(self.names) == str:
                self.names = [self.names]

            if self.plot_residues:
                # create parameters
                self.parameters = eos.Parameters.Defaults()
                if 'parameters' in item and 'parameters-from-mode' in item:
                    eos.warn('    overriding values read from \'parameters-from-mode\' with explicit values in \'parameters\'')

                if 'parameters-from-mode' in item and type(item['parameters-from-mode']) is str:
                    eos.warn('    overriding parameters from mode')
                    mode = eos.Mode(item['parameters-from-mode'])
                    for p, v in zip(mode.varied_parameters, mode.mode):
                        self.parameters.set(p['name'], v)

                if 'parameters' in item and type(item['parameters']) is dict:
                    for key, value in item['parameters'].items():
                        self.parameters.set(key, value)

                # create options
                self.options = eos.Options()
                if 'options' in item and type(item['options']) is dict:
                    for key, value in item['options'].items():
                        self.options.declare(key, value)


        def plot(self):
            import yaml

            constraints = eos.Constraints()
            obs_kinematics = eos.Kinematics()

            for name in self.names:
                entry = constraints[name]
                if not entry:
                    raise ValueError(f'unknown constraint {name}')

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
                        xerrors = None
                        if self.plot_residues:
                            obs_values = []
                            for x in xvalues:
                                obs_kinematics.declare(self.variable, x)
                                obs_values.append(eos.Observable.make(self.observable, self.parameters, obs_kinematics, self.options).evaluate())

                    elif (self.variable + '_min' in kinematics) and (self.variable + '_max' in kinematics):
                        xvalues = [(kinematics[self.variable + '_max'] + kinematics[self.variable + '_min']) / 2]
                        xerrors = [(kinematics[self.variable + '_max'] - kinematics[self.variable + '_min']) / 2]
                        if self.rescale_by_width:
                            width = (kinematics[self.variable + '_max'] - kinematics[self.variable + '_min'])
                        if self.plot_residues:
                            obs_values = []
                            for xmax, xmin in zip(kinematics[self.variable + '_max'], kinematics[self.variable + '_min']):
                                obs_kinematics.declare(self.variable + '_max', xmax)
                                obs_kinematics.declare(self.variable + '_min', xmin)
                                obs_values.append(eos.Observable.make(self.observable, self.parameters, obs_kinematics, self.options).evaluate() / width)

                    yvalues = [float(constraint['mean']) / width]
                    sigma_hi = np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2) / width
                    sigma_lo = np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2) / width
                    yerrors = [[sigma_hi, sigma_lo]]
                elif constraint['type'] == 'MultivariateGaussian(Covariance)':
                    if not self.observable:
                        raise KeyError('observable needs to be specified for MultivariateGaussian(Covariance) constraints')
                    covariance = np.array(constraint['covariance'])
                    observables = constraint['observables']
                    means = constraint['means']
                    dim = len(means)
                    kinematics = constraint['kinematics']

                    xvalues = []
                    xerrors = []
                    yvalues = []
                    yerrors = []
                    obs_values = []
                    for i in range(0, dim):
                        width = 1

                        if not observables[i] == self.observable:
                            continue
                        _kinematics = kinematics[i]
                        if self.variable in _kinematics:
                            xvalues.append(_kinematics[self.variable])
                            xerrors = None
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
                            if self.plot_residues:
                                obs_kinematics.declare(self.variable + '_max', var_max)
                                obs_kinematics.declare(self.variable + '_min', var_min)
                                obs_values.append(eos.Observable.make(self.observable, self.parameters, obs_kinematics, self.options).evaluate() / width)

                        yvalues.append(np.double(means[i]) / width)
                        yerrors.append(np.sqrt(np.double(covariance[i, i])) / width)
                elif constraint['type'] == 'MultivariateGaussian':
                    if not self.observable:
                        raise KeyError('observable needs to be specified for MultivariateGaussian constraints')
                    sigma_stat_hi = np.array(constraint['sigma-stat-hi'])
                    sigma_stat_lo = np.array(constraint['sigma-stat-lo'])
                    sigma_sys = np.array(constraint['sigma-sys'])
                    sigma = np.sqrt(np.power(sigma_sys, 2) + 0.25 * np.power(sigma_stat_hi + sigma_stat_lo, 2))
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
                            xerrors = None
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
                    raise ValueError('type of constraint presently not supported')

                xvalues = np.array(xvalues)
                yvalues = np.array(yvalues)
                yerrors = np.array(yerrors)

                if len(xvalues) == 0:
                    eos.info(f'   skipping plot for constraint {name} since it does not contain the requested observable')
                    return

                if self.plot_residues:
                    yvalues -= obs_values

                if self.xrange:
                    mask = np.logical_and(xvalues > min(self.xrange), xvalues < max(self.xrange))
                else:
                    mask = np.array([True] * len(xvalues))

                if xerrors:
                    xerrors = np.array(xerrors)
                    self.plotter.ax.errorbar(x=xvalues[mask], y=yvalues[mask], xerr=xerrors[mask], yerr=yerrors[mask].T,
                        color=self.color, elinewidth=1.0, fmt='_', linestyle='none', label=self.label)
                else:
                    self.plotter.ax.errorbar(x=xvalues[mask], y=yvalues[mask], yerr=yerrors[mask].T,
                        color=self.color, elinewidth=1.0, fmt='_', linestyle='none', label=self.label)
                # disable the label for subsequent plots
                self.label = None


    class Constraint2D(BasePlot):
        """Plot 2D contours for two correlated observables from a single constraint"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            # mandatory keys
            if 'constraint' not in item:
                raise KeyError('no constraint specified')
            self.constraint = item['constraint']

            self.sigmas = item['sigmas'] if 'sigmas' in item else [1.0]
            self.sigmas.sort(reverse=True)
            self.alphas = np.linspace(0.0, self.alpha, len(self.sigmas) + 1)[1:]

        def plot(self):
            item = self.item
            import yaml

            constraints = eos.Constraints()
            entry = constraints[self.constraint]
            if not entry:
                raise ValueError(f'unknown constraint {constraint}')

            constraint = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)
            if constraint['type'] not in ['Gaussian', 'MultivariateGaussian', 'MultivariateGaussian(Covariance)']:
                raise ValueError('constraint type \'{}\' presently not supported'.format(constraint['type']))

            if constraint['type'] in ['MultivariateGaussian(Covariance)', 'MultivariateGaussian']:
                xmean  = None
                xsigma = None
                ymean  = None
                ysigma = None
                rho    = None

                if 'x' not in item:
                    raise KeyError('no x-axis specification for 2D constraint')

                xitem = item['x']
                if 'observable' not in xitem:
                    raise KeyError('no observable in 2D constraint\'s x-axis specification')
                self.xobservable = xitem['observable']
                self.xkinematics = xitem['kinematics'] if 'kinematics' in xitem else {}
                self.xoptions    = xitem['options']    if 'options'    in xitem else {}

                if 'y' not in item:
                    raise KeyError('no y-axis specification for 2D constraint')

                yitem = item['y']
                if 'observable' not in yitem:
                    raise KeyError('no observable in 2D constraint\'s y-axis specification')
                self.yobservable = yitem['observable']
                self.ykinematics = yitem['kinematics'] if 'kinematics' in yitem else {}
                self.yoptions    = yitem['options']    if 'options'    in yitem else {}

                if constraint['type'] == 'MultivariateGaussian(Covariance)':
                    observables = list(zip(constraint['observables'], constraint['kinematics'], constraint['options']))

                    xobservable = (self.xobservable, self.xkinematics, self.xoptions)
                    yobservable = (self.yobservable, self.ykinematics, self.yoptions)

                    if xobservable not in observables:
                        raise ValueError(f'x-axis observable {self.xobservable} not contained in constraint {self.constraint}')

                    if yobservable not in observables:
                        raise ValueError(f'y-axis observable {self.yobservable} not contained in constraint {self.constraint}')

                    xidx = observables.index(xobservable)
                    yidx = observables.index(yobservable)

                    means  = np.array(constraint['means'])
                    xmean  = means[xidx]
                    ymean  = means[yidx]

                    cov    = np.array(constraint['covariance'])
                    xsigma = np.sqrt(covariance[xidx, xidx])
                    ysigma = np.sqrt(covariance[yidx, yidx])
                    rho    = cov[xidx, yidx] / (xsigma * ysigma)
                elif constraint['type'] == 'MultivariateGaussian':
                    observables = list(zip(constraint['observables'], constraint['kinematics'], constraint['options']))

                    xobservable = (self.xobservable, self.xkinematics, self.xoptions)
                    yobservable = (self.yobservable, self.ykinematics, self.yoptions)

                    if xobservable not in observables:
                        raise ValueError(f'x-axis observable {self.xobservable} not contained in constraint {self.constraint}')

                    if yobservable not in observables:
                        raise ValueError(f'y-axis observable {self.yobservable} not contained in constraint {self.constraint}')

                    xidx = observables.index(xobservable)
                    yidx = observables.index(yobservable)

                    means  = np.array(constraint['means'])
                    xmean  = means[xidx]
                    ymean  = means[yidx]

                    sigmas = np.sqrt(
                          (np.abs(np.array(constraint['sigma-stat-hi'])) + np.abs(np.array(constraint['sigma-stat-lo'])))**2 / 4.0
                        +  np.array(constraint['sigma-sys'])**2
                    )
                    xsigma = sigmas[xidx]
                    ysigma = sigmas[yidx]

                    rho    = np.array(constraint['correlations'])[xidx, yidx]

                for sigma, alpha in zip(self.sigmas, self.alphas):
                    xwidth = 2.0 * np.sqrt(1.0 + rho)
                    ywidth = 2.0 * np.sqrt(1.0 - rho)

                    ellipse = matplotlib.patches.Ellipse((0.0, 0.0), width=xwidth, height=ywidth,
                                          alpha=alpha, color=self.color, linewidth=1, fill=True)
                    transf = matplotlib.transforms.Affine2D() \
                        .rotate_deg(45) \
                        .scale(xsigma * sigma, ysigma * sigma) \
                        .translate(xmean, ymean)
                    ellipse.set_transform(transf + self.plotter.ax.transData)
                    self.plotter.ax.add_patch(ellipse)
            elif constraint['type'] == 'Gaussian':
                if 'x' in self.item and 'y' in self.item:
                    raise ValueError('both \'x\' and \'y\' specified for univariate constraint')

                if 'x' in self.item:
                    xitem = self.item['x']

                    if 'observable' not in xitem:
                        raise KeyError('no observable in 2D constraint\'s x-axis specification')

                    self.xobservable = xitem['observable']
                    self.xkinematics = xitem['kinematics'] if 'kinematics' in xitem else {}
                    self.xoptions    = xitem['options']    if 'options'    in xitem else {}

                    observable  = (constraint['observable'], constraint['kinematics'], constraint['options'])
                    xobservable = (self.xobservable, self.xkinematics, self.xoptions)

                    if xobservable != observable:
                        raise ValueError(f'x-axis observable {self.xobservable} not contained in constraint {self.constraint}')

                    mean     = float(constraint['mean'])
                    sigma_hi = np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2)
                    sigma_lo = np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2)

                    for sigma, alpha in zip(self.sigmas, self.alphas):
                        xlo, xhi = (mean - sigma * sigma_lo, mean + sigma * sigma_hi)
                        ylo, yhi = self.plotter.yrange

                        rect = matplotlib.patches.Rectangle((xlo, ylo), width=(xhi - xlo) * sigma, height=yhi - ylo,
                                                       alpha=alpha, color=self.color)
                        self.plotter.ax.add_patch(rect)
                elif 'y' in self.item:
                    yitem = self.item['y']

                    if 'observable' not in yitem:
                        raise KeyError('no observable in 2D constraint\'s y-axis specification')

                    self.yobservable = yitem['observable']
                    self.ykinematics = yitem['kinematics'] if 'kinematics' in yitem else {}
                    self.yoptions    = yitem['options']    if 'options'    in yitem else {}

                    observable  = (constraint['observable'], constraint['kinematics'], constraint['options'])
                    yobservable = (self.yobservable, self.ykinematics, self.yoptions)

                    if yobservable != observable:
                        raise ValueError(f'y-axis observable {self.yobservable} not contained in constraint {self.constraint}')

                    mean     = float(constraint['mean'])
                    sigma_hi = np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2)
                    sigma_lo = np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2)

                    for sigma, alpha in zip(self.sigmas, self.alphas):
                        ylo, yhi = (mean - sigma * sigma_lo, mean + sigma * sigma_hi)
                        xlo, xhi = self.plotter.xrange

                        rect = matplotlib.patches.Rectangle((xlo, ylo), width=xhi - xlo, height=(yhi - ylo) * sigma,
                                                       alpha=alpha, color=self.color)
                        self.plotter.ax.add_patch(rect)


        def handles_labels(self):
            if self.label:
                handle = plt.Rectangle((0,0),1,1, color=self.color, alpha=self.alpha)
                label  = self.label
                return ([handle], [label])
            else:
                return ([], [])


    class ConstraintOverview(BasePlot):
        """Plots overview of several constraints from the EOS library of experimental and theoretical likelihoods"""
        def __init__(self, plotter, item):
            import yaml

            super().__init__(plotter, item)

            if 'constraints' not in item:
                raise KeyError('no constraints specified')

            # extract information
            self.names            = item['constraints']
            self.rotation         = 'vertical' if 'rotation' not in item else item['rotation']
            self.constraints      = []

            if type(self.names) == str:
                self.names = [self.names]

            constraints = eos.Constraints()

            for name in self.names:
                entry = constraints[name]
                if not entry:
                    raise ValueError(f'unknown constraint {name}')

                constraint = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)
                self.constraints.append(constraint)

        def plot(self):
            xvalues     = []
            xticklabels = []
            yvalues     = []
            yerrors     = []
            idx = 0
            for constraint in self.constraints:
                if constraint['type'] == 'Gaussian':
                    sigma_hi = np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2) / width
                    sigma_lo = np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2) / width
                    latex = '$' + eos.Observables()[constraint['observable']].latex() + '$'
                    xvalues.append(idx)
                    xticklabels.append(latex)
                    yvalues.append(constraint['mean'])
                    yerrors.append([sigma_hi, sigma_lo])
                    idx = idx + 1
                elif constraint['type'] == 'MultivariateGaussian(Covariance)':
                    observables = constraint['observables']
                    kinematics  = constraint['kinematics']
                    options     = constraint['options']
                    covariance  = np.array(constraint['covariance'])
                    means       = constraint['means']
                    dim         = len(means)

                    for i in range(0, dim):
                        latex = '$' + eos.Observables()[observables[i]].latex() + '$'
                        sigma = np.sqrt(covariance[i, i])
                        xvalues.append(idx)
                        xticklabels.append(latex)
                        yvalues.append(means[i])
                        yerrors.append([sigma, sigma])
                        idx = idx + 1
                elif constraint['type'] == 'MultivariateGaussian':
                    observables = constraint['observables']
                    sigma_stat_hi = np.array(constraint['sigma-stat-hi'])
                    sigma_stat_lo = np.array(constraint['sigma-stat-lo'])
                    sigma_sys = np.array(constraint['sigma-sys'])
                    sigma = np.sqrt(np.power(sigma_sys, 2) + 0.25 * np.power(sigma_stat_hi + sigma_stat_lo, 2))
                    means = constraint['means']
                    dim = len(means)

                    for i in range(0, dim):
                        latex = '$' + eos.Observables()[observables[i]].latex() + '$'
                        xvalues.append(idx)
                        xticklabels.append(latex)
                        yvalues.append(means[i])
                        yerrors.append(sigma[i])
                        idx = idx + 1
                else:
                    raise ValueError('type of constraint presently not supported')

            xvalues = np.array(xvalues)
            yvalues = np.array(yvalues)
            yerrors = np.array(yerrors)

            self.plotter.ax.tick_params(axis='x', which='minor', bottom=False)
            self.plotter.ax.xticks(xvalues, xticklabels, rotation=self.rotation)
            self.plotter.ax.errorbar(x=xvalues, y=yvalues, xerr=None, yerr=yerrors.T,
                color=self.color, elinewidth=1.0, fmt='_', linestyle='none', label=self.label)
            self.plotter.ax.margins(0.2)
            # Tweak spacing to prevent clipping of tick-labels
            self.plotter.ax.subplots_adjust(bottom=0.15)


    class Contours2D(BasePlot):
        """Plots 2D contours of a pair of parameters based on pre-existing random samples"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            if 'hdf5-file' not in item:
                raise KeyError('no hdf5-file specified')

            h5fname = item['hdf5-file']
            eos.info(f'   plotting 2D contours from file "{h5fname}"')
            datafile = eos.data.load_data_file(h5fname)

            if 'variables' not in item:
                raise KeyError('no variables specificed')

            xvariable, yvariable = item['variables']

            data   = datafile.data()
            xindex = datafile.variable_indices[xvariable]
            xdata  = data[:, xindex]
            yindex = datafile.variable_indices[yvariable]
            ydata  = data[:, yindex]

            if not np.array(self.plotter.xrange).any():
                self.plotter.xrange = [np.amin(xdata), np.amax(xdata)]
                self.plotter.ax.set_xlim(tuple(self.plotter.xrange))
            if not np.array(self.plotter.yrange).any():
                self.plotter.yrange = [np.amin(ydata), np.amax(ydata)]
                self.plotter.ax.set_ylim(tuple(self.plotter.yrange))
            plt.show()

            xbins = 100
            ybins = 100

            H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(xbins, ybins), normed=True)
            x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,xbins))
            y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((ybins,1))
            pdf = (H*(x_bin_sizes*y_bin_sizes))

            # find the PDF value corresponding to a given cummulative probability
            plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
            pone_sigma = scipy.optimize.brentq(plevel, 0., 1., args=(pdf, 0.68))
            ptwo_sigma = scipy.optimize.brentq(plevel, 0., 1., args=(pdf, 0.95))
            pthree_sigma = scipy.optimize.brentq(plevel, 0., 1., args=(pdf, 0.99))
            levels = [pone_sigma, ptwo_sigma, pthree_sigma]
            labels = ['68%', '95%', '99%']

            CS = self.plotter.ax.contour(pdf.transpose(),
                             colors='OrangeRed',
                             extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                             levels=levels[::-1])

            fmt = {}
            for level, label in zip(CS.levels, labels[::-1]):
                fmt[level] = label

            self.plotter.ax.clabel(CS, inline=1, fmt=fmt, fontsize=10)


    class KernelDensityEstimate1D(BasePlot):
        """Plots a 1D Kernel Density Estimate (KDE) of pre-existing random samples"""

        _api_doc = inspect.cleandoc("""\
        Plotting Kernel Density Estimates
        ---------------------------------

        Contents items of type ``kde`` are used to display a kernel density estimate (a smooth histogram) of samples of a probability density,
        be it a prior, a posterior, or a signal PDF.

        The following key is mandatory:

         * ``data`` (*dict*, see below) -- The data on the probability density that will be histogramed.

           Within the data object, the following keys are understood.

            * ``samples`` (*list* of *float*) -- The samples that will be histogramed. Mandatory.
            * ``weights`` or ``log_weights`` (*list* of *float*, optional) -- The weights of the samples, on a linear or logarithmic scale.
              Defaults to uniform weights.


        The following keys are optional:

         * ``bandwidth`` (*float*) -- The factor by which the automatically determined kernel bandwidth is scaled. See the SciPy documentation
           for ``gaussian_kde``, ``bw_method='silverman'``. Defaults to 1.
         * ``range`` (*tuple* of two *float*) -- The minimum and maximum value of the x coordinate for which the smooth histogram is plotted.
         * ``level`` (*float*) -- The probability level of the contour. Defaults to ``68.27``.

        Example:

        .. code-block::

           analysis = ... # eos.Analysis object as discussed in the example notebook `inference.ipynb`
           parameter_samples, _, = analysis.sample(N=5000, pre_N=1000)
           plot_args = {
               'plot': {
                   'x': { 'label': r'$|V_{cb}|$', 'range':[38e-3, 45e-3] }
               },
               'contents': [
                   {
                       'type': 'kde', 'color': C0, 'label': 'posterior', 'bandwidth': 1,
                       'range': [40e-3, 45e-3],
                       'data': { 'samples': parameter_samples[:, 0] }
                   }
               ]
           }

        """)
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'data-file' not in item:
                raise KeyError('neither data nor data-file specified')

            if 'data' in item and 'data-file' in item:
                eos.warn('   both data and data-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples']

                if 'weights' in item['data'] and 'log_weights' in item['data']:
                    raise KeyError("Only one of 'weights' and 'log_weights' must be specified")
                elif 'weights' in item['data']:
                    self.weights = item['data']['weights']
                elif 'log_weights' in item['data']:
                    self.weights = np.exp(item['data']['log_weights'])
                else:
                    self.weights = None
            else:
                if 'variable' not in item:
                    raise KeyError('no variable specificed')

                dfname = item['data-file']
                eos.info(f'   plotting KDE for "{dfname}"')
                prefix = os.path.split(dfname)[-1]
                eos.info(f'   prefix = {prefix}')
                if prefix.startswith('mcmc-'):
                    df  = eos.data.MarkovChain(dfname)
                    idx = df.lookup_table[item['variable']]
                    self.samples = df.samples[:, idx]
                    self.weights = None
                elif prefix.startswith('samples'):
                    df = eos.data.ImportanceSamples(dfname)
                    idx = df.lookup_table[item['variable']]
                    self.samples = df.samples[:, idx]
                    self.weights = df.weights
                elif prefix.startswith('pred-'):
                    df = eos.data.Prediction(dfname)
                    idx = df.lookup_table[item['variable']]
                    self.samples = df.samples[:, idx]
                    self.weights = df.weights
                else:
                    raise ValueError(f'Do not recognize data-file prefix: {dfname}')

            self.bw       = item['bandwidth'] if 'bandwidth' in item else None
            self.level    = item['level']     if 'level'     in item else 68.27
            self.xrange   = item['range']     if 'range'     in item else plotter.xrange

        def plot(self):
            kde = gaussian_kde(self.samples, weights=self.weights)
            kde.set_bandwidth(bw_method='silverman')
            if self.bw:
                kde.set_bandwidth(bw_method=kde.factor * self.bw)

            x = np.linspace(self.xrange[0], self.xrange[1], self.xsamples)
            pdf = kde(x)
            pdf_norm = pdf.sum()

            # find the PDF value corresponding to a given cummulative probability
            if self.level:
                plevelf = lambda x, pdf, P: pdf[pdf > x * pdf_norm].sum() - P * pdf_norm
                plevel = scipy.optimize.brentq(plevelf, 0., 1., args=(pdf, self.level / 100.0 ))
                self.plotter.ax.fill_between(np.ma.masked_array(x, mask=pdf < plevel * pdf_norm),
                                             np.ma.masked_array(pdf, mask=pdf < plevel * pdf_norm, fill_value=np.nan),
                                             facecolor=self.color, alpha=self.alpha)

            self.plotter.ax.plot(x, pdf, color=self.color, linestyle=self.style, label=self.label)


    class KernelDensityEstimate2D(BasePlot):
        """Plots contours of a 2D Kernel Density Estimate (KDE) of pre-existing random samples"""

        _api_doc = inspect.cleandoc("""\
        Plotting 2D Kernel Density Estimates
        ------------------------------------

        Contents items of type ``kde2D`` are used to display contours of a two-dimensional
        kernel density estimate (a 2D smooth histogram) of samples of a probability density,
        be it a prior, a posterior, or a signal PDF.

        The following key is mandatory:

         * ``data`` (*dict*, see below) -- The data on the probability density that will be histogramed.

           Within the data object, the following keys are understood.

            * ``samples`` (*list* of *float* with shape (N, 2)) -- The samples that will be histogramed. Mandatory.
            * ``weights`` or ``log_weights`` (*list* of *float*, optional) -- The weights of the samples, on a linear or logarithmic scale.
              Defaults to uniform weights.


        The following keys are optional:

         * ``bandwidth`` (*float*) -- The factor by which the automatically determined kernel bandwidth is scaled. See the SciPy documentation
           for ``gaussian_kde``, ``bw_method='silverman'``. Defaults to 1.
         * ``contours`` (a *list* containing ``'lines'``, ``'areas'``, or both) -- The setting for the illustration of the contours.
           If ``'lines'`` is provided, the contour lines are drawn.
           If ``'areas'`` is provided, the contour areas are filled.
           Defaults to ``['lines']``.
         * ``levels`` (*list* of *float*) -- The probability levels of the contours. Defaults to ``[68, 95, 99]``.

        Example:

        .. code-block::

           analysis = ... # eos.Analysis object as discussed in the example notebook `inference.ipynb`
           parameter_samples, _, = analysis.sample(N=5000, pre_N=1000)
           plot_args = {
               'plot': {
                   'x': { 'label': r'$|V_{cb}|$', 'range': [38e-3, 47e-3] },
                   'y': { 'label': r'$f_+(0)$',   'range': [0.6, 0.75] },
               },
               'contents': [
                   {
                       'type': 'kde2D', 'color': 'C1', 'label': 'posterior',
                       'levels': [68, 95], 'contours': ['lines','areas'], 'bandwidth':3,
                       'data': { 'samples': parameter_samples[:, (0,1)] }
                   }
               ]
           }
           eos.plot.Plotter(plot_args).plot()

        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'data-file' not in item:
                raise KeyError('neither data nor data-file specified')

            if 'data' in item and 'data-file' in item:
                eos.warn('   both data and data-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples'].T

                if 'weights' in item['data'] and 'log_weights' in item['data']:
                    raise KeyError("Only one of 'weights' and 'log_weights' must be specified")
                elif 'weights' in item['data']:
                    self.weights = item['data']['weights']
                elif 'log_weights' in item['data']:
                    self.weights = np.exp(item['data']['log_weights'])
                else:
                    self.weights = None

            self.bw       = item['bandwidth'] if 'bandwidth' in item else None
            self.levels   = item['levels']    if 'levels'    in item else [0, 68, 95, 99]
            if 0 not in self.levels:
                self.levels = [0] + self.levels
            self.contours = item['contours']  if 'contours'  in item else ['lines']
            self.xrange   = plotter.xrange    if plotter.xrange      else (np.amin(self.samples[:, 0]), np.amax(self.samples[:, 0]))
            self.yrange   = plotter.yrange    if plotter.yrange      else (np.amin(self.samples[:, 1]), np.amax(self.samples[:, 1]))
            if type(self.style) is not list:
                self.style= [self.style]

        def plot(self):
            kde = gaussian_kde(self.samples, weights=self.weights)
            kde.set_bandwidth(bw_method='silverman')
            if self.bw:
                kde.set_bandwidth(bw_method=kde.factor * self.bw)

            xx, yy = np.mgrid[self.xrange[0]:self.xrange[1]:100j, self.yrange[0]:self.yrange[1]:100j]
            positions = np.vstack([xx.ravel(), yy.ravel()])
            pdf = np.reshape(kde(positions).T, xx.shape)
            pdf /= pdf.sum()

            # find the PDF value corresponding to a given cummulative probability
            plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
            plevels = []
            labels = []
            for level in self.levels:
                plevels.append(scipy.optimize.brentq(plevel, 0., 1., args=(pdf, level / 100.0)))
                labels.append(f'{level}%')

            if 'areas' in self.contours:
                colors = [matplotlib.colors.to_rgba(self.color, alpha) for alpha in np.linspace(0.50, 1.00, len(self.levels))]
                self.plotter.ax.contourf(pdf.transpose(),
                             colors=colors,
                             extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                             levels=plevels[::-1])

            CS = self.plotter.ax.contour(pdf.transpose(),
                             colors=self.color,
                             extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                             levels=plevels[::-1],
                             linestyles=self.style[::-1])

            if 'labels' in self.contours:
                fmt = {}
                for level, label in zip(CS.levels, labels[::-1]):
                    fmt[level] = label

                self.plotter.ax.clabel(CS, inline=1, fmt=fmt, fontsize=10)

        def handles_labels(self):
            if self.label:
                handle = None
                if 'areas' in self.contours:
                    handle = plt.Rectangle((0,0),1,1, color=self.color)
                else:
                    handle = plt.Line2D((0,1),(0.5,0.), color=self.color, linestyle=self.style[0])

                return ([handle], [self.label])
            else:
                return ([], [])


    class ExternalLikelihood2D(BasePlot):
        """Plots contours of a user-provided function"""

        _api_doc = inspect.cleandoc("""\
        Plotting a user-provided Likelihood
        -----------------------------------

        The following key is mandatory:

         * ``likelihood`` (*evaluable*, see below) -- The 2D probability density that will be plotted.

        The following keys are optional:

         * ``contours`` (a *list* containing ``'lines'``, ``'areas'``, or both) -- The setting for the illustration of the contours.
           If ``'lines'`` is provided, the contour lines are drawn.
           If ``'areas'`` is provided, the contour areas are filled.
           Defaults to ``['lines']``.
         * ``levels`` (*list* of *float*) -- The probability levels of the contours. Defaults to ``[68, 95, 99]``.
         * ``xrange`` (*list* of *float*) -- The x-axis range where the likelihood is evaluated, defaults to the plot ranges if provided.
         * ``yrange`` (*list* of *float*) -- The y-axis range where the likelihood is evaluated, defaults to the plot ranges if provided.

        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if not self.plotter.xrange and 'xrange' not in item:
                raise KeyError('xrange was not specified')
            if not self.plotter.yrange and 'yrange' not in item:
                raise KeyError('yrange was not specified')
            if 'likelihood' not in item:
                raise KeyError('likelihood was not specified')

            self.likelihood = item['likelihood']
            self.levels   = item['levels']    if 'levels'    in item else [0, 68, 95, 99]
            if 0 not in self.levels:
                self.levels = [0] + self.levels
            self.contours = item['contours']  if 'contours'  in item else ['lines']
            self.xrange   = item['xrange']    if 'xrange'    in item else self.plotter.xrange
            self.yrange   = item['yrange']    if 'yrange'    in item else self.plotter.yrange
            if type(self.style) is not list:
                self.style= [self.style]

        def plot(self):
            xx, yy = np.mgrid[self.xrange[0]:self.xrange[1]:100j, self.yrange[0]:self.yrange[1]:100j]
            positions = np.vstack([xx.ravel(), yy.ravel()])
            pdf = np.reshape([self.likelihood(point) for point in positions.T], xx.shape)
            pdf /= pdf.sum()

            # find the PDF value corresponding to a given cummulative probability
            plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
            plevels = []
            labels = []
            for level in self.levels:
                plevels.append(scipy.optimize.brentq(plevel, 0., 1., args=(pdf, level / 100.0)))
                labels.append(f'{level}%')

            if 'areas' in self.contours:
                colors = [matplotlib.colors.to_rgba(self.color, alpha) for alpha in np.linspace(0.50, 1.00, len(self.levels))]
                self.plotter.ax.contourf(pdf.transpose(),
                             colors=colors,
                             extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                             levels=plevels[::-1])

            CS = self.plotter.ax.contour(pdf.transpose(),
                             colors=self.color,
                             extent=[self.xrange[0], self.xrange[1], self.yrange[0], self.yrange[1]],
                             levels=plevels[::-1],
                             linestyles=self.style[::-1])

            if 'labels' in self.contours:
                fmt = {}
                for level, label in zip(CS.levels, labels[::-1]):
                    fmt[level] = label

                self.plotter.ax.clabel(CS, inline=1, fmt=fmt, fontsize=10)

        def handles_labels(self):
            if self.label:
                handle = None
                if 'areas' in self.contours:
                    handle = plt.Rectangle((0,0),1,1, color=self.color)
                else:
                    handle = plt.Line2D((0,1),(0.5,0.), color=self.color, linestyle=self.style[0])

                return ([handle], [self.label])
            else:
                return ([], [])


    class Histogram1D(BasePlot):
        """Plots a 1D histogram of pre-existing random samples"""

        _api_doc = inspect.cleandoc("""\
        Plotting Histograms
        -------------------

        Contents items of type ``histogram`` are used to display samples of a probability density, be it a prior, a posterior, or a signal PDF.
        The following key is mandatory:

         * ``data`` (*dict*, see below) -- The data on probability density that will be histogramed.

        Within the data object, the following keys are understood.

         * ``samples`` (*list* of *float*) -- The samples that will be histogramed. Mandatory.
         * ``weights`` or ``log_weights`` (*list* of *float*, optional) -- The weights of the samples, on a linear or logarithmic scale.
           Defaults to uniform weights.

        Example:

        .. code-block::

           analysis = ... # eos.Analysis object as discussed in the example notebook `inference.ipynb`
           parameter_samples, _, = analysis.sample(N=5000, pre_N=1000)
           plot_args = {
               'plot': {
                   'x': { 'label': r'$|V_{cb}|$' },
                   'y': { 'label': r'$d\\mathcal{B}/dq^2$' }
               },
               'contents': [
                   {
                       'type': 'histogram',
                       'data': { 'samples': parameter_samples[:, 0] },
                   },
               ]
           }

        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'data-file' not in item:
                raise KeyError('neither data nor hdf5-file specified')

            if 'data' in item and 'data-file' in item:
                eos.warn('   both data and data-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples']

                if 'weights' in item['data'] and 'log_weights' in item['data']:
                    raise KeyError("Only one of 'weights' and 'log_weights' must be specified")
                elif 'weights' in item['data']:
                    self.weights = item['data']['weights']
                elif 'log_weights' in item['data']:
                    self.weights = np.exp(item['data']['log_weights'])
                else:
                    self.weights = None
            else:
                if 'variable' not in item:
                    raise KeyError('no variable specificed')

                dfname = item['data-file']
                eos.info(f'   plotting histogram from "{dfname}"')
                prefix = os.path.split(dfname)[-1]
                eos.info(f'   prefix = {prefix}')
                if prefix.startswith('mcmc-'):
                    df  = eos.data.MarkovChain(dfname)
                    idx = df.lookup_table[item['variable']]
                    self.samples = df.samples[:, idx]
                    self.weights = None
                elif prefix.startswith('samples'):
                    df = eos.data.ImportanceSamples(dfname)
                    idx = df.lookup_table[item['variable']]
                    self.samples = df.samples[:, idx]
                    self.weights = df.weights
                elif prefix.startswith('pred-'):
                    df = eos.data.Prediction(dfname)
                    idx = df.lookup_table[item['variable']]
                    self.samples = df.samples[:, idx]
                    self.weights = df.weights
                else:
                    raise ValueError(f'Do not recognize data-file prefix: {dfname}')

            if 'bins' in item and type(item['bins']) is list:
                self.bins = item['bins']
            else:
                nbins = int(item['bins']) if 'bins' in item else 100
                xmin, xmax = self.xlo, self.xhi
                if xmin >= xmax:
                    eos.error(f'1D histogram expects the xrange to be non-zero; using empirical xrange instead')
                    xmin, xmax = np.min(self.samples), np.max(self.samples)
                self.bins = np.linspace(xmin, xmax, nbins)
            self.histtype = item['histtype']  if 'histtype'  in item     else 'bar'
            self.lw       = 3                 if self.histtype == 'step' else 1

        def plot(self):
            self.plotter.ax.hist(self.samples, weights=self.weights, density=True,
                    alpha=self.alpha, bins=self.bins, color=self.color,
                    histtype=self.histtype, label=self.label, lw=self.lw)


    class Histogram2D(BasePlot):
        """Plots a 2D histogram of pre-existing random samples"""

        _api_doc = inspect.cleandoc("""\
        Plotting 2D Histograms
        ----------------------

        Contents items of type ``histogram2D`` are used to display samples of a two-dimension probability density,
        be it a prior, a posterior, or a signal PDF.
        The following key is mandatory:

         * ``data`` (*dict*, see below) -- The data on probability density that will be histogramed.

        Within the data object, the following keys are understood.

         * ``samples`` (*list* of *float* with shape (N, 2)) -- The samples that will be histogramed. Mandatory.
         * ``weights`` or ``log_weights`` (*list* of *float* with length N, optional) -- The weights of the samples,
           on a linear or logarithmic scale, respectively. Defaults to uniform weights.

        Additional optional keys are:

         * ``bins`` (*int* or *list* of *int*, optional) -- The number of bins for both dimensions together  or in each dimension seperately. Defaults to 100 for both dimensions.

        Example:

        .. code-block::

           analysis = ... # eos.Analysis object as discussed in the example notebook `inference.ipynb`
           dstarlnu_kinematics = ... # create kineamtics and SignalPDF as discussed in the example notebook `simulation.ipynb`
           dstarlnu_pdf        = ...
           dstarlnu_samples, _ = dstarlnu_pdf.sample_mcmc(N=50000, stride=5, pre_N=1000, preruns=3, rng=rng)
           plot_args = {
               'plot': {
                   'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [ 0.0, 10.50] },
                   'y': { 'label': r'$cos(\theta_\\ell)$',                     'range': [-1.0,  +1.0] },
               },
               'contents': [
                   {
                       'label': r'samples ($\\ell=\\mu$)',
                       'type': 'histogram2D',
                       'data': {
                           'samples': dstarlnu_samples[:, (0, 1)]
                       },
                       'bins': 40
                   },
               ]
           }
           eos.plot.Plotter(plot_args).plot()

        """)

        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'data-file' not in item:
                raise KeyError('neither data not data-file specified')

            if 'data' in item and 'data-file' in item:
                eos.warn('   both data and hdf5-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples']

                if 'weights' in item['data'] and 'log_weights' in item['data']:
                    raise KeyError("Only one of 'weights' and 'log_weights' must be specified")
                elif 'weights' in item['data']:
                    self.weights = item['data']['weights']
                elif 'log_weights' in item['data']:
                    self.weights = np.exp(item['data']['log_weights'])
                else:
                    self.weights = None

            self.bins  = item['bins']    if 'bins'    in item else 100

        def plot(self):
            #cmap = plt.get_cmap('viridis')
            #cmap.set_under('w', 1)

            self.plotter.ax.hist2d(self.samples[:, 0], self.samples[:, 1], bins=self.bins, cmin=1, cmap=plt.get_cmap('viridis'),
                       label=self.label)


    class SignalPDF(BasePlot):
        """Plots a single EOS signal PDF w/o uncertainties as a function of one kinemtic variable"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            pname = item['pdf']
            eos.info(f'   plotting EOS PDF "{pname}"')

            # create parameters
            parameters = eos.Parameters.Defaults()
            if 'parameters' in item and 'parameters-from-file' in item:
                eos.warn('    overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            if 'parameters-from-file' in item and type(item['parameters-from-file']) is str:
                eos.warn('    overriding parameters from file')
                parameters.override_from_file(item['parameters-from-file'])

            if 'parameters' in item and type(item['parameters']) is dict:
                for key, value in item['parameters'].items():
                    parameters.set(key, value)

            # create kinematics
            kinematics = eos.Kinematics()
            if not ('kinematic' in item or 'variable' in item):
                raise KeyError('no kinematic variable specified; do not know how to map x to a variable')
            if 'kinematic' in item:
                var = kinematics.declare(item['kinematic'], np.nan)
            elif 'variable' in item:
                var = kinematics.declare(item['variable'], np.nan)

            if 'kinematics' in item:
                for k, v in item['kinematics'].items():
                    kinematics.declare(k, v)

            # create options
            options = eos.Options()
            if 'options' in item and type(item['options']) is dict:
                for key, value in item['options'].items():
                    options.declare(key, value)

            # create observable
            pdf  = eos.SignalPDF.make(pname, parameters, kinematics, options)
            norm = pdf.normalization()

            xvalues = np.linspace(self.xlo, self.xhi, self.xsamples + 1)
            pvalues = np.array([])
            for xvalue in xvalues:
                var.set(xvalue)
                pvalues = np.append(pvalues, pdf.evaluate())

            self.plotter.ax.plot(xvalues, np.exp(pvalues - norm), alpha=self.alpha, color=self.color, label=self.label, ls=self.style)


    class Expression(BasePlot):
        """Plots a given expression"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            if 'f' not in item:
                raise KeyError('no function specificied')
            f = item['f']
            alpha  = item['opacity'] if 'opacity' in item else 1.0
            color  = item['color']   if 'color'   in item else 'black'
            samples = item['samples'] if 'samples' in item else 100
            label  = item['label']   if 'label'   in item else None

            xmin, xmax = self.plotter.ax.get_xlim()
            x = np.linspace(xmin, xmax, samples)
            y = []

            for xvalue in x:
                y.append(eval(f, {}, {'x': xvalue}))

            self.plotter.ax.plot(x, y, color=color, alpha=alpha, linestyle=self.style, label=label)


    class Watermark(BasePlot):
        """Inserts an EOS watermark into the plots"""
        def __init__(self, plotter, item):
            super().__init__(plotter, item)
            self.z_order = sys.maxsize

        def plot(self):
            item = self.item
            xdelta, ydelta = (0.04, 0.04)

            hpos, vpos = item['position'] if 'position' in item else ['right', 'top']

            if hpos == 'right':
                x = 1 - xdelta
            elif hpos == 'left':
                x = xdelta
            elif hpos == 'center':
                x = 0.5
            else:
                raise ValueError(f'invalid horizontal position \'{hpos}\'')

            if vpos == 'bottom':
                y = 0 + ydelta
            elif vpos == 'top':
                y = 1 - ydelta
            elif vpos == 'center':
                y = 0.5
            else:
                raise ValueError(f'invalid vertical position \'{hpos}\'')

            ax = self.plotter.ax
            color = 'OrangeRed'
            version = f'v{eos.__version__}'
            if 'preliminary' in item and item['preliminary']:
                color = 'red'
                version = 'Preliminary'
            ax.text(x, y, fr'\textsf{{\textbf{{EOS {version}}}}}',
                    transform=ax.transAxes,
                    color=color, alpha=0.5, bbox=dict(facecolor='white', alpha=0.5, lw=0),
                    horizontalalignment=hpos, verticalalignment=vpos, zorder=+5)


    # This dict defines the user-facing keys for the plot items
    # after all classes are defined. Documentation can be generated
    # based on this class attribute.
    plot_types = {
        'band':                  Band,
        'constraint':            Constraint,
        'constraint2D':          Constraint2D,
        'constraint-overview':   ConstraintOverview,
        #'contours2D':            Contours2D,
        'graph':                 Graph,
        'expression':            Expression,
        'errorbar':              ErrorBar,
        'histogram':             Histogram1D,
        'histogram2D':           Histogram2D,
        'kde':                   KernelDensityEstimate1D,
        'kde2D':                 KernelDensityEstimate2D,
        'likelihood2D':          ExternalLikelihood2D,
        'observable':            Observable,
        'point':                 Point,
        'signal-pdf':            SignalPDF,
        'uncertainty':           Uncertainty,
        'uncertainty-binned':    UncertaintyBinned,
        #'uncertainty-overview':  UncertaintyOverview,
        'watermark':             Watermark,
    }


    def plot_contents(self):
        """Plots the contents specified in the instructions provided to Plotter"""
        if not 'contents' in self.instructions:
            return

        contents = self.instructions['contents']

        plots = []
        for item in contents:
            if not type(item) is dict:
                TypeError(f'wrong data type for content item {str(item)}')

            name = item['name'] if 'name' in item else 'unnamed'
            if not 'type' in item:
                raise KeyError(f'plot content "{name}" has no type')
            item_type = item['type']

            if 'name' not in item:
                name = None
                eos.debug(f'plotting anonymous contents of type \'{item_type}\'')
            else:
                name = item['name']
                eos.info(f'plotting "{name}"')

            if item_type not in self.plot_types:
                KeyError(f'unknown content type: "{item_type}"')

            plots.append(self.plot_types[item_type](self, item))

        # ensure watermarking
        if Plotter.Watermark not in [type(p) for p in plots]:
            plots.append(Plotter.Watermark(self, item))

        plots.sort()
        for plot in plots:
            plot.plot()

        if 'legend' in self.instructions['plot']:
            if 'location' in self.instructions['plot']['legend']:
                handles, labels = self.ax.get_legend_handles_labels()
                for plot in plots:
                    add_handles, add_labels = plot.handles_labels()
                    handles = handles + add_handles
                    labels  = labels  + add_labels
                loc  = self.instructions['plot']['legend']['location']
                ncol = self.instructions['plot']['legend']['columns'] if 'columns' in self.instructions['plot']['legend'] else 1
                self.ax.legend(handles, labels, loc=loc, ncol=ncol)


    def plot(self):
        """Produces the plot

        :returns:
            - fig (:py:class:`matplotlib.figure.Figure`) - the created figure
            - ax (:py:class:`matplotlib.axes.Axes`) - the created axes
        """
        self.setup_plot()
        self.plot_contents()

        plt.tight_layout()

        if self.output:
            plt.savefig(self.output, bbox_inches='tight', dpi=300)

        return self.fig, self.ax


def variable_to_latex(variable):
    p = eos.Parameters.Defaults()
    if variable in [pp.name() for pp in p]:
        return p[variable].latex()
    else:
        return r'\verb{' + variable + '}'
