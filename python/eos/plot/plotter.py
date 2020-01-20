# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2017, 2018 Danny van Dyk
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
from logging import debug, info, warn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import gaussian_kde
import sys

class Plotter:
    """ Produces publication-quality plots

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


    """ Returns the next available z-order value, incremented for each plot w/o pre-defined z-order value. """
    @property
    def next_z_order(self):
        result = self.__next_z_order
        self.__next_z_order += 1
        return(result)

    """ Returns the next available color. """
    @property
    def next_color(self):
        result = self.colors[self.__next_color]
        self.__next_color = (self.__next_color + 1) % len(self.colors)
        return(result)


    """ Setting up the plot based on the provided instruction. """
    def setup_plot(self):
        if not 'plot' in self.instructions:
            raise KeyError('no plot metadata specified')

        myplot = self.instructions['plot']

        self.fig, self.ax = plt.subplots()

        mytitle = ''
        myylabel = ''
        myxlabel = ''
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

            if 'format' in myx:
                self.xformat = myx['format']
            else:
                self.xformat = '${x}$'

            if 'scale' in myx:
                self.xscale = float(myx['scale'])
                self.xticks = matplotlib.ticker.FuncFormatter(lambda x, pos, xscale=self.xscale: self.xformat.format(x=x / xscale))
            else:
                self.xticks = matplotlib.ticker.StrMethodFormatter(self.xformat)

            self.ax.xaxis.set_major_formatter(self.xticks)

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

            if 'format' in myy:
                self.yformat = myy['format']
            else:
                self.yformat = '${x}$'

            if 'scale' in myy:
                self.yscale = float(myy['scale'])
                self.yticks = matplotlib.ticker.FuncFormatter(lambda y, pos, yscale=self.yscale: self.yformat.format(x=y / yscale))
            else:
                self.yticks = matplotlib.ticker.StrMethodFormatter(self.yformat)

            self.ax.yaxis.set_major_formatter(self.yticks)

        if 'grid' in myplot:
            self.ax.grid(b=True, which=myplot['grid'])

        self.ax.set(xlabel=myxlabel, ylabel=myylabel, title=mytitle)

    """ Base class for any of the plots supported by Plotter. """
    class BasePlot:
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
            self.xlo, self.xhi = item['range']   if 'range'   in item else self.plotter.xrange

        def __lt__(self, other):
            return(self.z_order < other.z_order)

        def handles_labels(self):
            return ([], [])


    """ Plots a single point. """
    class Point(BasePlot):
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
            plt.plot(self.x, self.y,
                alpha=self.alpha, color=self.color,
                label=self.label, linestyle='None',
                marker=self.marker, markeredgecolor=self.color, markerfacecolor='None',
                markersize=self.markersize)


    """ Plots a single errors bar. """
    class ErrorBar(BasePlot):
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
            plt.errorbar(x=self.x, y=self.y, xerr=self.xerr, yerr=self.yerr,
                color=self.color, elinewidth=self.elinewidth, fmt='_', linestyle='none', label=self.label)


    """ Plots a shaded band. """
    class Band(BasePlot):
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
            handle = plt.Rectangle((0,0),1,1, color=self.color)
            label  = self.label
            return ([handle], [label])


    """ Plots a single EOS observable w/o uncertainties as a function of one kinemtic variable or one parameter. """
    class Observable(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            oname = item['observable']
            info('   plotting EOS observable "{}"'.format(oname))

            # create parameters
            parameters = eos.Parameters.Defaults()
            if 'parameters' in item and 'parameters-from-file' in item:
                warn('    overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            if 'parameters-from-file' in item and type(item['parameters-from-file']) is str:
                warn('    overriding parameters from file')
                parameters.override_from_file(item['parameters-from-file'])

            if 'parameters' in item and type(item['parameters']) is dict:
                for key, value in item['parameters'].items():
                    parameters.set(key, value)

            # create kinematics
            kinematics = eos.Kinematics()
            if not ('kinematic' in item or 'variable' in item) and not 'parameter' in item:
                raise KeyError('neither kinematic variable nor parameter found; do not know how to map x to a variable')
            if ('kinematic' in item or 'variable' in item) and 'parameter' in item:
                raise KeyError('both kinematic variable and parameter found; do not know how to map x to a variable')
            if 'kinematic' in item:
                var = kinematics.declare(item['kinematic'], np.nan)
            elif 'variable' in item:
                var = kinematics.declare(item['variable'], np.nan)
            else:
                var = parameters.declare(item['parameter'], np.nan)

            if 'kinematics' in item:
                for k, v in item['kinematics'].items():
                    kinematics.declare(k, v)

            # create (empty) options
            options = eos.Options()

            # create observable
            observable = eos.Observable.make(oname, parameters, kinematics, options)

            xvalues = np.linspace(self.xlo, self.xhi, self.xsamples + 1)
            ovalues = np.array([])
            for xvalue in xvalues:
                var.set(xvalue)
                ovalues = np.append(ovalues, observable.evaluate())

            plt.plot(xvalues, ovalues, alpha=self.alpha, color=self.color, label=self.label, ls=self.style)

    """ Plots an uncertainty band as a function of one kinematic variable.

    This routine expects the uncertainty propagation to have produces an HDF5 file.
    """
    class Uncertainty(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'hdf5-file' not in item:
                raise KeyError('neither data nor hdf5-file specified')

            if 'data' in item and 'hdf5-file' in item:
                warn('   both data and hdf5-file specified; assuming interactive use is intended')

            self.xvalues = None
            self.samples = None
            if 'data' in item:
                self.xvalues = np.array(item['data']['xvalues'])
                self.samples = item['data']['samples']
            else:
                h5fname = item['hdf5-file']
                info('   plotting uncertainty propagation from file "{}"'.format(h5fname))
                uncfile = eos.data.UncertaintyDataFile(h5fname)

                _xvalues = []
                for o in uncfile.parameters:
                    kin = o[1].split(b',')
                    if len(kin) > 1:
                        raise ValueError('more than one kinematic variable specified')

                    name, value = kin[0].split(b'=')
                    value = float(value)
                    _xvalues.append(value)

                self.xvalues = np.array(_xvalues)
                self.samples = uncfile.data()

            self.xrange = item['range'] if 'range' in item else None

        def plot(self):
            _ovalues_lower   = []
            _ovalues_central = []
            _ovalues_higher  = []
            for i in range(len(self.samples[0])):
                lower   = np.percentile(self.samples[:, i], q=15.865)
                central = np.percentile(self.samples[:, i], q=50.000)
                higher  = np.percentile(self.samples[:, i], q=84.135)
                _ovalues_lower.append(lower)
                _ovalues_central.append(central)
                _ovalues_higher.append(higher)

            xvalues = np.linspace(np.min(self.xvalues),np.max(self.xvalues),100)
            if self.xrange:
                xvalues = np.ma.masked_outside(xvalues, float(self.xrange[0]), float(self.xrange[1]))

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

            self.plotter.ax.fill_between(xvalues, ovalues_lower, ovalues_higher, alpha=self.alpha, color=self.color, label=self.label, lw=0)
            self.plotter.ax.plot(xvalues, ovalues_lower,                         alpha=self.alpha, color=self.color)
            self.plotter.ax.plot(xvalues, ovalues_central,                       alpha=self.alpha, color=self.color)
            self.plotter.ax.plot(xvalues, ovalues_higher,                        alpha=self.alpha, color=self.color)


    """ Plots one or more uncertainty band integrated over one kinematic variable.

    This routine expects the uncertainty propagation to have produces an HDF5 file.
    """
    class UncertaintyBinned(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'hdf5-file' not in item and 'data' not in item:
                raise KeyError('neither hdf5-file nor data specified')

            if 'data' in item:
                self.xvalues = item['data']['xvalues']
                self.samples = item['data']['samples']
            else:
                h5fname = item['hdf5-file']
                info('   plotting uncertainty propagation from file "{}"'.format(h5fname))

                uncfile = eos.data.UncertaintyDataFile(h5fname)

                if 'kinematic' not in item:
                    raise KeyError('kinematic not found; do not know how to map x to a kinematic variable')

                xname = item['kinematic']

                self.xvalues = []
                for o in uncfile.parameters:
                    kin = o[1].decode('ascii').split(',')
                    if len(kin) != 2:
                        raise ValueError('expected exactly two kinematic variables, got {}'.format(len(kin)))

                    name,value = kin[0].strip().split('=')
                    if name == xname + '_min':
                        xmin = float(value)
                    elif name == xname + '_max':
                        xmax = float(value)
                    else:
                        raise ValueError('unexpected kinematic variable \'{}\''.format(name))

                    name,value = kin[1].strip().split('=')
                    if name == xname + '_min':
                        xmin = float(value)
                    elif name == xname + '_max':
                        xmax = float(value)
                    else:
                        raise ValueError('unexpected kinematic variable \'{}\''.format(name))

                    self.xvalues.append([xmin, xmax])

                self.samples = uncfile.data()

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
                lower   = np.percentile(self.samples[:, i], q=15.865)
                central = np.percentile(self.samples[:, i], q=50.000)
                higher  = np.percentile(self.samples[:, i], q=84.135)
                ovalues_lower.append(lower)
                ovalues_central.append(central)
                ovalues_higher.append(higher)

            for [xmin, xmax], olo, ocentral, ohi in zip(self.xvalues, ovalues_lower, ovalues_central, ovalues_higher):
                width = (xmax - xmin) if self.rescale_by_width else 1
                olo      /= width
                ocentral /= width
                ohi      /= width
                print("{xmin} ... {xmax} -> {ocentral} with interval {olo} .. {ohi}".format(xmin=xmin, xmax=xmax, olo=olo, ocentral=ocentral, ohi=ohi))
                plt.fill_between([xmin, xmax], [olo, olo], [ohi, ohi], lw=0, color=self.color, alpha=self.alpha, label=self.label)
                label = None
                plt.plot([xmin, xmax], [olo,      olo],      color=self.color, alpha=self.alpha)
                plt.plot([xmin, xmax], [ocentral, ocentral], color=self.color, alpha=self.alpha)
                plt.plot([xmin, xmax], [ohi,      ohi],      color=self.color, alpha=self.alpha)
                self.label = None


    """ Plots an overview of uncertainty estimates.

    This routine expects the uncertainty propagation to have produced an HDF5 file.
    """
    class UncertaintyOverview(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'hdf5-file' not in item:
                raise KeyError('hdf5-file not specified')

            if 'observables' not in item:
                raise KeyError('observables not specified')

            h5fname = item['hdf5-file']
            info('   plotting uncertainty propagation from file "{}"'.format(h5fname))
            uncfile = eos.data.UncertaintyDataFile(h5fname)

            self.observables    = item['observables']
            self.observable_map = {}
            for idx, p in enumerate(uncfile.parameters):
                name = p[0].decode('ascii')
                self.observable_map[name] = idx

            for observable in self.observables:
                if observable not in self.observable_map:
                    raise ValueError('observable \'{}\' not contained in HDF5 file'.format(observable))

            self.samples = uncfile.data()

        def plot(self):
            for xvalue, observable in enumerate(self.observables):
                data_idx = self.observable_map[observable]
                lower    = np.percentile(self.samples[:, data_idx], q=15.865)
                upper    = np.percentile(self.samples[:, data_idx], q=84.135)

                self.plotter.ax.fill_between([xvalue - 0.5, xvalue + 0.5], [lower, lower], [upper, upper], alpha=self.alpha, color=self.color, label=self.label, lw=0)
                self.label = None


    """ Plots constraints from the EOS library of experimental and theoretical likelihoods. """
    class Constraint(BasePlot):
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

            if type(self.names) == str:
                self.names = [self.names]

        def plot(self):
            import yaml

            constraints = eos.Constraints()

            for name in self.names:
                entry = constraints[name]
                if not entry:
                    raise ValueError('unknown constraint {}'.format(name))

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
                    elif (self.variable + '_min' in kinematics) and (self.variable + '_max' in kinematics):
                        xvalues = [(kinematics[self.variable + '_max'] + kinematics[self.variable + '_min']) / 2]
                        xerrors = [(kinematics[self.variable + '_max'] - kinematics[self.variable + '_min']) / 2]
                        if self.rescale_by_width:
                            width = (kinematics[self.variable + '_max'] - kinematics[self.variable + '_min'])

                    yvalues = [float(constraint['mean']) / width]
                    sigma_hi = np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2) / width
                    sigma_lo = np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2) / width
                    yerrors = [[sigma_hi, sigma_lo]]
                elif constraint['type'] == 'MultivariateGaussian(Covariance)':
                    if not self.observable:
                        raise KeyError('observable needs to be specified for MultivariateGaussian(Covariance) constraints')
                    dim = constraint['dim']
                    covariance = np.array(constraint['covariance'])
                    observables = constraint['observables']
                    means = constraint['means']
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
                        elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                            xvalues.append((_kinematics[self.variable + '_max'] + _kinematics[self.variable + '_min']) / 2)
                            xerrors.append((_kinematics[self.variable + '_max'] - _kinematics[self.variable + '_min']) / 2)
                            if self.rescale_by_width:
                                width = (_kinematics[self.variable + '_max'] - _kinematics[self.variable + '_min'])

                        yvalues.append(means[i] / width)
                        yerrors.append(np.sqrt(covariance[i, i]) / width)
                elif constraint['type'] == 'MultivariateGaussian':
                    if not self.observable:
                        raise KeyError('observable needs to be specified for MultivariateGaussian constraints')
                    dim = constraint['dim']
                    sigma_stat_hi = np.array(constraint['sigma-stat-hi'])
                    sigma_stat_lo = np.array(constraint['sigma-stat-lo'])
                    sigma_sys = np.array(constraint['sigma-sys'])
                    sigma = np.sqrt(np.power(sigma_sys, 2) + 0.25 * np.power(sigma_stat_hi + sigma_stat_lo, 2))
                    observables = constraint['observables']
                    means = constraint['means']
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
                        elif (self.variable + '_min' in _kinematics) and (self.variable + '_max' in _kinematics):
                            xvalues.append((_kinematics[self.variable + '_max'] + _kinematics[self.variable + '_min']) / 2)
                            xerrors.append((_kinematics[self.variable + '_max'] - _kinematics[self.variable + '_min']) / 2)
                            if self.rescale_by_width:
                                width = (_kinematics[self.variable + '_max'] - _kinematics[self.variable + '_min'])

                        yvalues.append(means[i] / width)
                        yerrors.append(sigma[i] / width)
                else:
                    raise ValueError('type of constraint presently not supported')

                xvalues = np.array(xvalues)
                if xerrors:
                    xerrors = np.array(xerrors)
                yvalues = np.array(yvalues)
                yerrors = np.array(yerrors)

                plt.errorbar(x=xvalues, y=yvalues, xerr=xerrors, yerr=yerrors.T,
                    color=self.color, elinewidth=1.0, fmt='_', linestyle='none', label=self.label)
                # disable the label for subsequent plots
                self.label = None


    """ Plots overview of several constraints from the EOS library of experimental and theoretical likelihoods. """
    class ConstraintOverview(BasePlot):
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
                    raise ValueError('unknown constraint {}'.format(name))

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
                    dim         = constraint['dim']
                    observables = constraint['observables']
                    kinematics  = constraint['kinematics']
                    options     = constraint['options']
                    covariance  = np.array(constraint['covariance'])
                    means       = constraint['means']

                    for i in range(0, dim):
                        latex = '$' + eos.Observables()[observables[i]].latex() + '$'
                        sigma = np.sqrt(covariance[i, i])
                        xvalues.append(idx)
                        xticklabels.append(latex)
                        yvalues.append(means[i])
                        yerrors.append([sigma, sigma])
                        idx = idx + 1
                elif constraint['type'] == 'MultivariateGaussian':
                    dim = constraint['dim']
                    observables = constraint['observables']
                    sigma_stat_hi = np.array(constraint['sigma-stat-hi'])
                    sigma_stat_lo = np.array(constraint['sigma-stat-lo'])
                    sigma_sys = np.array(constraint['sigma-sys'])
                    sigma = np.sqrt(np.power(sigma_sys, 2) + 0.25 * np.power(sigma_stat_hi + sigma_stat_lo, 2))
                    means = constraint['means']

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
            plt.xticks(xvalues, xticklabels, rotation=self.rotation)
            plt.errorbar(x=xvalues, y=yvalues, xerr=None, yerr=yerrors.T,
                color=self.color, elinewidth=1.0, fmt='_', linestyle='none', label=self.label)
            plt.margins(0.2)
            # Tweak spacing to prevent clipping of tick-labels
            plt.subplots_adjust(bottom=0.15)


    """ Plots 2D contours of a pair of parameters based on pre-existing random samples. """
    class Contours2D(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            if 'hdf5-file' not in item:
                raise KeyError('no hdf5-file specified')

            h5fname = item['hdf5-file']
            info('   plotting 2D contours from file "{}"'.format(h5fname))
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

            CS = plt.contour(pdf.transpose(),
                             colors='OrangeRed',
                             extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                             levels=levels[::-1])

            fmt = {}
            for level, label in zip(CS.levels, labels[::-1]):
                fmt[level] = label

            plt.clabel(CS, inline=1, fmt=fmt, fontsize=10)


    """ Plots a 1D Kernel Density Estimate (KDE) of pre-existing random samples. """
    class KernelDensityEstimate1D(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'hdf5-file' not in item:
                raise KeyError('neither data nor hdf5-file specified')

            if 'data' in item and 'hdf5-file' in item:
                warn('   both data and hdf5-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples']
                self.weights = np.exp(item['data']['log_weights']) if 'log_weights' in item['data'] else None
            else:
                h5fname = item['hdf5-file']
                info('   plotting histogram from file "{}"'.format(h5fname))
                datafile = eos.data.load_data_file(h5fname)

                if 'variable' not in item:
                    raise KeyError('no variable specificed')
                variable = item['variable']

                if variable not in datafile.variable_indices:
                    raise ValueError('variable {} not contained in data file'.format(variable))

                index        = datafile.variable_indices[variable]
                self.samples = datafile.data()[:, index]
                # TODO: use weights from data file
                self.weights = None

            self.bw       = item['bandwidth'] if 'bandwidth' in item else None
            self.level    = item['level']     if 'level'     in item else 68.29
            self.xrange   = plotter.xrange    if plotter.xrange      else (np.amin(self.samples), np.amax(self.samples))

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
                self.plotter.ax.fill_between(x[pdf >= plevel * pdf_norm], pdf[pdf >= plevel * pdf_norm], facecolor=self.color, alpha=self.alpha)

            plt.plot(x, pdf, color=self.color, label=self.label)


    """ Plots contours of a 2D Kernel Density Estimate (KDE) of pre-existing random samples. """
    class KernelDensityEstimate2D(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'hdf5-file' not in item:
                raise KeyError('neither data nor hdf5-file specified')

            if 'data' in item and 'hdf5-file' in item:
                warn('   both data and hdf5-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples'].T
                self.weights = np.exp(item['data']['log_weights']) if 'log_weights' in item['data'] else None
            else:
                h5fname = item['hdf5-file']
                info('   plotting 2D KDE from file "{}"'.format(h5fname))
                datafile = eos.data.load_data_file(h5fname)

                if 'variables' not in item:
                    raise KeyError('no variables specificed')

                xvariable, yvariable = item['variables']

                if xvariable not in datafile.variable_indices:
                    raise ValueError('x-axis variable {} not contained in data file'.format(variable))

                if yvariable not in datafile.variable_indices:
                    raise ValueError('y-axis variable {} not contained in data file'.format(variable))

                stride = item['stride'] if 'stride' in item else 1
                data   = datafile.data()
                xindex = datafile.variable_indices[xvariable]
                xdata  = data[::stride, xindex]
                yindex = datafile.variable_indices[yvariable]
                ydata  = data[::stride, yindex]

                self.samples = np.vstack([xdata, ydata])
                # TODO: use weights from data file
                self.weights = None

            self.bw       = item['bandwidth'] if 'bandwidth' in item else None
            self.levels   = item['levels']    if 'levels'    in item else [0, 68, 95, 99]
            if 0 not in self.levels:
                self.levels = [0] + self.levels
            self.contours = item['contours']  if 'contours'  in item else ['lines']
            self.xrange   = plotter.xrange    if plotter.xrange      else (np.amin(self.samples[:, 0]), np.amax(self.samples[:, 0]))
            self.yrange   = plotter.xrange    if plotter.yrange      else (np.amin(self.samples[:, 1]), np.amax(self.samples[:, 1]))
            if type(self.style) is not list:
                self.style= [self.style]

        def plot(self):
            kde = gaussian_kde(self.samples, weights=self.weights)
            kde.set_bandwidth(bw_method='silverman')
            if self.bw:
                kde.set_bandwidth(bw_method=kde.factor * self.bw)

            xx, yy = np.mgrid[self.plotter.xrange[0]:self.plotter.xrange[1]:100j, self.plotter.yrange[0]:self.plotter.yrange[1]:100j]
            positions = np.vstack([xx.ravel(), yy.ravel()])
            pdf = np.reshape(kde(positions).T, xx.shape)
            pdf /= pdf.sum()

            # find the PDF value corresponding to a given cummulative probability
            plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
            plevels = []
            labels = []
            for level in self.levels:
                plevels.append(scipy.optimize.brentq(plevel, 0., 1., args=(pdf, level / 100.0)))
                labels.append('{}%'.format(level))

            if 'areas' in self.contours:
                colors = [matplotlib.colors.to_rgba(self.color, alpha) for alpha in np.linspace(0.50, 1.00, len(self.levels))]
                plt.contourf(pdf.transpose(),
                             colors=colors,
                             extent=[self.plotter.xrange[0], self.plotter.xrange[1], self.plotter.yrange[0], self.plotter.yrange[1]],
                             levels=plevels[::-1])

            CS = plt.contour(pdf.transpose(),
                             colors=self.color,
                             extent=[self.plotter.xrange[0], self.plotter.xrange[1], self.plotter.yrange[0], self.plotter.yrange[1]],
                             levels=plevels[::-1],
                             linestyles=self.style[::-1])

            if 'labels' in self.contours:
                fmt = {}
                for level, label in zip(CS.levels, labels[::-1]):
                    fmt[level] = label

                plt.clabel(CS, inline=1, fmt=fmt, fontsize=10)

        def handles_labels(self):
            handle = None
            if 'areas' in self.contours:
                handle = plt.Rectangle((0,0),1,1, color=self.color)
            else:
                handle = plt.Line2D((0,1),(0.5,0.), color=self.color, linestyle=self.style[0])

            return ([handle], [self.label])


    """ Plots a 1D histogram of pre-existing random samples. """
    class Histogram1D(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

            if 'data' not in item and 'hdf5-file' not in item:
                raise KeyError('neither data nor hdf5-file specified')

            if 'data' in item and 'hdf5-file' in item:
                warn('   both data and hdf5-file specified; assuming interactive use is intended')

            self.samples = None
            self.weights = None
            if 'data' in item:
                self.samples = item['data']['samples']
                self.weights = np.exp(item['data']['log_weights']) if 'log_weights' in item['data'] else None
            else:
                h5fname = item['hdf5-file']
                info('   plotting histogram from file "{}"'.format(h5fname))
                datafile = eos.data.load_data_file(h5fname)

                if 'variable' not in item:
                    raise KeyError('no variable specificed')
                variable = item['variable']

                if variable not in datafile.variable_indices:
                    raise ValueError('variable {} not contained in data file'.format(variable))

                index        = datafile.variable_indices[variable]
                self.samples = datafile.data()[:, index]
                # TODO: use weights from data file
                self.weights = None

            self.bins     = item['bins']      if 'bins'      in item     else 100
            self.histtype = item['histtype']  if 'histtype'  in item     else 'bar'
            self.lw       = 3                 if self.histtype == 'step' else 1

        def plot(self):
            plt.hist(self.samples, weights=self.weights, density=True,
                    alpha=self.alpha, bins=self.bins, color=self.color,
                    histtype=self.histtype, label=self.label, lw=self.lw)


    """ Plots a 2D histogram of pre-existing random samples. """
    class Histogram2D(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            if 'hdf5-file' not in item:
                raise KeyError('no hdf5-file specified')

            h5fname = item['hdf5-file']
            info('   plotting 2D histogram from file "{}"'.format(h5fname))
            datafile = eos.data.load_data_file(h5fname)

            if 'variables' not in item:
                raise KeyError('no variables specificed')

            xvariable, yvariable = item['variables']

            xindex = datafile.variable_indices[xvariable]
            yindex = datafile.variable_indices[yvariable]
            data = datafile.data()

            cmap = plt.get_cmap('viridis')
            cmap.set_under('w', 1)

            xdata = data[:, xindex]
            ydata = data[:, yindex]
            bins  = item['bins']    if 'bins'    in item else 100
            plt.hist2d(xdata, ydata, bins=bins, cmin=1)


    """ Plots a given expression. """
    class Expression(BasePlot):
        def __init__(self, plotter, item):
            super().__init__(plotter, item)

        def plot(self):
            item = self.item
            if 'f' not in item:
                raise KeyError('no function specificied')
            f = item['f']
            alpha  = item['opacity'] if 'opacity' in item else 1.0
            color  = item['color']   if 'color'   in item else 'black'
            style  = item['style']   if 'style'   in item else '-'
            samples = item['samples'] if 'samples' in item else 100
            label  = item['label']   if 'label'   in item else None

            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, samples)
            y = []

            for xvalue in x:
                y.append(eval(f, {}, {'x': xvalue}))

            plt.plot(x, y, color=color, alpha=alpha, linestyle=style, label=label)


    """ Inserts an EOS watermark into the plots. """
    class Watermark(BasePlot):
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
                raise ValueError('invalid horizontal position \'{}\''.format(hpos))

            if vpos == 'bottom':
                y = 0 + ydelta
            elif vpos == 'top':
                y = 1 - ydelta
            elif vpos == 'center':
                y = 0.5
            else:
                raise ValueError('invalid vertical position \'{}\''.format(hpos))

            ax = plt.gca()
            color = 'OrangeRed'
            version = 'v{version}'.format(version=eos.version())
            if 'preliminary' in item and item['preliminary']:
                color = 'red'
                version = 'Preliminary'
            ax.text(x, y, r'\textsf{{\textbf{{EOS {version}}}}}'.format(version=version),
                    transform=ax.transAxes,
                    color=color, alpha=0.5, bbox=dict(facecolor='white', alpha=0.5, lw=0),
                    horizontalalignment=hpos, verticalalignment=vpos, zorder=+5)


    """ Plots the contents specified in the instructions provided to Plotter. """
    def plot_contents(self):
        if not 'contents' in self.instructions:
            return

        plot_types = {
            'band':                  Plotter.Band,
            'constraint':            Plotter.Constraint,
            'constraint-overview':   Plotter.ConstraintOverview,
            'contours2D':            Plotter.Contours2D,
            'expression':            Plotter.Expression,
            'errorbar':              Plotter.ErrorBar,
            'histogram':             Plotter.Histogram1D,
            'histogram2D':           Plotter.Histogram2D,
            'kde':                   Plotter.KernelDensityEstimate1D,
            'kde2D':                 Plotter.KernelDensityEstimate2D,
            'observable':            Plotter.Observable,
            'point':                 Plotter.Point,
            'uncertainty':           Plotter.Uncertainty,
            'uncertainty-binned':    Plotter.UncertaintyBinned,
            'uncertainty-overview':  Plotter.UncertaintyOverview,
            'watermark':             Plotter.Watermark,
        }

        contents = self.instructions['contents']

        plots = []
        for item in contents:
            if not type(item) is dict:
                TypeError('wrong data type for content item {}'.format(str(item)))

            name = item['name'] if 'name' in item else 'unnamed'
            if not 'type' in item:
                raise KeyError('plot content "{}" has no type'.format(name))
            item_type = item['type']

            if 'name' not in item:
                name = None
                debug('plotting anonymous contents of type \'{}\''.format(item_type))
            else:
                name = item['name']
                info('plotting "{}"'.format(name))

            if item_type not in plot_types:
                KeyError('unknown content type: "{}"'.format(item_type))

            plots.append(plot_types[item_type](self, item))

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


    """ Produces the plot. """
    def plot(self):
        self.setup_plot()
        self.plot_contents()

        plt.tight_layout()

        if self.output:
            plt.savefig(self.output, bbox_inches='tight', dpi=300)

        plt.show()


def variable_to_latex(variable):
    p = eos.Parameters.Defaults()
    if variable in [pp.name() for pp in p]:
        return p[variable].latex()
    else:
        return r'\verb{' + variable + '}'


