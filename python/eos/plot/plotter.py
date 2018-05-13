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
from scipy.interpolate import spline
import sys

class Plotter:
    def __init__(self, instructions, output):
        self.instructions = instructions
        self.output = output
        self.fig = None
        self.ax = None
        self.xrange = None
        self.yrange = None


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
            # convert from cm to inches
            xwidth /= 2.54 # cm / inch
            ywidth /= 2.54 # cm / inch
            plt.gcf().set_size_inches((xwidth, ywidth))

        plt.locator_params(axis='x', nbins=5)
        plt.locator_params(axis='y', nbins=5)

        if 'x' in myplot:
            myx = myplot['x']

            if 'label' in myx:
                myxlabel = myx['label']

            if 'unit' in myx:
                myxlabel += r'\,[' + myx['unit'] + r']'

            if 'range' in myx:
                self.xrange = myx['range']
                self.ax.set_xlim(tuple(self.xrange))

            if 'scale' in myx:
                self.xscale = float(myx['scale'])
                self.xticks = matplotlib.ticker.FuncFormatter(lambda x, pos, xscale=self.xscale: '${0:.2f}$'.format(x / xscale))
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

            if 'scale' in myy:
                self.yscale = float(myy['scale'])
                self.yticks = matplotlib.ticker.FuncFormatter(lambda y, pos, yscale=self.yscale: '${0:.2f}$'.format(y / yscale))
                self.ax.yaxis.set_major_formatter(self.yticks)

        self.ax.set(xlabel=myxlabel, ylabel=myylabel, title=mytitle)


    def plot_observable(self, item):
        oname = item['observable']
        info('   plotting EOS observable "{}"'.format(oname))

        # create parameters
        parameters = eos.Parameters.Defaults()
        if 'parameters' in item and 'parameters-from-file' in item:
            warn('    overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

        if 'parameters-from-file' in item and type(item['parameters-from-file']) is str:
            parameters.override_from_file(item['parameters-from-file'])

        if 'parameters' in item and type(item['parameters']) is dict:
            for key, value in item['parameters'].items():
                parameters.set(key, value)

        # create kinematics
        kinematics = eos.Kinematics()
        if not 'kinematic' in item:
            raise KeyError('kinematic not found; do not know how to map x to a kinematic variable')
        kvar = kinematics.declare(item['kinematic'], np.nan)

        # create (empty) options
        options = eos.Options()

        # create observable
        observable = eos.Observable.make(oname, parameters, kinematics, options)

        xlo, xhi = self.xrange
        if 'range' in item:
            xlo, xhi = item['range']

        samples = 100
        if 'samples' in item:
            samples = item['samples']

        xvalues = np.linspace(xlo, xhi, samples + 1)
        ovalues = np.array([])
        for xvalue in xvalues:
            kvar.set(xvalue)
            ovalues = np.append(ovalues, observable.evaluate())

        color = 'black'
        if 'color' in item:
            color = item['color']

        plt.plot(xvalues, ovalues, color=color)


    def plot_uncertainty(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')
            return

        h5fname = item['hdf5-file']
        info('   plotting uncertainty propagation from file "{}"'.format(h5fname))

        uncfile = eos.data.UncertaintyDataFile(h5fname)
        _xvalues = []
        for o in uncfile.parameters:
            kin = o[1].split(b',')
            if len(kin) > 1:
                raise ValueError('more than one kinematic variable specified')

            name,value = kin[0].split(b'=')
            _xvalues.append(float(value))

        data = uncfile.data()
        _ovalues_lower   = []
        _ovalues_central = []
        _ovalues_higher  = []
        for i in range(len(uncfile.parameters)):
            lower   = np.percentile(data[:, i], q=0.15865)
            central = np.percentile(data[:, i], q=0.5)
            higher  = np.percentile(data[:, i], q=0.84135)
            _ovalues_lower.append(lower)
            _ovalues_central.append(central)
            _ovalues_higher.append(higher)

        color = 'black'
        if 'color' in item:
            color = item['color']

        alpha = 1.0
        if 'opacity' in item:
            alpha = item['opacity']

        # TODO: replace scipy.interpolate.spline
        xvalues = np.linspace(np.min(_xvalues),np.max(_xvalues),100)
        ovalues_lower   = scipy.interpolate.spline(_xvalues, _ovalues_lower,   xvalues)
        ovalues_central = scipy.interpolate.spline(_xvalues, _ovalues_central, xvalues)
        ovalues_higher  = scipy.interpolate.spline(_xvalues, _ovalues_higher,  xvalues)

        plt.fill_between(xvalues, ovalues_lower, ovalues_higher, lw=0, color=color, alpha=alpha)
        plt.plot(xvalues, ovalues_lower,   color=color, alpha=alpha)
        plt.plot(xvalues, ovalues_central, color=color, alpha=alpha)
        plt.plot(xvalues, ovalues_higher,  color=color, alpha=alpha)


    def plot_contours2d(self, item):
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

        if not np.array(self.xrange).any():
            self.xrange = [np.amin(xdata), np.amax(xdata)]
            self.ax.set_xlim(tuple(self.xrange))
        if not np.array(self.yrange).any():
            self.yrange = [np.amin(ydata), np.amax(ydata)]
            self.ax.set_ylim(tuple(self.yrange))
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


    def plot_kde(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')

        h5fname = item['hdf5-file']
        info('   plotting histogram from file "{}"'.format(h5fname))
        datafile = eos.data.load_data_file(h5fname)

        if 'variable' not in item:
            raise KeyError('no variable specificed')
        variable = item['variable']

        print(datafile.variable_indices)
        if variable not in datafile.variable_indices:
            raise ValueError('variable {} not contained in data file'.format(variable))

        index = datafile.variable_indices[variable]
        data  = datafile.data()[:, index]
        alpha = item['opacity']   if 'opacity'   in item else 0.3
        color = item['color']     if 'color'     in item else 'blue'
        bw    = item['bandwidth'] if 'bandwidth' in item else None

        kde = gaussian_kde(data)
        kde.set_bandwidth(bw_method='silverman')
        if 'bandwidth' in item:
            kde.set_bandwidth(bw_method=kde.factor * item['bandwidth'])

        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 1000)

        plt.plot(x, kde(x), color=color)


    def plot_histogram(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')

        h5fname = item['hdf5-file']
        info('   plotting histogram from file "{}"'.format(h5fname))
        datafile = eos.data.load_data_file(h5fname)

        if 'variable' not in item:
            raise KeyError('no variable specificed')
        variable = item['variable']

        print(datafile.variable_indices)
        if variable not in datafile.variable_indices:
            raise ValueError('variable {} not contained in data file'.format(variable))

        index = datafile.variable_indices[variable]
        data  = datafile.data()[:, index]
        alpha = item['opacity'] if 'opacity' in item else 0.3
        color = item['color']   if 'color'   in item else 'blue'
        bins  = item['bins']    if 'bins'    in item else 100
        plt.hist(data, alpha=alpha, bins=bins, color=color, density=1)


    def plot_histogram2d(self, item):
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


    def plot_eos_watermark(self, item):
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

        logofont = matplotlib.font_manager.FontProperties(family='sans-serif', size='15')
        ax = plt.gca()
        ax.text(x, y, r'\textsf{{\textbf{{EOS v{version}}}}}'.format(version=eos.version()),
                transform=ax.transAxes, fontproperties=logofont,
                color='OrangeRed', alpha=0.5, bbox=dict(facecolor='white', alpha=0.5, lw=0),
                horizontalalignment=hpos, verticalalignment=vpos, zorder=+5)


    def plot_contents(self):
        if not 'contents' in self.instructions:
            return

        plot_functions = {
            'contours2D':  Plotter.plot_contours2d,
            'histogram':   Plotter.plot_histogram,
            'histogram2D': Plotter.plot_histogram2d,
            'kde':         Plotter.plot_kde,
            'observable':  Plotter.plot_observable,
            'uncertainty': Plotter.plot_uncertainty,
            'watermark':   Plotter.plot_eos_watermark,
        }

        anonymous_types = [
            'watermark',
        ]

        contents = self.instructions['contents']

        for item in contents:
            if not type(item) is dict:
                TypeError('wrong data type for content item {}'.format(str(item)))

            if not 'type' in item:
                raise KeyError('plot content "{}" has no type'.format(name))
            item_type = item['type']

            if item_type not in anonymous_types and 'name' not in item:
                raise KeyError('unnamed plot content')
            elif 'name' not in item:
                name = None
                debug('plotting anonymous contents of type \'{}\''.format(item_type))
            else:
                name = item['name']
                info('plotting "{}"'.format(name))

            if item_type not in plot_functions:
                KeyError('unknown content type: "{}"'.format(item_type))

            plot_functions[item_type](self, item)


    def plot(self):
        self.setup_plot()
        self.plot_contents()

        plt.show()
        plt.savefig(self.output)


def variable_to_latex(variable):
    p = eos.Parameters.Defaults()
    if variable in [pp.name() for pp in p]:
        return p[variable].latex()
    else:
        return r'\verb{' + variable + '}'


