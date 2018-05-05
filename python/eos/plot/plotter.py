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

        if 'x' in myplot:
            myx = myplot['x']

            if 'label' in myx:
                myxlabel = myx['label']

            if 'unit' in myx:
                myxlabel += r'\,[' + myx['unit'] + r']'

            if 'range' in myx:
                self.xrange = myx['range']
                self.ax.set_xlim(tuple(self.xrange))

        if 'y' in myplot:
            myy = myplot['y']

            if 'label' in myy:
                myylabel = myy['label']

            if 'unit' in myy:
                myylabel += r'\,[' + myy['unit'] + r']'

            if 'range' in myy:
                self.yrange = myy['range']
                self.ax.set_ylim(tuple(self.yrange))

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


    def plot_eos_watermark(self, item):
        width, height = (0.115, 0.035)

        hpos, vpos = item['position'] if 'position' in item else ['right', 'bottom']

        if hpos == 'right':
            x = 0 + width
        elif hpos == 'left':
            x = 1 - width
        elif hpos == 'center':
            x = 0.5
        else:
            raise ValueError('invalid horizontal position \'{}\''.format(hpos))

        if vpos == 'bottom':
            y = 0 + height
        elif vpos == 'top':
            y = 1 - height
        elif vpos == 'center':
            y = 0.5
        else:
            raise ValueError('invalid vertical position \'{}\''.format(hpos))

        logofont = matplotlib.font_manager.FontProperties(family='sans-serif', size='15')
        ax = plt.gca()
        ax.text(x, y, r'\textsf{{\textbf{{EOS v{version}}}}}'.format(version=eos.version()),
                transform=ax.transAxes, fontproperties=logofont,
                color='OrangeRed', alpha=0.5, bbox=dict(facecolor='white', alpha=0.5, lw=0),
                horizontalalignment='center', verticalalignment='center', zorder=+5)


    def plot_contents(self):
        if not 'contents' in self.instructions:
            return

        plot_functions = {
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



class Plotter1D:
    def __init__(self, datafile, pdffile):
        self.datafile = datafile
        self.pdffile = pdffile


    def histogram(self, index, **options):
        data = self.datafile.data()[:, index]
        # (name, min, max, nuisance, prior)
        parameter = self.datafile.parameters[index]

        plt.clf()
        plt.figure(figsize=(10, 10), dpi=80)

        # x axis
        plt.xlabel("{} [{}]".format(str(parameter[0], 'utf-8'), str(parameter[1], 'utf-8')))
        xmin = float(options['xmin'] if options['xmin'] != None else parameter[2])
        xmax = float(options['xmax'] if options['xmax'] != None else parameter[3])
        plt.xlim(xmin, xmax)

        # y axis
        plt.ylabel('frequency')

        # plot
        plt.hist(data, bins=100, normed=1, alpha=.3)
        if options['kde']:
            kde = gaussian_kde(data)
            kde.set_bandwidth(bw_method='silverman')
            kde.set_bandwidth(bw_method=kde.factor * options['kde_bandwidth'])
            x = numpy.linspace(xmin, xmax, 1000)
            plot.plot(x, kde(x), 'r')

        plt.tight_layout()

        # save figure
        plt.savefig(self.pdffile)



class Plotter2D:
    def __init__(self, datafile, pdffile):
        self.datafile = datafile
        self.pdffile = pdffile


    def histogram(self, index1, index2, **options):
        data = self.datafile.data()

        # param 1
        parameter1 = self.datafile.parameters[index1]
        data1 = data[:, index1]

        parameter2 = self.datafile.parameters[index2]
        data2 = data[:, index2]

        plt.clf()
        plt.figure(figsize=(10, 10), dpi=80)

        # x axis
        plt.xlabel("{} [{}]".format(str(parameter1[0], 'utf-8'), str(parameter1[1], 'utf-8')))
        xmin = options['xmin'] if options['xmin'] != None else parameter1[2]
        xmax = options['xmax'] if options['xmax'] != None else parameter1[3]
        plt.xlim(xmin, xmax)

        # y axis
        plt.ylabel("{} [{}]".format(str(parameter2[0], 'utf-8'), str(parameter2[1], 'utf-8')))
        ymin = options['ymin'] if options['ymin'] != None else parameter2[2]
        ymax = options['ymax'] if options['ymax'] != None else parameter2[3]
        plt.ylim(ymin, ymax)

        # plot
        plt.hist2d(data1, data2, bins=100)
        plt.tight_layout()

        plt.savefig(self.pdffile)
