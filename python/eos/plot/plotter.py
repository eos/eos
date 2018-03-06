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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

import scipy
from scipy.interpolate import spline

# TODO monkeypathc warn() etc. cleaner
import os
def error(message):
    print('%s: error: %s' % (os.path.basename(sys.argv[0]), message), file=sys.stderr)
    exit(1)

def warn(message):
    print('%s: warning: %s' % (os.path.basename(sys.argv[0]), message), file=sys.stderr)

def info(message):
    print('%s: info: %s' % (os.path.basename(sys.argv[0]), message), file=sys.stderr)

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
            error('no plot metadata specified')

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
        if 'parameters' in item and type(item['parameters']) is dict:
            for key, value in item['parameters'].items():
                parameters.set(key, value)

        # create kinematics
        kinematics = eos.Kinematics()
        if not 'kinematic' in item:
            error('kinematic not found; do not know how to map x to a kinematic variable')
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
            error('no hdf5-file specified')
            return

        h5fname = item['hdf5-file']
        info('   plotting uncertainty propagation from file "{}"'.format(h5fname))

        uncfile = eos.data.UncertaintyDataFile(h5fname)
        _xvalues = []
        for o in uncfile.parameters:
            kin = o[1].split(b',')
            if len(kin) > 1:
                error('more than one kinematic variable specified')
                return
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


    def plot_contents(self):
        if not 'contents' in self.instructions:
            return

        contents = self.instructions['contents']
        for item in contents:
            if not type(item) is dict:
                error('wrong data type for content item {}'.format(str(item)))

            if not 'name' in item:
                error('unnamed plot content')
            name = item['name']

            if not 'type' in item:
                error('plot content "{}" has no type'.format(name))
            item_type = item['type']

            info('plotting "{}"'.format(name))

            plot_functions = {
                'observable': Plotter.plot_observable,
                'uncertainty': Plotter.plot_uncertainty,
            }

            if not item_type in plot_functions:
                error('unknown content type: "{}"'.format(item_type))

            plot_functions[item_type](self, item)


    def plot(self):
        self.setup_plot()
        self.plot_contents()

        plt.show()
        plt.savefig(self.output)
