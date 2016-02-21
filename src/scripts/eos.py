# Copyright (c) 2016 Danny van Dyk
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

import argparse
import h5py
import matplotlib
# use PDF renderer by default
matplotlib.use('AGG')
import matplotlib.pyplot as plot
import numpy
import os
import sys

from scipy.stats.kde import gaussian_kde


def error(message):
    print >> sys.stderr, '%s: error: %s' % (os.path.basename(sys.argv[0]), message)
    exit(1)

def warn(message):
    print >> sys.stderr, '%s: warning: %s' % (os.path.basename(sys.argv[0]), message)

class FileFormatError(Exception):
    def __init__(self, expected, found):
        self.expected = expected
        self.found = found

    def __str__(self):
        return 'Expected file format %s, found %s instead' % (self.expected, self.found)


class PMCDataFile:
    def __init__(self, file):
        # open the input file for reading
        self.file = h5py.File(file, 'r')

        # check that the input file has format=PMC
        if 'format' not in self.file.attrs:
            warn('input file does not have attribute \'format\'; assuming format \'PMC\'')
        elif 'PMC' != self.file.attrs['format']:
            raise FileFormatError('PMC', self.file.attrs['format'])

        # extract parameter descriptions of the n-tuples
        self.parameters = None
        if '/descriptions/parameters' in self.file:
            self.parameters = self.file['/descriptions/parameters']
        else:
            error('input file has no valid parameter descriptions: is it corrupted?')

    def __del__(self):
        self.file.close()

    """Retrieve data"""
    def data(self):
        step = 'final'

        if '/data/%s/samples' % step not in self.file:
            error('input file does not contain stored data for step %s' % step)

        dataset = self.file['/data/%s/samples' % step]

        return numpy.array(dataset[:])


class MCMCDataFile:
    def __init__(self, file):
        # open the input file for reading
        self.file = h5py.File(file, 'r')

        # check that the input file has format=MCMC
        if 'format' not in self.file.attrs:
            warn('input file does not have attribute \'format\'; assuming format \'MCMC\'')
        elif 'MCMC' != self.file.attrs['format']:
            raise FileFormatError('MCMC', self.file.attrs['format'])

        # extract parameter descriptions of the n-tuples
        self.parameters = None
        if '/descriptions/main run/chain #0/parameters' in self.file:
            self.parameters = self.file['/descriptions/main run/chain #0/parameters']
        elif '/descriptions/prerun/chain #0/parameters' in self.file:
            self.parameters = self.file['/descriptions/prerun/chain #0/parameters']
        else:
            error('input file has no valid parameter descriptions: is it corrupted?')

    def __del__(self):
        self.file.close()

    """Retrieve data"""
    def data(self):
        groupname = 'main run'

        if 'main run' not in self.file:
            warn('input file does not contain results from a main run')
            groupname = 'prerun'

        group = self.file[groupname]

        # start with no data
        data = None

        # append each dataset to data
        for chainname in group:
            chain = group[chainname]
            dset = chain['samples']

            if data == None:
                data = numpy.array(dset[:])
            else:
                data = numpy.append(data, dset[:], axis=0)

        return data


class Plotter1D:
    def __init__(self, datafile, pdffile):
        self.datafile = datafile
        self.pdffile = pdffile

    def histogram(self, index, **options):
        data = self.datafile.data()[:, index]
        # (name, min, max, nuisance, prior)
        parameter = self.datafile.parameters[index]

        plot.clf()
        plot.figure(figsize=(10, 10), dpi=80)

        # x axis
        plot.xlabel(parameter[0])
        xmin = options['xmin'] if options['xmin'] != None else parameter[1]
        xmax = options['xmax'] if options['xmax'] != None else parameter[2]
        plot.xlim(xmin, xmax)

        # y axis
        plot.ylabel('frequency')

        # plot
        plot.hist(data, bins=100, normed=1, alpha=.3)
        if 'kde' in options:
            kde = gaussian_kde(data)
            kde.set_bandwidth(bw_method='silverman')
            kde.set_bandwidth(bw_method=kde.factor / 3)
            x = numpy.linspace(xmin, xmax, 1000)
            plot.plot(x, kde(x), 'r')
        plot.tight_layout()

        # save figure
        plot.savefig(self.pdffile)


class Plotter2D:
    def __init__(self, datafile, pdffile):
        self.datafile = datafile
        self.pdffile = pdffile

    def histogram(self, index1, index2, **options):
        data = self.datafile.data()
        parameter1 = self.datafile.parameters[index1]
        parameter2 = self.datafile.parameters[index2]

        # compute histogram
        heatmap, xedges, yedges = numpy.histogram2d(data[:, index1], data[:, index2], bins=100)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        plot.clf()
        plot.figure(figsize=(10, 10), dpi=80)

        # x axis
        plot.xlabel(parameter1[0])
        xmin = options['xmin'] if options['xmin'] != None else parameter1[1]
        xmax = options['xmax'] if options['xmax'] != None else parameter1[2]
        plot.xlim(xmin, xmax)

        # y axis
        plot.ylabel(parameter2[0])
        ymin = options['ymin'] if options['ymin'] != None else parameter2[1]
        ymax = options['ymax'] if options['ymax'] != None else parameter2[2]
        plot.ylim(ymin, ymax)

        # plot
        plot.imshow(heatmap, extent=extent, aspect='auto')
        plot.tight_layout()

        plot.savefig(self.pdffile)


# set some default values for plotting
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.serif'] = 'Computer Modern Sans serif'
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['font.weight'] = 400

matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['axes.linewidth'] = 1

matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
