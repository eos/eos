# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2016, 2017, 2018 Danny van Dyk
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
import h5py
import numpy
import os
import sys

class FileFormatError(Exception):
    def __init__(self, expected, found):
        self.expected = expected
        self.found = found

    def __str__(self):
        return 'Expected file format %s, found %s instead' % (self.expected, self.found)


class HDF5DataFile:
    def __init__(self):
        # generate map
        self.variable_indices = { self.parameters[i][0].decode('ascii'): i for i in range(0, len(self.parameters)) }


class PMCDataFile(HDF5DataFile):
    def __init__(self, file):
        # open the input file for reading
        self.file = h5py.File(file, 'r')

        # check that the input file has format=PMC
        if 'format' not in self.file.attrs:
            eos.warn('input file does not have attribute \'format\'; assuming format \'PMC\'')
        elif 'PMC' != self.file.attrs['format']:
            raise FileFormatError('PMC', self.file.attrs['format'])

        # extract parameter descriptions of the n-tuples
        self.parameters = None
        if '/descriptions/parameters' in self.file:
            self.parameters = self.file['/descriptions/parameters']
        else:
            RuntimeError('input file has no valid parameter descriptions: is it corrupted?')

        super().__init__()


    def __del__(self):
        self.file.close()


    """Retrieve data"""
    def data(self):
        step = 'final'

        if '/data/%s/samples' % step not in self.file:
            RuntimeError('input file does not contain stored data for step %s' % step)

        dataset = self.file['/data/%s/samples' % step]

        return numpy.array(dataset[:])


class MCMCDataFile(HDF5DataFile):
    def __init__(self, file):
        # open the input file for reading
        self.file = h5py.File(file, 'r')

        # check that the input file has format=MCMC
        if 'format' not in self.file.attrs:
            eos.warn('input file does not have attribute \'format\'; assuming format \'MCMC\'')
        elif 'MCMC' != self.file.attrs['format']:
            raise FileFormatError('MCMC', self.file.attrs['format'])

        # extract parameter descriptions of the n-tuples
        self.parameters = None
        if '/descriptions/main run/chain #0/parameters' in self.file:
            self.parameters = self.file['/descriptions/main run/chain #0/parameters']
        elif '/descriptions/prerun/chain #0/parameters' in self.file:
            self.parameters = self.file['/descriptions/prerun/chain #0/parameters']
        else:
            RuntimeError('input file has no valid parameter descriptions: is it corrupted?')

        super().__init__()


    def __del__(self):
        self.file.close()


    """Retrieve data"""
    def data(self):
        groupname = 'main run'

        if 'main run' not in self.file:
            eos.warn('input file does not contain results from a main run')
            groupname = 'prerun'

        group = self.file[groupname]

        # start with no data
        data = None

        # append each dataset to data
        for chainname in group:
            chain = group[chainname]
            dset = chain['samples']

            if data is None:
                data = numpy.array(dset[:])
            else:
                data = numpy.append(data, dset[:], axis=0)

        return data

    """Retrieve the modes of the chains"""
    def modes(self):
        groupname = 'main run'

        if 'main run' not in self.file:
            eos.warn('input file does not contain results from a main run')
            groupname = 'prerun'

        group = self.file[groupname]

        # start with no data
        result = []

        # append each dataset to data
        for chainname in group:
            chain = group[chainname]
            dset = chain['stats/mode']

            log_posterior = dset[-1][-1]
            mode          = dset[-1][0:-1]

            result.append((mode, log_posterior))

        return result

class UncertaintyDataFile(HDF5DataFile):
    def __init__(self, file):
        self.name = file
        # open the input file for reading
        self.file = h5py.File(file, 'r')

        # check that the input file has format=PMC
        if 'format' not in self.file.attrs:
            eos.warn('input file does not have attribute \'format\'; assuming format \'UNC\'')
        elif 'UNC' != self.file.attrs['format']:
            raise FileFormatError('UNC', self.file.attrs['format'])

        # extract parameter descriptions of the n-tuples
        self.parameters = []
        if '/descriptions/parameters' in self.file:
            for i in range(len(self.file['/descriptions/observables'])):
                desc = self.file['/descriptions/observables/%d' % i]
                name = desc.attrs.get("name")
                kinematics = desc.attrs.get("kinematics")
                self.parameters.append([name, kinematics, sys.float_info.min, sys.float_info.max])
        else:
            RuntimeError('input file has no valid parameter descriptions: is it corrupted?')

        super().__init__()


    def __del__(self):
        self.file.close()


    """Retrieve data"""
    def data(self):
        if '/data/observables' not in self.file:
            RuntimeError('input file does not contain stored observables' % step)

        dataset = self.file['/data/observables']
        data = numpy.array(dataset[:])

        # adjust min,max range for each parameter based on data
        for i in range(0, len(self.parameters)):
            self.parameters[i][2] = numpy.min(data[:, i])
            self.parameters[i][3] = numpy.max(data[:, i])

        return data


""" open HDF5 data file regardless of file type """
def load_data_file(name):
    basename = os.path.basename(name)
    if basename.startswith('mcmc'):
        return MCMCDataFile(name)
    elif basename.startswith('pmc'):
        return PMCDataFile(name)
    elif basename.startswith('unc'):
        return UncertaintyDataFile(name)
    else:
        raise RuntimeError('cannot determine HDF5 file type based on the file name')

