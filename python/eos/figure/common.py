# Copyright (c) 2024 Danny van Dyk
# Copyright (c) 2024 Méril Reboud
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

from dataclasses import dataclass, field
from eos.deserializable import Deserializable

import eos

@dataclass(kw_only=True)
class DataFile(Deserializable):
    r""" Collects the relevant information to load a data file, to extract data, variables, and labels from it.

    :param path: Path to the data file (either in eos.Prediction or eos.ImportanceSample format) that will be used.
    :type path: str
    :param label: Title to appear in a figure's legend in association with this data file.
    :type label: str
    :param color: Color to be used in the plot.
    :type color: str (optional)
    """
    path:str
    label:str
    color:str=None

    def prepare(self):
        if hasattr(self, '_datafile'):
            return self._datafile

        if not os.path.exists(self.path):
            eos.error(f"Data file '{self.path}' does not exist")

        abspath = os.path.abspath(self.datafilepath)
        name = os.path.split(abspath)[-1]

        if name == 'samples':
            # The datafile refers to importance (parameter) samples
            self._datafile = eos.data.ImportanceSamples(abspath)
            self._type = 'samples'
        elif name.startswith('pred-'):
            # The datafile refers to observable predictions
            self._datafile = eos.data.Prediction(abspath)
            self._type = 'prediction'
        else:
            eos.error(f"Data file '{self.datafilepath}' has an unsupported format")
            raise NotImplementedError

        self._variables = set(self._datafile.lookup_table.keys())
        return self._datafile

    @property
    def variables(self):
        self.prepare()
        return self._variables

    def labels(self, variables):
        self.load()

        if self._type == 'samples':
            p = eos.Parameters()
            return [p[dist].latex() for dist in variables]
        elif self._type == 'prediction':
            o = eos.Observables()
            return ["$" + o[dist].latex() + "$" for dist in variables]
        else:
            eos.error(f"Data file '{self.path}' has an unsupported format")
            raise NotImplementedError
