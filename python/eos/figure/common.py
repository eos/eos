# Copyright (c) 2024-2025 Danny van Dyk
# Copyright (c) 2024      Méril Reboud
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
from eos.analysis_file_description import AnalysisFileContext
from eos.deserializable import Deserializable
from eos.figure.item import ItemColorCycler

import eos
import os

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

    def __post_init__(self):
        if self.color is None:
            self.color = ItemColorCycler.next_color()

    def prepare(self, context:AnalysisFileContext=None):
        context = AnalysisFileContext() if context is None else context

        if hasattr(self, '_datafile'):
            return self._datafile

        path = context.data_path(self.path)
        os.path.exists(path) or eos.error(f"Data file '{path}' does not exist")
        name = os.path.split(path)[-1]

        if name == 'samples':
            # The datafile refers to importance (parameter) samples
            self._datafile = eos.data.ImportanceSamples(path)
            self._type = 'samples'
        elif name.startswith('pred-'):
            # The datafile refers to observable predictions
            self._datafile = eos.data.Prediction(path)
            self._type = 'prediction'
        else:
            eos.error(f"Data file '{path}' has an unsupported format")
            raise NotImplementedError

        self._variables = set(self._datafile.lookup_table.keys())
        return self._datafile

    @property
    def variables(self):
        self.prepare()
        return self._variables

    def labels(self, variables):
        self.prepare()

        if self._type == 'samples':
            p = eos.Parameters()
            return [p[dist].latex() for dist in variables]
        elif self._type == 'prediction':
            o = eos.Observables()
            return ["$" + o[dist].latex() + "$" for dist in variables]
        else:
            eos.error(f"Data file '{self.path}' has an unsupported format")
            raise NotImplementedError
