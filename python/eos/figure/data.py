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
import numpy as np

@dataclass(kw_only=True)
class DataFile(Deserializable):
    r""" Collects the relevant information to load a data file, to extract data, variables, and labels from it.

    :param path: Path to the data file (either in eos.Prediction or eos.ImportanceSample format) that will be used.
    :type path: str
    :param label: Title to appear in a figure's legend in association with this data file.
    :type label: str
    :param color: Color to be used in the plot.
    :type color: str (optional)
    :param kde: A boolean determining whether to use kernel density estimates (KDE) to visualize the distributions. Defaults to False, in which case histograms are used.
    :type kde: bool
    """
    path:str
    label:str
    color:str=None
    kde:bool=False

    def __post_init__(self):
        if self.color is None:
            self.color = ItemColorCycler.next_color()

    def prepare(self, context:AnalysisFileContext=None):
        if hasattr(self, '_datafile'):
            return self._datafile

        context = AnalysisFileContext() if context is None else context

        path = context.data_path(self.path)
        if not os.path.exists(path):
            raise ValueError(f"Data file '{path}' does not exist")
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

        self._variables = list(self._datafile.lookup_table.keys())
        return self._datafile

    @property
    def variables(self):
        self.prepare()
        return self._variables

    @property
    def empirical_range(self):
        self.prepare()
        lower  = np.quantile(self._datafile.samples, 0.01, weights=self._datafile.weights, method='inverted_cdf', axis=0)
        median = np.quantile(self._datafile.samples, 0.50, weights=self._datafile.weights, method='inverted_cdf', axis=0)
        upper  = np.quantile(self._datafile.samples, 0.99, weights=self._datafile.weights, method='inverted_cdf', axis=0)

        min = median - 1.30 * (median - lower)
        max = median + 1.30 * (upper - median)
        return (min, max)

    def labels(self, variables):
        self.prepare()

        if self._type == 'samples':
            p = eos.Parameters()
            return [p[dist].latex() for dist in variables]
        elif self._type == 'prediction':
            o = eos.Observables()
            return ["$" + o[dist.split(";")[0]].latex() + "$" for dist in variables]
        else:
            eos.error(f"Data file '{self.path}' has an unsupported format")
            raise NotImplementedError
