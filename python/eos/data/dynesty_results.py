# Copyright (c) 2022-2023 Filip Novak
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
import os
import numpy as _np
import dynesty
import yaml


class DynestyResults:
    r"""Represents the results of a nested-sampling run with dynesty, stored on disk.

    Wraps a :class:`dynesty.results.Results` object together with the descriptions of the varied
    parameters. Instances are created either by reading existing results from disk (passing their
    ``path`` to the constructor) or by writing new results with :meth:`create`.

    :ivar type: The type identifier of the data object, always ``'DynestyResults'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its column index in :attr:`samples`.
    :ivar results: The underlying :class:`dynesty.results.Results` object.
    :ivar samples: The samples in parameter space as a 2D array.
    :ivar weights: The importance weights on a linear scale, derived from the nested-sampling log-weights.
    """

    def __init__(self, path):
        """ Read Results object (in the dynesty.results module) from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Description file {f} does not exist or is not a file')

        with open(f) as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'DynestyResults':
            raise RuntimeError(f'Path {path} not pointing to a DynestyResults file')

        self.type = 'DynestyResults'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'dynesty_results.npy')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Dynesty results file {f} does not exist or is not a file')

        res_dict = _np.load(f, allow_pickle=True).item()
        if "blob" in res_dict:
            res_dict.pop('blob')
        self.results = dynesty.results.Results(res_dict)
        self.samples = self.results.samples
        self.weights = _np.exp(self.results.logwt - self.results.logz[-1])


    @staticmethod
    def create(path, parameters, results):
        """ Write a new Results object (in the dynesty.results module) to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param results: The results of a nested sampling run.
        :type results: dynesty.results.Results
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'DynestyResults'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)

        res_dict = results.asdict()
        _np.save(os.path.join(path, 'dynesty_results.npy'), res_dict)
