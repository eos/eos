# Copyright (c) 2022-2024 Méril Reboud
# Copyright (c) 2023 Danny van Dyk
# Copyright (c) 2023 Stephan Kürten
# Copyright (c) 2024 Lorenz Gärtner
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
import yaml


class Mode:
    def __init__(self, path):
        """ Read a posterior's (local) mode from a file.

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

        if not description['type'] == 'Mode':
            raise RuntimeError(f'Path {path} not pointing to a Mode file')

        self.type = 'Mode'
        self.varied_parameters = description['parameters']
        self.mode = description['mode']
        self.pvalue = description['pvalue']
        self.local_pvalues = description['local_pvalues']
        self.global_chi2 = None
        if 'global_chi2' in description:
            self.global_chi2 = description['global_chi2']
        self.dof = None
        if 'dof' in description:
            self.dof = description['dof']

    @staticmethod
    def create(path, parameters, mode, pvalue, local_pvalues, global_chi2, dof):
        """ Write a new Mode object to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param mode: The mode to be stored.
        :type mode: numpy.ndarray
        :param pvalue: The p-value of the mode.
        :type pvalue: float
        :param local_pvalues: The local p-values of the mode.
        :type local_pvalues: dict
        :param global_chi2: The global chi2 value of the mode.
        :type global_chi2: float
        :param dof: The degrees of freedom of the mode.
        :type dof: float
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'Mode'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]
        description['mode'] = mode.tolist()
        description['pvalue'] = float(pvalue) if pvalue is not None else None
        description['local_pvalues'] = local_pvalues
        description['global_chi2'] = global_chi2
        description['dof'] = dof

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)
