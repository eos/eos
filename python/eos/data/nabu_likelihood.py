# Copyright (c) 2024 Danny van Dyk
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


class NabuLikelihood:
    r"""Represents a likelihood serialized with the nabu package, stored on disk.

    Wraps a :class:`nabu.Likelihood` object together with the descriptions of the varied parameters.
    Instances are created either by reading an existing likelihood from disk (passing its ``path`` to
    the constructor) or by writing a new likelihood with :meth:`create`. Reading requires the optional
    nabu module.

    :ivar type: The type identifier of the data object, always ``'NabuLikelihood'``.
    :ivar varied_parameters: The descriptions (name, min, max) of the varied parameters.
    :ivar lookup_table: A mapping from each parameter name to its index in :attr:`varied_parameters`.
    :ivar likelihood: The underlying :class:`nabu.Likelihood` object.
    """

    def __init__(self, path):
        """ Read a nabu serialized likelihood from disk.

        :param path: Path to the storage location.
        :type path: str
        """
        import nabu

        if not os.path.exists(path) or not os.path.isdir(path):
            raise RuntimeError(f'Path {path} does not exist or is not a directory')

        f = os.path.join(path, 'description.yaml')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Description file {f} does not exist or is not a file')

        with open(f) as df:
            description = yaml.load(df, Loader=yaml.SafeLoader)

        if not description['type'] == 'NabuLikelihood':
            raise RuntimeError(f'Path {path} not pointing to a NabuLikelihood file')

        self.type = 'NabuLikelihood'
        self.varied_parameters = description['parameters']
        self.lookup_table = { item['name']: idx for idx, item in enumerate(self.varied_parameters) }

        f = os.path.join(path, 'likelihood.nabu')
        if not os.path.exists(f) or not os.path.isfile(f):
            raise RuntimeError(f'Nabu likelihood file {f} does not exist or is not a file')

        self.likelihood = nabu.Likelihood.load(f)


    @staticmethod
    def create(path, parameters, likelihood):
        """ Write a new Results object (in the dynesty.results module) to disk.

        :param path: Path to the storage location, which will be created as a directory.
        :type path: str
        :param parameters: Parameter descriptions as a 1D array of shape (N, ).
        :type parameters: list or iterable of eos.Parameter
        :param likelihood: The likelihood object created by nabu.
        :type results: nabu.Likelihood or descendant
        """
        description = {}
        description['version'] = eos.__version__
        description['type'] = 'NabuLikelihood'
        description['parameters'] = [{
            'name': p.name(),
            'min': p.min(),
            'max': p.max()
        } for p in parameters]

        os.makedirs(path, exist_ok=True)
        with open(os.path.join(path, 'description.yaml'), 'w') as description_file:
            yaml.dump(description, description_file, default_flow_style=False)

        likelihood.save(os.path.join(path, 'likelihood.nabu'))
