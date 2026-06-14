# Copyright (c) 2026 Danny van Dyk
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

import unittest

import eos
import os
from pathlib import Path

from eos.datasets_file_description import PublicLikelihoodDescription


class DataSetsTests(unittest.TestCase):

    # a prepared registry of data sets; pointing the storage directory here
    # avoids any download of the registry from the network
    _storage = str(Path(__file__).parent / 'datasets_TEST.d')

    def test_loading(self):
        "Load the registry of data sets from a prepared fixture."
        ds = eos.DataSets(storage_directory=self._storage)

        ids = {id for id, _ in ds.datasets()}
        self.assertEqual(ids, {'EXAMPLE:2024A', 'ANOTHER:2024B'})

        example = dict(ds.datasets())['EXAMPLE:2024A']
        self.assertEqual(example.authors, ['Jane Doe', 'John Roe'])
        self.assertEqual(example.title, 'Example data set')
        self.assertEqual(example.keywords, ['B decays', 'CKM'])
        self.assertEqual(example.eos_version, '1.0.0')
        self.assertEqual(example.doi, '10.5281/zenodo.0000000')
        # the likelihoods are deserialized into PublicLikelihoodDescription instances
        self.assertEqual(set(example.likelihoods.keys()), {'ckm', 'wet'})
        self.assertIsInstance(example.likelihoods['ckm'], PublicLikelihoodDescription)
        self.assertEqual(example.likelihoods['wet'].filetype, 'NabuLikelihood')

    def test_default_doi(self):
        "A data set without a DOI has 'doi' set to None."
        ds = eos.DataSets(storage_directory=self._storage)
        another = dict(ds.datasets())['ANOTHER:2024B']
        self.assertIsNone(another.doi)

    def test_path(self):
        "path() resolves an identifier against the storage directory."
        ds = eos.DataSets(storage_directory=self._storage)
        self.assertEqual(ds.path('EXAMPLE:2024A'), os.path.join(self._storage, 'EXAMPLE:2024A'))

    def test_exists(self):
        "exists() reports False for a data set that has not been downloaded."
        ds = eos.DataSets(storage_directory=self._storage)
        self.assertFalse(ds.exists('EXAMPLE:2024A'))

    def test_unknown_likelihood(self):
        "Requesting an unknown likelihood from a data set raises an error."
        ds = eos.DataSets(storage_directory=self._storage)
        with self.assertRaises(ValueError):
            ds.likelihood('EXAMPLE:2024A', 'does-not-exist')


if __name__ == '__main__':
    unittest.main(verbosity=5)
