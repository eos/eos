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

from eos.datasets_file_description import PublicLikelihoodDescription, DataSetDescription


class PublicLikelihoodDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = PublicLikelihoodDescription.from_dict(filename='ckm/mixture_density', filetype='MixtureDensity')
        self.assertEqual(desc.filename, 'ckm/mixture_density')
        self.assertEqual(desc.filetype, 'MixtureDensity')

        desc = PublicLikelihoodDescription.from_dict(filename='wet/nabu', filetype='NabuLikelihood')
        self.assertEqual(desc.filetype, 'NabuLikelihood')

    def test_invalid_filetype(self):
        with self.assertRaises(ValueError):
            PublicLikelihoodDescription.from_dict(filename='foo', filetype='NotAFileType')


class DataSetDescriptionTests(unittest.TestCase):

    def test_from_dict(self):
        desc = DataSetDescription.from_dict(
            authors=['Jane Doe', 'John Roe'],
            title='Example data set',
            keywords=['B decays', 'CKM'],
            eos_version='1.0.0',
            doi='10.5281/zenodo.0000000',
            likelihoods={
                'ckm': {'filename': 'ckm/mixture_density', 'filetype': 'MixtureDensity'},
                'wet': {'filename': 'wet/nabu',            'filetype': 'NabuLikelihood'},
            },
        )
        self.assertEqual(desc.authors, ['Jane Doe', 'John Roe'])
        self.assertEqual(desc.title, 'Example data set')
        self.assertEqual(desc.keywords, ['B decays', 'CKM'])
        self.assertEqual(desc.eos_version, '1.0.0')
        self.assertEqual(desc.doi, '10.5281/zenodo.0000000')
        # the likelihoods are deserialized into PublicLikelihoodDescription instances
        self.assertEqual(set(desc.likelihoods.keys()), {'ckm', 'wet'})
        self.assertIsInstance(desc.likelihoods['ckm'], PublicLikelihoodDescription)
        self.assertEqual(desc.likelihoods['wet'].filetype, 'NabuLikelihood')

    def test_default_doi(self):
        desc = DataSetDescription.from_dict(
            authors=['Alice Smith'],
            title='Another data set',
            keywords=['rare decays'],
            eos_version='1.0.1',
            likelihoods={'main': {'filename': 'main/mixture_density', 'filetype': 'MixtureDensity'}},
        )
        self.assertIsNone(desc.doi)

    def test_invalid_nested_likelihood(self):
        "An invalid likelihood file type is rejected during deserialization."
        with self.assertRaises(ValueError):
            DataSetDescription.from_dict(
                authors=['Alice Smith'],
                title='Broken data set',
                keywords=[],
                eos_version='1.0.1',
                likelihoods={'main': {'filename': 'main/data', 'filetype': 'NotAFileType'}},
            )


if __name__ == '__main__':
    unittest.main(verbosity=5)
