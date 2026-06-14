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
import eos.figure
import os

from eos.analysis_file_description import AnalysisFileContext

class DataFileTests(unittest.TestCase):

    def test_defaults(self):
        "Check the construction of a DataFile and its default values."
        from eos.figure.data import DataFile

        df = DataFile.from_dict(path='path/to/data', label='my label')
        self.assertEqual(df.path, 'path/to/data')
        self.assertEqual(df.label, 'my label')
        self.assertFalse(df.kde)
        # when no color is given, one is assigned from the color cycler
        self.assertIsNotNone(df.color)

    def test_custom(self):
        "Check the construction of a DataFile with explicit color and KDE flag."
        from eos.figure.data import DataFile

        df = DataFile.from_dict(path='path/to/data', label='my label', color='red', kde=True)
        self.assertEqual(df.color, 'red')
        self.assertTrue(df.kde)

    def test_prepare_samples(self):
        "Check that a DataFile can load a file of importance samples."
        from eos.figure.data import DataFile

        context = AnalysisFileContext(base_directory=os.environ['SOURCE_DIR'])
        df = DataFile.from_dict(path='eos/data/importance_samples_TEST.d/samples', label='samples')
        df.prepare(context=context)

        self.assertEqual(df._type, 'samples')
        # the variables are populated from the data file's lookup table
        self.assertGreater(len(df.variables), 0)

    def test_prepare_missing_file(self):
        "Check that loading a non-existent data file raises an error."
        from eos.figure.data import DataFile

        context = AnalysisFileContext(base_directory=os.environ['SOURCE_DIR'])
        df = DataFile.from_dict(path='eos/data/this_does_not_exist.d/samples', label='missing')
        with self.assertRaises(ValueError):
            df.prepare(context=context)


if __name__ == '__main__':
    unittest.main(verbosity=5)
