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
import os

from eos.analysis_file_context import AnalysisFileContext


class AnalysisFileContextTests(unittest.TestCase):

    def test_data_path(self):
        "Check that data_path resolves a relative path against the base directory."
        ctx = AnalysisFileContext()
        path = ctx.data_path('some/relative/path')
        self.assertTrue(os.path.isabs(path))
        self.assertTrue(path.endswith(os.path.join('some', 'relative', 'path')))

    def test_invalid_base_directory(self):
        "Check that a non-existent or non-directory base directory is rejected."
        with self.assertRaises(ValueError):
            AnalysisFileContext(base_directory='/this/does/not/exist/at/all')
        # a path that exists but is a file, not a directory
        with self.assertRaises(ValueError):
            AnalysisFileContext(base_directory=__file__)


if __name__ == '__main__':
    unittest.main(verbosity=5)
