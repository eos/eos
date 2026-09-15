# Copyright (c) 2026 Lorenz Gaertner
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

import os
import unittest
from pathlib import Path

os.environ.setdefault('EOS_SOURCE_DIR', str(Path(__file__).resolve().parents[3]))



def mcp_is_missing():
    try:
        import mcp
        return False
    except ImportError:
        return True


if not mcp_is_missing():
    from eos.mcp.server import observables, parameters, constraints, references, cli_reference_doc, defining_observables_doc
    from mcp.server.mcpserver.exceptions import ToolError


@unittest.skipIf(mcp_is_missing(), "Test requires the 'mcp' module")
class ObservablesToolTests(unittest.TestCase):

    def test_known_qualified_name(self):
        "A known qualified name is described in the result."

        result = observables(qualified_name='B->pilnu::BR')

        self.assertIn('B->pilnu::BR', result)

    def test_unknown_qualified_name(self):
        "An unknown qualified name raises a ToolError."

        with self.assertRaisesRegex(ToolError, "not known"):
            observables(qualified_name='prefix::Philipp')

    def test_prefix_filter(self):
        "The prefix filter narrows the returned observables."

        result = observables(prefix='B->pilnu')

        self.assertIn('B->pilnu::BR', result)
        self.assertNotIn('B->Dlnu::dBR/dq2', result)


@unittest.skipIf(mcp_is_missing(), "Test requires the 'mcp' module")
class ParametersToolTests(unittest.TestCase):

    def test_known_qualified_name(self):
        "A known qualified name is described in the result."

        result = parameters(qualified_name='life_time::tau')

        self.assertIn('life_time::tau', result)

    def test_unknown_qualified_name(self):
        "An unknown qualified name raises a ToolError."

        with self.assertRaisesRegex(ToolError, "Unknown parameter"):
            parameters(qualified_name='prefix::Philipp')

    def test_prefix_filter(self):
        "The prefix filter narrows the returned parameters."

        result = parameters(prefix='life_time')

        self.assertIn('life_time::tau', result)
        self.assertNotIn('mass::tau', result)


@unittest.skipIf(mcp_is_missing(), "Test requires the 'mcp' module")
class ConstraintsToolTests(unittest.TestCase):

    def test_known_qualified_name(self):
        "A known qualified name is described in the result."

        result = constraints(qualified_name='0->pipi::Abs{f_+}^2@Belle:2008C')

        self.assertIn('0->pipi::Abs{f_+}^2@Belle:2008C', result)

    def test_unknown_qualified_name(self):
        "An unknown qualified name raises a ToolError."

        with self.assertRaisesRegex(ToolError, "unknown"):
            constraints(qualified_name='prefix::Philipp')

    def test_prefix_filter(self):
        "The prefix filter narrows the returned constraints."

        result = constraints(prefix='0->pipi')

        self.assertIn('0->pipi::Abs{f_+}^2@Belle:2008C', result)
        self.assertTrue(all(qn.startswith('0->pipi') for qn in result))


@unittest.skipIf(mcp_is_missing(), "Test requires the 'mcp' module")
class ReferencesToolTests(unittest.TestCase):

    def test_known_key(self):
        "A known reference key is described in the result."

        result = references(key='AAGW:2001A')

        self.assertIn('AAGW:2001A', result)
        self.assertEqual(result['AAGW:2001A']['inspire_id'], 'Asatryan:2001zw')

    def test_unknown_key(self):
        "An unknown reference key raises a ToolError."

        with self.assertRaisesRegex(ToolError, "Unknown reference"):
            references(key='NoSuchRef:9999X')

    def test_year_filter(self):
        "The year filter narrows the returned references."

        result = references(year='2001')

        self.assertIn('AAGW:2001A', result)
        self.assertNotIn('ABHH:1999A', result)


@unittest.skipIf(mcp_is_missing(), "Test requires the 'mcp' module")
class DocResourceTests(unittest.TestCase):

    def test_cli_reference_doc(self):
        "The CLI reference doc is served with recognizable content."

        result = cli_reference_doc()

        self.assertIn('eos-analysis', result)

    def test_defining_observables_doc(self):
        "The defining-observables doc is served with recognizable content."

        result = defining_observables_doc()

        self.assertIn('observable', result.lower())


if __name__ == '__main__':
    unittest.main(verbosity=5)
