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

from matplotlib import pyplot as plt

class WatermarkTests(unittest.TestCase):

    def test_defaults(self):
        "Check the default position and the derived placement of a Watermark."
        from eos.figure.common import Watermark

        wm = Watermark.from_dict()
        self.assertEqual(wm.position, 'upper right')
        self.assertFalse(wm.preliminary)
        # 'upper right' places the watermark in the top-right corner
        self.assertAlmostEqual(wm._x, 0.96)
        self.assertAlmostEqual(wm._y, 0.96)
        self.assertEqual(wm._halign, 'right')
        self.assertEqual(wm._valign, 'top')
        # a non-preliminary watermark shows the EOS version
        self.assertEqual(wm._color, 'OrangeRed')
        self.assertTrue(wm._version.startswith('v'))

    def test_custom(self):
        "Check a custom position and the preliminary flag."
        from eos.figure.common import Watermark

        wm = Watermark.from_dict(position='lower left', preliminary=True)
        self.assertAlmostEqual(wm._x, 0.04)
        self.assertAlmostEqual(wm._y, 0.04)
        self.assertEqual(wm._halign, 'left')
        self.assertEqual(wm._valign, 'bottom')
        # a preliminary watermark is shown in red
        self.assertEqual(wm._color, 'red')
        self.assertEqual(wm._version, 'Preliminary')

        # a centered horizontal position is also valid
        self.assertAlmostEqual(Watermark.from_dict(position='upper center')._x, 0.5)

    def test_invalid_position(self):
        "Check that invalid positions are rejected."
        from eos.figure.common import Watermark

        # invalid horizontal position
        with self.assertRaises(ValueError):
            Watermark.from_dict(position='upper sideways')
        # invalid vertical position
        with self.assertRaises(ValueError):
            Watermark.from_dict(position='center right')
        # missing separator between the two positions
        with self.assertRaises(ValueError):
            Watermark.from_dict(position='upperright')

    def test_draw(self):
        "Check that a Watermark can be drawn onto a set of axes."
        from eos.figure.common import Watermark

        try:
            _, ax = plt.subplots()
            Watermark.from_dict().draw(ax)
        except Exception as e:
            self.fail(f"Error when drawing Watermark: {e}")


if __name__ == '__main__':
    unittest.main(verbosity=5)
