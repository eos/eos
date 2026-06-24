# Copyright (c) 2024-2026 Danny van Dyk
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
import yaml

from matplotlib import pyplot as plt

class PlotTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'e' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'mu' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
            """
            plot = eos.figure.PlotFactory.from_yaml(input)
            plot.prepare()
            fig, ax = plt.subplots()
            plot.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing plot of type '2D': {e}")

class EmptyPlotTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'empty'
            """
            plot = eos.figure.PlotFactory.from_yaml(input)
            plot.prepare()
            fig, ax = plt.subplots()
            plot.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing plot of type 'empty': {e}")

class LegendTests(unittest.TestCase):

    def test_defaults_and_custom(self):
        from eos.figure.plot import Legend

        self.assertEqual(Legend.from_dict().position, 'best')
        self.assertEqual(Legend.from_dict(position='lower left').position, 'lower left')

    def test_draw(self):
        from eos.figure.plot import Legend

        try:
            _, ax = plt.subplots()
            ax.plot([0, 1], [0, 1], label='line')
            Legend.from_dict(position='upper right').draw(ax)
            # also exercise the explicit-entries branch
            handle, = ax.plot([0, 1], [1, 0], label='other')
            Legend().draw(ax, entries=[(handle, 'other')])
        except Exception as e:
            self.fail(f"Error when drawing Legend: {e}")

    def test_draw_empty_entries(self):
        from eos.figure.plot import Legend

        # an empty list of entries (no labelled items) must not draw a legend and must not
        # raise: matplotlib >= 3.10 rejects empty handles/labels instead of warning.
        _, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        try:
            Legend().draw(ax, entries=[])
        except Exception as e:
            self.fail(f"Error when drawing Legend with empty entries: {e}")
        self.assertIsNone(ax.get_legend())


class GridTests(unittest.TestCase):

    def test_defaults_and_custom(self):
        from eos.figure.plot import Grid

        grid = Grid.from_dict()
        self.assertFalse(grid.visible)
        self.assertEqual(grid.axis, 'both')

        grid = Grid.from_dict(visible=True, axis='x')
        self.assertTrue(grid.visible)
        self.assertEqual(grid.axis, 'x')

    def test_invalid_axis(self):
        from eos.figure.plot import Grid

        with self.assertRaises(ValueError):
            Grid.from_dict(axis='diagonal')

    def test_draw(self):
        from eos.figure.plot import Grid

        try:
            _, ax = plt.subplots()
            Grid.from_dict(visible=True, axis='both').draw(ax)
        except Exception as e:
            self.fail(f"Error when drawing Grid: {e}")


class XTicksTests(unittest.TestCase):

    def test_defaults_and_custom(self):
        from eos.figure.plot import XTicks

        ticks = XTicks.from_dict()
        self.assertTrue(ticks.minor)
        self.assertEqual(ticks.position, 'bottom')
        self.assertTrue(ticks.visible)

        ticks = XTicks.from_dict(minor=False, position='both', visible=False)
        self.assertFalse(ticks.minor)
        self.assertEqual(ticks.position, 'both')
        self.assertFalse(ticks.visible)

    def test_invalid_position(self):
        from eos.figure.plot import XTicks

        with self.assertRaises(ValueError):
            XTicks.from_dict(position='left')

    def test_draw(self):
        from eos.figure.plot import XTicks

        try:
            _, ax = plt.subplots()
            XTicks.from_dict(position='top').draw(ax)
            _, ax = plt.subplots()
            XTicks.from_dict(visible=False).draw(ax)
        except Exception as e:
            self.fail(f"Error when drawing XTicks: {e}")

    def test_format(self):
        import matplotlib.ticker
        from eos.figure.plot import XTicks

        # defaults to None (matplotlib's default formatter)
        self.assertIsNone(XTicks.from_dict().format)
        self.assertEqual(XTicks.from_dict(format='%.2f').format, '%.2f')

        # a format string installs a FormatStrFormatter on the major ticks
        _, ax = plt.subplots()
        XTicks.from_dict(format='%.2f').draw(ax)
        formatter = ax.xaxis.get_major_formatter()
        self.assertIsInstance(formatter, matplotlib.ticker.FormatStrFormatter)
        self.assertEqual(formatter.fmt, '%.2f')

        # without a format, the major formatter is left untouched (not a FormatStrFormatter)
        _, ax = plt.subplots()
        XTicks.from_dict().draw(ax)
        self.assertNotIsInstance(ax.xaxis.get_major_formatter(), matplotlib.ticker.FormatStrFormatter)


class YTicksTests(unittest.TestCase):

    def test_defaults_and_custom(self):
        from eos.figure.plot import YTicks

        ticks = YTicks.from_dict()
        self.assertTrue(ticks.minor)
        self.assertEqual(ticks.position, 'left')
        self.assertTrue(ticks.visible)

        ticks = YTicks.from_dict(minor=False, position='both', visible=False)
        self.assertFalse(ticks.minor)
        self.assertEqual(ticks.position, 'both')
        self.assertFalse(ticks.visible)

    def test_invalid_position(self):
        from eos.figure.plot import YTicks

        with self.assertRaises(ValueError):
            YTicks.from_dict(position='bottom')

    def test_draw(self):
        from eos.figure.plot import YTicks

        try:
            _, ax = plt.subplots()
            YTicks.from_dict(position='right').draw(ax)
            _, ax = plt.subplots()
            YTicks.from_dict(visible=False).draw(ax)
        except Exception as e:
            self.fail(f"Error when drawing YTicks: {e}")

    def test_format(self):
        import matplotlib.ticker
        from eos.figure.plot import YTicks

        # defaults to None (matplotlib's default formatter)
        self.assertIsNone(YTicks.from_dict().format)
        self.assertEqual(YTicks.from_dict(format='%.2f').format, '%.2f')

        # a format string installs a FormatStrFormatter on the major ticks
        _, ax = plt.subplots()
        YTicks.from_dict(format='%.2f').draw(ax)
        formatter = ax.yaxis.get_major_formatter()
        self.assertIsInstance(formatter, matplotlib.ticker.FormatStrFormatter)
        self.assertEqual(formatter.fmt, '%.2f')

        # without a format, the major formatter is left untouched (not a FormatStrFormatter)
        _, ax = plt.subplots()
        YTicks.from_dict().draw(ax)
        self.assertNotIsInstance(ax.yaxis.get_major_formatter(), matplotlib.ticker.FormatStrFormatter)


class XAxisTests(unittest.TestCase):

    def test_defaults(self):
        from eos.figure.plot import XAxis, XTicks

        xaxis = XAxis.from_dict()
        self.assertIsNone(xaxis.label)
        self.assertIsNone(xaxis.range)
        self.assertEqual(xaxis.scale, 'linear')
        # the ticks default to an XTicks instance
        self.assertIsInstance(xaxis.ticks, XTicks)

    def test_nested_ticks_and_range(self):
        from eos.figure.plot import XAxis, XTicks

        xaxis = XAxis.from_dict(label='$q^2$', range=[1, 6], ticks={'position': 'both'})
        # nested ticks dict is deserialized into an XTicks instance
        self.assertIsInstance(xaxis.ticks, XTicks)
        self.assertEqual(xaxis.ticks.position, 'both')
        # the range is converted to a tuple of floats
        self.assertEqual(xaxis.range, (1.0, 6.0))
        self.assertTrue(all(isinstance(x, float) for x in xaxis.range))

    def test_invalid_range(self):
        from eos.figure.plot import XAxis

        with self.assertRaises(ValueError):
            XAxis.from_dict(range=[1, 2, 3])

    def test_draw(self):
        from eos.figure.plot import XAxis

        try:
            _, ax = plt.subplots()
            XAxis.from_dict(label='mass', unit='GeV', range=[0.0, 1.0], scale='linear').draw(ax)
        except Exception as e:
            self.fail(f"Error when drawing XAxis: {e}")


class YAxisTests(unittest.TestCase):

    def test_defaults(self):
        from eos.figure.plot import YAxis, YTicks

        yaxis = YAxis.from_dict()
        self.assertIsNone(yaxis.label)
        self.assertIsNone(yaxis.range)
        self.assertEqual(yaxis.scale, 'linear')
        # the ticks default to a YTicks instance
        self.assertIsInstance(yaxis.ticks, YTicks)

    def test_nested_ticks_and_range(self):
        from eos.figure.plot import YAxis, YTicks

        yaxis = YAxis.from_dict(label='$d\\mathcal{B}/dq^2$', range=[0, 5], ticks={'position': 'right'})
        # nested ticks dict is deserialized into a YTicks instance
        self.assertIsInstance(yaxis.ticks, YTicks)
        self.assertEqual(yaxis.ticks.position, 'right')
        # the range is converted to a tuple of floats
        self.assertEqual(yaxis.range, (0.0, 5.0))
        self.assertTrue(all(isinstance(y, float) for y in yaxis.range))

    def test_invalid_range(self):
        from eos.figure.plot import YAxis

        with self.assertRaises(ValueError):
            YAxis.from_dict(range=[1, 2, 3])

    def test_draw(self):
        from eos.figure.plot import YAxis

        try:
            _, ax = plt.subplots()
            YAxis.from_dict(label='rate', unit='GeV', range=[1.0e-3, 1.0], scale='log').draw(ax)
        except Exception as e:
            self.fail(f"Error when drawing YAxis: {e}")


if __name__ == '__main__':
    unittest.main(verbosity=5)
