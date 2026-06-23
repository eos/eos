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

class SingleFigureTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'single'
            plot:
              legend:
                position: 'lower left'
              xaxis:
                label: 'q^2'
              yaxis:
                label: 'dBR/dq^2'
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
            figure = eos.figure.FigureFactory.from_yaml(input)
            figure.draw()
        except Exception as e:
            self.fail(f"Error when testing figure of type 'single': {e}")


class InsetFigureTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: inset
            plot:
              xaxis: { label: '$q^2$',                range: [0.0, 11.6]   }
              yaxis: { label: '$d\\mathcal{B}/dq^2$', range: [0.0, 5.4e-3] }
              items:
                - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'mu' },  label: '$\\ell = \\mu$',
                    variable: 'q2', range: [0.02, 11.6], resolution: 5800
                  }
                - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'tau' }, label: '$\\ell = \\tau$',
                    variable: 'q2', range: [3.17, 11.6], resolution: 421
                  }
            inset:
              position: [0.5, 0.5]
              size: [0.485, 0.48]
              plot:
                xaxis: { range: [0.0, 0.25],   ticks: { visible: false } }
                yaxis: { range: [0.0, 4.4e-3], ticks: { visible: false } }
                items:
                  - { type: 'observable', observable: 'B->Dlnu::dBR/dq2', options: { 'l': 'mu' },  label: '$\\ell = \\mu$',
                      variable: 'q2', range: [0.01, 0.25], resolution: 230
                    }
            watermark:
              position: 'upper left'
            """
            figure = eos.figure.FigureFactory.from_yaml(input)
            figure.draw()
        except Exception as e:
            self.fail(f"Error when testing figure of type 'inset': {e}")


class GridFigureTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'grid'
            shape: [1, 2]
            plots:
              - xaxis: { label: '$q^2$' }
                yaxis: { label: '$d\\mathcal{B}/dq^2$' }
                items:
                  - type: 'observable'
                    observable: 'B->Dlnu::dBR/dq2'
                    options: { 'l': 'e' }
                    variable: 'q2'
                    range: [0.1, 1.0]
                    resolution: 100
              - xaxis: { label: '$q^2$' }
                yaxis: { label: '$d\\mathcal{B}/dq^2$' }
                items:
                  - type: 'observable'
                    observable: 'B->Dlnu::dBR/dq2'
                    options: { 'l': 'mu' }
                    variable: 'q2'
                    range: [0.1, 1.0]
                    resolution: 100
            """
            figure = eos.figure.FigureFactory.from_yaml(input)
            figure.draw()
        except Exception as e:
            self.fail(f"Error when testing figure of type 'grid': {e}")

    def test_size(self):

        input = """
        type: 'grid'
        shape: [1, 2]
        size: [8.0, 6.0]
        plots:
          - xaxis: { label: '$q^2$' }
            yaxis: { label: '$d\\mathcal{B}/dq^2$' }
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'e' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
          - xaxis: { label: '$q^2$' }
            yaxis: { label: '$d\\mathcal{B}/dq^2$' }
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'mu' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
        """
        figure = eos.figure.FigureFactory.from_yaml(input)
        size = figure._figure.get_size_inches()
        self.assertAlmostEqual(size[0], 8.0)
        self.assertAlmostEqual(size[1], 6.0)

    def test_default_size(self):

        # When 'size' is omitted, the figure falls back to (3.0 * ncol, 3.0 * nrow).
        input = """
        type: 'grid'
        shape: [1, 2]
        plots:
          - xaxis: { label: '$q^2$' }
            yaxis: { label: '$d\\mathcal{B}/dq^2$' }
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'e' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
          - xaxis: { label: '$q^2$' }
            yaxis: { label: '$d\\mathcal{B}/dq^2$' }
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'mu' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
        """
        figure = eos.figure.FigureFactory.from_yaml(input)
        size = figure._figure.get_size_inches()
        self.assertAlmostEqual(size[0], 6.0) # 3.0 * ncol = 3.0 * 2
        self.assertAlmostEqual(size[1], 3.0) # 3.0 * nrow = 3.0 * 1

    @staticmethod
    def _has_watermark(ax):
        return any('EOS' in t.get_text() for t in ax.texts)

    @staticmethod
    def _grid_yaml(shape, nplots, extra=''):
        plots = "".join("""
          - xaxis: { label: '$q^2$' }
            yaxis: { label: '$d\\mathcal{B}/dq^2$' }
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'e' }
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100""" for _ in range(nplots))
        return f"""
        type: 'grid'
        shape: {shape}
        {extra}
        plots:{plots}
        """

    def test_watermark_plot_single(self):

        # A flattened index selects a single panel; only that panel is stamped.
        input = self._grid_yaml('[1, 2]', 2, extra='watermark_plot: 1')
        figure = eos.figure.FigureFactory.from_yaml(input)
        figure.draw()
        self.assertFalse(self._has_watermark(figure._axes[0]))
        self.assertTrue(self._has_watermark(figure._axes[1]))

    def test_watermark_plot_tuple(self):

        # A 2D (row, col) address resolves to flat index row * ncol + col.
        # (1, 0) in a 2x2 grid -> flat index 2.
        input = self._grid_yaml('[2, 2]', 4, extra='watermark_plot: [1, 0]')
        figure = eos.figure.FigureFactory.from_yaml(input)
        figure.draw()
        for idx in range(4):
            self.assertEqual(self._has_watermark(figure._axes[idx]), idx == 2)

    def test_watermark_plot_default(self):

        # Without 'watermark_plot', every panel is stamped (backward-compatible).
        input = self._grid_yaml('[1, 2]', 2)
        figure = eos.figure.FigureFactory.from_yaml(input)
        figure.draw()
        self.assertTrue(self._has_watermark(figure._axes[0]))
        self.assertTrue(self._has_watermark(figure._axes[1]))

    def test_watermark_plot_out_of_range(self):

        # An out-of-range flattened index is rejected at construction.
        with self.assertRaises(Exception):
            eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 2]', 2, extra='watermark_plot: 5'))

        # An out-of-range 2D address is rejected as well.
        with self.assertRaises(Exception):
            eos.figure.FigureFactory.from_yaml(self._grid_yaml('[2, 2]', 4, extra='watermark_plot: [0, 9]'))

    @staticmethod
    def _two_range_grid(extra=''):
        # a single column with two panels carrying different x-ranges
        return f"""
        type: 'grid'
        shape: [2, 1]
        {extra}
        plots:
          - xaxis: {{ label: '$q^2$', range: [0.1, 1.0] }}
            yaxis: {{ label: '$d\\mathcal{{B}}/dq^2$' }}
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: {{ 'l': 'e' }}
                variable: 'q2'
                range: [0.1, 1.0]
                resolution: 100
          - xaxis: {{ label: '$q^2$', range: [2.0, 5.0] }}
            yaxis: {{ label: '$d\\mathcal{{B}}/dq^2$' }}
            items:
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: {{ 'l': 'mu' }}
                variable: 'q2'
                range: [2.0, 5.0]
                resolution: 100
        """

    def test_tight_layout_disabled(self):

        try:
            figure = eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 2]', 2, extra='tight_layout: false'))
            self.assertFalse(figure.tight_layout)
            figure.draw()
        except Exception as e:
            self.fail(f"Error when drawing grid figure with tight_layout disabled: {e}")

    def test_shared_axes_x(self):

        # With shared_axes 'x' the panels are joined in x and end up with one common x-range.
        figure = eos.figure.FigureFactory.from_yaml(self._two_range_grid(extra="shared_axes: 'x'"))
        self.assertTrue(figure._axes[0].get_shared_x_axes().joined(figure._axes[0], figure._axes[1]))
        self.assertFalse(figure._axes[0].get_shared_y_axes().joined(figure._axes[0], figure._axes[1]))
        figure.draw()
        self.assertEqual(figure._axes[0].get_xlim(), figure._axes[1].get_xlim())

    def test_shared_axes_y(self):

        # 'y' shares the y-axis per row. The 2x1 single-column grid has one panel
        # per row, so use a 1x2 grid to exercise row sharing.
        figure = eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 2]', 2, extra="shared_axes: 'y'"))
        self.assertTrue(figure._axes[0].get_shared_y_axes().joined(figure._axes[0], figure._axes[1]))
        self.assertFalse(figure._axes[0].get_shared_x_axes().joined(figure._axes[0], figure._axes[1]))

    def test_shared_axes_both(self):

        figure = eos.figure.FigureFactory.from_yaml(self._grid_yaml('[2, 2]', 4, extra="shared_axes: 'both'"))
        # column-shared x and row-shared y
        self.assertTrue(figure._axes[0].get_shared_x_axes().joined(figure._axes[0], figure._axes[2]))
        self.assertTrue(figure._axes[0].get_shared_y_axes().joined(figure._axes[0], figure._axes[1]))

    def test_shared_axes_default_independent(self):

        # Without shared_axes (the default) the panels keep independent axes.
        figure = eos.figure.FigureFactory.from_yaml(self._two_range_grid())
        self.assertFalse(figure._axes[0].get_shared_x_axes().joined(figure._axes[0], figure._axes[1]))
        self.assertFalse(figure._axes[0].get_shared_y_axes().joined(figure._axes[0], figure._axes[1]))

    def test_shared_axes_invalid(self):

        with self.assertRaises(Exception):
            eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 2]', 2, extra="shared_axes: 'diagonal'"))


class CornerFigureTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'corner'
            contents:
              - path: 'path/to/datafile'
                label: 'label 1'
                color: 'red'
              - path: 'path/to/anotherdatafile'
                label: 'label 2'
                color: 'blue'
            variables: ['var1', 'var2']
            """
            figure = eos.figure.FigureFactory.from_yaml(input)
        except Exception as e:
            self.fail(f"Error when testing figure of type 'corner': {e}")


if __name__ == '__main__':
    unittest.main(verbosity=5)
