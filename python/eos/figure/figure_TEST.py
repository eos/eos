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
import eos.data
import eos.figure
import numpy as np
import os
import shutil
import tempfile
import yaml

from eos.validation_context import ValidationContext
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

    def test_height_ratios(self):

        input = self._grid_yaml('[2, 1]', 2, extra='height_ratios: [4, 1]')
        figure = eos.figure.FigureFactory.from_yaml(input)
        self.assertEqual(list(figure._gridspec.get_height_ratios()), [4, 1])

    def test_width_ratios(self):

        input = self._grid_yaml('[1, 2]', 2, extra='width_ratios: [3, 1]')
        figure = eos.figure.FigureFactory.from_yaml(input)
        self.assertEqual(list(figure._gridspec.get_width_ratios()), [3, 1])

    def test_height_ratios_wrong_length(self):

        input = self._grid_yaml('[2, 1]', 2, extra='height_ratios: [1, 2, 3]')
        with self.assertRaises(Exception):
            eos.figure.FigureFactory.from_yaml(input)

    def test_width_ratios_wrong_length(self):

        input = self._grid_yaml('[1, 2]', 2, extra='width_ratios: [1]')
        with self.assertRaises(Exception):
            eos.figure.FigureFactory.from_yaml(input)

    def test_single_cell(self):

        # a 1x1 grid must work: subplots() squeezes to a bare Axes by default, which
        # would break the flatten() logic without squeeze=False.
        try:
            figure = eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 1]', 1))
            self.assertEqual(len(figure._axes), 1)
            figure.draw()
        except Exception as e:
            self.fail(f"Error when testing a 1x1 grid figure: {e}")

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

        # A boolean is rejected rather than silently treated as the int 0/1.
        with self.assertRaises(Exception):
            eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 2]', 2, extra='watermark_plot: true'))

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

    def test_draw_forwards_context(self):

        # Regression: GridFigure.draw must forward its context to each plot's prepare(), so that
        # items resolve relative paths against the analysis file's base directory. Previously the
        # context was dropped and every plot fell back to a fresh CWD-rooted context.
        from eos.analysis_file_context import AnalysisFileContext

        figure = eos.figure.FigureFactory.from_yaml(self._grid_yaml('[1, 2]', 2))
        context = AnalysisFileContext()

        seen = []
        for plot in figure.plots:
            original = plot.prepare
            def spy(ctx=None, _original=original):
                seen.append(ctx)
                return _original(ctx)
            plot.prepare = spy

        figure.draw(context=context)

        self.assertEqual(len(seen), len(figure.plots))
        for ctx in seen:
            self.assertIs(ctx, context)


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


class OverviewFigureTests(unittest.TestCase):

    # two observables constrained by 'B->KJpsi::BR@PDG:2020A' and 'B->K^*Jpsi::BR@PDG:2020A'
    OBSERVABLES = ['B->Kpsi::BR', 'B->K^*psi::BR']
    OPTIONS     = { 'psi': 'J/psi', 'q': 'u' }
    CONSTRAINTS = ['B->KJpsi::BR@PDG:2020A', 'B->K^*Jpsi::BR@PDG:2020A']

    def setUp(self):
        self._directory = tempfile.mkdtemp()
        self._prediction = os.path.join(self._directory, 'prediction')

        parameters  = eos.Parameters.Defaults()
        kinematics  = eos.Kinematics()
        observables = [
            eos.Observable.make(name, parameters, kinematics, eos.Options(self.OPTIONS))
            for name in self.OBSERVABLES
        ]

        rng     = np.random.default_rng(42)
        samples = np.column_stack([
            rng.normal(1.0e-3, 5.0e-5, 100),
            rng.normal(1.4e-3, 8.0e-5, 100),
        ])
        eos.data.Prediction.create(self._prediction, observables, samples, np.ones(100))

    def tearDown(self):
        shutil.rmtree(self._directory, ignore_errors=True)

    def _description(self, observables=None, sources=None, xaxis=None):
        return {
            'type':        'overview',
            'xaxis':       xaxis if xaxis is not None else { 'label': 'BR', 'range': [0.0, 2.0e-3] },
            'observables': observables if observables is not None else [
                { f'{name};psi=J/psi,q=u': name } for name in self.OBSERVABLES
            ],
            'sources':     sources if sources is not None else [
                { 'type': 'prediction', 'label': 'this work', 'datafiles': [self._prediction] },
                { 'type': 'constraint', 'label': 'PDG 2020',  'names': self.CONSTRAINTS },
            ],
        }

    def test_full(self):
        figure = eos.figure.FigureFactory.from_dict(**self._description())
        figure.draw()

        self.assertEqual(2, len(figure._sourced_entries))
        self.assertEqual(['this work', 'PDG 2020'], [entry['label'] for entry in figure._sourced_entries])
        # every source covers every observable, and the two rows are staggered symmetrically
        for entry in figure._sourced_entries:
            self.assertEqual(2, len(entry['positions']))
        first, second = figure._sourced_entries
        self.assertAlmostEqual(
            first['positions'][0][1] + second['positions'][0][1],
            2.0 * figure._yvals[0]
        )

    def test_constraint_uncertainties(self):
        "A serialized sigma such as '8e-05' is a str under YAML 1.1 and must still yield a number."
        description = self._description(sources=[
            { 'type': 'constraint', 'label': 'PDG 2020', 'names': self.CONSTRAINTS },
        ])
        figure = eos.figure.FigureFactory.from_dict(**description)
        figure.draw()

        (entry,) = figure._sourced_entries
        self.assertEqual([[2.7e-05, 2.7e-05], [8.0e-05, 8.0e-05]], entry['xerrors'])
        self.assertEqual([0.001006, 0.00143], [position[0] for position in entry['positions']])

    def test_grouped_observables(self):
        grouped = [
            [{ 'B->Kpsi::BR;psi=J/psi,q=u': 'A' }],
            [{ 'B->K^*psi::BR;psi=J/psi,q=u': 'B' }],
        ]
        flat, group_ids, labels = eos.figure.figure.OverviewFigure._group(grouped)

        self.assertEqual(['B->Kpsi::BR;psi=J/psi,q=u', 'B->K^*psi::BR;psi=J/psi,q=u'], flat)
        self.assertEqual([0, 1], group_ids)
        self.assertEqual(['A', 'B'], labels)

        figure = eos.figure.FigureFactory.from_dict(**self._description(observables=grouped))
        figure.draw()

    def test_normalize_to_constraint(self):
        xaxis = { 'label': 'BR / BR(exp)', 'range': [0.0, 2.0], 'normalize': 'constraint' }
        figure = eos.figure.FigureFactory.from_dict(**self._description(xaxis=xaxis))
        figure.draw()

        constraint = next(e for e in figure._sourced_entries if e['label'] == 'PDG 2020')
        for position in constraint['positions']:
            self.assertAlmostEqual(1.0, position[0])

    def test_invalid_description(self):
        for key, value in [
            ('xaxis',       None),
            ('observables', []),
            ('sources',     []),
        ]:
            description = self._description()
            description[key] = value
            with self.assertRaises(ValueError):
                eos.figure.FigureFactory.from_dict(**description)

        for xaxis in [
            { 'label': 'BR', 'range': [0.0, 1.0], 'normalize': 'unknown' },
            { 'label': 'BR', 'range': [0.0, 1.0], 'normalize': 'posterior' },
        ]:
            with self.assertRaises(ValueError):
                eos.figure.FigureFactory.from_dict(**self._description(xaxis=xaxis))

    def test_invalid_source(self):
        for source in [
            { 'type': 'unknown',    'label': 'S', 'names': self.CONSTRAINTS },
            { 'type': 'prediction', 'label': 'S' },
            { 'type': 'constraint', 'label': 'S' },
            { 'type': 'prediction', 'label': 'S', 'datafiles': [self._prediction], 'names': self.CONSTRAINTS },
            { 'type': 'constraint', 'label': 'S', 'names': self.CONSTRAINTS, 'datafiles': [self._prediction] },
        ]:
            with self.assertRaises(ValueError):
                eos.figure.FigureFactory.from_dict(**self._description(sources=[source]))

    def test_observable_missing_from_prediction(self):
        # the absent observable comes first, so that a check made only after the loop would miss it
        observables = [{ 'B->Dlnu::BR;l=mu,q=d': 'A' }, { 'B->Kpsi::BR;psi=J/psi,q=u': 'B' }]
        description = self._description(observables=observables, sources=[
            { 'type': 'prediction', 'label': 'this work', 'datafiles': [self._prediction] },
        ])
        figure = eos.figure.FigureFactory.from_dict(**description)

        with self.assertRaises(KeyError):
            figure.prepare()

    def test_observable_from_two_datafiles(self):
        description = self._description(sources=[
            { 'type': 'prediction', 'label': 'this work', 'datafiles': [self._prediction, self._prediction] },
        ])
        figure = eos.figure.FigureFactory.from_dict(**description)

        with self.assertRaises(ValueError):
            figure.prepare()

    def test_validate_semantics(self):
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()

        figure = eos.figure.FigureFactory.from_dict(**self._description())
        self.assertEqual([], list(figure.validate_semantics(ValidationContext(description))))

        figure = eos.figure.FigureFactory.from_dict(**self._description(
            observables=[{ 'test::unknown-observable': 'A' }],
            sources=[{ 'type': 'constraint', 'label': 'S', 'names': ['test::unknown-constraint@Unknown:2000A'] }],
        ))
        diagnostics = list(figure.validate_semantics(ValidationContext(description)))

        self.assertEqual(
            [('observables', 0), ('sources', 0, 'names')],
            [diagnostic.path for diagnostic in diagnostics]
        )

        xaxis = { 'label': 'BR', 'range': [0.0, 2.0], 'normalize': 'posterior', 'normalize-source': 'no such source' }
        figure = eos.figure.FigureFactory.from_dict(**self._description(xaxis=xaxis))
        diagnostics = list(figure.validate_semantics(ValidationContext(description)))

        self.assertEqual([('xaxis', 'normalize-source')], [diagnostic.path for diagnostic in diagnostics])


if __name__ == '__main__':
    unittest.main(verbosity=5)
