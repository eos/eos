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

    def test_outside(self):
        from eos.figure.plot import Legend

        legend = Legend.from_dict(position='upper center', outside=True, ncol=2)
        self.assertTrue(legend.outside)
        self.assertEqual(2, legend.ncol)

        _, ax = plt.subplots()
        handle, = ax.plot([0, 1], [0, 1], label='line')
        legend.draw(ax, entries=[(handle, 'line')])
        self.assertIsNotNone(ax.get_legend())

    def test_invalid_outside_position(self):
        from eos.figure.plot import Legend

        # 'best' has no edge for the legend to sit against
        with self.assertRaises(ValueError):
            Legend.from_dict(position='best', outside=True)

    def test_invalid_ncol(self):
        from eos.figure.plot import Legend

        with self.assertRaises(ValueError):
            Legend.from_dict(ncol=0)


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

    def test_scaling_factor(self):
        import matplotlib.ticker
        from eos.figure.plot import XTicks

        # defaults to None (no rescaling)
        self.assertIsNone(XTicks.from_dict().scaling_factor)

        # scaling only: a FuncFormatter divides the tick value and renders it with the '%g' default
        _, ax = plt.subplots()
        XTicks.from_dict(scaling_factor=1e-3).draw(ax)
        formatter = ax.xaxis.get_major_formatter()
        self.assertIsInstance(formatter, matplotlib.ticker.FuncFormatter)
        self.assertEqual(formatter(4.2e-3, 0), '4.2')

        # scaling and format: the rescaled value is rendered with the given printf format
        _, ax = plt.subplots()
        XTicks.from_dict(scaling_factor=1e-3, format='%.0f').draw(ax)
        formatter = ax.xaxis.get_major_formatter()
        self.assertIsInstance(formatter, matplotlib.ticker.FuncFormatter)
        self.assertEqual(formatter(4.2e-3, 0), '4')

        # a zero scaling factor is rejected
        with self.assertRaises(ValueError):
            XTicks.from_dict(scaling_factor=0.0)

        # combining the scaling factor with a logarithmic axis warns the user
        _, ax = plt.subplots()
        ax.set_xscale('log')
        with self.assertLogs('EOS', level='WARNING'):
            XTicks.from_dict(scaling_factor=1e-3).draw(ax)

    def test_locations(self):
        import matplotlib.ticker
        from eos.figure.plot import XTicks

        # defaults to None (automatic tick placement)
        self.assertIsNone(XTicks.from_dict().locations)

        # explicit locations install a FixedLocator, overriding the automatic placement
        _, ax = plt.subplots()
        ax.set_xlim(-2.5, 2.5)
        XTicks.from_dict(locations=[-2, 0, 2]).draw(ax)
        self.assertIsInstance(ax.xaxis.get_major_locator(), matplotlib.ticker.FixedLocator)
        self.assertEqual(list(ax.get_xticks()), [-2, 0, 2])


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

    def test_scaling_factor(self):
        import matplotlib.ticker
        from eos.figure.plot import YTicks

        # defaults to None (no rescaling)
        self.assertIsNone(YTicks.from_dict().scaling_factor)

        # scaling only: a FuncFormatter divides the tick value and renders it with the '%g' default
        _, ax = plt.subplots()
        YTicks.from_dict(scaling_factor=1e-3).draw(ax)
        formatter = ax.yaxis.get_major_formatter()
        self.assertIsInstance(formatter, matplotlib.ticker.FuncFormatter)
        self.assertEqual(formatter(4.2e-3, 0), '4.2')

        # scaling and format: the rescaled value is rendered with the given printf format
        _, ax = plt.subplots()
        YTicks.from_dict(scaling_factor=1e-3, format='%.0f').draw(ax)
        formatter = ax.yaxis.get_major_formatter()
        self.assertIsInstance(formatter, matplotlib.ticker.FuncFormatter)
        self.assertEqual(formatter(4.2e-3, 0), '4')

        # a zero scaling factor is rejected
        with self.assertRaises(ValueError):
            YTicks.from_dict(scaling_factor=0.0)

        # combining the scaling factor with a logarithmic axis warns the user
        _, ax = plt.subplots()
        ax.set_yscale('log')
        with self.assertLogs('EOS', level='WARNING'):
            YTicks.from_dict(scaling_factor=1e-3).draw(ax)

    def test_locations(self):
        import matplotlib.ticker
        from eos.figure.plot import YTicks

        # defaults to None (automatic tick placement)
        self.assertIsNone(YTicks.from_dict().locations)

        # explicit locations install a FixedLocator, overriding the automatic placement
        _, ax = plt.subplots()
        ax.set_ylim(-2.5, 2.5)
        YTicks.from_dict(locations=[-2, 0, 2]).draw(ax)
        self.assertIsInstance(ax.yaxis.get_major_locator(), matplotlib.ticker.FixedLocator)
        self.assertEqual(list(ax.get_yticks()), [-2, 0, 2])

    def test_labels(self):
        from eos.figure.plot import YTicks

        ticks = YTicks.from_dict(locations=[1.0, 2.0], labels=['$A$', '$B$'], minor=False)
        _, ax = plt.subplots()
        ax.set_ylim(0.0, 3.0)
        ticks.draw(ax)

        ax.figure.canvas.draw()
        self.assertEqual(['$A$', '$B$'], [label.get_text() for label in ax.get_yticklabels()])

    def test_invalid_labels(self):
        from eos.figure.plot import YTicks

        # labels without locations have nothing to attach to
        with self.assertRaises(ValueError):
            YTicks.from_dict(labels=['$A$'])

        with self.assertRaises(ValueError):
            YTicks.from_dict(locations=[1.0, 2.0], labels=['$A$'])


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


class OverviewPlotTests(unittest.TestCase):

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

    def _description(self, observables=None, sources=None, xaxis=None, **kwargs):
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
            **kwargs,
        }

    def _plot(self, **kwargs):
        "The overview plot alone, as deserialized by the plot factory."
        return eos.figure.PlotFactory.from_dict(**self._description(**kwargs))

    def _drawn(self, **kwargs):
        "The overview plot, prepared and drawn onto its own axes."
        plot = self._plot(**kwargs)
        _, ax = plt.subplots()
        plot.prepare()
        plot.draw(ax)
        return plot, ax

    def test_full(self):
        plot, _ = self._drawn()

        self.assertEqual(2, len(plot._sourced_entries))
        self.assertEqual(['this work', 'PDG 2020'], [entry['label'] for entry in plot._sourced_entries])
        # every source covers every observable, and the two rows are staggered symmetrically
        for entry in plot._sourced_entries:
            self.assertEqual(2, len(entry['positions']))
        first, second = plot._sourced_entries
        self.assertAlmostEqual(
            first['positions'][0][1] + second['positions'][0][1],
            2.0 * plot.yaxis.ticks.locations[0]
        )
        self.assertEqual(self.OBSERVABLES, plot.yaxis.ticks.labels)

    def test_within_figure(self):
        "The overview plot must work wherever a plot is accepted, here a single-plot figure."
        figure = eos.figure.FigureFactory.from_dict(**{
            'type': 'single',
            'size': [6.4, 2.4],
            'plot': self._description(),
        })
        figure.draw()

        # the enclosing figure validates through the plot, prefixing its diagnostics
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()
        figure = eos.figure.FigureFactory.from_dict(**{
            'type': 'single',
            'plot': self._description(observables=[{ 'test::unknown-observable': 'A' }]),
        })
        diagnostics = list(figure.validate_semantics(ValidationContext(description)))

        self.assertEqual([('plot', 'observables', 0)], [diagnostic.path for diagnostic in diagnostics])

    def test_constraint_uncertainties(self):
        "A serialized sigma such as '8e-05' is a str under YAML 1.1 and must still yield a number."
        plot, _ = self._drawn(sources=[
            { 'type': 'constraint', 'label': 'PDG 2020', 'names': self.CONSTRAINTS },
        ])

        (entry,) = plot._sourced_entries
        self.assertEqual([[2.7e-05, 2.7e-05], [8.0e-05, 8.0e-05]], entry['xerrors'])
        self.assertEqual([0.001006, 0.00143], [position[0] for position in entry['positions']])

    def test_grouped_observables(self):
        grouped = [
            [{ 'B->Kpsi::BR;psi=J/psi,q=u': 'A' }],
            [{ 'B->K^*psi::BR;psi=J/psi,q=u': 'B' }],
        ]
        flat, group_ids, labels = eos.figure.plot.OverviewPlot._group(grouped)

        self.assertEqual(['B->Kpsi::BR;psi=J/psi,q=u', 'B->K^*psi::BR;psi=J/psi,q=u'], flat)
        self.assertEqual([0, 1], group_ids)
        self.assertEqual(['A', 'B'], labels)

        self._drawn(observables=grouped)

    def test_normalize_to_constraint(self):
        xaxis = { 'label': 'BR / BR(exp)', 'range': [0.0, 2.0] }
        plot, _ = self._drawn(xaxis=xaxis, normalize='constraint')

        constraint = next(e for e in plot._sourced_entries if e['label'] == 'PDG 2020')
        for position in constraint['positions']:
            self.assertAlmostEqual(1.0, position[0])

    def test_invalid_description(self):
        for key, value in [
            ('xaxis',       { 'label': 'BR' }),
            ('observables', []),
            ('sources',     []),
        ]:
            description = self._description()
            description[key] = value
            with self.assertRaises(ValueError):
                eos.figure.PlotFactory.from_dict(**description)

        for kwargs in [
            { 'normalize': 'unknown' },
            { 'normalize': 'posterior' },
            { 'observable_offset': 0.5 },
        ]:
            with self.assertRaises(ValueError):
                self._plot(**kwargs)

    def test_invalid_source(self):
        for source in [
            { 'type': 'unknown',    'label': 'S', 'names': self.CONSTRAINTS },
            { 'type': 'prediction', 'label': 'S' },
            { 'type': 'constraint', 'label': 'S' },
            { 'type': 'prediction', 'label': 'S', 'datafiles': [self._prediction], 'names': self.CONSTRAINTS },
            { 'type': 'constraint', 'label': 'S', 'names': self.CONSTRAINTS, 'datafiles': [self._prediction] },
        ]:
            with self.assertRaises(ValueError):
                self._plot(sources=[source])

    def test_observable_missing_from_prediction(self):
        # the absent observable comes first, so that a check made only after the loop would miss it
        observables = [{ 'B->Dlnu::BR;l=mu,q=d': 'A' }, { 'B->Kpsi::BR;psi=J/psi,q=u': 'B' }]
        description = self._description(observables=observables, sources=[
            { 'type': 'prediction', 'label': 'this work', 'datafiles': [self._prediction] },
        ])
        plot = eos.figure.PlotFactory.from_dict(**description)

        with self.assertRaises(KeyError):
            plot.prepare()

    def test_observable_from_two_datafiles(self):
        description = self._description(sources=[
            { 'type': 'prediction', 'label': 'this work', 'datafiles': [self._prediction, self._prediction] },
        ])
        plot = eos.figure.PlotFactory.from_dict(**description)

        with self.assertRaises(ValueError):
            plot.prepare()

    def test_validate_semantics(self):
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()

        plot = self._plot()
        self.assertEqual([], list(plot.validate_semantics(ValidationContext(description))))

        plot = self._plot(
            observables=[{ 'test::unknown-observable': 'A' }],
            sources=[{ 'type': 'constraint', 'label': 'S', 'names': ['test::unknown-constraint@Unknown:2000A'] }],
        )
        diagnostics = list(plot.validate_semantics(ValidationContext(description)))

        self.assertEqual(
            [('observables', 0), ('sources', 0, 'names')],
            [diagnostic.path for diagnostic in diagnostics]
        )

        plot = self._plot(normalize='posterior', normalize_source='no such source')
        diagnostics = list(plot.validate_semantics(ValidationContext(description)))

        self.assertEqual([('normalize_source',)], [diagnostic.path for diagnostic in diagnostics])

if __name__ == '__main__':
    unittest.main(verbosity=5)
