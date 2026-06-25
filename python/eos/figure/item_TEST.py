# Copyright (c) 2025-2026 Danny van Dyk
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
from matplotlib import pyplot as plt
from matplotlib.container import ErrorbarContainer
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

class ItemColorCyclerTests(unittest.TestCase):

    def test_cycle(self):

        from eos.figure.item import ItemColorCycler

        # the cycler holds class-level state that other tests may have advanced
        ItemColorCycler.reset()

        colors = ItemColorCycler._colors

        # next_color() yields the colors in order, starting from the first
        for expected in colors:
            self.assertEqual(ItemColorCycler.next_color(), expected)

        # after a full cycle the index wraps around to the first color
        self.assertEqual(ItemColorCycler.next_color(), colors[0])

        # reset() returns the cycle to its starting point
        ItemColorCycler.reset()
        self.assertEqual(ItemColorCycler.next_color(), colors[0])

class ObservableItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: q2
            range: [0.1, 1.0]
            resolution: 100
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

    def test_legend(self):

        # a labelled observable contributes a single line entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::dBR/dq2'
        variable: q2
        range: [0.1, 1.0]
        resolution: 10
        label: 'foo'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][1], 'foo')
        # the observable is drawn opaque, so its swatch is opaque too
        self.assertIsNone(entries[0][0].get_alpha())

        # an unlabelled item contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::dBR/dq2'
        variable: q2
        range: [0.1, 1.0]
        resolution: 10
        """)
        self.assertEqual(list(item.legend()), [])

class ExpressionItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: expression
            expression: 'exp(-x**2) * sin(2 * pi * x)'
            range: [0.0, 6.28]
            resolution: 100
            label: 'foo'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'expression': {e}")

    def test_constant(self):

        # a constant expression is broadcast to the full grid of sample points
        item = eos.figure.ItemFactory.from_yaml("""
        type: expression
        expression: '1.5'
        range: [0.0, 1.0]
        resolution: 7
        """)
        item.prepare()
        self.assertEqual(item._yvalues.shape, item._xvalues.shape)
        self.assertTrue((item._yvalues == 1.5).all())

    def test_invalid(self):

        # a syntactically invalid expression is rejected at construction time
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: expression
            expression: 'sin(x'
            range: [0.0, 1.0]
            """)

        # an empty expression is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: expression
            expression: '   '
            range: [0.0, 1.0]
            """)

        # an inverted range is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: expression
            expression: 'x'
            range: [1.0, 0.0]
            """)

        # a name that is not exposed to the expression cannot be used
        with self.assertRaises(ValueError):
            item = eos.figure.ItemFactory.from_yaml("""
            type: expression
            expression: 'os.getcwd()'
            range: [0.0, 1.0]
            """)
            item.prepare()

        # a complex-valued result cannot be converted to float and is reported as a ValueError
        with self.assertRaises(ValueError):
            item = eos.figure.ItemFactory.from_yaml("""
            type: expression
            expression: '(-1.0)**0.5'
            range: [0.0, 1.0]
            """)
            item.prepare()

        # a result that cannot be broadcast onto the grid is reported as a ValueError
        with self.assertRaises(ValueError):
            item = eos.figure.ItemFactory.from_yaml("""
            type: expression
            expression: 'np.array([1.0, 2.0, 3.0])'
            range: [0.0, 1.0]
            resolution: 100
            """)
            item.prepare()

    def test_legend(self):

        # a labelled expression contributes a single line entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: expression
        expression: 'sin(x)'
        range: [0.0, 6.28]
        label: 'foo'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][1], 'foo')

        # an unlabelled expression contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: expression
        expression: 'sin(x)'
        range: [0.0, 6.28]
        """)
        self.assertEqual(list(item.legend()), [])

class UncertaintyBandItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: uncertainty
            label: '$\\ell=\\mu$'
            variable: 'q2'
            range: [0.02, 11.63]
            datafile: 'eos/data/prediction_TEST.d/predictions'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

class BinnedUncertaintyItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: uncertainty-binned
            label: '$\\ell=\\mu$'
            variable: 'q2'
            range: [0.02, 11.63]
            datafile: 'eos/data/prediction_TEST.d/predictions-binned'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

class ConstraintItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: constraint
            label: 'Belle 2015 $\\ell=e,\\, q=d$'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            rescale_by_width: true
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

    def test_legend(self):

        # a labelled constraint contributes a single error-bar entry (a capped bar rendered
        # by HandlerErrorbar), not just the central marker
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        label: 'Belle'
        constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
        observable: 'B->Dlnu::BR'
        variable: 'q2'
        rescale_by_width: true
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], ErrorbarContainer)
        self.assertTrue(entries[0][0].has_yerr)
        self.assertEqual(entries[0][1], 'Belle')

        # after prepare(), the binned constraint also carries an x error bar
        item.prepare()
        entry = item.legend()[0][0]
        self.assertIsInstance(entry, ErrorbarContainer)
        self.assertTrue(entry.has_xerr)
        self.assertTrue(entry.has_yerr)

class ConstraintResidueItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'constraint-residue'
            label: r'Belle 2015 $\\ell=e,\\, q=d$'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            parameters: {"mass::e": 1.0}
            rescale_by_width: true
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

    def test_legend(self):

        # a labelled residue item contributes a single error-bar entry (a capped bar
        # rendered by HandlerErrorbar), matching how the residues are drawn
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        label: 'Belle'
        constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
        observable: 'B->Dlnu::BR'
        variable: 'q2'
        rescale_by_width: true
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], ErrorbarContainer)
        self.assertTrue(entries[0][0].has_yerr)
        self.assertEqual(entries[0][1], 'Belle')

        # after prepare(), the binned constraint also carries an x error bar
        item.prepare()
        entry = item.legend()[0][0]
        self.assertTrue(entry.has_xerr)
        self.assertTrue(entry.has_yerr)

class TwoDimensionalConstraintItemTests(unittest.TestCase):

    def test_multivariate(self):

        # a bivariate constraint is drawn as one covariance ellipse per requested confidence level
        try:
            input = """
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K+C_K@BaBar:2008A'
            x: { observable: 'B->K^*gamma::S_K^*gamma' }
            y: { observable: 'B->K^*gamma::C_K^*gamma' }
            sigmas: [1.0, 2.0]
            color: 'C0'
            label: 'BaBar 2008'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint2D': {e}")

    def test_univariate(self):

        # a univariate (Gaussian) constraint with only 'x' is drawn as a vertical band
        try:
            input = """
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
            x: { observable: 'B->K^*gamma::S_K^*gamma' }
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing univariate item of type 'constraint2D': {e}")

        # the same constraint with only 'y' is drawn as a horizontal band
        try:
            input = """
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
            y: { observable: 'B->K^*gamma::S_K^*gamma' }
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing univariate item of type 'constraint2D': {e}")

        # the band must span the full orthogonal axis independently of the current limits:
        # it is anchored in axes-fraction coordinates, not in data coordinates. Non-default
        # limits are set first so that a data-coordinate band (the previous behaviour) would fail.
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
        x: { observable: 'B->K^*gamma::S_K^*gamma' }
        """)
        item.prepare()
        _, ax = plt.subplots()
        ax.set_ylim(5.0, 17.0)
        item.draw(ax)
        rect = ax.patches[-1]
        # the vertical (y) extent spans the whole axes: anchored at y=0 with height 1 in axes fraction
        self.assertEqual(rect.get_y(), 0.0)
        self.assertEqual(rect.get_height(), 1.0)

        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
        y: { observable: 'B->K^*gamma::S_K^*gamma' }
        """)
        item.prepare()
        _, ax = plt.subplots()
        ax.set_xlim(5.0, 17.0)
        item.draw(ax)
        rect = ax.patches[-1]
        # the horizontal (x) extent spans the whole axes: anchored at x=0 with width 1 in axes fraction
        self.assertEqual(rect.get_x(), 0.0)
        self.assertEqual(rect.get_width(), 1.0)

    def test_invalid(self):

        # neither 'x' nor 'y' specified is rejected at construction time
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
            """)

        # an axis specification without an 'observable' key is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
            x: {}
            """)

        # specifying both 'x' and 'y' for a univariate constraint is rejected during prepare()
        with self.assertRaises(ValueError):
            item = eos.figure.ItemFactory.from_yaml("""
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K@BaBar:2008A'
            x: { observable: 'B->K^*gamma::S_K^*gamma' }
            y: { observable: 'B->K^*gamma::C_K^*gamma' }
            """)
            item.prepare()

        # specifying only one axis for a multivariate constraint is rejected during prepare()
        with self.assertRaises(ValueError):
            item = eos.figure.ItemFactory.from_yaml("""
            type: 'constraint2D'
            constraint: 'B^0->K^*0gamma::S_K+C_K@BaBar:2008A'
            x: { observable: 'B->K^*gamma::S_K^*gamma' }
            """)
            item.prepare()

    def test_legend(self):

        # a labelled 2D constraint contributes a single patch entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B^0->K^*0gamma::S_K+C_K@BaBar:2008A'
        x: { observable: 'B->K^*gamma::S_K^*gamma' }
        y: { observable: 'B->K^*gamma::C_K^*gamma' }
        label: 'BaBar 2008'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Rectangle)
        self.assertEqual(entries[0][1], 'BaBar 2008')

        # an unlabelled 2D constraint contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B^0->K^*0gamma::S_K+C_K@BaBar:2008A'
        x: { observable: 'B->K^*gamma::S_K^*gamma' }
        y: { observable: 'B->K^*gamma::C_K^*gamma' }
        """)
        self.assertEqual(list(item.legend()), [])

class OneDimensionalHistogramItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'histogram1D'
            variable: 'CKM::abs(V_ub)'
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

class TwoDimensionalHistogramItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'histogram2D'
            variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

    def test_legend(self):

        # a 2D density has no faithful swatch and contributes no entry, even when labelled
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'histogram2D'
        variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
        datafile: 'eos/data/importance_samples_TEST.d/samples'
        label: 'should be ignored'
        """)
        self.assertEqual(list(item.legend()), [])

class OneDimensionalKernelDensityItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'kde1D'
            bandwidth: 1.3
            level: 2
            xsamples: 150
            variable: 'CKM::abs(V_ub)'
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

class TwoDimensionalKernelDensityItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'kde2D'
            bandwidth: 3
            contours: ['lines', 'areas', 'labels']
            levels: [1, 3]
            variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

class BandItemTests(unittest.TestCase):

    def test_full(self):

        # x values only
        try:
            input = """
            type: band
            x: [-0.1, +0.1]
            color: 'blue'
            label: 'foo'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'band': {e}")

        # y values only
        try:
            input = """
            type: band
            y: [-0.1, +0.1]
            color: 'orange'
            alpha: 0.5
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'band': {e}")

    def test_legend(self):

        # a labelled band contributes a single patch entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: band
        x: [-0.1, +0.1]
        color: 'blue'
        label: 'foo'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Rectangle)
        self.assertEqual(entries[0][1], 'foo')

        # an unlabelled band contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: band
        y: [-0.1, +0.1]
        color: 'orange'
        """)
        self.assertEqual(list(item.legend()), [])

class SignalPDFItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: signal-pdf
            label: 'PDF ($\\ell=\\mu$)'
            pdf: 'B->Dlnu::P(q2);l=mu'
            variable: 'q2'
            range: [0.02, 11.60]
            resolution: 100
            kinematics:
              q2_min:  0.02
              q2_max: 11.60
            color: 'C0'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'signal-pdf': {e}")

class ComplexPlaneItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'complex-plane'
            observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
            variables: ['Re{q2}', 'Im{q2}']
            ranges: [[-1.0, +1.0], [-1.0, +1.0]]
            resolution: 10
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'complex-plane': {e}")

    def test_legend(self):

        # a pseudocolor plot has no faithful swatch and contributes no entry, even when labelled
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'complex-plane'
        observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
        variables: ['Re{q2}', 'Im{q2}']
        ranges: [[-1.0, +1.0], [-1.0, +1.0]]
        resolution: 10
        label: 'should be ignored'
        """)
        self.assertEqual(list(item.legend()), [])

class ErrorBarsItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'errorbars'
            positions: [[1, 2], [2, 3], [3, 5]]
            xerrors: [0.5, 0.5, 0.5]
            yerrors: [0.2, [0.2, 0.3], 0.5]
            color: 'black'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'signal-pdf': {e}")

    def test_legend(self):

        # a labelled error-bar item contributes a single error-bar entry (a capped bar
        # rendered by HandlerErrorbar), not just the central marker
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'errorbars'
        positions: [[1, 2], [2, 3]]
        yerrors: [0.2, 0.3]
        marker: 'o'
        label: 'data'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], ErrorbarContainer)
        # only the y error is present, so the swatch carries the y bar but not the x bar
        self.assertFalse(entries[0][0].has_xerr)
        self.assertTrue(entries[0][0].has_yerr)
        self.assertEqual(entries[0][1], 'data')

        # an item with both x and y errors carries both bars
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'errorbars'
        positions: [[1, 2]]
        xerrors: [0.4]
        yerrors: [0.3]
        label: 'data'
        """)
        entries = item.legend()
        self.assertIsInstance(entries[0][0], ErrorbarContainer)
        self.assertTrue(entries[0][0].has_xerr)
        self.assertTrue(entries[0][0].has_yerr)

        # an unlabelled item contributes no legend entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'errorbars'
        positions: [[1, 2]]
        yerrors: [0.3]
        """)
        self.assertEqual(list(item.legend()), [])

class VerticalLineItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'vertical'
            x: 3.8
            color: 'gray'
            linestyle: 'dashed'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'vertical': {e}")

    def test_legend(self):

        # a labelled vertical line contributes a single line entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'vertical'
        x: 3.8
        label: 'threshold'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][1], 'threshold')
        # the line is drawn with the item's alpha, so the swatch matches it
        self.assertEqual(entries[0][0].get_alpha(), item.alpha)

        # an unlabelled vertical line contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'vertical'
        x: 3.8
        """)
        self.assertEqual(list(item.legend()), [])

if __name__ == '__main__':
    unittest.main(verbosity=5)
