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
import tempfile

from eos.analysis_file_context import AnalysisFileContext
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

    def test_levels(self):

        # the 0% level (the peak) is prepended by default; its threshold must be the maximum
        # density, and all thresholds must stay within the data range (i.e. not the ~1.0 that
        # solving for P=0 numerically would return for a normalized density)
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'kde2D'
        bandwidth: 3
        levels: [68, 95]
        variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
        datafile: 'eos/data/importance_samples_TEST.d/samples'
        """)
        item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))

        # the 0% level is present by construction
        self.assertIn(0, item.levels)

        plevels = item._plevels()
        pdf_max = item._pdf.max()
        # the 0% level maps to the peak density, which is the largest threshold
        self.assertEqual(max(plevels), pdf_max)
        # every threshold lies within the data range (0, pdf_max], never the out-of-range ~1.0
        for p in plevels:
            self.assertGreater(p, 0.0)
            self.assertLessEqual(p, pdf_max)

class TwoDimensionalContoursItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'contours2D'
            bins: 50
            contours: ['lines', 'areas', 'labels']
            levels: [68, 95, 99]
            variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            label: 'posterior'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'contours2D': {e}")

    def test_levels(self):

        # the 0% level (the peak) is prepended by default; its threshold must be the maximum
        # density, and all thresholds must stay within the data range (i.e. not the ~1.0 that
        # solving for P=0 numerically would return for a normalized histogram)
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'contours2D'
        bins: 50
        levels: [68, 95]
        variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
        datafile: 'eos/data/importance_samples_TEST.d/samples'
        """)
        item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))

        # the 0% level is present by construction
        self.assertIn(0, item.levels)

        plevels = item._plevels()
        pdf_max = item._pdf.max()
        # the 0% level maps to the peak density, which is the largest threshold
        self.assertEqual(max(plevels), pdf_max)
        # every threshold lies within the data range (0, pdf_max], never the out-of-range ~1.0
        for p in plevels:
            self.assertGreater(p, 0.0)
            self.assertLessEqual(p, pdf_max)

    def test_invalid(self):

        # fewer than two bins is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'contours2D'
            bins: 1
            variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """)

        # an out-of-range credibility level is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'contours2D'
            levels: [68, 100]
            variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """)

        # an unsupported contour type is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'contours2D'
            contours: ['lines', 'blobs']
            variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
            datafile: 'eos/data/importance_samples_TEST.d/samples'
            """)

    def test_legend(self):

        # with areas, the swatch is a filled rectangle
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'contours2D'
        contours: ['lines', 'areas']
        variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
        datafile: 'eos/data/importance_samples_TEST.d/samples'
        label: 'posterior'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Rectangle)
        self.assertEqual(entries[0][1], 'posterior')

        # without areas, the swatch is a line
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'contours2D'
        contours: ['lines']
        variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
        datafile: 'eos/data/importance_samples_TEST.d/samples'
        label: 'posterior'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)

        # an unlabelled item contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'contours2D'
        variables: ['CKM::abs(V_ub)', 'B->pi::f_+(0)@BCL2008']
        datafile: 'eos/data/importance_samples_TEST.d/samples'
        """)
        self.assertEqual(list(item.legend()), [])

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

class PointItemTests(unittest.TestCase):

    def test_full(self):

        try:
            input = """
            type: 'point'
            x: 0.0
            y: 0.261
            marker: 'o'
            markersize: 12
            color: 'C0'
            label: 'LCSR (Bharucha 2012)'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'point': {e}")

    def test_invalid(self):

        # a point requires both 'x' and 'y'
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'point'
            x: 0.0
            """)

        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: 'point'
            y: 0.0
            """)

    def test_legend(self):

        # a labelled point contributes a single marker entry, drawn as an open (unfilled) marker
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'point'
        x: 0.0
        y: 0.261
        marker: 'o'
        label: 'foo'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][1], 'foo')
        self.assertEqual(entries[0][0].get_marker(), 'o')
        self.assertEqual(entries[0][0].get_markerfacecolor(), 'none')

        # an unlabelled point contributes no entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'point'
        x: 0.0
        y: 0.261
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

# Fixtures and helpers shared by the coverage tests below.
_ITEM_TEST_D = 'eos/figure/item_TEST.d'
# The full lookup keys (name + options + kinematics) of the two observables stored in the
# 'pred-uniq' Prediction fixture; used to exercise the full-name lookup of the data-driven items.
_PRED_FP = 'B->D::f_+(q2);form-factors=BCL2008,model=CKM[q2=1]'
_PRED_F0 = 'B->D::f_0(q2);form-factors=BCL2008,model=CKM[q2=1]'

def _source_context():
    "An analysis file context rooted at the test source directory, where the fixtures live."
    return AnalysisFileContext(base_directory=os.environ['SOURCE_DIR'])

class ItemFactoryTests(unittest.TestCase):

    def test_invalid(self):

        # a description without a 'type' key is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("observable: 'B->Dlnu::dBR/dq2'")

        # a description naming an unknown item type is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: 'not-a-real-item-type'")

class ObservableItemValidationTests(unittest.TestCase):

    def test_fixed_kinematics_options(self):

        # explicit options and a fixed (spectator) kinematic variable are declared without error
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::BR'
        variable: q2_min
        range: [0.1, 1.0]
        resolution: 3
        options: { l: 'mu' }
        fixed_kinematics: { q2_max: 11.6 }
        """)
        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)

    def test_fixed_parameters(self):

        # explicit fixed parameters are applied to the parameter set
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::dBR/dq2'
        variable: q2
        range: [0.1, 1.0]
        resolution: 3
        fixed_parameters: { 'mass::mu': 0.10566 }
        """)
        item.prepare()

    def test_variable_is_parameter(self):

        # 'variable' may name an EOS parameter (rather than a kinematic variable); the observable is
        # then evaluated as a function of that parameter. This exercises the QualifiedName branch of
        # __post_init__, which previously referred to a non-existent 'self.parameters' attribute.
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::dBR/dq2'
        options: { model: 'CKM' }
        variable: 'CKM::abs(V_cb)'
        fixed_kinematics: { q2: 5.0 }
        range: [0.038, 0.043]
        resolution: 5
        """)
        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)

        # the x-axis parameter is resolved to an EOS Parameter and is actually varied: dBR/dq2 grows
        # with |V_cb| (it is proportional to |V_cb|^2), and the parameter is left at the last grid point
        yvalues = list(item._yvalues)
        self.assertTrue(all(lo < hi for lo, hi in zip(yvalues, yvalues[1:])),
                        f"observable did not vary monotonically with the parameter: {yvalues}")
        self.assertAlmostEqual(float(item._variable), 0.043)

    def test_variable_invalid_parameter(self):

        # a 'variable' that is a well-formed QualifiedName but not a registered parameter is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: 'Not::a-parameter'
            range: [0.1, 0.2]
            """)

    def test_invalid(self):

        # a fixed kinematic variable that the observable does not declare is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: q2
            range: [0.1, 1.0]
            fixed_kinematics: { not_a_kinematic: 1.0 }
            """)

        # a 'variable' that is neither a kinematic variable nor a valid parameter name is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("""
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: 'not a valid name'
            range: [0.1, 1.0]
            """)

    def test_missing_fixed_parameters_from_file(self):

        # Regression (issue #1181): an item referencing a fixed_parameters_from_file that does not
        # exist must still load, because the file is read only when the item is prepared for drawing
        # (it may be the output of an earlier task, e.g. find-mode). A file that is genuinely missing
        # at prepare() then raises.
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::dBR/dq2'
        variable: q2
        range: [0.1, 1.0]
        resolution: 3
        fixed_parameters_from_file: 'does-not-exist.yaml'
        """)
        with self.assertRaises(RuntimeError):
            item.prepare()

    def test_fixed_parameters_from_file(self):

        # The parameter file is resolved relative to the context's base directory and applied in
        # prepare(), not at construction.
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'params.yaml'), 'w') as f:
                f.write("'mass::mu':\n  central: 0.10566\n")
            item = eos.figure.ItemFactory.from_yaml("""
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: q2
            range: [0.1, 1.0]
            resolution: 3
            fixed_parameters_from_file: 'params.yaml'
            """)
            item.prepare(context=AnalysisFileContext(base_directory=tmp))
            _, ax = plt.subplots()
            item.draw(ax)

class ExpressionItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        # a null range is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: expression\nexpression: 'x'\nrange: null")

        # a non-positive resolution is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: expression\nexpression: 'x'\nrange: [0.0, 1.0]\nresolution: 0")

class UncertaintyBandItemValidationTests(unittest.TestCase):

    def test_band_normalization(self):

        # a single band type given as a string is normalized to a set
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nband: 'median'")
        self.assertEqual(item.band, {'median'})

        # a list of band types is normalized to a set
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nband: ['area', 'outer']")
        self.assertEqual(item.band, {'area', 'outer'})

    def test_invalid(self):

        # a 'band' that is neither a string, list, nor set is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nband: 5")

        # an unknown band type is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nband: 'nonsense'")

        # an unknown interpolation type is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\ninterpolation: 'quadratic'")

        # a range that is not a pair is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nrange: [1.0]")

        # an inverted range is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nrange: [2.0, 1.0]")

class BinnedUncertaintyItemValidationTests(unittest.TestCase):

    def test_missing_binned_kinematics(self):

        # the data file must expose '<variable>_min' and '<variable>_max' for each prediction
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'nonexistent'\n"
                                                "datafile: 'eos/data/prediction_TEST.d/predictions-binned'")
        with self.assertRaises(RuntimeError):
            item.prepare(context=_source_context())

class OneDimensionalHistogramItemValidationTests(unittest.TestCase):

    def test_invalid_bins(self):

        # fewer than two bins is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: histogram1D\nvariable: 'x'\ndatafile: 'x'\nbins: 1")

    def test_prepare_errors(self):

        ctx = _source_context()

        # a data file whose name is neither 'samples' nor 'pred-*' has an unsupported format
        item = eos.figure.ItemFactory.from_yaml("type: histogram1D\nvariable: 'x'\ndatafile: '" + _ITEM_TEST_D + "'")
        with self.assertRaises(NotImplementedError):
            item.prepare(context=ctx)

        # a variable absent from the data file is rejected
        item = eos.figure.ItemFactory.from_yaml("type: histogram1D\nvariable: 'not-a-variable'\n"
                                                "datafile: 'eos/data/importance_samples_TEST.d/samples'")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

    def test_prediction(self):

        # the Prediction ('pred-*') branch, addressed by the full observable name
        item = eos.figure.ItemFactory.from_yaml("type: histogram1D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variable: '" + _PRED_FP + "'")
        item.prepare(context=_source_context())
        _, ax = plt.subplots()
        item.draw(ax)

class TwoDimensionalHistogramItemValidationTests(unittest.TestCase):

    def test_invalid_bins(self):

        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: histogram2D\nvariables: ['a', 'b']\ndatafile: 'x'\nbins: 1")

    def test_prepare_errors(self):

        ctx = _source_context()

        item = eos.figure.ItemFactory.from_yaml("type: histogram2D\nvariables: ['a', 'b']\ndatafile: '" + _ITEM_TEST_D + "'")
        with self.assertRaises(NotImplementedError):
            item.prepare(context=ctx)

        item = eos.figure.ItemFactory.from_yaml("type: histogram2D\nvariables: ['nope1', 'nope2']\n"
                                                "datafile: 'eos/data/importance_samples_TEST.d/samples'")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

    def test_prediction(self):

        item = eos.figure.ItemFactory.from_yaml("type: histogram2D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['" + _PRED_FP + "', '" + _PRED_F0 + "']")
        item.prepare(context=_source_context())
        _, ax = plt.subplots()
        item.draw(ax)

class OneDimensionalKernelDensityItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        # a credibility level outside (0, 100) is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: kde1D\nvariable: 'x'\ndatafile: 'x'\nlevel: 150")

        # a non-positive bandwidth factor is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: kde1D\nvariable: 'x'\ndatafile: 'x'\nbandwidth: -1.0")

    def test_prepare_errors(self):

        ctx = _source_context()

        item = eos.figure.ItemFactory.from_yaml("type: kde1D\nvariable: 'x'\ndatafile: '" + _ITEM_TEST_D + "'")
        with self.assertRaises(NotImplementedError):
            item.prepare(context=ctx)

        item = eos.figure.ItemFactory.from_yaml("type: kde1D\nvariable: 'nope'\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

    def test_prediction(self):

        ctx = _source_context()

        # the stripped-name lookup (name without options/kinematics)
        item = eos.figure.ItemFactory.from_yaml("type: kde1D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variable: 'B->D::f_+(q2)'")
        item.prepare(context=ctx)
        _, ax = plt.subplots()
        item.draw(ax)

        # the full-name lookup
        item = eos.figure.ItemFactory.from_yaml("type: kde1D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variable: '" + _PRED_FP + "'")
        item.prepare(context=ctx)

        # a stripped name that is ambiguous across several predictions is rejected
        item = eos.figure.ItemFactory.from_yaml("type: kde1D\ndatafile: '" + _ITEM_TEST_D + "/pred-dup'\n"
                                                "variable: 'B->D::f_+(q2)'")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

class TwoDimensionalKernelDensityItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        with self.assertRaises(ValueError):  # credibility level outside (0, 100)
            eos.figure.ItemFactory.from_yaml("type: kde2D\nvariables: ['a', 'b']\ndatafile: 'x'\nlevels: [150]")

        with self.assertRaises(ValueError):  # non-positive bandwidth factor
            eos.figure.ItemFactory.from_yaml("type: kde2D\nvariables: ['a', 'b']\ndatafile: 'x'\nbandwidth: 0.0")

        with self.assertRaises(ValueError):  # unsupported contour type
            eos.figure.ItemFactory.from_yaml("type: kde2D\nvariables: ['a', 'b']\ndatafile: 'x'\ncontours: ['nonsense']")

    def test_prepare_errors(self):

        ctx = _source_context()

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\nvariables: ['a', 'b']\ndatafile: '" + _ITEM_TEST_D + "'")
        with self.assertRaises(NotImplementedError):
            item.prepare(context=ctx)

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\nvariables: ['nope', '" + _PRED_F0 + "']\n"
                                                "datafile: '" + _ITEM_TEST_D + "/pred-uniq'")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

    def test_prediction(self):

        ctx = _source_context()

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['B->D::f_+(q2)', 'B->D::f_0(q2)']")
        item.prepare(context=ctx)
        _, ax = plt.subplots()
        item.draw(ax)

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['" + _PRED_FP + "', '" + _PRED_F0 + "']")
        item.prepare(context=ctx)

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: '" + _ITEM_TEST_D + "/pred-dup'\n"
                                                "variables: ['B->D::f_+(q2)', 'B->D::f_+(q2)']")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

class TwoDimensionalContoursItemPredictionTests(unittest.TestCase):

    def test_prepare_errors(self):

        ctx = _source_context()

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nvariables: ['a', 'b']\ndatafile: '" + _ITEM_TEST_D + "'")
        with self.assertRaises(NotImplementedError):
            item.prepare(context=ctx)

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nvariables: ['nope', '" + _PRED_F0 + "']\n"
                                                "datafile: '" + _ITEM_TEST_D + "/pred-uniq'")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

    def test_prediction(self):

        ctx = _source_context()

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['B->D::f_+(q2)', 'B->D::f_0(q2)']")
        item.prepare(context=ctx)
        _, ax = plt.subplots()
        item.draw(ax)

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['" + _PRED_FP + "', '" + _PRED_F0 + "']")
        item.prepare(context=ctx)

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: '" + _ITEM_TEST_D + "/pred-dup'\n"
                                                "variables: ['B->D::f_+(q2)', 'B->D::f_+(q2)']")
        with self.assertRaises(ValueError):
            item.prepare(context=ctx)

class ConstraintItemCoverageTests(unittest.TestCase):

    def test_gaussian(self):

        # a univariate Gaussian constraint is prepared into a single (x, y) point and drawn as one
        # error bar with an asymmetric y error. Regression test: drawing such a constraint previously
        # failed because the asymmetric error (a pair) was passed to Axes.errorbar with a shape that
        # matplotlib interpreted as two data points.
        item = eos.figure.ItemFactory.from_yaml("type: constraint\nconstraints: 'B->D::f_+@FKKM:2008A'\nvariable: q2")
        item.prepare()
        self.assertEqual(len(item._yvalues), 1)
        _, ax = plt.subplots()
        item.draw(ax)
        # exactly one error bar (one ErrorbarContainer) is drawn
        self.assertEqual(len(ax.containers), 1)

    def test_multivariate_gaussian(self):

        # a MultivariateGaussian constraint (with a non-matching observable per entry that is skipped)
        item = eos.figure.ItemFactory.from_yaml("type: constraint\nconstraints: 'B->D::f_++f_0@HPQCD:2015A'\n"
                                                "observable: 'B->D::f_+(q2)'\nvariable: q2")
        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)

    def test_invalid(self):

        # 'constraints' of an unsupported type is rejected at construction
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: constraint\nconstraints: 123\nvariable: q2")

class ConstraintResidueItemCoverageTests(unittest.TestCase):

    def test_gaussian(self):

        # a univariate Gaussian constraint residue is drawn as one error bar with an asymmetric y
        # error (see the regression note in ConstraintItemCoverageTests.test_gaussian)
        item = eos.figure.ItemFactory.from_yaml("type: constraint-residue\nconstraints: 'B->D::f_+@FKKM:2008A'\n"
                                                "observable: 'B->D::f_+(q2)'\nvariable: q2")
        item.prepare()
        self.assertEqual(len(item._yvalues), 1)
        _, ax = plt.subplots()
        item.draw(ax)
        self.assertEqual(len(ax.containers), 1)

    def test_multivariate_gaussian(self):

        item = eos.figure.ItemFactory.from_yaml("type: constraint-residue\nconstraints: 'B->D::f_++f_0@HPQCD:2015A'\n"
                                                "observable: 'B->D::f_+(q2)'\nvariable: q2")
        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)

    def test_invalid(self):

        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: constraint-residue\nconstraints: 123\n"
                                             "observable: 'B->D::f_+(q2)'\nvariable: q2")

class TwoDimensionalConstraintItemCoverageTests(unittest.TestCase):

    def test_invalid_type(self):

        # a 'constraint' that is neither a string nor a QualifiedName is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: constraint2D\nconstraint: 123\nx: { observable: 'a' }")

    def test_prepare_errors(self):

        # a constraint of an unsupported type (here: Amoroso) is rejected during prepare()
        item = eos.figure.ItemFactory.from_yaml("type: constraint2D\n"
                                                "constraint: 'B^0_s->mu^+mu^-::BR@CMS+LHCb:2014A'\n"
                                                "x: { observable: 'B_q->ll::BR@Untagged' }")
        with self.assertRaises(ValueError):
            item.prepare()

class BandItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        with self.assertRaises(ValueError):  # neither 'x' nor 'y'
            eos.figure.ItemFactory.from_yaml("type: band")

        with self.assertRaises(ValueError):  # 'x' not a pair
            eos.figure.ItemFactory.from_yaml("type: band\nx: [1.0]")

        with self.assertRaises(ValueError):  # inverted 'x'
            eos.figure.ItemFactory.from_yaml("type: band\nx: [2.0, 1.0]")

        with self.assertRaises(ValueError):  # 'y' not a pair
            eos.figure.ItemFactory.from_yaml("type: band\ny: [1.0]")

        with self.assertRaises(ValueError):  # inverted 'y'
            eos.figure.ItemFactory.from_yaml("type: band\ny: [2.0, 1.0]")

class VerticalLineItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        # 'x' is mandatory
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: vertical")

class SignalPDFItemCoverageTests(unittest.TestCase):

    def test_invalid(self):

        base = "type: signal-pdf\npdf: 'B->Dlnu::P(q2);l=mu'\nvariable: q2\n"

        with self.assertRaises(ValueError):  # null range
            eos.figure.ItemFactory.from_yaml(base + "range: null")

        with self.assertRaises(ValueError):  # inverted range
            eos.figure.ItemFactory.from_yaml(base + "range: [2.0, 1.0]")

        with self.assertRaises(ValueError):  # non-positive resolution
            eos.figure.ItemFactory.from_yaml(base + "range: [0.02, 11.6]\nresolution: 0")

        # a null kinematic variable is rejected in prepare()
        item = eos.figure.ItemFactory.from_yaml("type: signal-pdf\npdf: 'B->Dlnu::P(q2);l=mu'\n"
                                                "variable: null\nrange: [0.02, 11.6]")
        with self.assertRaises(ValueError):
            item.prepare()

    def test_options_and_parameters(self):

        # explicit options and parameters are applied during prepare()
        item = eos.figure.ItemFactory.from_yaml("""
        type: signal-pdf
        pdf: 'B->Dlnu::P(q2);l=mu'
        variable: q2
        range: [0.02, 11.6]
        resolution: 10
        kinematics: { q2_min: 0.02, q2_max: 11.6 }
        options: { l: 'mu' }
        parameters: { 'mass::mu': 0.10566 }
        """)
        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)

class ComplexPlaneItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        base = "type: complex-plane\nobservable: 'b->s::Re{F17}(Re{q2},Im{q2})'\n"

        # the first variable is not a kinematic variable
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml(base + "variables: ['not_a_var', 'Im{q2}']\nranges: [[-1, 1], [-1, 1]]")

        # the second variable is not a kinematic variable
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml(base + "variables: ['Re{q2}', 'not_a_var']\nranges: [[-1, 1], [-1, 1]]")

        # a fixed kinematic variable that the observable does not declare is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml(base + "variables: ['Re{q2}', 'Im{q2}']\nranges: [[-1, 1], [-1, 1]]\n"
                                             "fixed_kinematics: { bad_kin: 1.0 }")

        # an inverted x range is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml(base + "variables: ['Re{q2}', 'Im{q2}']\nranges: [[1, -1], [-1, 1]]")

        # an inverted y range is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml(base + "variables: ['Re{q2}', 'Im{q2}']\nranges: [[-1, 1], [1, -1]]")

        # a non-positive resolution is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml(base + "variables: ['Re{q2}', 'Im{q2}']\nranges: [[-1, 1], [-1, 1]]\nresolution: 0")

    def test_missing_fixed_parameters_from_file(self):

        # Regression (issue #1181): loading must not require the parameter file to exist; it is read
        # only in prepare(), which raises if the file is genuinely missing at draw time.
        item = eos.figure.ItemFactory.from_yaml("""
        type: complex-plane
        observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
        variables: ['Re{q2}', 'Im{q2}']
        ranges: [[-1.0, 1.0], [-1.0, 1.0]]
        resolution: 5
        fixed_parameters_from_file: 'does-not-exist.yaml'
        """)
        with self.assertRaises(RuntimeError):
            item.prepare()

class ErrorBarsItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        with self.assertRaises(ValueError):  # no positions
            eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: []\nyerrors: []")

        with self.assertRaises(ValueError):  # neither x nor y errors
            eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]")

        with self.assertRaises(ValueError):  # x error count mismatch
            eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]\nxerrors: [0.1, 0.2]")

        with self.assertRaises(ValueError):  # y error count mismatch
            eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]\nyerrors: [0.1, 0.2]")

    def test_invalid_error_specs(self):

        # an x error tuple of the wrong length is rejected in prepare()
        item = eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]\nxerrors: [[0.1, 0.2, 0.3]]")
        with self.assertRaises(ValueError):
            item.prepare()

        # an x error of an unsupported type is rejected in prepare()
        item = eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]\nxerrors: ['bad']")
        with self.assertRaises(ValueError):
            item.prepare()

        # a y error tuple of the wrong length is rejected in prepare()
        item = eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]\nyerrors: [[0.1, 0.2, 0.3]]")
        with self.assertRaises(ValueError):
            item.prepare()

        # a y error of an unsupported type is rejected in prepare()
        item = eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2]]\nyerrors: ['bad']")
        with self.assertRaises(ValueError):
            item.prepare()

    def test_xerrors_only(self):

        # error bars with only x errors leave the y errors unset
        item = eos.figure.ItemFactory.from_yaml("type: errorbars\npositions: [[1, 2], [2, 3]]\nxerrors: [0.1, 0.2]")
        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)
        self.assertIsNone(item._yerr)

class PointItemValidationTests(unittest.TestCase):

    def test_invalid(self):

        # an explicitly null 'x' is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: point\nx: null\ny: 0.0")

        # an explicitly null 'y' is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: point\nx: 0.0\ny: null")

if __name__ == '__main__':
    unittest.main(verbosity=5)
