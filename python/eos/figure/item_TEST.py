# Copyright (c) 2025-2026 Danny van Dyk
# Copyright (c) 2026      Dominik Suelmann
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

from unittest import mock

import eos
import eos.data
import eos.figure
import math
import numpy as np
import os
import shutil
import tempfile

from eos.analysis_file_context import AnalysisFileContext
from eos.figure.item import BandHandle, BandHandler, CompositeRegionHandle, CompositeRegionHandler, ConstraintItem, ConstraintResidueItem, Item
from eos.validation_context import ValidationContext
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from matplotlib import transforms as mtransforms
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

    def test_validate_semantics_checks_fixed_parameters(self):
        item = eos.figure.ItemFactory.from_yaml("""
        type: observable
        observable: 'B->Dlnu::dBR/dq2'
        variable: q2
        range: [0.1, 1.0]
        fixed_parameters: { 'test::unknown-parameter': 1.0 }
        """)
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()
        diagnostics = list(item.validate_semantics(ValidationContext(description)))

        self.assertEqual(1, len(diagnostics))
        self.assertEqual(
            ('fixed_parameters', 'test::unknown-parameter'),
            diagnostics[0].path,
        )

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

    _DATAFILE = 'eos/data/prediction_TEST.d/predictions'

    def prepared(self, item):
        "Prepares an item against the checked-in prediction fixture."
        item.prepare(context=AnalysisFileContext(base_directory=os.environ['SOURCE_DIR']))
        return item

    def datafile(self, *observables):
        "Writes a prediction file holding the given observables into a temporary directory."
        base = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, base, ignore_errors=True)
        eos.data.Prediction.create(os.path.join(base, 'predictions'), observables,
                                   np.random.RandomState(1).rand(6, len(observables)), np.ones(6))
        return base

    @staticmethod
    def observable(name, **kinematics):
        return eos.Observable.make(name, eos.Parameters.Defaults(), eos.Kinematics(**kinematics), eos.Options())

    def test_no_predictions(self):

        base = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, base, ignore_errors=True)
        eos.data.Prediction.create(os.path.join(base, 'predictions'), [], np.zeros((6, 0)), np.ones(6))

        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'predictions'\nvariable: 'q2'")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=AnalysisFileContext(base_directory=base))

        self.assertIn('does not contain any predictions', str(cm.exception))

    def test_variable_is_inferred(self):

        # without 'variable', the sole kinematic variable of the first prediction is adopted
        item = self.prepared(eos.figure.ItemFactory.from_yaml(f"type: uncertainty\ndatafile: '{self._DATAFILE}'"))

        self.assertEqual(item._variable, 'q2')
        self.assertEqual(str(item._observable), 'B->pilnu::dBR/dq2')

    def test_ambiguous_kinematic_variable(self):

        # a prediction over a bin has two kinematic variables, so 'variable' must be given
        base = self.datafile(*[self.observable('B->pilnu::BR', q2_min=float(i), q2_max=float(i + 1)) for i in range(3)])

        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'predictions'")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=AnalysisFileContext(base_directory=base))

        self.assertIn('contains more than one kinematic variable', str(cm.exception))

    def test_ambiguous_observable(self):

        # with more than one predicted observable, 'observable' must be given
        base = self.datafile(self.observable('B->pilnu::dBR/dq2', q2=1.0), self.observable('B->Dlnu::dBR/dq2', q2=1.0))

        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'predictions'\nvariable: 'q2'")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=AnalysisFileContext(base_directory=base))

        self.assertIn('contains more than one predicted observable', str(cm.exception))

    def test_variable_not_predicted(self):

        item = eos.figure.ItemFactory.from_yaml(f"type: uncertainty\ndatafile: '{self._DATAFILE}'\nvariable: 'k2'")
        with self.assertRaises(ValueError) as cm:
            self.prepared(item)

        self.assertIn("does not depend on the chosen kinematic variable 'k2'", str(cm.exception))

    def test_observable_options_are_matched(self):

        # the fixture was predicted with 'form-factors=BCL2008', which the item may state explicitly
        item = self.prepared(eos.figure.ItemFactory.from_yaml(
            f"type: uncertainty\ndatafile: '{self._DATAFILE}'\nvariable: 'q2'\n"
            "observable: 'B->pilnu::dBR/dq2;form-factors=BCL2008'"))

        self.assertEqual(len(item._xvalues), item.resolution)

    def test_observable_options_mismatch(self):

        # an option that no prediction of the data file was computed with leaves nothing to plot
        item = eos.figure.ItemFactory.from_yaml(
            f"type: uncertainty\ndatafile: '{self._DATAFILE}'\nvariable: 'q2'\n"
            "observable: 'B->pilnu::dBR/dq2;form-factors=BSZ2015'")
        with self.assertRaises(ValueError) as cm:
            self.prepared(item)

        self.assertIn('No observable of the data file matches', str(cm.exception))

    def test_cubic_interpolation(self):

        # both interpolations yield a band over the same abscissae, but not the same band
        linear = self.prepared(eos.figure.ItemFactory.from_yaml(
            f"type: uncertainty\ndatafile: '{self._DATAFILE}'\nvariable: 'q2'\ninterpolation: 'linear'"))
        cubic = self.prepared(eos.figure.ItemFactory.from_yaml(
            f"type: uncertainty\ndatafile: '{self._DATAFILE}'\nvariable: 'q2'\ninterpolation: 'cubic'"))

        self.assertEqual(list(cubic._xvalues), list(linear._xvalues))
        self.assertEqual(len(cubic._ovalues_central[0]), len(linear._ovalues_central[0]))
        self.assertFalse(np.allclose(cubic._ovalues_central[0], linear._ovalues_central[0]))

        _, ax = plt.subplots()
        cubic.draw(ax)

    def test_full(self):

        try:
            input = """
            type: uncertainty
            label: '$\\ell=\\mu$'
            variable: 'q2'
            levels: [95]
            range: [0.02, 11.63]
            datafile: 'eos/data/prediction_TEST.d/predictions'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

        try:
            input = """
            type: uncertainty
            label: '$\\ell=\\mu$'
            variable: 'q2'
            levels: [95]
            range: [0.02, 11.63]
            band: ['area']
            datafile: 'eos/data/prediction_TEST.d/predictions'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

    def test_legend(self):

        # the key of a labelled band is a single entry that reflects the 'band' option; the item
        # need not be prepared, since the key is determined by the item's options alone
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlabel: 'band'")
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], BandHandle)
        self.assertEqual(entries[0][1], 'band')

        # by default every part of the band is drawn, so the key shows a swatch for the single
        # credibility level, the boundary formed by the outer lines, and the median line
        handle = entries[0][0]
        self.assertEqual(len(handle.facecolors), 1)
        self.assertTrue(handle.boundary)
        self.assertTrue(handle.median)

        # the filled area carries one swatch per drawn region, with the shades of the drawn regions,
        # i.e. with opacity decreasing from the innermost region at the left to the outermost at the
        # right; see CompositeRegionHandle
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlabel: 'band'\nlevels: [68.27, 95.45]")
        handle = item.legend()[0][0]
        self.assertEqual(len(handle.facecolors), 2)
        opacities = [mcolors.to_rgba(facecolor)[3] for facecolor in handle.facecolors]
        self.assertEqual(opacities, sorted(opacities, reverse=True))

        # a band drawn as a median alone is keyed by neither a swatch nor a boundary ...
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlabel: 'band'\nband: 'median'\n"
                                                "color: 'C0'\nlinestyle: 'dashed'")
        handle = item.legend()[0][0]
        self.assertEqual(handle.facecolors, [])
        self.assertFalse(handle.boundary)
        self.assertTrue(handle.median)

        # ... but by a single line across the middle of the key, spanning its full width, i.e. the
        # key an ordinary line handle would produce. The handler does not consult the legend it is
        # passed, so None suffices here.
        artists = BandHandler().create_artists(None, handle, 0.0, 0.0, 30.0, 10.0, 10.0,
                                              mtransforms.IdentityTransform())
        self.assertEqual(len(artists), 1)
        self.assertIsInstance(artists[0], Line2D)
        self.assertEqual(list(artists[0].get_xdata()), [0.0, 30.0])
        self.assertEqual(list(artists[0].get_ydata()), [5.0, 5.0])
        self.assertEqual(mcolors.to_rgba(artists[0].get_color()), mcolors.to_rgba(item.color))
        self.assertEqual(artists[0].get_linestyle(), '--')

        # the outer lines alone are keyed by an unfilled boundary and nothing else
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlabel: 'band'\nband: 'outer'")
        handle = item.legend()[0][0]
        self.assertEqual(handle.facecolors, [])
        self.assertTrue(handle.boundary)
        self.assertFalse(handle.median)
        artists = BandHandler().create_artists(None, handle, 0.0, 0.0, 30.0, 10.0, 10.0,
                                              mtransforms.IdentityTransform())
        self.assertEqual(len(artists), 1)
        self.assertIsInstance(artists[0], Rectangle)
        self.assertFalse(artists[0].get_fill())

        # the filled area alone is keyed by swatches, with neither a boundary nor a median line
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlabel: 'band'\nband: ['area']\n"
                                                "levels: [68.27, 95.45]")
        handle = item.legend()[0][0]
        self.assertEqual(len(handle.facecolors), 2)
        self.assertFalse(handle.boundary)
        self.assertFalse(handle.median)
        artists = BandHandler().create_artists(None, handle, 0.0, 0.0, 30.0, 10.0, 10.0,
                                              mtransforms.IdentityTransform())
        self.assertEqual(len(artists), 2)
        for swatch in artists:
            self.assertTrue(swatch.get_fill())
            self.assertEqual(swatch.get_linewidth(), 0.0)

        # with every part drawn, the key shows the swatches, a boundary around them, and the median
        # line across their middle, all in the item's line style
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlabel: 'band'\n"
                                                "band: ['area', 'outer', 'median']\nlevels: [68.27, 95.45]\n"
                                                "linestyle: 'dashed'")
        handle = item.legend()[0][0]
        artists = BandHandler().create_artists(None, handle, 0.0, 0.0, 30.0, 10.0, 10.0,
                                               mtransforms.IdentityTransform())
        self.assertEqual(len(artists), len(handle.facecolors) + 2)
        swatches, boundary, median = artists[:-2], artists[-2], artists[-1]
        # the swatches abut and together span exactly the width allotted to the key
        self.assertAlmostEqual(sum(swatch.get_width() for swatch in swatches), 30.0, places=12)
        for lower, upper in zip(swatches[:-1], swatches[1:]):
            self.assertAlmostEqual(lower.get_x() + lower.get_width(), upper.get_x(), places=12)
        # a single unfilled boundary spans the whole key
        self.assertFalse(boundary.get_fill())
        self.assertAlmostEqual(boundary.get_width(), 30.0, places=12)
        self.assertEqual(boundary.get_linestyle(), 'dashed')
        # the median line is vertically centered within the key and spans its full width
        self.assertIsInstance(median, Line2D)
        self.assertEqual(list(median.get_xdata()), [0.0, 30.0])
        self.assertEqual(list(median.get_ydata()), [5.0, 5.0])
        self.assertEqual(median.get_linestyle(), '--')

        # a key that shows none of the band's parts is rejected rather than rendered as an empty key
        with self.assertRaises(ValueError):
            BandHandle()

        # an unlabelled band contributes no entry
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'")
        self.assertEqual(list(item.legend()), [])

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

        try:
            input = """
            type: uncertainty-binned
            label: '$\\ell=\\mu$'
            variable: 'q2'
            range: [0.02, 11.63]
            band: ['area']
            datafile: 'eos/data/prediction_TEST.d/predictions-binned'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare(context=AnalysisFileContext(base_directory=os.path.join(os.environ['SOURCE_DIR'])))
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'observable': {e}")

    def test_legend(self):

        # as for an unbinned band, the key of a labelled item reflects the 'band' option; the item
        # need not be prepared, since the key is determined by the item's options alone
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'q2'\ndatafile: 'x'\nlabel: 'band'")
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], BandHandle)
        self.assertEqual(entries[0][1], 'band')

        # by default every part of the band is drawn, so the key shows a swatch for the single
        # credibility level, the boundary formed by the outer lines, and the median line
        handle = entries[0][0]
        self.assertEqual(len(handle.facecolors), 1)
        self.assertTrue(handle.boundary)
        self.assertTrue(handle.median)

        # the filled area carries one swatch per drawn region
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'q2'\ndatafile: 'x'\nlabel: 'band'\n"
                                                "levels: [68.27, 95.45]")
        self.assertEqual(len(item.legend()[0][0].facecolors), 2)

        # a band drawn as a median alone is keyed by neither a swatch nor a boundary, i.e. by the
        # single line that an ordinary line handle would produce
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'q2'\ndatafile: 'x'\nlabel: 'band'\n"
                                                "band: 'median'")
        handle = item.legend()[0][0]
        self.assertEqual(handle.facecolors, [])
        self.assertFalse(handle.boundary)
        self.assertTrue(handle.median)

        # the outer lines alone are keyed by an unfilled boundary
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'q2'\ndatafile: 'x'\nlabel: 'band'\n"
                                                "band: 'outer'")
        handle = item.legend()[0][0]
        self.assertEqual(handle.facecolors, [])
        self.assertTrue(handle.boundary)
        self.assertFalse(handle.median)

        # the filled area alone is keyed by swatches, with neither a boundary nor a median line
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'q2'\ndatafile: 'x'\nlabel: 'band'\n"
                                                "band: ['area']")
        handle = item.legend()[0][0]
        self.assertEqual(len(handle.facecolors), 1)
        self.assertFalse(handle.boundary)
        self.assertFalse(handle.median)

        # an unlabelled band contributes no entry
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'q2'\ndatafile: 'x'")
        self.assertEqual(list(item.legend()), [])

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

    def test_binned_gaussian(self):

        # a univariate constraint binned in the variable is placed at the centre of its bin, with
        # the half bin width as the x error
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar:2010A'
        observable: 'B->pilnu::BR'
        variable: 'q2'
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0])
        self.assertEqual(list(item._xerrors), [2.0])
        self.assertAlmostEqual(item._yvalues[0], 3.13e-05)
        # the statistical and systematic uncertainties are added in quadrature, per side
        sigma = math.sqrt(3.0e-06**2 + 2.5e-06**2)
        self.assertAlmostEqual(item._yerrors[0][0], sigma)
        self.assertAlmostEqual(item._yerrors[0][1], sigma)

        _, ax = plt.subplots()
        item.draw(ax)

    def test_binned_gaussian_rescaled_by_width(self):

        # 'rescale_by_width' divides the central value and its uncertainties by the bin width,
        # while leaving the abscissa and its error untouched
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar:2010A'
        observable: 'B->pilnu::BR'
        variable: 'q2'
        rescale_by_width: true
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0])
        self.assertEqual(list(item._xerrors), [2.0])
        self.assertAlmostEqual(item._yvalues[0], 3.13e-05 / 4.0)
        sigma = math.sqrt(3.0e-06**2 + 2.5e-06**2) / 4.0
        self.assertAlmostEqual(item._yerrors[0][0], sigma)
        self.assertAlmostEqual(item._yerrors[0][1], sigma)

    def test_binned_multivariate_gaussian(self):

        # a multivariate constraint contributes one point per bin, in the order of its entries
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR'
        variable: 'q2'
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0, 6.0, 9.0, 11.0, 16.5])
        self.assertEqual(list(item._xerrors), [2.0, 2.0, 1.0, 1.0, 4.5])
        self.assertEqual(list(item._yvalues), [1.506e-05, 1.815e-05, 1.471e-05, 1.715e-05, 4.981e-05])
        self.assertTrue(all(item._yerrors > 0.0))

        _, ax = plt.subplots()
        item.draw(ax)

    def test_binned_multivariate_gaussian_rescaled_by_width(self):

        # each bin of a multivariate constraint is rescaled by its own width
        def prepared(rescale):
            item = eos.figure.ItemFactory.from_yaml(f"""
            type: constraint
            constraints: 'B^+->omegalnu::BR@BPR:2021A'
            observable: 'B->omegalnu::BR'
            variable: 'q2'
            rescale_by_width: {rescale}
            """)
            item.prepare()
            return item

        plain, rescaled = prepared('false'), prepared('true')
        widths = [4.0, 4.0, 2.0, 2.0, 9.0]

        self.assertEqual(list(rescaled._xvalues), list(plain._xvalues))
        self.assertEqual(list(rescaled._xerrors), list(plain._xerrors))
        for i, width in enumerate(widths):
            self.assertAlmostEqual(rescaled._yvalues[i], plain._yvalues[i] / width)
            self.assertAlmostEqual(rescaled._yerrors[i], plain._yerrors[i] / width)

    def test_binned_multivariate_gaussian_observable_with_options(self):

        # the observable's options are matched against those recorded in the constraint
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR;l=e'
        variable: 'q2'
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0, 6.0, 9.0, 11.0, 16.5])

    def test_unknown_constraint(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B->pi::NOSUCH@Nobody:2000A'
        observable: 'B->pilnu::BR'
        variable: 'q2'
        """)
        with self.assertRaises(ValueError) as cm:
            item.prepare()

        self.assertIn('unknown constraint B->pi::NOSUCH@Nobody:2000A', str(cm.exception))

    def test_unsupported_constraint_type(self):

        # only the Gaussian and the two multivariate Gaussian types can be drawn
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^0_s->mu^+mu^-::BR@CMS+LHCb:2014A'
        observable: 'B_q->ll::BR@Untagged'
        variable: 'q2'
        """)
        with self.assertRaises(ValueError) as cm:
            item.prepare()

        self.assertIn('constraint type Amoroso presently not supported', str(cm.exception))

    def test_range_limits_the_drawn_points(self):

        # 'range' masks the points outside it, exclusive at both ends
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR'
        variable: 'q2'
        range: [5.0, 12.0]
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0, 6.0, 9.0, 11.0, 16.5])
        self.assertEqual(list(item._mask), [False, True, True, True, False])

        # only the three unmasked points are drawn, one error bar each
        _, ax = plt.subplots()
        item.draw(ax)
        self.assertEqual(len(ax.containers), 3)

    def test_multivariate_gaussian_covariance_requires_an_observable(self):

        # a multivariate constraint bundles several observables, so one of them must be named
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
        variable: 'q2'
        """)
        with self.assertRaises(KeyError) as cm:
            item.prepare()

        self.assertIn('MultivariateGaussian(Covariance)', cm.exception.args[0])

    def test_multivariate_gaussian_requires_an_observable(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        variable: 'q2'
        """)
        with self.assertRaises(KeyError) as cm:
            item.prepare()

        self.assertIn('MultivariateGaussian', cm.exception.args[0])

    def test_point_kinematics_multivariate_gaussian(self):

        # a constraint given at points rather than in bins carries no x error, and only the
        # entries naming the requested observable contribute
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'
        observable: 'B->D::f_+(q2)'
        variable: 'q2'
        """)
        with mock.patch.object(eos, 'warn') as warn:
            item.prepare()

        self.assertEqual(list(item._xvalues), [-15.0, -10.0, -5.0, 0.0])
        self.assertTrue(all(xerr is None for xerr in item._xerrors))
        self.assertEqual(list(item._yvalues), [0.414161569, 0.474556651, 0.550406251, 0.649944038])

        # the constraint holds 33 entries, of which the 29 for other form factors are skipped
        self.assertEqual(warn.call_count, 29)

        _, ax = plt.subplots()
        item.draw(ax)

    def test_binned_multivariate_gaussian_option_mismatch(self):

        # an option that no entry of the constraint provides leaves nothing to draw, and each
        # skipped entry is reported with its name and options
        item = eos.figure.ItemFactory.from_yaml("""
        type: constraint
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR;l=tau'
        variable: 'q2'
        """)
        with mock.patch.object(eos, 'warn') as warn:
            item.prepare()

        self.assertEqual(len(item._xvalues), 0)
        self.assertEqual(warn.call_count, 5)
        self.assertIn('B->omegalnu::BR', warn.call_args[0][0])
        self.assertIn("'l': 'e'", warn.call_args[0][0])

        # nothing is drawn, and no exception is raised
        _, ax = plt.subplots()
        item.draw(ax)

class ConstraintResidueItemTests(unittest.TestCase):

    class _Parameter:
        "Minimal stand-in for eos.Parameter, exposing the name/min/max accessors that Mode.create() uses."

        def __init__(self, name, minimum, maximum):
            self._name = name
            self._min  = minimum
            self._max  = maximum

        def name(self):
            return self._name

        def min(self):
            return self._min

        def max(self):
            return self._max

    def test_full(self):

        # the default style is 'pull', drawing a bar per residue
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
            self.assertEqual(item.style, 'pull')
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
            self.assertGreater(len(ax.patches), 0)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint': {e}")

    def test_full_delta_style(self):

        # 'delta' falls back to the original error-bar drawing
        try:
            input = """
            type: 'constraint-residue'
            label: r'Belle 2015 $\\ell=e,\\, q=d$'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            parameters: {"mass::e": 1.0}
            rescale_by_width: true
            style: 'delta'
            """
            item = eos.figure.ItemFactory.from_yaml(input)
            item.prepare()
            _, ax = plt.subplots()
            item.draw(ax)
        except Exception as e:
            self.fail(f"Error when testing item of type 'constraint-residue' with style 'delta': {e}")

    def test_binned_gaussian(self):

        # a univariate constraint binned in the variable is placed at the centre of its bin, and
        # the observable is evaluated over that same bin
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar:2010A'
        observable: 'B->pilnu::BR'
        variable: 'q2'
        style: 'delta'
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0])
        self.assertEqual(list(item._xerrors), [2.0])
        self.assertEqual(len(item._yvalues), 1)
        self.assertTrue(math.isfinite(item._yvalues[0]))
        # the residue is the measurement minus the prediction, both integrated over the bin
        self.assertLess(item._yvalues[0], 3.13e-05)

        _, ax = plt.subplots()
        item.draw(ax)

    def test_binned_gaussian_rescaled_by_width(self):

        # 'rescale_by_width' divides the residue and its uncertainties by the bin width
        def prepared(rescale):
            item = eos.figure.ItemFactory.from_yaml(f"""
            type: 'constraint-residue'
            constraints: 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar:2010A'
            observable: 'B->pilnu::BR'
            variable: 'q2'
            rescale_by_width: {rescale}
            """)
            item.prepare()
            return item

        plain, rescaled = prepared('false'), prepared('true')

        self.assertEqual(list(rescaled._xvalues), list(plain._xvalues))
        self.assertEqual(list(rescaled._xerrors), list(plain._xerrors))
        self.assertAlmostEqual(rescaled._yvalues[0], plain._yvalues[0] / 4.0)
        self.assertAlmostEqual(rescaled._yerrors[0][0], plain._yerrors[0][0] / 4.0)
        self.assertAlmostEqual(rescaled._yerrors[0][1], plain._yerrors[0][1] / 4.0)

    def test_binned_multivariate_gaussian(self):

        # a multivariate constraint contributes one residue per bin, in the order of its entries
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR'
        variable: 'q2'
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0, 6.0, 9.0, 11.0, 16.5])
        self.assertEqual(list(item._xerrors), [2.0, 2.0, 1.0, 1.0, 4.5])
        self.assertEqual(len(item._yvalues), 5)
        self.assertTrue(all(math.isfinite(y) for y in item._yvalues))

        _, ax = plt.subplots()
        item.draw(ax)
        self.assertGreater(len(ax.patches), 0)

    def test_unknown_constraint(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B->pi::NOSUCH@Nobody:2000A'
        observable: 'B->pilnu::BR'
        variable: 'q2'
        """)
        with self.assertRaises(ValueError) as cm:
            item.prepare()

        self.assertIn('unknown constraint B->pi::NOSUCH@Nobody:2000A', str(cm.exception))

    def test_unsupported_constraint_type(self):

        # only the Gaussian and the two multivariate Gaussian types can be drawn
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^0_s->mu^+mu^-::BR@CMS+LHCb:2014A'
        observable: 'B_q->ll::BR@Untagged'
        variable: 'q2'
        """)
        with self.assertRaises(ValueError) as cm:
            item.prepare()

        self.assertIn('constraint type Amoroso presently not supported', str(cm.exception))

    def test_range_limits_the_drawn_residues(self):

        # 'range' masks the residues outside it, exclusive at both ends
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR'
        variable: 'q2'
        range: [5.0, 12.0]
        """)
        item.prepare()

        self.assertEqual(list(item._xvalues), [2.0, 6.0, 9.0, 11.0, 16.5])
        self.assertEqual(list(item._mask), [False, True, True, True, False])

        # only the three unmasked residues are drawn, one bar each
        _, ax = plt.subplots()
        item.draw(ax)
        self.assertEqual(len(ax.patches), 3)

    def test_nan_observable_gaussian(self):

        # an observable that cannot be evaluated yields a NaN residue, which is reported but
        # does not abort the preparation
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^0->pi^+lnu::BR[0.0,4.0]@BaBar:2010A'
        observable: 'B->pilnu::BR'
        variable: 'q2'
        parameters: {'CKM::lambda': .nan}
        """)
        with mock.patch.object(eos, 'warn') as warn:
            item.prepare()

        self.assertEqual(len(item._yvalues), 1)
        self.assertTrue(math.isnan(item._yvalues[0]))
        self.assertEqual(warn.call_count, 1)
        self.assertIn('B->pilnu::BR', warn.call_args[0][0])
        self.assertIn('evaluated to NaN', warn.call_args[0][0])

        # the residues are still drawn, as NaN-valued bars
        _, ax = plt.subplots()
        item.draw(ax)
        self.assertEqual(len(ax.patches), 1)

    def test_nan_observable_multivariate_gaussian_covariance(self):

        # every entry of a multivariate constraint is reported separately
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'
        observable: 'B->D::f_+(q2)'
        variable: 'q2'
        parameters: {'B->D::alpha^f+_0@BSZ2015': .nan}
        """)
        with mock.patch.object(eos, 'warn') as warn:
            item.prepare()

        self.assertEqual(len(item._yvalues), 4)
        self.assertTrue(all(math.isnan(yvalue) for yvalue in item._yvalues))
        self.assertEqual(len([call for call in warn.call_args_list if 'evaluated to NaN' in call[0][0]]), 4)

    def test_nan_observable_multivariate_gaussian(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        observable: 'B->omegalnu::BR'
        variable: 'q2'
        parameters: {'CKM::lambda': .nan}
        """)
        with mock.patch.object(eos, 'warn') as warn:
            item.prepare()

        self.assertEqual(len(item._yvalues), 5)
        self.assertTrue(all(math.isnan(yvalue) for yvalue in item._yvalues))
        self.assertEqual(len([call for call in warn.call_args_list if 'evaluated to NaN' in call[0][0]]), 5)

    def test_observable_not_contained_in_the_constraint(self):

        # a constraint that holds none of the requested observable is reported and skipped
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'
        observable: 'B->pi::f_+(q2)'
        variable: 'q2'
        """)
        with mock.patch.object(eos, 'warn'), mock.patch.object(eos, 'info') as info:
            item.prepare()

        self.assertEqual(len(item._xvalues), 0)
        self.assertIn('B->D^(*)::FormFactors', info.call_args[0][0])

        # nothing is drawn, and no exception is raised
        _, ax = plt.subplots()
        item.draw(ax)
        self.assertEqual(len(ax.patches), 0)

    def test_multivariate_gaussian_covariance_requires_an_observable(self):

        # a multivariate constraint bundles several observables, so one of them must be named
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
        variable: 'q2'
        """)
        with self.assertRaises(KeyError) as cm:
            item.prepare()

        self.assertIn('MultivariateGaussian(Covariance)', cm.exception.args[0])

    def test_multivariate_gaussian_requires_an_observable(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B^+->omegalnu::BR@BPR:2021A'
        variable: 'q2'
        """)
        with self.assertRaises(KeyError) as cm:
            item.prepare()

        self.assertIn('MultivariateGaussian', cm.exception.args[0])

    def test_point_kinematics_multivariate_gaussian(self):

        # a constraint given at points rather than in bins carries no x error, and the observable
        # is evaluated at each of those points
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'
        observable: 'B->D::f_+(q2)'
        variable: 'q2'
        """)
        with mock.patch.object(eos, 'warn') as warn:
            item.prepare()

        self.assertEqual(list(item._xvalues), [-15.0, -10.0, -5.0, 0.0])
        self.assertTrue(all(xerr is None for xerr in item._xerrors))

        # the constraint holds 33 entries, of which the 29 for other form factors are skipped
        self.assertEqual(warn.call_count, 29)

        # each residue is the constraint's mean less the form factor at that very point
        means = [0.414161569, 0.474556651, 0.550406251, 0.649944038]
        parameters = eos.Parameters.Defaults()
        options = eos.Options({'form-factors': 'BSZ2015'})
        for i, (q2, mean) in enumerate(zip([-15.0, -10.0, -5.0, 0.0], means)):
            prediction = eos.Observable.make('B->D::f_+(q2)', parameters, eos.Kinematics(q2=q2), options).evaluate()
            self.assertAlmostEqual(item._yvalues[i], mean - prediction)

        _, ax = plt.subplots()
        item.draw(ax)

    def test_point_kinematics_delta_style(self):

        # without an x error, the 'delta' style draws a plain vertical error bar
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        constraints: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'
        observable: 'B->D::f_+(q2)'
        variable: 'q2'
        style: 'delta'
        """)
        with mock.patch.object(eos, 'warn'):
            item.prepare()

        _, ax = plt.subplots()
        item.draw(ax)
        self.assertEqual(len(ax.containers), 4)
        self.assertTrue(all(not container.has_xerr for container in ax.containers))

    def test_binned_multivariate_gaussian_rescaled_by_width(self):

        # each bin of a multivariate constraint is rescaled by its own width
        def prepared(rescale):
            item = eos.figure.ItemFactory.from_yaml(f"""
            type: 'constraint-residue'
            constraints: 'B^+->omegalnu::BR@BPR:2021A'
            observable: 'B->omegalnu::BR'
            variable: 'q2'
            rescale_by_width: {rescale}
            """)
            item.prepare()
            return item

        plain, rescaled = prepared('false'), prepared('true')
        widths = [4.0, 4.0, 2.0, 2.0, 9.0]

        self.assertEqual(list(rescaled._xvalues), list(plain._xvalues))
        for i, width in enumerate(widths):
            self.assertAlmostEqual(rescaled._yvalues[i], plain._yvalues[i] / width)
            self.assertAlmostEqual(rescaled._yerrors[i], plain._yerrors[i] / width)

    def test_invalid_style(self):

        with self.assertRaises(Exception):
            eos.figure.ItemFactory.from_yaml("""
            type: 'constraint-residue'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            style: 'unknown'
            """)

    def test_pull(self):

        # a symmetric error is a plain division
        self.assertAlmostEqual(ConstraintResidueItem._pull(2.0, 1.0), 2.0)
        self.assertAlmostEqual(ConstraintResidueItem._pull(-2.0, 1.0), -2.0)

        # an asymmetric error uses the side towards zero: sigma_lo for a positive residue
        # (measured above the prediction), sigma_hi for a negative one
        self.assertAlmostEqual(ConstraintResidueItem._pull(2.0, (1.0, 0.5)), 4.0)
        self.assertAlmostEqual(ConstraintResidueItem._pull(-2.0, (1.0, 0.5)), -2.0)

        # a vanishing standard deviation yields NaN rather than raising ZeroDivisionError
        self.assertTrue(math.isnan(ConstraintResidueItem._pull(1.0, 0.0)))

    def test_legend(self):

        # a labelled residue item in the 'pull' style contributes a single filled-bar entry
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
        self.assertIsInstance(entries[0][0], Rectangle)
        self.assertEqual(entries[0][1], 'Belle')

    def test_legend_delta_style(self):

        # in the 'delta' style, a labelled residue item contributes a single error-bar entry
        # (a capped bar rendered by HandlerErrorbar), matching how the residues are drawn
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint-residue'
        label: 'Belle'
        constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
        observable: 'B->Dlnu::BR'
        variable: 'q2'
        rescale_by_width: true
        style: 'delta'
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

    def test_mode_file(self):

        # 'mode_file' loads the best-fit parameter values from a stored eos.data.Mode, deferred to
        # prepare() so that the file need not exist at load time
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mode-default')
            eos.data.Mode.create(path, [self._Parameter('mass::e', 0.0, 1.0)], [0.5], None, [], None, None)

            item = eos.figure.ItemFactory.from_yaml(f"""
            type: 'constraint-residue'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            rescale_by_width: true
            mode_file: '{path}'
            """)
            item.prepare()
            self.assertEqual(item._parameters['mass::e'].evaluate(), 0.5)

    def test_mode_file_relative_to_base_directory(self):

        # a relative 'mode_file' is resolved against the analysis file context's base_directory,
        # exactly like 'fixed_parameters_from_file' and 'parameters_from_file' elsewhere
        with tempfile.TemporaryDirectory() as d:
            eos.data.Mode.create(os.path.join(d, 'mode-default'), [self._Parameter('mass::e', 0.0, 1.0)], [0.5], None, [], None, None)

            item = eos.figure.ItemFactory.from_yaml("""
            type: 'constraint-residue'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            rescale_by_width: true
            mode_file: 'mode-default'
            """)
            item.prepare(context=AnalysisFileContext(base_directory=d))
            self.assertEqual(item._parameters['mass::e'].evaluate(), 0.5)

    def test_mode_file_overridden_by_parameters(self):

        # an explicit 'parameters' entry takes precedence over the value loaded from 'mode_file'
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'mode-default')
            eos.data.Mode.create(path, [self._Parameter('mass::e', 0.0, 1.0)], [0.5], None, [], None, None)

            item = eos.figure.ItemFactory.from_yaml(f"""
            type: 'constraint-residue'
            constraints: 'B^0->D^+e^-nu::BRs@Belle:2015A'
            observable: 'B->Dlnu::BR'
            variable: 'q2'
            rescale_by_width: true
            mode_file: '{path}'
            parameters: {{"mass::e": 1.0}}
            """)
            item.prepare()
            self.assertEqual(item._parameters['mass::e'].evaluate(), 1.0)

class TwoDimensionalConstraintItemTests(unittest.TestCase):

    def test_unknown_constraint(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B->pi::NOSUCH@Nobody:2000A'
        x: { observable: 'B->K^*gamma::S_K^*gamma' }
        y: { observable: 'B->K^*gamma::C_K^*gamma' }
        """)
        with self.assertRaises(ValueError) as cm:
            item.prepare()

        self.assertIn('unknown constraint B->pi::NOSUCH@Nobody:2000A', str(cm.exception))

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

        # a labelled 2D constraint contributes a single composite entry
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B^0->K^*0gamma::S_K+C_K@BaBar:2008A'
        x: { observable: 'B->K^*gamma::S_K^*gamma' }
        y: { observable: 'B->K^*gamma::C_K^*gamma' }
        label: 'BaBar 2008'
        """)
        entries = item.legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], CompositeRegionHandle)
        self.assertEqual(entries[0][1], 'BaBar 2008')

        # a single confidence level yields a single swatch, i.e. the composite key degenerates to
        # the one-shade key it replaces
        self.assertEqual(len(entries[0][0].facecolors), 1)

        # the key carries one swatch per drawn region ...
        item = eos.figure.ItemFactory.from_yaml("""
        type: 'constraint2D'
        constraint: 'B^0->K^*0gamma::S_K+C_K@BaBar:2008A'
        x: { observable: 'B->K^*gamma::S_K^*gamma' }
        y: { observable: 'B->K^*gamma::C_K^*gamma' }
        sigmas: [1.0, 2.0, 3.0]
        alpha: 0.5
        color: 'C0'
        linestyle: 'dashed'
        label: 'BaBar 2008'
        """)
        handle = item.legend()[0][0]
        self.assertEqual(len(handle.facecolors), 3)

        # ... whose shades are those of the drawn regions, i.e. each region's own opacity
        # composited onto the opacities of the regions enclosing it. The item draws its regions from
        # the outermost inwards, whereas the swatches run from the innermost region at the left to
        # the outermost at the right, so that opacity decreases from left to right.
        accumulated = []
        opacity = 0.0
        for alpha in item._alphas:
            opacity = opacity + alpha * (1.0 - opacity)
            accumulated.append(opacity)
        accumulated.reverse()
        self.assertEqual(accumulated, sorted(accumulated, reverse=True))
        for facecolor, opacity in zip(handle.facecolors, accumulated):
            self.assertAlmostEqual(mcolors.to_rgba(facecolor)[3], opacity, places=12)
            self.assertEqual(mcolors.to_rgba(facecolor)[:3], mcolors.to_rgba(item.color)[:3])

        # the boundary of the key reflects the line style of the outermost region
        self.assertEqual(handle.linestyle, 'dashed')
        self.assertEqual(handle.linewidth, item.linewidth)

        # the key is rendered as one swatch per region plus a single enclosing boundary. The handler
        # does not consult the legend it is passed, so None suffices here.
        artists = CompositeRegionHandler().create_artists(None, handle, 0.0, 0.0, 30.0, 10.0, 10.0,
                                                          mtransforms.IdentityTransform())
        self.assertEqual(len(artists), len(handle.facecolors) + 1)
        swatches, boundary = artists[:-1], artists[-1]
        # the swatches abut and together span exactly the width allotted to the key
        self.assertAlmostEqual(sum(s.get_width() for s in swatches), 30.0, places=12)
        for lower, upper in zip(swatches[:-1], swatches[1:]):
            self.assertAlmostEqual(lower.get_x() + lower.get_width(), upper.get_x(), places=12)
        # the swatches are not outlined, so that no divider appears between adjacent swatches
        for swatch in swatches:
            self.assertEqual(swatch.get_linewidth(), 0.0)
        # a single unfilled boundary spans the whole key
        self.assertFalse(boundary.get_fill())
        self.assertAlmostEqual(boundary.get_width(), 30.0, places=12)
        self.assertEqual(boundary.get_linestyle(), 'dashed')

        # an empty composite key is rejected rather than rendered as a boundary with no swatches
        with self.assertRaises(ValueError):
            CompositeRegionHandle([])

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

class ItemBaseClassTests(unittest.TestCase):

    def test_prepare_and_draw_are_abstract(self):

        item = Item(label='base')

        with self.assertRaises(NotImplementedError):
            item.prepare()

        with self.assertRaises(NotImplementedError):
            item.draw(None)

    def test_legend_is_empty_by_default(self):

        # an item type that draws nothing keyable contributes no entry, even when labelled
        self.assertEqual(list(Item(label='base').legend()), [])

    def test_legend_marker(self):

        entries = Item(label='base', color='C0')._legend_marker('o')
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][0].get_marker(), 'o')
        self.assertEqual(entries[0][1], 'base')

        # an unlabelled item contributes no entry
        self.assertEqual(list(Item()._legend_marker('o')), [])

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

    def test_fixed_parameters_take_precedence_over_the_file(self):

        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'params.yaml'), 'w') as f:
                f.write("'mass::mu':\n  central: 0.5\n")

            item = eos.figure.ItemFactory.from_yaml("""
            type: observable
            observable: 'B->Dlnu::dBR/dq2'
            variable: q2
            range: [0.1, 1.0]
            resolution: 3
            fixed_parameters_from_file: 'params.yaml'
            fixed_parameters: { 'mass::mu': 0.2 }
            """)
            with mock.patch.object(eos, 'warn') as warn:
                item.prepare(context=AnalysisFileContext(base_directory=tmp))

            self.assertEqual(item._parameters['mass::mu'].evaluate(), 0.2)
            self.assertEqual(len([call for call in warn.call_args_list if 'with explicit values' in call[0][0]]), 1)

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

        # an empty band set is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nband: []")

        # an unknown interpolation type is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\ninterpolation: 'quadratic'")

        # a level value larger than 100 is rejected
        with self.assertRaisesRegex(ValueError, "not in the interval"):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nlevels: [150]")

        # a range that is not a pair is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nrange: [1.0]")

        # an inverted range is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty\ndatafile: 'x'\nrange: [2.0, 1.0]")

class BinnedUncertaintyItemValidationTests(unittest.TestCase):

    def test_band_normalization(self):

        # a single band type given as a string is normalized to a set
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'nonexistent'\n"
                                                "datafile: 'eos/data/prediction_TEST.d/predictions-binned'\nband: 'median'")
        self.assertEqual(item.band, {'median'})

        # a list of band types is normalized to a set
        item = eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'nonexistent'\n"
                                                "datafile: 'eos/data/prediction_TEST.d/predictions-binned'\nband: ['area', 'outer']")
        self.assertEqual(item.band, {'area', 'outer'})

    def test_invalid(self):

        # a 'band' that is neither a string, list, nor set is rejected
        with self.assertRaisesRegex(ValueError, "must be a string, list of strings"):
            eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'nonexistent'\n"
                                                "datafile: 'eos/data/prediction_TEST.d/predictions-binned'\nband: 5")

        # an unknown band type is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'nonexistent'\n"
                                                "datafile: 'eos/data/prediction_TEST.d/predictions-binned'\nband: 'nonsense'")

        # an empty band set is rejected
        with self.assertRaises(ValueError):
            eos.figure.ItemFactory.from_yaml("type: uncertainty-binned\nvariable: 'nonexistent'\n"
                                                "datafile: 'eos/data/prediction_TEST.d/predictions-binned'\nband: []")

    def test_invalid_credibility_level(self):

        base = "type: uncertainty-binned\nvariable: 'x'\ndatafile: 'x'\n"

        with self.assertRaisesRegex(ValueError, "not in the interval"):
            eos.figure.ItemFactory.from_yaml(base + "levels: [150]")

        with self.assertRaisesRegex(ValueError, "not in the interval"):
            eos.figure.ItemFactory.from_yaml(base + "levels: [0]")

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

    def test_legend(self):

        item = eos.figure.ItemFactory.from_yaml("type: histogram1D\nvariable: 'x'\ndatafile: 'x'\nlabel: 'hist'")
        entries = item.legend()

        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Rectangle)
        self.assertEqual(entries[0][1], 'hist')

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

    def test_legend(self):

        item = eos.figure.ItemFactory.from_yaml("type: kde1D\nvariable: 'x'\ndatafile: 'x'\nlabel: 'kde'")
        entries = item.legend()

        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][1], 'kde')

    def test_samples_missing_variable(self):

        item = eos.figure.ItemFactory.from_yaml("type: kde1D\nvariable: 'not-a-variable'\n"
                                                "datafile: 'eos/data/importance_samples_TEST.d/samples'")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain samples of variable 'not-a-variable'", str(cm.exception))

class TwoDimensionalKernelDensityItemValidationTests(unittest.TestCase):

    def test_samples_missing_x_variable(self):

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: 'eos/data/importance_samples_TEST.d/samples'\n"
                                                "variables: ['not-an-x-variable', 'B->pi::f_+(0)@BCL2008']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain samples of variable 'not-an-x-variable'", str(cm.exception))

    def test_samples_missing_y_variable(self):

        # the x variable is present, so only the second guard can reject this
        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: 'eos/data/importance_samples_TEST.d/samples'\n"
                                                "variables: ['CKM::abs(V_ub)', 'not-a-y-variable']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain samples of variable 'not-a-y-variable'", str(cm.exception))

    def test_prediction_ambiguous_y_variable(self):

        # the x variable is given by its full name, the y variable by one 'pred-dup' holds twice
        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: '" + _ITEM_TEST_D + "/pred-dup'\n"
                                                "variables: ['" + _PRED_FP + "', 'B->D::f_+(q2)']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("contains multiple predictions for variable 'B->D::f_+(q2)'", str(cm.exception))

    def test_prediction_missing_y_variable(self):

        item = eos.figure.ItemFactory.from_yaml("type: kde2D\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['" + _PRED_FP + "', 'not-a-y-variable']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain predictions for variable 'not-a-y-variable'", str(cm.exception))

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

    def test_legend(self):

        base = "type: kde2D\nvariables: ['a', 'b']\ndatafile: 'x'\n"

        # filled contours are keyed by a swatch, unfilled ones by a line
        entries = eos.figure.ItemFactory.from_yaml(base + "label: 'kde'\ncontours: ['areas']").legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Rectangle)
        self.assertEqual(entries[0][1], 'kde')

        entries = eos.figure.ItemFactory.from_yaml(base + "label: 'kde'\ncontours: ['lines']").legend()
        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)

        # an unlabelled item contributes no entry
        self.assertEqual(list(eos.figure.ItemFactory.from_yaml(base + "contours: ['areas']").legend()), [])

class TwoDimensionalContoursItemPredictionTests(unittest.TestCase):

    def test_samples_missing_x_variable(self):

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: 'eos/data/importance_samples_TEST.d/samples'\n"
                                                "variables: ['not-an-x-variable', 'B->pi::f_+(0)@BCL2008']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain samples of variable 'not-an-x-variable'", str(cm.exception))

    def test_samples_missing_y_variable(self):

        # the x variable is present, so only the second guard can reject this
        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: 'eos/data/importance_samples_TEST.d/samples'\n"
                                                "variables: ['CKM::abs(V_ub)', 'not-a-y-variable']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain samples of variable 'not-a-y-variable'", str(cm.exception))

    def test_prediction_ambiguous_y_variable(self):

        # the x variable is given by its full name, the y variable by one 'pred-dup' holds twice
        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: '" + _ITEM_TEST_D + "/pred-dup'\n"
                                                "variables: ['" + _PRED_FP + "', 'B->D::f_+(q2)']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("contains multiple predictions for variable 'B->D::f_+(q2)'", str(cm.exception))

    def test_prediction_missing_y_variable(self):

        item = eos.figure.ItemFactory.from_yaml("type: contours2D\nbins: 20\ndatafile: '" + _ITEM_TEST_D + "/pred-uniq'\n"
                                                "variables: ['" + _PRED_FP + "', 'not-a-y-variable']")
        with self.assertRaises(ValueError) as cm:
            item.prepare(context=_source_context())

        self.assertIn("does not contain predictions for variable 'not-a-y-variable'", str(cm.exception))

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

    def test_validate_semantics(self):

        item = eos.figure.ItemFactory.from_yaml("type: constraint\nconstraints: 'B->D::f_+@FKKM:2008A'\n"
                                                "variable: q2\nobservable: 'test::unknown-observable'")
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()
        diagnostics = list(item.validate_semantics(ValidationContext(description)))

        self.assertEqual(1, len(diagnostics))
        self.assertEqual(('observable',), diagnostics[0].path)

        # an item without an observable has nothing to check
        item = eos.figure.ItemFactory.from_yaml("type: constraint\nconstraints: 'B->D::f_+@FKKM:2008A'\nvariable: q2")
        self.assertEqual([], list(item.validate_semantics(ValidationContext(description))))

    def test_constraints_as_a_qualified_name(self):

        # a single QualifiedName, as passed by programmatic construction, is wrapped in a list
        item = ConstraintItem(constraints=eos.QualifiedName('B->D::f_+@FKKM:2008A'), variable='q2')

        self.assertEqual([str(constraint) for constraint in item.constraints], ['B->D::f_+@FKKM:2008A'])

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

    def test_validate_semantics(self):

        item = eos.figure.ItemFactory.from_yaml("type: constraint-residue\nconstraints: 'B->D::f_+@FKKM:2008A'\n"
                                                "variable: q2\nobservable: 'test::unknown-observable'")
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()
        diagnostics = list(item.validate_semantics(ValidationContext(description)))

        self.assertEqual(1, len(diagnostics))
        self.assertEqual(('observable',), diagnostics[0].path)

        # an item without an observable has nothing to check
        item = eos.figure.ItemFactory.from_yaml("type: constraint-residue\nconstraints: 'B->D::f_+@FKKM:2008A'\nvariable: q2")
        self.assertEqual([], list(item.validate_semantics(ValidationContext(description))))

    def test_constraints_as_a_qualified_name(self):

        # a single QualifiedName, as passed by programmatic construction, is wrapped in a list
        item = ConstraintResidueItem(constraints=eos.QualifiedName('B->D::f_+@FKKM:2008A'), variable='q2')

        self.assertEqual([str(constraint) for constraint in item.constraints], ['B->D::f_+@FKKM:2008A'])

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

    def test_validate_semantics(self):

        item = eos.figure.ItemFactory.from_yaml("type: constraint2D\nconstraint: 'B->D::f_+@FKKM:2008A'\n"
                                                "x: { observable: 'test::unknown-observable' }")
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()
        diagnostics = list(item.validate_semantics(ValidationContext(description)))

        self.assertEqual(1, len(diagnostics))
        self.assertEqual(('x', 'observable'), diagnostics[0].path)

    def test_covariance_ellipse(self):

        # a MultivariateGaussian(Covariance) constraint yields the ellipse's widths and correlation
        item = eos.figure.ItemFactory.from_yaml("type: constraint2D\nconstraint: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'\n"
                                                "x: { observable: 'B->D::f_+(q2)' }\ny: { observable: 'B->D::f_0(q2)' }")
        item.prepare()

        self.assertEqual(item._shape, 'ellipse')
        self.assertGreater(item._xsigma, 0.0)
        self.assertGreater(item._ysigma, 0.0)
        self.assertGreater(abs(item._rho), 0.0)
        self.assertLessEqual(abs(item._rho), 1.0)

        _, ax = plt.subplots()
        item.draw(ax)

    def test_observable_not_contained(self):

        base = "type: constraint2D\nconstraint: 'B->D^(*)::FormFactors[f_+,f_0,A_0,A_1,A_2,V,T_1,T_2,T_23]@GKvD:2018A'\n"

        item = eos.figure.ItemFactory.from_yaml(base + "x: { observable: 'B->pi::f_+(q2)' }\n"
                                                       "y: { observable: 'B->D::f_0(q2)' }")
        with self.assertRaises(ValueError) as cm:
            item.prepare()
        self.assertIn('x-axis observable B->pi::f_+(q2) not contained', str(cm.exception))

        item = eos.figure.ItemFactory.from_yaml(base + "x: { observable: 'B->D::f_+(q2)' }\n"
                                                       "y: { observable: 'B->pi::f_0(q2)' }")
        with self.assertRaises(ValueError) as cm:
            item.prepare()
        self.assertIn('y-axis observable B->pi::f_0(q2) not contained', str(cm.exception))

        # the univariate (Gaussian) branch matches its single observable in the same way
        item = eos.figure.ItemFactory.from_yaml("type: constraint2D\nconstraint: 'B->D::f_+@FKKM:2008A'\n"
                                                "x: { observable: 'B->pi::f_+(q2)' }")
        with self.assertRaises(ValueError) as cm:
            item.prepare()
        self.assertIn('x-axis observable B->pi::f_+(q2) not contained', str(cm.exception))

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

    def test_legend(self):

        item = eos.figure.ItemFactory.from_yaml("type: signal-pdf\npdf: 'B->Dlnu::P(q2);l=mu'\n"
                                                "variable: q2\nrange: [0.02, 11.6]\nlabel: 'pdf'")
        entries = item.legend()

        self.assertEqual(len(entries), 1)
        self.assertIsInstance(entries[0][0], Line2D)
        self.assertEqual(entries[0][1], 'pdf')

    def test_unknown_pdf(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: signal-pdf
        pdf: 'B->Dlnu::NOSUCH(q2)'
        variable: q2
        range: [0.02, 11.6]
        resolution: 5
        """)
        with self.assertRaises(ValueError) as cm:
            item.prepare()

        self.assertIn('could not be created', str(cm.exception))

    def test_parameters_from_file(self):

        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'params.yaml'), 'w') as f:
                f.write("'mass::mu':\n  central: 0.5\n")

            item = eos.figure.ItemFactory.from_yaml("""
            type: signal-pdf
            pdf: 'B->Dlnu::P(q2);l=mu'
            variable: q2
            range: [0.02, 11.6]
            resolution: 5
            kinematics: { q2_min: 0.02, q2_max: 11.6 }
            parameters_from_file: 'params.yaml'
            """)
            with mock.patch.object(eos, 'warn') as warn:
                item.prepare(context=AnalysisFileContext(base_directory=tmp))

            self.assertEqual(item._parameters['mass::mu'].evaluate(), 0.5)
            self.assertEqual(len([call for call in warn.call_args_list if 'from file' in call[0][0]]), 1)

    def test_parameters_take_precedence_over_the_file(self):

        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'params.yaml'), 'w') as f:
                f.write("'mass::mu':\n  central: 0.5\n")

            item = eos.figure.ItemFactory.from_yaml("""
            type: signal-pdf
            pdf: 'B->Dlnu::P(q2);l=mu'
            variable: q2
            range: [0.02, 11.6]
            resolution: 5
            kinematics: { q2_min: 0.02, q2_max: 11.6 }
            parameters_from_file: 'params.yaml'
            parameters: { 'mass::mu': 0.2 }
            """)
            with mock.patch.object(eos, 'warn') as warn:
                item.prepare(context=AnalysisFileContext(base_directory=tmp))

            self.assertEqual(item._parameters['mass::mu'].evaluate(), 0.2)
            self.assertEqual(len([call for call in warn.call_args_list if 'with explicit values' in call[0][0]]), 1)

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

    def test_validate_semantics(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: complex-plane
        observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
        variables: ['Re{q2}', 'Im{q2}']
        ranges: [[0.1, 1.0], [0.1, 1.0]]
        fixed_parameters: { 'test::unknown-parameter': 1.0 }
        """)
        description = eos.analysis_file_description.AnalysisFileDescription.from_dict()
        diagnostics = list(item.validate_semantics(ValidationContext(description)))

        # the observable is known, the fixed parameter is not
        self.assertEqual(1, len(diagnostics))
        self.assertEqual(('fixed_parameters', 'test::unknown-parameter'), diagnostics[0].path)

    def test_fixed_kinematics(self):

        # kinematic variables other than the two swept ones are held at the given values
        item = eos.figure.ItemFactory.from_yaml("""
        type: complex-plane
        observable: 'B->D^*lnu::BRbar'
        variables: ['q2_e_min', 'q2_e_max']
        ranges: [[0.1, 0.2], [1.0, 1.1]]
        resolution: 2
        fixed_kinematics: { q2_mu_min: 0.02, q2_mu_max: 10.0 }
        """)

        self.assertIn('q2_mu_min', item._kinematics)
        self.assertIn('q2_mu_max', item._kinematics)
        self.assertEqual(float(item._kinematics['q2_mu_min']), 0.02)
        self.assertEqual(float(item._kinematics['q2_mu_max']), 10.0)

        item.prepare()

    def test_options(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: complex-plane
        observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
        variables: ['Re{q2}', 'Im{q2}']
        ranges: [[0.1, 1.0], [0.1, 1.0]]
        resolution: 3
        options: { contribution: 'Qc' }
        """)

        self.assertIn('contribution', item._options)
        self.assertEqual(str(item._options['contribution']), 'Qc')

        item.prepare()
        _, ax = plt.subplots()
        item.draw(ax)

    def test_fixed_parameters(self):

        item = eos.figure.ItemFactory.from_yaml("""
        type: complex-plane
        observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
        variables: ['Re{q2}', 'Im{q2}']
        ranges: [[0.1, 1.0], [0.1, 1.0]]
        resolution: 3
        fixed_parameters: { 'mass::mu': 0.2 }
        """)
        item.prepare()

        self.assertEqual(item._parameters['mass::mu'].evaluate(), 0.2)

    def test_fixed_parameters_take_precedence_over_the_file(self):

        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'params.yaml'), 'w') as f:
                f.write("'mass::mu':\n  central: 0.10566\n")

            item = eos.figure.ItemFactory.from_yaml("""
            type: complex-plane
            observable: 'b->s::Re{F17}(Re{q2},Im{q2})'
            variables: ['Re{q2}', 'Im{q2}']
            ranges: [[0.1, 1.0], [0.1, 1.0]]
            resolution: 3
            fixed_parameters_from_file: 'params.yaml'
            fixed_parameters: { 'mass::mu': 0.2 }
            """)
            with mock.patch.object(eos, 'warn') as warn:
                item.prepare(context=AnalysisFileContext(base_directory=tmp))

            self.assertEqual(item._parameters['mass::mu'].evaluate(), 0.2)
            self.assertEqual(len([call for call in warn.call_args_list if 'with explicit values' in call[0][0]]), 1)

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

    def test_asymmetric_x_errors(self):

        # an x error given as a pair is stored as (minus, plus), a scalar one on both sides
        item = eos.figure.ItemFactory.from_yaml("type: 'errorbars'\npositions: [[1, 2], [2, 3]]\n"
                                                "xerrors: [[0.1, 0.2], 0.5]\nyerrors: [0.2, 0.3]")
        item.prepare()

        self.assertEqual(list(item._xerr[0]), [0.1, 0.5])
        self.assertEqual(list(item._xerr[1]), [0.2, 0.5])

        _, ax = plt.subplots()
        item.draw(ax)

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
