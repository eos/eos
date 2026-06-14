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

if __name__ == '__main__':
    unittest.main(verbosity=5)
