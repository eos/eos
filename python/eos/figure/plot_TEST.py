# Copyright (c) 2024 Danny van Dyk
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

if __name__ == '__main__':
    unittest.main(verbosity=5)
