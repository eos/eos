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
