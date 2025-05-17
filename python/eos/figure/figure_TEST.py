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
                xvalues: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
              - type: 'observable'
                observable: 'B->Dlnu::dBR/dq2'
                options: { 'l': 'mu' }
                variable: 'q2'
                xvalues: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
            """
            figure = eos.figure.FigureFactory.from_yaml(input)
            figure.draw()
        except Exception as e:
            self.fail(f"Error when testing figure of type 'single': {e}")


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
