import unittest

import eos
from eos import figure
import os
import pypmc
import numpy as np

class FigureTests(unittest.TestCase):

    def test_single_figure(self):
        "Test creating a single figure."

        DATA = os.path.join(os.environ['SOURCE_DIR'], "eos/data/native_TEST.d/samples")
        data = eos.data.ImportanceSamples(DATA)
        kwargs = {
            'title': 'My first figure',
            'plot': {
                'description': {
                    'x_label': data.varied_parameters[0]['name'],
                },
                'items': [
                    {
                        'type': 'histogram',
                        'description': {
                            'data_file': DATA,
                            'variable': 'CKM::abs(V_ub)'
                        }
                    }
                ]
            }
        }
        f = eos.figure.SingleFigure(**kwargs)
        f.draw()


    def test_grid_figure(self):
        "Test creating a grid figure."

        DATA = os.path.join(os.environ['SOURCE_DIR'], "eos/data/native_TEST.d/samples")
        data = eos.data.ImportanceSamples(DATA)
        kwargs = {
            'title': 'My second figure',
            'ncols': 2,
            'nrows': 2,
            'plots': [
                {
                    'description': {
                        'x_label': data.varied_parameters[0]['name'],
                    },
                    'items': [
                        {
                            'type': 'histogram',
                            'description': {
                                'data_file': DATA,
                                'variable': 'CKM::abs(V_ub)'
                            }
                        }
                    ]
                },
                {
                },
                {

                },
                {
                    'description': {
                        'x_label': data.varied_parameters[1]['name'],
                    },
                    'items': [
                        {
                            'type': 'histogram',
                            'description': {
                                'data_file': DATA,
                                'variable': 'B->pi::f_+(0)@BCL2008'
                            }
                        }
                    ]
                }
            ]
        }
        f = GridFigure(**kwargs)
        f.draw()


if __name__ == '__main__':
    unittest.main(verbosity=5)
