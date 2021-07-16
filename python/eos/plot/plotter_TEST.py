import unittest
import eos

class PlotterObservableVariableTests(unittest.TestCase):
    "Test handling of argument 'variable' in Plotter.Observable"

    def setUp(self):
        self.plot_args = {
                'plot': {
                    'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [0.0, 11.63] },
                    'y': { 'label': r'$d\mathcal{B}/dq^2$',                    'range': [0.0,  5e-3] },
                    'legend': { 'location': 'lower left' }
                },
                'contents': [
                    {
                        'label': r'$\ell=e$',
                        'type': 'observable',
                        'observable': 'B->Dlnu::dBR/dq2;l=e,q=d',
                        'variable': 'q2',
                        'parameters': {
                            'mass::mu': 1.0,
                            'mass::tau': 1.0,
                        },
                        'color': 'black',
                        'range': [0.02, 11.63],
                    },
                ]
            }


    def test_variabe_key_0(self):
        "Check valid case: variable is kinematic"

        eos.plot.Plotter(self.plot_args).plot()


    def test_variabe_key_1(self):
        "Check valid case: variable is parameter"

        self.plot_args['contents'][0]['variable'] = 'mass::tau'
        self.plot_args['contents'][0]['kinematics'] = {'q2': 0.1}
        eos.plot.Plotter(self.plot_args).plot()


    def test_variable_key_2(self):
        "Handle 'variable' key not given"

        self.plot_args['contents'][0].pop('variable')
        with self.assertRaisesRegex(ValueError,
            "Missing key for plot of observable 'B->Dlnu::dBR/dq2;l=e,q=d': 'variable'"):
            eos.plot.Plotter(self.plot_args).plot()


    def test_variable_key_3(self):
        "Handle value of 'variable' not a kinematic var. nor a parameter"

        self.plot_args['contents'][0]['variable'] = 'neitherKinNorPar'

        with self.assertRaisesRegex(ValueError,
            "Value of 'variable' for observable 'B->Dlnu::dBR/dq2;l=e,q=d' is neither " +
            "a valid kinematic variable nor parameter: 'neitherKinNorPar'"):
            eos.plot.Plotter(self.plot_args).plot()


if __name__ == '__main__':
    unittest.main(verbosity=5)

