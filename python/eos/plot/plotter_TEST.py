import unittest
import copy
import eos

class PlotterObservableVariableTests(unittest.TestCase):
    "Test handling of argument 'variable' in Plotter.Observable"

    @classmethod
    def setUpClass(cls):

        cls.plot_args_prototype = {
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

        # All tests call Plotter.plot(), which relies on the presence of a
        # LaTeX installation. The tests of this class should be skipped if
        # that is not given.
        try:
            eos.plot.Plotter(cls.plot_args_prototype).plot()
        except Exception as e:
            msg = "Plotting fails: " + repr(e)
            raise unittest.SkipTest(msg)


    def setUp(self):
        self.plot_args = copy.deepcopy(self.plot_args_prototype)


    def test_variable_key_0(self):
        "Check valid case: variable is kinematic"

        eos.plot.Plotter(self.plot_args).plot()


    def test_variable_key_1(self):
        "Check valid case: variable is parameter"

        self.plot_args['contents'][0]['variable'] = 'mass::tau'
        self.plot_args['contents'][0]['parameters'].pop('mass::tau')
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


    def test_variable_key_4(self):
        "Pass the variable key again as a fix kinematic"

        self.plot_args['contents'][0]['kinematics'] = {'q2': 0.1}

        with self.assertRaisesRegex(ValueError,
            "Variable 'q2' of observable 'B->Dlnu::dBR/dq2;l=e,q=d' is also " +
            "specified as a fix kinematic with value 0.1"):
            eos.plot.Plotter(self.plot_args).plot()


    def test_variable_key_5(self):
        "Pass the variable key again as a fix parameter"

        self.plot_args['contents'][0]['variable'] = 'mass::tau'

        with self.assertRaisesRegex(ValueError,
            "Variable 'mass::tau' of observable 'B->Dlnu::dBR/dq2;l=e,q=d' is also " +
            "specified as a fix parameter with value 1.0"):
            eos.plot.Plotter(self.plot_args).plot()


    def test_kinematics_key_0(self):
        "Handle invalid kinematics key"

        self.plot_args['contents'][0]['kinematics'] = {'q3': 0.1}

        with self.assertRaisesRegex(ValueError,
            "Kinematic quantity 'q3' does not match known ones for observable " +
            "'B->Dlnu::dBR/dq2;l=e,q=d': \['q2'\]"):
            eos.plot.Plotter(self.plot_args).plot()


if __name__ == '__main__':
    unittest.main(verbosity=5)

