# vim: set sts=4 tw=120 :

import unittest
import os
import numpy as _np


class EnvironmentTests(unittest.TestCase):

    def test_PYTHONPATH(self):
        "Check that the first term in the PYTHONPATH environment variable ends in '.libs'."
        self.assertIn('PYTHONPATH', os.environ, 'PYTHONPATH not set')
        pythonpath = os.environ['PYTHONPATH'].split(os.pathsep)[0]
        self.assertTrue(pythonpath.endswith('python/.libs/'),
                'PYTHONPATH not correctly set for running test cases')

    def test_cpp_module(self):
        "Check if the module Python/C++ interface can be imported."
        try:
            import _eos
        except ImportError as e:
            self.fail(f'importing \'_eos\' failed: \'{str(e)}\'')

        try:
            import eos
        except ImportError as e:
            self.fail(f'importing \'eos\' failed: \'{str(e)}\'')


class QualifiedNameTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of QualifiedName can be created."
        from eos import QualifiedName

        try:
            QualifiedName('B->K^*ll::A_FB(s)@LargeRecoil;foo=bar,baz=harr')
        except Exception as e:
            self.fail(f'cannot initialize QualifiedName: {e}')


class KinematicsTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of Kinematics can be created."
        from eos import Kinematics

        try:
            Kinematics()
        except Exception as e:
            self.fail(f'cannot initialize empty Kinematics: {e}')

        try:
            Kinematics(s_min=1.0, s_max=6.0)
        except Exception as e:
            self.fail(f'cannot initialize Kinematics with kwargs: {e}')

        try:
            Kinematics(**{'s_min': 1.0, 's_max': 6.0})
        except Exception as e:
            self.fail(f'cannot initialize Kinematics with explicit kwargs: {e}')


class OptionsTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of Options can be created."
        from eos import Options

        try:
            Options()
        except Exception as e:
            self.fail(f'cannot initialize empty Options: {e}')

        try:
            Options(foo='bar', baz='harr')
        except Exception as e:
            self.fail(f'cannot initialize Options with kwargs: {e}')

        try:
            Options(**{'foo': 'bar', 'baz': 'harr'})
        except Exception as e:
            self.fail(f'cannot initialize Options with explicit kwargs: {e}')


class ParametersTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of Parameters can be created."
        from eos import Parameters

        try:
            par = Parameters.Defaults()
        except Exception as e:
            self.fail(f'cannot initialize default Parameters: {e}')

        self.assertFalse(_np.any([p.unit().__str__() == 'undefined' for p in par]),
                'some parameters have undefined units')

        try:
            par['mass::e']
        except Exception as e:
            self.fail(f'cannot lookup existing Parameter \'mass::e\': {e}')


class ObservableTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of Observable can be created."
        from eos import Observable, Parameters, Kinematics, Options

        try:
            obs = Observable.make('B->Dlnu::BR', Parameters.Defaults(), Kinematics(q2_min=0.02, q2_max=10),
                    Options(model='SM'))
        except Exception as e:
            self.fail(f'cannot create Observable: {e}')

        self.assertIsNotNone(obs, 'lookup of Observable failed')
        self.assertEqual(obs.name(), 'B->Dlnu::BR', 'cannot obtain Observable name')

        try:
            obs.evaluate()
        except Exception as e:
            self.fail(f'cannot evaluate Observable: {e}')


class ConstraintTests(unittest.TestCase):

    def test_creation(self):
        """
        Check if an instance of Constraint can be created, and if we can iterate over
        its LogLikelihoodBlock range.
        """
        from eos import Constraint, Options

        try:
            con = Constraint.make('B->D::f_++f_0@HPQCD:2015A', Options())
        except Exception as e:
            self.fail(f'cannot create Constraint: {e}')

        self.assertEqual(con.name(), 'B->D::f_++f_0@HPQCD:2015A', 'cannot obtain Constraint name')

        for llhb in con.blocks():
            print(llhb)

        for obs in con.observables():
            print(obs.name())


class ModelTests(unittest.TestCase):

    def test_creation(self):
        """
        Check if an instance of Model can be created, and if we can compute the
        running MSbar mass of the b quark.
        """
        from eos import Model, Parameters, Options

        try:
            p = Parameters.Defaults()
            o = Options()
            m = Model.make('SM', p, o)
        except Exception as e:
            self.fail(f'cannot create Model: {e}')

        try:
            pvalue = p['mass::b(MSbar)'].evaluate()
            mvalue = m.m_b_msbar(pvalue)
        except Exception as e:
            self.fail(f'cannot determine running b quark mass: {e}')

        self.assertEqual(pvalue, mvalue, 'internal error')


class UnitTests(unittest.TestCase):

    def test_latex(self):
        "Check if the latex strings of the Units are valid."
        from eos import Unit
        import matplotlib

        for attr, value in Unit.__dict__.items():
            if attr[0] != "_" and attr != "latex":
                s = '<problem in test logic>'
                try:
                    s = value.__func__().latex()
                    matplotlib.texmanager.TexManager.get_text_width_height_descent(f"${s}$", 1)
                except Exception as e:
                    self.fail(f'invalid latex string "${s}$" for unit \'{attr}\': {e}')


class LoggingTests(unittest.TestCase):

    def test_log_levels(self):
        "Check if the set log level is respected"

        import eos, logging, _eos
        from eos import _NativeLogLevel as ll
        levels = [ll.DEBUG, ll.INFO, ll.WARNING, ll.ERROR]

        for set_level in levels:
            _eos._set_native_log_level(set_level)

            for check_level in levels:

                if check_level <= set_level:
                    with self.assertLogs('EOS', level='DEBUG') as cm:
                        _eos._emit_native_log("myId", check_level, "msg")
                else:
                    # workaround for old unittest version, w/o assertNoLogs
                    with self.assertRaisesRegex(AssertionError,
                        "no logs of level DEBUG or higher triggered on EOS"):
                        with self.assertLogs('EOS', level='DEBUG') as cm:
                            _eos._emit_native_log("id", check_level, "msg")

    def test_log_000(self):
        "Computation of specific observable should log error"
        import eos

        with self.assertLogs('EOS', level='WARNING') as cm:
            eos.Observable.make(
                    "B->pilnu::BR",
                    eos.Parameters.Defaults(),
                    eos.Kinematics(q2_min=1, q2_max=6),
                    eos.Options(**{'P': 'D'})) # will  be overwritten by the
                                               # implementation of the observable
            self.assertEqual(cm.output,
                    [r"""ERROR:EOS:[ConcreteObservableEntry.make] Observable 'B->pilnu::BR' forces option key 'P' to value 'pi', overriding user-provided value 'D'"""])


if __name__ == '__main__':
    unittest.main(verbosity=5)
