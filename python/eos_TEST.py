# vim: set sts=4 tw=120 :

import unittest
import inspect
import os

class TestFailedError(BaseException):
    def __init__(self, msg):
        self.msg = msg

class PythonTests:
    def check_000_PYTHONPATH(self):
        """Check that the first term in the PYTHONPATH environment variable ends in '.libs'."""
        try:
            pythonpath = os.environ['PYTHONPATH'].split(os.pathsep)[0]
            print('PYTHONPATH={}'.format(os.environ['PYTHONPATH']))
            if not pythonpath.endswith('python/.libs/'):
                raise TestFailedError('PYTHONPATH not correctly set for running test cases')

        except KeyError:
            raise TestFailedError('PYTHONPATH not set')

    def check_001_cpp_module(self):
        """Check if the module Python/C++ interface can be imported."""
        try:
            import _eos
        except ImportError as e:
            raise TestFailedError(f'importing \'_eos\' failed: \'{str(e)}\'')

        try:
            import eos
        except ImportError as e:
            raise TestFailedError(f'importing \'eos\' failed: \'{str(e)}\'')

    def check_002_QualifiedName(self):
        """Check if an instance of QualifiedName can be created."""
        from eos import QualifiedName

        try:
            qn = QualifiedName('B->K^*ll::A_FB(s)@LargeRecoil;foo=bar,baz=harr')
        except:
            raise TestFailedError('cannot initialize QualifiedName')

    def check_003_Kinematics(self):
        """Check if an instance of Kinematics can be created."""
        from eos import Kinematics

        try:
            kin = Kinematics()
        except:
            raise TestFailedError('cannot initialize empty Kinematics')

        try:
            kin = Kinematics(s_min=1.0, s_max=6.0)
        except:
            raise TestFailedError('cannot initialize Kinematics with kwargs')

        try:
            kin = Kinematics(**{'s_min': 1.0, 's_max': 6.0})
        except:
            raise TestFailedError('cannot initialize Kinematics with explicit kwargs')

    def check_004_Options(self):
        """Check if an instance of Options can be created."""
        from eos import Options

        try:
            opt = Options()
        except:
            raise TestFailedError('cannot initialize empty Options')

        try:
            kin = Options(foo='bar', baz='harr')
        except:
            raise TestFailedError('cannot initialize Options with kwargs')

        try:
            kin = Options(**{'foo': 'bar', 'baz': 'harr'})
        except:
            raise TestFailedError('cannot initialize Options with explicit kwargs')

    def check_005_Parameters(self):
        """Check if an instance of Parameters can be created."""
        from eos import Parameters

        try:
            par = Parameters.Defaults()
        except:
            raise TestFailedError('cannot initialize default Parameters')

        try:
            electron_mass = par['mass::e']
        except:
            raise TestFailedError('cannot lookup existing Parameter \'mass::e\'')

    def check_006_Observable(self):
        """Check if an instance of Observable can be created."""
        from eos import Observable, Parameters, Kinematics, Options, QualifiedName

        obs = None
        try:
            obs = Observable.make('B->Dlnu::BR', Parameters.Defaults(), Kinematics(q2_min=0.02, q2_max=10),
                    Options(model='SM'))
        except:
            raise TestFailedError('cannot create Observable')

        if not obs:
            raise TestFailedError('lookup of Observable failed')

        if not obs.name() == 'B->Dlnu::BR':
            raise TestFailedError('cannot obtain Observable name')

        value = None
        try:
            value = obs.evaluate()
        except:
            raise TestFailedError('cannote evaluate Observable')

    def check_007_Constraint(self):
        """
        Check if an instance of Constraint can be created, and if we can iterate over
        its LogLikelihoodBlock range.
        """
        from eos import Constraint, Options

        con = None
        try:
            con = Constraint.make('B->D::f_++f_0@HPQCD:2015A', Options())
        except:
            raise TestFailedError('cannot create Constraint')

        if not con.name() == 'B->D::f_++f_0@HPQCD:2015A':
            raise TestFailedError('cannot obtain Constraint name')

        for llhb in con.blocks():
            print(llhb)

        for obs in con.observables():
            print(obs.name())


    """
    Check if an instance of Model can be created, and if we can compute the
    running MSbar mass of the b quark.
    """
    def check_008_Model(self):
        from eos import Model, Parameters, Options

        m = None
        p = None
        try:
            p = Parameters.Defaults()
            o = Options()
            m = Model.make('SM', p, o)
        except:
            raise TestFailedError('cannot create Model')

        try:
            pvalue = p['mass::b(MSbar)'].evaluate()
            mvalue = m.m_b_msbar(pvalue)
            if not pvalue == mvalue:
                raise TestFailedError('internal error')
        except:
            raise TestFailedError('cannot determine running b quark mass')


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
                    eos.Options(**{'U': 'c'})) # will  be overwritten by the
                                               # implementation of the observable
            self.assertEqual(cm.output,
                    [r"""ERROR:EOS:[ConcreteObservableEntry.make] Observable 'B->pilnu::BR' forces option key 'U' to value 'u', overriding user-provided value 'c'"""])


# Run legacy test cases
tests = PythonTests()
for (name, testcase) in inspect.getmembers(tests, predicate=inspect.ismethod):
    print("%s: " % name)
    try:
        testcase()
    except TestFailedError as e:
        print("    failed with error: %s" % e.msg)
        exit(1)
    print("    passed")

# Run new tests
unittest.main()
