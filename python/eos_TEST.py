# vim: set sts=4 tw=120 :

import inspect
import os

class TestFailedError(BaseException):
    def __init__(self, msg):
        self.msg = msg

class PythonTests:
    """
    Check that the first term in the PYTHONPATH environment variable ends in '.libs'.
    """
    def check_000_PYTHONPATH(self):
        try:
            pythonpath = os.environ['PYTHONPATH'].split(os.pathsep)[0]
            if not pythonpath.endswith('python/.libs/'):
                raise TestFailedError('PYTHONPATH not correctly set for running test cases')

        except KeyError:
            raise TestFailedError('PYTHONPATH not set')

    """
    Check if the module Python/C++ interface can be imported.
    """
    def check_001_cpp_module(self):
        try:
            import _eos
        except:
            raise TestFailedError('importing \'_eos\' failed.')

        try:
            import eos
        except:
            raise TestFailedError('importing \'eos\' failed.')

    """
    Check if an instance of QualifiedName can be created.
    """
    def check_002_QualifiedName(self):
        from eos import QualifiedName

        try:
            qn = QualifiedName('B->K^*ll::A_FB(s)@LargeRecoil;foo=bar,baz=harr')
        except:
            raise TestFailedError('cannot initialize QualifiedName')

    """
    Check if an instance of Kinematics can be created.
    """
    def check_003_Kinematics(self):
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

    """
    Check if an instance of Options can be created.
    """
    def check_004_Options(self):
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

    """
    Check if an instance of Parameters can be created.
    """
    def check_005_Parameters(self):
        from eos import Parameters

        try:
            par = Parameters.Defaults()
        except:
            raise TestFailedError('cannot initialize default Parameters')

        try:
            electron_mass = par['mass::e']
        except:
            raise TestFailedError('cannot lookup existing Parameter \'mass::e\'')

    """
    Check if an instance of Observable can be created.
    """
    def check_006_Observable(self):
        from eos import Observable, Parameters, Kinematics, Options, QualifiedName

        obs = None
        try:
            obs = Observable.make('B->K^*ll::A_FB(s)@LargeRecoil', Parameters.Defaults(), Kinematics(s=1.0),
                    Options(model='SM'))
        except:
            raise TestFailedError('cannot create Observable')

        if not obs:
            raise TestFailedError('lookup of Observable failed')

        if not obs.name() == 'B->K^*ll::A_FB(s)@LargeRecoil':
            raise TestFailedError('cannot obtain Observable name')

        value = None
        try:
            value = obs.evaluate()
        except:
            raise TestFailedError('cannote evaluate Observable')

    """
    Check if an instance of Constraint can be created, and if we can iterate over
    its LogLikelihoodBlock range.
    """
    def check_007_Constraint(self):
        from eos import Constraint, Options

        con = None
        try:
            con = Constraint.make('B->D::f_++f_0@HPQCD2015A', Options())
        except:
            raise TestFailedError('cannot create Constraint')

        if not con.name() == 'B->D::f_++f_0@HPQCD2015A':
            raise TestFailedError('cannot obtain Constraint name')

        for llhb in con.blocks():
            print(llhb.as_string())

        for obs in con.observables():
            print(obs.name())



"""
Run all test cases.
"""
tests = PythonTests()
for (name, testcase) in inspect.getmembers(tests, predicate=inspect.ismethod):
    print("%s: " % name)
    try:
        testcase()
    except TestFailedError as e:
        print("    failed with error: %s" % e.msg)
        exit(1)
    print("    passed")
