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

        # syntactically invalid names must raise
        with self.assertRaises(Exception):
            QualifiedName('')

    def test_parts(self):
        "Check if the individual parts of a QualifiedName are parsed correctly."
        from eos import QualifiedName

        qn = QualifiedName('B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010')

        self.assertEqual(str(qn.prefix_part()), 'B->K^*ll')
        self.assertEqual(str(qn.name_part()),   'A_FB(s)')
        self.assertEqual(str(qn.suffix_part()), 'LargeRecoil')

        # __str__ yields the short name (suffix included, options excluded)
        self.assertEqual(str(qn), 'B->K^*ll::A_FB(s)@LargeRecoil')
        # full() yields the complete name, including options
        self.assertEqual(qn.full(), 'B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010')

    def test_equality(self):
        "Check that QualifiedNames compare equal based on their short name only."
        from eos import QualifiedName

        # the two names share a short name and differ only in their options
        qn1 = QualifiedName('B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010')
        qn2 = QualifiedName('B->K^*ll::A_FB(s)@LargeRecoil;form-factors=BSZ2015')
        self.assertEqual(qn1, qn2)

        # a differing suffix changes the short name
        qn3 = QualifiedName('B->K^*ll::A_FB(s)@LowRecoil')
        self.assertNotEqual(qn1, qn3)


class ReferenceNameTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of ReferenceName can be created."
        from eos import ReferenceName

        try:
            ReferenceName('IKMvD:2014A')
        except Exception as e:
            self.fail(f'cannot initialize ReferenceName: {e}')

        # syntactically invalid names must raise
        for invalid in ['', 'A', 'A:199', 'A:1999', 'A:1999-', '[A:1999-B]']:
            with self.assertRaises(Exception):
                ReferenceName(invalid)

    def test_str_and_equality(self):
        "Check the string representation and equality of ReferenceName."
        from eos import ReferenceName

        rn = ReferenceName('LHCb:2010A')
        self.assertEqual(str(rn), 'LHCb:2010A')

        self.assertEqual(ReferenceName('BES2:2006A'), ReferenceName('BES2:2006A'))
        self.assertNotEqual(ReferenceName('BES2:2006A'), ReferenceName('LHCb:2010A'))


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


class SignalPDFTests(unittest.TestCase):

    def test_creation(self):
        "Check if an instance of SignalPDF can be created and evaluated."
        from eos import SignalPDF, Parameters, Kinematics, Options
        from math import log

        # 'TestLegendre1D::P(z)' is a PDF provided for testing purposes:
        # its unnormalized density is pdf(z) = z (4 - z) = 4 z - z^2, with zeros at z = 0
        # and z = 4 and a maximum at z = 2.
        try:
            pdf = SignalPDF.make('TestLegendre1D::P(z)', Parameters.Defaults(),
                    Kinematics(z=2.0, z_min=0.0, z_max=4.0), Options())
        except Exception as e:
            self.fail(f'cannot create SignalPDF: {e}')

        self.assertIsNotNone(pdf, 'lookup of SignalPDF failed')
        self.assertEqual(pdf.name(), 'TestLegendre1D::P(z)', 'cannot obtain SignalPDF name')

        try:
            value = pdf.evaluate()
        except Exception as e:
            self.fail(f'cannot evaluate SignalPDF: {e}')

        # evaluate() returns the logarithm of the unnormalized PDF;
        # pdf(z=2) = 4 * 2 - 2^2 = 4
        self.assertAlmostEqual(value, log(4.0), delta=1.0e-10)

        try:
            norm = pdf.normalization()
        except Exception as e:
            self.fail(f'cannot evaluate SignalPDF normalization: {e}')

        # normalization() returns the logarithm of the normalization integral;
        # int_{0}^{4} dz (4 z - z^2) = 2 * 16 - 64 / 3 = 32 / 3
        self.assertAlmostEqual(norm, log(32.0 / 3.0), delta=1.0e-10)

    def test_runtime_insertion(self):
        "Check that a new SignalPDF can be inserted at run time and evaluated."
        from eos import SignalPDF, SignalPDFs, Parameters, Kinematics, Options
        from math import log

        signal_pdfs = SignalPDFs()

        # reuse the observables backing the built-in 'TestLegendre1D::P(z)' PDF:
        # pdf(z) = z (4 - z), int_{0}^{4} dz pdf(z) = 32 / 3
        try:
            signal_pdfs.insert('TestLegendre1D::RuntimePyP(z)', 'runtime-inserted test PDF', Options(),
                    'TestLegendre1D::UnnormalizedPDF(z)', ['z'],
                    'TestLegendre1D::NormalizationPDF(z)', ['z_min', 'z_max'])
        except Exception as e:
            self.fail(f'cannot insert SignalPDF: {e}')

        # referencing an unknown observable must raise
        with self.assertRaisesRegex(RuntimeError, "is not a known observable"):
            signal_pdfs.insert('TestLegendre1D::RuntimePyBad(z)', 'bad', Options(),
                    'TestLegendre1D::Unknown(z)', ['z'],
                    'TestLegendre1D::NormalizationPDF(z)', ['z_min', 'z_max'])

        try:
            pdf = SignalPDF.make('TestLegendre1D::RuntimePyP(z)', Parameters.Defaults(),
                    Kinematics(z=2.0, z_min=0.0, z_max=4.0), Options())
        except Exception as e:
            self.fail(f'cannot create runtime-inserted SignalPDF: {e}')

        self.assertIsNotNone(pdf, 'lookup of runtime-inserted SignalPDF failed')
        self.assertEqual(pdf.name(), 'TestLegendre1D::RuntimePyP(z)', 'cannot obtain SignalPDF name')

        # evaluate() returns log(pdf(z = 2)) = log(4)
        self.assertAlmostEqual(pdf.evaluate(), log(4.0), delta=1.0e-10)

        # normalization() returns log(32 / 3)
        self.assertAlmostEqual(pdf.normalization(), log(32.0 / 3.0), delta=1.0e-10)


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


class ExternalLogPriorTests(unittest.TestCase):

    def test_creation(self):
        "Check if an external LogPrior can be created and evaluated."
        import eos
        from math import log, pi, sqrt

        parameters = eos.Parameters.Defaults()

        class _GaussianProvider:
            def __init__(self, parameters):
                self._p = parameters['mass::b(MSbar)']
                self._mu, self._sigma = 4.18, 0.02

            varied_parameters = ['mass::b(MSbar)']
            informative = True

            def evaluate(self):
                x = self._p.evaluate()
                return -0.5 * ((x - self._mu) / self._sigma) ** 2 - log(self._sigma * sqrt(2.0 * pi))

            def sample(self):
                self._p.set(self._mu)

            def compute_cdf(self):
                self._p.set_generator(0.5)

        try:
            prior = eos.LogPrior.External(parameters, _GaussianProvider)
        except Exception as e:
            self.fail(f'cannot create external LogPrior: {e}')

        # the prior declares the parameters it varies
        self.assertEqual([p.name() for p in prior.varied_parameters()], ['mass::b(MSbar)'])

        # evaluate() delegates to the provider
        parameters['mass::b(MSbar)'].set(4.18)
        self.assertAlmostEqual(prior.evaluate(), -log(0.02 * sqrt(2.0 * pi)), places=13)

        # sample() and compute_cdf() delegate to the provider
        parameters['mass::b(MSbar)'].set(0.0)
        prior.sample()
        self.assertAlmostEqual(parameters['mass::b(MSbar)'].evaluate(), 4.18, places=13)

        # the prior can be used as part of a posterior
        posterior = eos.LogPosterior(eos.LogLikelihood(parameters))
        self.assertTrue(posterior.add(prior, False))
        self.assertAlmostEqual(posterior.evaluate(), prior.evaluate(), places=13)

    def test_incomplete_provider(self):
        "Check that a provider lacking a required method is rejected at creation."
        import eos

        parameters = eos.Parameters.Defaults()

        class _IncompleteProvider:
            def __init__(self, parameters):
                pass

            varied_parameters = []
            informative = False

            def evaluate(self):
                return 0.0

            def compute_cdf(self):
                pass

        # the missing 'sample()' method must be caught when the prior is created
        with self.assertRaises(Exception):
            eos.LogPrior.External(parameters, _IncompleteProvider)


if __name__ == '__main__':
    unittest.main(verbosity=5)
