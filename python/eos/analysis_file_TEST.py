import unittest

import contextlib
import io
import eos
import os
from pathlib import Path

_TESTD = Path(__file__).parent / 'analysis_file_TEST.d'

class TestAnalysisFile(unittest.TestCase):

    def test_analysis_file(self):

        af = eos.AnalysisFile(_TESTD / 'valid-analysis-file.yaml')
        af.validate()

    def test_minimal_analysis_file(self):

        af = eos.AnalysisFile(_TESTD / 'minimal-analysis-file.yaml')
        af.validate()


class TestAnalysisFileConstructionErrors(unittest.TestCase):
    """The mandatory-structure and uniqueness checks in ``__init__`` must raise."""

    def test_nonexistent_file(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'this-file-does-not-exist.yaml')

    def test_not_a_file(self):
        # a directory is not a valid analysis file
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD)

    def test_missing_priors(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'no-priors.yaml')

    def test_missing_posteriors(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'no-posteriors.yaml')

    def test_unknown_prior_reference(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'unknown-prior-ref.yaml')

    def test_unknown_likelihood_reference(self):
        with self.assertRaises(RuntimeError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'unknown-likelihood-ref.yaml')

    def test_invalid_parameter_name(self):
        with self.assertRaises(ValueError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'bad-parameter.yaml')

    def test_duplicate_step_id(self):
        with self.assertRaises(ValueError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'duplicate-step-id.yaml')

    def test_duplicate_mask_name(self):
        with self.assertRaises(ValueError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'duplicate-mask-name.yaml')

    def test_unknown_mask_reference(self):
        with self.assertRaises(ValueError):
            eos.AnalysisFile(_TESTD / 'invalid' / 'unknown-mask-ref.yaml')


class TestAnalysisFileMethods(unittest.TestCase):
    """Exercises the accessors and factory methods on a fully-populated analysis file."""

    def setUp(self):
        # a fresh instance per test: analysis()/observables() read (and, for the latter, mutate in
        # place) the cached descriptions, so tests must not share a single instance
        self.af = eos.AnalysisFile(_TESTD / 'valid-analysis-file.yaml')

    def test_properties(self):
        self.assertIn('CKM', self.af.priors)
        self.assertIn('TH-pi', self.af.likelihoods)
        self.assertIn('CKM-all', self.af.posteriors)
        self.assertIn('leptonic-BR-CKM', self.af.predictions)
        self.assertIn('mask-A', self.af.masks)

    def test_repr_html(self):
        # the HTML representation renders every section present in the fixture
        html = self.af._repr_html_()
        for section in ('PRIORS', 'LIKELIHOODS', 'POSTERIORS', 'PREDICTIONS',
                        'OBSERVABLES', 'PARAMETERS', 'STEPS', 'MASKS'):
            self.assertIn(section, html)

    def test_dump(self):
        # dump() writes the input data back out as YAML
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            self.af.dump()
        self.assertIn('priors', buffer.getvalue())

    def test_analysis(self):
        analysis = self.af.analysis('CKM-all')
        self.assertIsNotNone(analysis)

    def test_analysis_unknown_posterior(self):
        with self.assertRaises(RuntimeError):
            self.af.analysis('does-not-exist')

    def test_observables(self):
        observables = self.af.observables('CKM-all', 'leptonic-BR-CKM', eos.Parameters())
        self.assertGreater(len(observables), 0)

    def test_observables_unknown_posterior(self):
        with self.assertRaises(RuntimeError):
            self.af.observables('does-not-exist', 'leptonic-BR-CKM', eos.Parameters())

    def test_observables_unknown_prediction(self):
        with self.assertRaises(RuntimeError):
            self.af.observables('CKM-all', 'does-not-exist', eos.Parameters())

    def test_observable(self):
        observable = self.af.observable('CKM-all', 'B_u->lnu::BR;l=e', eos.Parameters())
        self.assertIsNotNone(observable)

    def test_observable_applies_fixed_parameters(self):
        # 'WET-all' fixes 'CKM::abs(V_ub)', so observable() writes that value into the parameter set
        parameters = eos.Parameters()
        self.af.observable('WET-all', 'B_u->lnu::BR;l=e', parameters)
        self.assertAlmostEqual(float(parameters['CKM::abs(V_ub)']), 3.67e-3)

    def test_observable_with_option_part_in_name(self):
        # an observable name may carry an option-part; it is created under the posterior's global
        # options merged with that name
        observable = self.af.observable('WET-all', 'B_u->lnu::BR;l=e,model=CKM', eos.Parameters())
        self.assertIsNotNone(observable)

    def test_observable_unknown_posterior(self):
        with self.assertRaises(RuntimeError):
            self.af.observable('does-not-exist', 'B_u->lnu::BR;l=e', eos.Parameters())


class TestAnalysisFileValidation(unittest.TestCase):
    """The reporting (non-raising) branches of ``validate``."""

    def test_semantic_errors_are_reported(self):
        af = eos.AnalysisFile(_TESTD / 'invalid' / 'semantic-errors.yaml')
        messages = af.validate()
        # validate() collects one message per problem rather than raising; the fixture contains
        # many independent problems, so it must return a non-empty list
        self.assertGreater(len(messages), 0)
        blob = '\n'.join(messages)
        # a representative selection of the distinct problem classes in the fixture
        self.assertIn('Not::a-parameter', blob)                 # unknown prior parameter
        self.assertIn('Not::a-constraint@Nowhere:2000A', blob)  # unknown constraint
        self.assertIn('matches an already defined constraint', blob)  # colliding manual constraint
        self.assertIn('B->pilnu::NOTREAL', blob)                # unknown prediction observable
        self.assertIn('Not::a-fixed-parameter', blob)           # unknown posterior fixed parameter
        self.assertIn('Not::a-prediction-parameter', blob)      # unknown prediction fixed parameter
        self.assertIn('unknown steps', blob)                    # step dependency
        self.assertIn('unknown tasks', blob)                    # step default arguments
        self.assertIn("Posterior 'nonexistent-posterior'", blob)  # step task posterior
        self.assertIn('repeatedly', blob)                       # repeated mask expression name


if __name__ == '__main__':
    unittest.main(verbosity=5)
