# Copyright (c) 2026 Danny van Dyk
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

import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import yaml

from eos.cli.data import create_check_factory
from eos.cli.data_checks import (
    CheckContext,
    CheckFactory,
    CheckResult,
    CheckRunner,
    CheckScope,
    ExitStatus,
    Finding,
    PlainTextRenderer,
    Severity,
    _normalize_name,
    check_analysis_metadata,
    run_checks,
)


FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/data_checks_TEST.d'


def analysis_document(*, analysis_id='2026-01', title='Title', description='Description', authors=None):
    return {
        'metadata': {
            'id': analysis_id,
            'title': title,
            'description': description,
            'authors': [{'name': name} for name in (authors or ['Jane Doe'])],
        },
    }


def zenodo_document(*, dataset_id='2026-01', analysis_id='2026-01', description='Description', creators=None):
    return {
        'title': f'EOS/DATA-{dataset_id}: Supplementary material for EOS/ANALYSIS-{analysis_id}',
        'description': description,
        'creators': [{'name': name} for name in (creators or ['Jane Doe'])],
    }


class Dataset:
    def __init__(self):
        self._temporary = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary.name)

    def close(self):
        self._temporary.cleanup()

    def yaml(self, relative, document):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(document, sort_keys=False), encoding='utf-8')
        return path

    def text(self, relative, text):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf-8')
        return path

    def json(self, document):
        (self.root / '.zenodo.json').write_text(json.dumps(document), encoding='utf-8')

    def valid(self, **kwargs):
        document = analysis_document(**kwargs)
        self.yaml('analysis.yaml', document)
        metadata = document['metadata']
        self.json(zenodo_document(
            analysis_id=metadata['id'],
            description=metadata['description'],
            creators=[author['name'] for author in metadata['authors']],
        ))

    def run(self, *, files=(), main=None):
        context = CheckContext(self.root, tuple(Path(path) for path in files), Path(main) if main else None)
        return run_checks(create_check_factory(), context, scope=CheckScope.BASIC), context


class DatasetTestCase(unittest.TestCase):
    def setUp(self):
        self.dataset = Dataset()

    def tearDown(self):
        self.dataset.close()


class InfrastructureTests(unittest.TestCase):
    def test_factory_rejects_duplicate_and_unknown_registration(self):
        factory = CheckFactory()
        factory.declare('one', scope=CheckScope.BASIC)
        with self.assertRaises(ValueError):
            factory.declare('one', scope=CheckScope.BASIC)
        with self.assertRaises(KeyError):
            factory.register('unknown', lambda _context: ())
        factory.register('one', lambda _context: ())
        with self.assertRaises(ValueError):
            factory.register('one', lambda _context: ())

    def test_runner_continues_after_check_error_and_reports_infrastructure_failure(self):
        factory = CheckFactory()
        factory.declare('broken', scope=CheckScope.BASIC)
        factory.register('broken', lambda _context: (_ for _ in ()).throw(RuntimeError('boom')))
        factory.declare('independent', scope=CheckScope.BASIC)
        factory.register('independent', lambda _context: [Finding('independent', Severity.INFO, 'ran')])
        result = CheckRunner(factory, scope=CheckScope.BASIC).run(CheckContext(Path('.')))
        self.assertEqual([finding.message for finding in result.findings], ['ran'])
        self.assertEqual(result.exit_status, ExitStatus.USAGE_OR_INFRASTRUCTURE)
        self.assertIn("check 'broken' failed", result.failure)

    def test_runner_rejects_missing_checks_and_non_findings(self):
        missing = CheckFactory()
        missing.declare('missing', scope=CheckScope.BASIC)
        self.assertEqual(CheckRunner(missing, scope=CheckScope.BASIC).run(CheckContext('.')).exit_status, 2)
        invalid = CheckFactory()
        invalid.declare('invalid', scope=CheckScope.BASIC)
        invalid.register('invalid', lambda _context: ['not a finding'])
        self.assertEqual(CheckRunner(invalid, scope=CheckScope.BASIC).run(CheckContext('.')).exit_status, 2)

    def test_factory_orders_prerequisites_and_rejects_invalid_graphs(self):
        calls = []
        factory = CheckFactory()
        factory.declare('dependent', scope=CheckScope.BASIC, needs=('prerequisite',))
        factory.register('dependent', lambda _context: calls.append('dependent') or ())
        factory.declare('prerequisite', scope=CheckScope.BASIC)
        factory.register('prerequisite', lambda _context: calls.append('prerequisite') or ())
        self.assertEqual(CheckRunner(factory, scope=CheckScope.BASIC).run(CheckContext('.')).exit_status, 0)
        self.assertEqual(calls, ['prerequisite', 'dependent'])

        unknown = CheckFactory()
        unknown.declare('check', scope=CheckScope.BASIC, needs=('missing',))
        unknown.register('check', lambda _context: ())
        with self.assertRaisesRegex(ValueError, 'unknown prerequisite'):
            unknown.checks(scope=CheckScope.BASIC)

        cyclic = CheckFactory()
        cyclic.declare('one', scope=CheckScope.BASIC, needs=('two',))
        cyclic.declare('two', scope=CheckScope.BASIC, needs=('one',))
        cyclic.register('one', lambda _context: ())
        cyclic.register('two', lambda _context: ())
        with self.assertRaisesRegex(ValueError, 'cyclic check prerequisites'):
            cyclic.checks(scope=CheckScope.BASIC)

    def test_runner_does_not_run_dependents_after_prerequisite_failure(self):
        calls = []
        factory = CheckFactory()
        factory.declare('broken', scope=CheckScope.BASIC)
        factory.register('broken', lambda _context: 1 / 0)
        factory.declare('dependent', scope=CheckScope.BASIC, needs=('broken',))
        factory.register('dependent', lambda _context: calls.append('dependent') or ())
        factory.declare('independent', scope=CheckScope.BASIC)
        factory.register('independent', lambda _context: calls.append('independent') or ())
        result = CheckRunner(factory, scope=CheckScope.BASIC).run(CheckContext('.'))
        self.assertEqual(result.exit_status, 2)
        self.assertEqual(calls, ['independent'])
        self.assertIn('did not run because prerequisite', result.failure)

    def test_missing_cache_state_is_an_infrastructure_failure(self):
        factory = CheckFactory()
        factory.declare('prerequisite', scope=CheckScope.BASIC)
        factory.register('prerequisite', lambda _context: ())
        factory.declare('metadata', scope=CheckScope.BASIC, needs=('prerequisite',))
        factory.register('metadata', check_analysis_metadata)
        result = CheckRunner(factory, scope=CheckScope.BASIC).run(CheckContext('.'))
        self.assertEqual(result.exit_status, 2)
        self.assertIn('KeyError', result.failure)

    def test_renderer_order_summary_and_exit_outcomes(self):
        findings = (
            Finding('first', Severity.WARNING, 'warning'),
            Finding('second', Severity.INFO, 'information'),
        )
        result = CheckResult(findings)
        rendered = PlainTextRenderer().render(result)
        self.assertLess(rendered.index('[first]'), rendered.index('[second]'))
        self.assertIn('Summary: 0 error(s), 1 warning(s), 1 info finding(s).', rendered)
        self.assertEqual(result.exit_status, 0)
        self.assertEqual(CheckResult((Finding('bad', Severity.ERROR, 'error'),)).exit_status, 1)
        self.assertEqual(CheckResult(failure='failed').exit_status, 2)

    def test_name_normalization(self):
        self.assertEqual(_normalize_name('  DOE,  Jane-Mary '), 'jane mary doe')
        self.assertEqual(_normalize_name('Ａlice O\N{COMBINING DIAERESIS}zil'), _normalize_name('Alice \N{LATIN SMALL LETTER O WITH DIAERESIS}zil'))
        self.assertNotEqual(_normalize_name('Jose'), _normalize_name('Jos\N{LATIN SMALL LETTER E WITH ACUTE}'))


class SelectionTests(DatasetTestCase):
    def test_missing_root_analysis(self):
        result, _ = self.dataset.run()
        self.assertEqual(result.exit_status, 1)
        self.assertIn('supply every analysis file explicitly', result.findings[0].message)

    def test_valid_sole_root_is_automatically_main_and_does_not_mutate(self):
        source = FIXTURES / 'valid'
        before = {path.name: path.read_bytes() for path in source.iterdir()}
        context = CheckContext(source)
        first = run_checks(create_check_factory(), context, scope=CheckScope.BASIC)
        second = run_checks(create_check_factory(), CheckContext(source), scope=CheckScope.BASIC)
        self.assertEqual(first.exit_status, 0)
        self.assertEqual(context.main_analysis_path, (source / 'analysis.yaml').resolve())
        self.assertEqual(first.findings, second.findings)
        self.assertEqual(before, {path.name: path.read_bytes() for path in source.iterdir()})

    def test_invalid_findings_have_stable_check_order_and_summary_counts(self):
        self.dataset.text('analysis.yaml', 'metadata: [')
        self.dataset.text('.zenodo.json', '{')
        first, _ = self.dataset.run()
        second, _ = self.dataset.run()
        self.assertEqual(first.findings, second.findings)
        self.assertEqual(
            [finding.check_id for finding in first.findings],
            [
                'analysis-metadata.parse',
                'analysis-metadata.consistency',
                'zenodo-metadata.basic',
                'authors.creator-consistency',
            ],
        )
        rendered = PlainTextRenderer().render(first)
        self.assertIn('Summary: 2 error(s), 0 warning(s), 2 info finding(s).', rendered)

    def test_explicit_one_and_multiple_file_selection(self):
        self.dataset.yaml('one.yaml', analysis_document())
        self.dataset.yaml('nested/two.yml', analysis_document())
        self.dataset.json(zenodo_document())
        one, one_context = self.dataset.run(files=['one.yaml'])
        many, many_context = self.dataset.run(files=['one.yaml', 'nested/two.yml'], main='nested/two.yml')
        self.assertEqual(one.exit_status, 0)
        self.assertEqual(one_context.main_analysis_path, (self.dataset.root / 'one.yaml').resolve())
        self.assertEqual(many.exit_status, 0)
        self.assertEqual(many_context.main_analysis_path, (self.dataset.root / 'nested/two.yml').resolve())

    def test_multiple_requires_included_main(self):
        self.dataset.yaml('one.yaml', analysis_document())
        self.dataset.yaml('two.yaml', analysis_document())
        self.dataset.json(zenodo_document())
        missing, _ = self.dataset.run(files=['one.yaml', 'two.yaml'])
        outside, _ = self.dataset.run(files=['one.yaml', 'two.yaml'], main='three.yaml')
        self.assertTrue(any('require --main-analysis-file' in f.message for f in missing.findings))
        self.assertTrue(any('must identify one' in f.message for f in outside.findings))

    def test_duplicate_outside_non_file_and_wrong_extension(self):
        self.dataset.valid()
        self.dataset.text('directory.yaml/child', 'x')
        self.dataset.text('analysis.txt', 'x')
        outside = self.dataset.root.parent / 'outside.yaml'
        cases = [
            (['analysis.yaml', './analysis.yaml'], 'duplicated'),
            ([outside], 'outside'),
            (['directory.yaml'], 'regular file'),
            (['analysis.txt'], '.yaml or .yml'),
        ]
        for paths, message in cases:
            with self.subTest(paths=paths):
                result, _ = self.dataset.run(files=paths)
                self.assertEqual(result.exit_status, 1)
                self.assertTrue(any(message in finding.message for finding in result.findings))

    def test_conflicting_main_for_single_file(self):
        self.dataset.valid()
        self.dataset.yaml('other.yaml', analysis_document())
        result, _ = self.dataset.run(files=['analysis.yaml'], main='other.yaml')
        self.assertTrue(any('conflicts' in finding.message for finding in result.findings))


class AnalysisMetadataTests(DatasetTestCase):
    def test_empty_malformed_and_non_mapping_yaml(self):
        self.dataset.json(zenodo_document())
        for text, message in [('', 'top-level mapping'), ('metadata: [', 'Cannot parse YAML'), ('- item\n', 'top-level mapping')]:
            with self.subTest(text=text):
                self.dataset.text('analysis.yaml', text)
                result, _ = self.dataset.run()
                self.assertTrue(any(message in finding.message for finding in result.findings))
                self.assertTrue(any('Skipped' in finding.message for finding in result.findings))

    def test_absent_and_malformed_metadata_fields(self):
        self.dataset.json(zenodo_document())
        documents = [
            ({}, 'metadata must be a mapping'),
            ({'metadata': {'authors': [{'name': 'Jane Doe'}]}}, 'metadata.id'),
            ({'metadata': {'id': '2026-01', 'authors': 'Jane Doe'}}, 'authors must be a list'),
            ({'metadata': {'id': '2026-01', 'authors': [{}]}}, 'authors[0].name'),
            (analysis_document(title=''), 'metadata.title'),
            (analysis_document(description=' '), 'metadata.description'),
        ]
        for document, message in documents:
            with self.subTest(message=message):
                self.dataset.yaml('analysis.yaml', document)
                result, _ = self.dataset.run()
                self.assertEqual(result.exit_status, 1)
                self.assertTrue(any(message in finding.message for finding in result.findings))

    def test_author_entries_match_runtime_metadata_schema(self):
        self.dataset.json(zenodo_document())
        invalid_authors = [
            ({'name': 'Jane Doe', 'orcid': '0000-0000-0000-0000'}, "Unknown key 'orcid'"),
            ({'name': 'Jane Doe', 'affiliation': 7}, 'affiliation must be a string'),
            ({'name': 'Jane Doe', 'email': []}, 'email must be a string'),
        ]
        for author, expected in invalid_authors:
            with self.subTest(author=author):
                document = analysis_document()
                document['metadata']['authors'] = [author]
                self.dataset.yaml('analysis.yaml', document)
                result, _ = self.dataset.run()
                self.assertEqual(result.exit_status, 1)
                self.assertTrue(any(expected in finding.message for finding in result.findings))

    def test_analysis_id_boundaries_and_revisions(self):
        for analysis_id, valid in [('0000-00', True), ('9999-99', True), ('2026-01v2', False), ('26-01', False), ('2026-1', False)]:
            with self.subTest(analysis_id=analysis_id):
                self.dataset.valid(analysis_id=analysis_id)
                result, _ = self.dataset.run()
                self.assertEqual(result.exit_status == 0, valid)

    def test_inconsistent_ids_and_author_sets(self):
        self.dataset.yaml('one.yaml', analysis_document(authors=['Jane Doe', 'John Roe']))
        self.dataset.yaml('two.yaml', analysis_document(analysis_id='2026-02', authors=['Jane Doe']))
        self.dataset.json(zenodo_document(creators=['Jane Doe', 'John Roe']))
        result, _ = self.dataset.run(files=['one.yaml', 'two.yaml'], main='one.yaml')
        messages = [finding.message for finding in result.findings]
        self.assertTrue(any('differs from common ID' in message for message in messages))
        self.assertTrue(any('author set differs' in message for message in messages))

    def test_author_order_is_ignored_but_duplicates_and_ambiguity_are_errors(self):
        self.dataset.json(zenodo_document(creators=['Jane Doe', 'John Roe']))
        self.dataset.yaml('one.yaml', analysis_document(authors=['Jane Doe', 'John Roe']))
        self.dataset.yaml('two.yaml', analysis_document(authors=['John Roe', 'Jane Doe']))
        result, _ = self.dataset.run(files=['one.yaml', 'two.yaml'], main='one.yaml')
        self.assertEqual(result.exit_status, 0)
        self.dataset.yaml('two.yaml', analysis_document(authors=['Jane Doe', 'Jane-Doe']))
        duplicate, _ = self.dataset.run(files=['one.yaml', 'two.yaml'], main='one.yaml')
        self.assertEqual(duplicate.exit_status, 1)
        self.dataset.yaml('two.yaml', analysis_document(authors=['J. Doe', 'John Roe']))
        ambiguous, _ = self.dataset.run(files=['one.yaml', 'two.yaml'], main='one.yaml')
        self.assertEqual(ambiguous.exit_status, 1)


class ZenodoTests(DatasetTestCase):
    def test_missing_malformed_and_non_object_json(self):
        self.dataset.yaml('analysis.yaml', analysis_document())
        missing, _ = self.dataset.run()
        self.assertTrue(any('missing' in finding.message for finding in missing.findings))
        self.dataset.text('.zenodo.json', '{')
        malformed, _ = self.dataset.run()
        self.assertTrue(any('Cannot parse JSON' in finding.message for finding in malformed.findings))
        self.dataset.json([])
        non_object, _ = self.dataset.run()
        self.assertTrue(any('JSON object' in finding.message for finding in non_object.findings))

    def test_initial_and_revision_dataset_ids(self):
        self.dataset.yaml('analysis.yaml', analysis_document())
        for dataset_id, valid in [('2026-01', True), ('2026-01v2', True), ('2026-01v19', True), ('2026-01v1', False), ('2026-1', False)]:
            with self.subTest(dataset_id=dataset_id):
                self.dataset.json(zenodo_document(dataset_id=dataset_id))
                result, _ = self.dataset.run()
                self.assertEqual(result.exit_status == 0, valid)

    def test_title_wrapper_and_embedded_analysis_id(self):
        self.dataset.valid()
        document = zenodo_document()
        document['title'] = 'EOS/DATA-2026-01 - Supplementary material for EOS/ANALYSIS-2026-01'
        self.dataset.json(document)
        wrapper, _ = self.dataset.run()
        self.assertTrue(any('exactly match' in finding.message for finding in wrapper.findings))
        self.dataset.json(zenodo_document(analysis_id='2026-02'))
        mismatch, _ = self.dataset.run()
        self.assertTrue(any("expected '2026-01'" in finding.message for finding in mismatch.findings))

    def test_description_equality_is_exact(self):
        self.dataset.valid(description='Line one\nLine two')
        equal, _ = self.dataset.run()
        self.assertEqual(equal.exit_status, 0)
        self.dataset.json(zenodo_document(description='Line one\nLine two '))
        different, _ = self.dataset.run()
        self.assertTrue(any('does not exactly equal' in finding.message for finding in different.findings))

    def test_missing_or_invalid_creators(self):
        self.dataset.yaml('analysis.yaml', analysis_document())
        for creators in [None, [], [{}], [{'name': '---'}]]:
            document = zenodo_document()
            document['creators'] = creators
            self.dataset.json(document)
            result, _ = self.dataset.run()
            self.assertTrue(any('creators' in finding.message for finding in result.findings))


class CreatorConsistencyTests(DatasetTestCase):
    def test_normalized_names_are_reused_between_checks(self):
        self.dataset.valid()
        with mock.patch(
            'eos.cli.data_checks._normalize_name',
            wraps=_normalize_name,
        ) as normalize:
            result, _ = self.dataset.run()
        self.assertEqual(result.exit_status, 0)
        self.assertEqual(normalize.call_count, 2)

    def test_exact_family_given_case_whitespace_punctuation_and_unicode(self):
        authors = ['Jane Mary Doe', 'Alice \N{LATIN SMALL LETTER O WITH DIAERESIS}zil']
        self.dataset.valid(authors=authors)
        self.dataset.json(zenodo_document(creators=['DOE, jane-mary', '  Alice  O\N{COMBINING DIAERESIS}zil ']))
        result, _ = self.dataset.run()
        self.assertEqual(result.exit_status, 0)

    def test_diacritics_initials_order_and_compound_surnames_are_ambiguous(self):
        cases = [
            ('Jos\N{LATIN SMALL LETTER E WITH ACUTE} Perez', 'Jose Perez'),
            ('Jane Doe', 'J Doe'),
            ('Jane Mary Doe', 'Doe Jane Mary'),
            ('Jean Van Dyke', 'J Van Dyke'),
        ]
        for author, creator in cases:
            with self.subTest(author=author, creator=creator):
                self.dataset.valid(authors=[author])
                self.dataset.json(zenodo_document(creators=[creator]))
                result, _ = self.dataset.run()
                warnings = [finding for finding in result.findings if finding.severity is Severity.WARNING]
                self.assertEqual(result.exit_status, 0)
                self.assertTrue(warnings)
                self.assertIsNotNone(warnings[0].question)

    def test_non_unique_ambiguous_candidates_remain_warnings(self):
        self.dataset.valid(authors=['Jane Doe', 'Janet Doe'])
        self.dataset.json(zenodo_document(creators=['J Doe']))
        result, _ = self.dataset.run()
        warnings = [finding for finding in result.findings if finding.severity is Severity.WARNING]
        self.assertEqual(len(warnings), 2)
        self.assertEqual(result.exit_status, 0)

    def test_definite_missing_and_additional_creators(self):
        self.dataset.valid(authors=['Jane Doe', 'John Roe'])
        self.dataset.json(zenodo_document(creators=['Jane Doe', 'Alice Poe']))
        result, _ = self.dataset.run()
        errors = [finding.message for finding in result.findings if finding.severity is Severity.ERROR]
        self.assertTrue(any('missing from Zenodo' in message for message in errors))
        self.assertTrue(any('not an analysis author' in message for message in errors))


if __name__ == '__main__':
    unittest.main()
