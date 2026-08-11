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

import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import yaml

from eos.cli.data import create_check_factory
from eos.cli.data_checks import CacheKey, CheckContext, CheckScope, Severity, run_checks
from eos.cli.data_checks_dataset import (
    FileInventoryEntry,
    MAX_DATASET_SIZE,
    MAX_FILE_SIZE,
    _valid_orcid,
    check_dataset_sizes,
)
from eos.tasks import task_output_templates


def analysis_document(*, figures=None, steps=None):
    result = {
        'metadata': {
            'id': '2026-01',
            'title': 'Title',
            'description': 'Description',
            'authors': [{'name': 'Jane Doe'}],
        },
    }
    if figures is not None:
        result['figures'] = figures
    if steps is not None:
        result['steps'] = steps
    return result


def zenodo_document(**updates):
    result = {
        'title': 'EOS/DATA-2026-01: Supplementary material for EOS/ANALYSIS-2026-01',
        'description': 'Description',
        'upload_type': 'dataset',
        'license': 'CC-BY-4.0',
        'creators': [{'name': 'Jane Doe', 'affiliation': 'Example University'}],
        'grants': [],
    }
    result.update(updates)
    return result


class Dataset:
    def __init__(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def close(self):
        self.temporary.cleanup()

    def write(self, relative, content=''):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding='utf-8')
        return path

    def yaml(self, relative, document):
        return self.write(relative, yaml.safe_dump(document, sort_keys=False))

    def valid(self, *, analysis=None, zenodo=None, readme='# Dataset\n'):
        self.yaml('analysis.yaml', analysis or analysis_document())
        self.write('.zenodo.json', json.dumps(zenodo or zenodo_document()))
        self.write('README.md', readme)

    def run(self, *, paths=(), main=None):
        context = CheckContext(
            self.root,
            tuple(Path(path) for path in paths),
            Path(main) if main else None,
        )
        return run_checks(create_check_factory(), context, scope=CheckScope.COMPLETE), context


class DatasetCheckTestCase(unittest.TestCase):
    def setUp(self):
        self.dataset = Dataset()

    def tearDown(self):
        self.dataset.close()


class ZenodoCompletenessTests(DatasetCheckTestCase):
    def test_complete_metadata_warns_only_for_unconfirmed_grants(self):
        self.dataset.valid(zenodo=zenodo_document(grants=[{'id': 'EU-123'}], extra='kept'))
        result, _ = self.dataset.run()
        warnings = [finding for finding in result.findings if finding.check_id == 'zenodo-metadata.completeness']
        self.assertEqual(result.exit_status, 0)
        self.assertEqual(len(warnings), 1)
        self.assertEqual(warnings[0].severity, Severity.WARNING)
        self.assertEqual(warnings[0].details['grant_ids'], ('EU-123',))

    def test_required_fields_affiliation_and_grant_structure(self):
        documents = (
            (zenodo_document(upload_type='publication'), 'upload_type'),
            (zenodo_document(license=' '), 'license'),
            (zenodo_document(creators=[{'name': 'Jane Doe'}]), 'affiliation'),
            (zenodo_document(grants={}), 'grants must be a list'),
            (zenodo_document(grants=[{}]), 'grants[0].id'),
        )
        for document, expected in documents:
            with self.subTest(expected=expected):
                self.dataset.valid(zenodo=document)
                result, _ = self.dataset.run()
                self.assertEqual(result.exit_status, 1)
                self.assertTrue(any(expected in finding.message for finding in result.findings))

    def test_orcid_shape_checksum_and_final_x(self):
        self.assertTrue(_valid_orcid('0000-0002-1825-0097'))
        self.assertTrue(_valid_orcid('0000-0002-1694-233X'))
        self.assertFalse(_valid_orcid('0000-0002-1825-0098'))
        self.assertFalse(_valid_orcid('0000000218250097'))
        for orcid, valid in (
            ('0000-0002-1825-0097', True),
            ('0000-0002-1694-233X', True),
            ('0000-0002-1825-0098', False),
            ('bad', False),
        ):
            document = zenodo_document(creators=[{
                'name': 'Jane Doe', 'affiliation': 'Institute', 'orcid': orcid,
            }])
            self.dataset.valid(zenodo=document)
            result, _ = self.dataset.run()
            self.assertEqual(not any('orcid' in finding.message for finding in result.findings), valid)


class ReadmeFigureTests(DatasetCheckTestCase):
    def test_markdown_html_nested_query_fragment_and_external_references(self):
        readme = '''
![one](figures/nested/one.pdf?download=1)
<img src="/figures/nested/one.png#preview">
[external](https://example.test/figures/ignored.pdf)
'''
        self.dataset.valid(readme=readme)
        self.dataset.write('figures/nested/one.pdf')
        self.dataset.write('figures/nested/one.png')
        result, _ = self.dataset.run()
        messages = [finding.message for finding in result.findings if finding.check_id == 'figures.readme']
        self.assertFalse(any('not listed' in message for message in messages))
        self.assertFalse(any('redundantly' in message for message in messages))

    def test_missing_readme_pdf_png_other_format_duplicates_and_traversal(self):
        self.dataset.valid(readme='![one](figures/one.pdf)\n![again](figures/one.pdf)\n[x](../outside.pdf)')
        self.dataset.write('figures/one.png')
        self.dataset.write('figures/one.svg')
        result, _ = self.dataset.run()
        messages = [finding.message for finding in result.findings if finding.check_id == 'figures.readme']
        for expected in ('no PDF', 'additional format', 'redundantly', 'traverses outside'):
            self.assertTrue(any(expected in message for message in messages), expected)
        (self.dataset.root / 'README.md').unlink()
        missing, _ = self.dataset.run()
        self.assertTrue(any('README.md is missing' in finding.message for finding in missing.findings))

    def test_unlisted_and_missing_png(self):
        self.dataset.valid()
        self.dataset.write('figures/one.pdf')
        result, _ = self.dataset.run()
        messages = [finding.message for finding in result.findings if finding.check_id == 'figures.readme']
        self.assertTrue(any('no PNG' in message for message in messages))
        self.assertTrue(any('not listed' in message for message in messages))

    def test_declared_figures_across_selected_analysis_files(self):
        one = analysis_document(figures=[{'name': 'one'}])
        two = analysis_document(figures=[{'name': 'nested/two'}])
        self.dataset.valid(analysis=one, readme='![one](figures/one.pdf)\n![two](figures/nested/two.pdf)')
        self.dataset.yaml('two.yaml', two)
        self.dataset.write('figures/one.pdf')
        self.dataset.write('figures/one.png')
        self.dataset.write('figures/nested/two.pdf')
        result, _ = self.dataset.run(paths=('analysis.yaml', 'two.yaml'), main='analysis.yaml')
        findings = [finding for finding in result.findings if finding.check_id == 'figures.analysis-declarations']
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].severity, Severity.WARNING)
        self.assertEqual(findings[0].details['analysis_path'], 'two.yaml')
        self.assertEqual(findings[0].details['output_path'], 'figures/nested/two.png')


class SizeTests(unittest.TestCase):
    def result(self, sizes):
        context = CheckContext(Path('.'))
        context.cache[CacheKey.FILE_INVENTORY] = tuple(
            FileInventoryEntry(Path(f'file-{index}'), size)
            for index, size in enumerate(sizes)
        )
        return tuple(check_dataset_sizes(context))

    def test_file_threshold_boundaries(self):
        for size, error in ((MAX_FILE_SIZE - 1, False), (MAX_FILE_SIZE, False), (MAX_FILE_SIZE + 1, True)):
            with self.subTest(size=size):
                findings = self.result((size,))
                self.assertEqual(any(finding.severity is Severity.ERROR for finding in findings), error)

    def test_total_threshold_boundaries_and_exact_info(self):
        for total, warning in ((MAX_DATASET_SIZE - 1, False), (MAX_DATASET_SIZE, False), (MAX_DATASET_SIZE + 1, True)):
            with self.subTest(total=total):
                findings = self.result((MAX_FILE_SIZE,) * (total // MAX_FILE_SIZE) + (total % MAX_FILE_SIZE,))
                total_finding = next(finding for finding in findings if finding.severity is Severity.INFO)
                self.assertEqual(total_finding.details['total_bytes'], total)
                oversized = [finding for finding in findings if finding.severity is Severity.WARNING]
                self.assertEqual(bool(oversized), warning)
                if oversized:
                    self.assertIn('importance samples and posterior predictions', oversized[0].message)


class ReproducibleOutputTests(DatasetCheckTestCase):
    def output(self, relative, output_type='MarkovChain'):
        self.dataset.yaml(f'{relative}/description.yaml', {'type': output_type})

    def step(self, task='sample-mcmc', arguments=None, defaults=None, step_id='sample'):
        return {
            'title': 'Sample', 'id': step_id,
            'default_arguments': defaults or {},
            'tasks': [{'task': task, 'arguments': arguments or {}}],
        }

    def test_task_templates_are_introspected_from_decorator_registry(self):
        templates = task_output_templates()
        self.assertEqual(templates['sample-mcmc'], 'data/{posterior}/mcmc-{chain:04}')
        self.assertEqual(templates['predict-observables'], 'data/{posterior}/pred-{prediction}')

    def test_valid_claimed_output_and_step_defaults(self):
        step = self.step(arguments={'posterior': 'p'}, defaults={'sample-mcmc': {'chain': 7}})
        self.dataset.valid(analysis=analysis_document(steps=[step]))
        self.output('data/p/mcmc-0007')
        result, _ = self.dataset.run()
        findings = [finding for finding in result.findings if finding.check_id == 'outputs.reproducibility']
        self.assertEqual(findings, [])

    def test_missing_expected_unclaimed_and_missing_steps(self):
        self.dataset.valid(analysis=analysis_document(steps=[self.step(arguments={'posterior': 'p', 'chain': 1})]))
        self.output('data/p/mcmc-0002')
        result, _ = self.dataset.run()
        messages = [finding.message for finding in result.findings if finding.check_id == 'outputs.reproducibility']
        self.assertTrue(any('Expected task output is absent' in message for message in messages))
        self.assertTrue(any('not claimed' in message for message in messages))
        self.dataset.valid(analysis=analysis_document())
        missing_steps, _ = self.dataset.run()
        self.assertTrue(any('not claimed' in finding.message for finding in missing_steps.findings))

    def test_unexpandable_template_and_cross_analysis_collision(self):
        self.dataset.valid(analysis=analysis_document(steps=[self.step(arguments={'posterior': 'p'})]))
        invalid, _ = self.dataset.run()
        self.assertTrue(any('Cannot expand output template' in finding.message for finding in invalid.findings))

        with mock.patch('eos.tasks.task_output_templates', return_value={'sample-mcmc': 'data/{posterior!bad}'}):
            malformed, _ = self.dataset.run()
        self.assertTrue(any('Cannot expand output template' in finding.message for finding in malformed.findings))

        shared_step = self.step(arguments={'posterior': 'p', 'chain': 1})
        self.dataset.yaml('analysis.yaml', analysis_document(steps=[shared_step]))
        self.dataset.yaml('two.yaml', analysis_document(steps=[shared_step]))
        self.output('data/p/mcmc-0001')
        collision, _ = self.dataset.run(paths=('analysis.yaml', 'two.yaml'), main='analysis.yaml')
        finding = next(finding for finding in collision.findings if 'claim the same output' in finding.message)
        self.assertEqual(finding.details['output_path'], 'data/p/mcmc-0001')

    def test_legacy_argument_alias_and_description_without_type(self):
        step = self.step(arguments={'posterior': 'p', 'CHAIN-IDX': 3})
        self.dataset.valid(analysis=analysis_document(steps=[step]))
        self.dataset.yaml('data/p/mcmc-0003/description.yaml', {'parameters': []})
        result, _ = self.dataset.run()
        findings = [finding for finding in result.findings if finding.check_id == 'outputs.reproducibility']
        self.assertEqual(findings, [])

    def test_step_default_argument_alias_is_normalized(self):
        step = self.step(
            arguments={'posterior': 'p'},
            defaults={'sample-mcmc': {'CHAIN-IDX': 4}},
        )
        self.dataset.valid(analysis=analysis_document(steps=[step]))
        self.dataset.yaml('data/p/mcmc-0004/description.yaml', {'type': 'MarkovChain'})
        result, _ = self.dataset.run()
        findings = [
            finding for finding in result.findings
            if finding.check_id == 'outputs.reproducibility'
        ]
        self.assertEqual(findings, [])


class IntegrationTests(DatasetCheckTestCase):
    def test_analysis_dependent_checks_skip_when_any_selected_file_is_invalid(self):
        self.dataset.valid()
        self.dataset.write('broken.yaml', 'metadata: [')
        result, _ = self.dataset.run(
            paths=('analysis.yaml', 'broken.yaml'), main='analysis.yaml',
        )
        self.assertEqual(result.exit_status, 1)
        self.assertIsNone(result.failure)
        skipped = {
            finding.check_id: finding
            for finding in result.findings
            if finding.check_id in (
                'figures.analysis-declarations', 'outputs.reproducibility',
            )
        }
        self.assertEqual(set(skipped), {
            'figures.analysis-declarations', 'outputs.reproducibility',
        })
        self.assertTrue(all(finding.severity is Severity.INFO for finding in skipped.values()))

        inconsistent = analysis_document()
        inconsistent['metadata']['authors'] = [{'name': 'John Roe'}]
        self.dataset.yaml('broken.yaml', inconsistent)
        result, _ = self.dataset.run(
            paths=('analysis.yaml', 'broken.yaml'), main='analysis.yaml',
        )
        self.assertEqual(result.exit_status, 1)
        self.assertIsNone(result.failure)
        skipped_ids = {
            finding.check_id for finding in result.findings
            if 'not all selected analysis files passed validation' in finding.message
        }
        self.assertEqual(skipped_ids, {
            'figures.analysis-declarations', 'outputs.reproducibility',
        })

    def test_registered_order_warning_exit_and_no_mutation(self):
        self.dataset.valid()
        before = {path.relative_to(self.dataset.root): path.read_bytes() for path in self.dataset.root.rglob('*') if path.is_file()}
        result, _ = self.dataset.run()
        after = {path.relative_to(self.dataset.root): path.read_bytes() for path in self.dataset.root.rglob('*') if path.is_file()}
        self.assertEqual(result.exit_status, 0)
        self.assertEqual(before, after)
        ids = [registered.declaration.name for registered in create_check_factory().checks(scope=CheckScope.COMPLETE)]
        self.assertLess(ids.index('zenodo-metadata.completeness'), ids.index('figures.readme'))
        self.assertLess(ids.index('figures.readme'), ids.index('dataset.size'))
        self.assertTrue(any(finding.severity is Severity.WARNING for finding in result.findings))

    def test_inventory_does_not_follow_symlinks_or_count_git(self):
        self.dataset.valid()
        self.dataset.write('.git/large', 'x')
        outside = self.dataset.root.parent / 'outside-dataset-checks'
        outside.write_text('secret', encoding='utf-8')
        try:
            (self.dataset.root / 'linked').symlink_to(outside)
            result, _ = self.dataset.run()
            total = next(finding.details['total_bytes'] for finding in result.findings if finding.check_id == 'dataset.size' and finding.severity is Severity.INFO)
            expected = sum(path.stat().st_size for path in self.dataset.root.iterdir() if path.is_file() and not path.is_symlink())
            self.assertEqual(total, expected)
        finally:
            outside.unlink()


if __name__ == '__main__':
    unittest.main()
