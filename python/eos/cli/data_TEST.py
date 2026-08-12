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
import os
from pathlib import Path
import sys
import unittest
from unittest import mock

import eos

from eos.cli import data
from eos.cli.data_checks import CheckResult, Finding, Severity
from eos.cli.data_release import SourceCheckError


FIXTURE = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/data_checks_TEST.d/valid'


class FakeDataSets:
    def __init__(self):
        self.downloaded = []
        self.updated = False

    def download(self, dataset_id):
        self.downloaded.append(dataset_id)

    def update(self):
        self.updated = True

    def datasets(self):
        item = type('DataSet', (), {'title': 'A title', 'authors': ['Jane Doe']})()
        return [('2026-01', item)]


class ParserCompatibilityTests(unittest.TestCase):
    def test_legacy_download_list_and_update_dispatch(self):
        fake = FakeDataSets()
        with mock.patch.object(data, '_datasets', return_value=fake):
            self.assertEqual(data.main(['download', '2026-01']), 0)
            output = io.StringIO()
            self.assertEqual(data.main(['list'], stdout=output), 0)
            self.assertEqual(data.main(['update']), 0)
        self.assertEqual(fake.downloaded, ['2026-01'])
        self.assertTrue(fake.updated)
        self.assertIn('2026-01', output.getvalue())

    def test_check_defaults_and_repeatable_options(self):
        args = data._parser().parse_args(['check'])
        self.assertEqual(args.directory, '.')
        self.assertIsNone(args.analysis_file)
        args = data._parser().parse_args([
            'check', '--analysis-file', 'one.yaml', '--analysis-file', 'two.yaml',
            '--main-analysis-file', 'two.yaml', 'dataset',
        ])
        self.assertEqual(args.analysis_file, ['one.yaml', 'two.yaml'])
        self.assertEqual(args.main_analysis_file, 'two.yaml')
        self.assertEqual(args.directory, 'dataset')

    def test_create_and_publish_interfaces(self):
        args = data._parser().parse_args([
            'create', '2026-01', '/analysis', '--analysis-file', 'one.yaml',
            '--analysis-file', 'two.yaml', '--main-analysis-file', 'two.yaml',
            '--revision', '--dry-run',
        ])
        self.assertEqual((args.base_id, args.source), ('2026-01', '/analysis'))
        self.assertEqual(args.analysis_file, ['one.yaml', 'two.yaml'])
        self.assertTrue(args.revision)
        self.assertTrue(args.dry_run)
        args = data._parser().parse_args(['publish', '2026-01v2', '--dry-run'])
        self.assertEqual(args.data_id, '2026-01v2')
        self.assertTrue(args.dry_run)

        error = io.StringIO()
        self.assertEqual(
            data.main(['create', '2026-01', '/analysis', '--revision', '--replace'], stderr=error),
            2,
        )

    def test_interactive_aliases_are_rejected(self):
        for option in ('-i', '--interactive'):
            with self.subTest(option=option):
                error = io.StringIO()
                self.assertEqual(data.main(['check', option], stderr=error), 2)
                self.assertIn('not yet available', error.getvalue())
        for option in ('-i', '--interactive'):
            with self.subTest(command='create', option=option):
                error = io.StringIO()
                self.assertEqual(
                    data.main(['create', '2026-01', '/analysis', option], stderr=error), 2,
                )
                self.assertIn('not yet available', error.getvalue())

    def test_help_and_missing_command(self):
        output = io.StringIO()
        self.assertEqual(data.main(['--help'], stdout=output), 0)
        self.assertIn('check', output.getvalue())
        output = io.StringIO()
        self.assertEqual(data.main([], stdout=output), 0)
        self.assertIn('usage:', output.getvalue())

    def test_default_streams_are_resolved_at_call_time(self):
        output = io.StringIO()
        with mock.patch.object(sys, 'stdout', output):
            self.assertEqual(data.main(['--help']), 0)
        self.assertIn('usage:', output.getvalue())

        error = io.StringIO()
        with mock.patch.object(sys, 'stderr', error):
            self.assertEqual(data.main(['check', '--interactive']), 2)
        self.assertIn('not yet available', error.getvalue())

    def test_logging_uses_injected_stderr_and_reconfigures_each_call(self):
        fake = FakeDataSets()

        def update():
            eos.info('update log message')
            fake.updated = True

        fake.update = update
        verbose_error = io.StringIO()
        quiet_error = io.StringIO()
        with mock.patch.object(data, '_datasets', return_value=fake):
            self.assertEqual(data.main(['update', '-vv'], stderr=verbose_error), 0)
            self.assertEqual(data.main(['update'], stderr=quiet_error), 0)
        self.assertIn('update log message', verbose_error.getvalue())
        self.assertNotIn('update log message', quiet_error.getvalue())

    def test_wheel_package_configuration_includes_cli(self):
        project = Path(os.environ['SOURCE_DIR']).parent / 'pyproject.toml.in'
        package_line = next(
            line for line in project.read_text(encoding='utf-8').splitlines()
            if line.startswith('packages = ')
        )
        self.assertIn('"eos.cli"', package_line)

    def test_main_propagates_check_exit_status(self):
        for fixture, expected in ((FIXTURE, 0), (FIXTURE / 'missing', 1)):
            with self.subTest(expected=expected):
                self.assertEqual(data.main(['check', str(fixture)], stdout=io.StringIO()), expected)
        incomplete_factory = data.CheckFactory()
        incomplete_factory.declare('not-registered', scope=data.CheckScope.BASIC)
        self.assertEqual(
            data.main(
                ['check', str(FIXTURE)],
                stdout=io.StringIO(),
                check_factory=incomplete_factory,
            ),
            2,
        )

    def test_create_preserves_source_check_exit_status(self):
        validation = SourceCheckError(CheckResult(findings=(
            Finding('analysis.metadata', Severity.ERROR, 'invalid'),
        )))
        infrastructure = SourceCheckError(CheckResult(failure='runner unavailable'))
        for error, expected in ((validation, 1), (infrastructure, 2)):
            with self.subTest(expected=expected), mock.patch.object(
                data, 'plan_create', side_effect=error,
            ):
                self.assertEqual(
                    data.main(
                        ['create', '2026-01', '/analysis'],
                        stdout=io.StringIO(), stderr=io.StringIO(),
                    ),
                    expected,
                )


class IsolationTests(unittest.TestCase):
    def test_noninteractive_check_imports_neither_rich_nor_gh_and_uses_no_network(self):
        removed = {name: module for name, module in sys.modules.items() if name == 'rich' or name.startswith('rich.')}
        for name in removed:
            sys.modules.pop(name)
        old_path = os.environ.get('PATH')
        try:
            os.environ['PATH'] = ''
            with mock.patch('socket.create_connection', side_effect=AssertionError('network access')):
                self.assertEqual(data.main(['check', str(FIXTURE)], stdout=io.StringIO()), 0)
            self.assertFalse(any(name == 'rich' or name.startswith('rich.') for name in sys.modules))
        finally:
            if old_path is None:
                os.environ.pop('PATH', None)
            else:
                os.environ['PATH'] = old_path
            sys.modules.update(removed)


if __name__ == '__main__':
    unittest.main()
