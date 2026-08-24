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

# These tests cover the 'eos-analysis run' command, which carries out the tasks of a single
# step recorded in an analysis file.

import io
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import eos

from eos.cli import analysis


_FIXTURES = Path(os.environ.get('SOURCE_DIR', Path(__file__).parents[2])) / 'eos/cli/analysis_TEST.d'

# 'steps.yaml' holds steps that are cheap enough to be run for real; 'analysis.yaml' holds the
# steps of a realistic analysis, which are only ever dry-run here.
STEPS_FILE    = str(_FIXTURES / 'steps.yaml')
ANALYSIS_FILE = str(_FIXTURES / 'analysis.yaml')

POSTERIOR = 'FF'

_BASE = None


def run(*argv, base, analysis_file=STEPS_FILE):
    """Runs one command line and returns its exit status together with the captured output."""
    stdout, stderr = io.StringIO(), io.StringIO()
    status = analysis.main([*argv, '-f', analysis_file, '-b', base], stdout=stdout, stderr=stderr)
    return status, stdout.getvalue(), stderr.getvalue()


def setUpModule():
    """Runs both steps of 'steps.yaml' once, in their intended order."""
    global _BASE
    _BASE = tempfile.mkdtemp()
    for step in ('FF.sample', 'FF.mode+predict'):
        status, stdout, _ = run('run', step, base=_BASE)
        assert status == 0, f'failed to run step \'{step}\': {stdout}'


def tearDownModule():
    shutil.rmtree(_BASE, ignore_errors=True)


def data_path(*parts, base=None):
    return os.path.join(_BASE if base is None else base, 'data', POSTERIOR, *parts)


class DryRunTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # The 'run' task writes its log file to the base directory even for a dry run.
        cls._base = tempfile.mkdtemp()
        cls.addClassCleanup(shutil.rmtree, cls._base, ignore_errors=True)

    def dry_run(self, step, analysis_file=STEPS_FILE):
        "Returns the task invocations that running the given step would perform, one per line."
        status, stdout, _ = self.run_step(step, analysis_file)

        self.assertEqual(status, 0)
        return stdout.strip().splitlines()

    def run_step(self, step, analysis_file=STEPS_FILE):
        return run('run', '-d', step, base=self._base, analysis_file=analysis_file)

    def test_single_task_step(self):
        "A step with a single task reports that one task."

        lines = self.dry_run('FF.sample')

        self.assertEqual(len(lines), 1)
        self.assertTrue(lines[0].startswith('sample-nested('), lines[0])
        self.assertIn('posterior=FF', lines[0])

    def test_tasks_are_reported_in_order(self):
        "The tasks of a step are reported in the order in which they are recorded."

        lines = self.dry_run('FF.mode+predict')

        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('find-mode('),            lines[0])
        self.assertTrue(lines[1].startswith('predict-observables('),  lines[1])

    def test_default_arguments_are_applied(self):
        "The step's default arguments are applied to its task."

        lines = self.dry_run('FF.sample')

        for argument in ('nlive=20', 'dlogz=1.0', 'maxiter=100', 'seed=10'):
            self.assertIn(argument, lines[0])

    def test_default_arguments_apply_to_every_task(self):
        "The step's default arguments are applied to each of its tasks."

        lines = self.dry_run('both.sample', analysis_file=ANALYSIS_FILE)

        self.assertEqual(len(lines), 2)
        for line in lines:
            for argument in ('nlive=100', 'dlogz=1.0', 'maxiter=4000'):
                self.assertIn(argument, line)
        self.assertIn('posterior=FF',     lines[0])
        self.assertIn('posterior=CKM+FF', lines[1])

    def test_task_arguments_override_default_arguments(self):
        "A task's own arguments take precedence over the step's default arguments."

        lines = self.dry_run('FF.mode+predict')

        self.assertIn('label=step',      lines[0])
        self.assertNotIn('label=unused', lines[0])

    def test_task_arguments_are_reported(self):
        "The arguments of each task are reported."

        lines = self.dry_run('CKM+FF.sample', analysis_file=ANALYSIS_FILE)

        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith('sample-mcmc('),   lines[0])
        self.assertTrue(lines[1].startswith('find-clusters('), lines[1])
        self.assertTrue(lines[2].startswith('sample-pmc('),    lines[2])
        self.assertIn('chain=0',      lines[0])
        self.assertIn('pre_N=500',    lines[0])
        self.assertIn('step_N=1000',  lines[2])

    def test_analysis_file_and_base_directory_are_added(self):
        "Each task is passed the analysis file and the base directory."

        lines = self.dry_run('FF.mode+predict')

        for line in lines:
            self.assertIn(f'analysis_file={STEPS_FILE}',   line)
            self.assertIn(f'base_directory={self._base}',  line)

    def test_dependencies_are_not_reported(self):
        "The tasks of a step the given step depends on are not reported."

        lines = self.dry_run('FF.mode+predict')

        self.assertNotIn('sample-nested', '\n'.join(lines))

    def test_no_data_is_written(self):
        "A dry run does not write any data."

        self.dry_run('FF.sample')

        self.assertFalse(os.path.exists(os.path.join(self._base, 'data')))

    def test_unknown_step(self):
        "An unknown step id is reported as an error."

        status, stdout, _ = self.run_step('NOSUCH')

        self.assertEqual(status, 1)
        self.assertIn('Step with id \'NOSUCH\' not found in analysis file', stdout)


class RunTests(unittest.TestCase):

    def test_single_task_step(self):
        "Running a step carries out its task."

        self.assertTrue(os.path.isfile(data_path('nested', 'description.yaml')))
        self.assertTrue(os.path.isfile(data_path('nested', 'dynesty_results.npy')))
        self.assertTrue(os.path.isfile(data_path('samples', 'samples.npy')))

    def test_multi_task_step(self):
        "Running a step carries out each of its tasks."

        self.assertTrue(os.path.isfile(data_path('mode-step',      'description.yaml')))
        self.assertTrue(os.path.isfile(data_path('pred-dBRdq2',    'samples.npy')))
        self.assertTrue(os.path.isfile(data_path('pred-dBRdq2',    'weights.npy')))

    def test_task_arguments_override_default_arguments(self):
        "The task's label, not the step's default label, determines the output directory."

        self.assertTrue(os.path.isdir(data_path('mode-step')))
        self.assertFalse(os.path.exists(data_path('mode-unused')))

    def test_each_task_writes_a_log_file(self):
        "Each task records its log output below its own output directory."

        for output in ('nested', 'mode-step', 'pred-dBRdq2'):
            self.assertTrue(os.path.isfile(data_path(output, 'log')), output)

    def test_dependencies_are_not_run(self):
        "A step's dependencies are not run; only its own tasks are."

        with tempfile.TemporaryDirectory() as base:
            status, stdout, _ = run('run', 'FF.mode+predict', base=base)

            # 'find-mode' runs and succeeds, 'predict-observables' then fails for want of the
            # samples that the step 'FF.sample' would have produced.
            self.assertEqual(status, 1)
            self.assertTrue(os.path.isdir(data_path('mode-step', base=base)))
            self.assertIn(data_path('samples', base=base), stdout)

    def test_unknown_step(self):
        "An unknown step id is reported as an error."

        with tempfile.TemporaryDirectory() as base:
            status, stdout, _ = run('run', 'NOSUCH', base=base)

            self.assertEqual(status, 1)
            self.assertIn('Step with id \'NOSUCH\' not found in analysis file', stdout)


class ParserTests(unittest.TestCase):

    def test_defaults(self):
        "The command's arguments carry the documented defaults."

        args = analysis._parser().parse_args(['run', 'FF.sample'])

        self.assertEqual(args.id, 'FF.sample')
        self.assertFalse(args.dry_run)
        self.assertEqual(args.analysis_file, 'analysis.yaml')

    def test_arguments_are_forwarded(self):
        "Every option reaches the task under its documented name."

        with mock.patch.object(eos, 'run') as task:
            status, _, _ = run('run', '-d', 'FF.sample', base='/base')

        self.assertEqual(status, 0)
        task.assert_called_once_with(
            analysis_file=STEPS_FILE, id='FF.sample', base_directory='/base', dry_run=True,
        )


if __name__ == '__main__':
    unittest.main(verbosity=5)
