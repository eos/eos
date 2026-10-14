import os
import shutil
import sys
import tempfile
import types
import unittest
import unittest.mock
import eos


class _SuppressingOutput:
    "Stands in for ipywidgets.Output, which displays and then suppresses what is raised within it."

    def __init__(self, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, etype, evalue, tb):
        return True


class _Accordion:
    "Stands in for ipywidgets.Accordion."

    def __init__(self, children=None):
        self.selected_index = 0

    def set_title(self, index, title):
        pass


class TaskFailureTests(unittest.TestCase):

    def test_failure_escapes_the_output_widget(self):
        "A task that raises under IPython fails its caller, although the output widget suppresses."
        widgets         = types.ModuleType('ipywidgets')
        widgets.Output  = _SuppressingOutput
        widgets.Accordion = _Accordion

        ipython         = types.ModuleType('IPython')
        ipython.display = types.ModuleType('IPython.display')
        ipython.display.display = lambda *args, **kwargs: None

        modules = { 'ipywidgets': widgets, 'IPython': ipython, 'IPython.display': ipython.display }

        @eos.tasks.task('test-failure', '', logfile=False)
        def failing_task():
            raise RuntimeError('the task failed')

        self.addCleanup(eos.tasks._tasks.pop, 'test-failure', None)
        self.addCleanup(eos.tasks._task_outputs.pop, 'test-failure', None)
        self.addCleanup(eos.tasks._task_output_fields.pop, 'test-failure', None)

        with unittest.mock.patch.dict(sys.modules, modules), \
             unittest.mock.patch.object(eos.tasks._framework, '__ipython__', True):
            with self.assertRaises(RuntimeError):
                failing_task()



class TaskOutputTests(unittest.TestCase):
    "Tests of the output-directory handling of the task decorator (issue #1291)."

    def setUp(self):
        self.base = tempfile.mkdtemp(prefix='eos-task-output-')
        self.addCleanup(shutil.rmtree, self.base, ignore_errors=True)

        @eos.tasks.task('test-output', 'data/{posterior}/out-{label}', load_analysis_file=False)
        def output_task(posterior, label, base_directory, content, fail=False):
            path = eos.tasks._staged(os.path.join(base_directory, 'data', posterior, f'out-{label}'))
            with open(os.path.join(path, 'result'), 'w') as f:
                f.write(content)
            if fail:
                raise RuntimeError('the task failed')

        self.addCleanup(eos.tasks._tasks.pop, 'test-output', None)
        self.addCleanup(eos.tasks._task_outputs.pop, 'test-output', None)
        self.addCleanup(eos.tasks._task_output_fields.pop, 'test-output', None)
        self.task = output_task

    def _register(self, name, *args, **kwargs):
        self.addCleanup(eos.tasks._tasks.pop, name, None)
        self.addCleanup(eos.tasks._task_outputs.pop, name, None)
        self.addCleanup(eos.tasks._task_output_fields.pop, name, None)
        return eos.tasks.task(name, *args, load_analysis_file=False, **kwargs)

    def _read(self, label):
        with open(os.path.join(self.base, 'data', 'P', f'out-{label}', 'result')) as f:
            return f.read()

    def test_invalid_names_are_rejected(self):
        for label in ('foo/bar', 'foo\\bar', '../x', 'a b', '.', '..', ''):
            with self.subTest(label=label), self.assertRaises(ValueError):
                self.task(posterior='P', label=label, base_directory=self.base, content='x')
        for posterior in ('.', '..'):
            with self.subTest(posterior=posterior), self.assertRaises(ValueError):
                self.task(posterior=posterior, label='l', base_directory=self.base, content='x')
        self.assertEqual([], os.listdir(self.base))

    def test_success_creates_directory(self):
        self.task(posterior='P', label='l', base_directory=self.base, content='v1')
        self.assertEqual('v1', self._read('l'))
        self.assertEqual(['out-l'], os.listdir(os.path.join(self.base, 'data', 'P')))
        self.assertTrue(os.path.isfile(os.path.join(self.base, 'data', 'P', 'out-l', 'log')))

    def test_failure_leaves_no_directory(self):
        with self.assertRaises(RuntimeError):
            self.task(posterior='P', label='l', base_directory=self.base, content='v1', fail=True)
        self.assertEqual(['.out-l.failed.log'], os.listdir(os.path.join(self.base, 'data', 'P')))

    def test_failed_log_is_kept_until_success(self):
        failedlog = os.path.join(self.base, 'data', 'P', '.out-l.failed.log')
        with self.assertRaises(RuntimeError):
            self.task(posterior='P', label='l', base_directory=self.base, content='v1', fail=True)
        self.assertTrue(os.path.isfile(failedlog))
        self.task(posterior='P', label='l', base_directory=self.base, content='v1')
        self.assertFalse(os.path.exists(failedlog))

    def test_failure_keeps_previous_result(self):
        self.task(posterior='P', label='l', base_directory=self.base, content='v1')
        with self.assertRaises(RuntimeError):
            self.task(posterior='P', label='l', base_directory=self.base, content='v2', fail=True)
        self.assertEqual('v1', self._read('l'))
        self.assertEqual(['.out-l.failed.log', 'out-l'], sorted(os.listdir(os.path.join(self.base, 'data', 'P'))))

    def test_success_replaces_previous_result(self):
        self.task(posterior='P', label='l', base_directory=self.base, content='v1')
        self.task(posterior='P', label='l', base_directory=self.base, content='v2')
        self.assertEqual('v2', self._read('l'))
        self.assertEqual(['out-l'], os.listdir(os.path.join(self.base, 'data', 'P')))

    def test_nested_task_commits_with_outermost_task(self):
        @self._register('test-outer', 'data/{posterior}/outer')
        def outer_task(posterior, base_directory, fail=False):
            self.task(posterior=posterior, label='inner', base_directory=base_directory, content='inner')
            if fail:
                raise RuntimeError('the task failed')

        with self.assertRaises(RuntimeError):
            outer_task(posterior='P', base_directory=self.base, fail=True)
        self.assertEqual(['.outer.failed.log'], os.listdir(os.path.join(self.base, 'data', 'P')))

        outer_task(posterior='P', base_directory=self.base)
        self.assertEqual('inner', self._read('inner'))
        self.assertEqual(['out-inner', 'outer'], sorted(os.listdir(os.path.join(self.base, 'data', 'P'))))

    def test_derived_output_fields(self):
        @self._register('test-masked', 'data/{posterior}/out-{label}{mask}', output_fields=eos.tasks._mask_output_fields)
        def masked_task(posterior, label, base_directory, content, mask_name=None):
            path = eos.tasks._staged(os.path.join(base_directory, 'data', posterior, f'out-{label}{eos.tasks._mask_output_fields(mask_name)["mask"]}'))
            with open(os.path.join(path, 'result'), 'w') as f:
                f.write(content)

        masked_task(posterior='P', label='l', base_directory=self.base, content='unmasked')
        masked_task(posterior='P', label='l', base_directory=self.base, content='masked', mask_name='m')
        self.assertEqual('unmasked', self._read('l'))
        self.assertEqual('masked', self._read('l_mask-m'))
        self.assertEqual(['out-l', 'out-l_mask-m'], sorted(os.listdir(os.path.join(self.base, 'data', 'P'))))
        with self.assertRaises(ValueError):
            masked_task(posterior='P', label='l', base_directory=self.base, content='x', mask_name='a/b')

    def test_staged_files_keep_previous_files_on_failure(self):
        @self._register('test-files', 'files', staged=False)
        def files_task(base_directory, content, fail=False):
            for name in ('a.txt', 'b.txt'):
                with open(eos.tasks._staged_file(os.path.join(base_directory, 'shared', name)), 'w') as f:
                    f.write(content)
            if fail:
                raise RuntimeError('the task failed')

        shared = os.path.join(self.base, 'shared')
        files_task(base_directory=self.base, content='v1')
        with open(os.path.join(shared, 'other.txt'), 'w') as f:
            f.write('other')
        with self.assertRaises(RuntimeError):
            files_task(base_directory=self.base, content='v2', fail=True)
        self.assertEqual(['a.txt', 'b.txt', 'other.txt'], sorted(os.listdir(shared)))
        for name in ('a.txt', 'b.txt'):
            with open(os.path.join(shared, name)) as f:
                self.assertEqual('v1', f.read())

        files_task(base_directory=self.base, content='v2')
        self.assertEqual(['a.txt', 'b.txt', 'other.txt'], sorted(os.listdir(shared)))
        with open(os.path.join(shared, 'a.txt')) as f:
            self.assertEqual('v2', f.read())

    def test_staged_file_below_staged_directory(self):
        @self._register('test-below', 'data/{posterior}/out-{label}')
        def below_task(posterior, label, base_directory, fail=False):
            with open(eos.tasks._staged_file(os.path.join(base_directory, 'data', posterior, f'out-{label}', 'extra.txt')), 'w') as f:
                f.write('extra')
            if fail:
                raise RuntimeError('the task failed')

        with self.assertRaises(RuntimeError):
            below_task(posterior='P', label='l', base_directory=self.base, fail=True)
        self.assertEqual(['.out-l.failed.log'], os.listdir(os.path.join(self.base, 'data', 'P')))

        below_task(posterior='P', label='l', base_directory=self.base)
        with open(os.path.join(self.base, 'data', 'P', 'out-l', 'extra.txt')) as f:
            self.assertEqual('extra', f.read())
        self.assertEqual(['out-l'], os.listdir(os.path.join(self.base, 'data', 'P')))

    def test_failed_commit_restores_previous_outputs(self):
        @self._register('test-two', 'data/{posterior}/out-a')
        def two_task(posterior, base_directory, content):
            for label in ('a', 'b'):
                path = eos.tasks._staged(os.path.join(base_directory, 'data', posterior, f'out-{label}'))
                with open(os.path.join(path, 'result'), 'w') as f:
                    f.write(content)

        two_task(posterior='P', base_directory=self.base, content='v1')

        rename = os.rename
        def failing_rename(source, target):
            if os.path.basename(target) == 'out-b' and '.tmp' in os.path.basename(source):
                raise OSError('the commit failed')
            rename(source, target)

        with unittest.mock.patch('os.rename', failing_rename), self.assertRaises(OSError):
            two_task(posterior='P', base_directory=self.base, content='v2')
        self.assertEqual('v1', self._read('a'))
        self.assertEqual('v1', self._read('b'))
        self.assertEqual(['out-a', 'out-b'], sorted(os.listdir(os.path.join(self.base, 'data', 'P'))))

    def test_backed_up_files_are_restored_on_failure(self):
        @self._register('test-backup', 'files', staged=False)
        def backup_task(base_directory, content, fail=False):
            for name in ('a.txt', 'b.txt'):
                with open(eos.tasks._backed_up(os.path.join(base_directory, 'shared', name)), 'w') as f:
                    f.write(content)
            if fail:
                raise RuntimeError('the task failed')

        shared = os.path.join(self.base, 'shared')
        os.makedirs(shared)
        with open(os.path.join(shared, 'a.txt'), 'w') as f:
            f.write('v1')
        with self.assertRaises(RuntimeError):
            backup_task(base_directory=self.base, content='v2', fail=True)
        self.assertEqual(['a.txt'], os.listdir(shared))
        with open(os.path.join(shared, 'a.txt')) as f:
            self.assertEqual('v1', f.read())

        backup_task(base_directory=self.base, content='v2')
        self.assertEqual(['a.txt', 'b.txt'], sorted(os.listdir(shared)))
        with open(os.path.join(shared, 'a.txt')) as f:
            self.assertEqual('v2', f.read())

    def test_staging_parent_after_child_is_rejected(self):
        @self._register('test-parent', 'data/{posterior}/out')
        def parent_task(posterior, base_directory):
            eos.tasks._staged(os.path.join(base_directory, 'data', posterior, 'out', 'child'))
            eos.tasks._staged(os.path.join(base_directory, 'data', posterior))

        with self.assertRaises(RuntimeError):
            parent_task(posterior='P', base_directory=self.base)
        self.assertEqual(['.out.failed.log'], os.listdir(os.path.join(self.base, 'data', 'P')))


if __name__ == '__main__':
    unittest.main(verbosity=5)
