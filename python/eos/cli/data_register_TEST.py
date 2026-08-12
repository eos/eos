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

from contextlib import contextmanager
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

from eos.cli import data
from eos.cli.data_github import (
    GitHubClient, GitHubCommandError, GitHubRef, GitHubRelease, GitHubTag,
)
from eos.cli.data_register import (
    RegistrationError, execute_register, parse_doi, parse_eos_version,
    parse_keywords, parse_likelihoods, parse_repository, plan_register,
    render_register_plan,
)
from eos.cli.data_release import LocalGitClient, SourceCheckError, make_annotation


def _dataset(root: Path, *, authors=('Jane Doe', 'John Roe')):
    (root / 'README.md').write_text('# Dataset\n', encoding='utf-8')
    author_yaml = ''.join(f'    - name: {author}\n' for author in authors)
    (root / 'analysis.yaml').write_text(
        'metadata:\n  id: 2025-03\n  title: Main title\n'
        '  description: Release notes.\n  authors:\n' + author_yaml,
        encoding='utf-8',
    )
    (root / '.zenodo.json').write_text(json.dumps({
        'title': 'EOS/DATA-2026-01: Supplementary material for EOS/ANALYSIS-2025-03',
        'description': 'Release notes.',
        'creators': [
            {'name': author, 'affiliation': 'Lab'} for author in authors
        ],
        'upload_type': 'dataset', 'license': 'CC-BY-4.0', 'grants': [],
    }), encoding='utf-8')


REGISTRY = """2025-01:
  authors:
  - Existing Author
  title: Existing title
  likelihoods: {}
  keywords:
  - existing
  eos_version: 1.0.0
  future-field: retained
"""


class FakeGitHub:
    def __init__(self, tagged: Path, main: Path):
        annotation = make_annotation('2026-01', '2025-03', ('analysis.yaml',), 'analysis.yaml')
        self.tag = GitHubTag('2026-01', 'tag-object', 'tag-commit', annotation.serialize())
        self.roots = {'tag-commit': tagged, 'main-commit': main}
        self.branches = {'main': 'main-commit'}
        self.release_url = 'https://example.invalid/release'
        self.before_push = None
        self.after_push = None
        self.calls = []

    def list_tag_refs(self):
        self.calls.append(('list-tag-refs',))
        return (GitHubRef('2026-01', 'tag-object', 'tag'),)

    def get_tag(self, ref):
        self.calls.append(('get-tag', ref.name))
        return self.tag

    def get_release(self, name):
        self.calls.append(('get-release', name))
        return GitHubRelease(name, self.release_url)

    def get_branch(self, name):
        self.calls.append(('get-branch', name))
        return self.branches.get(name)

    @contextmanager
    def checkout_tree(self, commit):
        self.calls.append(('checkout', commit))
        yield self.roots[commit]

    def transfer_local_branch(self, root, branch, commit, main, tag_name, tag_object):
        self.calls.append(('push-branch', branch, commit, main, tag_name, tag_object))
        if self.before_push is not None:
            self.before_push()
        if (
            branch in self.branches or self.branches.get('main') != main
            or self.tag.name != tag_name or self.tag.object_id != tag_object
        ):
            raise RuntimeError('atomic ref lease rejected')
        self.branches[branch] = commit
        if self.after_push is not None:
            self.after_push()


class FakeGit:
    def __init__(self):
        self.calls = []

    def create_registration_commit(self, root, branch, parent, text, message):
        self.calls.append(('create-commit', branch, parent, message, text))
        return 'registration-commit'


class RegisterTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        base = Path(temporary.name)
        self.tagged, self.main = base / 'tagged', base / 'main'
        self.tagged.mkdir()
        self.main.mkdir()
        _dataset(self.tagged)
        (self.main / 'datasets.yaml').write_text(REGISTRY, encoding='utf-8')
        self.github = FakeGitHub(self.tagged, self.main)

    def plan(self, **kwargs):
        values = {
            'doi': '10.5281/zenodo.12345', 'keywords': ('flavour', 'fits'),
            'eos_version': '1.0.20', 'github': self.github,
            'check_factory': data.create_check_factory(),
        }
        values.update(kwargs)
        return plan_register('2026-01', **values)

    def test_input_parsing(self):
        self.assertEqual(parse_doi('10.5281/zenodo.123'), '10.5281/zenodo.123')
        self.assertEqual(parse_keywords(('one', 'two')), ('one', 'two'))
        self.assertEqual(parse_eos_version('1.0.20'), '1.0.20')
        self.assertEqual(parse_repository('other/data'), 'other/data')
        self.assertEqual(
            parse_likelihoods(('fit:path/to/file:MixtureDensity',)),
            (('fit', 'path/to/file', 'MixtureDensity'),),
        )
        for function, value in (
            (parse_doi, 'https://doi.org/10.5281/zenodo.1'),
            (parse_keywords, ('one', 'one')),
            (parse_eos_version, ' '),
            (parse_repository, 'data'),
            (parse_likelihoods, ('fit:file',)),
            (parse_likelihoods, ('fit:a:MixtureDensity', 'fit:b:NabuLikelihood')),
        ):
            with self.subTest(value=value), self.assertRaises(RegistrationError):
                function(value)
        self.assertEqual(parse_likelihoods(()), ())

    def test_cli_parser_requires_registry_fields_and_accepts_repeats(self):
        arguments = [
            'register', '2026-01', '--doi', '10.5281/zenodo.123',
            '--keyword', 'one', '--keyword', 'two', '--eos-version', '1.2.3',
            '--likelihood', 'fit:file:MixtureDensity', '--repo', 'other/data', '--dry-run',
        ]
        parsed = data._parser().parse_args(arguments)
        self.assertEqual(parsed.keyword, ['one', 'two'])
        self.assertEqual(parsed.likelihood, ['fit:file:MixtureDensity'])
        self.assertEqual(parsed.repo, 'other/data')
        self.assertTrue(parsed.dry_run)
        for missing in ('--doi', '--keyword', '--eos-version'):
            reduced = [
                'register', '2026-01', '--doi', '10.5281/zenodo.123',
                '--keyword', 'one', '--eos-version', '1.2.3',
            ]
            index = reduced.index(missing)
            del reduced[index:index + 2]
            self.assertEqual(data.main(reduced, stderr=io.StringIO()), 2)

    def test_plan_derives_metadata_and_preserves_registry_bytes(self):
        plan = self.plan(likelihoods=('fit:result.yaml:MixtureDensity',))
        self.assertEqual(plan.entry.authors, ('Jane Doe', 'John Roe'))
        self.assertEqual(plan.entry.title, 'Main title')
        self.assertEqual(plan.branch_name, 'register-2026-01')
        self.assertEqual(plan.commit_message, 'Add 2026-01')
        self.assertTrue(plan.updated_registry.endswith(REGISTRY))
        self.assertIn('future-field: retained', plan.updated_registry)
        rendered = render_register_plan(plan)
        for expected in ('Source main commit: main-commit', 'Prospective registry entry:',
                         'YAML diff:', 'push-new-branch'):
            self.assertIn(expected, rendered)

    def test_common_author_set_uses_main_file_order(self):
        (self.tagged / 'secondary.yaml').write_text(
            'metadata:\n  id: 2025-03\n  title: Secondary\n'
            '  description: Release notes.\n  authors:\n'
            '    - name: John Roe\n    - name: Jane Doe\n',
            encoding='utf-8',
        )
        annotation = make_annotation(
            '2026-01', '2025-03', ('secondary.yaml', 'analysis.yaml'), 'analysis.yaml',
        )
        self.github.tag = GitHubTag('2026-01', 'tag-object', 'tag-commit', annotation.serialize())
        self.assertEqual(self.plan().entry.authors, ('Jane Doe', 'John Roe'))
        (self.tagged / 'secondary.yaml').write_text(
            'metadata:\n  id: 2025-03\n  title: Secondary\n'
            '  description: Release notes.\n  authors:\n    - name: Jane Doe\n',
            encoding='utf-8',
        )
        with self.assertRaises(SourceCheckError):
            self.plan()

    def test_rejects_existing_id_malformed_registry_and_branch_conflict(self):
        valid = self.plan()
        for contents, message in (
            (valid.updated_registry, 'already contains'),
            ('- not-a-mapping\n', 'string-keyed mapping'),
        ):
            self.main.joinpath('datasets.yaml').write_text(contents, encoding='utf-8')
            with self.assertRaisesRegex(RegistrationError, message):
                self.plan()
        self.main.joinpath('datasets.yaml').write_text(REGISTRY, encoding='utf-8')
        self.github.branches['register-2026-01'] = 'other'
        with self.assertRaisesRegex(RegistrationError, 'already exists'):
            self.plan()

    def test_release_and_tag_are_mandatory(self):
        with mock.patch.object(self.github, 'get_release', return_value=None):
            with self.assertRaisesRegex(RegistrationError, 'release'):
                self.plan()
        with mock.patch.object(
            self.github, 'list_tag_refs',
            return_value=(GitHubRef('2026-01', 'commit', 'commit'),),
        ):
            with self.assertRaisesRegex(RegistrationError, 'lightweight'):
                self.plan()

    def test_execute_replans_commits_once_pushes_and_verifies(self):
        plan = self.plan()
        git = FakeGit()
        result = execute_register(
            plan, git=git, github=self.github, check_factory=data.create_check_factory(),
        )
        self.assertEqual(result.commit_id, 'registration-commit')
        self.assertEqual(len(git.calls), 1)
        writes = [call[0] for call in self.github.calls if call[0] == 'push-branch']
        self.assertEqual(writes, ['push-branch'])
        self.assertEqual(self.github.branches['register-2026-01'], 'registration-commit')

    def test_execute_rejects_stale_main_and_branch_race_without_push(self):
        plan = self.plan()
        git = FakeGit()
        self.github.roots['advanced-main'] = self.main
        self.github.branches['main'] = 'advanced-main'
        with self.assertRaisesRegex(RegistrationError, 'preconditions changed'):
            execute_register(
                plan, git=git, github=self.github, check_factory=data.create_check_factory(),
            )
        self.assertEqual(git.calls, [])
        self.github.branches['main'] = 'main-commit'
        original = git.create_registration_commit
        def raced(*arguments):
            commit = original(*arguments)
            self.github.branches['register-2026-01'] = 'racing-commit'
            return commit
        git.create_registration_commit = raced
        with self.assertRaisesRegex(RegistrationError, 'appeared before transfer'):
            execute_register(
                plan, git=git, github=self.github, check_factory=data.create_check_factory(),
            )
        self.assertNotIn('push-branch', [call[0] for call in self.github.calls])

    def test_atomic_push_rejects_main_or_tag_race_after_final_validation(self):
        for kind in ('main', 'tag'):
            with self.subTest(kind=kind):
                self.setUp()
                plan = self.plan()
                if kind == 'main':
                    self.github.before_push = lambda: self.github.branches.__setitem__(
                        'main', 'racing-main',
                    )
                else:
                    def move_tag():
                        self.github.tag = GitHubTag(
                            '2026-01', 'racing-tag', 'tag-commit', self.github.tag.message,
                        )
                    self.github.before_push = move_tag
                with self.assertRaisesRegex(RegistrationError, 'atomic ref lease rejected'):
                    execute_register(
                        plan, git=FakeGit(), github=self.github,
                        check_factory=data.create_check_factory(),
                    )
                self.assertNotIn('register-2026-01', self.github.branches)

    def test_post_push_release_race_is_reported_as_partial_failure(self):
        plan = self.plan()
        self.github.after_push = lambda: setattr(
            self.github, 'release_url', 'https://example.invalid/changed-release',
        )
        with self.assertRaisesRegex(
            RegistrationError, 'post-push publication verification failed',
        ) as caught:
            execute_register(
                plan, git=FakeGit(), github=self.github,
                check_factory=data.create_check_factory(),
            )
        self.assertEqual(
            self.github.branches.get('register-2026-01'), 'registration-commit',
        )
        self.assertIn('Completed operations:', str(caught.exception))

    def test_local_git_creates_one_commit_changing_only_registry(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        environment = {
            **os.environ, 'GIT_AUTHOR_NAME': 'Test', 'GIT_AUTHOR_EMAIL': 'test@example.com',
            'GIT_COMMITTER_NAME': 'Test', 'GIT_COMMITTER_EMAIL': 'test@example.com',
        }
        def git(*arguments):
            return subprocess.run(
                ['git', *arguments], cwd=root, env=environment, check=True,
                capture_output=True, text=True,
            ).stdout.strip()
        git('init', '-q')
        git('config', 'user.name', 'Test')
        git('config', 'user.email', 'test@example.com')
        (root / 'datasets.yaml').write_text(REGISTRY, encoding='utf-8')
        (root / 'untouched').write_text('same\n', encoding='utf-8')
        git('add', 'datasets.yaml', 'untouched')
        git('commit', '-q', '-m', 'main')
        parent = git('rev-parse', 'HEAD')
        commit = LocalGitClient().create_registration_commit(
            root, 'register-2026-01', parent, 'new registry\n', 'Add 2026-01',
        )
        self.assertEqual(git('show', '-s', '--format=%P', commit), parent)
        self.assertEqual(git('show', '-s', '--format=%s', commit), 'Add 2026-01')
        self.assertEqual(git('diff-tree', '--no-commit-id', '--name-only', '-r', commit), 'datasets.yaml')

    def test_cli_dry_run_is_read_only_and_real_run_prints_pr_reminder(self):
        arguments = [
            'register', '2026-01', '--doi', '10.5281/zenodo.12345',
            '--keyword', 'flavour', '--eos-version', '1.0.20',
        ]
        git = FakeGit()
        output = io.StringIO()
        files_before = self.main.joinpath('datasets.yaml').read_bytes()
        refs_before = dict(self.github.branches)
        self.assertEqual(data.main(
            [*arguments, '--dry-run'], stdout=output, stderr=io.StringIO(),
            git_client=git, github_client=self.github,
        ), 0)
        self.assertEqual(git.calls, [])
        self.assertNotIn('push-branch', [call[0] for call in self.github.calls])
        self.assertEqual(self.main.joinpath('datasets.yaml').read_bytes(), files_before)
        self.assertEqual(self.github.branches, refs_before)
        output = io.StringIO()
        self.assertEqual(data.main(
            arguments, stdout=output, stderr=io.StringIO(),
            git_client=git, github_client=self.github,
        ), 0)
        self.assertIn('Create a pull request', output.getvalue())

    def test_no_pr_capability_exists(self):
        self.assertFalse(hasattr(self.github, 'create_pull_request'))
        source = Path(data.__file__).with_name('data_register.py').read_text(encoding='utf-8')
        self.assertNotIn('gh pr', source)
        self.assertNotIn('create_pull_request', source)

    def test_github_branch_transfer_uses_guarded_argument_array(self):
        result = subprocess.CompletedProcess([], 0, stdout='', stderr='')
        runner = mock.Mock(return_value=result)
        GitHubClient(executable='unused', git_runner=runner).transfer_local_branch(
            Path('/checkout'), 'register-2026-01', 'commit-id', 'main-id',
            '2026-01', 'tag-id',
        )
        command = runner.call_args.args[0]
        self.assertEqual(command[:4], ['git', 'push', '--porcelain', '--atomic'])
        self.assertIn('--force-with-lease=refs/heads/main:main-id', command)
        self.assertIn('--force-with-lease=refs/tags/2026-01:tag-id', command)
        self.assertIn('--force-with-lease=refs/heads/register-2026-01:', command)
        self.assertEqual(command[-4:], [
            'origin', 'main-id:refs/heads/main', 'tag-id:refs/tags/2026-01',
            'commit-id:refs/heads/register-2026-01',
        ])

    def test_real_atomic_push_does_not_create_branch_after_main_or_tag_race(self):
        environment = {
            **os.environ, 'GIT_AUTHOR_NAME': 'Test', 'GIT_AUTHOR_EMAIL': 'test@example.com',
            'GIT_COMMITTER_NAME': 'Test', 'GIT_COMMITTER_EMAIL': 'test@example.com',
        }
        for raced_ref in ('main', 'tag'):
            with self.subTest(raced_ref=raced_ref):
                temporary = tempfile.TemporaryDirectory()
                self.addCleanup(temporary.cleanup)
                base = Path(temporary.name)
                remote, seed, client, racer = (
                    base / 'remote.git', base / 'seed', base / 'client', base / 'racer',
                )
                def run(root, *arguments):
                    return subprocess.run(
                        ['git', *arguments], cwd=root, env=environment, check=True,
                        capture_output=True, text=True,
                    ).stdout.strip()
                subprocess.run(
                    ['git', 'init', '--bare', '-q', str(remote)], check=True,
                    capture_output=True, text=True,
                )
                seed.mkdir()
                run(seed, 'init', '-q')
                run(seed, 'config', 'user.name', 'Test')
                run(seed, 'config', 'user.email', 'test@example.com')
                seed.joinpath('datasets.yaml').write_text(REGISTRY, encoding='utf-8')
                run(seed, 'add', 'datasets.yaml')
                run(seed, 'commit', '-q', '-m', 'main')
                run(seed, 'branch', '-M', 'main')
                run(seed, 'tag', '-a', '2026-01', '-m', 'dataset')
                run(seed, 'remote', 'add', 'origin', str(remote))
                run(seed, 'push', '-q', 'origin', 'main', 'refs/tags/2026-01')
                subprocess.run(
                    ['git', 'clone', '-q', '-b', 'main', str(remote), str(client)],
                    check=True, capture_output=True, text=True,
                )
                subprocess.run(
                    ['git', 'clone', '-q', '-b', 'main', str(remote), str(racer)],
                    check=True, capture_output=True, text=True,
                )
                for checkout in (client, racer):
                    run(checkout, 'config', 'user.name', 'Test')
                    run(checkout, 'config', 'user.email', 'test@example.com')
                expected_main = run(client, 'rev-parse', 'origin/main^{commit}')
                expected_tag = run(client, 'rev-parse', 'refs/tags/2026-01^{tag}')
                run(client, 'checkout', '-q', '-b', 'register-2026-01')
                client.joinpath('datasets.yaml').write_text('registration\n', encoding='utf-8')
                run(client, 'add', 'datasets.yaml')
                run(client, 'commit', '-q', '-m', 'Add 2026-01')
                registration_commit = run(client, 'rev-parse', 'HEAD^{commit}')
                if raced_ref == 'main':
                    racer.joinpath('race').write_text('main moved\n', encoding='utf-8')
                    run(racer, 'add', 'race')
                    run(racer, 'commit', '-q', '-m', 'race main')
                    run(racer, 'push', '-q', 'origin', 'HEAD:refs/heads/main')
                else:
                    run(racer, 'tag', '-f', '-a', '2026-01', '-m', 'moved dataset')
                    run(racer, 'push', '-q', '--force', 'origin', 'refs/tags/2026-01')
                with self.assertRaises(GitHubCommandError):
                    GitHubClient(executable='unused').transfer_local_branch(
                        client, 'register-2026-01', registration_commit,
                        expected_main, '2026-01', expected_tag,
                    )
                refs = subprocess.run(
                    ['git', 'ls-remote', '--heads', str(remote), 'register-2026-01'],
                    check=True, capture_output=True, text=True,
                ).stdout
                self.assertEqual(refs, '')

    def test_missing_gh_is_focused_register_failure(self):
        arguments = [
            'register', '2026-01', '--doi', '10.5281/zenodo.12345',
            '--keyword', 'flavour', '--eos-version', '1.0.20', '--dry-run',
        ]
        error = io.StringIO()
        with mock.patch('shutil.which', return_value=None):
            self.assertEqual(data.main(arguments, stderr=error), 2)
        self.assertIn("GitHub CLI 'gh' is required", error.getvalue())


if __name__ == '__main__':
    unittest.main(verbosity=5)
