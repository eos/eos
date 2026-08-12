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
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from eos.cli.data import create_check_factory
from eos.cli.data_github import GitHubClient, GitHubCommandError, GitHubRef, GitHubRelease, GitHubTag
from eos.cli.data_release import (
    LocalGitClient, ReleasePlanningError, ResolvedSource, TargetCheckout,
    make_annotation, plan_create, plan_publish, render_create_plan, render_publish_plan,
)
from eos.cli.data_release_execution import (
    PartialPublishError, ReleaseExecutionError, execute_create, execute_publish,
)


def _dataset(root: Path, data_id: str) -> None:
    (root / 'README.md').write_text('# Dataset\n', encoding='utf-8')
    (root / 'analysis.yaml').write_text(
        'metadata:\n  id: 2025-03\n  title: Main title\n'
        '  description: Release notes.\n  authors:\n    - name: Jane Doe\n',
        encoding='utf-8',
    )
    (root / '.zenodo.json').write_text(json.dumps({
        'title': f'EOS/DATA-{data_id}: Supplementary material for EOS/ANALYSIS-2025-03',
        'description': 'Release notes.',
        'creators': [{'name': 'Doe, Jane', 'affiliation': 'Lab'}],
        'upload_type': 'dataset', 'license': 'CC-BY-4.0', 'grants': [],
    }), encoding='utf-8')


def _tag(name: str, commit: str) -> GitHubTag:
    annotation = make_annotation(name, '2025-03', ('analysis.yaml',), 'analysis.yaml')
    return GitHubTag(name, f'tag-{name}', commit, annotation.serialize())


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ['git', *arguments], cwd=root, check=True, capture_output=True, text=True,
        env={
            **os.environ,
            'GIT_AUTHOR_NAME': 'Test', 'GIT_AUTHOR_EMAIL': 'test@example.com',
            'GIT_COMMITTER_NAME': 'Test', 'GIT_COMMITTER_EMAIL': 'test@example.com',
        },
    ).stdout.strip()


class ScriptedGit:
    def __init__(self, root: Path, tags=(), branch=None):
        self.root = root
        self.tags = {tag.name: (tag.object_id, tag.commit_id) for tag in tags}
        self.messages = {tag.name: tag.message for tag in tags}
        self.branch = branch
        self.parents = {tag.commit_id: (() if index == 0 else (tags[index - 1].commit_id,))
                        for index, tag in enumerate(tags)}
        self.calls = []
        self.next_commit = 'new-commit'

    def target_checkout(self, directory, *, require_clean=True):
        self.calls.append(('target', require_clean))
        return TargetCheckout(Path(directory), 'origin', 'git@github.com:eos/data.git')

    def local_tag_refs(self, target):
        self.calls.append(('local-tags',))
        return dict(self.tags)

    def local_branch(self, target, name):
        self.calls.append(('local-branch', name))
        return self.branch

    @contextmanager
    def checkout_source(self, source):
        self.calls.append(('check-source', str(source)))
        yield ResolvedSource(str(source), 'source-commit', self.root)

    @contextmanager
    def checkout_commit(self, source, commit):
        self.calls.append(('checkout-commit', str(source), commit))
        yield self.root

    def create_dataset_commit(self, target, source, message, parent):
        self.calls.append(('create-commit', message, parent))
        self.parents[self.next_commit] = () if parent is None else (parent,)
        return self.next_commit

    def prepare_local_branch(self, target, name, current, parent):
        self.calls.append(('prepare-branch', name, current, parent))

    def commit_parents(self, target, commit):
        self.calls.append(('local-parents', commit))
        return self.parents[commit]

    def update_local_branch(self, target, name, commit, expected_old):
        self.calls.append(('move-branch', name, commit, expected_old))
        if self.branch != expected_old:
            raise RuntimeError('local branch race')
        self.branch = commit

    def create_local_tag(self, target, name, commit, message, expected_old):
        self.calls.append(('create-tag', name, commit, expected_old))
        current = self.tags.get(name)
        if (current[0] if current else None) != expected_old:
            raise RuntimeError('local tag race')
        object_id = f'new-tag-{name}'
        self.tags[name] = (object_id, commit)
        self.messages[name] = message
        return self.tags[name]


class ScriptedGitHub:
    def __init__(self, root: Path, tags=(), releases=(), branch=None, parents=None):
        self.root = root
        self.tags = {tag.name: tag for tag in tags}
        self.releases = set(releases)
        self.branch = branch
        self.parents = dict(parents or {})
        self.calls = []
        self.git = None
        self.fail_release = False
        self.fail_branch = False
        self.wrong_transfer = False
        self.move_tag_after_release = False

    def list_tag_refs(self):
        self.calls.append(('list-tags',))
        return tuple(GitHubRef(tag.name, tag.object_id, 'tag') for tag in self.tags.values())

    def get_tag(self, ref):
        self.calls.append(('get-tag', ref.name))
        return self.tags[ref.name]

    def get_release(self, name):
        self.calls.append(('get-release', name))
        return GitHubRelease(name, f'https://example.invalid/{name}') if name in self.releases else None

    def get_branch(self, name):
        self.calls.append(('get-branch', name))
        return self.branch

    def commit_parents(self, commit):
        self.calls.append(('parents', commit))
        return tuple(self.parents.get(commit, ()))

    def is_ancestor(self, ancestor, descendant):
        self.calls.append(('ancestor', ancestor, descendant))
        current = descendant
        while self.parents.get(current):
            current = self.parents[current][0]
            if current == ancestor:
                return True
        return ancestor == descendant

    @contextmanager
    def checkout_tree(self, commit):
        self.calls.append(('checkout-tree', commit))
        yield self.root

    def transfer_local_tag(self, root, remote, name, expected_old):
        self.calls.append(('transfer-tag', remote, name, expected_old))
        current = self.tags.get(name)
        if (current.object_id if current else None) != expected_old:
            raise GitHubCommandError('remote tag race')
        object_id, commit_id = self.git.tags[name]
        if self.wrong_transfer:
            object_id = 'wrong-object'
        self.tags[name] = GitHubTag(name, object_id, commit_id, self.git.messages[name])
        self.parents[commit_id] = self.git.parents[commit_id]

    def create_release(self, name, title, notes):
        self.calls.append(('create-release', name, title, notes))
        if self.fail_release:
            raise GitHubCommandError('scripted release failure')
        self.releases.add(name)
        if self.move_tag_after_release:
            tag = self.tags[name]
            self.tags[name] = GitHubTag(name, 'moved-tag-object', 'moved-commit', tag.message)
        return GitHubRelease(name, f'https://example.invalid/{name}')

    def create_branch(self, name, commit):
        self.calls.append(('create-branch', name, commit))
        if self.fail_branch:
            raise GitHubCommandError('scripted branch failure')
        self.branch = commit

    def update_branch(self, name, commit, expected_old):
        self.calls.append(('update-branch', name, commit, expected_old))
        if self.fail_branch or self.branch != expected_old:
            raise GitHubCommandError('scripted guarded branch failure')
        self.branch = commit


class CreateExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        _dataset(self.root, '2026-01')
        self.factory = create_check_factory()

    def tearDown(self):
        self.temporary.cleanup()

    def _execute(self, git, github, **options):
        github.git = git
        plan = plan_create(
            '2026-01', self.root, Path('/target'), git=git, github=github,
            check_factory=self.factory, **options,
        )
        before = (dict(git.tags), git.branch, dict(github.tags), github.branch)
        self.assertIn('read-only', render_create_plan(plan))
        self.assertEqual(before, (dict(git.tags), git.branch, dict(github.tags), github.branch))
        result = execute_create(plan, git=git, github=github, check_factory=self.factory)
        return result

    def test_initial_creation_exact_write_order(self):
        git, github = ScriptedGit(self.root), ScriptedGitHub(self.root)
        result = self._execute(git, github)
        self.assertEqual(result.commit_id, 'new-commit')
        writes = [call[0] for call in git.calls + github.calls if call[0] in {
            'prepare-branch', 'create-commit', 'move-branch', 'create-tag', 'transfer-tag',
        }]
        self.assertEqual(writes, [
            'prepare-branch', 'create-commit', 'move-branch', 'create-tag', 'transfer-tag',
        ])
        self.assertIn('move-local-branch:2026-01:new-commit', result.completed)

    def test_replacement_uses_old_ids_and_guarded_transfer(self):
        base = _tag('2026-01', 'c1')
        git = ScriptedGit(self.root, (base,), 'c1')
        github = ScriptedGitHub(self.root, (base,), parents={'c1': ()})
        result = self._execute(git, github, replace=True)
        self.assertEqual((result.old_tag_object_id, result.old_commit_id), ('tag-2026-01', 'c1'))
        self.assertIn(('move-branch', '2026-01', 'new-commit', 'c1'), git.calls)
        self.assertIn(('transfer-tag', 'origin', '2026-01', 'tag-2026-01'), github.calls)

    def test_v2_and_v3_extend_latest_commit(self):
        for tags, expected_id, parent in (
            ((_tag('2026-01', 'c1'),), '2026-01v2', 'c1'),
            ((_tag('2026-01', 'c1'), _tag('2026-01v2', 'c2')), '2026-01v3', 'c2'),
        ):
            with self.subTest(expected_id=expected_id):
                _dataset(self.root, expected_id)
                parents = {'c1': (), 'c2': ('c1',)}
                git = ScriptedGit(self.root, tags, parent)
                github = ScriptedGitHub(self.root, tags, parents=parents)
                result = self._execute(git, github, revision=True)
                self.assertEqual(result.data_id, expected_id)
                subject = make_annotation(
                    expected_id, '2025-03', ('analysis.yaml',), 'analysis.yaml',
                ).subject
                self.assertIn(('create-commit', subject, parent), git.calls)

    def test_wrong_remote_object_is_reported_with_completed_state(self):
        git, github = ScriptedGit(self.root), ScriptedGitHub(self.root)
        github.wrong_transfer = True
        github.git = git
        plan = plan_create('2026-01', self.root, Path('/target'), git=git, github=github,
                           check_factory=self.factory)
        with self.assertRaises(ReleaseExecutionError) as caught:
            execute_create(plan, git=git, github=github, check_factory=self.factory)
        self.assertIn('remote tag ref failed', str(caught.exception))
        self.assertIn('move-local-branch:2026-01:new-commit', str(caught.exception))
        self.assertIn('create-annotated-tag', str(caught.exception))


class LocalMutationTests(unittest.TestCase):
    def test_commit_branch_and_annotated_tag_have_exact_tree_and_parent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target_root, source_root = root / 'target', root / 'source'
            remote_root = root / 'remote.git'
            target_root.mkdir()
            source_root.mkdir()
            _git(root, 'init', '--bare', str(remote_root))
            _git(target_root, 'init', '-b', 'main')
            _git(target_root, 'config', 'user.name', 'Test')
            _git(target_root, 'config', 'user.email', 'test@example.com')
            (target_root / 'old').write_text('old', encoding='utf-8')
            _git(target_root, 'add', 'old')
            _git(target_root, 'commit', '-m', 'unrelated main')
            _git(target_root, 'remote', 'add', 'origin', str(remote_root))
            _git(source_root, 'init', '-b', 'main')
            (source_root / 'new').write_text('new', encoding='utf-8')
            _git(source_root, 'add', 'new')
            _git(source_root, 'commit', '-m', 'source')
            source_commit = _git(source_root, 'rev-parse', 'HEAD')

            client = LocalGitClient()
            target = TargetCheckout(target_root, 'origin', 'git@github.com:eos/data.git')
            client.prepare_local_branch(target, '2026-01', None, None)
            with client.checkout_commit(source_root, source_commit) as checked_out:
                commit = client.create_dataset_commit(target, checked_out, 'dataset', None)
            self.assertEqual(client.commit_parents(target, commit), ())
            self.assertEqual(
                _git(target_root, 'rev-parse', f'{commit}^{{tree}}'),
                _git(source_root, 'rev-parse', f'{source_commit}^{{tree}}'),
            )
            client.update_local_branch(target, '2026-01', commit, None)
            self.assertEqual(sorted(path.name for path in target_root.iterdir() if path.name != '.git'), ['new'])
            message = make_annotation(
                '2026-01', '2025-03', ('analysis.yaml',), 'analysis.yaml',
            ).serialize()
            object_id, peeled = client.create_local_tag(target, '2026-01', commit, message, None)
            self.assertEqual(peeled, commit)
            GitHubClient(executable='unused').transfer_local_tag(
                target_root, 'origin', '2026-01', None,
            )
            self.assertEqual(
                _git(remote_root, 'rev-parse', 'refs/tags/2026-01'), object_id,
            )
            self.assertEqual(_git(remote_root, 'for-each-ref', '--format=%(refname)', 'refs/heads'), '')

    def test_local_tag_compare_and_swap_does_not_overwrite_a_race(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _git(root, 'init', '-b', 'main')
            _git(root, 'config', 'user.name', 'Test')
            _git(root, 'config', 'user.email', 'test@example.com')
            (root / 'value').write_text('one', encoding='utf-8')
            _git(root, 'add', 'value')
            _git(root, 'commit', '-m', 'one')
            first = _git(root, 'rev-parse', 'HEAD')
            (root / 'value').write_text('two', encoding='utf-8')
            _git(root, 'commit', '-am', 'two')
            second = _git(root, 'rev-parse', 'HEAD')
            _git(root, 'tag', '-a', '-m', 'old', '2026-01', first)
            old_object = _git(root, 'rev-parse', 'refs/tags/2026-01')

            class RacingGit(LocalGitClient):
                injected = False

                def _run(inner_self, arguments, **kwargs):
                    if (
                        arguments[:2] == ['update-ref', 'refs/tags/2026-01']
                        and not inner_self.injected
                    ):
                        inner_self.injected = True
                        _git(root, 'update-ref', 'refs/tags/2026-01', second, old_object)
                    return super()._run(arguments, **kwargs)

            target = TargetCheckout(root, 'origin', 'git@github.com:eos/data.git')
            with self.assertRaises(ReleasePlanningError):
                RacingGit().create_local_tag(
                    target, '2026-01', second, 'replacement metadata', old_object,
                )
            self.assertEqual(_git(root, 'rev-parse', 'refs/tags/2026-01'), second)


class PublishExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        _dataset(self.root, '2026-01')
        self.factory = create_check_factory()

    def tearDown(self):
        self.temporary.cleanup()

    def _plan(self, tags, *, data_id='2026-01', branch=None, parents=None):
        if isinstance(tags, GitHubTag):
            tags = (tags,)
        git = ScriptedGit(self.root)
        github = ScriptedGitHub(self.root, tags, branch=branch, parents=parents)
        plan = plan_publish(data_id, Path('/target'), git=git, github=github,
                            check_factory=self.factory)
        before = (set(github.releases), github.branch, list(github.calls))
        self.assertIn('read-only', render_publish_plan(plan))
        self.assertEqual((set(github.releases), github.branch), before[:2])
        return git, github, plan

    def test_first_publish_releases_verifies_then_creates_branch(self):
        git, github, plan = self._plan(_tag('2026-01', 'c1'), parents={'c1': ()})
        result = execute_publish(plan, git=git, github=github, check_factory=self.factory)
        writes = [call[0] for call in github.calls if call[0] in {'create-release', 'create-branch'}]
        self.assertEqual(writes, ['create-release', 'create-branch'])
        self.assertEqual(result.branch_commit_id, 'c1')

    def test_later_publish_advances_only_an_ancestor_branch(self):
        tags = (_tag('2026-01', 'c1'), _tag('2026-01v2', 'c2'))
        git, github, plan = self._plan(
            tags, data_id='2026-01v2', branch='c1', parents={'c1': (), 'c2': ('c1',)},
        )
        execute_publish(plan, git=git, github=github, check_factory=self.factory)
        self.assertIn(('update-branch', '2026-01', 'c2', 'c1'), github.calls)

    def test_release_failure_never_calls_branch_write(self):
        git, github, plan = self._plan(_tag('2026-01', 'c1'), parents={'c1': ()})
        github.fail_release = True
        with self.assertRaises(ReleaseExecutionError) as caught:
            execute_publish(plan, git=git, github=github, check_factory=self.factory)
        self.assertIn('branch was not touched', str(caught.exception))
        self.assertFalse(any(call[0] in {'create-branch', 'update-branch'} for call in github.calls))

    def test_post_release_tag_movement_is_partial_failure_without_branch_write(self):
        git, github, plan = self._plan(_tag('2026-01', 'c1'), parents={'c1': ()})
        github.move_tag_after_release = True
        with self.assertRaises(PartialPublishError) as caught:
            execute_publish(plan, git=git, github=github, check_factory=self.factory)
        self.assertIn("release '2026-01' exists", str(caught.exception))
        self.assertIn('tag ref mismatch', str(caught.exception))
        self.assertFalse(any(call[0] in {'create-branch', 'update-branch'} for call in github.calls))

    def test_branch_failure_reports_release_observed_target_and_recovery(self):
        git, github, plan = self._plan(_tag('2026-01', 'c1'), parents={'c1': ()})
        github.fail_branch = True
        with self.assertRaises(PartialPublishError) as caught:
            execute_publish(plan, git=git, github=github, check_factory=self.factory)
        message = str(caught.exception)
        for text in ("release '2026-01' exists", 'refs/heads/2026-01', 'c1', 'observed absent',
                     'set that ref explicitly'):
            self.assertIn(text, message)


if __name__ == '__main__':
    unittest.main()
