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
from unittest import mock

from eos.cli.data import create_check_factory
from eos.cli.data_github import (
    GitHubCLIUnavailable, GitHubClient, GitHubCommandError, GitHubRef,
    GitHubRelease, GitHubResponseError, GitHubTag,
)
from eos.cli.data_release import (
    LocalGitClient, OperationType, ReleasePlanningError, ResolvedSource,
    TargetCheckout, github_repository, make_annotation, parse_annotation,
    parse_base_id, parse_data_id, plan_create, plan_publish,
    render_create_plan, render_publish_plan, validate_family,
)


def _run_git(root, *arguments):
    return subprocess.run(
        ['git', *arguments], cwd=root, check=True, capture_output=True, text=True,
        env={**os.environ, 'GIT_AUTHOR_NAME': 'Test', 'GIT_AUTHOR_EMAIL': 'test@example.com',
             'GIT_COMMITTER_NAME': 'Test', 'GIT_COMMITTER_EMAIL': 'test@example.com'},
    ).stdout.strip()


def _dataset(root: Path, data_id='2026-01'):
    (root / 'README.md').write_text('# Dataset\n', encoding='utf-8')
    (root / 'analysis.yaml').write_text(
        'metadata:\n  id: 2025-03\n  title: Main title\n'
        '  description: Release notes.\n  authors:\n    - name: Jane Doe\n',
        encoding='utf-8',
    )
    (root / '.zenodo.json').write_text(json.dumps({
        'title': f'EOS/DATA-{data_id}: Supplementary material for EOS/ANALYSIS-2025-03',
        'description': 'Release notes.', 'creators': [{'name': 'Doe, Jane', 'affiliation': 'Lab'}],
        'upload_type': 'dataset', 'license': 'CC-BY-4.0', 'grants': [],
    }), encoding='utf-8')


class FakeGit:
    def __init__(self, source: Path, local_tags=None, local_branch=None):
        self.source = source
        self.tags = local_tags or {}
        self.branch = local_branch

    def target_checkout(self, directory, *, require_clean=True):
        return TargetCheckout(Path(directory), 'origin', 'git@github.com:eos/data.git')

    def local_tag_refs(self, target):
        return self.tags

    def local_branch(self, target, name):
        return self.branch

    @contextmanager
    def checkout_source(self, source):
        yield ResolvedSource(str(source), 'source-commit', self.source)


class FakeGitHub:
    def __init__(self, root: Path, tags=(), releases=(), branch=None, parents=None):
        self.root = root
        self.tags = {tag.name: tag for tag in tags}
        self.releases = set(releases)
        self.branch = branch
        self.parents = parents or {}
        self.calls = []

    def list_tag_refs(self):
        self.calls.append(('list-tag-refs',))
        return tuple(GitHubRef(tag.name, tag.object_id, 'tag') for tag in self.tags.values())

    def get_tag(self, ref):
        self.calls.append(('get-tag', ref.name))
        return self.tags[ref.name]

    def get_release(self, name):
        self.calls.append(('get-release', name))
        return GitHubRelease(name, 'https://example.invalid') if name in self.releases else None

    def get_branch(self, name):
        self.calls.append(('get-branch', name))
        return self.branch

    def commit_parents(self, commit):
        self.calls.append(('parents', commit))
        return tuple(self.parents.get(commit, ()))

    def is_ancestor(self, ancestor, descendant):
        self.calls.append(('is-ancestor', ancestor, descendant))
        return ancestor == descendant or ancestor in self.parents.get(descendant, ())

    @contextmanager
    def checkout_tree(self, commit):
        self.calls.append(('checkout', commit))
        yield self.root


def _tag(name, commit, analysis_id='2025-03', files=('analysis.yaml',), main='analysis.yaml'):
    annotation = make_annotation(name, analysis_id, files, main)
    return GitHubTag(name, f'tag-{name}', commit, annotation.serialize())


class IDAndAnnotationTests(unittest.TestCase):
    def test_ids(self):
        self.assertEqual(parse_base_id('2026-01'), '2026-01')
        self.assertEqual(parse_data_id('2026-01'), ('2026-01', 1))
        self.assertEqual(parse_data_id('2026-01v10'), ('2026-01', 10))
        for value in ('2026-1', '2026-01v2'):
            with self.assertRaises(ReleasePlanningError):
                parse_base_id(value)
        for value in ('2026-01v1', '2026-01v01', '2026-01v02', 'x'):
            with self.assertRaises(ReleasePlanningError):
                parse_data_id(value)

    def test_annotation_round_trip_is_sorted_and_retains_unknown_trailers(self):
        value = make_annotation('2026-01v2', '2025-03', ['z.yaml', 'a.yaml'], 'z.yaml')
        parsed = parse_annotation(value.serialize() + '\nFuture-Key: retained', expected_data_id='2026-01v2')
        self.assertEqual(parsed.analysis_files, (Path('a.yaml'), Path('z.yaml')))
        self.assertEqual(parsed.unknown_trailers, (('Future-Key', 'retained'),))

    def test_annotation_rejects_malformed_metadata(self):
        valid = make_annotation('2026-01', '2025-03', ['analysis.yaml'], 'analysis.yaml').serialize()
        cases = (
            valid.replace('Version: 1', 'Version: 2'),
            valid.replace('EOS-Analysis-ID: 2025-03\n', ''),
            valid + '\nEOS-Analysis-ID: 2025-03',
            valid + '\nEOS-Analysis-File: analysis.yaml',
            valid.replace('analysis.yaml', '../analysis.yaml'),
            valid.replace('EOS-Main-Analysis-File: analysis.yaml', 'EOS-Main-Analysis-File: other.yaml'),
            valid.replace('DATA-2026-01', 'DATA-2026-02'),
            valid + '\nnot a trailer',
        )
        for message in cases:
            with self.subTest(message=message[-50:]), self.assertRaises(ReleasePlanningError):
                parse_annotation(message, expected_data_id='2026-01')

    def test_supported_github_urls(self):
        for url in (
            'git@github.com:eos/data.git', 'ssh://git@github.com/eos/data.git',
            'ssh://github.com/eos/data', 'https://github.com/eos/data.git',
        ):
            self.assertEqual(github_repository(url), 'eos/data')
        for url in (
            'https://gitlab.com/eos/data', 'ssh://github/eos/data',
            'ssh://attacker@github.com/eos/data',
            'https://attacker@github.com/eos/data',
            'https://attacker:secret@github.com/eos/data',
        ):
            self.assertIsNone(github_repository(url))


class FamilyTests(unittest.TestCase):
    def test_consecutive_lineage(self):
        remote = (_tag('2026-01', 'c1'), _tag('2026-01v2', 'c2'))
        family = validate_family('2026-01', {}, remote, lambda commit: {'c1': (), 'c2': ('c1',)}[commit])
        self.assertEqual(family.next_id, '2026-01v3')

    def test_gaps_malformed_conflicts_and_lineage_fail(self):
        cases = (
            ({}, (_tag('2026-01', 'c1'), _tag('2026-01v3', 'c3')), {'c1': (), 'c3': ('c1',)}),
            ({}, (_tag('2026-01', 'c1'), GitHubTag('2026-01v1', 't1', 'c2', 'x')), {'c1': (), 'c2': ('c1',)}),
            ({'2026-01': ('different', 'c1')}, (_tag('2026-01', 'c1'),), {'c1': ()}),
            ({}, (_tag('2026-01', 'c1'), _tag('2026-01v2', 'c2')), {'c1': (), 'c2': ('other',)}),
        )
        for local, remote, parents in cases:
            with self.subTest(remote=[tag.name for tag in remote]), self.assertRaises(ReleasePlanningError):
                validate_family('2026-01', local, remote, lambda commit: parents[commit])


class LocalGitTests(unittest.TestCase):
    def test_local_source_uses_committed_main_and_excludes_dirty_content(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _run_git(root, 'init', '-b', 'main')
            (root / 'value').write_text('committed', encoding='utf-8')
            _run_git(root, 'add', 'value')
            _run_git(root, 'commit', '-m', 'initial')
            commit = _run_git(root, 'rev-parse', 'HEAD')
            (root / 'value').write_text('dirty', encoding='utf-8')
            (root / 'untracked').write_text('no', encoding='utf-8')
            client = LocalGitClient()
            with client.checkout_source(root) as source:
                self.assertEqual(source.commit_id, commit)
                self.assertEqual((source.root / 'value').read_text(), 'committed')
                self.assertFalse((source.root / 'untracked').exists())

    def test_remote_source_resolves_remote_main(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'source'
            root.mkdir()
            _run_git(root, 'init', '-b', 'main')
            (root / 'value').write_text('remote committed', encoding='utf-8')
            _run_git(root, 'add', 'value')
            _run_git(root, 'commit', '-m', 'initial')
            commit = _run_git(root, 'rev-parse', 'HEAD')
            url = root.as_uri()
            with LocalGitClient().checkout_source(url) as source:
                self.assertEqual(source.description, url)
                self.assertEqual(source.commit_id, commit)
                self.assertEqual((source.root / 'value').read_text(), 'remote committed')

    def test_target_requires_root_identity_and_clean_state(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _run_git(root, 'init', '-b', 'main')
            _run_git(root, 'remote', 'add', 'origin', 'https://github.com/eos/data.git')
            target = LocalGitClient().target_checkout(root)
            self.assertEqual(target.remote_name, 'origin')
            child = root / 'child'
            child.mkdir()
            with self.assertRaises(ReleasePlanningError):
                LocalGitClient().target_checkout(child)
            (root / 'untracked').write_text('unsafe', encoding='utf-8')
            with self.assertRaises(ReleasePlanningError):
                LocalGitClient().target_checkout(root)

    def test_target_requires_push_url_to_identify_eos_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _run_git(root, 'init', '-b', 'main')
            _run_git(root, 'remote', 'add', 'origin', 'https://github.com/eos/data.git')
            _run_git(root, 'remote', 'set-url', '--push', 'origin', 'https://github.com/attacker/data.git')
            with self.assertRaises(ReleasePlanningError):
                LocalGitClient().target_checkout(root)

            _run_git(root, 'remote', 'set-url', '--push', 'origin', 'https://github.com/eos/data.git')
            _run_git(
                root, 'remote', 'set-url', '--add', '--push', 'origin',
                'https://github.com/attacker/data.git',
            )
            with self.assertRaises(ReleasePlanningError):
                LocalGitClient().target_checkout(root)

    def test_remote_source_is_terminated_before_user_input(self):
        runner = mock.Mock(return_value=subprocess.CompletedProcess([], 2, '', 'expected failure'))
        with self.assertRaises(ReleasePlanningError):
            with LocalGitClient(runner=runner).checkout_source('-malicious-option'):
                pass
        command = runner.call_args.args[0]
        self.assertEqual(command[1:5], ['clone', '--no-checkout', '--', '-malicious-option'])


class PlannerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        _dataset(self.root)
        self.factory = create_check_factory()

    def tearDown(self):
        self.temporary.cleanup()

    def test_initial_plan_is_deterministic_and_read_only(self):
        github = FakeGitHub(self.root)
        plan = plan_create(
            '2026-01', self.root, Path('/target'), git=FakeGit(self.root), github=github,
            check_factory=self.factory,
        )
        self.assertEqual(plan.data_id, '2026-01')
        self.assertIsNone(plan.parent_commit_id)
        self.assertEqual([operation.kind for operation in plan.operations], [
            OperationType.PREPARE_LOCAL_BRANCH, OperationType.CREATE_COMMIT,
            OperationType.MOVE_LOCAL_BRANCH, OperationType.CREATE_ANNOTATED_TAG,
            OperationType.TRANSFER_TAG,
        ])
        self.assertIsNone(plan.operations[2].expected_old)
        self.assertEqual(plan.operations[2].target, '<new-commit>')
        rendered = render_create_plan(plan)
        for text in ('Create dataset/tag plan', 'Source commit: source-commit', 'Parent/ancestry: orphan',
                     'Analysis files: analysis.yaml', 'EOS-Data-Metadata-Version: 1'):
            self.assertIn(text, rendered)
        self.assertTrue(all(call[0] in {'list-tag-refs', 'get-branch'} for call in github.calls))

    def test_revision_and_replacement_policy(self):
        base = _tag('2026-01', 'c1')
        parents = {'c1': ()}
        github = FakeGitHub(self.root, (base,), parents=parents)
        _dataset(self.root, '2026-01v2')
        revision = plan_create(
            '2026-01', self.root, Path('/target'), revision=True,
            git=FakeGit(self.root, {'2026-01': ('tag-2026-01', 'c1')}, 'c1'), github=github,
            check_factory=self.factory,
        )
        self.assertEqual((revision.data_id, revision.parent_commit_id), ('2026-01v2', 'c1'))

        _dataset(self.root, '2026-01')
        replacement = plan_create(
            '2026-01', self.root, Path('/target'), replace=True,
            git=FakeGit(self.root, {'2026-01': ('tag-2026-01', 'c1')}, 'c1'), github=github,
            check_factory=self.factory,
        )
        self.assertEqual(replacement.old_tag_object_id, 'tag-2026-01')
        self.assertTrue(replacement.operations[-1].force)

    def test_replacement_disqualifying_states(self):
        base = _tag('2026-01', 'c1')
        revision = _tag('2026-01v2', 'c2')
        cases = (
            FakeGitHub(self.root, (base,), releases=('2026-01',), parents={'c1': ()}),
            FakeGitHub(self.root, (base,), branch='c1', parents={'c1': ()}),
            FakeGitHub(self.root, (base, revision), parents={'c1': (), 'c2': ('c1',)}),
        )
        for github in cases:
            with self.subTest(calls=github.calls), self.assertRaises(ReleasePlanningError):
                plan_create(
                    '2026-01', self.root, Path('/target'), replace=True,
                    git=FakeGit(self.root), github=github, check_factory=self.factory,
                )

    def test_publish_recovers_selection_and_orders_release_before_branch(self):
        tag = _tag('2026-01', 'c1')
        github = FakeGitHub(self.root, (tag,), parents={'c1': ()})
        plan = plan_publish(
            '2026-01', Path('/target'), git=FakeGit(self.root), github=github,
            check_factory=self.factory,
        )
        self.assertEqual((plan.release_title, plan.release_notes), ('Main title', 'Release notes.'))
        self.assertEqual([operation.kind for operation in plan.operations], [
            OperationType.CREATE_RELEASE, OperationType.CREATE_REMOTE_BRANCH,
        ])
        rendered = render_publish_plan(plan)
        self.assertLess(rendered.index('create-release'), rendered.index('create-remote-branch'))
        self.assertEqual(github.calls.count(('list-tag-refs',)), 1)
        self.assertEqual(github.calls.count(('get-tag', '2026-01')), 1)

    def test_publish_rejects_existing_release_and_unrelated_or_backward_branch(self):
        base = _tag('2026-01', 'c1')
        revision = _tag('2026-01v2', 'c2')
        cases = (
            ('2026-01', FakeGitHub(
                self.root, (base,), releases=('2026-01',), parents={'c1': ()},
            )),
            ('2026-01', FakeGitHub(
                self.root, (base,), branch='unrelated', parents={'c1': (), 'unrelated': ()},
            )),
            ('2026-01', FakeGitHub(
                self.root, (base, revision), branch='c2', parents={'c1': (), 'c2': ('c1',)},
            )),
        )
        for data_id, github in cases:
            with self.subTest(data_id=data_id, branch=github.branch), self.assertRaises(ReleasePlanningError):
                plan_publish(
                    data_id, Path('/target'), git=FakeGit(self.root), github=github,
                    check_factory=self.factory,
                )

    def test_planning_failure_does_not_mutate(self):
        _dataset(self.root, '2026-99')
        github = FakeGitHub(self.root)
        with self.assertRaises(ReleasePlanningError):
            plan_create(
                '2026-01', self.root, Path('/target'), git=FakeGit(self.root), github=github,
                check_factory=self.factory,
            )
        self.assertTrue(all(call[0] in {'list-tag-refs', 'get-branch'} for call in github.calls))

    def test_create_checks_source_before_github_preflight(self):
        events = []

        class OrderingGit(FakeGit):
            @contextmanager
            def checkout_source(inner_self, source):
                events.append('check-source')
                yield ResolvedSource(str(source), 'source-commit', inner_self.source)

        class OrderingGitHub(FakeGitHub):
            def list_tag_refs(inner_self):
                self.assertEqual(events, ['check-source'])
                events.append('github-preflight')
                return super().list_tag_refs()

        plan_create(
            '2026-01', self.root, Path('/target'), git=OrderingGit(self.root),
            github=OrderingGitHub(self.root), check_factory=self.factory,
        )
        self.assertEqual(events, ['check-source', 'github-preflight'])


class GitHubClientTests(unittest.TestCase):
    def test_missing_gh(self):
        with mock.patch('shutil.which', return_value=None):
            with self.assertRaises(GitHubCLIUnavailable):
                GitHubClient().list_tag_refs()

    def test_lightweight_tag_is_rejected_without_an_api_call(self):
        client = GitHubClient(executable='gh', runner=mock.Mock())
        with self.assertRaises(GitHubResponseError):
            client.get_tag(GitHubRef('2026-01', 'commit-id', 'commit'))
        client._runner.assert_not_called()

    def test_malformed_json_command_failure_and_unexpected_response(self):
        def result(returncode=0, stdout='', stderr=''):
            return subprocess.CompletedProcess([], returncode, stdout, stderr)
        for response, error in (
            (result(stdout='not json'), GitHubResponseError),
            (result(returncode=2, stderr='denied'), GitHubCommandError),
            (result(stdout='{}'), GitHubResponseError),
        ):
            with self.subTest(error=error):
                client = GitHubClient(executable='gh', runner=mock.Mock(return_value=response))
                with self.assertRaises(error):
                    client.list_tag_refs()

    def test_tag_ref_enumeration_paginates_and_flattens_all_pages(self):
        pages = [[
            {'ref': 'refs/tags/2026-01', 'object': {'sha': 'tag-1', 'type': 'tag'}},
        ], [
            {'ref': 'refs/tags/2026-01v2', 'object': {'sha': 'tag-2', 'type': 'tag'}},
        ]]
        runner = mock.Mock(return_value=subprocess.CompletedProcess([], 0, json.dumps(pages), ''))
        refs = GitHubClient(executable='gh', runner=runner).list_tag_refs()
        self.assertEqual([ref.name for ref in refs], ['2026-01', '2026-01v2'])
        self.assertEqual(runner.call_args.args[0][-2:], ['--paginate', '--slurp'])

    def test_tag_ref_response_requires_nonempty_string_fields(self):
        valid = {'ref': 'refs/tags/2026-01', 'object': {'sha': 'tag-1', 'type': 'tag'}}
        cases = []
        for field in ('ref',):
            for value in (None, ''):
                item = {**valid, field: value}
                cases.append(item)
        for field in ('sha', 'type'):
            for value in (None, ''):
                item = {**valid, 'object': {**valid['object'], field: value}}
                cases.append(item)
        cases.append({'object': valid['object']})
        cases.append({'ref': valid['ref'], 'object': {'type': 'tag'}})
        for item in cases:
            with self.subTest(item=item):
                response = subprocess.CompletedProcess([], 0, json.dumps([[item]]), '')
                client = GitHubClient(executable='gh', runner=mock.Mock(return_value=response))
                with self.assertRaises(GitHubResponseError):
                    client.list_tag_refs()

    def test_release_and_guarded_branch_update_use_exact_api_arguments(self):
        release = {
            'tag_name': '2026-01', 'html_url': 'https://example.invalid/release',
        }
        responses = [
            subprocess.CompletedProcess([], 0, json.dumps(release), ''),
            subprocess.CompletedProcess([], 0, '', ''),
        ]
        runner = mock.Mock(side_effect=responses)
        git_runner = mock.Mock(return_value=subprocess.CompletedProcess([], 0, '', ''))
        client = GitHubClient(executable='gh', runner=runner, git_runner=git_runner)
        result = client.create_release('2026-01', 'Title', 'Notes')
        self.assertEqual(result.tag_name, '2026-01')
        client.update_branch('2026-01', 'new', 'old')
        release_command = runner.call_args_list[0].args[0]
        self.assertEqual(release_command[:5], [
            'gh', 'api', 'repos/eos/data/releases', '--method', 'POST',
        ])
        self.assertIn('tag_name=2026-01', release_command)
        self.assertEqual(runner.call_args_list[1].args[0][:4], [
            'gh', 'repo', 'clone', 'eos/data',
        ])
        update_command = git_runner.call_args.args[0]
        self.assertEqual(update_command[:3], ['git', 'push', '--porcelain'])
        self.assertIn('--force-with-lease=refs/heads/2026-01:old', update_command)
        self.assertIn('new:refs/heads/2026-01', update_command)


if __name__ == '__main__':
    unittest.main()
