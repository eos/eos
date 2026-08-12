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

"""A thin, replaceable wrapper around the GitHub CLI for EOS dataset releases."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any


class GitHubError(RuntimeError):
    """A focused failure while reading GitHub state."""


class GitHubCLIUnavailable(GitHubError):
    pass


class GitHubCommandError(GitHubError):
    pass


class GitHubResponseError(GitHubError):
    pass


@dataclass(frozen=True)
class GitHubRef:
    name: str
    object_id: str
    object_type: str


@dataclass(frozen=True)
class GitHubTag:
    name: str
    object_id: str
    commit_id: str
    message: str


@dataclass(frozen=True)
class GitHubRelease:
    tag_name: str
    url: str


Runner = Callable[..., subprocess.CompletedProcess]


class GitHubClient:
    """Typed read operations implemented with ``gh`` argument arrays."""

    def __init__(
        self,
        repository: str = 'eos/data',
        *,
        executable: str | None = None,
        runner: Runner = subprocess.run,
        git_runner: Runner = subprocess.run,
    ) -> None:
        self.repository = repository
        self._executable = executable
        self._runner = runner
        self._git_runner = git_runner

    def _gh(self) -> str:
        executable = self._executable or shutil.which('gh')
        if executable is None:
            raise GitHubCLIUnavailable(
                "GitHub CLI 'gh' is required for release planning but was not found on PATH"
            )
        return executable

    def _run(self, arguments: Sequence[str]) -> subprocess.CompletedProcess:
        command = [self._gh(), *arguments]
        try:
            result = self._runner(command, capture_output=True, text=True, check=False)
        except OSError as error:
            raise GitHubCommandError(f"cannot execute GitHub CLI: {error}") from error
        if result.returncode != 0:
            diagnostic = (result.stderr or result.stdout or '').strip()
            suffix = f': {diagnostic}' if diagnostic else ''
            raise GitHubCommandError(f"GitHub CLI command failed ({result.returncode}){suffix}")
        return result

    def _api_json(
        self,
        endpoint: str,
        *,
        allow_missing: bool = False,
        paginate: bool = False,
        method: str | None = None,
        fields: Mapping[str, Any] | None = None,
    ) -> Any:
        command = [self._gh(), 'api', endpoint]
        if method is not None:
            command.extend(('--method', method))
        if paginate:
            command.extend(('--paginate', '--slurp'))
        for key, value in (fields or {}).items():
            option = '-F' if isinstance(value, (bool, int)) else '-f'
            rendered = str(value).lower() if isinstance(value, bool) else str(value)
            command.extend((option, f'{key}={rendered}'))
        try:
            result = self._runner(command, capture_output=True, text=True, check=False)
        except OSError as error:
            raise GitHubCommandError(f"cannot execute GitHub CLI: {error}") from error
        if allow_missing and result.returncode == 1 and '404' in (result.stderr or ''):
            return None
        if result.returncode != 0:
            diagnostic = (result.stderr or result.stdout or '').strip()
            suffix = f': {diagnostic}' if diagnostic else ''
            raise GitHubCommandError(f"GitHub API read failed ({result.returncode}){suffix}")
        try:
            return json.loads(result.stdout)
        except (TypeError, json.JSONDecodeError) as error:
            raise GitHubResponseError(f'malformed JSON from GitHub for {endpoint}: {error}') from error

    @staticmethod
    def _mapping(value: Any, description: str) -> Mapping[str, Any]:
        if not isinstance(value, Mapping):
            raise GitHubResponseError(f'unexpected GitHub response for {description}: expected object')
        return value

    @staticmethod
    def _string(mapping: Mapping[str, Any], field: str, description: str) -> str:
        value = mapping.get(field)
        if not isinstance(value, str) or not value:
            raise GitHubResponseError(
                f'unexpected GitHub response for {description}: '
                f'{field} must be a nonempty string'
            )
        return value

    def list_tag_refs(self) -> tuple[GitHubRef, ...]:
        document = self._api_json(
            f'repos/{self.repository}/git/matching-refs/tags/', paginate=True,
        )
        if not isinstance(document, list):
            raise GitHubResponseError('unexpected paginated GitHub tag-ref response: expected list')
        refs = []
        for page in document:
            if not isinstance(page, list):
                raise GitHubResponseError(
                    'unexpected paginated GitHub tag-ref response: expected each page to be a list'
                )
            for item in page:
                entry = self._mapping(item, 'tag ref')
                target = self._mapping(entry.get('object'), 'tag ref object')
                full_name = self._string(entry, 'ref', 'tag ref')
                object_id = self._string(target, 'sha', 'tag ref object')
                object_type = self._string(target, 'type', 'tag ref object')
                if not full_name.startswith('refs/tags/') or full_name == 'refs/tags/':
                    raise GitHubResponseError(
                        f"unexpected GitHub tag ref name: '{full_name}'"
                    )
                refs.append(GitHubRef(
                    full_name.removeprefix('refs/tags/'),
                    object_id,
                    object_type,
                ))
        return tuple(sorted(refs, key=lambda ref: ref.name))

    def get_tag(self, ref: GitHubRef) -> GitHubTag:
        if ref.object_type != 'tag':
            raise GitHubResponseError(
                f"tag '{ref.name}' is lightweight or has unexpected type '{ref.object_type}'"
            )
        document = self._mapping(
            self._api_json(f'repos/{self.repository}/git/tags/{ref.object_id}'),
            'annotated tag',
        )
        target = self._mapping(document.get('object'), 'annotated tag target')
        if target.get('type') != 'commit':
            raise GitHubResponseError(f"annotated tag '{ref.name}' does not point to a commit")
        commit_id = self._string(target, 'sha', 'annotated tag target')
        message = document.get('message')
        if not isinstance(message, str):
            raise GitHubResponseError('annotated-tag response is missing a string message')
        return GitHubTag(ref.name, ref.object_id, commit_id, message)

    def get_release(self, tag_name: str) -> GitHubRelease | None:
        document = self._api_json(f'repos/{self.repository}/releases/tags/{tag_name}', allow_missing=True)
        if document is None:
            return None
        entry = self._mapping(document, 'release')
        actual_tag = entry.get('tag_name')
        url = entry.get('html_url')
        if actual_tag != tag_name or not isinstance(url, str):
            raise GitHubResponseError(f"unexpected release response for tag '{tag_name}'")
        return GitHubRelease(tag_name, url)

    def get_branch(self, branch_name: str) -> str | None:
        document = self._api_json(
            f'repos/{self.repository}/git/ref/heads/{branch_name}', allow_missing=True,
        )
        if document is None:
            return None
        entry = self._mapping(document, 'branch ref')
        target = self._mapping(entry.get('object'), 'branch ref object')
        if target.get('type') != 'commit' or not isinstance(target.get('sha'), str):
            raise GitHubResponseError(f"unexpected branch response for '{branch_name}'")
        return target['sha']

    def commit_parents(self, commit_id: str) -> tuple[str, ...]:
        document = self._mapping(self._api_json(f'repos/{self.repository}/git/commits/{commit_id}'), 'commit')
        parents = document.get('parents')
        if not isinstance(parents, list):
            raise GitHubResponseError(f"commit response for '{commit_id}' has no parent list")
        result = []
        for parent in parents:
            entry = self._mapping(parent, 'commit parent')
            if not isinstance(entry.get('sha'), str):
                raise GitHubResponseError('commit parent is missing sha')
            result.append(entry['sha'])
        return tuple(result)

    def is_ancestor(self, ancestor: str, descendant: str) -> bool:
        document = self._mapping(
            self._api_json(f'repos/{self.repository}/compare/{ancestor}...{descendant}'),
            'commit comparison',
        )
        status = document.get('status')
        if status not in ('ahead', 'identical', 'behind', 'diverged'):
            raise GitHubResponseError('commit-comparison response has an unexpected status')
        return status in ('ahead', 'identical')

    @contextmanager
    def checkout_tree(self, commit_id: str):
        """Check out a remote committed tree using read-only ``gh repo clone``."""
        temporary = Path(tempfile.mkdtemp(prefix='eos-data-publish-'))
        try:
            self._run(['repo', 'clone', self.repository, str(temporary), '--', '--no-checkout'])
            try:
                result = self._runner(
                    ['git', 'checkout', '--detach', commit_id], cwd=temporary,
                    capture_output=True, text=True, check=False,
                )
            except OSError as error:
                raise GitHubCommandError(f'cannot execute Git while checking out tag: {error}') from error
            if result.returncode != 0:
                diagnostic = (result.stderr or result.stdout or '').strip()
                raise GitHubCommandError(f'cannot check out tagged commit {commit_id}: {diagnostic}')
            yield temporary
        finally:
            shutil.rmtree(temporary, ignore_errors=True)

    def transfer_local_tag(
        self,
        root: Path,
        remote_name: str,
        tag_name: str,
        expected_old: str | None,
    ) -> None:
        arguments = ['git', 'push', '--porcelain']
        if expected_old is not None:
            arguments.append(f'--force-with-lease=refs/tags/{tag_name}:{expected_old}')
        arguments.extend((remote_name, f'refs/tags/{tag_name}:refs/tags/{tag_name}'))
        try:
            result = self._git_runner(
                arguments, cwd=root, capture_output=True, text=True, check=False,
            )
        except OSError as error:
            raise GitHubCommandError(f'cannot execute Git while transferring tag: {error}') from error
        if result.returncode != 0:
            diagnostic = (result.stderr or result.stdout or '').strip()
            suffix = f': {diagnostic}' if diagnostic else ''
            raise GitHubCommandError(f'tag transfer failed ({result.returncode}){suffix}')

    def create_release(self, tag_name: str, title: str, notes: str) -> GitHubRelease:
        document = self._mapping(self._api_json(
            f'repos/{self.repository}/releases',
            method='POST',
            fields={'tag_name': tag_name, 'name': title, 'body': notes},
        ), 'created release')
        actual_tag = document.get('tag_name')
        url = document.get('html_url')
        if actual_tag != tag_name or not isinstance(url, str) or not url:
            raise GitHubResponseError('created release does not target the requested tag')
        return GitHubRelease(tag_name, url)

    def create_branch(self, branch_name: str, commit_id: str) -> None:
        self._create_ref(f'refs/heads/{branch_name}', commit_id)

    def update_branch(self, branch_name: str, commit_id: str, expected_old: str) -> None:
        temporary = Path(tempfile.mkdtemp(prefix='eos-data-publish-branch-'))
        try:
            self._run(['repo', 'clone', self.repository, str(temporary), '--', '--no-checkout'])
            command = [
                'git', 'push', '--porcelain',
                f'--force-with-lease=refs/heads/{branch_name}:{expected_old}',
                'origin', f'{commit_id}:refs/heads/{branch_name}',
            ]
            try:
                result = self._git_runner(
                    command, cwd=temporary, capture_output=True, text=True, check=False,
                )
            except OSError as error:
                raise GitHubCommandError(
                    f'cannot execute Git while advancing branch: {error}'
                ) from error
            if result.returncode != 0:
                diagnostic = (result.stderr or result.stdout or '').strip()
                suffix = f': {diagnostic}' if diagnostic else ''
                raise GitHubCommandError(
                    f'guarded branch update failed ({result.returncode}){suffix}'
                )
        finally:
            shutil.rmtree(temporary, ignore_errors=True)

    def transfer_local_branch(
        self,
        root: Path,
        branch_name: str,
        commit_id: str,
        main_commit_id: str,
        tag_name: str,
        tag_object_id: str,
    ) -> None:
        """Atomically push a new branch after confirming the main and tag refs are unchanged."""
        command = [
            'git', 'push', '--porcelain', '--atomic',
            f'--force-with-lease=refs/heads/main:{main_commit_id}',
            f'--force-with-lease=refs/tags/{tag_name}:{tag_object_id}',
            f'--force-with-lease=refs/heads/{branch_name}:',
            'origin',
            f'{main_commit_id}:refs/heads/main',
            f'{tag_object_id}:refs/tags/{tag_name}',
            f'{commit_id}:refs/heads/{branch_name}',
        ]
        try:
            result = self._git_runner(
                command, cwd=root, capture_output=True, text=True, check=False,
            )
        except OSError as error:
            raise GitHubCommandError(
                f'cannot execute Git while transferring registration branch: {error}'
            ) from error
        if result.returncode != 0:
            diagnostic = (result.stderr or result.stdout or '').strip()
            suffix = f': {diagnostic}' if diagnostic else ''
            raise GitHubCommandError(
                f'registration branch transfer failed ({result.returncode}){suffix}'
            )

    def _create_ref(self, ref_name: str, object_id: str) -> None:
        document = self._mapping(self._api_json(
            f'repos/{self.repository}/git/refs',
            method='POST', fields={'ref': ref_name, 'sha': object_id},
        ), 'created ref')
        target = self._mapping(document.get('object'), 'created ref object')
        if document.get('ref') != ref_name or target.get('sha') != object_id:
            raise GitHubResponseError(f"created ref '{ref_name}' has an unexpected target")
