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

"""Pure release planning for ``eos-data create`` and ``publish``."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Enum
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tempfile
from types import MappingProxyType
from typing import Any, Protocol
from urllib.parse import urlparse

import yaml

from .data_checks import (
    CacheKey, CheckContext, CheckResult, CheckScope, ExitStatus, PlainTextRenderer, run_checks,
)
from .data_github import GitHubRef, GitHubRelease, GitHubTag


_BASE_ID = re.compile(r'^[0-9]{4}-[0-9]{2}$')
_DATA_ID = re.compile(r'^(?P<base>[0-9]{4}-[0-9]{2})(?:v(?P<revision>[1-9][0-9]*))?$')
_SUBJECT = re.compile(
    r'^EOS/DATA-(?P<data_id>[0-9]{4}-[0-9]{2}(?:v(?:[2-9]|[1-9][0-9]+))?): '
    r'Supplementary material for EOS/ANALYSIS-(?P<analysis_id>[0-9]{4}-[0-9]{2})$'
)
_TRAILER = re.compile(r'^(?P<key>[A-Za-z0-9-]+): (?P<value>.*)$')


class ReleasePlanningError(RuntimeError):
    pass


class SourceCheckError(ReleasePlanningError):
    def __init__(self, result: CheckResult):
        self.result = result
        diagnostics = PlainTextRenderer().render(result).strip()
        super().__init__(
            f'source dataset check failed with exit status {int(result.exit_status)}:\n{diagnostics}'
        )


class GitError(ReleasePlanningError):
    pass


@dataclass(frozen=True)
class ResolvedSource:
    description: str
    commit_id: str
    root: Path


@dataclass(frozen=True)
class TargetCheckout:
    root: Path
    remote_name: str
    remote_url: str


@dataclass(frozen=True)
class TagAnnotation:
    data_id: str
    analysis_id: str
    main_analysis_file: PurePosixPath
    analysis_files: tuple[PurePosixPath, ...]
    unknown_trailers: tuple[tuple[str, str], ...] = ()

    @property
    def subject(self) -> str:
        return (
            f'EOS/DATA-{self.data_id}: Supplementary material for '
            f'EOS/ANALYSIS-{self.analysis_id}'
        )

    def serialize(self) -> str:
        lines = [
            self.subject,
            '',
            'EOS-Data-Metadata-Version: 1',
            f'EOS-Analysis-ID: {self.analysis_id}',
            f'EOS-Main-Analysis-File: {self.main_analysis_file.as_posix()}',
        ]
        lines.extend(f'EOS-Analysis-File: {path.as_posix()}' for path in self.analysis_files)
        lines.extend(f'{key}: {value}' for key, value in self.unknown_trailers)
        return '\n'.join(lines)


@dataclass(frozen=True)
class FamilyTag:
    name: str
    object_id: str
    commit_id: str
    parent_commit_id: str | None


@dataclass(frozen=True)
class TagFamily:
    base_id: str
    tags: tuple[FamilyTag, ...]

    @property
    def latest(self) -> FamilyTag | None:
        return self.tags[-1] if self.tags else None

    @property
    def next_id(self) -> str:
        return self.base_id if not self.tags else f'{self.base_id}v{len(self.tags) + 1}'


@dataclass(frozen=True)
class RemoteFamilyState:
    family: TagFamily
    remote_tags: tuple[GitHubTag, ...]
    releases: frozenset[str]
    branch_commit_id: str | None


class OperationType(Enum):
    PREPARE_LOCAL_BRANCH = 'prepare-local-branch'
    CREATE_COMMIT = 'create-commit'
    MOVE_LOCAL_BRANCH = 'move-local-branch'
    CREATE_ANNOTATED_TAG = 'create-annotated-tag'
    TRANSFER_TAG = 'transfer-tag'
    CREATE_RELEASE = 'create-release'
    CREATE_REMOTE_BRANCH = 'create-remote-branch'
    ADVANCE_REMOTE_BRANCH = 'advance-remote-branch'


@dataclass(frozen=True)
class PlannedOperation:
    kind: OperationType
    ref: str
    expected_old: str | None = None
    target: str | None = None
    force: bool = False


@dataclass(frozen=True)
class CreatePlan:
    source_description: str
    source_commit_id: str
    target: TargetCheckout
    data_id: str
    base_branch: str
    parent_commit_id: str | None
    analysis_files: tuple[PurePosixPath, ...]
    main_analysis_file: PurePosixPath
    commit_message: str
    annotation: TagAnnotation
    old_tag_object_id: str | None
    old_commit_id: str | None
    check_result: CheckResult
    operations: tuple[PlannedOperation, ...]


@dataclass(frozen=True)
class PublishPlan:
    target: TargetCheckout
    data_id: str
    base_branch: str
    tag_object_id: str
    tag_commit_id: str
    annotation: TagAnnotation
    analysis_files: tuple[PurePosixPath, ...]
    main_analysis_file: PurePosixPath
    release_title: str
    release_notes: str
    branch_commit_id: str | None
    check_result: CheckResult
    operations: tuple[PlannedOperation, ...]


class GitHubReads(Protocol):
    def list_tag_refs(self) -> tuple[GitHubRef, ...]: ...
    def get_tag(self, ref: GitHubRef) -> GitHubTag: ...
    def get_release(self, tag_name: str) -> GitHubRelease | None: ...
    def get_branch(self, branch_name: str) -> str | None: ...
    def commit_parents(self, commit_id: str) -> tuple[str, ...]: ...
    def is_ancestor(self, ancestor: str, descendant: str) -> bool: ...


Runner = Callable[..., subprocess.CompletedProcess]


class LocalGitClient:
    """Read-only local Git operations, including temporary checkouts of a source."""

    def __init__(self, *, executable: str = 'git', runner: Runner = subprocess.run) -> None:
        self.executable = executable
        self._runner = runner

    def _run(
        self,
        arguments: Sequence[str],
        *,
        cwd: Path | None = None,
        env: Mapping[str, str] | None = None,
        input_text: str | None = None,
    ) -> str:
        command = [self.executable, *arguments]
        try:
            keywords = {'cwd': cwd, 'capture_output': True, 'text': True, 'check': False}
            if env is not None:
                keywords['env'] = env
            if input_text is not None:
                keywords['input'] = input_text
            result = self._runner(command, **keywords)
        except OSError as error:
            raise GitError(f'cannot execute Git: {error}') from error
        if result.returncode != 0:
            diagnostic = (result.stderr or result.stdout or '').strip()
            suffix = f': {diagnostic}' if diagnostic else ''
            raise GitError(f"Git command failed ({' '.join(arguments[:2])}){suffix}")
        return result.stdout.strip()

    def target_checkout(self, directory: Path, *, require_clean: bool = True) -> TargetCheckout:
        requested = directory.resolve()
        root = Path(self._run(['rev-parse', '--show-toplevel'], cwd=requested)).resolve()
        if root != requested:
            raise ReleasePlanningError(f'target must be the Git checkout root: {root}')
        remotes = self._run(['remote'], cwd=root).splitlines()
        matches = []
        for name in remotes:
            fetch_urls = self._run(['remote', 'get-url', '--all', name], cwd=root).splitlines()
            push_urls = self._run([
                'remote', 'get-url', '--push', '--all', name,
            ], cwd=root).splitlines()
            if (
                fetch_urls
                and push_urls
                and all(github_repository(url) == 'eos/data' for url in fetch_urls)
                and all(github_repository(url) == 'eos/data' for url in push_urls)
            ):
                matches.append((name, fetch_urls[0]))
        if not matches:
            raise ReleasePlanningError("target Git checkout does not identify GitHub repository 'eos/data'")
        if require_clean and self._run(['status', '--porcelain=v1', '--untracked-files=all'], cwd=root):
            raise ReleasePlanningError('target Git worktree has modifications or untracked files')
        name, url = sorted(matches)[0]
        return TargetCheckout(root, name, url)

    def local_tag_refs(self, target: TargetCheckout) -> Mapping[str, tuple[str, str]]:
        output = self._run([
            'for-each-ref', '--format=%(refname) %(objectname) %(*objectname)',
        ], cwd=target.root)
        result = {}
        for line in output.splitlines():
            fields = line.split()
            if not fields or not fields[0].startswith('refs/tags/'):
                continue
            fields[0] = fields[0].removeprefix('refs/tags/')
            if len(fields) == 2:
                result[fields[0]] = (fields[1], '')
                continue
            if len(fields) != 3:
                raise ReleasePlanningError(f'unexpected local tag response: {line}')
            result[fields[0]] = (fields[1], fields[2])
        return MappingProxyType(result)

    def local_branch(self, target: TargetCheckout, branch_name: str) -> str | None:
        output = self._run([
            'for-each-ref', '--format=%(objectname)', f'refs/heads/{branch_name}',
        ], cwd=target.root)
        lines = output.splitlines()
        if not lines:
            return None
        if len(lines) != 1:
            raise ReleasePlanningError(f"unexpected local branch state for '{branch_name}'")
        return lines[0]

    @contextmanager
    def checkout_source(self, source: str | Path) -> Iterator[ResolvedSource]:
        value = str(source)
        local = Path(value).expanduser()
        if local.exists():
            source_description = str(local.resolve())
            repository = source_description
            commit = self._run(['-C', repository, 'rev-parse', '--verify', 'refs/heads/main^{commit}'])
        else:
            source_description = value
            repository = value
            commit = ''
        temporary = Path(tempfile.mkdtemp(prefix='eos-data-source-'))
        try:
            self._run(['clone', '--no-checkout', '--', repository, str(temporary)])
            if not commit:
                commit = self._run(
                    ['rev-parse', '--verify', 'refs/remotes/origin/main^{commit}'],
                    cwd=temporary,
                )
            self._run(['checkout', '--detach', commit], cwd=temporary)
            actual = self._run(['rev-parse', 'HEAD^{commit}'], cwd=temporary)
            yield ResolvedSource(source_description, actual, temporary)
        finally:
            shutil.rmtree(temporary, ignore_errors=True)

    @contextmanager
    def checkout_commit(self, source: str | Path, commit_id: str) -> Iterator[Path]:
        value = str(source)
        local = Path(value).expanduser()
        repository = str(local.resolve()) if local.exists() else value
        temporary = Path(tempfile.mkdtemp(prefix='eos-data-source-'))
        try:
            self._run(['clone', '--no-checkout', '--', repository, str(temporary)])
            self._run(['checkout', '--detach', commit_id], cwd=temporary)
            actual = self._run(['rev-parse', 'HEAD^{commit}'], cwd=temporary)
            if actual != commit_id:
                raise ReleasePlanningError(
                    f'source moved while executing: expected {commit_id}, observed {actual}'
                )
            yield temporary
        finally:
            shutil.rmtree(temporary, ignore_errors=True)

    def create_dataset_commit(
        self,
        target: TargetCheckout,
        source_root: Path,
        message: str,
        parent_commit_id: str | None,
    ) -> str:
        import os

        descriptor, index_name = tempfile.mkstemp(prefix='eos-data-index-')
        os.close(descriptor)
        Path(index_name).unlink()
        environment = dict(os.environ)
        environment['GIT_INDEX_FILE'] = index_name
        git_dir = self._run(['rev-parse', '--git-dir'], cwd=target.root)
        git_dir_path = (target.root / git_dir).resolve()
        prefix = [f'--git-dir={git_dir_path}', f'--work-tree={source_root}']
        try:
            self._run([*prefix, 'read-tree', '--empty'], cwd=source_root, env=environment)
            self._run([*prefix, 'add', '--all'], cwd=source_root, env=environment)
            tree_id = self._run([*prefix, 'write-tree'], cwd=source_root, env=environment)
        finally:
            Path(index_name).unlink(missing_ok=True)
        arguments = ['commit-tree', tree_id]
        if parent_commit_id is not None:
            arguments.extend(('-p', parent_commit_id))
        arguments.extend(('-m', message))
        return self._run(arguments, cwd=target.root)

    def prepare_local_branch(
        self,
        target: TargetCheckout,
        branch_name: str,
        current_commit_id: str | None,
        parent_commit_id: str | None,
    ) -> None:
        if current_commit_id is not None:
            self._run(['checkout', '--force', branch_name], cwd=target.root)
        elif parent_commit_id is not None:
            self._run(['checkout', '--detach', '--force', parent_commit_id], cwd=target.root)
        else:
            self._run(['checkout', '--orphan', branch_name], cwd=target.root)

    def commit_parents(self, target: TargetCheckout, commit_id: str) -> tuple[str, ...]:
        return tuple(self._run(
            ['show', '-s', '--format=%P', commit_id], cwd=target.root,
        ).split())

    def update_local_branch(
        self,
        target: TargetCheckout,
        branch_name: str,
        commit_id: str,
        expected_old: str | None,
    ) -> None:
        arguments = ['update-ref', f'refs/heads/{branch_name}', commit_id]
        arguments.append(expected_old if expected_old is not None else '0' * 40)
        self._run(arguments, cwd=target.root)
        self._run(['checkout', '--force', branch_name], cwd=target.root)

    def create_local_tag(
        self,
        target: TargetCheckout,
        tag_name: str,
        commit_id: str,
        message: str,
        expected_old: str | None,
    ) -> tuple[str, str]:
        tagger = self._run(['var', 'GIT_COMMITTER_IDENT'], cwd=target.root)
        tag_contents = (
            f'object {commit_id}\n'
            'type commit\n'
            f'tag {tag_name}\n'
            f'tagger {tagger}\n\n'
            f'{message}\n'
        )
        object_id = self._run(
            ['mktag'], cwd=target.root, input_text=tag_contents,
        )
        old_object_id = expected_old if expected_old is not None else '0' * 40
        self._run([
            'update-ref', f'refs/tags/{tag_name}', object_id, old_object_id,
        ], cwd=target.root)
        result = self.local_tag_refs(target).get(tag_name)
        if result is None or result != (object_id, commit_id):
            raise ReleasePlanningError(
                f"local annotated tag '{tag_name}' failed verification: "
                f"expected {(object_id, commit_id)}, observed {result}"
            )
        actual_message = self._run(['for-each-ref', '--format=%(contents)', f'refs/tags/{tag_name}'], cwd=target.root)
        if actual_message != message:
            raise ReleasePlanningError(f"local annotated tag '{tag_name}' has unexpected metadata")
        return result


def github_repository(url: str) -> str | None:
    """Return ``owner/repository`` for supported GitHub SSH/HTTPS URLs."""
    value = url.strip()
    if value.startswith('git@github.com:'):
        path = value[len('git@github.com:'):]
    else:
        parsed = urlparse(value)
        if parsed.scheme not in ('https', 'ssh') or parsed.hostname != 'github.com':
            return None
        if parsed.password is not None:
            return None
        if parsed.scheme == 'ssh' and parsed.username not in (None, 'git'):
            return None
        if parsed.scheme == 'https' and parsed.username is not None:
            return None
        path = parsed.path.lstrip('/')
    path = path.removesuffix('.git').strip('/')
    parts = path.split('/')
    return '/'.join(parts) if len(parts) == 2 and all(parts) else None


def parse_base_id(value: str) -> str:
    if _BASE_ID.fullmatch(value) is None:
        raise ReleasePlanningError("base ID must have suffix-less 'YYYY-NN' form")
    return value


def parse_data_id(value: str) -> tuple[str, int]:
    match = _DATA_ID.fullmatch(value)
    if match is None:
        raise ReleasePlanningError("data ID must have 'YYYY-NN' or 'YYYY-NNvN' form")
    revision = int(match.group('revision') or 1)
    if match.group('revision') is not None and revision < 2:
        raise ReleasePlanningError('dataset revision v1 is invalid; the initial tag has no suffix')
    return match.group('base'), revision


def _release_path(value: str | Path) -> PurePosixPath:
    raw = str(value).replace('\\', '/')
    path = PurePosixPath(raw)
    if path.is_absolute() or not path.parts or any(part in ('', '.', '..') for part in path.parts):
        raise ReleasePlanningError(f"analysis path must be normalized and repository-relative: '{value}'")
    return path


def make_annotation(
    data_id: str,
    analysis_id: str,
    analysis_files: Sequence[str | Path],
    main_analysis_file: str | Path,
) -> TagAnnotation:
    parse_data_id(data_id)
    if _BASE_ID.fullmatch(analysis_id) is None:
        raise ReleasePlanningError('analysis ID must have YYYY-NN form')
    files = tuple(sorted((_release_path(path) for path in analysis_files), key=lambda path: path.as_posix()))
    if len(set(files)) != len(files):
        raise ReleasePlanningError('annotated-tag analysis paths must be unique')
    main = _release_path(main_analysis_file)
    if main not in files:
        raise ReleasePlanningError('main analysis path is not in the complete analysis selection')
    return TagAnnotation(data_id, analysis_id, main, files)


def parse_annotation(message: str, *, expected_data_id: str | None = None) -> TagAnnotation:
    lines = message.splitlines()
    if len(lines) < 6 or lines[1] != '':
        raise ReleasePlanningError('annotated-tag message must contain a subject and trailer block')
    subject = _SUBJECT.fullmatch(lines[0])
    if subject is None:
        raise ReleasePlanningError('annotated-tag subject is malformed')
    if expected_data_id is not None and subject.group('data_id') != expected_data_id:
        raise ReleasePlanningError('annotated-tag subject disagrees with the tag name')
    trailers: dict[str, list[str]] = {}
    unknown = []
    known = {
        'EOS-Data-Metadata-Version', 'EOS-Analysis-ID',
        'EOS-Main-Analysis-File', 'EOS-Analysis-File',
    }
    for line in lines[2:]:
        match = _TRAILER.fullmatch(line)
        if match is None:
            raise ReleasePlanningError(f'malformed annotated-tag trailer: {line!r}')
        key, value = match.group('key'), match.group('value')
        trailers.setdefault(key, []).append(value)
        if key not in known:
            unknown.append((key, value))
    for key in ('EOS-Data-Metadata-Version', 'EOS-Analysis-ID', 'EOS-Main-Analysis-File'):
        if len(trailers.get(key, ())) != 1:
            raise ReleasePlanningError(f'annotated-tag trailer {key} must occur exactly once')
    if trailers['EOS-Data-Metadata-Version'][0] != '1':
        raise ReleasePlanningError('unsupported EOS data metadata version')
    if trailers['EOS-Analysis-ID'][0] != subject.group('analysis_id'):
        raise ReleasePlanningError('analysis ID trailer disagrees with the annotated-tag subject')
    files = trailers.get('EOS-Analysis-File', [])
    annotation = make_annotation(
        subject.group('data_id'), trailers['EOS-Analysis-ID'][0], files,
        trailers['EOS-Main-Analysis-File'][0],
    )
    return TagAnnotation(
        annotation.data_id, annotation.analysis_id, annotation.main_analysis_file,
        annotation.analysis_files, tuple(unknown),
    )


def validate_annotation_documents(annotation: TagAnnotation, root: Path) -> None:
    for path in annotation.analysis_files:
        candidate = root / Path(path.as_posix())
        try:
            document = yaml.safe_load(candidate.read_text(encoding='utf-8'))
            actual = document['metadata']['id']
        except (OSError, UnicodeError, yaml.YAMLError, TypeError, KeyError) as error:
            raise ReleasePlanningError(f"cannot validate tagged analysis '{path}': {error}") from error
        if actual != annotation.analysis_id:
            raise ReleasePlanningError(
                f"tag annotation analysis ID disagrees with tagged file '{path}'"
            )


def validate_family(
    base_id: str,
    local_tags: Mapping[str, tuple[str, str]],
    remote_tags: Sequence[GitHubTag],
    parent_lookup: Callable[[str], Sequence[str]],
) -> TagFamily:
    parse_base_id(base_id)
    combined: dict[str, tuple[str, str]] = {}
    prefix = re.compile(rf'^{re.escape(base_id)}(?:v.*)?$')
    for name, target in local_tags.items():
        if prefix.fullmatch(name):
            try:
                tag_base, _ = parse_data_id(name)
            except ReleasePlanningError as error:
                raise ReleasePlanningError(f"malformed member of tag family '{base_id}': {name}") from error
            if tag_base == base_id:
                if not target[1]:
                    raise ReleasePlanningError(
                        f"local tag '{name}' is lightweight or has an unexpected ref type"
                    )
                combined[name] = target
    for tag in remote_tags:
        if not prefix.fullmatch(tag.name):
            continue
        try:
            tag_base, _ = parse_data_id(tag.name)
        except ReleasePlanningError as error:
            raise ReleasePlanningError(f"malformed member of tag family '{base_id}': {tag.name}") from error
        if tag_base != base_id:
            continue
        target = (tag.object_id, tag.commit_id)
        if tag.name in combined and combined[tag.name] != target:
            raise ReleasePlanningError(f"local and remote tag '{tag.name}' have conflicting targets")
        combined[tag.name] = target
    revisions = {}
    for name in combined:
        _, revision = parse_data_id(name)
        revisions[revision] = name
    if revisions and sorted(revisions) != list(range(1, max(revisions) + 1)):
        raise ReleasePlanningError(f"tag family '{base_id}' has a revision gap")
    tags = []
    previous = None
    for revision in sorted(revisions):
        name = revisions[revision]
        object_id, commit_id = combined[name]
        parents = tuple(parent_lookup(commit_id))
        if revision == 1:
            if parents:
                raise ReleasePlanningError(f"initial tag '{name}' is not an orphan commit")
            parent = None
        else:
            if parents != (previous,):
                raise ReleasePlanningError(f"tag '{name}' does not directly extend the validated tag lineage")
            parent = previous
        tags.append(FamilyTag(name, object_id, commit_id, parent))
        previous = commit_id
    return TagFamily(base_id, tuple(tags))


def _remote_tags(github: GitHubReads, base_id: str) -> tuple[GitHubTag, ...]:
    result = []
    for ref in github.list_tag_refs():
        if ref.name == base_id or ref.name.startswith(f'{base_id}v'):
            result.append(github.get_tag(ref))
    return tuple(result)


def read_remote_family(
    github: GitHubReads,
    base_id: str,
    local_tags: Mapping[str, tuple[str, str]],
) -> RemoteFamilyState:
    tags = _remote_tags(github, base_id)
    family = validate_family(base_id, local_tags, tags, github.commit_parents)
    remote_names = {tag.name for tag in tags}
    local_family_names = {
        tag.name for tag in family.tags if tag.name in local_tags
    }
    missing_remote = sorted(local_family_names - remote_names)
    if missing_remote:
        raise ReleasePlanningError(
            'local tag family contains tag(s) absent from GitHub: ' + ', '.join(missing_remote)
        )
    releases = frozenset(tag.name for tag in family.tags if github.get_release(tag.name) is not None)
    branch = github.get_branch(base_id)
    return RemoteFamilyState(family, tags, releases, branch)


def _run_source_checks(
    root: Path,
    analysis_files: Sequence[str | Path],
    main_analysis_file: str | Path | None,
    factory,
) -> tuple[CheckResult, CheckContext]:
    context = CheckContext(
        root,
        tuple(Path(path) for path in analysis_files),
        Path(main_analysis_file) if main_analysis_file is not None else None,
    )
    result = run_checks(factory, context, scope=CheckScope.COMPLETE)
    if result.exit_status is not ExitStatus.SUCCESS:
        raise SourceCheckError(result)
    return result, context


def _exact_zenodo_data_id(context: CheckContext, expected: str) -> None:
    document = context.cache.get(CacheKey.ZENODO_METADATA)
    title = document.get('title') if isinstance(document, Mapping) else None
    prefix = f'EOS/DATA-{expected}: '
    if not isinstance(title, str) or not title.startswith(prefix):
        raise ReleasePlanningError(f'Zenodo dataset ID must equal the derived tag ID {expected}')


def plan_create(
    base_id: str,
    source: str | Path,
    target_directory: Path,
    *,
    analysis_files: Sequence[str | Path] = (),
    main_analysis_file: str | Path | None = None,
    revision: bool = False,
    replace: bool = False,
    git: LocalGitClient,
    github: GitHubReads,
    check_factory,
) -> CreatePlan:
    base_id = parse_base_id(base_id)
    if revision and replace:
        raise ReleasePlanningError('--revision and --replace are mutually exclusive')
    target = git.target_checkout(target_directory, require_clean=True)

    with git.checkout_source(source) as resolved:
        result, context = _run_source_checks(
            resolved.root, analysis_files, main_analysis_file, check_factory,
        )
        analysis_id = context.cache[CacheKey.COMMON_ANALYSIS_ID]
        selected = tuple(
            PurePosixPath(path.relative_to(context.dataset_root).as_posix())
            for path in context.analysis_paths
        )
        main = PurePosixPath(context.main_analysis_path.relative_to(context.dataset_root).as_posix())
        source_description = resolved.description
        source_commit_id = resolved.commit_id

    local_tags = git.local_tag_refs(target)
    local_branch = git.local_branch(target, base_id)
    remote = read_remote_family(github, base_id, local_tags)
    family = remote.family
    if not family.tags:
        if revision or replace:
            raise ReleasePlanningError('initial creation does not accept --revision or --replace')
        data_id = base_id
        mode = 'initial'
    else:
        if revision == replace:
            raise ReleasePlanningError('an existing base tag requires exactly one of --revision or --replace')
        if replace:
            if len(family.tags) != 1:
                raise ReleasePlanningError('replacement is forbidden when revision tags exist')
            if base_id in remote.releases:
                raise ReleasePlanningError('replacement is forbidden after GitHub release publication')
            if remote.branch_commit_id is not None:
                raise ReleasePlanningError('replacement is forbidden after the remote base branch exists')
            data_id = base_id
            mode = 'replacement'
        else:
            data_id = family.next_id
            mode = 'revision'
    if mode == 'initial':
        if local_branch is not None:
            raise ReleasePlanningError(f"local base branch '{base_id}' already exists without a tag family")
        if remote.branch_commit_id is not None:
            raise ReleasePlanningError(f"remote base branch '{base_id}' already exists without a tag family")
    elif mode == 'revision' and local_branch not in (None, family.latest.commit_id):
        raise ReleasePlanningError(f"local base branch '{base_id}' is not at the latest family commit")
    elif mode == 'replacement' and local_branch not in (None, family.latest.commit_id):
        raise ReleasePlanningError(f"local base branch '{base_id}' is inconsistent with the base tag")
    if mode == 'revision' and remote.branch_commit_id is not None:
        family_commits = {tag.commit_id for tag in family.tags}
        if remote.branch_commit_id not in family_commits:
            raise ReleasePlanningError('remote base branch is outside the validated tag lineage')
    if data_id in remote.releases:
        raise ReleasePlanningError(f"GitHub release already exists for '{data_id}'")

    _exact_zenodo_data_id(context, data_id)
    annotation = make_annotation(data_id, analysis_id, selected, main)

    parent = family.latest.commit_id if mode == 'revision' else None
    old = family.latest if mode == 'replacement' else None
    message = annotation.subject
    operations = [
        PlannedOperation(OperationType.PREPARE_LOCAL_BRANCH, base_id, old.commit_id if old else None, parent),
        PlannedOperation(OperationType.CREATE_COMMIT, base_id, parent, source_commit_id),
        PlannedOperation(
            OperationType.MOVE_LOCAL_BRANCH,
            base_id,
            local_branch,
            '<new-commit>',
            mode == 'replacement',
        ),
    ]
    operations.extend((
        PlannedOperation(
            OperationType.CREATE_ANNOTATED_TAG,
            data_id,
            old.object_id if old else None,
            '<new-tag-object>',
            mode == 'replacement',
        ),
        PlannedOperation(
            OperationType.TRANSFER_TAG,
            f'{target.remote_name}:{data_id}',
            old.object_id if old else None,
            '<new-tag-object>',
            mode == 'replacement',
        ),
    ))
    return CreatePlan(
        source_description, source_commit_id, target, data_id, base_id, parent,
        annotation.analysis_files, annotation.main_analysis_file, message, annotation,
        old.object_id if old else None, old.commit_id if old else None,
        result, tuple(operations),
    )


def plan_publish(
    data_id: str,
    target_directory: Path,
    *,
    git: LocalGitClient,
    github: GitHubReads,
    check_factory,
) -> PublishPlan:
    base_id, _ = parse_data_id(data_id)
    target = git.target_checkout(target_directory, require_clean=False)
    local_tags = git.local_tag_refs(target)
    remote = read_remote_family(github, base_id, local_tags)
    remote_tag = next((tag for tag in remote.remote_tags if tag.name == data_id), None)
    if remote_tag is None:
        raise ReleasePlanningError(f"remote annotated tag '{data_id}' does not exist")
    if data_id in remote.releases:
        raise ReleasePlanningError(f"GitHub release already exists for '{data_id}'")
    annotation = parse_annotation(remote_tag.message, expected_data_id=data_id)

    checkout = getattr(github, 'checkout_tree', None)
    if not callable(checkout):
        raise ReleasePlanningError('GitHub read client cannot check out the tagged commit tree')
    with checkout(remote_tag.commit_id) as root:
        validate_annotation_documents(annotation, root)
        result, context = _run_source_checks(
            root, annotation.analysis_files, annotation.main_analysis_file, check_factory,
        )
        if context.cache[CacheKey.COMMON_ANALYSIS_ID] != annotation.analysis_id:
            raise ReleasePlanningError('tag annotation analysis ID disagrees with tagged analyses')
        metadata = context.cache[CacheKey.ANALYSIS_METADATA][context.main_analysis_path]
        title = metadata['title']
        notes = metadata['description']

    branch = remote.branch_commit_id
    if branch is not None and not github.is_ancestor(branch, remote_tag.commit_id):
        raise ReleasePlanningError('remote base branch is not an ancestor of the tag commit')
    operations = [PlannedOperation(
        OperationType.CREATE_RELEASE, data_id, None, remote_tag.commit_id,
    )]
    if branch is None:
        operations.append(PlannedOperation(
            OperationType.CREATE_REMOTE_BRANCH, base_id, None, remote_tag.commit_id,
        ))
    elif branch != remote_tag.commit_id:
        operations.append(PlannedOperation(
            OperationType.ADVANCE_REMOTE_BRANCH, base_id, branch, remote_tag.commit_id,
        ))
    return PublishPlan(
        target, data_id, base_id, remote_tag.object_id, remote_tag.commit_id, annotation,
        annotation.analysis_files, annotation.main_analysis_file, title, notes, branch,
        result, tuple(operations),
    )


def render_create_plan(plan: CreatePlan) -> str:
    lines = [
        'Create dataset/tag plan (read-only)',
        f'Source repository: {plan.source_description}',
        f'Source commit: {plan.source_commit_id}',
        f'Target repository: eos/data ({plan.target.remote_name}: {plan.target.remote_url})',
        f'Derived tag: {plan.data_id}',
        f'Local base branch: {plan.base_branch}',
        f'Parent/ancestry: {plan.parent_commit_id or "orphan"}',
        'Analysis files: ' + ', '.join(path.as_posix() for path in plan.analysis_files),
        f'Main analysis file: {plan.main_analysis_file}',
        f'Commit message: {plan.commit_message}',
        f'Old tag object: {plan.old_tag_object_id or "none"}',
        f'Old commit: {plan.old_commit_id or "none"}',
        'Check findings:',
        PlainTextRenderer().render(plan.check_result).rstrip(),
        'Annotated tag message:',
        plan.annotation.serialize(),
        'Operations:',
    ]
    for index, operation in enumerate(plan.operations, 1):
        lines.append(
            f'{index}. {operation.kind.value}: {operation.ref}; '
            f'old={operation.expected_old or "none"}; target={operation.target or "execution-derived"}; '
            f'force={str(operation.force).lower()}'
        )
    return '\n'.join(lines) + '\n'


def render_publish_plan(plan: PublishPlan) -> str:
    lines = [
        'Publish release plan (read-only)',
        f'Target repository: eos/data ({plan.target.remote_name}: {plan.target.remote_url})',
        f'Tag: {plan.data_id}',
        f'Tag object: {plan.tag_object_id}',
        f'Tag commit: {plan.tag_commit_id}',
        f'Base branch: {plan.base_branch}',
        f'Current branch commit: {plan.branch_commit_id or "absent"}',
        'Analysis files: ' + ', '.join(path.as_posix() for path in plan.analysis_files),
        f'Main analysis file: {plan.main_analysis_file}',
        f'Release title: {plan.release_title}',
        'Release notes:',
        plan.release_notes,
        'Check findings:',
        PlainTextRenderer().render(plan.check_result).rstrip(),
        'Operations:',
    ]
    for index, operation in enumerate(plan.operations, 1):
        lines.append(
            f'{index}. {operation.kind.value}: {operation.ref}; '
            f'old={operation.expected_old or "none"}; target={operation.target or "none"}'
        )
    return '\n'.join(lines) + '\n'
