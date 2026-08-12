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

"""Planning and execution for ``eos-data register``."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import difflib
from pathlib import Path, PurePosixPath
import re
from typing import Any, Protocol

import yaml

from ..datasets_file_description import DataSetDescription, PublicLikelihoodDescription
from .data_checks import CacheKey
from .data_github import GitHubRef, GitHubRelease, GitHubTag
from .data_release import (
    LocalGitClient, ReleasePlanningError, TagAnnotation, parse_annotation,
    parse_data_id, validate_annotation_documents, _run_source_checks,
)


_DOI = re.compile(r'^10\.5281/zenodo\.[0-9]+$')
_REPOSITORY = re.compile(r'^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$')


class _RegistryDumper(yaml.SafeDumper):
    def increase_indent(self, flow: bool = False, indentless: bool = False):
        return super().increase_indent(flow, False)


class RegistrationError(RuntimeError):
    pass


@dataclass(frozen=True)
class RegistryEntry:
    data_id: str
    authors: tuple[str, ...]
    title: str
    doi: str
    likelihoods: tuple[tuple[str, str, str], ...]
    keywords: tuple[str, ...]
    eos_version: str

    def as_mapping(self) -> dict[str, Any]:
        return {
            'authors': list(self.authors),
            'title': self.title,
            'doi': self.doi,
            'likelihoods': {
                name: {'filename': filename, 'filetype': filetype}
                for name, filename, filetype in self.likelihoods
            },
            'keywords': list(self.keywords),
            'eos_version': self.eos_version,
        }


@dataclass(frozen=True)
class RegisterPlan:
    repository: str
    data_id: str
    tag_object_id: str
    tag_commit_id: str
    annotation: TagAnnotation
    release_url: str
    main_commit_id: str
    branch_name: str
    commit_message: str
    entry: RegistryEntry
    original_registry: str
    updated_registry: str
    diff: str
    operations: tuple[str, ...]


@dataclass(frozen=True)
class RegisterResult:
    data_id: str
    branch_name: str
    commit_id: str
    completed: tuple[str, ...]


class RegisterGitHub(Protocol):
    def list_tag_refs(self) -> tuple[GitHubRef, ...]: ...
    def get_tag(self, ref: GitHubRef) -> GitHubTag: ...
    def get_release(self, tag_name: str) -> GitHubRelease | None: ...
    def get_branch(self, branch_name: str) -> str | None: ...
    def checkout_tree(self, commit_id: str): ...
    def transfer_local_branch(
        self, root: Path, branch_name: str, commit_id: str, main_commit_id: str,
        tag_name: str, tag_object_id: str,
    ) -> None: ...


def parse_repository(value: str) -> str:
    if _REPOSITORY.fullmatch(value) is None:
        raise RegistrationError("repository must have 'OWNER/REPO' form")
    return value


def parse_doi(value: str) -> str:
    if _DOI.fullmatch(value) is None:
        raise RegistrationError("DOI must have '10.5281/zenodo.<record>' form")
    return value


def parse_keywords(values: Sequence[str]) -> tuple[str, ...]:
    keywords = tuple(value.strip() for value in values)
    if not keywords or any(not value for value in keywords):
        raise RegistrationError('at least one nonempty keyword is required')
    if len(set(keywords)) != len(keywords):
        raise RegistrationError('keywords must not contain duplicates')
    return keywords


def parse_eos_version(value: str) -> str:
    result = value.strip()
    if not result or any(character.isspace() for character in result):
        raise RegistrationError('EOS version must be a nonempty value without whitespace')
    return result


def parse_likelihoods(values: Sequence[str]) -> tuple[tuple[str, str, str], ...]:
    result = []
    names = set()
    for value in values:
        fields = value.split(':')
        if len(fields) != 3 or any(not field.strip() for field in fields):
            raise RegistrationError(
                "likelihood must have nonempty 'NAME:FILENAME:FILETYPE' fields"
            )
        name, filename, filetype = (field.strip() for field in fields)
        if name in names:
            raise RegistrationError(f"duplicate likelihood name '{name}'")
        try:
            PublicLikelihoodDescription(filename, filetype)
        except ValueError as error:
            raise RegistrationError(str(error)) from error
        names.add(name)
        result.append((name, filename, filetype))
    return tuple(result)


def _validate_registry(text: str) -> Mapping[str, Any]:
    try:
        document = yaml.safe_load(text)
    except yaml.YAMLError as error:
        raise RegistrationError(f'cannot parse remote datasets.yaml: {error}') from error
    if not isinstance(document, Mapping) or not all(isinstance(key, str) for key in document):
        raise RegistrationError('remote datasets.yaml must contain a string-keyed mapping')
    for data_id, raw in document.items():
        try:
            parse_data_id(data_id)
        except ReleasePlanningError as error:
            raise RegistrationError(f"registry entry ID '{data_id}' is malformed") from error
        if not isinstance(raw, Mapping):
            raise RegistrationError(f"registry entry '{data_id}' must be a mapping")
        required = {'authors', 'title', 'keywords', 'likelihoods'}
        missing = required - raw.keys()
        if missing:
            raise RegistrationError(
                f"registry entry '{data_id}' is missing {', '.join(sorted(missing))}"
            )
        authors = raw['authors']
        keywords = raw['keywords']
        likelihoods = raw['likelihoods']
        if (
            not isinstance(authors, list) or not authors
            or not all(isinstance(author, str) and author.strip() for author in authors)
        ):
            raise RegistrationError(f"registry entry '{data_id}' has invalid authors")
        if not isinstance(raw['title'], str) or not raw['title'].strip():
            raise RegistrationError(f"registry entry '{data_id}' has invalid title")
        if not isinstance(keywords, list) or not all(
            isinstance(keyword, str) and keyword.strip() for keyword in keywords
        ):
            raise RegistrationError(f"registry entry '{data_id}' has invalid keywords")
        if not isinstance(likelihoods, Mapping) or not all(
            isinstance(name, str) and name and isinstance(value, Mapping)
            for name, value in likelihoods.items()
        ):
            raise RegistrationError(f"registry entry '{data_id}' has invalid likelihoods")
        try:
            for value in likelihoods.values():
                PublicLikelihoodDescription.from_dict(**value)
        except (TypeError, ValueError) as error:
            raise RegistrationError(f"invalid registry entry '{data_id}': {error}") from error
        if 'doi' in raw:
            try:
                parse_doi(raw['doi'])
            except (TypeError, RegistrationError) as error:
                raise RegistrationError(f"invalid registry entry '{data_id}': {error}") from error
        if 'eos_version' in raw:
            try:
                parse_eos_version(raw['eos_version'])
            except (AttributeError, RegistrationError) as error:
                raise RegistrationError(f"invalid registry entry '{data_id}': {error}") from error
        known = {key: raw[key] for key in DataSetDescription.__dataclass_fields__ if key in raw}
        if 'eos_version' in known:
            try:
                DataSetDescription.from_dict(**known)
            except (TypeError, ValueError) as error:
                raise RegistrationError(f"invalid registry entry '{data_id}': {error}") from error
    return document


def _updated_registry(text: str, entry: RegistryEntry) -> tuple[str, str]:
    document = _validate_registry(text)
    if entry.data_id in document:
        adjective = 'identical' if document[entry.data_id] == entry.as_mapping() else 'conflicting'
        raise RegistrationError(f"registry already contains {adjective} ID '{entry.data_id}'")
    rendered = yaml.dump(
        {entry.data_id: entry.as_mapping()}, Dumper=_RegistryDumper,
        sort_keys=False, allow_unicode=True, width=4096,
    )
    updated = rendered + ('\n' if text and not text.startswith('\n') else '') + text
    difference = ''.join(difflib.unified_diff(
        text.splitlines(keepends=True), updated.splitlines(keepends=True),
        fromfile='main/datasets.yaml', tofile=f'register-{entry.data_id}/datasets.yaml',
    ))
    return updated, difference


def _exact_tag(github: RegisterGitHub, data_id: str) -> GitHubTag:
    refs = [ref for ref in github.list_tag_refs() if ref.name == data_id]
    if len(refs) != 1:
        raise RegistrationError(f"exact remote tag '{data_id}' does not exist")
    if refs[0].object_type != 'tag':
        raise RegistrationError(f"remote tag '{data_id}' is lightweight")
    tag = github.get_tag(refs[0])
    if (
        tag.name != data_id or tag.object_id != refs[0].object_id
        or not tag.commit_id or not isinstance(tag.message, str)
    ):
        raise RegistrationError(f"remote annotated tag '{data_id}' has inconsistent metadata")
    return tag


def plan_register(
    data_id: str,
    *,
    doi: str,
    keywords: Sequence[str],
    eos_version: str,
    likelihoods: Sequence[str] = (),
    repository: str = 'eos/data',
    github: RegisterGitHub,
    check_factory,
) -> RegisterPlan:
    parse_data_id(data_id)
    repository = parse_repository(repository)
    doi = parse_doi(doi)
    keywords = parse_keywords(keywords)
    eos_version = parse_eos_version(eos_version)
    likelihoods = parse_likelihoods(likelihoods)
    tag = _exact_tag(github, data_id)
    annotation = parse_annotation(tag.message, expected_data_id=data_id)
    release = github.get_release(data_id)
    if release is None or release.tag_name != data_id:
        raise RegistrationError(f"GitHub release for exact tag '{data_id}' does not exist")

    with github.checkout_tree(tag.commit_id) as root:
        validate_annotation_documents(annotation, root)
        _result, context = _run_source_checks(
            root, annotation.analysis_files, annotation.main_analysis_file, check_factory,
        )
        metadata = context.cache[CacheKey.ANALYSIS_METADATA][context.main_analysis_path]
        authors = context.cache[CacheKey.ANALYSIS_AUTHOR_NAMES][context.main_analysis_path]
        title = metadata['title']

    main_commit = github.get_branch('main')
    if main_commit is None:
        raise RegistrationError("remote branch 'main' does not exist")
    branch_name = f'register-{data_id}'
    if github.get_branch(branch_name) is not None:
        raise RegistrationError(f"remote branch '{branch_name}' already exists")
    with github.checkout_tree(main_commit) as root:
        registry_path = root / 'datasets.yaml'
        try:
            original = registry_path.read_text(encoding='utf-8')
        except (OSError, UnicodeError) as error:
            raise RegistrationError(f'cannot read remote main datasets.yaml: {error}') from error
    entry = RegistryEntry(
        data_id, tuple(authors), title, doi, tuple(likelihoods), tuple(keywords), eos_version,
    )
    updated, difference = _updated_registry(original, entry)
    message = f'Add {data_id}'
    return RegisterPlan(
        repository, data_id, tag.object_id, tag.commit_id, annotation, release.url,
        main_commit, branch_name, message, entry, original, updated, difference,
        ('checkout-remote-main', 'create-local-branch', 'update-datasets.yaml',
         'create-single-commit', 'push-new-branch', 'verify-remote-branch'),
    )


def render_register_plan(plan: RegisterPlan) -> str:
    entry = yaml.dump(
        {plan.data_id: plan.entry.as_mapping()}, Dumper=_RegistryDumper,
        sort_keys=False, allow_unicode=True, width=4096,
    ).rstrip()
    lines = [
        'Register dataset plan (read-only)',
        f'Repository: {plan.repository}',
        f'Tag: {plan.data_id}',
        f'Tag object: {plan.tag_object_id}',
        f'Tag commit: {plan.tag_commit_id}',
        f'GitHub release: {plan.release_url}',
        f'Source main commit: {plan.main_commit_id}',
        f'Branch: {plan.branch_name}',
        f'Commit message: {plan.commit_message}',
        'Prospective registry entry:', entry,
        'YAML diff:', plan.diff.rstrip(),
        'Operations:',
    ]
    lines.extend(f'{index}. {operation}' for index, operation in enumerate(plan.operations, 1))
    return '\n'.join(lines) + '\n'


def _signature(plan: RegisterPlan) -> tuple[Any, ...]:
    return (
        plan.repository, plan.data_id, plan.tag_object_id, plan.tag_commit_id,
        plan.annotation, plan.release_url, plan.main_commit_id, plan.branch_name,
        plan.commit_message, plan.entry, plan.original_registry, plan.updated_registry,
        plan.diff, plan.operations,
    )


def execute_register(
    plan: RegisterPlan,
    *,
    git: LocalGitClient,
    github: RegisterGitHub,
    check_factory,
) -> RegisterResult:
    fresh = plan_register(
        plan.data_id, doi=plan.entry.doi, keywords=plan.entry.keywords,
        eos_version=plan.entry.eos_version,
        likelihoods=tuple(':'.join(value) for value in plan.entry.likelihoods),
        repository=plan.repository, github=github, check_factory=check_factory,
    )
    if _signature(fresh) != _signature(plan):
        raise RegistrationError('registration preconditions changed; no writes performed')
    completed = []
    try:
        with github.checkout_tree(plan.main_commit_id) as root:
            completed.append(f'checkout-main:{plan.main_commit_id}')
            if github.get_branch('main') != plan.main_commit_id:
                raise RegistrationError('remote main advanced before registration commit creation')
            commit_id = git.create_registration_commit(
                root, plan.branch_name, plan.main_commit_id,
                plan.updated_registry, plan.commit_message,
            )
            completed.append(f'create-commit:{commit_id}')
            observed_tag = _exact_tag(github, plan.data_id)
            observed_release = github.get_release(plan.data_id)
            if (
                observed_tag.object_id != plan.tag_object_id
                or observed_tag.commit_id != plan.tag_commit_id
                or observed_tag.message.rstrip('\n') != plan.annotation.serialize()
                or observed_release is None
                or observed_release.tag_name != plan.data_id
                or observed_release.url != plan.release_url
            ):
                raise RegistrationError('published tag or release changed before branch transfer')
            if github.get_branch('main') != plan.main_commit_id:
                raise RegistrationError('remote main advanced before registration branch transfer')
            if github.get_branch(plan.branch_name) is not None:
                raise RegistrationError(f"remote branch '{plan.branch_name}' appeared before transfer")
            github.transfer_local_branch(
                root, plan.branch_name, commit_id, plan.main_commit_id,
                plan.data_id, plan.tag_object_id,
            )
            completed.append(f'push-branch:{plan.branch_name}:{commit_id}')
        observed_main = github.get_branch('main')
        observed_tag = _exact_tag(github, plan.data_id)
        observed_release = github.get_release(plan.data_id)
        observed_branch = github.get_branch(plan.branch_name)
        if (
            observed_main != plan.main_commit_id
            or observed_tag.object_id != plan.tag_object_id
            or observed_tag.commit_id != plan.tag_commit_id
            or observed_tag.message.rstrip('\n') != plan.annotation.serialize()
            or observed_release is None
            or observed_release.tag_name != plan.data_id
            or observed_release.url != plan.release_url
            or observed_branch != commit_id
        ):
            raise RegistrationError(
                'post-push publication verification failed: '
                f'main={observed_main or "absent"}, tag={observed_tag.object_id}, '
                f'release={observed_release.url if observed_release else "absent"}, '
                f'branch={observed_branch or "absent"}'
            )
        completed.append(f'verify-branch:{plan.branch_name}:{commit_id}')
        return RegisterResult(plan.data_id, plan.branch_name, commit_id, tuple(completed))
    except Exception as error:
        state = ', '.join(completed) if completed else 'none'
        raise RegistrationError(
            f'registration execution failed: {error}. Completed operations: {state}. '
            f"Safe recovery: inspect remote branch '{plan.branch_name}' and rerun only after "
            'reconciling its target; do not delete published release or tag state.'
        ) from error
