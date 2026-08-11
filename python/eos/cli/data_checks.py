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

from __future__ import annotations

from collections.abc import Iterable, Mapping, MutableMapping
from dataclasses import dataclass, field
from enum import Enum, IntEnum
import json
from pathlib import Path
import re
from types import MappingProxyType
from typing import Any, Protocol, TextIO
import unicodedata

import yaml

from ..analysis_file_description import MetadataAuthorDescription
from ..deserializable import load_yaml_file
from ..diagnostic import Severity


_ANALYSIS_ID = re.compile(r'^[0-9]{4}-[0-9]{2}$')
_DATASET_ID = re.compile(r'^[0-9]{4}-[0-9]{2}(?:v([2-9]|[1-9][0-9]+))?$')
_ZENODO_TITLE = re.compile(
    r'^EOS/DATA-(?P<dataset_id>[^:]+): Supplementary material for '
    r'EOS/ANALYSIS-(?P<analysis_id>[^\s]+)$'
)
_WHITESPACE = re.compile(r'\s+')


class ExitStatus(IntEnum):
    SUCCESS = 0
    CHECK_ERRORS = 1
    USAGE_OR_INFRASTRUCTURE = 2


class CheckScope(IntEnum):
    BASIC = 1


class CacheKey(Enum):
    ANALYSIS_SELECTION_VALID = 'analysis-selection-valid'
    ANALYSIS_METADATA_VALID = 'analysis-metadata-valid'
    ANALYSIS_METADATA = 'analysis-metadata'
    ANALYSIS_AUTHOR_NAMES = 'analysis-author-names'
    ANALYSIS_CONSISTENCY_VALID = 'analysis-consistency-valid'
    COMMON_ANALYSIS_ID = 'common-analysis-id'
    COMMON_AUTHORS_NORMALIZED = 'common-authors-normalized'
    ZENODO_CREATORS_VALID = 'zenodo-creators-valid'
    ZENODO_CREATORS_NORMALIZED = 'zenodo-creators-normalized'


def _immutable_mapping(values: Mapping[str, Any] | None) -> Mapping[str, Any]:
    return MappingProxyType(dict(values or {}))


@dataclass(frozen=True)
class InteractiveQuestion:
    prompt: str
    details: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, 'details', _immutable_mapping(self.details))


@dataclass(frozen=True)
class Finding:
    check_id: str
    severity: Severity
    message: str
    path: Path | None = None
    details: Mapping[str, Any] = field(default_factory=dict)
    question: InteractiveQuestion | None = None

    def __post_init__(self) -> None:
        if not self.check_id:
            raise ValueError('a finding requires a nonempty check identifier')
        if not self.message:
            raise ValueError('a finding requires a nonempty message')
        if self.path is not None:
            object.__setattr__(self, 'path', Path(self.path))
        object.__setattr__(self, 'details', _immutable_mapping(self.details))


@dataclass
class CheckContext:
    dataset_root: Path
    analysis_paths: tuple[Path, ...] = ()
    main_analysis_path: Path | None = None
    cache: MutableMapping[CacheKey, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.dataset_root = Path(self.dataset_root)
        self.analysis_paths = tuple(Path(path) for path in self.analysis_paths)
        if self.main_analysis_path is not None:
            self.main_analysis_path = Path(self.main_analysis_path)


class CheckFunction(Protocol):
    def __call__(self, context: CheckContext) -> Iterable[Finding]:
        ...


@dataclass(frozen=True)
class CheckDeclaration:
    name: str
    scope: CheckScope
    needs: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError('a check declaration requires a nonempty name')
        if not isinstance(self.scope, CheckScope):
            raise TypeError('a check declaration requires a valid scope')
        object.__setattr__(self, 'needs', frozenset(self.needs))


@dataclass(frozen=True)
class RegisteredCheck:
    declaration: CheckDeclaration
    function: CheckFunction


class CheckFactory:
    """An ordered, closed registry of dataset checks."""

    def __init__(self) -> None:
        self._declarations: dict[str, CheckDeclaration] = {}
        self._checks: dict[str, RegisteredCheck] = {}

    def declare(
        self, name: str, *, scope: CheckScope, needs: Iterable[str] = (),
    ) -> CheckDeclaration:
        if name in self._declarations:
            raise ValueError(f"check '{name}' is already declared")

        declaration = CheckDeclaration(name=name, scope=scope, needs=frozenset(needs))
        self._declarations[name] = declaration
        return declaration

    def register(self, name: str, function: CheckFunction) -> None:
        if name not in self._declarations:
            raise KeyError(f"unknown check registration: '{name}'")
        if name in self._checks:
            raise ValueError(f"check '{name}' is already registered")
        if not callable(function):
            raise TypeError(f"check '{name}' is not callable")

        self._checks[name] = RegisteredCheck(self._declarations[name], function)

    def checks(self, *, scope: CheckScope) -> tuple[RegisteredCheck, ...]:
        enabled = {
            name: declaration
            for name, declaration in self._declarations.items()
            if declaration.scope <= scope
        }
        for declaration in enabled.values():
            unknown = declaration.needs - self._declarations.keys()
            if unknown:
                names = ', '.join(sorted(unknown))
                raise ValueError(f"check '{declaration.name}' has unknown prerequisite(s): {names}")
            unavailable = declaration.needs - enabled.keys()
            if unavailable:
                names = ', '.join(sorted(unavailable))
                raise ValueError(
                    f"check '{declaration.name}' depends on check(s) outside the requested "
                    f'scope: {names}'
                )

        ordered: list[RegisteredCheck] = []
        remaining = dict(enabled)
        completed: set[str] = set()
        while remaining:
            ready = [
                name
                for name, declaration in remaining.items()
                if declaration.needs <= completed
            ]
            if not ready:
                names = ', '.join(remaining)
                raise ValueError(f'cyclic check prerequisites: {names}')
            for name in ready:
                if name in self._checks:
                    ordered.append(self._checks[name])
                completed.add(name)
                del remaining[name]

        return tuple(ordered)

    def missing(self, *, scope: CheckScope) -> tuple[CheckDeclaration, ...]:
        return tuple(
            declaration
            for name, declaration in self._declarations.items()
            if declaration.scope <= scope and name not in self._checks
        )

    def __len__(self) -> int:
        return len(self._checks)


@dataclass(frozen=True)
class CheckResult:
    findings: tuple[Finding, ...] = ()
    failure: str | None = None

    @property
    def completed(self) -> bool:
        return self.failure is None

    @property
    def exit_status(self) -> ExitStatus:
        if not self.completed:
            return ExitStatus.USAGE_OR_INFRASTRUCTURE
        if any(finding.severity is Severity.ERROR for finding in self.findings):
            return ExitStatus.CHECK_ERRORS
        return ExitStatus.SUCCESS


class CheckRunner:
    def __init__(self, factory: CheckFactory, *, scope: CheckScope) -> None:
        self._checks = factory.checks(scope=scope)
        self._missing = factory.missing(scope=scope)

    def run(self, context: CheckContext) -> CheckResult:
        if self._missing:
            names = ', '.join(declaration.name for declaration in self._missing)
            return CheckResult(failure=f'unregistered built-in checks: {names}')
        if not self._checks:
            return CheckResult(failure='built-in dataset checks are not installed')

        findings: list[Finding] = []
        failures: list[str] = []
        failed_checks: set[str] = set()
        for registered in self._checks:
            failed_needs = registered.declaration.needs & failed_checks
            if failed_needs:
                names = ', '.join(sorted(failed_needs))
                failures.append(
                    f"check '{registered.declaration.name}' did not run because prerequisite "
                    f'check(s) failed: {names}'
                )
                failed_checks.add(registered.declaration.name)
                continue
            try:
                produced = registered.function(context)
                for finding in produced:
                    if not isinstance(finding, Finding):
                        raise TypeError(
                            f"check '{registered.declaration.name}' returned a non-Finding value"
                        )
                    findings.append(finding)
            except SystemExit as error:
                failed_checks.add(registered.declaration.name)
                failures.append(
                    f"check '{registered.declaration.name}' attempted to exit with status {error.code}"
                )
            except Exception as error:
                failed_checks.add(registered.declaration.name)
                failures.append(
                    f"check '{registered.declaration.name}' failed: {type(error).__name__}: {error}"
                )

        return CheckResult(tuple(findings), '; '.join(failures) or None)


class PlainTextRenderer:
    def render(self, result: CheckResult) -> str:
        lines = []
        for finding in result.findings:
            location = f' {finding.path}' if finding.path is not None else ''
            lines.append(
                f'{finding.severity.value.upper()} [{finding.check_id}]{location}: {finding.message}'
            )

        counts = {
            severity: sum(finding.severity is severity for finding in result.findings)
            for severity in Severity
        }
        if result.failure is not None:
            lines.append(f'ERROR [infrastructure]: {result.failure}')
            lines.append('Summary: checks did not complete.')
        else:
            lines.append(
                'Summary: '
                f'{counts[Severity.ERROR]} error(s), '
                f'{counts[Severity.WARNING]} warning(s), '
                f'{counts[Severity.INFO]} info finding(s).'
            )

        return '\n'.join(lines) + '\n'

    def write(self, result: CheckResult, stream: TextIO) -> None:
        stream.write(self.render(result))


def run_checks(
    factory: CheckFactory,
    context: CheckContext,
    *,
    scope: CheckScope,
) -> CheckResult:
    return CheckRunner(factory, scope=scope).run(context)


def _finding(
    check_id: str,
    severity: Severity,
    message: str,
    *,
    path: Path | None = None,
    details: Mapping[str, Any] | None = None,
) -> Finding:
    return Finding(check_id, severity, message, path, details or {})


def _relative_path(root: Path, path: Path) -> Path:
    try:
        return path.relative_to(root)
    except ValueError:
        return path


def _resolve_inside(root: Path, value: Path) -> Path | None:
    candidate = value if value.is_absolute() else root / value
    candidate = candidate.resolve()
    try:
        candidate.relative_to(root)
    except ValueError:
        return None
    return candidate


def check_analysis_file_selection(context: CheckContext) -> Iterable[Finding]:
    check_id = 'analysis-files.selection'
    root = context.dataset_root.resolve()
    context.dataset_root = root
    context.cache[CacheKey.ANALYSIS_SELECTION_VALID] = False

    if not root.is_dir():
        yield _finding(check_id, Severity.ERROR, 'Dataset root is not a directory.', path=root)
        return

    requested = context.analysis_paths
    if not requested:
        default = root / 'analysis.yaml'
        if not default.is_file():
            yield _finding(
                check_id,
                Severity.ERROR,
                'No root analysis.yaml exists; supply every analysis file explicitly with '
                'repeatable --analysis-file options.',
                path=default,
            )
            return
        selected = [default.resolve()]
    else:
        selected = []
        seen: set[Path] = set()
        invalid = False
        for requested_path in requested:
            resolved = _resolve_inside(root, requested_path)
            display_path = requested_path if requested_path.is_absolute() else root / requested_path
            if resolved is None:
                invalid = True
                yield _finding(
                    check_id,
                    Severity.ERROR,
                    'Selected analysis file is outside the dataset root.',
                    path=display_path,
                )
                continue
            if resolved in seen:
                invalid = True
                yield _finding(
                    check_id,
                    Severity.ERROR,
                    'Selected analysis file is duplicated.',
                    path=_relative_path(root, resolved),
                )
                continue
            seen.add(resolved)
            if resolved.suffix.lower() not in ('.yaml', '.yml'):
                invalid = True
                yield _finding(
                    check_id,
                    Severity.ERROR,
                    'Selected analysis file must have a .yaml or .yml suffix.',
                    path=_relative_path(root, resolved),
                )
                continue
            if not resolved.is_file():
                invalid = True
                yield _finding(
                    check_id,
                    Severity.ERROR,
                    'Selected analysis path is not a regular file.',
                    path=_relative_path(root, resolved),
                )
                continue
            selected.append(resolved)
        if invalid:
            return

    requested_main = context.main_analysis_path
    if len(selected) == 1:
        main = selected[0]
        if requested_main is not None:
            resolved_main = _resolve_inside(root, requested_main)
            if resolved_main != main:
                yield _finding(
                    check_id,
                    Severity.ERROR,
                    '--main-analysis-file conflicts with the automatically selected analysis file.',
                    path=requested_main,
                )
                return
    else:
        if requested_main is None:
            yield _finding(
                check_id,
                Severity.ERROR,
                'Multiple analysis files require --main-analysis-file.',
            )
            return
        main = _resolve_inside(root, requested_main)
        if main is None or main not in selected:
            yield _finding(
                check_id,
                Severity.ERROR,
                '--main-analysis-file must identify one of the selected analysis files.',
                path=requested_main,
            )
            return

    context.analysis_paths = tuple(selected)
    context.main_analysis_path = main
    context.cache[CacheKey.ANALYSIS_SELECTION_VALID] = True


def _usable_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _usable_name(value: Any) -> bool:
    return _usable_string(value) and any(character.isalnum() for character in value)


def _load_yaml(path: Path) -> tuple[Any, str | None]:
    try:
        return load_yaml_file(path), None
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        return None, str(error)


def _author_names(metadata: Mapping[str, Any]) -> tuple[list[str], list[str]]:
    authors = metadata.get('authors')
    if not isinstance(authors, list):
        return [], ['metadata.authors must be a list']
    if not authors:
        return [], ['metadata.authors must not be empty']

    names: list[str] = []
    errors: list[str] = []
    for index, author in enumerate(authors):
        if not isinstance(author, Mapping):
            errors.append(f'metadata.authors[{index}].name must be a nonempty string')
            continue
        if not all(isinstance(key, str) for key in author):
            errors.append(f'metadata.authors[{index}] must use string keys')
            continue

        description = MetadataAuthorDescription.from_dict(**author)
        for diagnostic in description.validate_structure():
            errors.append(
                f'metadata.authors[{index}].{diagnostic.location()}: {diagnostic.message}'
            )
        if not _usable_name(author.get('name')):
            errors.append(f'metadata.authors[{index}].name must be a nonempty string')
        else:
            names.append(author['name'].strip())
    return names, errors


def check_analysis_metadata(context: CheckContext) -> Iterable[Finding]:
    check_id = 'analysis-metadata.parse'
    context.cache[CacheKey.ANALYSIS_METADATA_VALID] = False
    if not context.cache[CacheKey.ANALYSIS_SELECTION_VALID]:
        yield _finding(check_id, Severity.INFO, 'Skipped because analysis-file selection failed.')
        return

    metadata_by_path: dict[Path, Mapping[str, Any]] = {}
    author_names: dict[Path, tuple[str, ...]] = {}
    valid = True
    for path in context.analysis_paths:
        display = _relative_path(context.dataset_root, path)
        document, error = _load_yaml(path)
        if error is not None:
            valid = False
            yield _finding(
                check_id, Severity.ERROR, f'Cannot parse YAML: {error}', path=display,
                details={'field': '$'},
            )
            continue
        if not isinstance(document, Mapping):
            valid = False
            yield _finding(
                check_id, Severity.ERROR, 'Analysis YAML must be a top-level mapping.', path=display,
                details={'field': '$'},
            )
            continue
        metadata = document.get('metadata')
        if not isinstance(metadata, Mapping):
            valid = False
            yield _finding(
                check_id, Severity.ERROR, 'metadata must be a mapping.', path=display,
                details={'field': 'metadata'},
            )
            continue
        metadata_by_path[path] = metadata

        analysis_id = metadata.get('id')
        if not _usable_string(analysis_id):
            valid = False
            yield _finding(
                check_id, Severity.ERROR, 'metadata.id must be a nonempty string.', path=display,
                details={'field': 'metadata.id'},
            )
        elif _ANALYSIS_ID.fullmatch(analysis_id) is None:
            valid = False
            yield _finding(
                check_id, Severity.ERROR,
                'metadata.id must have exactly the YYYY-NN form without a revision suffix.',
                path=display, details={'field': 'metadata.id', 'value': analysis_id},
            )

        names, author_errors = _author_names(metadata)
        if author_errors:
            valid = False
            for message in author_errors:
                yield _finding(
                    check_id, Severity.ERROR, message + '.', path=display,
                    details={'field': message.split(' must', 1)[0]},
                )
        else:
            author_names[path] = tuple(names)

        if path == context.main_analysis_path:
            # 'description' is Optional on MetadataDescription -- an ongoing or unofficial
            # analysis need not have one -- but the release process requires it.
            for field_name in ('title', 'description'):
                if not _usable_string(metadata.get(field_name)):
                    valid = False
                    yield _finding(
                        check_id, Severity.ERROR,
                        f'metadata.{field_name} must be a nonempty string in the main analysis file.',
                        path=display, details={'field': f'metadata.{field_name}'},
                    )

    context.cache[CacheKey.ANALYSIS_METADATA] = metadata_by_path
    context.cache[CacheKey.ANALYSIS_AUTHOR_NAMES] = author_names
    context.cache[CacheKey.ANALYSIS_METADATA_VALID] = valid


def _normalize_name(name: str) -> str:
    value = unicodedata.normalize('NFKC', name).casefold().strip()
    if value.count(',') == 1:
        family, given = value.split(',', 1)
        if family.strip() and given.strip():
            value = f'{given} {family}'
    value = ''.join(
        character
        if unicodedata.category(character)[0] in ('L', 'M', 'N')
        else ' '
        for character in value
    )
    return _WHITESPACE.sub(' ', value).strip()


def _unique_normalized(names: Iterable[str]) -> tuple[dict[str, str], list[str]]:
    normalized: dict[str, str] = {}
    duplicates: list[str] = []
    for name in names:
        key = _normalize_name(name)
        if key in normalized:
            duplicates.append(name)
        else:
            normalized[key] = name
    return normalized, duplicates


def check_analysis_consistency(context: CheckContext) -> Iterable[Finding]:
    check_id = 'analysis-metadata.consistency'
    context.cache[CacheKey.ANALYSIS_CONSISTENCY_VALID] = False
    if not context.cache[CacheKey.ANALYSIS_METADATA_VALID]:
        yield _finding(check_id, Severity.INFO, 'Skipped because analysis metadata is invalid.')
        return

    metadata_by_path = context.cache[CacheKey.ANALYSIS_METADATA]
    authors_by_path = context.cache[CacheKey.ANALYSIS_AUTHOR_NAMES]
    ids = {path: metadata_by_path[path]['id'] for path in context.analysis_paths}
    common_id = ids[context.analysis_paths[0]]
    valid = True
    for path, analysis_id in ids.items():
        if analysis_id != common_id:
            valid = False
            yield _finding(
                check_id, Severity.ERROR,
                f"Analysis ID '{analysis_id}' differs from common ID '{common_id}'.",
                path=_relative_path(context.dataset_root, path),
                details={'field': 'metadata.id', 'expected': common_id, 'actual': analysis_id},
            )

    normalized_authors: dict[Path, dict[str, str]] = {}
    for path, names in authors_by_path.items():
        normalized, duplicates = _unique_normalized(names)
        normalized_authors[path] = normalized
        for duplicate in duplicates:
            valid = False
            yield _finding(
                check_id, Severity.ERROR,
                f"Duplicate or unresolved analysis author '{duplicate}'.",
                path=_relative_path(context.dataset_root, path),
                details={'field': 'metadata.authors', 'name': duplicate},
            )

    common_authors = set(normalized_authors[context.analysis_paths[0]])
    for path in context.analysis_paths[1:]:
        current_authors = set(normalized_authors[path])
        if current_authors != common_authors:
            valid = False
            yield _finding(
                check_id, Severity.ERROR,
                'Analysis author set differs from the common author set.',
                path=_relative_path(context.dataset_root, path),
                details={
                    'field': 'metadata.authors',
                    'missing': sorted(common_authors - current_authors),
                    'additional': sorted(current_authors - common_authors),
                },
            )

    if valid:
        context.cache[CacheKey.COMMON_ANALYSIS_ID] = common_id
        context.cache[CacheKey.COMMON_AUTHORS_NORMALIZED] = normalized_authors[
            context.analysis_paths[0]
        ]
        context.cache[CacheKey.ANALYSIS_CONSISTENCY_VALID] = True


def check_zenodo_metadata(context: CheckContext) -> Iterable[Finding]:
    check_id = 'zenodo-metadata.basic'
    context.cache[CacheKey.ZENODO_CREATORS_VALID] = False
    path = context.dataset_root / '.zenodo.json'
    if not path.is_file():
        yield _finding(check_id, Severity.ERROR, 'Required root .zenodo.json is missing.', path=path)
        return
    try:
        with path.open('r', encoding='utf-8') as stream:
            document = json.load(stream)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        yield _finding(check_id, Severity.ERROR, f'Cannot parse JSON: {error}', path=path)
        return
    if not isinstance(document, Mapping):
        yield _finding(check_id, Severity.ERROR, '.zenodo.json must contain a JSON object.', path=path)
        return
    creators = document.get('creators')
    creator_names: list[str] = []
    creators_valid = True
    if not isinstance(creators, list) or not creators:
        creators_valid = False
        yield _finding(
            check_id, Severity.ERROR, 'creators must be a nonempty list.', path=path,
            details={'field': 'creators'},
        )
    else:
        for index, creator in enumerate(creators):
            if not isinstance(creator, Mapping) or not _usable_name(creator.get('name')):
                creators_valid = False
                yield _finding(
                    check_id, Severity.ERROR,
                    f'creators[{index}].name must be a nonempty string.', path=path,
                    details={'field': f'creators[{index}].name'},
                )
            else:
                creator_names.append(creator['name'].strip())
    if creators_valid:
        normalized, duplicates = _unique_normalized(creator_names)
        if duplicates:
            creators_valid = False
            for duplicate in duplicates:
                yield _finding(
                    check_id, Severity.ERROR, f"Duplicate creator '{duplicate}'.", path=path,
                    details={'field': 'creators', 'name': duplicate},
                )
        else:
            context.cache[CacheKey.ZENODO_CREATORS_NORMALIZED] = normalized
            context.cache[CacheKey.ZENODO_CREATORS_VALID] = True

    title = document.get('title')
    if not isinstance(title, str):
        yield _finding(
            check_id, Severity.ERROR, 'title must be a string.', path=path,
            details={'field': 'title'},
        )
    else:
        match = _ZENODO_TITLE.fullmatch(title)
        if match is None:
            yield _finding(
                check_id, Severity.ERROR,
                'title must exactly match "EOS/DATA-{dataset-id}: Supplementary material for '
                'EOS/ANALYSIS-{analysis-id}".',
                path=path, details={'field': 'title', 'value': title},
            )
        else:
            dataset_id = match.group('dataset_id')
            if _DATASET_ID.fullmatch(dataset_id) is None:
                yield _finding(
                    check_id, Severity.ERROR,
                    'The dataset ID embedded in title must be YYYY-NN or YYYY-NNvN with N >= 2.',
                    path=path, details={'field': 'title', 'dataset_id': dataset_id},
                )
            if context.cache[CacheKey.ANALYSIS_CONSISTENCY_VALID]:
                expected = context.cache[CacheKey.COMMON_ANALYSIS_ID]
                actual = match.group('analysis_id')
                if actual != expected:
                    yield _finding(
                        check_id, Severity.ERROR,
                        f"The analysis ID embedded in title is '{actual}', expected '{expected}'.",
                        path=path,
                        details={'field': 'title', 'expected': expected, 'actual': actual},
                    )
            else:
                yield _finding(
                    check_id, Severity.INFO,
                    'Skipped embedded analysis-ID comparison because analysis metadata is invalid.',
                    path=path,
                )

    description = document.get('description')
    if not isinstance(description, str):
        yield _finding(
            check_id, Severity.ERROR, 'description must be a string.', path=path,
            details={'field': 'description'},
        )
    elif context.cache[CacheKey.ANALYSIS_METADATA_VALID]:
        main_metadata = context.cache[CacheKey.ANALYSIS_METADATA][context.main_analysis_path]
        expected = main_metadata['description']
        if description != expected:
            yield _finding(
                check_id, Severity.ERROR,
                'description does not exactly equal main analysis metadata.description.',
                path=path, details={'field': 'description'},
            )
    else:
        yield _finding(
            check_id, Severity.INFO,
            'Skipped description comparison because main analysis metadata is invalid.', path=path,
        )


def _without_diacritics(value: str) -> str:
    return ''.join(
        character
        for character in unicodedata.normalize('NFKD', value)
        if not unicodedata.combining(character)
    )


def _initial_parts_match(left_parts: list[str], right_parts: list[str]) -> bool:
    if len(left_parts) != len(right_parts) or left_parts == right_parts:
        return False
    for left_part, right_part in zip(left_parts, right_parts):
        if left_part == right_part:
            continue
        if len(left_part) == 1 and right_part.startswith(left_part):
            continue
        if len(right_part) == 1 and left_part.startswith(right_part):
            continue
        return False
    return True


def _ambiguous_match(left: str, right: str) -> tuple[str, ...]:
    left_plain = _without_diacritics(left)
    right_plain = _without_diacritics(right)
    left_parts = left_plain.split()
    right_parts = right_plain.split()
    reasons: list[str] = []

    if left_plain == right_plain and left != right:
        reasons.append('diacritics')
    elif sorted(left_parts) == sorted(right_parts) and left_parts != right_parts:
        reasons.append('name-order')
        if left_plain != left or right_plain != right:
            reasons.append('diacritics')
    elif _initial_parts_match(left_parts, right_parts):
        reasons.append('initials')
        if left_plain != left or right_plain != right:
            reasons.append('diacritics')
    elif _initial_parts_match(left_parts, list(reversed(right_parts))):
        reasons.extend(('initials', 'name-order'))
        if left_plain != left or right_plain != right:
            reasons.append('diacritics')

    return tuple(reasons)


def check_creator_consistency(context: CheckContext) -> Iterable[Finding]:
    check_id = 'authors.creator-consistency'
    if not context.cache[CacheKey.ANALYSIS_CONSISTENCY_VALID]:
        yield _finding(check_id, Severity.INFO, 'Skipped because analysis author metadata is invalid.')
        return
    if not context.cache[CacheKey.ZENODO_CREATORS_VALID]:
        yield _finding(check_id, Severity.INFO, 'Skipped because Zenodo creator metadata is invalid.')
        return

    authors = context.cache[CacheKey.COMMON_AUTHORS_NORMALIZED]
    creators = context.cache[CacheKey.ZENODO_CREATORS_NORMALIZED]
    unmatched_authors = set(authors) - set(creators)
    unmatched_creators = set(creators) - set(authors)
    candidates = [
        (author, creator, reasons)
        for author in sorted(unmatched_authors)
        for creator in sorted(unmatched_creators)
        if (reasons := _ambiguous_match(author, creator))
    ]
    ambiguous_authors = {author for author, _creator, _reasons in candidates}
    ambiguous_creators = {creator for _author, creator, _reasons in candidates}

    for author, creator, reasons in candidates:
        yield Finding(
            check_id,
            Severity.WARNING,
            f"Plausible but ambiguous author/creator match: '{authors[author]}' and "
            f"'{creators[creator]}'.",
            context.dataset_root / '.zenodo.json',
            {
                'analysis_author': authors[author],
                'zenodo_creator': creators[creator],
                'analysis_normalized': author,
                'creator_normalized': creator,
                'reasons': reasons,
            },
            InteractiveQuestion(
                'Do these names identify the same person?',
                {
                    'analysis_author': authors[author],
                    'zenodo_creator': creators[creator],
                    'reasons': reasons,
                },
            ),
        )

    for author in sorted(unmatched_authors - ambiguous_authors):
        yield _finding(
            check_id, Severity.ERROR,
            f"Analysis author '{authors[author]}' is missing from Zenodo creators.",
            path=context.dataset_root / '.zenodo.json',
            details={'analysis_author': authors[author]},
        )
    for creator in sorted(unmatched_creators - ambiguous_creators):
        yield _finding(
            check_id, Severity.ERROR,
            f"Zenodo creator '{creators[creator]}' is not an analysis author.",
            path=context.dataset_root / '.zenodo.json',
            details={'zenodo_creator': creators[creator]},
        )


def register_basic_checks(factory: CheckFactory) -> None:
    checks = (
        ('analysis-files.selection', (), check_analysis_file_selection),
        ('analysis-metadata.parse', ('analysis-files.selection',), check_analysis_metadata),
        (
            'analysis-metadata.consistency',
            ('analysis-metadata.parse',),
            check_analysis_consistency,
        ),
        (
            'zenodo-metadata.basic',
            ('analysis-metadata.consistency',),
            check_zenodo_metadata,
        ),
        (
            'authors.creator-consistency',
            ('analysis-metadata.consistency', 'zenodo-metadata.basic'),
            check_creator_consistency,
        ),
    )
    for name, needs, function in checks:
        factory.declare(name, scope=CheckScope.BASIC, needs=needs)
        factory.register(name, function)
