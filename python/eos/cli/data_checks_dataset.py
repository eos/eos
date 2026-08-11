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

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import html
import inspect
import os
from pathlib import Path, PurePosixPath
import re
import string
from typing import Any
from urllib.parse import unquote, urlsplit

import yaml

from ..diagnostic import Severity
from .data_checks import (
    CacheKey, CheckContext, CheckFactory, CheckScope, Finding,
    _finding, _relative_path, _usable_string,
)


MAX_FILE_SIZE = 100 * 1024 * 1024
MAX_DATASET_SIZE = 1024 * 1024 * 1024
_ORCID = re.compile(r'^(\d{4})-(\d{4})-(\d{4})-(\d{3}[\dX])$')
_MARKDOWN_TARGET = re.compile(r'!?\[[^\]]*\]\((?:<([^>]+)>|([^\s)]+))(?:\s+[^)]*)?\)')
_HTML_TARGET = re.compile(r'''(?:href|src)\s*=\s*["']([^"']+)["']''', re.IGNORECASE)
_RECOGNIZED_OUTPUT_TYPES = frozenset({
    'DynestyResults', 'ImportanceSamples', 'MarkovChain', 'MixtureDensity', 'Mode',
    'NabuLikelihood', 'PMCSampler', 'Prediction', 'SampleMask',
})


@dataclass(frozen=True)
class FileInventoryEntry:
    path: Path
    size: int


def _inventory(context: CheckContext) -> tuple[FileInventoryEntry, ...]:
    cached = context.cache.get(CacheKey.FILE_INVENTORY)
    if cached is not None:
        return tuple(cached)

    root = context.dataset_root
    entries: list[FileInventoryEntry] = []
    for current, directories, files in os.walk(root, followlinks=False):
        current_path = Path(current)
        directories[:] = sorted(
            directory for directory in directories
            if not (current_path == root and directory == '.git')
            and not (current_path / directory).is_symlink()
        )
        for filename in sorted(files):
            path = current_path / filename
            try:
                if not path.is_symlink() and path.is_file():
                    entries.append(FileInventoryEntry(path.relative_to(root), path.stat().st_size))
            except OSError:
                continue
    result = tuple(sorted(entries, key=lambda entry: entry.path.as_posix()))
    context.cache[CacheKey.FILE_INVENTORY] = result
    return result


def _valid_orcid(value: str) -> bool:
    match = _ORCID.fullmatch(value)
    if match is None:
        return False
    digits = ''.join(match.groups())
    total = 0
    for digit in digits[:15]:
        total = (total + int(digit)) * 2
    remainder = (12 - total % 11) % 11
    expected = 'X' if remainder == 10 else str(remainder)
    return digits[-1] == expected


def check_zenodo_completeness(context: CheckContext) -> Iterable[Finding]:
    check_id = 'zenodo-metadata.completeness'
    document = context.cache.get(CacheKey.ZENODO_METADATA)
    if not isinstance(document, Mapping):
        yield _finding(check_id, Severity.INFO, 'Skipped because Zenodo metadata is invalid.')
        return
    path = context.dataset_root / '.zenodo.json'
    if document.get('upload_type') != 'dataset':
        yield _finding(
            check_id, Severity.ERROR, 'upload_type must be exactly "dataset".',
            path=path, details={'field': 'upload_type'},
        )
    for field_name in ('title', 'description', 'license'):
        if not _usable_string(document.get(field_name)):
            yield _finding(
                check_id, Severity.ERROR, f'{field_name} must be a nonempty string.',
                path=path, details={'field': field_name},
            )

    creators = document.get('creators')
    if isinstance(creators, list):
        for index, creator in enumerate(creators):
            if not isinstance(creator, Mapping):
                continue
            if not _usable_string(creator.get('affiliation')):
                yield _finding(
                    check_id, Severity.ERROR,
                    f'creators[{index}].affiliation must be a nonempty string.',
                    path=path, details={'field': f'creators[{index}].affiliation'},
                )
            if 'orcid' in creator and (not isinstance(creator['orcid'], str) or not _valid_orcid(creator['orcid'])):
                yield _finding(
                    check_id, Severity.ERROR, f'creators[{index}].orcid must be a valid ORCID.',
                    path=path,
                    details={
                        'field': f'creators[{index}].orcid', 'value': creator.get('orcid'),
                    },
                )

    grants = document.get('grants')
    if not isinstance(grants, list):
        yield _finding(
            check_id, Severity.ERROR, 'grants must be a list (an empty list is permitted).',
            path=path, details={'field': 'grants'},
        )
        return
    grant_ids: list[str] = []
    malformed = False
    for index, grant in enumerate(grants):
        if not isinstance(grant, Mapping) or not _usable_string(grant.get('id')):
            malformed = True
            yield _finding(
                check_id, Severity.ERROR, f'grants[{index}].id must be a nonempty string.',
                path=path, details={'field': f'grants[{index}].id'},
            )
        else:
            grant_ids.append(grant['id'].strip())
    if not malformed:
        yield _finding(
            check_id, Severity.WARNING, 'Grant completeness has not been confirmed.',
            path=path, details={'grant_ids': tuple(grant_ids)},
        )


def _figure_files(context: CheckContext) -> dict[PurePosixPath, list[str]]:
    figures: dict[PurePosixPath, list[str]] = defaultdict(list)
    for entry in _inventory(context):
        parts = entry.path.parts
        if len(parts) < 2 or parts[0] != 'figures':
            continue
        suffix = entry.path.suffix.lower()
        if not suffix:
            continue
        figures[PurePosixPath(entry.path.with_suffix('').as_posix())].append(suffix)
    return dict(sorted(figures.items(), key=lambda item: item[0].as_posix()))


def _readme_targets(text: str) -> list[str]:
    markdown = [first or second for first, second in _MARKDOWN_TARGET.findall(text)]
    return markdown + _HTML_TARGET.findall(text)


def _local_target(value: str) -> tuple[PurePosixPath | None, bool]:
    value = html.unescape(value.strip())
    split = urlsplit(value)
    if split.scheme or split.netloc or not split.path:
        return None, False
    raw = unquote(split.path).replace('\\', '/')
    path = PurePosixPath(raw.lstrip('/'))
    normalized: list[str] = []
    for part in path.parts:
        if part in ('', '.'):
            continue
        if part == '..':
            if not normalized:
                return None, True
            normalized.pop()
        else:
            normalized.append(part)
    return PurePosixPath(*normalized), False


def check_readme_figures(context: CheckContext) -> Iterable[Finding]:
    check_id = 'figures.readme'
    path = context.dataset_root / 'README.md'
    figures = _figure_files(context)
    if not path.is_file() or path.is_symlink():
        yield _finding(check_id, Severity.ERROR, 'Required root README.md is missing.', path=path)
        targets: list[str] = []
    else:
        try:
            targets = _readme_targets(path.read_text(encoding='utf-8'))
        except (OSError, UnicodeError) as error:
            yield _finding(check_id, Severity.ERROR, f'Cannot read README.md: {error}', path=path)
            targets = []

    references: Counter[tuple[PurePosixPath, str]] = Counter()
    for target in targets:
        local, traversal = _local_target(target)
        if traversal:
            yield _finding(
                check_id, Severity.ERROR,
                f'README reference traverses outside the dataset root: {target}',
                path=Path('README.md'), details={'target': target},
            )
            continue
        if local is None or len(local.parts) < 2 or local.parts[0] != 'figures' or not local.suffix:
            continue
        references[(local.with_suffix(''), local.suffix.lower())] += 1

    for stem, suffixes in figures.items():
        suffix_set = set(suffixes)
        if '.pdf' not in suffix_set:
            yield _finding(
                check_id, Severity.ERROR, 'Physical figure has no PDF representation.',
                path=Path(stem.as_posix()), details={'figure': stem.as_posix()},
            )
        if '.png' not in suffix_set:
            yield _finding(
                check_id, Severity.WARNING, 'Physical figure has no PNG representation.',
                path=Path(stem.as_posix()), details={'figure': stem.as_posix()},
            )
        for suffix in sorted(suffix_set - {'.pdf', '.png'}):
            yield _finding(
                check_id, Severity.WARNING, f'Physical figure uses additional format {suffix}.',
                path=Path(stem.as_posix() + suffix),
                details={'figure': stem.as_posix(), 'format': suffix},
            )
        if not any(reference_stem == stem for reference_stem, _suffix in references):
            yield _finding(
                check_id, Severity.ERROR, 'Physical figure is not listed in README.md.',
                path=Path(stem.as_posix()), details={'figure': stem.as_posix()},
            )
        duplicate_count = sum(
            max(0, count - 1)
            for (reference_stem, _suffix), count in references.items()
            if reference_stem == stem
        )
        if duplicate_count:
            yield _finding(
                check_id, Severity.WARNING,
                'Logical figure is listed redundantly in README.md.', path=Path('README.md'),
                details={'figure': stem.as_posix(), 'duplicate_count': duplicate_count},
            )


def check_analysis_figures(context: CheckContext) -> Iterable[Finding]:
    check_id = 'figures.analysis-declarations'
    if (
        not context.cache.get(CacheKey.ANALYSIS_METADATA_VALID, False)
        or not context.cache.get(CacheKey.ANALYSIS_CONSISTENCY_VALID, False)
    ):
        yield _finding(
            check_id, Severity.INFO,
            'Skipped because not all selected analysis files passed validation.',
        )
        return
    documents = context.cache.get(CacheKey.ANALYSIS_DOCUMENTS)
    if not isinstance(documents, Mapping):
        yield _finding(check_id, Severity.INFO, 'Skipped because analysis documents are invalid.')
        return
    available = {entry.path.as_posix() for entry in _inventory(context)}
    for analysis_path in context.analysis_paths:
        document = documents[analysis_path]
        figures = document.get('figures', [])
        if not isinstance(figures, list):
            yield _finding(
                check_id, Severity.ERROR, 'figures must be a list.',
                path=_relative_path(context.dataset_root, analysis_path),
                details={'field': 'figures'},
            )
            continue
        for index, figure in enumerate(figures):
            name = figure.get('name') if isinstance(figure, Mapping) else None
            if not _usable_string(name):
                yield _finding(
                    check_id, Severity.ERROR,
                    f'figures[{index}].name must be a nonempty string.',
                    path=_relative_path(context.dataset_root, analysis_path),
                    details={'field': f'figures[{index}].name'},
                )
                continue
            stem, traversal = _local_target(f'figures/{name}')
            details = {'analysis_path': _relative_path(context.dataset_root, analysis_path).as_posix(), 'figure': name}
            if traversal or stem is None or stem.parts[0] != 'figures':
                yield _finding(
                    check_id, Severity.ERROR, f'Invalid declared figure name: {name}',
                    path=_relative_path(context.dataset_root, analysis_path), details=details,
                )
                continue
            for suffix, severity in (('.pdf', Severity.ERROR), ('.png', Severity.WARNING)):
                expected = stem.as_posix() + suffix
                if expected not in available:
                    yield _finding(
                        check_id, severity,
                        f'Declared figure is missing its {suffix[1:].upper()} file.',
                        path=Path(expected), details={**details, 'output_path': expected},
                    )


def check_dataset_sizes(context: CheckContext) -> Iterable[Finding]:
    check_id = 'dataset.size'
    inventory = _inventory(context)
    total = sum(entry.size for entry in inventory)
    for entry in inventory:
        if entry.size > MAX_FILE_SIZE:
            yield _finding(
                check_id, Severity.ERROR, f'File exceeds 100 MiB ({entry.size} bytes).',
                path=entry.path,
                details={'size_bytes': entry.size, 'limit_bytes': MAX_FILE_SIZE},
            )
    yield _finding(check_id, Severity.INFO, f'Dataset contains exactly {total} bytes.', details={'total_bytes': total})
    if total > MAX_DATASET_SIZE:
        yield _finding(
            check_id, Severity.WARNING,
            'Dataset exceeds 1 GiB; remove unnecessary importance samples and posterior predictions.',
            details={'total_bytes': total, 'limit_bytes': MAX_DATASET_SIZE},
        )


def _normalize_output_path(template: str, arguments: Mapping[str, Any]) -> tuple[PurePosixPath | None, str | None]:
    try:
        fields = [field for _literal, field, _spec, _conversion in string.Formatter().parse(template) if field]
        if any('.' in field or '[' in field for field in fields):
            return None, 'attribute and item access are not permitted in output templates'
        rendered = template.format(**arguments).replace('\\', '/')
    except (KeyError, ValueError, TypeError, IndexError) as error:
        return None, str(error)
    path = PurePosixPath(rendered)
    if path.is_absolute() or any(part in ('', '..') for part in path.parts):
        return None, 'output path is absolute or traverses outside the dataset root'
    return path, None


def _output_template_patterns(templates: Mapping[str, str]) -> tuple[re.Pattern, ...]:
    patterns = []
    for template in templates.values():
        if not template.startswith('data/'):
            continue
        try:
            parts = string.Formatter().parse(template)
            expression = ''.join(
                re.escape(literal) + ('[^/]+' if field is not None else '')
                for literal, field, _spec, _conversion in parts
            )
        except ValueError:
            continue
        patterns.append(re.compile(f'^{expression}$'))
    return tuple(patterns)


def _recognized_outputs(
    context: CheckContext,
    templates: Mapping[str, str],
) -> tuple[PurePosixPath, ...]:
    results: list[PurePosixPath] = []
    layout_patterns = _output_template_patterns(templates)
    for entry in _inventory(context):
        if entry.path.name != 'description.yaml' or not entry.path.parts or entry.path.parts[0] != 'data':
            continue
        try:
            document = yaml.safe_load((context.dataset_root / entry.path).read_text(encoding='utf-8'))
        except (OSError, UnicodeError, yaml.YAMLError):
            continue
        if not isinstance(document, Mapping):
            continue
        output = PurePosixPath(entry.path.parent.as_posix())
        if (
            document.get('type') in _RECOGNIZED_OUTPUT_TYPES
            or any(pattern.fullmatch(output.as_posix()) for pattern in layout_patterns)
        ):
            results.append(output)
    return tuple(sorted(set(results), key=lambda path: path.as_posix()))


def check_reproducible_outputs(context: CheckContext) -> Iterable[Finding]:
    check_id = 'outputs.reproducibility'
    if (
        not context.cache.get(CacheKey.ANALYSIS_METADATA_VALID, False)
        or not context.cache.get(CacheKey.ANALYSIS_CONSISTENCY_VALID, False)
    ):
        yield _finding(
            check_id, Severity.INFO,
            'Skipped because not all selected analysis files passed validation.',
        )
        return
    documents = context.cache.get(CacheKey.ANALYSIS_DOCUMENTS)
    if not isinstance(documents, Mapping):
        yield _finding(check_id, Severity.INFO, 'Skipped because analysis documents are invalid.')
        return
    from .. import tasks as task_registry
    from ..analysis_file_description import InvalidComponent, TaskComponent

    templates = task_registry.task_output_templates()
    available_files = {entry.path.as_posix() for entry in _inventory(context)}
    available_directories = {
        PurePosixPath(*PurePosixPath(path).parts[:index]).as_posix()
        for path in available_files
        for index in range(1, len(PurePosixPath(path).parts))
    }
    claims: dict[PurePosixPath, list[dict[str, str]]] = defaultdict(list)
    for analysis_path in context.analysis_paths:
        document = documents[analysis_path]
        steps = document.get('steps', [])
        if steps is None:
            steps = []
        if not isinstance(steps, list):
            yield _finding(
                check_id, Severity.ERROR, 'steps must be a list.',
                path=_relative_path(context.dataset_root, analysis_path),
                details={'field': 'steps'},
            )
            continue
        for step_index, step in enumerate(steps):
            if not isinstance(step, Mapping):
                yield _finding(
                    check_id, Severity.ERROR, f'steps[{step_index}] must be a mapping.',
                    path=_relative_path(context.dataset_root, analysis_path),
                )
                continue
            step_id = step.get('id', str(step_index))
            defaults = step.get('default_arguments', {})
            task_entries = step.get('tasks', [])
            if not isinstance(defaults, Mapping) or not isinstance(task_entries, list):
                yield _finding(
                    check_id, Severity.ERROR,
                    f'Step {step_id} has malformed defaults or tasks.',
                    path=_relative_path(context.dataset_root, analysis_path),
                    details={'step': step_id},
                )
                continue
            for task in task_entries:
                if not isinstance(task, Mapping):
                    yield _finding(
                        check_id, Severity.ERROR, 'Task invocation must be a mapping.',
                        path=_relative_path(context.dataset_root, analysis_path),
                        details={'analysis_path': _relative_path(
                            context.dataset_root, analysis_path,
                        ).as_posix(), 'step': str(step_id)},
                    )
                    continue
                raw_task_name = task.get('task')
                raw_arguments = task.get('arguments', {})
                if (
                    not isinstance(raw_task_name, str)
                    or not isinstance(raw_arguments, Mapping)
                    or not all(isinstance(key, str) for key in raw_arguments)
                ):
                    yield _finding(
                        check_id, Severity.ERROR, 'Task name or arguments are malformed.',
                        path=_relative_path(context.dataset_root, analysis_path),
                        details={'analysis_path': _relative_path(
                            context.dataset_root, analysis_path,
                        ).as_posix(), 'step': str(step_id)},
                    )
                    continue
                task_component = TaskComponent.from_dict(**task)
                task_name = getattr(task_component, 'task', None)
                details = {
                    'analysis_path': _relative_path(
                        context.dataset_root, analysis_path,
                    ).as_posix(),
                    'step': str(step_id),
                    'task': str(task_name),
                }
                if task_name not in templates:
                    yield _finding(
                        check_id, Severity.ERROR, f'Unknown task {task_name!r}.',
                        path=_relative_path(context.dataset_root, analysis_path),
                        details=details,
                    )
                    continue
                template = templates[task_name]
                if not template or not template.startswith('data/'):
                    continue
                function = task_registry._tasks[task_name]
                arguments = {
                    name: parameter.default
                    for name, parameter in inspect.signature(function).parameters.items()
                    if parameter.default is not inspect.Parameter.empty
                }
                task_defaults = defaults.get(task_name, {})
                task_arguments = getattr(task_component, 'arguments', {})
                if isinstance(task_component, InvalidComponent):
                    task_arguments = {}
                if (
                    not isinstance(task_defaults, Mapping)
                    or not all(isinstance(key, str) for key in task_defaults)
                    or not isinstance(task_arguments, Mapping)
                ):
                    yield _finding(
                        check_id, Severity.ERROR,
                        f'Task {task_name!r} has malformed arguments.',
                        path=_relative_path(context.dataset_root, analysis_path), details=details,
                    )
                    continue
                default_component = TaskComponent.from_dict(
                    task=task_name, arguments=dict(task_defaults),
                )
                if isinstance(default_component, InvalidComponent):
                    yield _finding(
                        check_id, Severity.ERROR,
                        f'Task {task_name!r} has malformed default arguments.',
                        path=_relative_path(context.dataset_root, analysis_path), details=details,
                    )
                    continue
                task_defaults = default_component.arguments
                arguments.update(task_defaults)
                arguments.update(task_arguments)
                output, error = _normalize_output_path(template, arguments)
                if error is not None:
                    yield _finding(
                        check_id, Severity.ERROR,
                        f'Cannot expand output template for task {task_name!r}: {error}.',
                        path=_relative_path(context.dataset_root, analysis_path),
                        details={**details, 'template': template},
                    )
                    continue
                details['output_path'] = output.as_posix()
                claims[output].append(details)
                if output.as_posix() not in available_directories and output.as_posix() not in available_files:
                    yield _finding(
                        check_id, Severity.ERROR, 'Expected task output is absent.',
                        path=Path(output.as_posix()), details=details,
                    )

    for output, output_claims in sorted(claims.items(), key=lambda item: item[0].as_posix()):
        analysis_paths = {claim['analysis_path'] for claim in output_claims}
        if len(analysis_paths) > 1:
            yield _finding(
                check_id, Severity.ERROR, 'Separate analysis files claim the same output.',
                path=Path(output.as_posix()),
                details={'output_path': output.as_posix(), 'claims': tuple(output_claims)},
            )
    claimed = set(claims)
    for output in _recognized_outputs(context, templates):
        if output not in claimed:
            yield _finding(
                check_id, Severity.ERROR,
                'Recognized EOS output object is not claimed by any analysis step.',
                path=Path(output.as_posix()), details={'output_path': output.as_posix()},
            )


def register_dataset_checks(factory: CheckFactory) -> None:
    checks = (
        ('zenodo-metadata.completeness', ('zenodo-metadata.basic',), check_zenodo_completeness),
        ('figures.readme', ('analysis-files.selection',), check_readme_figures),
        ('figures.analysis-declarations', ('analysis-metadata.parse',), check_analysis_figures),
        ('dataset.size', ('analysis-files.selection',), check_dataset_sizes),
        ('outputs.reproducibility', ('analysis-metadata.parse',), check_reproducible_outputs),
    )
    for name, needs, function in checks:
        factory.declare(name, scope=CheckScope.COMPLETE, needs=needs)
        factory.register(name, function)
