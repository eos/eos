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

from collections.abc import Callable, Iterable, Mapping, MutableMapping
from dataclasses import dataclass, field
from enum import Enum, IntEnum
from pathlib import Path
from types import MappingProxyType
from typing import Any, Protocol, TextIO


class Severity(Enum):
    ERROR = 'error'
    WARNING = 'warning'
    INFO = 'info'


class ExitStatus(IntEnum):
    SUCCESS = 0
    CHECK_ERRORS = 1
    USAGE_OR_INFRASTRUCTURE = 2


class CheckScope(IntEnum):
    BASIC = 1


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
    cache: MutableMapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.dataset_root = Path(self.dataset_root)
        self.analysis_paths = tuple(Path(path) for path in self.analysis_paths)
        if self.main_analysis_path is not None:
            self.main_analysis_path = Path(self.main_analysis_path)

    def cached(self, key: str, loader: Callable[[], Any]) -> Any:
        if key not in self.cache:
            self.cache[key] = loader()
        return self.cache[key]


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
        return tuple(
            self._checks[name]
            for name, declaration in self._declarations.items()
            if declaration.scope <= scope and name in self._checks
        )

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
        try:
            for registered in self._checks:
                produced = registered.function(context)
                for finding in produced:
                    if not isinstance(finding, Finding):
                        raise TypeError(
                            f"check '{registered.declaration.name}' returned a non-Finding value"
                        )
                    findings.append(finding)
        except SystemExit as error:
            return CheckResult(tuple(findings), f'a check attempted to exit with status {error.code}')
        except Exception as error:
            return CheckResult(tuple(findings), f'{type(error).__name__}: {error}')

        return CheckResult(tuple(findings))


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
