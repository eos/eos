#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

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

from dataclasses import dataclass
import enum


class Severity(enum.Enum):
    ERROR = 'error'
    WARNING = 'warning'
    INFO = 'info'


@dataclass(frozen=True)
class Diagnostic:
    path: tuple
    severity: Severity
    message: str

    def prefixed(self, *segments) -> 'Diagnostic':
        return Diagnostic(tuple(segments) + self.path, self.severity, self.message)

    def location(self) -> str:
        result = ''
        for segment in self.path:
            if isinstance(segment, int):
                result += f'[{segment}]'
            else:
                if result:
                    result += '/'
                result += segment

        return result

    def _sort_key(self):
        path = tuple((0, segment) if isinstance(segment, str) else (1, segment) for segment in self.path)
        return path, self.severity.value

    def __lt__(self, other):
        if not isinstance(other, Diagnostic):
            return NotImplemented

        return self._sort_key() < other._sort_key()

    def __str__(self) -> str:
        return f'{self.location()}: {self.message}'


def _check_qualified(context, value, kind, path):
    import eos

    try:
        qn = eos.QualifiedName(value)
    except RuntimeError as e:
        yield Diagnostic(path, Severity.ERROR, f"'{value}' is not a valid qualified name: {e}")
        return
    if not context.lookup(kind, qn):
        yield Diagnostic(path, Severity.ERROR, f"{kind} '{value}' is unknown to EOS")
