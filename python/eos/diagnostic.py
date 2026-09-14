#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2026 Danny van Dyk
# Copyright (c) 2026 Lorenz Gärtner
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

from dataclasses import dataclass, replace
import enum
import yaml

_STR_TAG = 'tag:yaml.org,2002:str'


class Severity(enum.Enum):
    ERROR = 'error'
    WARNING = 'warning'
    INFO = 'info'


@dataclass(frozen=True)
class Diagnostic:
    path: tuple
    severity: Severity
    message: str
    # The 1-based line number in the source YAML file that 'path' resolves to, or None if it has
    # not been resolved.
    line: int | None = None

    def prefixed(self, *segments) -> 'Diagnostic':
        return Diagnostic(tuple(segments) + self.path, self.severity, self.message, self.line)

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
        if self.line is not None:
            return f'line {self.line}: {self.location()}: {self.message}'
        return f'{self.location()}: {self.message}'


# The one _segments() call site (analysis_file_description.py:1651) that overrides the default
# identifier='name' with 'id', keyed by the path of the sequence it applies to.
_SEQUENCE_IDENTIFIER_OVERRIDES = {('steps',): 'id'}

# Sections whose YAML value is a mapping (name -> description) that analysis_file_description.py
# converts to a [{'name': name, **description}, ...] list before calling _segments() on it
# (observables, parameters: analysis_file_description.py:1633-1646; a likelihood's
# manual_constraints: analysis_file_description.py:604-608), keyed by their own last path segment.
# Each entry's path segment is therefore its key when that is a clean string, or otherwise its
# position among its siblings -- the same rule _sequence_child_segment applies to actual sequences.
_NAME_KEYED_MAPPING_SECTIONS = {'observables', 'parameters', 'manual_constraints'}


def _is_clean_identifier(node):
    return (
        node.tag == _STR_TAG
        and '/' not in node.value
        and not any(character.isspace() for character in node.value)
    )


def _sequence_child_segment(child, index, identifier):
    # Mirror eos.analysis_file_description._segments(): use the child's 'identifier' field as the
    # path segment when it is a clean string (no '/', no whitespace); otherwise fall back to 'index'.
    if isinstance(child, yaml.MappingNode):
        for key_node, value_node in child.value:
            if key_node.tag == _STR_TAG and key_node.value == identifier and _is_clean_identifier(value_node):
                return value_node.value
    return index


def _named_mapping_child_segment(key_node, index):
    # Mirror the {'name': key, **description} / _segments() construction used for the sections in
    # _NAME_KEYED_MAPPING_SECTIONS: a clean string key is used verbatim, otherwise the fallback is
    # the entry's position among its siblings (enumeration order matches dict/YAML source order).
    return key_node.value if _is_clean_identifier(key_node) else index


def _walk_yaml_nodes(node, path, index):
    if isinstance(node, yaml.MappingNode):
        name_keyed = bool(path) and path[-1] in _NAME_KEYED_MAPPING_SECTIONS
        for i, (key_node, value_node) in enumerate(node.value):
            if key_node.tag != _STR_TAG:
                continue
            segment = _named_mapping_child_segment(key_node, i) if name_keyed else key_node.value
            child_path = path + (segment,)
            index[child_path] = key_node.start_mark.line + 1
            _walk_yaml_nodes(value_node, child_path, index)
    elif isinstance(node, yaml.SequenceNode):
        identifier = _SEQUENCE_IDENTIFIER_OVERRIDES.get(path, 'name')
        for i, child in enumerate(node.value):
            child_path = path + (_sequence_child_segment(child, i, identifier),)
            index[child_path] = child.start_mark.line + 1
            _walk_yaml_nodes(child, child_path, index)


def _build_line_index(path: str) -> dict:
    """Map each eos.diagnostic.Diagnostic path that may occur for the YAML file at 'path' to its
    1-based source line, by walking the raw parse tree (yaml.compose() retains source positions;
    yaml.safe_load()'s plain dict/list result does not)."""
    with open(path, encoding='utf-8') as stream:
        root = yaml.compose(stream)

    index = {(): root.start_mark.line + 1 if root is not None else 1}
    if root is not None:
        _walk_yaml_nodes(root, (), index)
    return index


def _line_for_path(index: dict, path: tuple) -> int:
    # Fall back to the closest ancestor path present in the index, e.g. for a missing-key
    # diagnostic whose path does not itself appear in the YAML source.
    for depth in range(len(path), -1, -1):
        prefix = path[:depth]
        if prefix in index:
            return index[prefix]
    return index[()]


def attach_line_numbers(diagnostics, path: str) -> list:
    """Return 'diagnostics' with each entry's 'line' resolved against the YAML source at 'path'."""
    diagnostics = list(diagnostics)
    if not diagnostics:
        # Avoid re-parsing the whole YAML document just to build an index nothing will look up.
        return diagnostics
    index = _build_line_index(path)
    return [replace(diagnostic, line=_line_for_path(index, diagnostic.path)) for diagnostic in diagnostics]


def _check_qualified(context, value, kind, path):
    import eos

    try:
        qn = eos.QualifiedName(value)
    except RuntimeError as e:
        yield Diagnostic(path, Severity.ERROR, f"'{value}' is not a valid qualified name: {e}")
        return
    if not context.lookup(kind, qn):
        yield Diagnostic(path, Severity.ERROR, f"{kind} '{value}' is unknown to EOS")
