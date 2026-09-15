# Copyright (c) 2026 Lorenz Gaertner
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

import os
from pathlib import Path

import eos
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ResourceError, ToolError

mcp = MCPServer('eos')


def _read_doc(relative_path: str) -> str:
    source_dir = Path(os.environ.get('EOS_SOURCE_DIR', Path.cwd()))
    path = source_dir / 'doc' / 'reference' / relative_path
    try:
        return path.read_text()
    except OSError as e:
        raise ResourceError(
            f"Could not read {path}: {e}. Set the EOS_SOURCE_DIR environment "
            "variable to the root of an eos checkout, or run eos-mcp-server "
            "from within one."
        ) from e


def _filter_kwargs(prefix, name, suffix):
    kwargs = {}
    if prefix is not None:
        kwargs['prefix'] = prefix
    if name is not None:
        kwargs['name'] = name
    if suffix is not None:
        kwargs['suffix'] = suffix
    return kwargs


def _describe_observable(entry):
    return {
        'latex': entry.latex(),
        'unit': str(entry.unit()),
        'kinematic_variables': [str(kv) for kv in entry.kinematic_variables()],
    }


@mcp.tool()
def observables(prefix: str | None = None, name: str | None = None,
                suffix: str | None = None, qualified_name: str | None = None) -> dict:
    "List EOS observables matching the given filters, or describe a single one by qualified name."

    if qualified_name is not None:
        try:
            entry = eos.Observables()[qualified_name]
        except RuntimeError as e:
            raise ToolError(str(e)) from e
        return {qualified_name: _describe_observable(entry)}

    container = eos.Observables(**_filter_kwargs(prefix, name, suffix))
    return {
        str(qn): _describe_observable(entry)
        for qn, entry in container
        if container.filter_entry(qn)
    }


def _describe_parameter(p):
    return {
        'latex': p.latex(),
        'unit': str(p.unit()),
        'central': p.central(),
        'min': p.min(),
        'max': p.max(),
    }


@mcp.tool()
def parameters(prefix: str | None = None, name: str | None = None,
               suffix: str | None = None, qualified_name: str | None = None) -> dict:
    "List EOS parameters matching the given filters, or describe a single one by qualified name."

    if qualified_name is not None:
        try:
            p = eos.Parameters()[qualified_name]
        except RuntimeError as e:
            raise ToolError(str(e)) from e
        return {qualified_name: _describe_parameter(p)}

    container = eos.Parameters(**_filter_kwargs(prefix, name, suffix))
    return {
        str(p.name()): _describe_parameter(p)
        for p in container
        if container.filter_entry(p.name())
    }


def _describe_constraint(entry):
    return {
        'type': str(entry.type()),
        'observables': [str(o) for o in entry.observables()],
        'references': [str(r) for r in entry.references()],
    }


@mcp.tool()
def constraints(prefix: str | None = None, name: str | None = None,
                suffix: str | None = None, qualified_name: str | None = None) -> dict:
    "List EOS constraints matching the given filters, or describe a single one by qualified name."

    if qualified_name is not None:
        try:
            entry = eos.Constraints()[qualified_name]
        except RuntimeError as e:
            raise ToolError(str(e)) from e
        return {qualified_name: _describe_constraint(entry)}

    container = eos.Constraints(**_filter_kwargs(prefix, name, suffix))
    return {
        str(qn): _describe_constraint(entry)
        for qn, entry in container
        if container.filter_entry(qn)
    }


def _describe_reference(entry):
    return {
        'title': entry.title(),
        'authors': entry.authors(),
        'eprint_archive': entry.eprint_archive(),
        'eprint_id': entry.eprint_id(),
    }


@mcp.tool()
def references(year: str | None = None, index: str | None = None, key: str | None = None) -> dict:
    "List EOS references matching the given filters, or describe a single one by key."

    if key is not None:
        container = eos.References()
        if key not in container:
            raise ToolError(f"Unknown reference: '{key}' not known")
        return {key: _describe_reference(container[key])}

    container = eos.References(year=year, index=index)
    return {
        str(rn): _describe_reference(entry)
        for rn, entry in container
        if container.filter_entry(rn)
    }


@mcp.resource('eos://doc/cli-reference')
def cli_reference_doc() -> str:
    "The eos command-line interface reference documentation."
    return _read_doc('command-line-interface.rst')


@mcp.resource('eos://doc/defining-observables')
def defining_observables_doc() -> str:
    "How to define custom observables in an analysis file."
    return _read_doc('defining-observables.rst')
