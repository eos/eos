**********
MCP Server
**********

.. _mcp-server:

.. note::

   The EOS MCP server is completely optional and does not provide any functionality beyond the interactive Python interface.

``eos-mcp-server`` exposes EOS's observables, parameters, constraints, and references as read-only
lookup tools for an `MCP <https://modelcontextprotocol.io/>`_ client (Claude Code, Claude Desktop,
Codex, ...), so an LLM coding agent can look up correct qualified names, units, and kinematics
instead of guessing. It runs locally over stdio; it performs no computation and never modifies
EOS state.

Installing
==========

After building EOS with ``./configure --enable-python --prefix=<venv>`` and ``make install`` (see
:ref:`installation-from-source`), install the ``mcp`` SDK into the same environment::

    <venv>/bin/pip install mcp

``<venv>/bin/eos-mcp-server`` is then ready to use.

Connecting a client
====================

Point the client at ``<venv>/bin/eos-mcp-server``, with ``EOS_SOURCE_DIR`` set to your checkout
root (needed only to serve the two doc resources below; the lookup tools work without it).

Claude Code (project ``.mcp.json``, or ``claude mcp add``)::

    { "mcpServers": { "eos": {
      "command": "<venv>/bin/eos-mcp-server",
      "env": { "EOS_SOURCE_DIR": "/path/to/eos" }
    } } }

Codex CLI (``~/.codex/config.toml``, or ``codex mcp add``)::

    [mcp_servers.eos]
    command = "<venv>/bin/eos-mcp-server"
    env = { EOS_SOURCE_DIR = "/path/to/eos" }

Tools and resources
====================

  * ``observables``, ``parameters``, ``constraints`` -- list entries matching a ``prefix``/``name``/``suffix``
    filter, or describe a single one given its ``qualified_name``.
  * ``references`` -- list entries matching a ``year``/``index`` filter, or describe a single one given its ``key``.
  * ``eos://doc/cli-reference``, ``eos://doc/defining-observables`` -- the corresponding reference documents.
