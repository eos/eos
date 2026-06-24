#!/usr/bin/env bash
set -euxo pipefail

# Smoke-test the relocated, installed package in conda's clean test environment.
# (The exhaustive suite has already run via ``make check`` at build time.)

# The compiled _eos extension, the eos Python package and the bundled data
# must all import cleanly from ${PREFIX}.
python -c "import eos; print('eos', eos.__version__)"

# The command-line entry points must be on PATH and runnable. ``--help`` avoids
# emitting the full (very large) observable/reference listings.
eos-list-observables --help
eos-list-references  --help
eos-plot             --help
