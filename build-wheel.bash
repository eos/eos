#!/usr/bin/env bash
# vim: set sw=4 sts=4 et tw=80 :
#
# Build an EOS wheel locally, on Linux or on macOS (arm64), using exactly the
# cibuildwheel settings that .github/workflows/pypi-build+check+deploy.yaml
# uses in CI -- the matrix row, the CIBW_* environment, and the before-build/
# repair/test commands are all extracted from that workflow file at run time
# rather than copied into this script by hand, so the two cannot silently
# drift apart. See doc/installation.rst.in, section "Building a Wheel
# Manually", for prerequisites and a walkthrough.
#
# Usage:
#   ./build-wheel.bash [--python cp312,cp313,cp314] [--arch ARCH]
#                       [--quick] [--output-dir DIR] [--print-env] [-h]
#
# Two constraints inherited from the setuptools/PEP 517 build (see
# python/setup.py.in and configure.ac):
#   - EOS_STAGING_PREFIX must be set to a 'make install' prefix; this script
#     gets that from the workflow's own CIBW_ENVIRONMENT, unmodified.
#   - AM_CXXFLAGS freezes an absolute -I<source path> at 'configure' time, so
#     the build must run against the same tree that was configured. Run this
#     script in place; do not copy or relocate the checkout mid-build.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

WORKFLOW="$ROOT/.github/workflows/pypi-build+check+deploy.yaml"
PYTHON_BIN="${EOS_BUILD_WHEEL_PYTHON:-$(command -v python3 || true)}"
CIBUILDWHEEL_VERSION="3.2.1" # matches the workflow's 'uses: pypa/cibuildwheel@v3.2.1'

case "$(uname -s)" in
    Linux)  HOST_OS=linux ;;
    Darwin) HOST_OS=macos ;;
    *) echo "error: unsupported host OS '$(uname -s)'" >&2; exit 1 ;;
esac

case "$(uname -m)" in
    x86_64)  HOST_ARCH=x86_64 ;;
    aarch64) HOST_ARCH=aarch64 ;;
    arm64)   HOST_ARCH=arm64 ;;
    *) echo "error: unsupported host architecture '$(uname -m)'" >&2; exit 1 ;;
esac

PY_VERSIONS=""
ARCH="$HOST_ARCH"
QUICK=0
OUTPUT_DIR="$ROOT/wheelhouse"
PRINT_ENV=0

usage() {
    cat <<USAGE
Usage: $(basename "$0") [options]

Build an EOS wheel locally with the same cibuildwheel settings CI uses.

Options:
  --python LIST     Comma-separated CPython tags to build, e.g. cp312,cp313.
                     Default: every tag the workflow ships for this OS/arch
                     (currently cp312,cp313,cp314).
  --arch ARCH       Target architecture (x86_64, aarch64 on Linux; arm64 on
                     macOS). Default: the host architecture ($HOST_ARCH).
                     A non-host Linux arch needs QEMU registered with Docker:
                       docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
  --quick           Skip 'make check' inside CIBW_BEFORE_BUILD for faster
                     iteration. Off by default, matching CI.
  --output-dir DIR  Where to place the built wheel(s). Default: ./wheelhouse
  --print-env       Print the resolved CIBW_* environment for each requested
                     version and exit, without building anything.
  -h, --help        Show this message.

Prerequisites:
  Linux: Docker or Podman (cibuildwheel builds inside the manylinux
         container referenced by the workflow -- no host toolchain needed).
  macOS: run natively on arm64 hardware, with Homebrew installed; the
         dependency bundle built by eos/docker.io is fetched via 'oras'
         (installed through Homebrew as part of the fetch step).
  Both:  the venv at /mt/home/dvan-dyk/.eos/master (or \$EOS_BUILD_WHEEL_PYTHON)
         is used to run cibuildwheel and the workflow-parsing helper;
         cibuildwheel is installed into it automatically if missing.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python) PY_VERSIONS="$2"; shift 2 ;;
        --arch) ARCH="$2"; shift 2 ;;
        --quick) QUICK=1; shift ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --print-env) PRINT_ENV=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "error: unrecognized argument '$1'" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ "$HOST_OS" == "macos" && "$ARCH" != "arm64" ]]; then
    echo "error: macOS wheels are arm64-only (see issues/setuptools+macos/DESIGN.md, decision D4); got --arch $ARCH" >&2
    exit 1
fi

if [[ -z "$PY_VERSIONS" ]]; then
    PY_VERSIONS="cp312,cp313,cp314"
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "error: python interpreter '$PYTHON_BIN' not found or not executable (set EOS_BUILD_WHEEL_PYTHON to override)" >&2
    exit 1
fi

if [[ "$HOST_OS" == "linux" ]]; then
    if ! command -v docker >/dev/null 2>&1 && ! command -v podman >/dev/null 2>&1; then
        echo "error: cibuildwheel needs Docker or Podman on Linux; neither was found on PATH" >&2
        exit 1
    fi
fi

if ! "$PYTHON_BIN" -c 'import cibuildwheel' >/dev/null 2>&1; then
    echo "==> installing cibuildwheel==$CIBUILDWHEEL_VERSION into $PYTHON_BIN's venv" >&2
    "$PYTHON_BIN" -m pip install --quiet "cibuildwheel==$CIBUILDWHEEL_VERSION"
fi
CIBUILDWHEEL="$("$PYTHON_BIN" -c 'import sysconfig, os; print(os.path.join(sysconfig.get_path("scripts"), "cibuildwheel"))')"
if [[ ! -x "$CIBUILDWHEEL" ]]; then
    echo "error: cibuildwheel not found at expected path '$CIBUILDWHEEL' after install" >&2
    exit 1
fi

# Renders the CIBW_* environment for one matrix row straight out of the
# workflow file: loads the YAML, picks the row matching (host_os, arch,
# version), substitutes the handful of '${{ matrix.* }}' /
# '${{ steps.prerelease.outputs.option }}' expressions those strings use
# (a local build is never a GitHub release, so the prerelease option is
# always empty), and fails loudly if the workflow uses an expression this
# script does not know how to substitute -- that failure IS the drift check
# called for in issues/setuptools+macos/STEP-7.md task 2.
render_env() {
    local version="$1"
    "$PYTHON_BIN" - "$WORKFLOW" "$HOST_OS" "$ARCH" "$version" "$QUICK" <<'PY'
import re
import shlex
import sys

import yaml

workflow_path, host_os, arch, version, quick = sys.argv[1:6]
quick = quick == "1"

with open(workflow_path) as f:
    workflow = yaml.safe_load(f)

job = workflow["jobs"]["build_wheels"]
rows = job["strategy"]["matrix"]["include"]

runner_prefix = "macos" if host_os == "macos" else "ubuntu"
candidates = [
    r for r in rows
    if r["arch"] == arch and r["version"] == version and r["runner"].startswith(runner_prefix)
]
if len(candidates) != 1:
    sys.exit(
        f"error: expected exactly one matrix row for arch={arch} version={version} "
        f"runner~={runner_prefix} in {workflow_path}, found {len(candidates)} -- "
        "the workflow's matrix changed shape; update build-wheel.bash's row selector"
    )
row = candidates[0]

try:
    build_step = next(
        s for s in job["steps"] if s.get("uses", "").startswith("pypa/cibuildwheel")
    )
except StopIteration:
    sys.exit(f"error: no 'pypa/cibuildwheel' step found in {workflow_path}'s build_wheels job")
cibw_env = build_step.get("env", {})

fetch_step = next(
    (s for s in job["steps"] if s.get("name") == "Fetch and unpack the macOS dependency bundle"),
    None,
)

substitutions = {
    "matrix.arch": row["arch"],
    "matrix.cxxflags": row["cxxflags"],
    "matrix.version": row["version"],
    "matrix.boost_python_suffix": str(row["boost_python_suffix"]),
    "matrix.prefix": row.get("prefix", ""),
    "matrix.deployment_target": row.get("deployment_target", ""),
    "matrix.dep_tarball_tag": row.get("dep_tarball_tag", ""),
    "matrix.lto_option": row.get("lto_option", ""),
    "steps.prerelease.outputs.option": "",
}

def substitute(text):
    def repl(m):
        key = m.group(1).strip()
        if key not in substitutions:
            sys.exit(
                f"error: workflow expression '{{{{ {key} }}}}' has no local "
                "substitution in build-wheel.bash -- the workflow changed; "
                "teach the script about the new expression"
            )
        return substitutions[key]
    rendered = re.sub(r"\$\{\{\s*(.+?)\s*\}\}", repl, text)
    if "${{" in rendered:
        sys.exit(f"error: unsubstituted GitHub Actions expression left in: {rendered!r}")
    return rendered

for key, value in cibw_env.items():
    if not isinstance(value, str):
        continue
    rendered = substitute(value)
    if quick and key == "CIBW_BEFORE_BUILD":
        new_rendered = re.sub(r"^\s*make -j4 check VERBOSE=1\s*\n", "", rendered, flags=re.MULTILINE)
        if new_rendered == rendered:
            print(
                "warning: --quick requested but no 'make -j4 check' line found "
                "to remove from CIBW_BEFORE_BUILD -- printing it unmodified",
                file=sys.stderr,
            )
        rendered = new_rendered
    print(f"export {key}={shlex.quote(rendered)}")

if host_os == "macos" and fetch_step is not None:
    print(f"export EOS_MACOS_FETCH_CMD={shlex.quote(substitute(fetch_step['run']))}")
PY
}

mkdir -p "$OUTPUT_DIR"
# The macOS repair command (CIBW_REPAIR_WHEEL_COMMAND_MACOS, rendered
# verbatim from the workflow) writes debug symbols under
# "$GITHUB_WORKSPACE/macos-debug"; GITHUB_WORKSPACE is normally an
# ambient GitHub Actions variable, so a local build must provide it itself.
export GITHUB_WORKSPACE="$ROOT"

# cibuildwheel checks for {setup.py, setup.cfg, pyproject.toml} at the package
# root before it does anything else -- before CIBW_BEFORE_BUILD ever runs --
# but pyproject.toml is generated from pyproject.toml.in by './configure'
# (D2), which is itself part of CIBW_BEFORE_BUILD. Mirrors the workflow step
# "Materialize a placeholder pyproject.toml for cibuildwheel's own package
# check" in pypi-build+check+deploy.yaml; './configure' overwrites this with
# the real, version-substituted file before 'pip wheel' reads it.
cp "$ROOT/pyproject.toml.in" "$ROOT/pyproject.toml"

IFS=',' read -r -a versions <<< "$PY_VERSIONS"
fetched_macos_bundle=0

for version in "${versions[@]}"; do
    env_file="$(mktemp)"
    render_env "$version" > "$env_file"

    if [[ "$PRINT_ENV" -eq 1 ]]; then
        echo "# --- resolved environment for $version-$ARCH ---"
        cat "$env_file"
        rm -f "$env_file"
        continue
    fi

    # shellcheck disable=SC1090
    source "$env_file"
    rm -f "$env_file"

    if [[ "$HOST_OS" == "macos" && "$fetched_macos_bundle" -eq 0 ]]; then
        echo "==> fetching macOS dependency bundle" >&2
        bash -euo pipefail -c "$EOS_MACOS_FETCH_CMD"
        fetched_macos_bundle=1
    fi

    echo "==> building $version-$ARCH" >&2
    "$CIBUILDWHEEL" --platform "$HOST_OS" --output-dir "$OUTPUT_DIR" "$ROOT"
done

if [[ "$PRINT_ENV" -eq 0 ]]; then
    echo "==> wheel(s) written to $OUTPUT_DIR" >&2
fi
