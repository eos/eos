# AGENTS.md

Working notes for agents operating on this EOS checkout. Scope: build, test, and
contribution conventions. Physics/analysis notes are deliberately not included here.

## Local instructions

If `.agents-local/AGENTS.md` or `.agents-local/CLAUDE.md` exist, read them at the start of
every session and follow them in addition to this file. They are personal to each
developer and never committed. If they conflict with this file, stop and ask the user
which rule applies.

## Building the C++ code

EOS is written in C++20 (`configure.ac` sets `-std=c++20`) and requires g++ 10.1 or
higher (see `doc/installation.rst.in`). Some g++ installations that otherwise claim
C++20 support still lack `<format>`, which EOS uses widely; if a bare `make` fails with
`fatal error: format: No such file or directory`, select or install a newer g++ (e.g. a
`gcc-toolset`/`devtoolset` on RHEL-family systems, or an updated Homebrew/apt package
elsewhere) that fully implements it, and put it on `PATH` ahead of the default compiler.

After editing a `Makefile.am`, regenerate with `automake eos/<subdir>/Makefile`, then
`./config.status eos/<subdir>/Makefile`, then build. `Makefile.in` is generated
(`MAINTAINERCLEANFILES`) and not tracked.

## Running tests

Run tests from the repository root through the automake `check` target with `TESTS=`
set to the test binary name (no path) — this gets the `.trs`/`.log` generation and the
proper environment:

```
make -C eos/maths check VERBOSE=1 TESTS="dft-plan_TEST"
make -C python    check VERBOSE=1 TESTS="eos/data/markov_chain_TEST.py"
```

Add `-j<N>` (e.g. `-j8`) to parallelize the build step of a slow `check` invocation:

```
make -C eos/form-factors -j8 VERBOSE=1 TESTS="parametric-bsz2015_TEST" check
```

Never invoke `./<name>_TEST` directly. The `check` target sets environment variables
that the tests require. Moreover, C++ test binaries are libtool wrapper scripts whose
rpath includes the configured installation prefix; if a prior `make install` populated
that prefix, it can **shadow** the freshly built in-tree `.libs/*.so`, so source edits
silently appear to have no effect.

## Python

Build and use a Python virtual environment with EOS installed via `./configure
--enable-python` (see `doc/installation.rst.in`) for anything importing this repo's
modules (the compiled `eos` extension, `eos.figure`, tests, coverage). Do not `pip
install` this project ad hoc into a different interpreter or mix in user-site packages
alongside a venv install — mixing installs of the compiled `_eos` extension across
locations triggers a numpy "cannot load module more than once per process"
`ImportError`. Install extra tooling (e.g. `coverage`) into that same venv.

Register new Python modules in **both** `EXTRA_DIST` and the relevant `*_SCRIPTS` list
in `python/Makefile.am`.

## Verbosity

Keep changes, commit messages, and code comments concise. Prefer the smallest diff
that solves the task, a short commit message (see Commits below), and few or no
comments (see Code style below) — do not pad any of these.

## GitHub

`gh` may not be on the default `PATH` on a given checkout; locate it (`which gh`, or
check an active project virtual environment's `bin/`) before relying on it.

Treat GitHub access as read-only by default. Use `gh`/`git fetch`/`pull` only when the
task requires it — inspecting a branch or commit range that already exists locally does
not, so work from local refs instead of fetching first. Never `push`, open a GitHub
issue or pull request, or modify an existing one (editing, commenting, labelling,
closing, reopening, merging, reviewing) without the user's explicit approval for that
specific action.

Before opening a pull request, confirm that the build succeeds and all tests pass
at the tip of the branch. Never build with all available cores; keep `-j<N>` below
`nproc`.

## Commits

- **Copyright headers:** before every commit, check each modified file's existing
  `Copyright (c) <BEGIN>-<END> <Name>` header lines (e.g. `eos/observable.cc`) and, for
  the line matching the change's author, extend `<END>` to the current year, using
  range notation (`2023-2026`), not a comma-separated list. Only extend a line that
  already names that author — files with no header, or no line for that author, are
  left alone; this is about updating existing headers, not adding new ones. Fold the
  header fix into the same commit as the change.
- **No co-author trailer:** never add `Co-Authored-By: Claude ...` or any other Claude
  co-authorship acknowledgement.
- Follow the repo's subject-prefix convention, e.g. `[eos]`, `[python]`, `[infra]` —
  match the subdirectory/module being changed (see `git log --oneline` for existing
  prefixes).
- **Default to a one-line commit message** (title only). Add a body only in
  extraordinary circumstances: to warn about a trap for future changes, to flag a break
  with forward compatibility, or similar. Keep bodies to a minimum — typically a single
  paragraph.

## Code style

New file names use hyphens, never underscores (`observable-cache_TEST.cc`, not
`observable_cache_TEST.cc`); underscored names are legacy and are being converted, so never
rename in that direction. Python modules are exempt, since an importable module name cannot
contain a hyphen; their `_TEST.py` scripts follow the module's name.

`clang-format` is the canonical formatter. In `.cc` files, it correctly emits one
leading space before an opening `namespace eos {` — do not flag that leading space as a
style issue in review.

**Format through `pre-commit`, not through a local `clang-format` binary.** The canonical
version is pinned in `.pre-commit-config.yaml`; the `clang-format`
on `PATH` may well be newer and format some constructs differently, rewriting lines that
are already correct — which looks like pre-existing drift but is not. Use
`pre-commit run clang-format --files <path>`. If a newer binary has already been run,
discard the lines it touched outside the change at hand.

Comment sparingly — code says *what*, comments say *why*. Keep in-code comments to a
minimum. When a comment is necessary, aim for a single line. In exceptional cases, use
two lines; comments must never exceed four lines.
**Do not use** the following software development/engineering jargon: boundary, materialize,
provenance, 'single truth' (or 'sole truth').
**Do not refer to development phases or steps** when writing comments.

## References and bibliography

The bibliography is generated automatically from `eos/references.yaml`; do not
hand-maintain bibliography output. Each entry carries an `inspire-id` (e.g.
`Asatryan:2001zw`) alongside `authors`, `title`, and the `eprint` archive/id.

An agent must **never** guess or invent an INSPIRE ID. If a new reference is needed and
the ID has not been provided, ask the user for it and wait — do not fabricate a
plausible-looking `Author:YYnnn` key.

## Changelog

All modifications to `Changelog.md` are collected in a single commit at the **tip of the
branch**, separate from the code commits. Its message reads:

```
[eos] Update changelog for PR #ABCD
```

The PR number `ABCD` is to be provided by the user — do not guess it. If it has not been
given and the branch has no PR open yet, determine `N`, the largest number among all
open **and** closed PRs, issues, and discussions in the repository, and use `N+1` (the
number GitHub will assign to the next-created item of any of these three kinds, since
they share one counter). Use `gh` for the lookup, e.g.:

```
gh pr list    --repo eos/eos --state all --limit 1 --json number --jq '.[0].number'
gh issue list --repo eos/eos --state all --limit 1 --json number --jq '.[0].number'
gh api graphql -f query='{ repository(owner:"eos", name:"eos") { discussions(first:1, orderBy:{field:CREATED_AT, direction:DESC}) { nodes { number } } } }'
```

If the user has already provided a PR number, use that instead and skip this lookup.

## Ancillary files for issues, PRs, and features

Scratch and working material produced while addressing an item goes into a
per-item directory at the repository root, named by its number alone (no `#`) or, for
an untracked feature, by its name:

- working on issue #ABCD → `issues/ABCD/`
- resolving problems with pull request #ABCD → `prs/ABCD/`
- working on an untracked feature `foo` → `features/foo/`

Put notes, reproducers, logs, scripts, and intermediate data there rather than at the
repository root or in `/tmp`. These directories are untracked working space; nothing in
them belongs in a commit.
