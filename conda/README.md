# Conda packaging for EOS

This directory contains a [conda-build](https://docs.conda.io/projects/conda-build/)
recipe that builds, tests, and installs EOS from its autotools sources. The
autotools build system remains the authorative installation system; the recipe is a thin
driver around `autogen.bash` / `configure` / `make` / `make install`.

## Building locally

`conda/recipe/meta.yaml` is generated from `conda/recipe/meta.yaml.in` by
`configure`, which substitutes the package version from `configure.ac`
(`@PACKAGE_VERSION@`). Generate it once before building:

```bash
./autogen.bash
./configure            # writes conda/recipe/meta.yaml
```

Then, from the repository root, with `conda-build` (or `boa`) available:

```bash
conda build conda/recipe
```

To build against a specific Python version:

```bash
conda build conda/recipe --python 3.12
```

## What the recipe does

- **`meta.yaml`** declares the build/host/run dependencies. The C/C++
  dependencies (Boost, GSL, yaml-cpp, FFTW3) are located via `pkg-config`; the
  Python runtime dependencies mirror `python/setup.py.in`.
- **`build.sh`** regenerates the build system, configures with `--prefix=$PREFIX`,
  builds, runs the full `make check` test suite in the build tree, and installs.
  It deliberately avoids the PyPI-wheel `eoshep-before` target and its
  `chrpath`/`$ORIGIN` surgery -- conda-build relocates RPATHs itself.
- **`run_test.sh`** smoke-tests the relocated, installed package (imports `eos`
  and runs the CLI entry points) in conda's clean test environment.

## Source selection

By default the recipe builds the committed state of the current checkout
(`git_url: ../..`, `git_rev: HEAD`). For a released build, edit the `source:`
section of `meta.yaml` to point at a tagged sdist `url` + `sha256`.
