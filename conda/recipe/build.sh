#!/usr/bin/env bash
set -euxo pipefail

# Make the pkg-config (.pc) files and aclocal macros shipped by the host
# dependencies discoverable. conda's activation usually sets these, but we
# export them defensively so the build does not depend on activation order.
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export ACLOCAL_PATH="${PREFIX}/share/aclocal:${ACLOCAL_PATH:-}"

# Regenerate the build system from the autotools sources (configure et al.
# are not exported by ``git_rev``).
./autogen.bash

# conda-forge's boost ships libboost_python<major><minor>.so, which is exactly
# the suffix configure.ac now defaults to. We pass it explicitly so the build
# is self-documenting and independent of that default.
BOOST_PY_SUFFIX=$(${PYTHON} -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")

mkdir -p _build
cd _build

../configure \
    --prefix="${PREFIX}" \
    --enable-python \
    --enable-cli \
    --with-boost-python-suffix="${BOOST_PY_SUFFIX}" \
    --with-custom-pythondir="site-packages" \
    PYTHON="${PYTHON}"
# gsl, yaml-cpp and fftw are located via pkg-config (PKG_CONFIG_PATH above);
# the compiler activation packages export CXX/CXXFLAGS/LDFLAGS pointing at
# ${PREFIX}, so no --with-* path overrides are required.

make -j"${CPU_COUNT}"

# Run the full autotools test suite in the build tree (C++ TESTS + Python
# *_TEST.py). The tests rely on EOS_TESTS_* env vars that the Makefiles set
# relative to the source/build tree, so they must run here, not in conda's
# separate (installed-only) test phase.
make -j"${CPU_COUNT}" check VERBOSE=1

make install

# We intentionally do NOT invoke the python/ 'eoshep-before' target: that path
# performs chrpath/$ORIGIN surgery to vendor shared objects into a PyPI wheel.
# conda-build relocates RPATHs to ${PREFIX}/lib on its own, so a plain
# ``make install`` produces a correctly relocatable package.
