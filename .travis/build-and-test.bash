#/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

function build_and_test_ubuntu() {
    pushd /src
    ./autogen.bash
    popd
    pushd /build
    export CXXFLAGS="-O2 -g -march=x86-64" # build for generic x86-64 architecture
    CONFIGURE_FLAGS="--enable-pmc --enable-python --prefix=/usr"
    /src/configure ${CONFIGURE_FLAGS}
    make distcheck -j2 DISTCHECK_CONFIGURE_FLAGS="${CONFIGURE_FLAGS}" VERBOSE=1
    make install
    export PYTHONPATH+=":$(make print-pythondir)"
    make -C /src/manual/examples examples
    if [[ -n ${TAG} ]] && [[ "g++" == ${CXX} ]]; then
        echo Building debian package for ${OS}
        export DESTDIR=/tmp/eos-${TAG#v}
        make deb DESTDIR=${DESTDIR} OS=${OS}
        package_cloud push eos/eos/ubuntu/${OS} /tmp/eos-${TAG#v}.deb
    fi
    popd
}

function build_and_coverage_ubuntu() {
    pushd /build
    make distclean
    export CXXFLAGS="-O2 -g --coverage"
    /src/configure \
        --enable-pmc \
        --enable-python \
        --prefix=/usr
    make all -j2
    make check -j2
    lcov --directory . --capture --output-file /tmp/coverage.info
    popd
    lcov \
        --remove /tmp/coverage.info '/usr/*' '/test/*' \
        --output-file /tmp/coverage.info
    pushd /src
    coveralls-lcov --repo-token ${COVERALLS_TOKEN} /tmp/coverage.info
    popd
}

export OS=${1}
shift 1

export CXX=${1}
shift 1

export TAG=${1}
shift 1

echo "==========="
${CXX} --version

[[ -n ${COVERALLS_TOKEN} ]] || echo 'Skipping coverage report since $COVERALLS_TOKEN is empty'
echo "==========="

[[ -n ${PACKAGECLOUD_TOKEN} ]] || echo 'Skipping packagecloud.io deployment since $PACKAGECLOUD_TOKEN is empty'
echo "==========="

if [[ "xenial" == ${OS} ]] || [[ "bionic" == ${OS} ]]; then
    build_and_test_ubuntu $@
fi

if [[ -n ${COVERALLS_TOKEN} ]] && [[ "xenial" == ${OS} ]] && [[ "g++" == ${CXX} ]] ; then
    build_and_coverage_ubuntu $@
fi
