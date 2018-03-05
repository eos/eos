#/bin/bash

function build_and_test_ubuntu() {
    pushd /src
    ./autogen.bash || exit 1
    popd
    pushd /build
    export CXXFLAGS="-O2 -g"
    /src/configure \
        --enable-pmc \
        --enable-python \
        --prefix=/usr \
        || exit 1
    make distcheck -j2 DISTCHECK_CONFIGURE_FLAGS="--enable-pmc --enable-python --prefix=/usr" VERBOSE=1 || exit 1
    make install || exit 1
    make -C /src/manual/examples examples || exit 1
    popd
}

function build_and_coverage_ubuntu() {
    pushd /build
    make distclean
    export CXXFLAGS="-O2 -g --coverage"
    /src/configure \
        --enable-pmc \
        --enable-python \
        --prefix=/usr \
        || exit 1
    make all -j2 || exit 1
    make check -j2 || exit 1
    lcov --directory . --capture --output-file /tmp/coverage.info || exit 1
    popd
    lcov \
        --remove /tmp/coverage.info '/usr/*' '/test/*' \
        --output-file /tmp/coverage.info || exit 1
    pushd /src
    coveralls-lcov --repo-token ${COVERALLS_TOKEN} /tmp/coverage.info || exit 1
    popd
}

export OS=${1}
shift 1

export CXX=${1}
shift 1

echo "==========="
${CXX} --version

[[ -n ${COVERALLS_TOKEN} ]] || echo 'Skipping coverage report since $COVERALLS_TOKEN is empty'
echo "==========="

if [[ "xenial" == ${OS} ]] || [[ "artful" == ${OS} ]]; then
    build_and_test_ubuntu $@
fi

if [[ -n ${COVERALLS_TOKEN} ]] && [[ "xenial" == ${OS} ]] && [[ "g++" == ${CXX} ]] ; then
    build_and_coverage_ubuntu $@
fi
