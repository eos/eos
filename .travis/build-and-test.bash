#/bin/bash

function build_and_test_ubuntu() {
    pushd /src
    ./autogen.bash || exit 1
    popd
    pushd /build
    /src/configure \
        --enable-pmc \
        --enable-python \
        --prefix=/usr \
        || exit 1
    make all || exit 1
    make check VERBOSE=1 || exit 1
    popd
}

name=${1}
shift 1

if [[ "xenial" == ${name} ]] ; then
    build_and_test_ubuntu $@
fi
