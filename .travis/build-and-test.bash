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

export OS=${1}
shift 1

export CXX=${1}
shift 1

echo "==========="
${CXX} --version
echo "==========="

if [[ "xenial" == ${OS} ]] || [[ "artful" == ${OS} ]]; then
    build_and_test_ubuntu $@
fi
