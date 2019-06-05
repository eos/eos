#/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

#############################
# Ubuntu specific functions #
#############################

function ubuntu_coverage() {
    [[ -z ${COVERALLS_TOKEN} ]] && exit 0
    pushd /src
    ./autogen.bash
    popd
    pushd /build
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
    coveralls-lcov --source-encoding=utf-8 --repo-token ${COVERALLS_TOKEN} /tmp/coverage.info
    popd
}

function ubuntu_deploy() {
    pushd /src
    export EOS_VERSION=${TAG:-$(git describe --abbrev=0 --tags)}
    export EOS_VERSION=${EOS_VERSION#v}
    ./autogen.bash
    popd
    pushd /build
    export CXXFLAGS="-O2 -g -march=x86-64" # build for generic x86-64 architecture
    CONFIGURE_FLAGS="--enable-pmc --enable-python --prefix=/usr"
    /src/configure ${CONFIGURE_FLAGS}
    make -j2 all
    make -j2 check VERBOSE=1
    make install
    export PYTHONPATH+=":$(make print-pythondir)"
    make -C /src/manual/examples examples
    echo Building debian package for ${OS}
    export DESTDIR=/tmp/eos-${EOS_VERSION}
    make deb DESTDIR=${DESTDIR} OS=${OS}
    dpkg -i /tmp/eos-${EOS_VERSION}.deb
    if [[ -n ${TAG} ]] && [[ "g++" == ${CXX} ]]; then
        echo Deploying debian package to package cloud ...
        package_cloud push eos/eos/ubuntu/${OS} /tmp/eos-${EOS_VERSION}.deb
    fi
    popd
}

function ubuntu_distcheck() {
    pushd /src
    export EOS_VERSION=${TAG:-$(git describe --abbrev=0 --tags)}
    export EOS_VERSION=${EOS_VERSION#v}
    ./autogen.bash
    popd
    pushd /build
    export CXXFLAGS="-O2 -g -march=x86-64" # build for generic x86-64 architecture
    CONFIGURE_FLAGS="--enable-pmc --enable-python --prefix=/usr"
    /src/configure ${CONFIGURE_FLAGS}
    make distcheck -j2 DISTCHECK_CONFIGURE_FLAGS="${CONFIGURE_FLAGS}" VERBOSE=1
    popd
}

##########################
# OSX specific functions #
##########################
function osx_build_and_test() {
    SUFFIX=$(python3 -c "import sys; print('{0}{1}'.format(sys.version_info[0], sys.version_info[1]))")
    echo using boost-python suffix ${SUFFIX}
    ./autogen.bash
    ./configure \
        --enable-pmc \
        --enable-python \
        --with-boost-python-suffix=${SUFFIX} \
        --prefix=/usr/local
    make all -j2
    make install
    make check -j2 VERBOSE=1
    export PYTHONPATH+=":$(make print-pythondir)"
    make -C manual/examples examples
}

export OS=${1}
shift 1

export CXX=${1}
shift 1

export TAG=${1}
shift 1

export USE=${1}
shift 1

echo "======================="
${CXX} --version
echo "Running with arguments:"
echo "    OS  = ${OS}"
echo "    CXX = ${CXX}"
echo "    TAG = ${TAG}"
echo "    USE = ${USE}"
echo "======================="

[[ -n ${COVERALLS_TOKEN} ]] || echo 'Skipping coverage report since $COVERALLS_TOKEN is empty'
echo "==========="

[[ -n ${PACKAGECLOUD_TOKEN} ]] || echo 'Skipping packagecloud.io deployment since $PACKAGECLOUD_TOKEN is empty'
echo "==========="

if [[ "xenial" == ${OS} ]] || [[ "bionic" == ${OS} ]]; then
    case ${USE} in
        "distcheck")
            ubuntu_distcheck $@
            ;;
        "coverage")
            ubuntu_coverage $@
            ;;
        "deploy")
            ubuntu_deploy $@
            ;;
        *)
            echo "ERROR: Unknown USE in Ubuntu tests!"
            ;;
    esac
elif [[ "osx" == ${OS} ]] ; then
    osx_build_and_test $@
fi

