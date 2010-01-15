#!/usr/bin/env bash
# vim: set sw=4 sts=4 et tw=80 :

run() {
    echo ">>> $@"
    if ! $@ ; then
        echo "oops!" 1>&2
        exit 127
    fi
}

get() {
    type ${1}-${2}    &>/dev/null && echo ${1}-${2}    && return
    type ${1}${2//.}  &>/dev/null && echo ${1}${2//.}  && return
    type ${1}         &>/dev/null && echo ${1}         && return
    echo "Could not find ${1} ${2}" 1>&2
    exit 127
}

run mkdir -p config
run $(get libtoolize 1.5 ) --copy --force --automake
rm -f config.cache
run $(get aclocal 1.9 )
run $(get autoheader 2.59 )
run $(get autoconf 2.59 )
run $(get automake 1.9 ) -a --copy
