#!/bin/bash

mkdir /tmp/build-minuit2
cd /tmp/build-minuit2
wget http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz
tar zxf Minuit2-5.34.14.tar.gz
cd Minuit2-5.34.14
./configure \
    --disable-openmp \
    --prefix=/tmp/Minuit2-5.34.14
make all
make install
