#!/bin/bash
# Install prerequisites for TravisCI builds

set -ve

if [[ "$TRAVIS_OS_NAME" != "osx" ]]; then
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
    sudo apt-get update
    sudo apt-get install -y --allow-unauthenticated g++-6 clang-3.9
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 60 --slave /usr/bin/g++ g++ /usr/bin/g++-6
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 40 --slave /usr/bin/g++ g++ /usr/bin/g++-4.8
    sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.9 60 --slave /usr/bin/clang++ clang++ /usr/bin/clang++-3.9
    sudo apt-get install -y cmake git libgfortran3 gfortran openmpi-bin openmpi-common openmpi-doc libopenmpi-dev libblas-dev liblapack-dev libfftw3-dev libgmp-dev hdf5-tools libhdf5-serial-dev python-h5py python-dev python-numpy python-scipy python-jinja2 python-virtualenv python-matplotlib python-tornado python-zmq python-mpi4py python-mako clang-format-3.9 libclang-3.9-dev python-clang-3.9 python-sphinx libjs-mathjax valgrind
  - export OMPI_CC=${CC}
  - export OMPI_CXX=${CXX}
else
    brew cask uninstall oclint # TravisCI bug, see https://github.com/travis-ci/travis-ci/issues/8826#issuecomment-350103392
    brew install hdf5
    brew list
fi
