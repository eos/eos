############
Installation
############

If you intend to use EOS without any modifications, then you should follow the steps outlined in `Installing with a Package Manager`_
or `Installing with the Python Package Index (PyPI)`_.
If you intend to modify EOS for your personal purpose, then you must follow the steps outlined in `Installing from Source`_.


*********************************
Installing with a Package Manager
*********************************


Ubuntu Linux
============

Pre-built binary packages of the most recent release of EOS for Ubuntu Linux LTS release 18.04, are available through the Packagecloud web service.
They can be installed through the following steps:

1. Create the file ``/etc/apt/sources.list.d/eos.list`` by running the following command in a shell

::

  DIST=$(lsb_release -sc)
  echo deb https://packagecloud.io/eos/eos/ubuntu/ $DIST main \
      | sudo tee /etc/apt/sources.list.d/eos.list

2. Add the official EOS repository's GPG key to your package manager's keychain by running:

::

  curl -L "https://packagecloud.io/eos/eos/gpgkey" 2> /dev/null | apt-key add - &>/dev/null

You might need to install the ``curl`` program first, if not already installed. This can be done by running the following command

::

  apt-get install curl

3. Install the pre-built binary packages and its dependencies by running the following commands

::

  apt-get update
  apt-get install eos

4. On older version of Ubuntu, you might need to update your installation of the Matplotlib Python package.
   We recommend to install the most up-to-date Python package with from the Python Package Index, using

::

  pip3 install --user h5py matplotlib numpy pypmc pyyaml scipy


MacOS X
=======

The most recent development version of EOS (and all of its dependencies) can be installed automatically using `Homebrew <https://brew.sh/>`_:

::

  brew tap eos/eos
  brew install --HEAD eos

We recommend to install the most up-to-date Python package with from the Python Package Index, using

::

  pip3 install --user h5py matplotlib numpy pypmc pyyaml scipy

***********************************************
Installing with the Python Package Index (PyPI)
***********************************************

Linux
=====

Pre-built binary packages of the most recent release of EOS are available through the `Python Package Index <https://pypi.org/>`_.
We are presently building these binary packages, also known as wheels, for the ``manylinux2014`` platform.
We assume that you have access to the ``pip3`` command.
If not, you need to install ``pip3`` through via your package manager, or ask your system administrator to install it.

EOS can be installed by running the following command:

::

  pip3 install --user eoshep


MacOS X
=======

We are presently not building Python wheels for the MacOS platform.

**********************
Installing from Source
**********************

The following instructions explain how to install EOS from source on a Linux-based operating system,
such as Debian or Ubuntu, or MacOS X.
Other flavors of Linux or UNIX-like operating systems will likely work as well.
However, note that we will exclusively use Debian/Ubuntu-specific package names when listing
the neccessary softwares that EOS depends upon.

Installing the dependencies on Linux
====================================

The dependencies can be roughly categorized as either system software, Python-related software, or scientific software.
EOS requires the following system software:

g++
  the GNU C++ compiler, in version 4.8.1 or higher;

autoconf
  the GNU tool for creating configure scripts, in version 2.69 or higher;

automake
  the GNU tool for creating makefiles, in version 1.14.1 or higher;

libtool
  the GNU tool for generic library support, in version 2.4.2 or higher;

pkg-config
  the freedesktop.org library helper, in version 0.28 or higher;

BOOST
  the BOOST C++ libraries (sub-libraries ``boost-filesystem`` and ``boost-system``);

yaml-cpp
  a C++ YAML parser and interpreter library.


The Python interface to EOS requires the additional software:

python3
  the Python interpreter in version 3.5 or higher, and required header files;

BOOST
  the BOOST C++ library ``boost-python`` for interfacing Python and C++;

h5py
  the Python interface to HDF5;

matplotlib
  the Python plotting library in version 2.0 or higher;

scipy
  the Python scientific library;

PyYAML
  the Python YAML parser and emitter library.

We recommend you install the above packages via your system's software management system.


EOS requires the following scientific software:

GSL
  the GNU Scientific Library \cite{GSL}, in version 1.16 or higher;
HDF5
  the \gls{HDF5} \cite{HDF5}, in version 1.8.11 or higher;
FFTW3
  the C subroutine library for computing the discrete Fourier transform;
Minuit2
  the physics analysis tool for function minimization, in version 5.28.00 or higher.


If you have administrator access to the computers on which you use EOS,
we recommend you install the above packages via your system's software management system.
The ``Minuit2`` is excempted from this recommendation, since its Debian/Ubuntu packages are affected by an unresolved bug.
Installing a prebuilt version of ``Minuit2`` is discussed in `Installing Minuit2 via APT`_.
Installing ``Minuit2`` from source is discussed in `Installing Minuit2 from Source`_.

On a Debian/Ubuntu based operating system you can install the prerequisite software with the ``APT`` package management system,
by running the following commands

::

  # for the 'System Software'
  sudo apt-get install g++ autoconf automake libtool pkg-config libboost-filesystem-dev libboost-system-dev libyaml-cpp-dev
  # for the 'Python Software'
  sudo apt-get install python3-dev libboost-python-dev python3-h5py python3-matplotlib python3-scipy python3-yaml
  # for the 'Scientific Software'
  sudo apt-get install libgsl0-dev libhdf5-serial-dev libfftw3-dev

Do not install the ``Minuit2`` software via ``APT``, since there is presently a bug in the Debian/Ubuntu packages,
which prevents EOS linking against the needed libraries.
We recommend that you upgrade ``matplotlib`` to the latest available version by running the following command:

::

  # for the 'pip3' command
  apt-get install python3-pip
  pip3 install matplotlib --user --upgrade


Installing Minuit2 via APT
--------------------------

There are pre-built binary package files for the ``Minuit2`` software available for the Ubuntu long-term-support releases 16.04 and 18.04 via the Packagecloud web service.
To use the EOS third-party repository, create a new file ``eos.list`` within the directory ``/etc/apt/sources.list.d`` with the following contents:

::

  deb https://packagecloud.io/eos/eos/ubuntu/ DIST main
  deb-src https://packagecloud.io/eos/eos/ubuntu/ DIST main

where you must replace the metavariable ``DIST`` with either ``xenial`` or ``bionic``, depending on your version of Ubuntu.
Add our repository's GPG key by making sure that the ``curl`` command line utility is installed, and

::

  curl -L "https://packagecloud.io/eos/eos/gpgkey" 2> /dev/null | sudo apt-key add -

You can then install the binary package through

::

  apt-get update
  apt-get install minuit2

You can then proceed with the EOS installation from source in `Installing EOS`_.


Installing Minuit2 from Source
------------------------------

We recommend to Minuit2 to be installed below ``/usr/local``, if you have administrator access to your computer.
In preparation for the installation, run the following command

::

  export PREFIX=/usr/local

Instead, if you do not have administrator access to your computer, we recommend to ``Minuit2`` to be installed below ``$HOME/.local``.
In preparation for the installation, run the following command

::

  export PREFIX=${HOME}/.local

When installing ``Minuit2`` from source, you need to disable the automatic support for parallelization with the OpenMP framework.
To this end, run the following commands

::

  mkdir /tmp/Minuit2
  pushd /tmp/Minuit2
  wget http://www.cern.ch/mathlibs/sw/5_28_00/Minuit2/Minuit2-5.28.00.tar.gz
  tar zxf Minuit2-5.28.00.tar.gz
  pushd Minuit2-5.28.00
  ./configure --prefix=$PREFIX --disable-openmp
  make all
  make install # Use 'sudo make install' if you install e.g. to 'PREFIX=/usr/local'
               # or a similarly privileged directory
  popd
  popd
  rm -R /tmp/Minuit2


Installing the Dependencies with Homebrew and PyPi
==================================================

You can install most of the prerequisite software via ``Homebrew``.
You will need to make ``Homebrew`` aware of the EOS third-party repository by running the following command in a shell

::

  brew tap eos/eos

To install the packages, run the following commands in a shell:

::

  # for the 'System Software'
  brew install autoconf automake libtool pkg-config boost yaml-cpp
  # for the 'Python Software'
  brew install python3 boost-python3
  # for the 'Scientific Software'
  brew install gsl hdf5 minuit2

You can now use the ``pip3`` command to install the remaining packages from the \package{PyPi} package index.

.. note::
    Due to problems with the Python 3 installation provided by Mac OS X, we strongly recommend to use instead the ``pip3`` programm
    provided by Homebrew, which should be available as ``/usr/local/bin/pip3``.

To install the remaining packages, run the following command in a shell

::

  pip3 install h5py matplotlib scipy PyYAML


Installing EOS
==============

You can obtain the EOS source code from the public Github repository.
To download it for the first time, clone the repository by running the following command:

::

  git clone -o eos -b master https://github.com/eos/eos.git

To install from the source code repository, you must first create all the necessary build scripts by running the following commands:

::

  cd eos
  ./autogen.bash

You must now decide where EOS will be installed.
To proceed we require you to set the environment variable ``PREFIX``.
We recommend to install to your home directory.
To do this, run the following command

::

  export PREFIX=${HOME}/.local

Next, you must configure the EOS build using the ``configure`` script.
The fo

1. To use the EOS Python interface you must pass ``--enable-python`` to the call ``configure``.
   The default is \cli{--disable-python}.

2. To use the ROOT software's internal copy of Minuit2 you must pass ``--with-minuit2=root`` to the call to ``configure``.

3. Otherwise, if you have installed ``Minuit2`` from source to a non-standard location you must specify its installation
   prefix by passing on ``--with-minuit2=MINUIT2-INSTALLATION-PREFIX``.
   If you followed the instructions in this manual, then the metavariable ``MINUIT2-INSTALLATION-PREFIX``
   corresponds to either ``/usr/local`` or ``$HOME/.local``.

The recommended configuration is achieved by running the following command

::

  ./configure \
      --prefix=$PREFIX \
      --enable-python \
      --enable-pmc

If the ``configure`` script finds any problems with your system, it will complain loudly.


After successful configuration, build EOS by running the following command

::

  make -j all

The ``-j`` option instructs the ``make`` programm to use all available processors to parallelize the build process.
#We strongly recommend testing the build by running the command

::

  make -j check VERBOSE=1

#within the build directory.
Please contact the authors if any test fails by opening an issue in the official `EOS Github repository <https://github.com/eos/eos>`_.
If all tests pass, install EOS by running the command

::

  make install # Use 'sudo make install' if you install e.g. to 'PREFIX=/usr/local'
               # or a similarly privileged directory

If you installd EOS to a non-standard location (i.e. not `/usr/local``),
to use it from the command line you must set up some environment variable.
For ``BASH``, which is the default Debian/Ubuntu shell, add the following lines to ``\$HOME/.bash_profile``:

::

  export PATH+=":$PREFIX/bin"
  export PYTHONPATH+=":$PREFIX/lib/python3.6/site-packages"

Note that in the above ``python3.6`` piece must be replaced by the appropriate Python version with which EOS was built.
You can determine the correct value by running the following command

::

  python3 -c "import sys; print('python{0}.{1}'.format(sys.version_info[0], sys.version_info[1]))"
