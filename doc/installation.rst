############
Installation
############

We distribute EOS as packaged software for end-users and also as source code for developers.


*****
Linux
*****

Installation with `pip3`
========================

This is the recommended installation method using the 'package installer for Python' `pip3`:

::

  pip3 install --user eoshep

This way, you will install the latest regular release of EOS.
You can now use the EOS Python module via ``import eos``. In addition, some command line scripts have become available.
To update EOS, run ``pip3 install --user --upgrade eoshep``.

This installs the most recent development version of EOS:

::

  pip3 install --user --pre eoshep

The development version corresponds to the most recent commit in the Github EOS repository.
To update, run ``pip3 install --user --pre --upgrade eoshep``.


Installation with `apt` (Ubuntu, Debian and Siblings)
=====================================================

EOS is available through the Packagecloud web service.
The installation requires several steps:


1. Create the file ``/etc/apt/sources.list.d/eos.list`` by running the following command in a shell

::

  DIST=$(lsb_release -sc)
  echo deb https://packagecloud.io/eos/eos/ubuntu/ $DIST main \
      | sudo tee /etc/apt/sources.list.d/eos.list

2. Add the official EOS repository's GPG key to your package manager's keychain by running:

::

  curl -L "https://packagecloud.io/eos/eos/gpgkey" 2> /dev/null | apt-key add - &>/dev/null

You might need to install the ``curl`` program first, if not already installed. This can be done by running the following command:

::

  apt-get install curl

3. Install the pre-built binary package and its dependencies by running the following commands:

::

  apt-get update
  apt-get install eos

4. On older versions of Ubuntu, you might need to update your installation of the Matplotlib Python package.
   We recommend installing the most up-to-date Python package from the Python Package Index, using

::

  pip3 install --user --upgrade h5py matplotlib numpy pypmc pyyaml scipy



*****
macOS
*****

Installation with `brew`
========================

EOS (and all of its dependencies) can be installed automatically using `Homebrew <https://brew.sh/>`_:

::

  brew tap eos/eos
  brew install --HEAD eos

We recommend installing the most up-to-date Python package from the Python Package Index, using

::

  pip3 install --user --upgrade h5py matplotlib numpy pypmc pyyaml scipy



*****************************
Build and Install from Source
*****************************

The following instructions explain how to install EOS from source on a Linux-based operating system,
such as Debian or Ubuntu, and macOS.
Other flavours of Linux or UNIX-like operating systems will likely work as well.
However, note that we will exclusively use Debian/Ubuntu-specific package names when listing
software dependencies.

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


If you have administrator access to the computers on which you use EOS,
we recommend you install the above packages via your system's software management system.

On a Debian/Ubuntu-based operating system you can install the prerequisite software with the ``apt`` package management system,
by running the following commands:

::

  # for the 'System Software'
  sudo apt-get install g++ autoconf automake libtool pkg-config libboost-filesystem-dev libboost-system-dev libyaml-cpp-dev
  # for the 'Python Software'
  sudo apt-get install python3-dev libboost-python-dev python3-h5py python3-matplotlib python3-scipy python3-yaml
  # for the 'Scientific Software'
  sudo apt-get install libgsl0-dev

We recommend that you upgrade ``matplotlib`` to the latest available version by running the following command:

::

  # for the 'pip3' command
  apt-get install python3-pip
  pip3 install matplotlib --user --upgrade


Installing the dependencies on macOS with Homebrew and PyPi
===========================================================

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
  brew install gsl

You can now use the ``pip3`` command to install the remaining packages from the \package{PyPi} package index.

.. note::
    Due to problems with the Python 3 installation provided by macOS, we strongly recommend using instead the ``pip3`` program
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
We recommend installing to your home directory.
To do this, run the following command:

::

  export PREFIX=${HOME}/.local

Next, you must configure the EOS build using the ``configure`` script.
To use the EOS Python interface you must pass ``--enable-python`` to the call ``configure``.
The default is ``--disable-python``.

The recommended configuration is achieved by running the following command:

::

  ./configure \
      --prefix=$PREFIX \
      --enable-python

If the ``configure`` script finds any problems with your system, it will complain loudly.

The flag ``with-boost-python-suffix`` might be necessary, depending on the installation of Python and BOOST.
For example, when `boost-python3` is installed on macOS via `brew`, you can find the suffix by inspecting the installed libraries:
``ls /usr/local/lib/libboost_python*``
might yield
``/usr/local/lib/libboost_python39.a``.
Here, the flag ``--with-boost-python-suffix=39`` is required.

After successful configuration, build EOS by running the following command:

::

  make -j all

The ``-j`` option instructs the ``make`` program to use all available processors to parallelize the build process.
#We strongly recommend testing the build by running the command

::

  make -j check VERBOSE=1

#within the build directory.
Please contact the authors if any test fails by opening an issue in the official `EOS Github repository <https://github.com/eos/eos>`_.
If all tests pass, install EOS by running the command

::

  make install # Use 'sudo make install' if you install e.g. to 'PREFIX=/usr/local'
               # or a similarly privileged directory

If you installed EOS to a non-standard location (i.e. not ``/usr/local``),
to use it from the command line you must set up some environment variable.
For ``BASH``, which is the default Debian/Ubuntu shell, add the following lines to ``\$HOME/.bash_profile``:

::

  export PATH+=":$PREFIX/bin"
  export PYTHONPATH+=":$PREFIX/lib/python3.6/site-packages"

Note that in the above the ``python3.6`` piece must be replaced by the appropriate Python version with which EOS was built.
You can determine the correct value by running the following command:

::

  python3 -c "import sys; print('python{0}.{1}'.format(sys.version_info[0], sys.version_info[1]))"
