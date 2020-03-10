############
Installation
############

***********************************
Installation with a Package Manager
***********************************

Ubuntu Linux
============

Pre-built binary package of the most recent release of EOS for Ubuntu Linux LTS releases 16.04 and 18.04, as well as the current release 19.04, are available through the Packagecloud web service.
They can be installed through the following steps:

1. Create the file ``/etc/apt/sources.list.d/eos.list`` by running the following command in a shell

::

  DIST=$(lsb_release -sc)
  echo deb https://packagecloud.io/eos/eos/ubuntu/ $DIST main \
      | sudo tee /etc/apt/sources.list.d/eos.list

2. Add the official EOS repository's GPG key to your package manager's keychain by running:

::

  curl -L "https://packagecloud.io/eos/eos/gpgkey" 2> /dev/null | apt-key add - &>/dev/null

You might need to install the ``curl`` program first, if not already installed.

3. Install the pre-built binary packages and its dependencies by running:

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

************************
Installation from Source
************************

.. todo::

  Write section on installation from source
