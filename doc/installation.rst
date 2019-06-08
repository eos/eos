############
Installation
############

***********************************
Installation with a Package Manager
***********************************

Ubuntu Linux
============

Pre-built binary package of the most recent release of EOS for Ubuntu Linux LTS 16.04 and 18.04 are available through the Packagecloud web service. They can be installed through the following steps:

1. Create the file ``/etc/apt/sources.list.d/eos.list`` with the following contents

::

  deb https://packagecloud.io/eos/eos/ubuntu/ DIST main
  deb-src https://packagecloud.io/eos/eos/ubuntu/ DIST main

where you must replace ``DIST`` with either ``xenial`` (Ubuntu 16.04) or ``bionic`` (Ubuntu 18.04), depending on your Ubuntu installation.

2. Add the official EOS repository's GPG key by running:

::

  curl -L "https://packagecloud.io/eos/eos/gpgkey" 2> /dev/null | apt-key add - &>/dev/null

3. Install the pre-built binary packages and its dependencies by running:

::

  apt-get update
  apt-get install eos


MacOS X
=======

The most recent release of EOS (and all of its dependencies) can be installed automatically using `Homebrew <https://brew.sh/>`_:

::

  brew tap eos/eos
  brew install eos


Similarly, the most recent development version of EOS (and all of its dependencies) can be installed automatically using:

::

  brew tap eos/eos
  brew install --HEAD eos


***********************
Installation from Souce
***********************

.. todo::

  Write section on installation from source
