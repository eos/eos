[![PyPi version](https://img.shields.io/pypi/v/eoshep)](https://img.shields.io/pypi/v/eoshep)
[![Build Status](https://github.com/eos/eos/actions/workflows/pypi-build+check+deploy.yaml/badge.svg)](https://github.com/eos/eos/actions/workflows/pypi-build+check+deploy.yaml)
[![Build Status](https://github.com/eos/eos/actions/workflows/ubuntu-build+check+deploy.yaml/badge.svg)](https://github.com/eos/eos/actions/workflows/ubuntu-build+check+deploy.yaml)
[![Discord](https://img.shields.io/discord/808999754989961236.svg?label=&logo=discord&logoColor=ffffff&color=7389D8&labelColor=6A7EC2)](https://discord.gg/hyPu7f7K6W)


![EOS logo](https://eos.github.io/images/github-eos-logo.png)

EOS - A software for Flavor Physics Phenomenology
=================================================

EOS is a software package that addresses several use cases in the field of
high-energy flavor physics:

 1. [theory predictions of and uncertainty estimation for flavor observables](https://eos.github.io/doc/use-cases.html#theory-predictions-and-their-uncertainties)
   within the Standard Model or within the Weak Effective Theory;
 2. [Bayesian parameter inference](https://eos.github.io/doc/use-cases.html#parameter-inference)
    from both experimental and theoretical constraints; and
 3. [Monte Carlo simulation of pseudo events](https://eos.github.io/doc/use-cases.html#pseudo-event-simulation) for flavor processes.

An up-to-date list of publications that use EOS can be found [here](https://eos.github.io/publications/).

EOS is written in C++20 and designed to be used through its Python 3 interface,
ideally within a Jupyter notebook environment.
It depends on as a small set of external software:

 - the GNU Scientific Library (libgsl),
 - a subset of the BOOST C++ libraries,
 - the Python 3 interpreter.

For details on these dependencies we refer to the [online documentation](https://eos.github.io/doc/installation.html#installing-the-dependencies-on-linux).

Installation
------------

EOS supports several methods of installation. For Linux users, the recommended method
is installation via PyPI:
```
pip3 install eoshep
```
Development versions tracking the master branch are also available via PyPi:
```
pip3 install --pre eoshep
```

For instructions on how to build and install EOS on your computer please have a
look at the [online documentation](https://eos.github.io/doc/installation.html).

Contact
-------

If you want to report an error or file a request, please file an issue [here](https://github.com/eos/eos/issues).
For additional information, please contact any of the main authors, e.g. via our [Discord server](https://discord.gg/hyPu7f7K6W).

Authors and Contributors
------------------------

The main authors are:

 * Danny van Dyk <danny.van.dyk@gmail.com>,
 * Frederik Beaujean,
 * Christoph Bobeth,
 * Carolina Bolognani <carolinabolognani@gmail.com>,
 * Nico Gubernari <nicogubernari@gmail.com>,
 * Meril Reboud <merilreboud@gmail.com>,

with further code contributions by:

 * Marzia Bordone,
 * Thomas Blake,
 * Lorenz Gaertner,
 * Elena Graverini,
 * Stephan Jahn,
 * Ahmet Kokulu,
 * Viktor Kuschke,
 * Stephan Kürten,
 * Philip Lüghausen,
 * Bastian Müller,
 * Filip Novak,
 * Stefanie Reichert,
 * Eduardo Romero,
 * Rafael Silva Coutinho,
 * Ismo Tojiala,
 * K. Keri Vos,
 * Christian Wacker.

We would like to extend our thanks to the following people whose input and
support were most helpful in either the development or the maintenance of EOS:

 * Gudrun Hiller
 * Gino Isidori
 * David Leverton
 * Thomas Mannel
 * Ciaran McCreesh
 * Hideki Miyake
 * Konstantinos Petridis
 * Nicola Serra
 * Alexander Shires
