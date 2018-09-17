EOS - A HEP Program for Flavor Observables
==========================================

EOS is a software package that addresses several use cases in the field of
high-energy flavor physics (HEP):

 1. calculation and uncertainty estimation of flavor observables within
   various models,
 2. Bayesian inference of parameters from experimental and/or theoretical
   constraints, and
 3. sampling process-specific probability density functions.

An up-to-date list of EOS related publications can be found [here](https://eos.github.io/publications/).

EOS is written in C++14, with an optional interface to Python, and depends on
as a small set of external software libraries:

 - the GNU Scientific Library (libgsl),
 - a subset of the BOOST C++ libraries,
 - the Hierarchical Data Format v5 library (libdf5),
 - the minimizer Minuit2 (as of ROOT version 5.14.00 or later),
 - the Population Monte Carlo (PMC) library pmclib (optional),
 - the Python interpreter (optional).

For details on these dependencies we refer to the [user manual](https://eos.github.io/manual/manual.pdf).

Installation
------------

Presently EOS supports installation from source only. For Ubuntu users, two of the external software
dependencies that are not available from the main repositories are provided in the official
[EOS repository](https://packagecloud.io/eos/eos).

For instructions on how to build and install EOS on your computer please have a
look at the [user manual](https://eos.github.io/manual/manual.pdf).

Authors and Contributors
------------------------

The main authors are:

 * Danny van Dyk <danny.van.dyk@gmail.com>,
 * Frederik Beaujean <frederik.beaujean@lmu.de>,
 * Christoph Bobeth <christoph.bobeth@gmail.com>,

with code contributions by:

 * Marzia Bordone,
 * Thomas Blake,
 * Elena Graverini,
 * Nico Gubernari,
 * Stephan Jahn,
 * Ahmet Kokulu,
 * Bastian Müller,
 * Stefanie Reichert,
 * Eduardo Romero,
 * Rafael Silva Coutinho,
 * Ismo Tojiala,
 * Keri Vos,
 * Christian Wacker.

We would like to extend our thanks to the following people whose input and
support were most helpful in either the development or the maintenance of EOS:

 * Gudrun Hiller
 * David Leverton
 * Ciaran McCreesh
 * Hideki Miyake
 * Konstantinos Petridis
 * Alexander Shires

Contact
-------

For additional information, please contact any of the main authors. If you want to report an
error or file a request, please file an issue [here](https://github.com/eos/eos/issues).
