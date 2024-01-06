EOS: A software for Flavor Physics Phenomenology
================================================

EOS is a software framework for applications in high-energy physics; in particular flavor physics.
It is written in C++20, and provides both a C++20 and a Python3 API. The Python3 API is the recommended
interface.

EOS has been authored with three use cases in mind.

1. `Theory predictions and uncertainty estimates <use-cases.html#theory-predictions-and-their-uncertainties>`_
   of observables and further theoretical quantities in the field of flavor physics.
   EOS aspires to produce theory estimates and their inherent uncertainties of publication quality,
   and has produced such estimates in the past.

2. `Parameter Inference <use-cases.html#parameter-inference>`_ based on experimental measurements and/or theoretical constraints.
   For this use case, EOS defaults to the Bayesian framework of parameter inference.
   Moreover, EOS provides a large database of experimental measurements and theoretical constraints for immediate use.

3. `Production of pseudo events <use-cases.html#production-of-pseudo-events>`_ for a variety of flavor-physics-related processes using Monte Carlo methods.

Installation instructions can be found `here <installation.html>`_.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   user-guide/index
   faq
   reference/index
