EOS: A software for Flavor Physics Phenomenology
================================================

EOS is a software framework for applications in high-energy physics; in particular flavor physics.
It is written in C++20, and provides both a C++20 and a Python3 API. The Python3 API is the recommended
interface. Binary packages for ``x86_64`` and ``ARMv8`` instruction sets are available and
installation instructions can be found `here <installation.html>`_.

EOS has been authored with three use cases in mind. The following three documents provide a low-level
introduction to the respective use case, and how Python code can use EOS to address it.

1. `Theory predictions and uncertainty estimates <user-guide/predictions.html>`_
   of observables and further theoretical quantities in the field of flavor physics.
   EOS aspires to produce theory estimates and their inherent uncertainties of publication quality,
   and has produced such estimates in the past.

2. `Parameter Inference <user-guide/inference.html>`_ based on experimental measurements and/or theoretical constraints.
   For this use case, EOS defaults to the Bayesian framework of parameter inference.
   Moreover, EOS provides a large database of experimental measurements and theoretical constraints for immediate use.

3. `Production of pseudo events <user-guide/simulation.html>`_ for a variety of flavor-physics-related processes using Monte Carlo methods.

All three use cases can be addressed with EOS' high-level interface, which utilizes analysis files written in YAML;
their power to organise one or more analyses is demonstrated in the `analysis files <user-guide/analysis-organisation.html>`_ section.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   user-guide/index
   faq
   reference/index
