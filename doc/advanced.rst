##############
Advanced Usage
##############

*********************
The Expression Parser
*********************

A user can define new observables for their analyses.
These new observables must be arithmetic expressions of already defined observables and/or parameters.
Observables defined in this way benefit from the same optimizations as all built-in observables, including multi-threaded evaluation and caching.
This paragraph describes and exemplifies the syntax of the expression parser that interprets user defined python strings as new EOS observables.

The Construction of Expressions
===============================

The Rules
~~~~~~~~~

Basic rules are used to parse the input string

  * Spaces are ignored and can be added arbitrarily for readability.
  * The parser supports usual arithmetic operations ``+``, ``-``, ``*``, ``/`` and ``^`` and parenthesized expressions; usual precedence rules of arithmetics apply.
  * ``<<...>>`` encapsulates the name of an EOS object, which must either be the name of a parameter or of an observable.
    Any such must therefore adhere to the restrictions of ``eos.QualifiedName``, for example ``<<mass::mu>>`` or ``<<B_u->lnu::BR>>``.
  * ``{...}`` enscapsulates the name of a kinematic variable, e.g. ``{q2}``.

The following strings are valid observable expressions

  * ``"(<<mass::B_d>>^2 - 4 * <<mass::mu>>^2) ^ 0.5"``
  * ``"1.0 / <<B_u->lnu::BR@l=mu>>"``
  * ``"<<B->K::f_+(q2)>> * (1 - {q2} / <<mass::B_d>>^2)"``

Aliasing
~~~~~~~~

By default, the kinematic arguments of the observables are transferred to the expression.
For example, the expression ``1.0 / <<B->pilnu::dBR/dq2>>`` will expect a kinematic specification for ``q2``
(either through an ``eos.Kinematic`` object or indirectly in a plotting routine).
When more than one observable appears in the expression, it is useful to rename the kinematic variables.
This can be done via an alias specification ``<<...>>[...]``.
Two types of specification are supported

  * the ``=`` operator fixes a kinematic variable to a given value. E.g. ``<<B->pilnu::BR>>[q2_min=0.1]`` only expects a specification for ``q2_max``.
  * the ``=>`` operator renames the kinematic variable on its left-hand side to the name on the right-hand side.
    E.g. ``<<B_u->lnu::BR@l=mu>>[q2=>q2_mu] / <<B_u->lnu::BR@l=e>>[q2=>q2_e]`` requires two kinematic specifications, ``q2_mu`` and ``q2_e``.

Note that these specifications can be combined in comma-separated list, for example ``[q2_min=1.0, q2_max=>q2_mu]``.

The Insert Method
=================

Once a new observable is defined via its expression string, it can be added to the list of observables via the ``insert`` method

.. code:: ipython3

    eos.Observables().insert(name, latex, unit, options, expression)

where ``name``, ``latex`` and ``unit`` are the (Qualified)name, the latex representation and the unit of the new observable;
``options`` takes an ``eos.Options`` object and allows to specify `global` options (i.e. applied to all observables in the expression);
and ``expression`` is the expression string to be parsed.

We conclude with a concrete example

.. code:: ipython3

    eos.Observables().insert('B->Kll::R_K_example', R'(R_K)', eos.Unit.Unity(), eos.Options(),
                             '( <<B->Kll::BR;l=mu>>[q2_max=6, q2_min=>q2_mu_min] / <<B->Kll::BR;l=e>>[q2_max=6,q2_min=>q2_e_min] )')

    R_K = eos.Observable.make('B->Kll::R_K_example', eos.Parameters.Defaults(), eos.Kinematics(q2_e_min=1.1, q2_mu_min=1.1), eos.Options(**{'tag':'BFS2004'}))

    R_K.evaluate()   # should be ~1

******************************
The EOS Command-Line Interface
******************************

Although using EOS within an interactive Jupyter notebook on your personal computer or laptop
is useful to prototype an analysis, this approach sometimes suffers from limited computing power.
To circumvent this problem, you can alternatively

  * use EOS in Jupyter interactively on a remote workstation computer via an SSH tunnel (see the `FAQ <faq>`_);
  * use EOS on remote workstations or compute clusters via the command-line interface.

In the following we document the command-line interface and the file format used in conjunction with it.

.. note::

   The EOS command-line interface is completely optional and does not provide any means beyond the
   interactive Python interface.

The Analysis Description Format
===============================

EOS uses a YAML file to describe the individual steps of one or more statistical analyses.
At the top level, the format includes the following YAML keys:

 - ``priors`` (**mandatory**) --- The list of priors within the analysis.
 - ``likelihoods`` (**mandatory**) --- The list of likelihoods within the analysis..
 - ``posteriors`` (**mandatory**) --- The list of posteriors within the analysis.
 - ``predictions`` (**optional**) --- The list of theory predictions within the analysis.

Describing Priors
~~~~~~~~~~~~~~~~~

The ``priors`` key contains a list of *named* priors. Each prior has two mandatory keys:

  - ``name`` (**mandatory**) --- The unique name of this prior.
  - ``parameters`` (**mandatory**) --- The ordered list of parameters described by this prior.

The description of each individual parameter follows the prior description used in the
:class:`Analysis <eos.Analysis>` constructor.


Describing Likelihoods
~~~~~~~~~~~~~~~~~~~~~~

The ``likelihoods`` key contains a list of *named* likelihoods. Each likelihood has two mandatory keys:

  - ``name`` (**mandatory**) --- The unique name of this likelihood.
  - ``constraints`` (**mandatory**) --- The ordered list of EOS constraint names that comprise this likelihood.

Describing Posteriors
~~~~~~~~~~~~~~~~~~~~~

The ``posteriors`` key contains a list of *named* posteriors. Each posterior can contain several keys:

  - ``name`` (**mandatory**) --- The unique name of this posterior.
  - ``global_options`` (**optional**) --- A key/value map providing global options, i.e., options that apply to all observables used by this posterior.
  - ``prior`` (**mandatory**) --- The ordered list of named priors that are used as part of this posterior.
  - ``likelihood`` (**optional**) --- The ordered list of named likelihoods that are used as part of this posterior.
  - ``fixed_parameter`` (**optional**) --- A key/value map providing values for parameters that deviate from the default values.

Example
~~~~~~~

.. toggle-header::
   :header: Example `examples/cli/btopilnu.analysis <https://github.com/eos/eos/tree/master/examples/cli/btopilnu.analysis>`_

   .. literalinclude:: ../examples/cli/btopilnu.analysis
      :language: YAML

|

The Command-Line Interface
==========================

.. argparse::
   :filename: ../src/scripts/eos-analysis
   :func: _parser
   :prog: eos-analysis
   :nodescription:
   :nodefault:

   The ``eos-analysis`` script provides several subcommands that

    * inspect the analysis file;
    * sample from a posterior density with Monte Carlo methods;
    * perform auxiliary tasks on intermediate results.

   The output of these commands are stored on disk as directories filled with YAML files
   (for descriptions and small numerical datasets) and Numpy datafiles (for samples).
   The datafiles can be access with the classes documented as part of the :obj:`eos.data` module.
