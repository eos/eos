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
``options`` takes an ``eos.Options`` object and allows to specify *global* options (i.e. applied to all observables in the expression);
and ``expression`` is the expression string to be parsed.

We conclude with a concrete example

.. code:: ipython3

    eos.Observables().insert('B->Kll::R_K_example', R'(R_K)', eos.Unit.Unity(), eos.Options(),
                             '( <<B->Kll::BR;l=mu>>[q2_max=6, q2_min=>q2_mu_min] / <<B->Kll::BR;l=e>>[q2_max=6,q2_min=>q2_e_min] )')

    R_K = eos.Observable.make('B->Kll::R_K_example', eos.Parameters.Defaults(), eos.Kinematics(q2_e_min=1.1, q2_mu_min=1.1), eos.Options(**{'tag':'BFS2004'}))

    R_K.evaluate()   # should be ~1
