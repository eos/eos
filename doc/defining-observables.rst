********************
Defining Observables
********************

.. _defining_observables:

A user can define new observables for their analyses by combining existing observables, parameters, and kinematic variables within an arithmetic expression.
Observables defined in this way benefit from the same optimizations as all built-in observables, including multi-threaded evaluation and caching.
This section describes and exemplifies the syntax of the expression parser that interprets user defined python strings as new EOS observables.

Syntax
======

Arithmetic Expressions
~~~~~~~~~~~~~~~~~~~~~~

Basic rules are used to parse the input string

  * Whitespaces such as spaces and tabs do not have semantic meaning and are ignored when parsing the expression.
    Hence, they can be added arbitrarily for readability.

  * Expressions are built upon the arithmetic operations

     - addition ``+``,
     - subtraction ``-``,
     - multiplication ``*``,
     - division ``/``,
     - and exponentiation ``^``.

    Parentheses are recognized in expressions, and the usual precedence rules of arithmetics apply.

  * Expression typically reference one or more EOS objects (either observables or a parameters) by name, it must be encapsulated by two chevrons: ``<<...>>``.
    Names must be fulfill the requirements of :code:`eos.QualifiedName`. All features of a qualified name are supported,
    including specifying option.
    Names are assumed to point to an observable by default. If no observable of that name exists, a parameter of the same name is used.
    If neither an observable nor a parameter of that name exists, an exception is raised.
    Examples for a reference include ``<<mass::mu>>`` or ``<<B_u->lnu::BR>>``.

  * Expressions can reference a kinematic variable. Its name must be encapsulated by curly braces: ``{...}``.
    The kinematic variable can already be used by an observable within the expression. If no observable uses the kinematic variable,
    a new kinematic variable is automatically created. Examples for a kinematic variable include ``{q2}`` or ``{E_gamma}``.

The following strings are valid observable expressions

  * ``(<<mass::B_d>>^2 - 4 * <<mass::mu>>^2) ^ 0.5``

  * ``1.0 / <<B_u->lnu::BR;l=mu>>``

  * ``<<B->K::f_+(q2)>> * (1 - {q2} / <<mass::B_d^*>>^2)``

Dealing with Kinematic Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, a referenced observable's dependence on a kinematic variables is transferred to the newly defined observable.
For example, the observable defined by the expression ``1.0 / <<B->pilnu::dBR/dq2>>`` depends on the kinematic variable
``q2``, since the reference observable ``B->pilnu::dBR/dq2`` does.
When more than one observable appears in the expression, it is sometimes useful or even essential to rename kinematic variables,
i.e., create an alias for the name of the kinematic variable.
For example, when defining a lepton-flavour universality observable such as the ratio of branching fractions ``<<B->Dlnu::BR;l=tau>> / <<B->Dlnu::BR;l=mu>>``,
the kinematic variable ``q2`` should be renamed to ``q2_tau`` and ``q2_mu`` for the two observables, respectively, to avoid any ambiguity.
Moreover, phenomenological application sometimes involve the evaluation of an observable at a fixed value of a kinematic variable.
Both types of modification to the kinematic variables of a referenced observable are supported.
They are achieved by appending a kinematic specification to the reference of an EOS observable: ``<<...>>[...]``.

Two types of kinematics specification are supported

  * the ``=>`` operator renames the kinematic variable on its left-hand side to the name on its right-hand side.
    For example, the expression ``<<B_u->lnu::BR@l=mu>>[q2=>q2_mu] / <<B_u->lnu::BR@l=e>>[q2=>q2_e]`` expects
    two kinematic variable: ``q2_mu`` and ``q2_e``. The value of the former is only used in the evaluation of the numerator;
    the value of the latter is only used in the evaluation of the denominator.

  * the ``=`` operator fixes a kinematic variable to a given value.
    For example, the expression ``<<B->pilnu::BR>>[q2_min=0.1]`` expects only one kinematic variable ``q2_max``.
    The kinematic variable ``q2_min`` that the references observable ``B->pilnu::BR`` also depends on is
    fixed to the value ``0.1``.

Note that both types of specifications can be combined in comma-separated list in an arbitrary order, for example ``[q2_min=1.0,q2_max=>q2_max_mu]``.

Making a new Observable known
=============================

To make a new observable known to EOS, it must be added to the list of observables.
This can be done using the :method:`eos.Observables.insert` method.

.. code:: ipython3

    eos.Observables().insert(name, latex, unit, options, expression)

Here the arguments ``name``, ``latex`` and ``unit`` are the qualified name, the latex representation, and the unit of the new observable, respectively.
Known units include:


The argument ``options`` takes an :class:`eos.Options` object and allows to specify *global* options (i.e. applied to all observables in the expression),
for the newly defined observable.
The argument ``expression`` is Python string that is parsed as described above.

We conclude with a concrete example for the definition of the lepton-flavour universality observable :math:`R_K`:

.. code:: ipython3

    eos.Observables().insert('B->Kll::R_K_example', R'(R_K)', eos.Unit.Unity(), eos.Options(),
                             '( <<B->Kll::BR;l=mu>>[q2_max=6, q2_min=>q2_mu_min] / <<B->Kll::BR;l=e>>[q2_max=6,q2_min=>q2_e_min] )')

    R_K = eos.Observable.make('B->Kll::R_K_example', eos.Parameters.Defaults(), eos.Kinematics(q2_e_min=1.1, q2_mu_min=1.1), eos.Options(**{'tag':'BFS2004'}))

    R_K.evaluate()   # should yield ~1
