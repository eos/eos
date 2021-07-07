#########
Basics
#########


.. note::

   The following examples are also available as an interactive `Jupyter <https://jupyter.org/>`_ notebook `from here <https://github.com/eos/eos/blob/master/examples/basics.ipynb>`_.

********************
Importing the module
********************

Before you can begin using EOS, you need to load it into your Python
interpreter. You achieve this with the ``import`` keyword.

.. code:: ipython3

    import eos

The import statement should go through without any error messages. If
there are errors, they might be related to the following issues: 

- missing Python packages (try to run ``!pip3 install -U matplotlib numpy scipy``);
- you might have forgotten to tell the Python interpreter where you installed EOS to ( this should only happen if you have **not** installed EOS using ``pip``)

**************************************
(Almost) Everything can be a Parameter
**************************************

In EOS we make a lot of use of the ``eos.Parameters`` class. It provides
access to one real-valued and scalar parameter, among a diverse and
large set of parameters. You cannot create an object of type
``eos.Parameters`` yourself. Instead, you can create a new set of
parameters, from which your favourite parameter can be extracted and
looked at.

.. code:: ipython3

    parameters = eos.Parameters()

The variable ``parameters`` now contains all of your parameters. You can
use ``display(parameters)`` to see the full list, which is rather
lengthy. Let’s have a look at something more manageable: the parameter
for the muon mass:

.. code:: ipython3

    parameters['mass::mu']

.. list-table::
   :widths: 25, 25

   * - math:m_{\mu}
     - eos.Parameter
   * - current value
     - 0.105658
   * - default value
     - 0.1050658



Let’s highlight some aspects of this output:

- Parameters have a rigid naming scheme, and must contain at least two parts: the *prefix* part, which is separated by ``::`` from the *name part*.
- The masses of all fields, elementary or composite, are collected in the ``mass`` namespace.
- EOS uses Giga electron Volt (GeV) for all masses, momenta and energies.
- A parameter has a current value, and a default value. Both are shown, which makes it easy to see if a parameter has been modified.

Let’s do just that: we can modify the value of the muon mass. This is
unlikely to be necessary, since it has been measured to such a high
precision, but useful to illustrate how to work with parameters.

.. code:: ipython3

    m_mu = parameters['mass::mu']
    m_mu.set(1.779)  # we just made the muon as heavy as the tauon!
    display(m_mu)



.. list-table::
   :widths: 25, 25

   * - math:m_{\mu}
     - eos.Parameter
   * - current value
     - 1.779
   * - default value
     - 0.1050658


It’s important to understand that the object is not living in
“isolation”, it is still part of the same parameter set
(``parameters``). Accessing the muon mass through parameters again will
therefore show the same information as ``m_mu``, but this time the
modified mass value will be included as “current value” .

.. code:: ipython3

    parameters['mass::mu']



.. list-table::
   :widths: 25, 25

   * - math:m_{\mu}
     - eos.Parameter
   * - current value
     - 1.779
   * - default value
     - 0.1050658



Brilliant! But what about if we want or need to have two independent
sets of parameters? This could be needed, for example, to compare theory
predictions for two different choices of parameters.

.. code:: ipython3

    more_parameters = eos.Parameters()
    more_parameters['mass::mu']




.. list-table::
   :widths: 25, 25

   * - math:m_{\mu}
     - eos.Parameter
   * - current value
     - 0.1050658
   * - default value
     - 0.1050658



Very good! We see that each call to ``eos.Parameters`` creates an
independent set of parameters, which start out with the defaults. You
can browse the full list of parameters known to EOS by running
``display(eos.Parameters())`` or by going online to the `EOS
documentation <https://eos.github.io/doc/parameters>`__.

You can access a parameter’s name and LaTeX representation using the
``name`` or ``latex`` method seen below. The value is obtained using the
``evaluate`` method.

.. code:: ipython3

    display(m_mu.name())
    display(m_mu.latex())
    display(m_mu.evaluate())



.. parsed-literal::

    'mass::mu'



.. parsed-literal::

    ':math:m_{\mu}'



.. parsed-literal::

    1.779


You can also treat a Parameter object just like any other Python object.
It can be part of a ``list``, a ``dict``, or a ``tuple``:

.. code:: ipython3

    lepton_masses = [parameters['mass::' + l] for l in ['e', 'mu', 'tau']]
    [display(p) for p in lepton_masses]
    translation = { p.name(): p.latex() for p in lepton_masses}
    display(translation)



.. list-table::
   :widths: 25, 25

   * - math:m_e
     - eos.Parameter
   * - current value
     - 0.000510999
   * - default value
     - 0.000510999



.. list-table::
   :widths: 25, 25

   * - math:m_{\mu}
     - eos.Parameter
   * - current value
     - 1.779
   * - default value
     - 0.105658



.. list-table::
   :widths: 25, 25

   * - math:m_{\tau}
     - eos.Parameter
   * - current value
     - 1.77682
   * - default value
     - 1.77682



.. parsed-literal::

    {'mass::e': '$m_e$', 'mass::mu': '$m_\\mu$', 'mass::tau': '$m_\\tau$'}


These properties make it possible to bind a function to an arbitrary
number of parameters, let the function evaluate these parameters in a
computationally efficient way, and let the user change these parameters
at a whim. Parameters are meant to be **shared**, i.e., a single set of
parameters is meant to be used by arbitrary number of functions. The
sharing of parameters is cause for the versatility in the EOS use cases.

**********************************
Kinematics and kinematic variables
**********************************

EOS makes also plentiful use of the ``eos.Kinematics`` class. Objects of
this class are used to store a set of real-valued and scalar kinematic
variables by name. Contrary to the parameters, there are no default
variables or values. This would not make much sense: kinematic variables
pertain to a single process. Consequently, naming of kinematic variables
does not require a distinction by any sort of ``prefix``. We will see
what kinematic variables a process requires in a later section. You can
create an empty set of kinematic variables as follows:

.. code:: ipython3

    kinematics = eos.Kinematics()

We can populate this object with a choice of kinematic variables: 

- We generally use ``q2`` or ``p2`` to denote the square of four momentum ``q`` or ``p``. - We generally use ``E_pi`` or ``E_gamma`` for the energy of a final state in the rest frame of the respective initial state.
- We usually parametrize helicity angles via their cosines, e.g., ``cos(theta_l)`` or ``cos(theta_pi)``. As for parameters, we use powers of GeV as units.

.. code:: ipython3

    k1 = kinematics.declare('q2',             1.0)   # 1 GeV^2
    k2 = kinematics.declare('E_pi',           0.139) # 139 MeV, a pion at rest!
    k3 = kinematics.declare('cos(theta_pi)', -1.0)   # negative values are OK!

.. code:: ipython3

    display(kinematics)

.. list-table::
   :widths: 25, 25

   * - q2
     - 1.0
   * - E_pi
     - 0.139
   * - cos(theta_pi)
     - -1.0



There are two ways to directly create the ``eos.Kinematics`` object
populated with the variables you need. The first way works if all your
kinematic variables names are also valid Python identifiers. (In the
above example, ``cos(...)`` is **not** a valid identifier.)

.. code:: ipython3

    kinematics = eos.Kinematics(q2=1.0, E_pi=0.139)
    display(kinematics)



.. list-table::
   :widths: 25, 25

   * - q2
     - 1.0
   * - E_pi
     - 0.139



The second way works also for variable names that are not valid
identifiers. It uses Python *keyword arguments*:

.. code:: ipython3

    kinematics = eos.Kinematics(**{
        'q2': 1.0, 'E_pi': 0.139, 'cos(theta_l)': -1.0
    })
    display(kinematics)



.. list-table::
   :widths: 25, 25

   * - q2
     - 1.0
   * - E_pi
     - 0.139
   * - cos(theta_l)
     - -1.0



We can extract a kinematic variable from the set using the ``[...]``
operator:

.. code:: ipython3

    k1 = kinematics['q2']
    display(k1)


.. list-table::
   :widths: 25, 25

   * - q2
     - (eos.KinematicVariable)
   * - current value
     - 1.0



We can also modify the set by setting the value of an individual
kinematic variable:

.. code:: ipython3

    k1.set(16.0)
    display(kinematics)



.. list-table::
   :widths: 25, 25

   * - q2
     - 16.0
   * - E_pi
     - 0.139
   * - cos(theta_l)
     - -1.0



Kinematic variables and their naming *usually* pertain to only a single
observable. When creating observables, we therefore *usually* create an
independet set of kinematic variables per observable. Nevertheless, it
is possible to create observables that have a common set of kinematic
variables. This makes it possible to investigate correlations among
observables that share a kinematic variable (e.g., LFU ratios such as
:math:`R_K` as a functions of the lower dilepton momentum cut-off).

****************************
Options and what they impact
****************************

EOS allows us to modify the behaviour of processes through objects of
the ``eos.Option`` class. In many cases, the processes have a default
set of options, e.g., the process’ default choice of hadronic form
factors, its default choice of BSM model, and so on. In some cases, it
does not make sense to have a default choice, e.g., when evaluating a
single hadronic form factor. In such cases, you will see an error
expressed through a Python exception if the mandatory option is not
specified by you!

Contrary to parameters and kinematic variables, EOS does not permit to
change a process’ options after the creation. To make this abundantly
clear: if you change an ``eos.Option`` object after it has been used
(in, e.g. an observable), none of your modifications to the options will
have any effect on the user (again, e.g. an observable)

You can create a new and empty set of options as follows:

.. code:: ipython3

    options = eos.Options()

We can now populate this object with individual options. Options are
pairs of strings. Within each pair, we refer to the first element as the
**key** and to the second element as the **value**. Typical keys
include:
- ``model``, to select a BSM model;
- ``form-factors``, to select the parametrization of the hadronic form factors in semileptonic decays;
- ``l``, to select a lepton flavour;
- ``q``, to select a quark flavour (typically for a spectator quark).

Option values are specific to both the process and the option key: -
``model`` can typically take values such as ``SM`` (for the Standard
Model), ``CKMScan`` (to parametrize each CKM matrix element and fit for
absolute value or complex phase), and ``WilsonScan`` (to parametrize the
Wilson coefficients of the Weak Effective Theory); - ``form-factors``
can typically take values that identify a single paper (e.g. ``BSZ2015``
for a parametrization used in Bharucha, Straub, Zwicky 2015); - ``l`` can
typically take values ``e``, ``mu`` and ``tau``; - ``q`` can typically
take values ``u``, ``d``, ``s``, and ``c``.

Option keys are specific to a process. Presently, there is no way to list
all option keys that a process understands, or to see their possible
values. We are working on that!

Adding options to the set can be achieved by the following:

.. code:: ipython3

    options.set('model', 'CKMScan')
    options.set('form-factors', 'BSZ2015')
    options.set('l', 'mu')                 # Since we are all so "cautiously excited"!
    options.set('q', 's')

.. code:: ipython3

    options



.. list-table::
   :widths: 25, 25

   * - form-factor
     - BSZ2015
   * - l
     - mu
   * - model
     - CKMScan
   * - q
     - s



There are two ways to directly create an ``eos.Options`` object
populated with the options you need. The first way works if all your
option keys are also valid Python identifiers. (In the above example
``form-factors`` is not a valid identifier.)

.. code:: ipython3

    options = eos.Options(l='mu', q='s', model='CKMScan')
    display(options)



.. list-table::
   :widths: 25, 25

   * - l
     - mu
   * - model
     - CKMScan
   * - q
     - s



The second way works also for option keys that are not valid
identifiers. It uses Python “keyword arguments”:

.. code:: ipython3

    options = eos.Options(**{
        'form-factors': 'BSZ2015',
        'model': 'WilsonScan',
        'l': 'tau',
        'q': 's'
    })
    display(options)



.. list-table::
   :widths: 25, 25

   * - form-factors
     - BSZ2015
   * - l
     - mu
   * - model
     - WilsonScan
   * - q
     - s



*******************************************************
Applying what we learned: creating our first Observable
*******************************************************

EOS uses the class ``eos.Observable`` to provide theory predictions for
a variety of processes and their associated observables. A list of all
observables known to EOS is available online in the documentation, or
via:

.. code:: ipython3

     eos.Observables()

To create an observable, we require: 

- its name, 
- a set of parameters that it will be bound to, 
- a set of kinematic variables that it will be bound to, 
- a set of options.

.. code:: ipython3

    observable1 = eos.Observable.make('B_q->ll::BR@Untagged',
            eos.Parameters(),
            eos.Kinematics(),
            eos.Options(model='WilsonScan', q='s', l='mu')
    )
    display(observable1)
    print('----')
    observable2 = eos.Observable.make('B->D^*lnu::A_FB',
            eos.Parameters(),
            eos.Kinematics(q2_min=0.02, q2_max=10.67),
            eos.Options(l='mu')
    )
    display(observable2)


.. list-table::
   :widths: 25, 25, 25

   * - B_q->ll::BR@Untagged
     - (eos.Observable)
     -
   * - kinematics
     - none
     -
   * -
     - l
     - mu
   * - options
     - model
     - WilsonScan
   * -
     - q
     - s
   * - current value
     - 3.567e-09
     -




.. parsed-literal::

    ----



.. list-table::
   :widths: 25, 25, 25

   * - B->D^*lnu::A_FB
     - (eos.Observable)
     -
   * - kinematics
     - q2_min
     - 0.02
   * -
     - q2_mmax
     - 10.67
   * - options
     - l
     - mu
   * - current value
     - 0.1674
     -



You can access an observable’s set of parameters through the
``parameters`` method:

.. code:: ipython3

    observable1.parameters() == observable2.parameters()




.. parsed-literal::

    False



As you can see, the two observables do not share a common set of
parameters, even though all their parameter values are identical.
Changes to the parameters of ``observable1`` do not affect
``observable2``.

You can also access the set of kinematic variables through the
``kinematics`` method:

.. code:: ipython3

    observable2.kinematics()



.. list-table::
   :widths: 25, 25

   * - q2_min
     - 0.02
   * - q2_mmax
     - 10.67




Finally, you can access the set of options with which an
``eos.Observable`` has been constructed through the ``options`` method:

.. code:: ipython3

    observable1.options()



.. list-table::
   :widths: 25, 25

   * - l
     - mu
   * - model
     - WilsonScan
   * - q
     - s



An observable can be handled just like any other Python object. For
example, we can readily create a list of observables that differ only by
one of their kinematic variables, e.g. to plot an observables as a
function of the kinematic variable.

.. code:: ipython3

    import numpy

    parameters  = eos.Parameters()
    observables = [
        eos.Observable.make('B->D^*lnu::A_FB(q2)', parameters, eos.Kinematics(q2=q2), eos.Options())
        for q2 in numpy.linspace(1.00, 10.67, 10)
    ]

    values = [o.evaluate() for o in observables]
    display(values)



.. parsed-literal::

    [0.07568167619600394,
     0.14540687361030025,
     0.18550212944207703,
     0.20752414012171227,
     0.2160779168296866,
     0.21324089622890918,
     0.19946563045552954,
     0.17338309618702236,
     0.12940814736260464,
     0.017723854042111166]


With this knowledge, you can now proceed to look at the example use
cases in the documentation and the other interactive notebook examples.
