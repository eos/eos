====================
Analysis File Format
====================

The `user guide <../user-guide/index.html>`_ illustrates various use cases for EOS that mostly cover relatively simple phenomenological analyses.
The number of observables, parameters, and statistical constraints is rather small.
Typical analyses carried out using EOS require a high degree of organisation,
as discussed in the `example on analysis organisation <../user-guide/analysis-organisation.html>`_.
This example highlights that EOS provides the means to store one or more analysis within an external text file file,
a so-called analysis file. Using this feature is beneficial in particular when the several analyses share common elements.
Every EOS analysis file is a YAML file that define the individual steps of one or more Bayesian analysis.
At the top level, the format includes the following YAML keys:

 - ``priors`` (**mandatory**) --- The list of priors within the analysis.

 - ``likelihoods`` (**mandatory**) --- The list of likelihoods within the analysis..

 - ``posteriors`` (**mandatory**) --- The list of posteriors within the analysis.

 - ``observables`` (**optional**) --- The list of custom observables defined for the scope of the analysis.

 - ``parameters`` (**optional**) --- The list of new parameters defined for the scope of the analysis.

 - ``predictions`` (**optional**) --- The list of theory predictions within the analysis.

The following example illustrates the analysis file format at the hand of a real-world example.

.. toggle-header::
   :header: Example `examples/b-to-u-l-nu.yaml <https://github.com/eos/eos/tree/master/examples/b-to-u-l-nu.yaml>`_

   .. literalinclude:: ../../examples/b-to-u-l-nu.yaml
      :language: YAML


Priors
~~~~~~

The ``priors`` key contains a list of *named* priors. Each prior has one mandatory key:

  - ``name`` (**mandatory**) --- The unique name of this prior.
  - ``descriptions`` (**mandatory**) --- The ordered list of parameters described by this prior.

For a univariate prior, ``descriptions`` contains only a single prior description, which is a key/value map describing the prior.
For an uncorrelated multivariate prior, ``descriptions`` contains more than one prior description.
The format for a prior description is the same as used in the :class:`Analysis <eos.Analysis>` constructor.
For example, the following code create two named priors ``CKM`` (univariate) and ``FF-pi`` (multivariate but uncorrelated).

.. code-block:: yaml

   priors:
     - name: CKM
       descriptions:
        - { 'parameter': 'CKM::abs(V_ub)', 'min': 3.0e-3, 'max': 4.0e-3, 'type': 'uniform' }

     - name: FF-pi
       descriptions:
         - { 'parameter':  'B->pi::f_+(0)@BCL2008' , 'min':   0.21 , 'max':   0.32 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_+^1@BCL2008'  , 'min':  -2.96 , 'max':  -0.60 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_+^2@BCL2008'  , 'min':  -3.98 , 'max':   4.38 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_+^3@BCL2008'  , 'min': -18.30 , 'max':   9.27 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_0^1@BCL2008'  , 'min':  -0.10 , 'max':   1.35 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_0^2@BCL2008'  , 'min':  -2.08 , 'max':   4.65 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_0^3@BCL2008'  , 'min':  -4.73 , 'max':   9.07 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_0^4@BCL2008'  , 'min': -60.00 , 'max':  38.00 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::f_T(0)@BCL2008' , 'min':   0.18 , 'max':   0.32 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_T^1@BCL2008'  , 'min':  -3.91 , 'max':  -0.33 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_T^2@BCL2008'  , 'min':  -4.32 , 'max':   2.00 , 'type': 'uniform' }
         - { 'parameter':  'B->pi::b_T^3@BCL2008'  , 'min':  -7.39 , 'max':  10.60 , 'type': 'uniform' }


For a correlated multivariate prior, ``descriptions`` contains a single element consisting of a key/value pair.
The mandatorys key is ``constraint`` and the value is the qualified name of an EOS constraint.
The following example illustrates the organisation of a correlated multivariate prior.

.. code-block:: yaml

   priors:
     - name: FF-rho
       descriptions:
         - { 'constraint': 'B->rho::FormFactors[parametric,LCSR]@BSZ:2015A' }


Likelihoods
~~~~~~~~~~~

The ``likelihoods`` key contains a list of *named* likelihoods. Each likelihood has two mandatory keys:

  - ``name`` (**mandatory**) --- The unique name of this likelihood.
  - ``constraints`` or ``manual_constraints`` or ``pyhf`` (**mandatory**) --- The ordered list of EOS objects that comprise this likelihood.

The likelihoods can be of three types:

  - ``constraints`` The following example illustrates the organisation of a likelihood for a simple constraint.

  .. code-block:: yaml

    - name: EXP-pi
      constraints:
        - 'B^0->pi^-l^+nu::BR@HFLAV:2019A;form-factors=BCL2008-4'

  - ``manual_constraints`` Additional manually-specified constraints can also be added.
    The syntax needs to follow the syntax of the usual ``EOS constraints``, as in the following example.

  .. code-block:: yaml

    - name: manual-TH-pi
      manual_constraints:
        "B->pi::form-factor-ratio":
          type: "Gaussian"
          observable: "B->pi::f_0(q2)/f_+(q2)"
          kinematics: {'q2': 0}
          options: {'form-factors': 'BSZ2015'}
          mean: 1
          sigma-stat: {"hi": 0., "lo": 0.}
          sigma-sys:  {"hi": 0.1, "lo": 0.1}

  - ``pyhf``

Posteriors
~~~~~~~~~~

The ``posteriors`` key contains a list of *named* posteriors. Each posterior contains two mandatory and various optional keys:

  - ``name`` (**mandatory**) --- The unique name of this posterior.
  - ``prior`` (**mandatory**) --- The ordered list of named priors that are used as part of this posterior.
  - ``likelihood`` (**mandatory**) --- The ordered list of named likelihoods that are used as part of this posterior.
  - ``global_options`` (**optional**) --- A key/value map providing global options, i.e., options that apply to all observables used within this posterior.
  - ``fixed_parameters`` (**optional**) --- A key/value map providing values for parameters that deviate from the default values.

The following example illustrates the organisation of a posterior.

.. code-block:: yaml

   posteriors:
     - name: CKM-pi
       global_options:
         l: e
         model: CKM
       prior:
         - CKM
         - FF-pi
       likelihood:
         - TH-pi
         - EXP-pi


Observables
~~~~~~~~~~~

New observables can be defined and used in the analysis description by following the syntax described in :ref:`the corresponding section <defining_observables>`.

For example, the following code defines the ratio of two :math:`B \to \pi` form-factors as a new observable.

.. code-block:: yaml

  observables:
    'B->pi::f_+(q2)/f_0(q2)':
      latex: '$\frac{f_+}{f_0}$'
      unit: '1'
      options: {}
      expression:
        '<<B->pi::f_+(q2)>> / <<B->pi::f_0(q2)>>'


Parameters
~~~~~~~~~~~

New parameters can also be defined in the analysis description. This can be useful in two cases:

  1. The new parameter(s) can be directly used in a custom observable and added to the analysis priors.
  The combination of new observables, parameters and manual constraints make ``EOS`` extremely flexible.
  The syntax for a new parameter follows:

  .. code-block:: yaml

    parameters:
      'prefix::name' :
          central: +1.0
          min:     +0.0
          max:     +2.0
          unit:     '1'
          latex:    '$p_\mathrm{user}$'

  2. New parameters can also be used as aliases for existing parameters. Varying the alias will then vary all the aliased parameters.
  This is particularly useful in analyses that assumes some symmetry amongst the parameters.
  E.g. for a fit to Wilson coefficients under the assumption of lepton flavor universality, we can use

  .. code-block:: yaml

    parameters:
      'ublnul::Re{cVL}' :
        alias_of: [ 'ubenue::Re{cVL}', 'ubmunumu::Re{cVL}', 'ubtaunutau::Re{cVL}' ]
        central: +1.0
        min:     +0.0
        max:     +2.0
        unit:     '1'
        latex:    '$\mathrm{Re}\, \mathcal{C}^{\bar{u}b\bar{\ell}\nu_\ell}_{V_L}$'


Predictions
~~~~~~~~~~~

The last step of an analysis usually consists in the prediction of a set of observables based on previously obtained importance samples.
The recognized ``predictions`` keys are:

  - ``name`` (**mandatory**) The name of the set of predictions.
  - ``observables`` (**mandatory**) The list of observables that need to be predicted. This should contain valid existing or manually-specified observables.
  - ``global_options`` (**optional**) The global options that should be used in the evaluation of the observables.
  - ``fixed_parameters`` (**optional**) A dictionary of parameters and their values that will be fixed in the evaluation of the observables.

The observables accept two keys:
  - ``name`` (**mandatory**) The qualified name of the observable.
    Options can be specified in the observable name following the syntax of :class:`QualifiedName <eos.QualifiedName>`.
    A warning will be raised if the observable option override the global options defined above.
  - ``kinematics`` (**optional**) The dictionary of kinematics specifications for the observables.
    For brevety, a list of kinematic specification can be provided. In this case, one observable per specification will be created.

The following code provides a valid example of predictions.

.. code-block:: yaml

  predictions:
  - name: BR
    global_options:
      model: CKM
    observables:
      - name: B_u->lnu::BR;l=e
      - name: B_u->lnu::BR;l=mu
      - name: B_u->lnu::BR;l=tau

  - name: dBR
    global_options:
      l: e
      q: d
      model: CKM
      form-factors: BCL2008
    observables:
      - name: B->pilnu::dBR/dq2
        kinematics: [ { q2:  1.0 }, { q2:  2.0 }, { q2:  3.0 }, { q2:  4.0 }, { q2:  5.0 },
                      { q2:  6.0 }, { q2:  7.0 }, { q2:  8.0 }, { q2:  9.0 }, { q2: 10.0 },
                      { q2: 11.0 }, { q2: 12.0 }, { q2: 13.0 }, { q2: 14.0 }, { q2: 15.0 },
                      { q2: 16.0 }, { q2: 17.0 }, { q2: 18.0 }, { q2: 19.0 }, { q2: 20.0 },
                      { q2: 21.0 }, { q2: 22.0 }, { q2: 23.0 }, { q2: 24.0 }, { q2: 25.0 },
                      { q2: 26.0 }, { q2: 27.0 } ]
