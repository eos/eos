====================
Analysis File Format
====================

The `user guide <../user-guide/index.html>`_ illustrates various use cases for EOS that mostly cover relatively simple phenomenological analyses.
The number of observables, parameters, and statistical constraints is rather small.
Typical analyses carried out using EOS require a high degree of organisation,
as discussed in the `example on analysis organisation <../user-guide/analysis-organisation.html>`_.
This example highlights that EOS provides the means to store one or more analysis within an external text file file,
a so-called analysis file. Using this feature is beneficial in particular when the several analyses share common elements.
Every EOS analysis file is a YAML file that define the individual steps of one or more Bayesina analysis.
At the top level, the format includes the following YAML keys:

 - ``priors`` (**mandatory**) --- The list of priors within the analysis.

 - ``likelihoods`` (**mandatory**) --- The list of likelihoods within the analysis..

 - ``posteriors`` (**mandatory**) --- The list of posteriors within the analysis.

 - ``predictions`` (**optional**) --- The list of theory predictions within the analysis.

 - ``observables`` (**optional**) --- The list of custom observables defined for the scope of the analysis.

The following example illustrates the analysis file format at the hand of a real-world example.

.. toggle-header::
   :header: Example `examples/b-to-u-l-nu.yaml <https://github.com/eos/eos/tree/master/examples/b-to-u-l-nu.yaml>`_

   .. literalinclude:: ../../examples/b-to-u-l-nu.yaml
      :language: YAML


Priors
~~~~~~

The ``priors`` key contains a list of *named* priors. Each prior has one mandatory key:

  - ``name`` (**mandatory**) --- The unique name of this prior.
  - ``parameters`` (**mandatory**) --- The ordered list of parameters described by this prior.

For a univariate prior, ``parameters`` contains only a single prior description, which is a key/value map describing the prior.
For an uncorrelated multivariate prior, ``parameters`` contains more than one prior description.
The format for a prior description is the same as used in the :class:`Analysis <eos.Analysis>` constructor.
For example, the following code create two named priors ``CKM`` (univariate) and ``FF-pi`` (multivariate but uncorrelated).

.. code-block:: yaml

   priors:
     - name: CKM
       parameters:
        - { 'parameter': 'CKM::abs(V_ub)', 'min': 3.0e-3, 'max': 4.0e-3, 'type': 'uniform' }

     - name: FF-pi
       parameters:
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


For a correlated multivariate prior, ``parameters`` contains a single element consisting of a key/value pair.
The mandatorys key is ``constraint`` and the value is the qualified name of an EOS constraint.
The following example illustrates the organisation of a correlated multivariate prior.

.. code-block:: yaml

   priors:
     - name: FF-rho
       parameters:
         - { 'constraint': 'B->rho::FormFactors[parametric,LCSR]@BSZ:2015A' }


Likelihoods
~~~~~~~~~~~

The ``likelihoods`` key contains a list of *named* likelihoods. Each likelihood has two mandatory keys:

  - ``name`` (**mandatory**) --- The unique name of this likelihood.
  - ``constraints`` (**mandatory**) --- The ordered list of EOS constraint names that comprise this likelihood.

The following example illustrates the organisation of a likelihood.

.. code-block:: yaml

  - name: EXP-pi
    constraints:
      - 'B^0->pi^-l^+nu::BR@HFLAV:2019A;form-factors=BCL2008-4'

Posteriors
~~~~~~~~~~

The ``posteriors`` key contains a list of *named* posteriors. Each posterior contains two mandator and various optional keys:

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
