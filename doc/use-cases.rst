#########
Use Cases
#########


******************************************
Theory Predictions and their Uncertainties
******************************************

EOS can produce theory predictions for any of its built-in observables.
The examples following in this section illustrate how to find a specific observable from the list of all built-in observables,
construct an :class:`Observable <eos.Observable>` object and evaluate it,
and estimate the theoretical uncertainties associated with it.

.. note::

   The following examples are also available as an interactive `Jupyter <https://jupyter.org/>`_ notebook `from here <https://github.com/eos/eos/blob/master/examples/predictions.ipynb>`_.


Listing the built-in Observables
================================

The full list of built-in observables for the most-recent EOS release is available online `here <https://eos.github.io/doc/observables>`_.
For an interactive apporach to see the list of observables available to you, run the following in a Jupyter notebook:

.. code-block::

   import eos
   display(eos.Observables())

Searching for a specific observable is possible filtering by the observable name's `prefix`, `name`, or `suffix`, e.g.:

.. code-block::

   display(eos.Observables(name='A_FB'))
   display(eos.Observables(prefix='B->pi'))
   display(eos.Observables(suffix='Untagged'))



Constructing and Evaluating an Observable
=========================================

To make theory predictions of any observable, EOS requires its full name, its :class:`Parameters <eos.Parameters>` object,
its :class:`Kinematics <eos.Kinematics>` object, and its :class:`Options <eos.Options>` object.
As an example, we will use the integrated branching ratio of :math:`\bar{B}\to D\ell^-\bar\nu`,
which is represented by the :class:`QualifiedName <eos.QualifiedName>` ``B->Dlnu::BR``.
Additional information about any given observable can be obtained by displaying the full database entry,
which also contains information about the kinematic variables required:

.. code-block::

  display(eos.Observables()['B->pilnu::BR'])

which outputs the following table:

.. list-table::
   :widths: 25, 25

   * - Qualified Name
     - ``B->Dlnu::BR``
   * - Description
     - :math:`\mathcal{B}(\bar{B}\to D\ell^-\bar\nu)`
   * - Kinematic Variable
     - ``q2_min``, ``q2_max``

From this we understand that ``B->Dlnu::BR`` expects two kinematic variables,
corresponding here to the lower and upper integration boundaries of the dilepton invariant mass :math:`q^2`.

To create an :class:`Obervable <eos.Observable>` object with the default set of parameters and options through:
.. code-block::

   params = eos.Parameters.Defaults()
   kinematics = eos.Kinematics(q2_min=0.02, q2_max=10)
   obs = eos.Observable.make('B->Dlnu::BR', params, kinematics, eos.Options())
   display(obs)

The default options select a spectator :math:`\ell=\mu` yielding a value of :math:`2.3\%`,
which is compatible with the current world average for the :math:`\bar{B}^-\to D^0\mu^-\bar\nu` branching ratio.

.. code-block::

   parameters = eos.Parameters.Defaults()
   kinematics = eos.Kinematics(q2_min=3.17, q2_max=11.60)
   obs = eos.Observable.make('B->Dlnu::BR', parameters, kinematics, eos.Options(l='tau'))
   display(obs)

By setting the ``l`` option to the value ``'tau'``, we have create a different observable representing the :math:`\bar{B}^-\to D^0\tau^-\bar\nu` branching ratio.
The new observable yields a value of :math:`0.69\%`.

So far we evaluated the integrated branching ratio. EOS also provides the corresponding differential branching ratio as a function of :math:`q^2`.
It is accessible through the name ``B->Dlnu::dBR/dq2``. To illustrate it, we use EOS's plot functions:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [0.0, 11.60] },
           'y': { 'label': r'$d\mathcal{B}/dq^2$',                    'range': [0.0,  5e-3] },
           'legend': { 'location': 'upper center' }
       },
       'contents': [
           {
               'label': r'$\ell=\mu$',
               'type': 'observable',
               'observable': 'B->Dlnu::dBR/dq2;l=mu',
               'kinematic': 'q2',
               'range': [0.02, 11.60],
           },
           {
               'label': r'$\ell=\tau$',
               'type': 'observable',
               'observable': 'B->Dlnu::dBR/dq2;l=tau',
               'kinematic': 'q2',
               'range': [3.17, 11.60],
           }
       ]
   }
   eos.plot.Plotter(plot_args).plot()

which yields:

.. image:: /images/use-cases_prediction_plot-example.png
   :width: 600


Estimating Theoretical Uncertainties
====================================

To estimate theoretical uncertainties of the observables, EOS uses Bayesian statistics.
The latter interprets the theory parameters as random variables and assigns *a priori* probability density functions (prior PDFs) for each parameter.

.. note::

  For technical reasons, EOS can only use uncorrelated prior PDFs.
  The same effects a having correlated prior PDFs can be achieved by using a correlated likelihood and uniform prior PDFs.

We carry on using the integrated branching ratios of :math:`\bar{B}^-\to D^0\left\lbrace\mu^-, \tau^-\right\rbrace\bar\nu` decays as examples.
The largest source of theoretical uncertainty in these decays arises from the hadronic matrix elements, i.e.,
from the form factors :math:`f^{B\to \bar{D}}_+(q^2)` and :math:`f^{B\to \bar{D}}_0(q^2)`.
Both form factors have been obtained independently using lattice QCD simulations by the HPQCD and Fermilab/MILC (FNALMILC) collaborations.
The joint likelihoods for both form factors at different :math:`q^2` values of each experiment are available in EOS under the names ``B->D::f_++f_0@HPQCD2015A`` and ``B->D::f_++f_0@FNALMILC2015A``.
For this example, we will use both results and create a combined likelihood:

.. code-block::

   analysis_args = {
       'global_options': None,
       'priors': [
           { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0, 'max':  1.0, 'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f+_1@BSZ2015', 'min': -5.0, 'max': +5.0, 'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f+_2@BSZ2015', 'min': -5.0, 'max': +5.0, 'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f0_1@BSZ2015', 'min': -5.0, 'max': +5.0, 'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f0_2@BSZ2015', 'min': -5.0, 'max': +5.0, 'type': 'uniform' }
       ],
       'likelihood': [
           'B->D::f_++f_0@HPQCD2015A'
           'B->D::f_++f_0@FNALMILC2015A'
       ]
   }
   analysis = eos.Analysis(**analysis_args)

Next we create three observables: the semi-muonic branching ratio, the semi-tauonic branching ratio, and the ratio of the former two.
By using :code:`analysis.parameter` we ensure that all observables and the Analysis object share the same parameter set.


.. code-block::

   obs_mu  = eos.Observable.make(
       'B->Dlnu::BR',
       analysis.parameters,
       eos.Kinematics(q2_min=0.02, q2_max=11.60),
       eos.Options(**{'l':'mu', 'form-factors':'BSZ2015'})
   )
   obs_tau = eos.Observable.make(
       'B->Dlnu::BR',
       analysis.parameters,
       eos.Kinematics(q2_min=3.17, q2_max=11.60),
       eos.Options(**{'l':'tau','form-factors':'BSZ2015'})
   )
   obs_R_D = eos.Observable.make(
       'B->Dlnu::R_D',
       analysis.parameters,
       eos.Kinematics(q2_mu_min=0.02, q2_mu_max=11.60, q2_tau_min=3.17, q2_tau_max=11.60),
       eos.Options(**{'form-factors':'BSZ2015'})
   )
   observables=(obs_mu, obs_tau, obs_R_D)

In the above, we made sure to use :code:`form-factors=BSZ2015` to ensure that the right form factor plugin is used.


Sampling from the log(posterior) and -- at the same time -- producing posterior-predictive samples of the :code:`observables` is achieved by running:

.. code-block:

   parameter_samples, log_weights, observable_samples = analysis.sample(N=5000, pre_N=1000, observables=observables)

Here :code:`N=5000` samples are produced. To illustrate these samples we use EOS' plotting framework:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$d\mathcal{B}/dq^2$',  'range': [0.0,  3e-2] },
           'legend': { 'location': 'upper center' }
       },
       'contents': [
           { 'label': r'$\ell=\mu$', 'type': 'histogram', 'bins': 30, 'data': { 'samples': observable_samples[:, 0], 'log_weights': log_weights }},
           { 'label': r'$\ell=\tau$','type': 'histogram', 'bins': 30, 'data': { 'samples': observable_samples[:, 1], 'log_weights': log_weights }},
       ]
   }
   eos.plot.Plotter(plot_args).plot()

.. image:: /images/use-cases_prediction_hist-b-to-d-l-nu.png
   :width: 600

and

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$d\mathcal{B}/dq^2$',  'range': [0.28,  0.32] },
           'legend': { 'location': 'upper left' }
       },
       'contents': [
           { 'label': r'$R_D$ (EOS)',     'type': 'histogram', 'bins': 30, 'color': 'C3', 'data': { 'samples': observable_samples[:, 2], 'log_weights': log_weights }},
           { 'label': r'$R_D$ (manually)','type': 'histogram', 'bins': 30, 'color': 'C4', 'data': { 'samples': [o[1] / o[0] for o in observable_samples[:]], 'log_weights': log_weights },
             'histtype': 'step'},
       ]
   }
   eos.plot.Plotter(plot_args).plot()

.. image:: /images/use-cases_prediction_hist-r-d.png
   :width: 600

Using the Numpy routines :code:`numpy.average` and :code:`numpy.var` we can produce numerical estimates
of the weighted mean and its standard deviation:

.. code-block::

   import numpy as np
   print('{obs};{opt}  = {mean:.4f} +/- {std:.4f}'.format(
       obs=obs_mu.name(), opt=obs_mu.options(),
       mean=np.average(observable_samples[:,0], weights=np.exp(log_weights)),
       std=np.sqrt(np.var(observable_samples[:, 0]))
   ))
   print('{obs};{opt} = {mean:.4f} +/- {std:.4f}'.format(
       obs=obs_tau.name(), opt=obs_tau.options(),
       mean=np.average(observable_samples[:,1], weights=np.exp(log_weights)),
       std=np.sqrt(np.var(observable_samples[:, 1]))
   ))
   print('{obs};{opt}      = {mean:.4f} +/- {std:.4f}'.format(
       obs=obs_R_D.name(), opt=obs_R_D.options(),
       mean=np.average(observable_samples[:,2], weights=np.exp(log_weights)),
       std=np.sqrt(np.var(observable_samples[:, 1]))
   ))

From the above we obtain:

.. code-block::

   B->Dlnu::BR;form-factors=BSZ2015,l=mu  = 0.0231 +/- 0.0007
   B->Dlnu::BR;form-factors=BSZ2015,l=tau = 0.0070 +/- 0.0001
   B->Dlnu::R_D;form-factors=BSZ2015      = 0.3020 +/- 0.0001

To obtain uncertainty bands for a plot of the differential branching ratios, we can now produce a
sequence of observables at different points in phase space. We then pass these observables on to
:method:`analysis.sample <eos.Analysis.sample>`, to obtain posterior-predictive samples:

.. code-block:

   mu_q2values  = np.unique(np.concatenate((np.linspace(0.02,  1.00, 20), np.linspace(1.00, 11.60, 20))))
   mu_obs       = [eos.Observable.make(
                      'B->Dlnu::dBR/dq2', analysis.parameters, eos.Kinematics(q2=q2),
                      eos.Options(**{'form-factors': 'BSZ2015', 'l': 'mu'}))
                  for q2 in mu_q2values]
   tau_q2values = np.linspace(3.17, 11.60, 40)
   tau_obs      = [eos.Observable.make(
                      'B->Dlnu::dBR/dq2', analysis.parameters, eos.Kinematics(q2=q2),
                      eos.Options(**{'form-factors': 'BSZ2015', 'l': 'tau'}))
                  for q2 in tau_q2values]
   _, log_weights, mu_samples  = analysis.sample(N=5000, pre_N=1000, observables=mu_obs)
   _, log_weights, tau_samples = analysis.sample(N=5000, pre_N=1000, observables=tau_obs)

We can plot the so-obtained posterior-predictive samples with EOS' plotting framework by running:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [0.0, 11.60] },
           'y': { 'label': r'$d\mathcal{B}/dq^2$',                    'range': [0.0,  5e-3] },
           'legend': { 'location': 'upper center' }
       },
       'contents': [
           {
             'label': r'$\ell=\mu$', 'type': 'uncertainty', 'range': [0.02, 11.60],
             'data': { 'samples': mu_samples, 'xvalues': mu_q2values }
           },
           {
             'label': r'$\ell=\tau$','type': 'uncertainty', 'range': [3.17, 11.60],
             'data': { 'samples': tau_samples, 'xvalues': tau_q2values }
           },
       ]
   }
   eos.plot.Plotter(plot_args).plot()

.. image:: /images/use-cases_prediction_plot-uncertainty.png
   :width: 600


*******************
Parameter Inference
*******************

EOS can infer parameters based on a database of experimental or theoetical constraints and its built-in observables.
The examples following in this section illustrate how to find a specific constraint from the list of all built-in observables,
construct an :class:`Analysis <eos.Analysis>` object that represents the statistical analysis,
and infer mean value and standard deviation of a list of parameters.

.. note::

   The following examples are also available as an interactive `Jupyter <https://jupyter.org/>`_ notebook `from here <https://github.com/eos/eos/blob/master/examples/inference.ipynb>`_.


Listing the built-in Constraints
================================

The full list of built-in constraints for the most-recent EOS release is available online `here <https://eos.github.io/doc/constraints>`_.
For an interactive apporach to see the list of constraints available to you, run the following in a Jupyter notebook:

.. code-block::

   import eos
   display(eos.Constraints())

Searching for a specific observable is possible filtering by the constraint name's `prefix`, `name`, or `suffix`, e.g.:

.. code-block::

   display(eos.Constraints(prefix='B^0->D^+'))


Visualizing the built-in Constraints
====================================

For what follows we will use the two experimental constraints ``B^0->D^+e^-nu::BRs@Belle-2015A`` and ``B^0->D^+mu^-nu::BRs@Belle-2015A``,
to infer the CKM matrix element :math:`|V_{cb}|`. We can readily display these two constraints, along with the default theory prediction,
using the following code:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [0.0, 11.63] },
           'y': { 'label': r'$d\mathcal{B}/dq^2$',                    'range': [0.0,  5e-3] },
           'legend': { 'location': 'lower left' }
       },
       'contents': [
           {
               'label': r'$\ell=e$',
               'type': 'observable',
               'observable': 'B->Dlnu::dBR/dq2;l=e,q=d',
               'kinematic': 'q2',
               'color': 'black',
               'range': [0.02, 11.63],
           },
           {
               'label': r'Belle 2015 $\ell=e,\, q=d$',
               'type': 'constraint',
               'color': 'C0',
               'constraints': 'B^0->D^+e^-nu::BRs@Belle-2015A',
               'observable': 'B->Dlnu::BR',
               'variable': 'q2',
               'rescale-by-width': False
           },
           {
               'label': r'Belle 2015 $\ell=\mu,\,q=d$',
               'type': 'constraint',
               'color': 'C1',
               'constraints': 'B^0->D^+mu^-nu::BRs@Belle-2015A',
               'observable': 'B->Dlnu::BR',
               'variable': 'q2',
               'rescale-by-width': False
           },
       ]
   }
   eos.plot.Plotter(plot_args).plot()

The resulting plot looks like this:

.. image:: /images/use-cases_inference_plot-a-priori.png
   :width: 600


Defining the statistical Analysis
=================================

To define our statistical analysis for the inference of :math:`|V_{cb}|` from :math:`\bar{B}\to D\ell^-\bar\nu` branching ratios,
some decisions are needed.
First, we must decide how to parametrize the hadronic form factors that emerge in semileptonic :math:`\bar{B}\to D` transitions.
For what follows we will use the [BSZ2015]_ parametrization.
Next, we must decide the theory input for the form factors. For what follows we will combine the correlated lattice QCD results published by the Fermilab/MILC and HPQCD collaborations in 2015.

We then create an :class:`Analysis <eos.Analysis>` object as follows:

.. code-block::

   analysis_args = {
       'global_options': { 'form-factors': 'BSZ2015', 'model': 'CKMScan' },
       'priors': [
           { 'parameter': 'CKM::abs(V_cb)',           'min':  38e-3, 'max':  45e-3, },
           { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0,   'max':  1.0,   },
           { 'parameter': 'B->D::alpha^f+_1@BSZ2015', 'min': -4.0,   'max': -1.0,   },
           { 'parameter': 'B->D::alpha^f+_2@BSZ2015', 'min': +4.0,   'max': +6.0,   },
           { 'parameter': 'B->D::alpha^f0_1@BSZ2015', 'min': -1.0,   'max': +2.0,   },
           { 'parameter': 'B->D::alpha^f0_2@BSZ2015', 'min': -2.0,   'max':  0.0,   }
       ],
       'likelihood': [
           'B->D::f_++f_0@HPQCD2015A',
           'B->D::f_++f_0@FNALMILC2015A',
           'B^0->D^+e^-nu::BRs@Belle-2015A',
           'B^0->D^+mu^-nu::BRs@Belle-2015A'
       ]
   }
   analysis = eos.Analysis(**analysis_args)
   analysis.parameters['CKM::abs(V_cb)'].set(42.0e-3)

In the above, the global options ensure that our choice of form factor parametrization is used throughout,
and that for CKM matrix elements the `CKMScan` model is used. The latter provides access to :math:`V_{cb}` matrix element through two :class:`parameters <eos.Parameter>`:
the absolute value ``CKM::abs(V_cb)`` and the complex phase ``CKM::arg(V_cb)``.
We provide the parameters in our analysis through the specifications of the Bayesian priors.
In the above, each prior is by default a uniform prior that covers the range from ``min`` to ``max``.
The likelihood is defined through a list constraints, which in the above includes both the experimental measurements by the Belle collaboration as well as the theoretical lattice QCD results.
Finally, we set the starting value of ``CKM::abs(V_cb)`` to a sensible value of :math:`42\cdot 10^{-3}`.

We can now proceed to optimize the log(posterior) through a call to ``analysis.optimize``.
In a Jupyter notebook, it is useful to display the return value of this method, which illustrates the best-fit point.
We can further display a summary of the goodness-of-fit information.

.. code-block::

   bfp = analysis.optimize()
   display(bfp)
   display(analysis.goodness_of_fit())

The resulting best-fit point looks like this:

.. list-table::
   :widths: 25, 25

   * - parameter
     - value
   * - :math:`|V_{cb}|`
     - 0.0422
   * - ``B->D::alpha^f+_0@BSZ2015``
     - 0.6671
   * - ``B->D::alpha^f+_1@BSZ2015``
     - -2.5314
   * - ``B->D::alpha^f+_2@BSZ2015``
     - 4.8813
   * - ``B->D::alpha^f0_1@BSZ2015``
     - 0.2660
   * - ``B->D::alpha^f0_2@BSZ2015``
     - -0.8410

The goodness-of-fit summary consists of a table listing all constraints,

.. list-table::
   :widths: 25, 25, 25

   * - constraint
     - :math:`\chi^2`
     - d.o.f.
   * - ``B->D::f_++f_0@FNALMILC2015A``
     - 3.4847
     - 7
   * - ``B->D::f_++f_0@HPQCD2015A``
     - 3.1016
     - 5
   * - ``B^0->D^+e^-nu::BRs@Belle-2015A``
     - 11.8206
     - 10
   * - ``B^0->D^+mu^-nu::BRs@Belle-2015A``
     - 5.2242
     - 10

and the overall information including the p value:

.. list-table::
   :widths: 25, 25

   * - total :math:`\chi^2`
     - 23.6310
   * - total degrees of freedom
     - 26
   * - p-value
     - 59.7053%


Sampling from the posterior
===========================

To sample from the posterior, EOS provides the :meth:`sample <eos.Analysis.sample>` method.
Optionally, it can also produce posterior-predictive samples for a list of observables.
For this example, we produce such samples for the differential :math:`\bar{B}\to D^+e^-\bar\nu` branching ratio in 40 points in the phase space.
We can use this to illustrate the results of our fit in relation to the constraints.
The call is:

.. code-block::

   e_q2values  = np.unique(np.concatenate((np.linspace(0.02,  1.00, 20), np.linspace(1.00, 11.60, 20))))
   e_obs       = [eos.Observable.make(
                     'B->Dlnu::dBR/dq2', analysis.parameters, eos.Kinematics(q2=q2),
                     eos.Options(**{'form-factors': 'BSZ2015', 'l': 'e', 'q': 'd'}))
                 for q2 in e_q2values]
   parameter_samples, log_weights, e_samples  = analysis.sample(N=20000, stride=5, pre_N=1000, preruns=5, cov_scale=0.05, start_point=bfp.point, observables=e_obs)

In the above we start sampling at the best-fit point, and carry out :code:`preruns = 5` burn-in runs/preruns of :code:`pre_N = 1000` samples each,
and the main run with a total of :code:`N * stride = 100000` random Markov Chain samples.
The latter are thinned down by a factor of :code:`stride = 5` to obtain :code:`N = 20000` samples in :code:`parameter_samples`.
The values of the log(posterior) are stored in :code:`log_weights`.
The posterior-preditive samples for the observables are stored in :code:`e_samples`, and are only returned if the :code:`observables` keyword argument is provided.


We can illustrate the posterior samples either as a histogram or as a kernel density estimate (KDE) using the built-in plotting functions:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$|V_{cb}|$', 'range': [38e-3, 45e-3] },
           'legend': { 'location': 'upper left' }
       },
       'contents': [
           {
               'type': 'histogram',
               'data': { 'samples': parameter_samples[:, 0], 'log_weights': log_weights }
           },
           {
               'type': 'kde', 'color': 'C0', 'label': 'posterior', 'bandwidth': 2,
               'range': [40e-3, 45e-3],
               'data': { 'samples': parameter_samples[:, 0], 'log_weights': log_weights }
           }
       ]
   }
   eos.plot.Plotter(plot_args).plot()

The result looks like this:

.. image:: /images/use-cases_inference_hist-Vcb.png
   :width: 600

Contours at given levels of posterior probability can be obtained for any pair of parameters using:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$|V_{cb}|$', 'range': [38e-3, 47e-3] },
           'y': { 'label': r'$f_+(0)$', 'range': [0.6, 0.75] },
       },
       'contents': [
           {
               'type': 'kde2D', 'color': 'C1', 'label': 'posterior',
               'range': [40e-3, 45e-3], 'levels': [68, 99], 'bandwidth': 3,
               'data': { 'samples': parameter_samples[:, (0,1)], 'log_weights': log_weights }
           }
       ]
   }
   eos.plot.Plotter(plot_args).plot()

The result looks like this:

.. image:: /images/use-cases_inference_hist-Vcb-f_+-zero.png
   :width: 600

We can visualize the posterior-predictive samples using:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [0.0, 11.63] },
           'y': { 'label': r'$d\mathcal{B}/dq^2$',                    'range': [0.0,  5e-3] },
           'legend': { 'location': 'lower left' }
       },
       'contents': [
           {
             'label': r'$\ell=\mu$', 'type': 'uncertainty', 'range': [0.02, 11.60],
             'data': { 'samples': e_samples, 'xvalues': e_q2values }
           },
           {
               'label': r'Belle 2015 $\ell=e,\, q=d$',
               'type': 'constraint',
               'color': 'C0',
               'constraints': 'B^0->D^+e^-nu::BRs@Belle-2015A',
               'observable': 'B->Dlnu::BR',
               'variable': 'q2',
               'rescale-by-width': False
           },
           {
               'label': r'Belle 2015 $\ell=\mu,\,q=d$',
               'type': 'constraint',
               'color': 'C1',
               'constraints': 'B^0->D^+mu^-nu::BRs@Belle-2015A',
               'observable': 'B->Dlnu::BR',
               'variable': 'q2',
               'rescale-by-width': False
           },
       ]
   }
   eos.plot.Plotter(plot_args).plot()

The result looks like this:

.. image:: /images/use-cases_inference_plot-a-posteriori.png
   :width: 600


.. image:: /images/use-cases_inference_hist-Vcb-f_+-zero.png
   :width: 600

**************************
Production of Peudo Events
**************************

.. todo::

  Write section on production of pseudo events
