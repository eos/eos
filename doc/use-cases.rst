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

.. image:: /images/use-cases_plot-example.png
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
   eos.plot.Plotter(plot_args, '/home/dvandyk/Repositories/eos/master/doc/images/use-cases_hist-b-to-d-l-nu.png').plot()

.. image:: /images/use-cases_hist-b-to-d-l-nu.png
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

.. image:: /images/use-cases_hist-r-d.png
   :width: 600

Numerically, we obtain:

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

To obtain uncertainty bands for a plot of the differential branching ratios, we can now run:

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

We can plot these samples with EOS' plotting framework by running:

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

.. image:: /images/use-cases_plot-uncertainty.png
   :width: 600

*******************
Parameter Inference
*******************

.. todo::

  Write section on parameter inference

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

.. image:: /images/use-cases_inference_plot-a-priori.png
   :width: 600

.. code-block::

   analysis_args = {
       'global_options': { 'form-factors': 'BSZ2015', 'model': 'CKMScan' },
       'priors': [
           { 'parameter': 'CKM::abs(V_cb)',           'min':  38e-3, 'max':  45e-3,  'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f+_0@BSZ2015', 'min':  0.0,   'max':  1.0,    'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f+_1@BSZ2015', 'min': -4.0,   'max': -1.0,    'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f+_2@BSZ2015', 'min': +4.0,   'max': +6.0,    'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f0_1@BSZ2015', 'min': -1.0,   'max': +2.0,    'type': 'uniform' },
           { 'parameter': 'B->D::alpha^f0_2@BSZ2015', 'min': -2.0,   'max':  0.0,    'type': 'uniform' }
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

.. code-block::

   e_q2values  = np.unique(np.concatenate((np.linspace(0.02,  1.00, 20), np.linspace(1.00, 11.60, 20))))
   e_obs       = [eos.Observable.make(
                     'B->Dlnu::dBR/dq2', analysis.parameters, eos.Kinematics(q2=q2),
                     eos.Options(**{'form-factors': 'BSZ2015', 'l': 'e', 'q': 'd'}))
                 for q2 in e_q2values]
   parameter_samples, log_weights, e_samples  = analysis.sample(N=20000, stride=5, pre_N=1000, preruns=5, cov_scale=0.05, start_point=bfp.point, observables=e_obs)


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

.. image:: /images/use-cases_inference_plot-a-posteriori.png
   :width: 600

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

.. image:: /images/use-cases_inference_hist-Vcb.png
   :width: 600

**************************
Production of Peudo Events
**************************

.. todo::

  Write section on production of pseudo events
