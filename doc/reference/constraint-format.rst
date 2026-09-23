*****************
Constraint Format
*****************

Every entry of :class:`eos.Constraints` is described on disk by a YAML map whose ``type`` key selects one of the concrete field schemas below.
The full list of constraints currently shipped with EOS is documented separately on the :doc:`constraints` page.


``Amoroso``
~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observable``
     - scalar
     - yes
     - The qualified name of the observable this constraint applies to.
   * - ``kinematics``
     - map
     - yes
     - The kinematic variables at which the observable is evaluated.
   * - ``options``
     - map
     - yes
     - The options with which the observable is evaluated.
   * - ``physical-limit``
     - scalar
     - yes
     - The distribution's physical limit.
   * - ``theta``
     - scalar
     - yes
     - The distribution's scale parameter theta.
   * - ``alpha``
     - scalar
     - yes
     - The distribution's shape parameter alpha.
   * - ``beta``
     - scalar
     - yes
     - The distribution's shape parameter beta.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   B^0_s->mu^+mu^-::BR@CMS:2013B:
     type: Amoroso
     observable: B_q->ll::BR@Untagged
     kinematics: {}
     options: {l: mu, q: s}
     physical-limit: 0
     theta: 1.98596335e-09
     alpha: 2.7971996
     beta: 2.02188458
     references:
       []


``Gaussian``
~~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observable``
     - scalar
     - yes
     - The qualified name of the observable this constraint applies to.
   * - ``kinematics``
     - map
     - yes
     - The kinematic variables at which the observable is evaluated.
   * - ``options``
     - map
     - yes
     - The options with which the observable is evaluated.
   * - ``mean``
     - scalar
     - yes
     - The measurement's central value.
   * - ``sigma-stat``
     - map
     - yes
     - The measurement's statistical uncertainty, as 'hi'/'lo' one-sigma values.
   * - ``sigma-sys``
     - map
     - yes
     - The measurement's systematic uncertainty, as 'hi'/'lo' one-sigma values.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   B^0_s->mu^+mu^-::BR@CMS:2016A:
     type: Gaussian
     observable: B_q->ll::BR@Untagged
     kinematics: {}
     options: {l: mu, q: s}
     mean: 3e-09
     sigma-stat: {hi: 1e-09, lo: 9e-10}
     sigma-sys: {hi: 0, lo: 0}
     references:
       []


``LogGamma``
~~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observable``
     - scalar
     - yes
     - The qualified name of the observable this constraint applies to.
   * - ``kinematics``
     - map
     - yes
     - The kinematic variables at which the observable is evaluated.
   * - ``options``
     - map
     - yes
     - The options with which the observable is evaluated.
   * - ``mode``
     - scalar
     - yes
     - The distribution's mode.
   * - ``sigma``
     - map
     - yes
     - The distribution's width, as 'hi'/'lo' one-sigma values.
   * - ``alpha``
     - scalar
     - yes
     - The distribution's shape parameter alpha.
   * - ``lambda``
     - scalar
     - yes
     - The distribution's shape parameter lambda.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   # No constraint of type LogGamma currently ships with EOS; the following illustrates the format.
   Example::Constraint@LogGamma:2026A:
       type: LogGamma
       observable: mass::b(MSbar)
       kinematics: {}
       options: {}
       mode: 0.53
       sigma: { hi: 0.1, lo: 0.19 }
       alpha: 0.383056
       lambda: 0.0687907



``Mixture``
~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observables``
     - sequence
     - yes
     - The qualified names of the observables this constraint applies to.
   * - ``kinematics``
     - sequence
     - yes
     - One kinematics map per observable, in the same order as 'observables'.
   * - ``options``
     - sequence
     - yes
     - One options map per observable, in the same order as 'observables'.
   * - ``components``
     - sequence
     - yes
     - The mixture's components, each a map with keys 'means' and 'covariance'.
   * - ``weights``
     - sequence
     - yes
     - The mixture weight of each component, in the same order as 'components'.
   * - ``test statistics``
     - map
     - yes
     - The test-statistic map, with keys 'sigma' and 'densities'.
   * - ``dof``
     - scalar
     - no
     - The number of degrees of freedom; defaults to the dimension of the first component.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   ublnu::P(WET)@LMNRvD:2023A:
     type: Mixture
     observables:
       - ubenue::Re{cVL}
       - ubenue::Re{cVR}
       - ubenue::Re{cSL}
       - ubenue::Re{cSR}
       - ubenue::Re{cT}
     kinematics:
       - {}
       - {}
       - {}
       - {}
       - {}
     options:
       - {}
       - {}
       - {}
       - {}
       - {}
     components:
       - means: [0.239210427, 0.745753385, 0.0134320629, 0.00712165087, 0.07476569]
         covariance:
           - [0.00413002622, -0.00332498388, -0.000848109332, -0.000406691892, 0.00149330229]
           - [-0.00332498388, 0.00423865185, 0.000609995065, 0.000447884848, -0.00198342431]
           - [-0.000848109332, 0.000609995065, 0.0179180233, -0.0139526675, -0.000419354195]
           - [-0.000406691892, 0.000447884848, -0.0139526675, 0.0176090923, -0.000427770435]
           - [0.00149330229, -0.00198342431, -0.000419354195, -0.000427770435, 0.00201519131]
       - means: [0.757398556, 0.228773992, 0.00110913152, -0.000873656294, 0.0677685293]
         covariance:
           - [0.00372855018, -0.00283559277, 7.81683737e-05, -5.53574956e-05, -0.00175599948]
           - [-0.00283559277, 0.00370910325, -0.00013336092, 0.00012761393, 0.00128510668]
           - [7.81683737e-05, -0.00013336092, 0.0177483896, -0.0130322168, 6.4885982e-05]
           - [-5.53574956e-05, 0.00012761393, -0.0130322168, 0.0169730206, -4.18798439e-05]
           - [-0.00175599948, 0.00128510668, 6.4885982e-05, -4.18798439e-05, 0.00195438375]
       - means: [0.198516518, 0.792744105, -0.0472557982, -0.0333198714, -0.00746000472]
         covariance:
           - [0.00233305594, -0.00112909094, 0.000503448472, 0.000338975635, -0.000152968596]
           - [-0.00112909094, 0.00192412867, 0.000601037901, 1.78341407e-05, 0.000272259322]
           - [0.000503448472, 0.000601037901, 0.0148029835, -0.0123754828, 0.000114676622]
           - [0.000338975635, 1.78341407e-05, -0.0123754828, 0.0141841404, -9.50121964e-05]
           - [-0.000152968596, 0.000272259322, 0.000114676622, -9.50121964e-05, 0.00307078657]
       - means: [0.55162879, 0.421223066, -0.00326568492, 0.0073358451, -0.136095987]
         covariance:
           - [0.0136726272, -0.0124196149, -6.77457689e-05, -0.000186221927, 0.00113994882]
           - [-0.0124196149, 0.0131990105, -0.0001525186, 0.000187918672, -0.00072572775]
           - [-6.77457689e-05, -0.0001525186, 0.0231210904, -0.0194229642, -0.000140483516]
           - [-0.000186221927, 0.000187918672, -0.0194229642, 0.0217373176, 7.90708268e-05]
           - [0.00113994882, -0.00072572775, -0.000140483516, 7.90708268e-05, 0.000671170943]
       - means: [0.5034653, 0.46872954, 0.00419595692, -0.00736434769, 0.138663722]
         covariance:
           - [0.0143097582, -0.0131789312, 5.84795881e-05, -0.000179578405, -0.000449735291]
           - [-0.0131789312, 0.0141062439, 0.000283478865, 0.00014711228, 5.31608868e-05]
           - [5.84795881e-05, 0.000283478865, 0.0235784881, -0.0198343226, -0.0001527965]
           - [-0.000179578405, 0.00014711228, -0.0198343226, 0.0220192126, 0.000118910175]
           - [-0.000449735291, 5.31608868e-05, -0.0001527965, 0.000118910175, 0.00060768473]
       - means: [0.767568627, 0.220278888, -0.00097655811, 0.000878976359, -0.057475313]
         covariance:
           - [0.00317257454, -0.00234444581, -6.06395131e-05, 3.61794251e-05, 0.00156711916]
           - [-0.00234444581, 0.00326005989, 0.00010849244, -9.86664846e-05, -0.00111646512]
           - [-6.06395131e-05, 0.00010849244, 0.0174746435, -0.012628733, 7.86544453e-05]
           - [3.61794251e-05, -9.86664846e-05, -0.012628733, 0.0166955983, -5.66194913e-05]
           - [0.00156711916, -0.00111646512, 7.86544453e-05, -5.66194913e-05, 0.00210933949]
       - means: [0.250297463, 0.733580913, 0.0125323274, 0.0115232794, -0.0826979928]
         covariance:
           - [0.0048791311, -0.00412718468, -0.0007456679, -0.000709134398, -0.00172573604]
           - [-0.00412718468, 0.00513665626, 0.000498154877, 0.000598726604, 0.00223295979]
           - [-0.0007456679, 0.000498154877, 0.0183029766, -0.0144801713, 0.000482037239]
           - [-0.000709134398, 0.000598726604, -0.0144801713, 0.0178354501, 0.000332905413]
           - [-0.00172573604, 0.00223295979, 0.000482037239, 0.000332905413, 0.00193181161]
     weights: [0.171118343, 0.190789926, 0.0943923462, 0.0940510903, 0.0906591629, 0.194029573, 0.16490176]
     test statistics:
       sigma: []
       densities: []
     references:
       []
     dof: 5


``MultivariateGaussian``
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observables``
     - sequence
     - yes
     - The qualified names of the observables this constraint applies to.
   * - ``kinematics``
     - sequence
     - yes
     - One kinematics map per observable, in the same order as 'observables'.
   * - ``options``
     - sequence
     - yes
     - One options map per observable, in the same order as 'observables'.
   * - ``means``
     - sequence
     - yes
     - The measurements' central values, in the same order as 'observables'.
   * - ``sigma-stat-hi``
     - sequence
     - yes
     - The measurements' upper statistical one-sigma uncertainties.
   * - ``sigma-stat-lo``
     - sequence
     - yes
     - The measurements' lower statistical one-sigma uncertainties.
   * - ``sigma-sys``
     - sequence
     - yes
     - The measurements' (symmetric) systematic one-sigma uncertainties.
   * - ``correlations``
     - sequence
     - yes
     - The measurements' correlation matrix, as a sequence of rows.
   * - ``dof``
     - scalar
     - no
     - The number of degrees of freedom; defaults to the number of measurements.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   B^0->K^*0gamma::S_K+C_K@BaBar:2008A:
     type: MultivariateGaussian
     observables:
       - B->K^*gamma::S_K^*gamma
       - B->K^*gamma::C_K^*gamma
     kinematics:
       - {}
       - {}
     options:
       - {q: d, tag: BFS2004}
       - {q: d, tag: BFS2004}
     means: [-0.03, -0.14]
     sigma-stat-hi: [0.29, 0.16]
     sigma-stat-lo: [0.29, 0.16]
     sigma-sys: [0.03, 0.03]
     correlations:
       - [1, 0.05]
       - [0.05, 1]
     references:
       []
     dof: 2


``MultivariateGaussian(Covariance)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observables``
     - sequence
     - yes
     - The qualified names of the observables (predictions) this constraint applies to.
   * - ``kinematics``
     - sequence
     - yes
     - One kinematics map per observable, in the same order as 'observables'.
   * - ``options``
     - sequence
     - yes
     - One options map per observable, in the same order as 'observables'.
   * - ``means``
     - sequence
     - yes
     - The measurements' central values.
   * - ``covariance``
     - sequence
     - yes
     - The measurements' covariance matrix, as a sequence of rows.
   * - ``response``
     - sequence
     - no
     - The response matrix mapping predictions to measurements, if they are not in 1:1 correspondence.
   * - ``dof``
     - scalar
     - no
     - The number of degrees of freedom; defaults to the number of measurements.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   B->pi::f_0+f_+@FLAG:2024A:
     type: MultivariateGaussian(Covariance)
     observables:
       - B->pi::f_+(q2)
       - B->pi::f_+(q2)
       - B->pi::f_+(q2)
       - B->pi::f_0(q2)
       - B->pi::f_0(q2)
     kinematics:
       - {q2: 18}
       - {q2: 22}
       - {q2: 26}
       - {q2: 18}
       - {q2: 26}
     options:
       - {form-factors: G2026}
       - {form-factors: G2026}
       - {form-factors: G2026}
       - {form-factors: G2026}
       - {form-factors: G2026}
     means: [1.07854, 2.0037, 6.11813, 0.48709, 0.95403]
     covariance:
       - [0.0034726032, 0.00499931193, 0.0114332044, -6.23611537e-06, -1.11172232e-05]
       - [0.00499931193, 0.00911756815, 0.0263339334, -1.1469038e-05, -2.1550285e-05]
       - [0.0114332044, 0.0263339334, 0.118033742, -1.94665745e-05, -3.47890026e-05]
       - [-6.23611537e-06, -1.1469038e-05, -1.94665745e-05, 0.000460187487, 0.000709222282]
       - [-1.11172232e-05, -2.1550285e-05, -3.47890026e-05, 0.000709222282, 0.00146191433]
     references:
       []
     dof: 5


``UniformBound``
~~~~~~~~~~~~~~~~

.. list-table::
   :widths: auto
   :header-rows: 1

   * - key
     - kind
     - required
     - description
   * - ``observables``
     - sequence
     - yes
     - The qualified names of the observables this constraint applies to.
   * - ``kinematics``
     - sequence
     - yes
     - One kinematics map per observable, in the same order as 'observables'.
   * - ``options``
     - sequence
     - yes
     - One options map per observable, in the same order as 'observables'.
   * - ``bound``
     - scalar
     - yes
     - The shared upper bound on all observables.
   * - ``uncertainty``
     - scalar
     - yes
     - The (relative) uncertainty on the bound.
   * - ``references``
     - sequence
     - no
     - The list of reference names this constraint originates from.

Example:

.. code-block:: yaml

   b->u::MesonicBound[1^+_A]:
     type: UniformBound
     observables:
       - B->pipi::Saturation[1^+_A]
     kinematics:
       - {}
     options:
       - {C: +-}
     bound: 1
     uncertainty: 0.05
     references:
       []
