===========================================
Hadronic Matrix Elements of Local Operators
===========================================

EOS uses the following conventions for the hadronic matrix elements of local operators,
i.e. the decay constants and the semileptonic form factors.

Definition of Decay Constants
-----------------------------

EOS uses meson to vacuum elements for the decay of both pseudoscalar (:math:`P`) and vector (:math:`V`) mesons.
The decay constants are defined as:

.. math::
   :nowrap:

   \begin{align*}
    \braket{0 | \bar{q}_1 \gamma^\mu \gamma_5 q_2 | P(p)} & = i f_P p^\mu, \\
    \braket{0 | \bar{q}_1 \gamma^\mu q_2 | V(p, \epsilon)} & = f_V m_V \epsilon^\mu, &
    \braket{0 | \bar{q}_1 \sigma^{\mu\nu} q_2 | V(p, \epsilon)} & = i f_V^T \left( \epsilon^\mu p^\nu - \epsilon^\nu p^\mu \right),
   \end{align*}

Here :math:`q_1` and :math:`q_2` are the quark fields, :math:`p` is the momentum of the initial-state meson,
and :math:`\epsilon` is the polarization vector of the vector meson.

The parameters describing these decay constants use the prefix ``decay-constant``.

Definition of Semileptonic Form Factors
---------------------------------------

EOS uses the following definitions for semileptonic form factors.

In the case of the decay of a pseudoscalar meson (:math:`P_1`) to another pseudoscalar meson (:math:`P_2`), the form factors are defined as:

.. math::
   :nowrap:

   \begin{align*}
      \braket{P_2(k) | \bar{q}_1 \gamma^\mu q_2 | P_1(p)}
         & = f_+^{P_1 \to P_2}(q^2) \left[ \left(p + k\right)^\mu - \frac{M_{P_1}^2 - M_{P_2}^2}{q^2} q^\mu\right]
           + f_0^{P_1 \to P_2}(q^2) \frac{M_{P_1}^2 - M_{P_2}^2}{q^2} q^\mu\,, \\
      \braket{P_2(k) | \bar{q}_1 \sigma^{\mu\nu} q_2 | P_1(p)}
         & = \frac{i f_T^{P_1 \to P_2}(q^2)}{M_{P_1} + M_{P_2}} \left[ (p + k)^\mu q^\nu - q^\mu (p + k)^\nu\right]\,.
   \end{align*}

Here, :math:`q = p - k` is the momentum transfer.

In the case of the decay of a pseudoscalar meson (:math:`P`) to a vector meson (:math:`V`), the form factors are defined as:

.. math::
   :nowrap:

   \begin{align*}
      \braket{V(k, \eta) | \bar{q}_1 \gamma^\mu q_2 | P(p)}
         & = \frac{2 V^{P \to V}(q^2)}{M_P + M_V} \varepsilon^{\mu\nu\alpha\beta} \eta^*_\nu p_\alpha k_\beta\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \gamma^\mu \gamma_5 q_2 | P(p)}
         & = i \eta_\nu^* \left[ A_1^{P \to V}(q^2) (M_P + M_V) g^{\mu\nu} - A_2^{P \to V}(q^2) \frac{(p + k)^\mu q_\nu}{M_P + M_V} - (A_3 - A_0) \frac{2 M_V q^\mu q^\nu}{q^2}\right]\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \sigma^{\mu\nu} q_2 | P(p)}
         & = 2 T_1^{P \to V} \varepsilon^{\mu\nu\alpha\beta} \eta^*_\nu p_\alpha k_\beta\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \sigma^{\mu\nu} \gamma_5 q_2 | P(p)}
         & = i T_2^{P \to V} \left[ \eta^{\mu*} (p + k)^\nu - \eta^{\nu*} (p + k)^\mu \right] + \frac{i T_3^{P \to V}}{M_P^2 - M_V^2} (\eta^* \cdot q) \left[q^\mu (p+k)^\nu - q^\nu (p+k)^\mu\right]\,,
   \end{align*}

where we abbreviate:

.. math::
   :nowrap:

   \begin{equation*}
      A_3^{P \to V}(q^2) = A_1^{P \to V}(q^2) \frac{M_P + M_V}{2 M_V} - A_2^{P \to V}(q^2) \frac{M_P - M_V}{2 M_V}\,.
   \end{equation*}

Here, :math:`\eta` is the polarization vector of the vector meson, and again :math:`q = p - k` is the momentum transfer.

In the case of the decay of a spin-parity :math:`J^P = 1/2^+` baryon (:math:`\mathcal{B}_1`) to another :math:`1/2^+` baryon (:math:`\mathcal{B}_2`),
EOS uses the following helicity-based definitions.
Abbreviating :math:`s_\pm \equiv (M_{\mathcal{B}_1} \pm M_{\mathcal{B}_2})^2 - q^2`,
the vector and axialvector matrix elements are defined as:

.. math::
   :nowrap:

   \begin{align*}
      \braket{\mathcal{B}_2(k, s') | \bar{q}_1 \gamma^\mu q_2 | \mathcal{B}_1(p, s)}
         & = \bar{u}(k, s') \bigg[
              f_t^V(q^2)\, (M_{\mathcal{B}_1} - M_{\mathcal{B}_2}) \frac{q^\mu}{q^2}
            + f_0^V(q^2)\, \frac{M_{\mathcal{B}_1} + M_{\mathcal{B}_2}}{s_+} \left( (p + k)^\mu - \frac{M_{\mathcal{B}_1}^2 - M_{\mathcal{B}_2}^2}{q^2} q^\mu \right) \\
         & \qquad\qquad + f_\perp^V(q^2) \left( \gamma^\mu - \frac{2 M_{\mathcal{B}_2}}{s_+} p^\mu - \frac{2 M_{\mathcal{B}_1}}{s_+} k^\mu \right)
              \bigg] u(p, s)\,, \\
      \braket{\mathcal{B}_2(k, s') | \bar{q}_1 \gamma^\mu \gamma_5 q_2 | \mathcal{B}_1(p, s)}
         & = -\bar{u}(k, s')\, \gamma_5 \bigg[
              f_t^A(q^2)\, (M_{\mathcal{B}_1} + M_{\mathcal{B}_2}) \frac{q^\mu}{q^2}
            + f_0^A(q^2)\, \frac{M_{\mathcal{B}_1} - M_{\mathcal{B}_2}}{s_-} \left( (p + k)^\mu - \frac{M_{\mathcal{B}_1}^2 - M_{\mathcal{B}_2}^2}{q^2} q^\mu \right) \\
         & \qquad\qquad + f_\perp^A(q^2) \left( \gamma^\mu + \frac{2 M_{\mathcal{B}_2}}{s_-} p^\mu - \frac{2 M_{\mathcal{B}_1}}{s_-} k^\mu \right)
              \bigg] u(p, s)\,,
   \end{align*}

while the tensor and axialtensor matrix elements are defined as:

.. math::
   :nowrap:

   \begin{align*}
      \braket{\mathcal{B}_2(k, s') | \bar{q}_1\, i \sigma^{\mu\nu} q_\nu q_2 | \mathcal{B}_1(p, s)}
         & = -\bar{u}(k, s') \bigg[
              f_0^T(q^2)\, \frac{q^2}{s_+} \left( (p + k)^\mu - \frac{M_{\mathcal{B}_1}^2 - M_{\mathcal{B}_2}^2}{q^2} q^\mu \right) \\
         & \qquad\qquad + f_\perp^T(q^2)\, (M_{\mathcal{B}_1} + M_{\mathcal{B}_2}) \left( \gamma^\mu - \frac{2 M_{\mathcal{B}_2}}{s_+} p^\mu - \frac{2 M_{\mathcal{B}_1}}{s_+} k^\mu \right)
              \bigg] u(p, s)\,, \\
      \braket{\mathcal{B}_2(k, s') | \bar{q}_1\, i \sigma^{\mu\nu} q_\nu \gamma_5 q_2 | \mathcal{B}_1(p, s)}
         & = -\bar{u}(k, s')\, \gamma_5 \bigg[
              f_0^{T5}(q^2)\, \frac{q^2}{s_-} \left( (p + k)^\mu - \frac{M_{\mathcal{B}_1}^2 - M_{\mathcal{B}_2}^2}{q^2} q^\mu \right) \\
         & \qquad\qquad + f_\perp^{T5}(q^2)\, (M_{\mathcal{B}_1} - M_{\mathcal{B}_2}) \left( \gamma^\mu + \frac{2 M_{\mathcal{B}_2}}{s_-} p^\mu - \frac{2 M_{\mathcal{B}_1}}{s_-} k^\mu \right)
              \bigg] u(p, s)\,.
   \end{align*}

Here, :math:`u(p, s)` and :math:`u(k, s')` are the Dirac spinors of the initial- and final-state baryon,
and again :math:`q = p - k` is the momentum transfer.
The ten form factors are labelled by the current (:math:`V`, :math:`A`, :math:`T`, :math:`T5`) and by their helicity:
:math:`t` (timelike), :math:`0` (longitudinal), and :math:`\perp` (transverse); the tensor currents have no timelike component.
Within EOS the corresponding parameters carry the labels ``time``, ``long``, and ``perp``.


Parametrisation of Semileptonic Form Factors
--------------------------------------------

By default, EOS uses the general-purpose ``SSE`` (Simplified Series Expansion) parametrisation for all semileptonic form factors, which is described in the subsection below.
The remaining parametrisations require domain-specific knowledge and are discussed in their respective subsections.

General Parametrisation: ``SSE``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EOS provides one general-purpose parametrisation of all semileptonic form factors. It is referred to as the ``SSE`` parametrisation.
For a generic form factor :math:`F(q^2)`, it reads:

.. math::
   :nowrap:

   \begin{equation*}
      F(q^2) = \frac{1}{1 - q^2 / M_R^2} \left[\sum_{i=0}^N \alpha^{(F)}_{i}\, z(q^2)^i \right]\,,
   \end{equation*}

with :math:`N = 2` in the current implementation.
Here :math:`M_R` corresponds to the mass of the first resonance seen by that form factor and :math:`\alpha_i^{(F)}` are free parameters.

This parametrisation uses the conformal mapping :math:`q^2 \mapsto z(q^2) = z(q^2; t_+, t_0)`, which reads

.. math::
   :nowrap:

   \begin{equation*}
      z(q^2; t_+, t_0) = \frac{\sqrt{t_+ - q^2} - \sqrt{t_+ - t_0}}{\sqrt{t_+ - q^2} + \sqrt{t_+ - t_0}}\,.
   \end{equation*}

In the above, :math:`t_+` is the pair-production threshold of the respective form factor, chosen as the *lowest-lying* two- or three-particle
threshold accessible in the isospin-symmetry limit. The expansion point :math:`t_0` is fixed to its optimised value
:math:`t_0 = t_+ \left(1 - \sqrt{1 - t_- / t_+}\right)`, where :math:`t_- \equiv (M_1 - M_2)^2` for hadron masses :math:`M_{1,2}`.

Since the vector and scalar currents couple to a lower-lying threshold than the axial and pseudoscalar currents, EOS distinguishes two thresholds:
:math:`t_+^V` for the vector/scalar-current form factors (e.g. :math:`V`, :math:`T_1`, and all :math:`P\to P` form factors),
and :math:`t_+^A` for the (pseudo)scalar/axial-current form factors (e.g. :math:`A_0`, :math:`A_1`, :math:`A_{12}`, :math:`T_2`, :math:`T_{23}`).
For a :math:`b\to s` transition, for instance, :math:`t_+^V = (M_B + M_K)^2` is the :math:`B K` threshold,
while :math:`t_+^A = (M_{B_s} + 2 M_\pi)^2` is the lowest isospin-allowed three-particle threshold.
Both thresholds are themselves accessible as parameters and are set to sensible default values.

The coefficients :math:`\alpha_i^{(F)}` are treated as unconstrained free parameters:
there is no expectation for the magnitude of any individual coefficient, nor is there an expectation that the series of coefficients converges.
In contrast to some other parametrisations, the series is expanded directly in :math:`z(q^2)`, and *not* in the shifted variable :math:`z(q^2) - z(0)`.

As a consequence, the kinematic constraints and equations of motion relating individual form factors
(for example :math:`f_0(0) = f_+(0)` for :math:`P\to P` transitions, or :math:`T_2(0) = T_1(0)` and :math:`A_{12}(0) = R\, A_0(0)` for :math:`P\to V` transitions,
with :math:`R \equiv (M_P^2 - M_V^2) / (8 M_P M_V)`)
fix the leading coefficient :math:`\alpha_0^{(F)}` of the constrained form factors. This coefficient is therefore removed from the set of free parameters.

Further relations hold at the kinematic endpoint :math:`t_- = (M_1 - M_2)^2`, where they ensure that the form factors obtained by dividing
by the Källén function :math:`\lambda(M_1^2, M_2^2, q^2)`, which vanishes at :math:`q^2 = t_-`, remain finite. For :math:`P\to V` transitions these are
:math:`A_{12}(t_-) = R\, A_1(t_-)`, which renders :math:`A_2` finite and fixes :math:`\alpha_0^{(A_1)}`, and
:math:`T_{23}(t_-) = K\, T_2(t_-)` with :math:`K \equiv (M_P + M_V)^2 / (4 M_P M_V)`, which renders :math:`T_3` finite and fixes :math:`\alpha_0^{(T_{23})}`.
For :math:`1/2^+ \to 1/2^+` transitions, the analogous relations are :math:`f_\perp^A(t_-) = f_\text{long}^A(t_-)` and
:math:`f_\text{long}^{T5}(t_-) = f_\perp^{T5}(t_-)`.

The ``SSE`` parametrisation covers :math:`P\to P` and :math:`P\to V` mesonic transitions as well as :math:`1/2^+ \to 1/2^+` baryonic transitions.
Examples for accessing the parameters :math:`\alpha_i^{(F)}` and the thresholds are:

+-----------------------------+-------------------+---------------------------------------------+
| Transition                  | Quantity          | Parameter(s)                                |
+=============================+===================+=============================================+
| :math:`B\to K`              | :math:`f_+`       | ``B->K::alpha^f+_0@SSE``,                   |
|                             |                   | ``B->K::alpha^f+_1@SSE``,                   |
|                             |                   | ...                                         |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`B\to K`              | :math:`t_+`       | ``B->K::tp@SSE``                            |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`B\to K`              | :math:`f_0`       | ``B->K::alpha^f0_1@SSE``,                   |
|                             |                   | ``B->K::alpha^f0_2@SSE``,                   |
|                             |                   | ...                                         |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`B\to K^*`            | :math:`V`         | ``B->K^*::alpha^V_0@SSE``,                  |
|                             |                   | ``B->K^*::alpha^V_1@SSE``,                  |
|                             |                   | ...                                         |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`B\to K^*`            | :math:`t_+^V`     | ``B->K^*::tp_v@SSE``                        |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`B\to K^*`            | :math:`A_1`       | ``B->K^*::alpha^A1_1@SSE``,                 |
|                             |                   | ``B->K^*::alpha^A1_2@SSE``,                 |
|                             |                   | ...                                         |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`B\to K^*`            | :math:`t_+^A`     | ``B->K^*::tp_a@SSE``                        |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`\Lambda_b\to\Lambda` | :math:`f_\perp^V` | ``Lambda_b->Lambda::alpha^(perp,V)_0@SSE``, |
|                             |                   | ``Lambda_b->Lambda::alpha^(perp,V)_1@SSE``, |
|                             |                   | ...                                         |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`\Lambda_b\to\Lambda` | :math:`t_+^V`     | ``Lambda_b->Lambda::tp_v@SSE``              |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`\Lambda_b\to\Lambda` | :math:`f_0^A`     | ``Lambda_b->Lambda::alpha^(0,A)_0@SSE``,    |
|                             |                   | ``Lambda_b->Lambda::alpha^(0,A)_1@SSE``,    |
|                             |                   | ...                                         |
+-----------------------------+-------------------+---------------------------------------------+
| :math:`\Lambda_b\to\Lambda` | :math:`t_+^A`     | ``Lambda_b->Lambda::tp_a@SSE``              |
+-----------------------------+-------------------+---------------------------------------------+


Dispersively Bounded Parametrisation: ``SE``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EOS further provides a ``Series Expansion`` (``SE``) parametrisation, in which the expansion coefficients are subject to a dispersive bound.
Using this parametrisation requires domain-specific knowledge and is intended for dedicated form-factor analyses.

The underlying idea goes back to [BGL:1997A]: the analyticity and unitarity of a suitable two-point correlator of the relevant quark current
imply an upper bound on a weighted integral of the form factors, and thereby on a sum of squared expansion coefficients.
A generic form factor is written as

.. math::
   :nowrap:

   \begin{equation*}
      F(q^2) = \frac{1}{P(q^2)\, \phi_F(q^2)} \sum_{k=0}^{K} a^{(F)}_{k}\, p_k\big(z(q^2)\big)\,,
   \end{equation*}

where :math:`P(q^2)` is a Blaschke product that removes all sub-threshold resonance poles, :math:`\phi_F` is a process- and form-factor-specific *outer function*,
and the :math:`p_k` are polynomials in the conformal variable :math:`z(q^2)`.
The outer functions are constructed such that the dispersive bound assumes the diagonal form

.. math::
   :nowrap:

   \begin{equation*}
      \sum_{F \in J^P} \sum_{k} \left(a^{(F)}_{k}\right)^2 \leq 1\,.
   \end{equation*}

This bound holds channel by channel: for each set of quantum numbers :math:`J^P` of the interpolating current, the first sum runs only over those form factors
:math:`F` that contribute to the corresponding two-point correlator. EOS therefore provides one saturation observable per channel and current, for example
``B->K::Saturation[1^-_V]`` and ``B->K::Saturation[0^+_V]``, rather than a single number for the entire transition.

The left-hand side is referred to as the *saturation* of the bound. It is obtained by integrating the squared, outer-function-weighted form factors
over the arc of the unit circle in the complex :math:`z` plane onto which the pair-production branch cut is mapped.
Choosing the :math:`p_k` to be orthonormal with respect to (the measure supported on) this arc is precisely what renders the bound diagonal.

Note that EOS absorbs the isospin-degeneracy factor :math:`n_I` of the transition into the normalisation of the outer functions, :math:`\phi_F \propto \sqrt{n_I}`.
Here :math:`n_I` counts the isospin-related exclusive channels that contribute to the same correlator; for example :math:`n_I = 2` for :math:`B\to K` (accounting for
both :math:`B^-\to K^-` and :math:`\bar{B}^0 \to \bar{K}^0`) and for :math:`B\to K^*`, while :math:`n_I = 1` for :math:`B_s\to \phi`, :math:`B_s\to K`, and all
baryonic transitions presently implemented.
As a consequence, the right-hand side of the bound above equals :math:`1` rather than :math:`1/n_I`, and the coefficients :math:`a^{(F)}_k` are larger by a factor
of :math:`\sqrt{n_I}` than in conventions that keep :math:`n_I` explicit on the right-hand side.

The concrete implementation in the traditional basis of mesonic form factors used within EOS follows [BFW:2010A].
For the baryonic :math:`1/2^+ \to 1/2^+` transitions, the dispersive bound, the integration over the arc of the unit circle, and the associated orthonormal polynomials :math:`p_k` follow [BMRvD:2022A].
The extension to :math:`1/2^+ \to 3/2^-` baryonic transitions follows [ABR:2022A].
Within EOS the corresponding coefficients carry the label ``SE``, for example ``B->K^*::a^V_0@SE`` or ``Lambda_b->Lambda::a^(0,V)_0@SE``.
