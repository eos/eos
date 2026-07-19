===================================
Observables in Semileptonic Decays
===================================

EOS predicts a large number of observables for exclusive semileptonic decays,
i.e. decays that proceed through a charged-current transition of the type described by the
`semileptonic charged-current operators <conventions.html#semileptonic-charged-current-operators>`_.
The available observables depend on the spin of the initial- and final-state hadrons.
This section summarises the observables that EOS provides for three common cases:
the decay of a pseudoscalar meson to a pseudoscalar meson (:math:`P`),
the decay of a pseudoscalar meson to a vector meson (:math:`V`), and
the decay of a spin-:math:`1/2` baryon to a spin-:math:`1/2` baryon.
All three cases are built from the corresponding
`semileptonic form factors <hadronic-matrix-elements.html#definition-of-semileptonic-form-factors>`_.

Throughout, :math:`q^2` denotes the invariant mass squared of the lepton-neutrino pair,
and the physical :math:`q^2` range extends from :math:`m_\ell^2` to :math:`(M_1 - M_2)^2`,
where :math:`M_1` and :math:`M_2` are the masses of the initial- and final-state hadron.
For each decay, EOS provides observables that are

* *fully differential* in :math:`q^2` and the decay angles;
* *single differential* in :math:`q^2`, obtained by integrating over the angles; and
* *integrated* over a :math:`q^2` interval :math:`[q^2_\text{min}, q^2_\text{max}]`.

In addition, EOS provides a set of *normalised* observables, in which the dependence on the CKM matrix
element :math:`|V_{UD}|` is removed (i.e. :math:`|V_{UD}| = 1`), and a set of *lepton-flavour-universality (LFU) ratios*.

Decays to a Pseudoscalar Meson
------------------------------

For a decay :math:`P_1 \to P_2\, \ell\, \bar\nu` to a pseudoscalar meson, the kinematics are
described by :math:`q^2` and a single angle: the helicity angle :math:`\theta_\ell` of the charged lepton,
defined in the rest frame of the lepton-neutrino pair.
The fully differential decay rate is a second-order polynomial in :math:`\cos\theta_\ell`:

.. math::
   :nowrap:

   \begin{equation*}
      \frac{d^2\Gamma}{dq^2\, d\cos\theta_\ell}
        = a_\ell(q^2) + b_\ell(q^2) \cos\theta_\ell + c_\ell(q^2) \cos^2\theta_\ell\,.
   \end{equation*}

From this distribution EOS derives, as functions of :math:`q^2` and as quantities integrated over a :math:`q^2` bin:

* the differential decay width :math:`d\Gamma/dq^2` and branching ratio :math:`d\mathcal{B}/dq^2`
  (e.g. ``B->pilnu::dBR/dq2``), and the integrated branching ratio (``B->pilnu::BR``);
* the leptonic forward-backward asymmetry :math:`A_\text{FB}(q^2)`
  (``B->pilnu::A_FB(q2)``), which is sensitive to :math:`b_\ell(q^2)`;
* the flat term :math:`F_H` (``B->pilnu::F_H``), which quantifies the angular-independent contribution
  to the distribution and vanishes in the massless-lepton limit; and
* the longitudinal lepton-polarisation asymmetry :math:`A_\lambda^\ell` (``B->pilnu::A_l``).

LFU ratios compare the rates into different lepton flavours, for example

.. math::
   :nowrap:

   \begin{equation*}
      R_{P} = \frac{\mathcal{B}(P_1 \to P_2\, \tau\, \bar\nu)}{\mathcal{B}(P_1 \to P_2\, \ell\, \bar\nu)}\,,
      \qquad \ell = e, \mu\,,
   \end{equation*}

which EOS exposes both integrated (``B->pilnu::R_pi``) and differential in :math:`q^2` (``B->pilnu::R_pi(q2)``).

Decays to a Vector Meson
------------------------

For a decay :math:`P \to V(\to P_1 P_2)\, \ell\, \bar\nu` to a vector meson, the subsequent decay of the
vector meson makes the full angular distribution accessible.
The kinematics are described by :math:`q^2` and three angles: the helicity angle :math:`\theta_\ell` of the
charged lepton, the helicity angle :math:`\theta_V` of the vector meson (defined through its decay products),
and the azimuthal angle :math:`\phi` between the two decay planes.
The fully differential decay rate reads

.. math::
   :nowrap:

   \begin{align*}
      \frac{d^4\Gamma}{dq^2\, d\cos\theta_\ell\, d\cos\theta_V\, d\phi}
        = \frac{9}{32\pi} \Big[
        & \left(J_{1s} \sin^2\theta_V + J_{1c} \cos^2\theta_V\right)
          + \left(J_{2s} \sin^2\theta_V + J_{2c} \cos^2\theta_V\right) \cos 2\theta_\ell \\
        & + J_3 \sin^2\theta_V \sin^2\theta_\ell \cos 2\phi
          + J_4 \sin 2\theta_V \sin 2\theta_\ell \cos\phi
          + J_5 \sin 2\theta_V \sin\theta_\ell \cos\phi \\
        & + \left(J_{6s} \sin^2\theta_V + J_{6c} \cos^2\theta_V\right) \cos\theta_\ell
          + J_7 \sin 2\theta_V \sin\theta_\ell \sin\phi \\
        & + J_8 \sin 2\theta_V \sin 2\theta_\ell \sin\phi
          + J_9 \sin^2\theta_V \sin^2\theta_\ell \sin 2\phi
        \Big]\,,
   \end{align*}

where the twelve angular coefficient functions :math:`J_i \equiv J_i(q^2)` encode the full information
of the decay.
EOS provides each of the :math:`J_i` both differentially in :math:`q^2` (e.g. ``B->D^*lnu::J_1c(q2)``)
and integrated over a :math:`q^2` bin (e.g. ``B->D^*lnu::J_1c``).

From the angular coefficients EOS derives further observables, as functions of :math:`q^2` and integrated
over a :math:`q^2` bin, including:

* the differential decay width :math:`d\Gamma/dq^2` and branching ratio :math:`d\mathcal{B}/dq^2`
  (e.g. ``B->D^*lnu::dBR/dq2``), and the integrated branching ratio (``B->D^*lnu::BR``);
* the leptonic forward-backward asymmetry :math:`A_\text{FB}(q^2)` (``B->D^*lnu::A_FB(q2)``);
* the longitudinal polarisation fraction :math:`F_\text{L}` of the vector meson (``B->D^*lnu::F_L``); and
* a set of angular and CP asymmetries :math:`A_\text{C}^{1,2,3}` and :math:`A_\text{T}^{1,2,3}`
  (e.g. ``B->D^*lnu::A_C^1``, ``B->D^*lnu::A_T^1``).

As in the pseudoscalar case, EOS provides LFU ratios, for example

.. math::
   :nowrap:

   \begin{equation*}
      R_{V} = \frac{\mathcal{B}(P \to V\, \tau\, \bar\nu)}{\mathcal{B}(P \to V\, \ell\, \bar\nu)}\,,
      \qquad \ell = e, \mu\,,
   \end{equation*}

available both integrated (``B->D^*lnu::R_D^*``) and differential in :math:`q^2` (``B->D^*lnu::R_{D^*}^{tau/mu}(q2)``).

Decays to a Spin-1/2 Baryon
---------------------------

For the decay of a spin-parity :math:`J^P = 1/2^+` baryon to another :math:`1/2^+` baryon,
:math:`\mathcal{B}_1 \to \mathcal{B}_2\, \ell\, \bar\nu`
(e.g. :math:`\Lambda_b \to \Lambda_c\, \ell\, \bar\nu`), the final-state baryon is itself unstable.
Its subsequent decay :math:`\mathcal{B}_2 \to \mathcal{B}_3\, \pi` (e.g. :math:`\Lambda_c \to \Lambda\, \pi`)
acts as a polarisation analyser and thereby makes the full angular distribution accessible.
The kinematics are described by :math:`q^2` and three angles: the helicity angle :math:`\theta_\ell`
of the charged lepton (in the lepton-neutrino rest frame), the helicity angle :math:`\theta_\mathcal{B}`
of the daughter baryon (defined through :math:`\mathcal{B}_2 \to \mathcal{B}_3\, \pi`), and the azimuthal
angle :math:`\phi` between the two decay planes.
Assuming the initial-state baryon :math:`\mathcal{B}_1` to be **unpolarised**, the fully differential decay rate reads

.. math::
   :nowrap:

   \begin{align*}
      \frac{d^4\Gamma}{dq^2\, d\cos\theta_\ell\, d\cos\theta_\mathcal{B}\, d\phi}
        = \frac{3}{8\pi} \Big[
        & \left(K_{1ss} \sin^2\theta_\ell + K_{1cc} \cos^2\theta_\ell + K_{1c} \cos\theta_\ell\right) \\
        & + \left(K_{2ss} \sin^2\theta_\ell + K_{2cc} \cos^2\theta_\ell + K_{2c} \cos\theta_\ell\right) \cos\theta_\mathcal{B} \\
        & + \left(K_{3sc} \sin\theta_\ell \cos\theta_\ell + K_{3s} \sin\theta_\ell\right) \sin\theta_\mathcal{B} \sin\phi \\
        & + \left(K_{4sc} \sin\theta_\ell \cos\theta_\ell + K_{4s} \sin\theta_\ell\right) \sin\theta_\mathcal{B} \cos\phi
        \Big]\,,
   \end{align*}

where the ten angular coefficient functions :math:`K_i \equiv K_i(q^2)` encode the full information
of the decay.
A non-zero polarisation of the initial-state baryon would introduce additional angular structures,
and hence further angular coefficients, beyond these ten.
EOS provides the :math:`K_i` integrated over a :math:`q^2` bin and normalised to the decay width
(e.g. ``Lambda_b->Lambda_clnu::K_1ss``).

From the angular coefficients EOS derives further observables, as functions of :math:`q^2` and integrated
over a :math:`q^2` bin, including:

* the differential decay width and branching ratio :math:`d\mathcal{B}/dq^2`
  (``Lambda_b->Lambda_clnu::dBR/dq2``), and the integrated branching ratio (``Lambda_b->Lambda_clnu::BR``);
* the leptonic forward-backward asymmetry :math:`A_\text{FB}^\ell` (``Lambda_b->Lambda_clnu::A_FB^l(q2)``);
* the hadronic forward-backward asymmetry :math:`A_\text{FB}^h` (``Lambda_b->Lambda_clnu::A_FB^h(q2)``),
  which arises from the polarisation of the daughter baryon;
* the combined lepton-hadron forward-backward asymmetry :math:`A_\text{FB}^{h\ell}`
  (``Lambda_b->Lambda_clnu::A_FB^c(q2)``); and
* the fraction :math:`F_0` (``Lambda_b->Lambda_clnu::F_0(q2)``).

As in the meson cases, EOS provides LFU ratios, for example

.. math::
   :nowrap:

   \begin{equation*}
      R_{\mathcal{B}_2} = \frac{\mathcal{B}(\mathcal{B}_1 \to \mathcal{B}_2\, \tau\, \bar\nu)}{\mathcal{B}(\mathcal{B}_1 \to \mathcal{B}_2\, \ell\, \bar\nu)}\,,
      \qquad \ell = e, \mu\,,
   \end{equation*}

available as an integrated observable (``Lambda_b->Lambda_clnu::R(Lambda_c)``).
