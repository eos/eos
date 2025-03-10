==================================
Summary of the Physics Conventions
==================================

EOS implicitly uses the following conventions in the predictions of flavour physics observables.

Units
-----

EOS uses natural units :math:`(\hbar = c =1)` and, with the exceptions of lifetimes that are expressed in seconds, all quantities are expressed in powers of GeV.
Units are specified in the definition of observables and parameters.
They can be found in this documentation (`here <constraints.html>`_ and `here <parameters.html>`_ respectively),
or directly from python, as described in the `basic examples <../user-guide/basics.html>`_.

Weak Effective Theory
---------------------

We follow the conventions of the WCxf format for the Wilson coefficients of the effective Hamiltonian at mass dimension six.
Each sector of the effective Hamiltonian is defined by a unique name, corresponding to the flavour quantum numbers
of the fermion fields.

Semileptonic Charged-Current Operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the Bern class-III notation for the semileptonic charged-current operators.
For down-type (:math:`D`) to up-type (:math:`U`) semileptonic decays, the sector of the effective theory is described by the Lagrangian:

.. math::
   :nowrap:

   \begin{align*}
      \mathcal{L}^{UD\ell\nu}
        = -\frac{4 G_\text{F}}{\sqrt{2}} V_{UD} \sum_{i} \mathcal{C}_i^{UD\ell\nu}(\mu) \, \mathcal{O}_i + \text{h.c.}\,,
   \end{align*}


where :math:`V_{UD}` is the CKM matrix element, :math:`G_\text{F}` is the Fermi constant. The operator basis :math:`\mathcal{O}_i` reads:

.. math::
   :nowrap:

   \begin{align*}
      \mathcal{O}_{VL}^{UD\ell\nu} & = [\bar{U} \gamma_\mu      P_L D] [\bar{\ell} \gamma^\mu      P_L \nu]\,, &
      \mathcal{O}_{VR}^{UD\ell\nu} & = [\bar{U} \gamma_\mu      P_R D] [\bar{\ell} \gamma^\mu      P_L \nu]\,, \\
      \mathcal{O}_{SL}^{UD\ell\nu} & = [\bar{U}                 P_L D] [\bar{\ell}                 P_L \nu]\,, &
      \mathcal{O}_{SR}^{UD\ell\nu} & = [\bar{U}                 P_R D] [\bar{\ell}                 P_L \nu]\,, \\
      \mathcal{O}_{T}^{UD\ell\nu}  & = [\bar{U} \sigma_{\mu\nu}     D] [\bar{\ell} \sigma^{\mu\nu} P_L \nu]\,.
   \end{align*}

The Wilson coefficients :math:`\mathcal{C}_i^{UD\ell\nu}(\mu)` are defined at the sector-specific scale :math:`\mu^{UD\ell\nu}` and are dimensionless.
Their parameters use the prefix ``UDlnul``, where ``U``, ``D``, and ``l`` are understood as metavariables.


Hadronic Matrix Elements
------------------------

Definition of Decay Constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
         & = i \eta_\nu^* \left[ A_1^{P \to V}(q^2) (M_B + M_V) g^{\mu\nu} - A_2^{P \to V}(q^2) \frac{(p + k)^\mu q_\nu}{M_B + M_V} - (A_3 - A_0) \frac{2 M_V q^\mu q^\nu}{q^2}\right]\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \sigma^{\mu\nu} q_2 | P(p)}
         & = 2 T_1^{P \to V} \varepsilon^{\mu\nu\alpha\beta} \eta^*_\nu p_\alpha k_\beta\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \sigma^{\mu\nu} \gamma_5 q_2 | P(p)}
         & = i \eta^*_\nu \left[ T_2^{P \to V} \left( (M_P^2 - M_V^2) g^{\mu\nu} - (p + k)^\mu q^\nu \right) + T_3^{P \to V} \left(q^\mu - \frac{q^2}{M_P^2 - M_V^2} (p + k)^\mu\right) q^\nu\right]\,,
   \end{align*}

where we abbreviate:

.. math::
   :nowrap:

   \begin{equation*}
      A_3^{P \to V}(q^2) = A_1^{P \to V}(q^2) \frac{M_B + M_V}{2 M_V} - A_2^{P \to V}(q^2) \frac{M_B - M_V}{2 M_V}\,.
   \end{equation*}

Here, :math:`\eta` is the polarization vector of the vector meson, and again :math:`q = p - k` is the momentum transfer.


Parametrisation of Semileptonic Form Factors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

EOS provides one general-purpose parametrisation of all semileptonic form factors. It is referred to as the ``BSZ2015`` parametrisation.
For a generic form factor :math:`F(q^2)`, it reads:

.. math::
   :nowrap:

   \begin{equation*}
      F(q^2) = \frac{1}{q^2 - M_F^2} \left[\sum_{i=0}^N \alpha^{(F)}_{i} z^i(q^2) \right]\,.
   \end{equation*}

Here :math:`M_R` correspond to the mass of the first resonance seen by that form factor and :math:`\alpha_i^{(F)}` are free parameters.
We use the conformal mapping :math:`q^2 \mapsto z(q^2) = z(q^2; t_+, t_0)`, which reads

.. math::
   :nowrap:

   \begin{equation*}
      z(q^2; t_+, t_0) = \frac{\sqrt{t_+ - q^2} - \sqrt{t_+ - t_0}}{\sqrt{t_+ - q^2} - \sqrt{t_+ - t_0}}\,.
   \end{equation*}

In the above, :math:`t_+ \equiv (M_1 + M_2)^2` represents the two-body threshold for the respective form factor for a reaction with hadron masses :math:`M_{1,2}`,
and :math:`t_0 < (M_1 - M_2)^2` is a free parameter.
Examples for accessing the parameters :math:`\alpha_i^{(F)}` are:

+-------------------+--------------+------------------------------+
| Transition        | Form Factor  | Parameter                    |
+-------------------+--------------+------------------------------+
| :math:`B\to K`    | :math:`f_+`  | ``B->K::alpha^f+_0@BSZ2015`` |
+-------------------+--------------+------------------------------+
| :math:`D\to K`    | :math:`f_0`  | ``D->K::alpha^f0_2@BSZ2015`` |
+-------------------+--------------+------------------------------+
| :math:`B\to K^*`  | :math:`V`    | ``B->K^*::alpha^V_1@BZS2015``|
+-------------------+--------------+------------------------------+
