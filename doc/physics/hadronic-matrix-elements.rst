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
         & = i \eta_\nu^* \left[ A_1^{P \to V}(q^2) (M_B + M_V) g^{\mu\nu} - A_2^{P \to V}(q^2) \frac{(p + k)^\mu q_\nu}{M_B + M_V} - (A_3 - A_0) \frac{2 M_V q^\mu q^\nu}{q^2}\right]\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \sigma^{\mu\nu} q_2 | P(p)}
         & = 2 T_1^{P \to V} \varepsilon^{\mu\nu\alpha\beta} \eta^*_\nu p_\alpha k_\beta\,, \\
      \braket{V(k, \eta) | \bar{q}_1 \sigma^{\mu\nu} \gamma_5 q_2 | P(p)}
         & = i T_2^{P \to V} \left[ \eta^{\mu*} (p + k)^\nu - \eta^{\nu*} (p + k)^\mu \right] + \frac{i T_3^{P \to V}}{M_P^2 - M_V^2} (\eta^* \cdot q) \left[q^\mu (p+k)^\nu - q^\nu (p+k)^\mu\right]\,,
   \end{align*}

where we abbreviate:

.. math::
   :nowrap:

   \begin{equation*}
      A_3^{P \to V}(q^2) = A_1^{P \to V}(q^2) \frac{M_B + M_V}{2 M_V} - A_2^{P \to V}(q^2) \frac{M_B - M_V}{2 M_V}\,.
   \end{equation*}

Here, :math:`\eta` is the polarization vector of the vector meson, and again :math:`q = p - k` is the momentum transfer.


Parametrisation of Semileptonic Form Factors
--------------------------------------------

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
