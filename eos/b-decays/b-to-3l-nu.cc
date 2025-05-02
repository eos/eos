/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2025 Danny van Dyk
 * Copyright (c) 2022      Stephan Kürten
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/b-decays/b-to-3l-nu.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::norm;
    using std::abs;
    using std::real;
    using std::imag;
    using std::sqrt;
    using std::conj;
    using std::log;
    using std::cos;
    using std::sin;

    /*
     * Decay: B_q^- -> lprime^+ lprime^- l^- nubar, cf. [KKvDZ:2022A]
     */
    template <>
    struct Implementation<BToThreeLeptonsNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToGammaOffShell>> form_factors;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        UsedParameter alpha_qed;

        LeptonFlavorOption opt_lprime;

        UsedParameter m_lprime;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            form_factors(FormFactorFactory<PToGammaOffShell>::create("B->gamma^*::" + o.get("form-factors"_ok, "KKvDZ2022"), p, o)),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_B(p["mass::B_u"], u),
            f_B(p["decay-constant::B_u"], u),
            tau_B(p["life_time::B_u"], u),
            alpha_qed(p["QED::alpha_e(m_b)"],u),
            opt_lprime(o, options, "lprime"_ok),
            m_lprime(p["mass::" + opt_lprime.str()], u),
            opt_l(o, options, "l"_ok),
            m_l(p["mass::" + opt_l.str()], u)
        {
            Context ctx("When constructing B->l'l'lnu observable");

            u.uses(*model);
            u.uses(*form_factors);

            if (opt_l.str() == opt_lprime.str())
            {
                throw InvalidOptionValueError("lprime"_ok, opt_l.str(), "e, mu, tau (but may not be the value of l)");
            }
        }

        /*!
         * differential decay width of 2 kinematic variables
         * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - m_l )^2
         * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B - sqrt(q2) )^2
         */
        double double_differential_decay_width(const double & q2, const double & k2) const
        {
            const double m_B_sq = m_B * m_B, m_B_3 = m_B_sq * m_B, m_B_4 = m_B_sq * m_B_sq;
            const double m_B_6 = m_B_3 * m_B_3, m_B_8 = m_B_4 * m_B_4;
            const double m_l_sq = m_l * m_l, m_l3 = m_l_sq * m_l, m_l_4 = m_l_sq * m_l_sq;
            const double m_l_6 = m_l3 * m_l3;
            const double m_lprime_sq = m_lprime * m_lprime;

            const complex<double> F_1 = form_factors->F_1(q2, k2);
            const complex<double> F_2 = form_factors->F_2(q2, k2);
            const complex<double> F_3 = form_factors->F_3(q2, k2);
            const complex<double> F_4 = form_factors->F_4(q2, k2);

            const double prefactor = (power_of<2>(g_fermi * abs(model->ckm_ub())
            * alpha_qed * 4 * M_PI) / (2.0 * power_of<2>(q2)) * (1.0 - m_l_sq / k2)
            * sqrt(1.0 - (4.0 * m_lprime_sq) / q2) * sqrt(power_of<2>(k2) - 2.0 * k2 * m_B_sq
            + m_B_4 - 2.0 * k2 * q2 - 2.0 * m_B_sq * q2 + power_of<2>(q2))
            / (32768.0 * m_B_3 * power_of<6>(M_PI)));

            const double f11 = (64 * (2 * power_of<2>(k2) - k2 * m_l_sq - m_l_4) * M_PI
            * (power_of<2>(k2) - 2 * k2 * (m_B_sq - 2 * q2) + power_of<2>(m_B_sq - q2))
            * (2 * m_lprime_sq + q2)) / (9. * k2 * m_B_sq);
            const double f22 = (32 * (2 * power_of<2>(k2) - k2 * m_l_sq - m_l_4) * M_PI
            * (power_of<2>(k2) - 2 * k2 * (m_B_sq - 5 * q2) + power_of<2>(m_B_sq - q2))
            * q2 * (2 * m_lprime_sq + q2)) / (9. * power_of<2>(k2) * m_B_sq);
            const double f33 = (32 * m_l_sq * (k2 - m_l_sq) * M_PI * q2 * (2 * m_lprime_sq + q2)
            * (power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2))) / (3.
            * power_of<2>(k2) * m_B_sq);
            const double f44 = (64 * (k2 - m_l_sq) * (2 * k2 + m_l_sq) * M_PI * (2 * m_lprime_sq
            + q2) * (power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2)))
            / (9. * k2 * m_B_sq);
            const double f12 = (128 * (2 * power_of<2>(k2) - k2 * m_l_sq - m_l_4) * M_PI * q2
            * (k2 - m_B_sq + q2) * (2 * m_lprime_sq + q2)) / (3. * k2 * m_B_sq);
            const double f15 = -(256 * m_l_sq * M_PI * (2 * m_lprime_sq + q2) * (- ((k2 - m_l_sq)
            * (2 * k2 - m_B_sq + m_l_sq + q2) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2)
            - 2 * k2 * (m_B_sq + q2))) + 2 * k2 * (k2 - m_B_sq) * (k2 - m_B_sq + 2 * m_l_sq + q2)
            * std::atanh(((k2 - m_l_sq) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2
            * (m_B_sq + q2)))/(power_of<2>(k2) + m_l_sq * (- m_B_sq + q2) - k2
            * (m_B_sq - m_l_sq + q2))))) / (3. * m_B * (k2 - m_l_sq) * sqrt(power_of<2>(k2)
            + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2)));
            const double f25 = -(- 256 * m_l_sq * M_PI * q2 * (2 * m_lprime_sq + q2)
            * (3 * (k2 - m_l_sq) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2
            * (m_B_sq + q2)) + (- 4 * power_of<2>(k2) + 6 * k2 * m_B_sq - 4 * k2 * m_l_sq
            + 2 * m_l_4) * std::atanh(((k2 - m_l_sq) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2)
            - 2 * k2 * (m_B_sq + q2)))/(power_of<2>(k2) + m_l_sq * (- m_B_sq + q2) - k2 * (m_B_sq
            - m_l_sq + q2))))) / (3. * m_B * (k2 - m_l_sq) * sqrt(power_of<2>(k2)
            + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2)));
            const double f35 = -(- 128 * m_l_sq * M_PI * q2 * (2 * m_lprime_sq + q2)
            * (- ((k2 - m_l_sq) * (power_of<2>(k2) + k2 * (3 * m_B_sq - 3 * m_l_sq - q2)
            + m_l_sq * (- m_B_sq + q2)) * sqrt(power_of<2>(k2)
            + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2))) + 4 * k2 * (k2 - m_B_sq)
            * (k2 * m_B_sq - m_l_4) * std::atanh(((k2 - m_l_sq) * sqrt(power_of<2>(k2)
            + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2)))/(power_of<2>(k2) + m_l_sq
            * (- m_B_sq + q2) - k2 * (m_B_sq - m_l_sq + q2))))) / (3. * k2 * m_B * (k2 - m_B_sq)
            * (k2 - m_l_sq) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2)
            - 2 * k2 * (m_B_sq + q2)));
            const double f45 = -(- 256 * m_l_sq * M_PI * (2 * m_lprime_sq + q2) * (- ((k2 - m_l_sq)
            * (k2 - m_B_sq + q2) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2
            * (m_B_sq + q2))) + 2 * k2 * (power_of<2>(k2) + m_B_4 - m_B_sq * q2 + 2 * m_l_sq
            * q2 - k2 * (2 * m_B_sq + q2)) * std::atanh(((k2 - m_l_sq) * sqrt(power_of<2>(k2)
            + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2)))/(power_of<2>(k2) + m_l_sq
            * (- m_B_sq + q2) - k2 * (m_B_sq - m_l_sq + q2))))) / (3. * m_B * (k2 - m_l_sq)
            * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2)));
            const double f55 = (128 * m_l_sq * M_PI * (2 * m_lprime_sq + q2) * (- ((k2 - m_l_sq)
            * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2))
            * (2 * power_of<4>(k2) * m_l_sq + power_of<3>(k2) * (- 8 * m_l_4 + 4 * m_B_sq * q2
            - 5 * m_l_sq * q2) + m_l_sq * (2 * m_B_8 - 2 * m_B_6 * q2
            + 3 * m_B_4 * m_l_sq * q2 - m_B_sq * m_l_sq * power_of<2>(q2) + m_l_4
            * power_of<2>(q2)) + k2 * (4 * m_B_6 * q2 + 2 * m_B_sq * m_l_sq * q2 * (m_l_sq + q2)
            - 2 * m_l_4 * q2 * (2 * m_l_sq + q2) - m_B_4 * (8 * m_l_4 + 5 * m_l_sq * q2))
            + power_of<2>(k2) * (- 4 * m_B_4 * (m_l_sq + q2) + m_l_sq * q2 * (7 * m_l_sq + q2)
            + m_B_sq * (16 * m_l_4 - power_of<2>(q2))))) + 4 * k2 * (k2 - m_B_sq)
            * (power_of<4>(k2) * m_l_sq + power_of<3>(k2) * (m_B_sq * (- 2 * m_l_sq + q2)
            - 2 * m_l_sq * (m_l_sq + q2)) + power_of<2>(k2) * (2 * m_B_4 * m_l_sq - 2 * m_l_6
            + 3 * m_l_4 * q2 + m_l_sq * power_of<2>(q2) + m_B_sq * (6 * m_l_4 - power_of<2>(q2)))
            + k2 * (- (m_l_4 * power_of<2>(q2)) + m_B_6 * (- 2 * m_l_sq + q2) + m_B_sq
            * (4 * m_l_6 - 2 * m_l_4 * q2) + m_B_4 * (- 6 * m_l_4 - 2 * m_l_sq * q2
            + power_of<2>(q2))) + m_l_sq * (m_B_8 + 2 * m_B_6 * m_l_sq - 2 * m_l_6 * q2 + m_B_sq
            * m_l_sq * q2 * (4 * m_l_sq + q2) - m_B_4 * (2 * m_l_4 + m_l_sq * q2
            + power_of<2>(q2)))) * std::atanh(((k2 - m_l_sq) * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq
            - q2) - 2 * k2 * (m_B_sq + q2)))/(power_of<2>(k2) + m_l_sq * (- m_B_sq + q2) - k2
            * (m_B_sq - m_l_sq + q2))))) / (3. * power_of<2>(k2 - m_B_sq) * (k2 - m_l_sq)
            * sqrt(power_of<2>(k2) + power_of<2>(m_B_sq - q2) - 2 * k2 * (m_B_sq + q2))
            * (power_of<2>(k2) * m_l_sq - k2 * m_l_sq * q2 + k2 * m_B_sq * (- 2 * m_l_sq + q2)
            + m_l_sq * (m_B_4 - m_B_sq * q2 + m_l_sq * q2)));

            const double amp = f11 * norm(F_1) + f22 * norm(F_2) + f33 * norm(F_3)
            + f44 * norm(F_4) + f12 * real(F_1 * conj(F_2)) + f15 * f_B * real(F_1)
            + f25 * f_B * real(F_2) + f35 * f_B * real(F_3) + f45 * f_B * real(F_4)
            + f55 * power_of<2>(f_B);

            if ( power_of<2>(m_B - m_l) < q2 || q2 < 4 * m_lprime_sq
                || power_of<2>(m_B - sqrt(q2)) < k2 || k2 < m_l_sq )
            {
                return 0.0;
            }

            return prefactor * amp;
        }

        /*!
         * differential decay width of 5 kinematic variables
         * quintuple_differential_decay_width(q2, k2, z_gamma, z_w, phi)
         * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - m_l )^2
         * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B - sqrt(q2) )^2
         * z_gamma is the angle between the negatively charged lepton l' and the negative z-axis
         * z_w is the angle between the charged lepton l and the positive z-axis
         * phi is the angle between the q2 plane and the k2 plane
         */
        double quintuple_differential_decay_width(const double & q2, const double & k2,
                const double & z_gamma, const double & z_w, const double & phi) const
        {
            const double m_B_sq = m_B * m_B, m_B_3 = m_B_sq * m_B, m_B_4 = m_B_sq * m_B_sq;
            const double m_B_6 = m_B_3 * m_B_3, m_B_8 = m_B_4 * m_B_4, m_l_sq = m_l * m_l;
            const double m_l_4 = m_l_sq * m_l_sq, m_lprime_sq = m_lprime * m_lprime;
            const double k4 = power_of<2>(k2), q4 = power_of<2>(q2), z_w_sq = power_of<2>(z_w);
            const double z_gamma_sq = power_of<2>(z_gamma);
            const double k6 = k4 * k2, q6 = q4 * q2, k8 = k4 * k4, q8 = q4 * q4, k10 = k8 * k2;
            const double m_B_sq_q2_diff2 = power_of<2>(m_B_sq - q2);
            const double m_l_sq_k2_diff2 = power_of<2>(m_l_sq - k2);
            const double m_B_sq_k2_diff2 = power_of<2>(m_B_sq - k2);

            const complex<double> F_1 = form_factors->F_1(q2, k2);
            const complex<double> F_2 = form_factors->F_2(q2, k2);
            const complex<double> F_3 = form_factors->F_3(q2, k2);
            const complex<double> F_4 = form_factors->F_4(q2, k2);

            const double prefactor = (power_of<2>(g_fermi * abs(model->ckm_ub())
            * alpha_qed * 4 * M_PI) / (2.0 * power_of<2>(q2)) * (1.0 - m_l_sq / k2)
            * sqrt(1.0 - (4.0 * m_lprime_sq) / q2) * sqrt(power_of<2>(k2) - 2.0 * k2 * m_B_sq
            + m_B_4 - 2.0 * k2 * q2 - 2.0 * m_B_sq * q2 + power_of<2>(q2))
            / (32768.0 * m_B_3 * power_of<6>(M_PI)));

            const double f11 = - (((k2 - m_l_sq) * (- 4 * k6 * m_lprime_sq + 8 * k4 * m_B_sq
            * m_lprime_sq - 4 * k2 * m_B_4 * m_lprime_sq - 4 * k4 * m_l_sq * m_lprime_sq + 8 * k2
            * m_B_sq * m_l_sq * m_lprime_sq - 4 * m_B_4 * m_l_sq * m_lprime_sq - k6 * q2 + 2 * k4
            * m_B_sq * q2 - k2 * m_B_4 * q2 - k4 * m_l_sq * q2 + 2 * k2 * m_B_sq * m_l_sq * q2
            - m_B_4 * m_l_sq * q2 - 8 * k4 * m_lprime_sq * q2 + 8 * k2 * m_B_sq * m_lprime_sq * q2
            - 8 * k2 * m_l_sq * m_lprime_sq * q2 + 8 * m_B_sq * m_l_sq * m_lprime_sq * q2 - 10 * k4
            * q4 + 2 * k2 * m_B_sq * q4 - 2 * k2 * m_l_sq * q4 + 2 * m_B_sq * m_l_sq * q4 - 4 * k2
            * m_lprime_sq * q4 - 4 * m_l_sq * m_lprime_sq * q4 - k2 * q6 - m_l_sq * q6 + 4 * k6
            * m_lprime_sq * z_gamma_sq - 8 * k4 * m_B_sq * m_lprime_sq * z_gamma_sq + 4 * k2 * m_B_4
            * m_lprime_sq * z_gamma_sq + 4 * k4 * m_l_sq * m_lprime_sq * z_gamma_sq - 8 * k2
            * m_B_sq * m_l_sq * m_lprime_sq * z_gamma_sq + 4 * m_B_4 * m_l_sq * m_lprime_sq
            * z_gamma_sq - k6 * q2 * z_gamma_sq + 2 * k4 * m_B_sq * q2 * z_gamma_sq - k2 * m_B_4
            * q2 * z_gamma_sq - k4 * m_l_sq * q2 * z_gamma_sq + 2 * k2 * m_B_sq * m_l_sq * q2
            * z_gamma_sq - m_B_4 * m_l_sq * q2 * z_gamma_sq - 24 * k4 * m_lprime_sq * q2
            * z_gamma_sq - 8 * k2 * m_B_sq * m_lprime_sq * q2 * z_gamma_sq + 8 * k2 * m_l_sq
            * m_lprime_sq * q2 * z_gamma_sq - 8 * m_B_sq * m_l_sq * m_lprime_sq * q2 * z_gamma_sq
            + 6 * k4 * q4 * z_gamma_sq + 2 * k2 * m_B_sq * q4 * z_gamma_sq - 2 * k2 * m_l_sq * q4
            * z_gamma_sq + 2 * m_B_sq * m_l_sq * q4 * z_gamma_sq + 4 * k2 * m_lprime_sq * q4
            * z_gamma_sq + 4 * m_l_sq * m_lprime_sq * q4 * z_gamma_sq - k2 * q6 * z_gamma_sq
            - m_l_sq * q6 * z_gamma_sq - 4 * k6 * m_lprime_sq * z_w_sq + 8 * k4 * m_B_sq
            * m_lprime_sq * z_w_sq - 4 * k2 * m_B_4 * m_lprime_sq * z_w_sq + 4 * k4 * m_l_sq
            * m_lprime_sq * z_w_sq - 8 * k2 * m_B_sq * m_l_sq * m_lprime_sq * z_w_sq + 4 * m_B_4
            * m_l_sq * m_lprime_sq * z_w_sq - k6 * q2 * z_w_sq + 2 * k4 * m_B_sq * q2 * z_w_sq - k2
            * m_B_4 * q2 * z_w_sq + k4 * m_l_sq * q2 * z_w_sq - 2 * k2 * m_B_sq * m_l_sq * q2
            * z_w_sq + m_B_4 * m_l_sq * q2 * z_w_sq - 8 * k4 * m_lprime_sq * q2 * z_w_sq + 8 * k2
            * m_B_sq * m_lprime_sq * q2 * z_w_sq + 8 * k2 * m_l_sq * m_lprime_sq * q2 * z_w_sq - 8
            * m_B_sq * m_l_sq * m_lprime_sq * q2 * z_w_sq + 6 * k4 * q4 * z_w_sq + 2 * k2 * m_B_sq
            * q4 * z_w_sq - 6 * k2 * m_l_sq * q4 * z_w_sq - 2 * m_B_sq * m_l_sq * q4 * z_w_sq - 4
            * k2 * m_lprime_sq * q4 * z_w_sq + 4 * m_l_sq * m_lprime_sq * q4 * z_w_sq - k2 * q6
            * z_w_sq + m_l_sq * q6 * z_w_sq + 4 * k6 * m_lprime_sq * z_gamma_sq * z_w_sq
            - 8 * k4 * m_B_sq * m_lprime_sq * z_gamma_sq * z_w_sq + 4 * k2 * m_B_4 * m_lprime_sq
            * z_gamma_sq * z_w_sq - 4 * k4 * m_l_sq * m_lprime_sq * z_gamma_sq * z_w_sq + 8 * k2
            * m_B_sq * m_l_sq * m_lprime_sq * z_gamma_sq * z_w_sq - 4 * m_B_4 * m_l_sq
            * m_lprime_sq * z_gamma_sq * z_w_sq - k6 * q2 * z_gamma_sq * z_w_sq + 2 * k4 * m_B_sq
            * q2 * z_gamma_sq * z_w_sq - k2 * m_B_4 * q2 * z_gamma_sq * z_w_sq + k4 * m_l_sq * q2
            * z_gamma_sq * z_w_sq - 2 * k2 * m_B_sq * m_l_sq * q2 * z_gamma_sq * z_w_sq + m_B_4
            * m_l_sq * q2 * z_gamma_sq * z_w_sq + 40 * k4 * m_lprime_sq * q2 * z_gamma_sq * z_w_sq
            - 8 * k2 * m_B_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 40 * k2 * m_l_sq
            * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + 8 * m_B_sq * m_l_sq * m_lprime_sq * q2
            * z_gamma_sq * z_w_sq - 10 * k4 * q4 * z_gamma_sq * z_w_sq + 2 * k2 * m_B_sq * q4
            * z_gamma_sq * z_w_sq + 10 * k2 * m_l_sq * q4 * z_gamma_sq * z_w_sq - 2 * m_B_sq
            * m_l_sq * q4 * z_gamma_sq * z_w_sq + 4 * k2 * m_lprime_sq * q4 * z_gamma_sq * z_w_sq
            - 4 * m_l_sq * m_lprime_sq * q4 * z_gamma_sq * z_w_sq - k2 * q6 * z_gamma_sq * z_w_sq
            + m_l_sq * q6 * z_gamma_sq * z_w_sq - 8 * sqrt(k2) * (k2 - m_l_sq) * sqrt(q2)
            * (k2 - m_B_sq + q2) * (- 4 * m_lprime_sq + q2) * z_gamma * sqrt(1 - z_gamma_sq) * z_w
            * sqrt(1 - z_w_sq) * cos(phi) + (k2 - m_l_sq) * (4 * m_lprime_sq - q2)
            * power_of<2>(k2 - m_B_sq + q2) * (- 1 + z_gamma_sq) * (- 1 + z_w_sq)
            * cos(2 * phi))) / (k2 * m_B_sq));
            const double f22 = (2 * (k2 - m_l_sq) * q2 * (k6 * q2 - 2 * k4 * m_B_sq * q2 + k2
            * m_B_4 * q2 + 8 * k4 * m_lprime_sq * q2 + 8 * k2 * m_l_sq * m_lprime_sq * q2 + 4 * k4
            * q4 - 2 * k2 * m_B_sq * q4 + 2 * k2 * m_l_sq * q4 + k2 * q6 + 4 * k6 * m_lprime_sq
            * z_gamma_sq - 8 * k4 * m_B_sq * m_lprime_sq * z_gamma_sq + 4 * k2 * m_B_4
            * m_lprime_sq * z_gamma_sq - k6 * q2 * z_gamma_sq + 2 * k4 * m_B_sq * q2 * z_gamma_sq
            - k2 * m_B_4 * q2 * z_gamma_sq - 8 * k2 * m_B_sq * m_lprime_sq * q2 * z_gamma_sq - 8
            * k2 * m_l_sq * m_lprime_sq * q2 * z_gamma_sq + 2 * k2 * m_B_sq * q4 * z_gamma_sq + 2
            * k2 * m_l_sq * q4 * z_gamma_sq + 4 * k2 * m_lprime_sq * q4 * z_gamma_sq
            - k2 * q6 * z_gamma_sq - k6 * q2 * z_w_sq + 2 * k4 * m_B_sq * q2 * z_w_sq - k2 * m_B_4
            * q2 * z_w_sq + k4 * m_l_sq * q2 * z_w_sq - 2 * k2 * m_B_sq * m_l_sq * q2 * z_w_sq
            + m_B_4 * m_l_sq * q2 * z_w_sq + 8 * k4 * m_lprime_sq * q2 * z_w_sq - 8 * k2 * m_l_sq
            * m_lprime_sq * q2 * z_w_sq + 2 * k2 * m_B_sq * q4 * z_w_sq - 2 * m_B_sq * m_l_sq * q4
            * z_w_sq - k2 * q6 * z_w_sq + m_l_sq * q6 * z_w_sq - 4 * k6 * m_lprime_sq * z_gamma_sq
            * z_w_sq + 8 * k4 * m_B_sq * m_lprime_sq * z_gamma_sq * z_w_sq - 4 * k2 * m_B_4
            * m_lprime_sq * z_gamma_sq * z_w_sq + 4 * k4 * m_l_sq * m_lprime_sq * z_gamma_sq
            * z_w_sq - 8 * k2 * m_B_sq * m_l_sq * m_lprime_sq * z_gamma_sq * z_w_sq + 4 * m_B_4
            * m_l_sq * m_lprime_sq * z_gamma_sq * z_w_sq + k6 * q2 * z_gamma_sq * z_w_sq - 2 * k4
            * m_B_sq * q2 * z_gamma_sq * z_w_sq + k2 * m_B_4 * q2 * z_gamma_sq * z_w_sq
            - k4 * m_l_sq * q2 * z_gamma_sq * z_w_sq + 2 * k2 * m_B_sq * m_l_sq * q2 * z_gamma_sq
            * z_w_sq - m_B_4 * m_l_sq * q2 * z_gamma_sq * z_w_sq - 16 * k4 * m_lprime_sq * q2
            * z_gamma_sq * z_w_sq + 8 * k2 * m_B_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + 16
            * k2 * m_l_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 8 * m_B_sq * m_l_sq
            * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + 4 * k4 * q4 * z_gamma_sq * z_w_sq - 2 * k2
            * m_B_sq * q4 * z_gamma_sq * z_w_sq - 4 * k2 * m_l_sq * q4 * z_gamma_sq * z_w_sq + 2
            * m_B_sq * m_l_sq * q4 * z_gamma_sq * z_w_sq - 4 * k2 * m_lprime_sq * q4 * z_gamma_sq
            * z_w_sq + 4 * m_l_sq * m_lprime_sq * q4 * z_gamma_sq * z_w_sq + k2 * q6 * z_gamma_sq
            * z_w_sq - m_l_sq * q6 * z_gamma_sq * z_w_sq +  4 * sqrt(k2) * (k2 - m_l_sq) * sqrt(q2)
            * (k2 - m_B_sq + q2) * (- 4 * m_lprime_sq + q2) * z_gamma * sqrt(1 - z_gamma_sq) * z_w
            * sqrt(1 - z_w_sq) * cos(phi) - 2 * k2 * (k2 - m_l_sq) * (4 * m_lprime_sq - q2)
            * q2 * (- 1 + z_gamma_sq) * (- 1 + z_w_sq) * cos(2 * phi))) / (k4 * m_B_sq);
            const double f33 = (2 * m_l_sq * (- k2 + m_l_sq) * q2 * (k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * (- 4 * m_lprime_sq * z_gamma_sq + q2 * (- 1 + z_gamma_sq)))
            / (k4 * m_B_sq);
            const double f44 = ((k2 - m_l_sq) * (k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) *
            (- ((4 * m_lprime_sq * (- 1 + z_gamma_sq) - q2 * (1 + z_gamma_sq)) * (k2 + m_l_sq + k2
            * z_w_sq - m_l_sq * z_w_sq)) + (k2 - m_l_sq) * (4 * m_lprime_sq - q2) * (- 1
            + z_gamma_sq) * (- 1 + z_w_sq) * cos(2 * phi))) / (k2 * m_B_sq);
            const double f2c1Re =(4 * (k2 - m_l_sq) * sqrt(q2) * ((k2 - m_l_sq) * (k4 - 2 * k2
            * (m_B_sq - 3 * q2) + m_B_sq_q2_diff2) * (- 4 * m_lprime_sq + q2) * z_gamma
            * sqrt(1 - z_gamma_sq) * z_w * sqrt(1 - z_w_sq) * cos(phi) + sqrt(k2) * sqrt(q2)
            * (k2 - m_B_sq + q2) * (m_l_sq * (q2 * (1 + z_w_sq + z_gamma_sq * (1 - 3 * z_w_sq))
            + 4 * m_lprime_sq * (1 - z_w_sq + z_gamma_sq * (- 1 + 3 * z_w_sq)))
            + k2 * (4 * m_lprime_sq * (1 + z_w_sq + z_gamma_sq * (1 - 3 * z_w_sq)) + q2 * (3
            - z_w_sq + z_gamma_sq * (- 1 + 3 * z_w_sq))) - (k2 - m_l_sq) * (4 * m_lprime_sq - q2)
            * (- 1 + z_gamma_sq) * (- 1 + z_w_sq) * cos(2 * phi)))) / (std::pow(k2,1.5) * m_B_sq);
            const double f3c1Re =(4 * m_l_sq * (k2 - m_l_sq) * sqrt(q2) * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * (2 * sqrt(k2) * sqrt(q2)
            * (- 4 * m_lprime_sq * z_gamma_sq + q2 * (- 1 + z_gamma_sq)) * z_w - (4 * m_lprime_sq
            - q2) * (k2 - m_B_sq + q2) * z_gamma * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * cos(phi))) / (std::pow(k2,1.5) * m_B_sq);
            const double f4c1Re = (4 * (k2 - m_l_sq) * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * ((k2 - m_B_sq + q2) * (- 4 * m_lprime_sq * (- 1 + z_gamma_sq) + q2
            * (1 + z_gamma_sq)) * z_w + 2 * sqrt(k2) * sqrt(q2) * (- 4 * m_lprime_sq + q2) * z_gamma
            * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * cos(phi))) / m_B_sq;
            const double f3c2Re = (4 * m_l_sq * (k2 - m_l_sq) * q2 * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * ((k2 - m_B_sq + q2) * (- 4 * m_lprime_sq * z_gamma_sq + q2
            * (- 1 + z_gamma_sq)) * z_w + 2 * sqrt(k2) * sqrt(q2) * (- 4 * m_lprime_sq + q2)
            * z_gamma * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * cos(phi))) / (k4 * m_B_sq);
            const double f4c2Re = (4 * (k2 - m_l_sq) * sqrt(q2) * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * (2 * sqrt(k2) * sqrt(q2) * (- 4 * m_lprime_sq
            * (- 1 + z_gamma_sq) + q2 * (1 + z_gamma_sq)) * z_w - (4 * m_lprime_sq - q2) * (k2
            - m_B_sq + q2) * z_gamma * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * cos(phi)))
            / (sqrt(k2) * m_B_sq);
            const double f2c1Im = (4.0 * (k2 - m_l_sq) * sqrt(q2) * (- 4 * m_lprime_sq + q2) * (k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * sin(phi)) / (sqrt(k2) * m_B_sq);
            const double f4c1Im = (- 4.0 * m_l_sq_k2_diff2 * (- 4 * m_lprime_sq + q2) * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) *
            (2 * sqrt(k2) * sqrt(q2) * z_gamma * z_w + (k2 - m_B_sq + q2) * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * cos(phi)) * sin(phi)) / (k2 * m_B_sq);
            const double f4c2Im = (- 4.0 * m_l_sq_k2_diff2 * sqrt(q2) * (- 4 * m_lprime_sq + q2)
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * ((k2 - m_B_sq + q2) * z_gamma * z_w + 2 * sqrt(k2) * sqrt(q2)
            * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * cos(phi)) * sin(phi)) / (std::pow(k2,1.5)
            * m_B_sq);
            const double f4c3Im =(4.0 * m_l_sq * (- k2 + m_l_sq) * sqrt(q2) * (- 4 * m_lprime_sq
            + q2) * (k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * sin(phi)) / (std::pow(k2,1.5) * m_B_sq);
            const double f51Re = (8 * m_l_sq * (-k2 + m_l_sq) * (-(sqrt(q2) * (-4 * m_lprime_sq
            + q2) * z_gamma * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * (-7 * k6 * z_w + k4
            * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + 6 * m_B_sq * z_w + 7 * m_l_sq
            * z_w - 2 * q2 * z_w) - m_l_sq * (m_B_sq - q2) * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) + m_B_sq * z_w - q2 * z_w) + k2 * (m_l_sq * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) - q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            + m_B_4 * z_w + 2 * m_l_sq * q2 * z_w + q4 * z_w - m_B_sq * (sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) + 6 * m_l_sq * z_w + 2 * q2 * z_w))) * cos(phi)) + sqrt(k2)
            * (-4 * m_B_6 * m_lprime_sq + 4 * m_B_4 * m_l_sq * m_lprime_sq - m_B_6 * q2 + m_B_4
            * m_l_sq * q2 + 8 * m_B_4 * m_lprime_sq * q2 - 4 * m_B_sq * m_l_sq * m_lprime_sq * q2
            + 2 * m_B_4 * q4 - m_B_sq * m_l_sq * q4 - 4 * m_B_sq * m_lprime_sq * q4 - m_B_sq * q6 +
            4 * m_B_6 * m_lprime_sq * z_gamma_sq - 4 * m_B_4 * m_l_sq * m_lprime_sq * z_gamma_sq
            - m_B_6 * q2 * z_gamma_sq + m_B_4 * m_l_sq * q2 * z_gamma_sq - 8 * m_B_4 * m_lprime_sq
            * q2 * z_gamma_sq + 4 * m_B_sq * m_l_sq * m_lprime_sq * q2 * z_gamma_sq + 2 * m_B_4
            * q4 * z_gamma_sq - m_B_sq * m_l_sq * q4 * z_gamma_sq + 4 * m_B_sq * m_lprime_sq * q4
            * z_gamma_sq - m_B_sq * q6 * z_gamma_sq - 4 * m_B_4 * m_lprime_sq
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w - m_B_4 * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 4 * m_B_sq * m_lprime_sq * q2
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + m_B_sq * q4
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 2 * m_l_sq * q4
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 4 * m_B_4 * m_lprime_sq
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w - m_B_4 * q2
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w - 4 * m_B_sq
            * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq
            * z_w + 8 * m_l_sq * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w + m_B_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w - 2 * m_l_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w - 4 * m_B_4 * m_l_sq * m_lprime_sq * z_w_sq - m_B_4
            * m_l_sq * q2 * z_w_sq + 4 * m_B_sq * m_l_sq * m_lprime_sq * q2 * z_w_sq + 3 * m_B_sq
            * m_l_sq * q4 * z_w_sq - 2 * m_l_sq * q6 * z_w_sq + 4 * m_B_4 * m_l_sq * m_lprime_sq
            * z_gamma_sq * z_w_sq - m_B_4 * m_l_sq * q2 * z_gamma_sq * z_w_sq + 4 * m_B_sq * m_l_sq
            * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - m_B_sq * m_l_sq * q4 * z_gamma_sq * z_w_sq -
            8 * m_l_sq * m_lprime_sq * q4 * z_gamma_sq * z_w_sq + 2 * m_l_sq * q6 * z_gamma_sq
            * z_w_sq + k6 * (-4 * m_lprime_sq * (-1 + z_gamma_sq) + q2 * (1 + z_gamma_sq)) * z_w_sq
            + k4 * (4 * m_lprime_sq * q2 + 5 * q4 + 12 * m_lprime_sq * q2 * z_gamma_sq - 3 * q4
            * z_gamma_sq - 4 * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w - q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 4 * m_lprime_sq
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w - q2
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 4
            * m_lprime_sq * q2 * z_w_sq - 5 * q4 * z_w_sq - 28 * m_lprime_sq * q2 * z_gamma_sq
            * z_w_sq + 7 * q4 * z_gamma_sq * z_w_sq + m_l_sq * (4 * m_lprime_sq * (-1 + z_gamma_sq)
            - q2 * (1 + z_gamma_sq)) * (-1 + z_w_sq) + m_B_sq * (4 * m_lprime_sq * (-1 + z_gamma_sq)
            - q2 * (1 + z_gamma_sq)) * (1 + 2 * z_w_sq)) + k2 * (-(m_B_4 * (4 * m_lprime_sq
            * (-1 + z_gamma_sq) - q2 * (1 + z_gamma_sq)) * (2 + z_w_sq)) + m_B_sq * (-2 * m_l_sq
            * (4 * m_lprime_sq * (-1 + z_gamma_sq) - q2 * (1 + z_gamma_sq)) * (-1 + z_w_sq) +
            q2 * (2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * (1 + z_gamma_sq) * z_w
            + q2 * (-7 - 3 * z_w_sq + z_gamma_sq * (1 + z_w_sq))) - 4 * m_lprime_sq * (2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * (-1 + z_gamma_sq) * z_w + q2 * (3 + z_w_sq
            + z_gamma_sq * (1 + z_w_sq)))) + q2 * (4 * m_lprime_sq * q2 + q4 - 4 * m_lprime_sq * q2
            * z_gamma_sq + q4 * z_gamma_sq - 4 * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_w - 3 * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w - 4 * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma_sq * z_w + q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma_sq * z_w + 2 * q4 * z_w_sq + 8 * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 2
            * q4 * z_gamma_sq * z_w_sq + m_l_sq * (q2 * (1 + 5 * z_w_sq + z_gamma_sq * (1 - 7
            * z_w_sq)) + 4 * m_lprime_sq * (1 - z_w_sq + z_gamma_sq * (-1 + 7 * z_w_sq))))) - (k2
            - m_B_sq) * (k2 - m_l_sq) * (4 * m_lprime_sq - q2) * (k2 - m_B_sq + q2) * (-1
            + z_gamma_sq) * (-1 + z_w_sq) * cos(2 * phi))))/(sqrt(k2) * m_B * (k2 - m_B_sq)
            * (k4 - m_l_sq * (m_B_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w) + k2 * (-m_B_sq + m_l_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_w)));
            const double f52Re = (-8 * m_l_sq * (-k2 + m_l_sq) * sqrt(q2) * (sqrt(k2) * (-4
            * m_lprime_sq + q2) * z_gamma * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) *
            (m_B_4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) - m_B_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + 2 * m_l_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) - 2 * k6 * z_w + 2 * m_B_4 * m_l_sq * z_w
            - 2 * m_l_sq * q4 * z_w + k4 * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            + 4 * m_B_sq * z_w + 2 * m_l_sq * z_w - 8 * q2 * z_w) - k2 * (2 * m_B_4 * z_w + 2
            * m_B_sq * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + 2 * m_l_sq * z_w) +
            q2 * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) - 8 * m_l_sq * z_w - 2 * q2
            * z_w))) * cos(phi) + sqrt(q2) * (-(m_l_sq * (m_B_sq - q2) * (-4 * m_lprime_sq
            * z_gamma_sq + q2 * (-1 + z_gamma_sq)) * z_w * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) + m_B_sq * z_w - q2 * z_w)) + k6 * (q2 * (-2 + z_w_sq + z_gamma_sq
            * (2 - 5 * z_w_sq)) + 4 * m_lprime_sq * (-2 * z_w_sq + z_gamma_sq * (-2 + 5 * z_w_sq)))
            + k4 * (-8 * m_l_sq * m_lprime_sq - 2 * m_l_sq * q2 - 8 * m_lprime_sq * q2 - 4 * q4 + 8
            * m_l_sq * m_lprime_sq * z_gamma_sq - 2 * m_l_sq * q2 * z_gamma_sq +
            8 * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 3 * q2
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w - 4 * m_lprime_sq * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 8 * m_l_sq
            * m_lprime_sq * z_w_sq - m_l_sq * q2 * z_w_sq + 2 * q4 * z_w_sq - 20 * m_l_sq
            * m_lprime_sq * z_gamma_sq * z_w_sq + 5 * m_l_sq * q2 * z_gamma_sq * z_w_sq + 8
            * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 2 * q4 * z_gamma_sq * z_w_sq + m_B_sq * (8
            * m_lprime_sq * (1 + z_w_sq + z_gamma_sq * (1 - 2 * z_w_sq)) + 2 * q2 * (3 + z_gamma_sq
            * (-1 + 2 * z_w_sq)))) - k2 * (-((-4 * m_lprime_sq * z_gamma_sq + q2 * (-1
            + z_gamma_sq)) * z_w * (q2 * (-sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + q2
            * z_w) + m_l_sq * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + 2 * q2 * z_w)))
            + m_B_4 * (4 * m_lprime_sq * (2 + z_gamma_sq * z_w_sq) + q2 * (4 - (-1 + z_gamma_sq)
            * z_w_sq)) + m_B_sq * (-8 * m_lprime_sq * q2 - 4 * q4 + 8 * m_lprime_sq * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 3 * q2 * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * z_w - 4 * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w + q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w - 2 * q4 * z_w_sq - 8 * m_lprime_sq * q2 * z_gamma_sq * z_w_sq
            + 2 * q4 * z_gamma_sq * z_w_sq - 2 * m_l_sq * (q2 + q2 * z_gamma_sq * (1 - 2 * z_w_sq)
            + m_lprime_sq * (4 - 4 * z_w_sq + z_gamma_sq * (-4 + 8 * z_w_sq))))) + 2 * k2 * (k2
            - m_B_sq) * (k2 - m_l_sq) * (4 * m_lprime_sq - q2) * (-1 + z_gamma_sq) * (-1 + z_w_sq)
            * cos(2 * phi))))/(k2 * m_B * (k2 - m_B_sq) * (k4 - m_l_sq * (m_B_sq - q2 + sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) + k2 * (-m_B_sq + m_l_sq - q2
            + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)));
            const double f53Re = (8 * m_l_sq * (-k2 + m_l_sq) * sqrt(q2) * (-(sqrt(q2) * (-4
            * m_lprime_sq * z_gamma_sq + q2 * (-1 + z_gamma_sq)) * (k6 - k4 * (2 * m_B_sq + m_l_sq
            + 2 * q2 - sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) -
            m_l_sq * (m_B_sq - q2) * (m_B_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_w) + k2 * (m_B_4 + 2 * m_l_sq * q2 + q4 - 3 * m_l_sq * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w - q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2
            * k2 * (m_B_sq + q2)) * z_w + m_B_sq * (2 * m_l_sq - 2 * q2 + 3 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)))) - sqrt(k2) * (k2 - m_B_sq)
            * (k2 + m_B_sq - 2 * m_l_sq - q2) * (-4 * m_lprime_sq + q2) * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * z_gamma * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq)
            * cos(phi)))/(k2 * m_B * (k2 - m_B_sq) * (k4 - m_l_sq * (m_B_sq - q2 + sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) + k2 * (-m_B_sq + m_l_sq - q2
            + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)));
            const double f54Re = (-8 * m_l_sq * (-k2 + m_l_sq) * ((4 * m_lprime_sq * (-1
            + z_gamma_sq) - q2 * (1 + z_gamma_sq)) * (k4 + (m_B_sq - q2) * (m_B_sq - q2 + sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) - k2 * (2 * m_B_sq + 2 * q2
            + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)) + 2 * sqrt(k2) * sqrt(q2)
            * (-4 * m_lprime_sq + q2) * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * cos(phi)))/(m_B * (-k4 + k2
            * (m_B_sq - m_l_sq + q2 - sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) +
            m_l_sq * (m_B_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)));
            const double f52Im = (-8.0 * m_l_sq * (k2 - m_l_sq) * sqrt(q2) * (-4 * m_lprime_sq + q2)
            * (k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * sin(phi))/(sqrt(k2) * m_B * (k4 - m_l_sq * (m_B_sq - q2 + sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) + k2 * (-m_B_sq + m_l_sq - q2
            + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)));
            const double f53Im = (-8.0 * m_l_sq * (-k2 + m_l_sq) * sqrt(q2) * (-4
            * m_lprime_sq + q2) * (k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma
            * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) * sin(phi))/(sqrt(k2) * m_B * (k4 - m_l_sq
            * (m_B_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w) + k2
            * (-m_B_sq + m_l_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)));
            const double f54Im = (8.0 * m_l_sq * m_l_sq_k2_diff2 * (-4 * m_lprime_sq + q2) * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * sqrt(1 - z_gamma_sq) * sqrt(1 - z_w_sq) *
            (sqrt(q2) * z_gamma * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + 3 * k2
            * z_w + m_B_sq * z_w - q2 * z_w) + 2 * sqrt(k2) * (k2 - m_B_sq) * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * cos(phi)) * sin(phi))/(sqrt(k2) * m_B * (k2 - m_B_sq) * (k4
            - m_l_sq * (m_B_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)
            + k2 * (-m_B_sq + m_l_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w)));
            const double f55 = (8 * m_l_sq * (k2 - m_l_sq) * (4 * k10 * m_lprime_sq - 8 * k8
            * m_B_sq * m_lprime_sq + 8 * k6 * m_B_4 * m_lprime_sq - 8 * k4 * m_B_6 * m_lprime_sq +
            4 * k2 * m_B_8 * m_lprime_sq - 4 * k8 * m_l_sq * m_lprime_sq - 8 * k6 * m_B_sq * m_l_sq
            * m_lprime_sq + 32 * k4 * m_B_4 * m_l_sq * m_lprime_sq - 24 * k2 * m_B_6 * m_l_sq
            * m_lprime_sq + 4 * m_B_8 * m_l_sq * m_lprime_sq + 8 * k6 * m_l_4 * m_lprime_sq - 16
            * k4 * m_B_sq * m_l_4 * m_lprime_sq + 8 * k2 * m_B_4 * m_l_4 * m_lprime_sq + k10 * q2
            - 2 * k8 * m_B_sq * q2 + 2 * k6 * m_B_4 * q2 - 2 * k4 * m_B_6 * q2 + k2 * m_B_8 * q2
            - k8 * m_l_sq * q2 - 2 * k6 * m_B_sq * m_l_sq * q2 + 8 * k4 * m_B_4 * m_l_sq * q2
            - 6 * k2 * m_B_6 * m_l_sq * q2 + m_B_8 * m_l_sq * q2 + 2 * k6 * m_l_4 * q2 - 4 * k4
            * m_B_sq * m_l_4 * q2 + 2 * k2 * m_B_4 * m_l_4 * q2 - 8 * k8 * m_lprime_sq * q2 + 8 * k6
            * m_B_sq * m_lprime_sq * q2 + 8 * k4 * m_B_4 * m_lprime_sq * q2 - 8 * k2 * m_B_6
            * m_lprime_sq * q2 + 8 * k6 * m_l_sq * m_lprime_sq * q2 - 24 * k4 * m_B_sq * m_l_sq
            * m_lprime_sq * q2 + 24 * k2 * m_B_4 * m_l_sq * m_lprime_sq * q2 - 8 * m_B_6 * m_l_sq
            * m_lprime_sq * q2 - k8 * q4 + 3 * k4 * m_B_4 * q4 - 2 * k2 * m_B_6 * q4 + 4 * k6
            * m_l_sq * q4 - 10 * k4 * m_B_sq * m_l_sq * q4 + 8 * k2 * m_B_4 * m_l_sq * q4 - 2
            * m_B_6 * m_l_sq * q4 + k4 * m_l_4 * q4 - 2 * k2 * m_B_sq * m_l_4 * q4 + m_B_4 * m_l_4
            * q4 + 4 * k6 * m_lprime_sq * q4 - 8 * k4 * m_B_sq * m_lprime_sq * q4 + 4 * k2 * m_B_4
            * m_lprime_sq * q4 + 4 * k4 * m_l_sq * m_lprime_sq * q4 - 8 * k2 * m_B_sq * m_l_sq
            * m_lprime_sq * q4 + 4 * m_B_4 * m_l_sq * m_lprime_sq * q4 - k6 * q6 - 4 * k4 * m_B_sq
            * q6 + k2 * m_B_4 * q6 + 5 * k4 * m_l_sq * q6 + 2 * k2 * m_B_sq * m_l_sq * q6 + m_B_4
            * m_l_sq * q6 - 2 * k2 * m_l_4 * q6 - 2 * m_B_sq * m_l_4 * q6 + k4 * q8 - 2 * k2
            * m_l_sq * q8 + m_l_4 * q8 - 4 * k10 * m_lprime_sq * z_gamma_sq + 8 * k8 * m_B_sq
            * m_lprime_sq * z_gamma_sq - 8 * k6 * m_B_4 * m_lprime_sq * z_gamma_sq +
            8 * k4 * m_B_6 * m_lprime_sq * z_gamma_sq - 4 * k2 * m_B_8 * m_lprime_sq * z_gamma_sq
            + 4 * k8 * m_l_sq * m_lprime_sq * z_gamma_sq + 8 * k6 * m_B_sq * m_l_sq * m_lprime_sq
            * z_gamma_sq - 32 * k4 * m_B_4 * m_l_sq * m_lprime_sq * z_gamma_sq + 24 * k2 * m_B_6
            * m_l_sq * m_lprime_sq * z_gamma_sq - 4 * m_B_8 * m_l_sq * m_lprime_sq * z_gamma_sq
            - 8 * k6 * m_l_4 * m_lprime_sq * z_gamma_sq + 16 * k4 * m_B_sq * m_l_4 * m_lprime_sq
            * z_gamma_sq - 8 * k2 * m_B_4 * m_l_4 * m_lprime_sq * z_gamma_sq + k10 * q2 * z_gamma_sq
            - 2 * k8 * m_B_sq * q2 * z_gamma_sq + 2 * k6 * m_B_4 * q2 * z_gamma_sq - 2 * k4 * m_B_6
            * q2 * z_gamma_sq + k2 * m_B_8 * q2 * z_gamma_sq - k8 * m_l_sq * q2 * z_gamma_sq - 2
            * k6 * m_B_sq * m_l_sq * q2 * z_gamma_sq + 8 * k4 * m_B_4 * m_l_sq * q2 * z_gamma_sq -
            6 * k2 * m_B_6 * m_l_sq * q2 * z_gamma_sq + m_B_8 * m_l_sq * q2 * z_gamma_sq + 2 * k6
            * m_l_4 * q2 * z_gamma_sq - 4 * k4 * m_B_sq * m_l_4 * q2 * z_gamma_sq +
            2 * k2 * m_B_4 * m_l_4 * q2 * z_gamma_sq + 12 * k8 * m_lprime_sq * q2 * z_gamma_sq - 16
            * k6 * m_B_sq * m_lprime_sq * q2 * z_gamma_sq - 4 * k4 * m_B_4 * m_lprime_sq * q2
            * z_gamma_sq + 8 * k2 * m_B_6 * m_lprime_sq * q2 * z_gamma_sq + 8 * k4 * m_B_sq * m_l_sq
            * m_lprime_sq * q2 * z_gamma_sq - 16 * k2 * m_B_4 * m_l_sq * m_lprime_sq * q2
            * z_gamma_sq + 8 * m_B_6 * m_l_sq * m_lprime_sq * q2 * z_gamma_sq + 4 * k4 * m_l_4
            * m_lprime_sq * q2 * z_gamma_sq - 8 * k2 * m_B_sq * m_l_4 * m_lprime_sq * q2
            * z_gamma_sq + 4 * m_B_4 * m_l_4 * m_lprime_sq * q2 * z_gamma_sq - 3 * k8 * q4
            * z_gamma_sq + 4 * k6 * m_B_sq * q4 * z_gamma_sq + k4 * m_B_4 * q4 * z_gamma_sq -
            2 * k2 * m_B_6 * q4 * z_gamma_sq - 2 * k4 * m_B_sq * m_l_sq * q4 * z_gamma_sq + 4 * k2
            * m_B_4 * m_l_sq * q4 * z_gamma_sq - 2 * m_B_6 * m_l_sq * q4 * z_gamma_sq -
            k4 * m_l_4 * q4 * z_gamma_sq + 2 * k2 * m_B_sq * m_l_4 * q4 * z_gamma_sq - m_B_4 * m_l_4
            * q4 * z_gamma_sq - 12 * k6 * m_lprime_sq * q4 * z_gamma_sq - 4 * k2 * m_B_4
            * m_lprime_sq * q4 * z_gamma_sq + 12 * k4 * m_l_sq * m_lprime_sq * q4 * z_gamma_sq + 24
            * k2 * m_B_sq * m_l_sq * m_lprime_sq * q4 * z_gamma_sq - 4 * m_B_4 * m_l_sq
            * m_lprime_sq * q4 * z_gamma_sq - 8 * k2 * m_l_4 * m_lprime_sq * q4 * z_gamma_sq - 8
            * m_B_sq * m_l_4 * m_lprime_sq * q4 * z_gamma_sq + 3 * k6 * q6 * z_gamma_sq + k2 * m_B_4
            * q6 * z_gamma_sq - 3 * k4 * m_l_sq * q6 * z_gamma_sq - 6 * k2 * m_B_sq * m_l_sq * q6
            * z_gamma_sq + m_B_4 * m_l_sq * q6 * z_gamma_sq + 2 * k2 * m_l_4 * q6 * z_gamma_sq + 2
            * m_B_sq * m_l_4 * q6 * z_gamma_sq + 4 * k4 * m_lprime_sq * q6 * z_gamma_sq - 8 * k2
            * m_l_sq * m_lprime_sq * q6 * z_gamma_sq + 4 * m_l_4 * m_lprime_sq * q6 * z_gamma_sq
            - k4 * q8 * z_gamma_sq + 2 * k2 * m_l_sq * q8 * z_gamma_sq - m_l_4 * q8 * z_gamma_sq
            - 8 * k6 * m_l_sq * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w + 24 * k4 * m_B_sq * m_l_sq * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_w - 24 * k2 * m_B_4 * m_l_sq * m_lprime_sq * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 8 * m_B_6 * m_l_sq * m_lprime_sq
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w -
            2 * k6 * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w +
            6 * k4 * m_B_sq * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w - 6 * k2 * m_B_4 * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_w + 2 * m_B_6 * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_w - 8 * k4 * m_l_sq * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_w + 16 * k2 * m_B_sq * m_l_sq * m_lprime_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w - 8 * m_B_4 * m_l_sq * m_lprime_sq
            * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w +
            2 * k6 * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w +
            6 * k4 * m_B_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w -
            10 * k4 * m_l_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w -
            4 * k2 * m_B_sq * m_l_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w - 2 * m_B_4 * m_l_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_w + 6 * k2 * m_l_4 * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w
            + 2 * m_B_sq * m_l_4 * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w -
            2 * k4 * q6 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w + 4 * k2
            * m_l_sq * q6 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w -
            2 * m_l_4 * q6 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w +
            8 * k6 * m_l_sq * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma_sq * z_w - 24 * k4 * m_B_sq * m_l_sq * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 24 * k2 * m_B_4 * m_l_sq * m_lprime_sq
            * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w -
            8 * m_B_6 * m_l_sq * m_lprime_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma_sq * z_w - 2 * k6 * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w + 6 * k4 * m_B_sq * m_l_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w -
            6 * k2 * m_B_4 * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma_sq * z_w + 2 * m_B_6 * m_l_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w + 8 * k6 * m_lprime_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w +
            24 * k4 * m_B_sq * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w - 24 * k4 * m_l_sq * m_lprime_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w -
            48 * k2 * m_B_sq * m_l_sq * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w + 8 * m_B_4 * m_l_sq * m_lprime_sq * q2 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 24 * k2 * m_l_4
            * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq
            * z_w + 8 * m_B_sq * m_l_4 * m_lprime_sq * q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w - 2 * k6 * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w - 6 * k4 * m_B_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 6 * k4 * m_l_sq * q4 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 12 * k2 * m_B_sq
            * m_l_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w -
            2 * m_B_4 * m_l_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            * z_gamma_sq * z_w - 6 * k2 * m_l_4 * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w - 2 * m_B_sq * m_l_4 * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w - 8 * k4 * m_lprime_sq * q4 * sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w + 16 * k2 * m_l_sq
            * m_lprime_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_gamma_sq
            * z_w - 8 * m_l_4 * m_lprime_sq * q4 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w + 2 * k4 * q6 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq
            + q2)) * z_gamma_sq * z_w - 4 * k2 * m_l_sq * q6 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) * z_gamma_sq * z_w + 2 * m_l_4 * q6 * sqrt(k4 + m_B_sq_q2_diff2 - 2
            * k2 * (m_B_sq + q2)) * z_gamma_sq * z_w - 4 * k10 * m_lprime_sq * z_w_sq + 8 * k8
            * m_B_sq * m_lprime_sq * z_w_sq - 8 * k6 * m_B_4 * m_lprime_sq * z_w_sq + 8 * k4 * m_B_6
            * m_lprime_sq * z_w_sq - 4 * k2 * m_B_8 * m_lprime_sq * z_w_sq + 12 * k8 * m_l_sq
            * m_lprime_sq * z_w_sq - 24 * k6 * m_B_sq * m_l_sq * m_lprime_sq * z_w_sq + 16 * k4
            * m_B_4 * m_l_sq * m_lprime_sq * z_w_sq - 8 * k2 * m_B_6 * m_l_sq * m_lprime_sq * z_w_sq
            + 4 * m_B_8 * m_l_sq * m_lprime_sq * z_w_sq - 8 * k6 * m_l_4 * m_lprime_sq * z_w_sq + 16
            * k4 * m_B_sq * m_l_4 * m_lprime_sq * z_w_sq - 8 * k2 * m_B_4 * m_l_4 * m_lprime_sq
            * z_w_sq - k10 * q2 * z_w_sq + 2 * k8 * m_B_sq * q2 * z_w_sq - 2 * k6 * m_B_4 * q2
            * z_w_sq + 2 * k4 * m_B_6 * q2 * z_w_sq - k2 * m_B_8 * q2 * z_w_sq + 3 * k8 * m_l_sq
            * q2 * z_w_sq - 6 * k6 * m_B_sq * m_l_sq * q2 * z_w_sq + 4 * k4 * m_B_4 * m_l_sq * q2
            * z_w_sq - 2 * k2 * m_B_6 * m_l_sq * q2 * z_w_sq + m_B_8 * m_l_sq * q2 * z_w_sq - 2 * k6
            * m_l_4 * q2 * z_w_sq + 4 * k4 * m_B_sq * m_l_4 * q2 * z_w_sq - 2 * k2 * m_B_4 * m_l_4
            * q2 * z_w_sq + 8 * k8 * m_lprime_sq * q2 * z_w_sq - 8 * k6 * m_B_sq * m_lprime_sq * q2
            * z_w_sq - 8 * k4 * m_B_4 * m_lprime_sq * q2 * z_w_sq + 8 * k2 * m_B_6 * m_lprime_sq
            * q2 * z_w_sq - 8 * k6 * m_l_sq * m_lprime_sq * q2 * z_w_sq + 8 * k4 * m_B_sq * m_l_sq
            * m_lprime_sq * q2 * z_w_sq + 8 * k2 * m_B_4 * m_l_sq * m_lprime_sq * q2 * z_w_sq - 8
            * m_B_6 * m_l_sq * m_lprime_sq * q2 * z_w_sq + 3 * k8 * q4 * z_w_sq + 4 * k6 * m_B_sq
            * q4 * z_w_sq + 7 * k4 * m_B_4 * q4 * z_w_sq + 2 * k2 * m_B_6 * q4 * z_w_sq - 12 * k6
            * m_l_sq * q4 * z_w_sq - 10 * k4 * m_B_sq * m_l_sq * q4 * z_w_sq -
            8 * k2 * m_B_4 * m_l_sq * q4 * z_w_sq - 2 * m_B_6 * m_l_sq * q4 * z_w_sq + 9 * k4
            * m_l_4 * q4 * z_w_sq + 6 * k2 * m_B_sq * m_l_4 * q4 * z_w_sq + m_B_4 * m_l_4 * q4
            * z_w_sq - 4 * k6 * m_lprime_sq * q4 * z_w_sq + 8 * k4 * m_B_sq * m_lprime_sq * q4
            * z_w_sq - 4 * k2 * m_B_4 * m_lprime_sq * q4 * z_w_sq + 4 * k4 * m_l_sq * m_lprime_sq
            * q4 * z_w_sq - 8 * k2 * m_B_sq * m_l_sq * m_lprime_sq * q4 * z_w_sq +
            4 * m_B_4 * m_l_sq * m_lprime_sq * q4 * z_w_sq - 3 * k6 * q6 * z_w_sq - 4 * k4 * m_B_sq
            * q6 * z_w_sq - k2 * m_B_4 * q6 * z_w_sq + 9 * k4 * m_l_sq * q6 * z_w_sq + 6 * k2
            * m_B_sq * m_l_sq * q6 * z_w_sq + m_B_4 * m_l_sq * q6 * z_w_sq - 6 * k2 * m_l_4 * q6
            * z_w_sq - 2 * m_B_sq * m_l_4 * q6 * z_w_sq + k4 * q8 * z_w_sq - 2 * k2 * m_l_sq * q8
            * z_w_sq + m_l_4 * q8 * z_w_sq + 4 * k10 * m_lprime_sq * z_gamma_sq * z_w_sq - 8 * k8
            * m_B_sq * m_lprime_sq * z_gamma_sq * z_w_sq + 8 * k6 * m_B_4 * m_lprime_sq * z_gamma_sq
            * z_w_sq - 8 * k4 * m_B_6 * m_lprime_sq * z_gamma_sq * z_w_sq + 4 * k2 * m_B_8
            * m_lprime_sq * z_gamma_sq * z_w_sq - 12 * k8 * m_l_sq * m_lprime_sq * z_gamma_sq
            * z_w_sq + 24 * k6 * m_B_sq * m_l_sq * m_lprime_sq * z_gamma_sq * z_w_sq - 16 * k4
            * m_B_4 * m_l_sq * m_lprime_sq * z_gamma_sq * z_w_sq + 8 * k2 * m_B_6 * m_l_sq
            * m_lprime_sq * z_gamma_sq * z_w_sq - 4 * m_B_8 * m_l_sq * m_lprime_sq * z_gamma_sq
            * z_w_sq + 8 * k6 * m_l_4 * m_lprime_sq * z_gamma_sq * z_w_sq - 16 * k4 * m_B_sq
            * m_l_4 * m_lprime_sq * z_gamma_sq * z_w_sq + 8 * k2 * m_B_4 * m_l_4 * m_lprime_sq
            * z_gamma_sq * z_w_sq - k10 * q2 * z_gamma_sq * z_w_sq + 2 * k8 * m_B_sq * q2
            * z_gamma_sq * z_w_sq - 2 * k6 * m_B_4 * q2 * z_gamma_sq * z_w_sq + 2 * k4 * m_B_6 * q2
            * z_gamma_sq * z_w_sq - k2 * m_B_8 * q2 * z_gamma_sq * z_w_sq + 3 * k8 * m_l_sq * q2
            * z_gamma_sq * z_w_sq - 6 * k6 * m_B_sq * m_l_sq * q2 * z_gamma_sq * z_w_sq + 4 * k4
            * m_B_4 * m_l_sq * q2 * z_gamma_sq * z_w_sq - 2 * k2 * m_B_6 * m_l_sq * q2 * z_gamma_sq
            * z_w_sq + m_B_8 * m_l_sq * q2 * z_gamma_sq * z_w_sq - 2 * k6 * m_l_4 * q2 * z_gamma_sq
            * z_w_sq + 4 * k4 * m_B_sq * m_l_4 * q2 * z_gamma_sq * z_w_sq - 2 * k2 * m_B_4 * m_l_4
            * q2 * z_gamma_sq * z_w_sq - 4 * k8 * m_lprime_sq * q2 * z_gamma_sq * z_w_sq +
            32 * k6 * m_B_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + 44 * k4 * m_B_4
            * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 8 * k2 * m_B_6 * m_lprime_sq * q2
            * z_gamma_sq * z_w_sq - 32 * k6 * m_l_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 56
            * k4 * m_B_sq * m_l_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq - 48 * k2 * m_B_4
            * m_l_sq * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + 8 * m_B_6 * m_l_sq * m_lprime_sq
            * q2 * z_gamma_sq * z_w_sq + 36 * k4 * m_l_4 * m_lprime_sq * q2 * z_gamma_sq * z_w_sq +
            24 * k2 * m_B_sq * m_l_4 * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + 4 * m_B_4 * m_l_4
            * m_lprime_sq * q2 * z_gamma_sq * z_w_sq + k8 * q4 * z_gamma_sq * z_w_sq -
            8 * k6 * m_B_sq * q4 * z_gamma_sq * z_w_sq - 11 * k4 * m_B_4 * q4 * z_gamma_sq * z_w_sq
            + 2 * k2 * m_B_6 * q4 * z_gamma_sq * z_w_sq + 8 * k6 * m_l_sq * q4 * z_gamma_sq * z_w_sq
            + 14 * k4 * m_B_sq * m_l_sq * q4 * z_gamma_sq * z_w_sq + 12 * k2 * m_B_4 * m_l_sq * q4
            * z_gamma_sq * z_w_sq - 2 * m_B_6 * m_l_sq * q4 * z_gamma_sq * z_w_sq - 9 * k4 * m_l_4
            * q4 * z_gamma_sq * z_w_sq - 6 * k2 * m_B_sq * m_l_4 * q4 * z_gamma_sq * z_w_sq -
            m_B_4 * m_l_4 * q4 * z_gamma_sq * z_w_sq - 4 * k6 * m_lprime_sq * q4 * z_gamma_sq
            * z_w_sq - 32 * k4 * m_B_sq * m_lprime_sq * q4 * z_gamma_sq * z_w_sq + 4 * k2 * m_B_4
            * m_lprime_sq * q4 * z_gamma_sq * z_w_sq + 28 * k4 * m_l_sq * m_lprime_sq * q4
            * z_gamma_sq * z_w_sq + 40 * k2 * m_B_sq * m_l_sq * m_lprime_sq * q4 * z_gamma_sq
            * z_w_sq - 4 * m_B_4 * m_l_sq * m_lprime_sq * q4 * z_gamma_sq * z_w_sq - 24 * k2 * m_l_4
            * m_lprime_sq * q4 * z_gamma_sq * z_w_sq - 8 * m_B_sq * m_l_4 * m_lprime_sq * q4
            * z_gamma_sq * z_w_sq + k6 * q6 * z_gamma_sq * z_w_sq + 8 * k4 * m_B_sq * q6
            * z_gamma_sq * z_w_sq - k2 * m_B_4 * q6 * z_gamma_sq * z_w_sq - 7 * k4 * m_l_sq * q6
            * z_gamma_sq * z_w_sq - 10 * k2 * m_B_sq * m_l_sq * q6 * z_gamma_sq * z_w_sq + m_B_4
            * m_l_sq * q6 * z_gamma_sq * z_w_sq + 6 * k2 * m_l_4 * q6 * z_gamma_sq * z_w_sq +
            2 * m_B_sq * m_l_4 * q6 * z_gamma_sq * z_w_sq + 4 * k4 * m_lprime_sq * q6 * z_gamma_sq
            * z_w_sq - 8 * k2 * m_l_sq * m_lprime_sq * q6 * z_gamma_sq * z_w_sq + 4 * m_l_4
            * m_lprime_sq * q6 * z_gamma_sq * z_w_sq - k4 * q8 * z_gamma_sq * z_w_sq + 2 * k2
            * m_l_sq * q8 * z_gamma_sq * z_w_sq - m_l_4 * q8 * z_gamma_sq * z_w_sq - 2 * sqrt(k2)
            * (k2 - m_B_sq) * sqrt(q2) * (-4 * m_lprime_sq + q2) * z_gamma * sqrt(1 - z_gamma_sq)
            * sqrt(1 - z_w_sq) * (k6 * z_w + k4 * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2)) + 4 * m_B_sq * z_w - 7 * m_l_sq * z_w - 2 * q2 * z_w) + m_l_sq * (-3
            * m_B_sq + 2 * m_l_sq + q2) * (sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2))
            + m_B_sq * z_w - q2 * z_w) + k2 * (-(m_l_sq * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2
            * (m_B_sq + q2))) - q2 * sqrt(k4 + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) + 3
            * m_B_4 * z_w + 6 * m_l_4 * z_w + 4 * m_l_sq * q2 * z_w + q4 * z_w + m_B_sq * (sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) - 6 * m_l_sq * z_w - 4 * q2 * z_w)))
            * cos(phi) + 2 * k2 * m_B_sq_k2_diff2 * (k2 - m_l_sq) * (m_B_sq - m_l_sq) * (4
            * m_lprime_sq - q2) * (-1 + z_gamma_sq) * (-1 + z_w_sq) * cos(2 * phi)))
            /(m_B_sq_k2_diff2 * power_of<2>(k4 - m_l_sq * (m_B_sq - q2 + sqrt(k4 + m_B_sq_q2_diff2
            - 2 * k2 * (m_B_sq + q2)) * z_w) + k2 * (-m_B_sq + m_l_sq - q2 + sqrt(k4
            + m_B_sq_q2_diff2 - 2 * k2 * (m_B_sq + q2)) * z_w)));

            if ( power_of<2>(m_B - m_l) < q2 || q2 < 4 * m_lprime_sq
                || power_of<2>(m_B - sqrt(q2)) < k2 || k2 < m_l_sq )
            {
                return 0.0;
            }

            return prefactor * (f11 * norm(F_1) + f22 * norm(F_2) + f33 * norm(F_3)
                + f44 * norm(F_4) + f2c1Re * real(conj(F_2) * F_1) + f3c1Re * real(conj(F_3) * F_1)
                + f4c1Re * real(conj(F_4) * F_1) + f3c2Re * real(conj(F_3) * F_2)
                + f4c2Re * real(conj(F_4) * F_2) + f2c1Im * imag(conj(F_2) * F_1)
                + f4c1Im * imag(conj(F_4) * F_1) + f4c2Im * imag(conj(F_4) * F_2)
                + f4c3Im * imag(conj(F_4) * F_3) + f51Re * f_B * real(F_1) + f52Re * f_B * real(F_2)
                + f53Re * f_B * real(F_3) + f54Re * f_B * real(F_4) + f52Im * f_B * imag(F_2)
                + f53Im * f_B * imag(F_3) + f54Im * f_B * imag(F_4) + f55 * power_of<2>(f_B));
        }

        /*!
         * differential Branching ratio of 5 kinematic variables
         * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - m_l )^2
         * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B - sqrt(q2) )^2
         * z_gamma is the angle between the negatively charged lepton l' and the negative z-axis
         * z_w is the angle between the charged lepton l and the positive z-axis
         * phi is the angle between the q2 plane and the k2 plane
         */
        double quintuple_differential_branching_ratio(const double & q2, const double & k2,
                const double & z_gamma, const double & z_w, const double & phi) const
        {
            return quintuple_differential_decay_width(q2, k2, z_gamma, z_w, phi) * tau_B / hbar;
        }

        /*!
         * differential branching ratio of 2 kinematic variables
         * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - m_l )^2
         * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B - sqrt(q2) )^2
         */
        double double_differential_branching_ratio(const double & q2, const double & k2) const
        {
            return double_differential_decay_width(q2, k2) * tau_B / hbar;
        }

        double integrated_branching_ratio(const double & q2_min, const double & q2_max,
            const double & k2_min, const double & k2_max) const
        {
            std::function<double(const std::array<double, 2u> &)>
                integrand = [&] (const std::array<double, 2u> & x)
            {
                return this->double_differential_branching_ratio(x[0], x[1]);
            };

            auto config_cubature = cubature::Config().epsrel(10e-5);

            std::array<double, 2> x_min{ q2_min, k2_min };
            std::array<double, 2> x_max{ q2_max, k2_max };

            return integrate(integrand, x_min, x_max, config_cubature);
        }

        double _asymmetry_numerator(const double & q2,const double & k2) const
        {
            const double m_B_sq = m_B * m_B, m_B_3 = m_B_sq * m_B, m_B_4 = m_B_sq * m_B_sq;
            const double k4 = power_of<2>(k2), q4 = power_of<2>(q2);
            const double m_l_sq = m_l * m_l, m_lprime_sq = m_lprime * m_lprime;
            const double m_l_sq_k2_diff2 = power_of<2>(m_l_sq - k2);

            const complex<double> F_1 = form_factors->F_1(q2, k2);
            const complex<double> F_2 = form_factors->F_2(q2, k2);
            const complex<double> F_3 = form_factors->F_3(q2, k2);
            const complex<double> F_4 = form_factors->F_4(q2, k2);

            const double prefactor = (power_of<2>(g_fermi * abs(model->ckm_ub())
            * alpha_qed * 4 * M_PI) / (2.0 * power_of<2>(q2)) * (1.0 - m_l_sq / k2)
            * sqrt(1.0 - (4.0 * m_lprime_sq) / q2) * sqrt(power_of<2>(k2) - 2.0 * k2 * m_B_sq
            + m_B_4 - 2.0 * k2 * q2 - 2.0 * m_B_sq * q2 + power_of<2>(q2))
            / (32768.0 * m_B_3 * power_of<6>(M_PI)));

            const double g13 = (64 * m_l_sq * (k2 - m_l_sq) * M_PI * q2 * (2 * m_lprime_sq + q2)
            * sqrt(k4 - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4))
            / (3. * k2 * m_B_sq);
            const double g14 = (64 * (k2 - m_l_sq) * M_PI * (- k2 + m_B_sq - q2) * (2 * m_lprime_sq
            + q2) * sqrt(k4 - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4))
            / (3. * m_B_sq);
            const double g23 = (- 32 * m_l_sq * (k2 - m_l_sq) * M_PI * (- k2 + m_B_sq - q2) * q2
            * (2 * m_lprime_sq + q2) * sqrt(k4 - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq
            * q2 + q4)) / (3. * k4 * m_B_sq);
            const double g24 = (- 128 * (k2 - m_l_sq) * M_PI * q2 * (2 * m_lprime_sq + q2) * sqrt(k4
            - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4)) / (3. * m_B_sq);
            const double g15 = -(64 * m_l_sq * M_PI * (2 * m_lprime_sq + q2) * (2 * m_l_sq_k2_diff2
            * sqrt(k4 - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4) - (4 * k2
            * power_of<2>(- k2 + m_B_sq) * (- k2 + m_B_sq - 2 * m_l_sq - q2) * log((4 * k2
            * power_of<2>(- k2 + m_B_sq) * m_l_sq + 4 * k2 * (k2 - m_l_sq) * (m_B_sq - m_l_sq) * q2)
            / power_of<2>((k2 + m_l_sq) * (- k2 + m_B_sq - q2) + 2 * k2 * q2))) / sqrt(k4 - 2 * k2
            * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4))) / (3. * m_B * (- k2 + m_B_sq)
            * (k2 - m_l_sq));
            const double g25 = -(64 * m_l_sq * M_PI * q2 * (2 * m_lprime_sq + q2) * (m_l_sq_k2_diff2
            * sqrt(k4 - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4) + (4 * k2
            * (- k2 + m_B_sq) * (3 * k2 * (- k2 + m_B_sq) + m_l_sq_k2_diff2) * log((4 * k2
            * power_of<2>(- k2 + m_B_sq) * m_l_sq + 4 * k2 * (k2 - m_l_sq) * (m_B_sq - m_l_sq) * q2)
            / power_of<2>((k2 + m_l_sq) * (- k2 + m_B_sq - q2) + 2 * k2 * q2))) / sqrt(k4 - 2 * k2
            * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4))) / (3. * k2 * m_B
            * (- k2 + m_B_sq) * (k2 - m_l_sq));
            const double g35 = -(256 * m_l_sq * (k2 * (- k2 + m_B_sq) + (k2 - m_l_sq)
            * (k2 + m_l_sq)) * M_PI * q2 * (2 * m_lprime_sq + q2) * log((4 * k2 * power_of<2>(- k2
            + m_B_sq) * m_l_sq + 4 * k2 * (k2 - m_l_sq) * (m_B_sq - m_l_sq) * q2) / power_of<2>((k2
            + m_l_sq) * (- k2 + m_B_sq - q2) + 2 * k2 * q2))) / (3. * m_B * (k2 - m_l_sq) * sqrt(k4
            - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4));
            const double g45 = -(256 * k2 * m_l_sq * M_PI * (2 * m_lprime_sq + q2) * ((- k2 + m_B_sq)
            * (- k2 + m_B_sq - q2) - 2 * (k2 - m_l_sq) * q2) * log((4 * k2 * power_of<2>(- k2
            + m_B_sq) * m_l_sq + 4 * k2 * (k2 - m_l_sq) * (m_B_sq - m_l_sq) * q2) / power_of<2>((k2
            + m_l_sq) * (- k2 + m_B_sq - q2) + 2 * k2 * q2))) / (3. * m_B * (k2 - m_l_sq) * sqrt(k4
            - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4));
            const double g55 = (256 * k2 * m_l_sq * M_PI * (2 * m_lprime_sq + q2) * ((- k2 + m_B_sq)
            * m_l_sq_k2_diff2 * (m_B_sq - m_l_sq) * (2 * m_l_sq + q2) * sqrt(k4 - 2 * k2 * m_B_sq
            + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4) + ((4 * m_l_sq * (m_B_sq - m_l_sq) + 2
            * (power_of<2>(- k2 + m_B_sq) + (- k2 + m_B_sq) * (k2 - m_l_sq) + m_l_sq_k2_diff2)
            - (- k2 + m_B_sq) * (- k2 + m_B_sq - q2)) * ((k2 + m_l_sq) * (- k2 + m_B_sq - q2) + 2
            * k2 * q2) * (power_of<2>(- k2 + m_B_sq) * m_l_sq + (k2 - m_l_sq) * (m_B_sq - m_l_sq)
            * q2) * log((4 * k2 * power_of<2>(- k2 + m_B_sq) * m_l_sq + 4 * k2 * (k2 - m_l_sq)
            * (m_B_sq - m_l_sq) * q2) / power_of<2>((k2 + m_l_sq) * (- k2 + m_B_sq - q2) + 2 * k2
            * q2))) / sqrt(k4 - 2 * k2 * m_B_sq + m_B_4 - 2 * k2 * q2 - 2 * m_B_sq * q2 + q4)))
            / (3. * (- k2 + m_B_sq) * (k2 - m_l_sq) * ((k2 + m_l_sq) * (- k2 + m_B_sq - q2) + 2 * k2
            * q2) * (power_of<2>(- k2 + m_B_sq) * m_l_sq + (k2 - m_l_sq) * (m_B_sq - m_l_sq) * q2));

            if ( power_of<2>(m_B - m_l) < q2 || q2 < 4 * m_lprime_sq
                || power_of<2>(m_B - sqrt(q2)) < k2 || k2 < m_l_sq )
            {
                return 0.0;
            }

            return prefactor * (g13 * real(conj(F_3) * F_1) + g14 * real(conj(F_4) * F_1)
                + g23 * real(conj(F_3) * F_2) + g24 * real(conj(F_4) * F_2)
                + g15 * real(F_1) * f_B + g25 * real(F_2) * f_B + g35 * real(F_3) * f_B
                + g45 * real(F_4) * f_B + g55 * power_of<2>(f_B));
        }

        double double_differential_forward_backward_asymmetry(const double & q2,const double & k2) const
        {
            const double m_l_sq = m_l * m_l, m_lprime_sq = m_lprime * m_lprime;

            if ( power_of<2>(m_B - m_l) < q2 || q2 < 4 * m_lprime_sq
                || power_of<2>(m_B - sqrt(q2)) < k2 || k2 < m_l_sq )
            {
                return 0.0;
            }

            return _asymmetry_numerator(q2, k2) / double_differential_decay_width(q2, k2);
        }

        double integrated_forward_backward_asymmetry(const double & q2_min,const double & q2_max,
                const double & k2_min, const double & k2_max) const
        {
            std::function<double(const std::array<double, 2u> &)>
                integrand = [&] (const std::array<double, 2u> & x)
            {
                return this->_asymmetry_numerator(x[0], x[1]);
            };

            auto config_cubature = cubature::Config().epsrel(10e-5);

            std::array<double, 2> x_min{ q2_min, k2_min };
            std::array<double, 2> x_max{ q2_max, k2_max };

            return integrate(integrand, x_min, x_max, config_cubature)
                / ((hbar/tau_B) * integrated_branching_ratio(q2_min, q2_max, k2_min, k2_max));
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToThreeLeptonsNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToGammaOffShell>::option_specification(),
        { "lprime"_ok, { "e", "mu", "tau" }, "mu" },
        { "l"_ok, { "e", "mu", "tau" }, "e" },
    };

    BToThreeLeptonsNeutrino::BToThreeLeptonsNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToThreeLeptonsNeutrino>(new Implementation<BToThreeLeptonsNeutrino>(parameters, options, *this))
    {
    }

    BToThreeLeptonsNeutrino::~BToThreeLeptonsNeutrino()
    {
    }

    double
    BToThreeLeptonsNeutrino::double_differential_branching_ratio(const double & q2,
        const double & k2) const
    {
        return _imp->double_differential_branching_ratio(q2, k2);
    }

    double
    BToThreeLeptonsNeutrino::quintuple_differential_branching_ratio(const double & q2,
        const double & k2, const double & z_gamma, const double & z_w, const double & phi) const
    {
        return _imp->quintuple_differential_branching_ratio(q2, k2, z_gamma, z_w, phi);
    }

    double
    BToThreeLeptonsNeutrino::integrated_branching_ratio(const double & q2_min,
        const double & q2_max, const double & k2_min, const double & k2_max) const
    {
        return _imp->integrated_branching_ratio(q2_min, q2_max, k2_min, k2_max);
    }

    double
    BToThreeLeptonsNeutrino::double_differential_forward_backward_asymmetry(const double & q2,
        const double & k2) const
    {
        return _imp->double_differential_forward_backward_asymmetry(q2, k2);
    }

    double
    BToThreeLeptonsNeutrino::integrated_forward_backward_asymmetry(const double & q2_min,
        const double & q2_max, const double & k2_min, const double & k2_max) const
    {
        return _imp->integrated_forward_backward_asymmetry(q2_min, q2_max, k2_min, k2_max);
    }

    const std::string
    BToThreeLeptonsNeutrino::description = "\
The decayB^- -> l^- nubar l_prime^+ l_prime^-, where l is either e, mu or tau and l_prime is either e or mu";

    const std::string
    BToThreeLeptonsNeutrino::kinematics_description_q2 = "\
The invariant mass of the l_prime pair in GeV^2.";

    const std::string
    BToThreeLeptonsNeutrino::kinematics_description_k2 = "\
The invariant mass of the l^- nubar pair in GeV^2.";

    const std::string
    BToThreeLeptonsNeutrino::kinematics_description_z_gamma = "\
The cosine of the angle between l_prime^- and the negative photon direction of flight.";

    const std::string
    BToThreeLeptonsNeutrino::kinematics_description_z_w = "\
The cosine of the angle between charged l^- and the negative W boson direction of flight.";

    const std::string
    BToThreeLeptonsNeutrino::kinematics_description_phi = "\
The angle between the q2 and the k2 plane.";

    const std::set<ReferenceName>
    BToThreeLeptonsNeutrino::references
    {
        "KKvDZ:2022A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToThreeLeptonsNeutrino::begin_options()
    {
        return Implementation<BToThreeLeptonsNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToThreeLeptonsNeutrino::end_options()
    {
        return Implementation<BToThreeLeptonsNeutrino>::options.cend();
    }
}
