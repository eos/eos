/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Méril Reboud
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

#include <eos/rare-b-decays/b-to-kstar-ll-gvdv2020.hh>
#include <eos/utils/kinematic.hh>

using namespace std;

namespace eos
{
    using namespace std::placeholders;
    using std::norm;
    using std::sqrt;

    BToKstarDileptonAmplitudes<tag::GvDV2020>::BToKstarDileptonAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        tau(p["life_time::B_" + o.get("q", "d")], *this),
        q(o, "q", { "d", "u" }, "d"),
        opt_nonlocal_formfactor(o, "nonlocal-formfactor", { "GvDV2020", "GRvDV2021", "naive" }, "GvDV2020"),
        nonlocal_formfactor(NonlocalFormFactor<nff::PToV>::make("B->K^*::" + opt_nonlocal_formfactor.value(), p, o))
    {
    }

    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::GvDV2020>::amplitudes(const double & s) const
    {
        BToKstarDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        // classic form factors
        const double
                ff_V  = form_factors->v(s),
                ff_A0 = form_factors->a_0(s),
                ff_A1 = form_factors->a_1(s),
                ff_A2 = form_factors->a_2(s),
                ff_T1 = form_factors->t_1(s),
                ff_T2 = form_factors->t_2(s),
                ff_T3 = form_factors->t_3(s);

        // kinematics
        const double
                sqrt_s      = std::sqrt(s),
                m_B         = this->m_B(),
                m_B2        = power_of<2>(m_B),
                m_V         = this->m_Kstar(),
                m_V2        = power_of<2>(m_V),
                lambda      = eos::lambda(m_B2, m_V2, s),
                sqrt_lambda = std::sqrt(lambda);

        // vectorial form factors, cf. [GvDV2020], eq. (A.11)
        const double
                calF_perp = sqrt(2.0) * sqrt_lambda / (m_B * (m_B + m_V)) * ff_V,
                calF_para = sqrt(2.0) * (m_B + m_V) / m_B * ff_A1,
                calF_long = ((m_B2 - m_V2 - s) * power_of<2>(m_B + m_V) * ff_A1 - lambda * ff_A2)
                          / (2.0 * m_V * m_B2 * (m_B + m_V)),
                calF_time = ff_A0;

        // tensorial form factors, cf. [GvDV2020], eq. (A.11)
        const double
                calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1,
                calF_T_para = sqrt(2.0) * (m_B2 - m_V2) / m_B2 * ff_T2,
                calF_T_long = s / (2.0 * power_of<3>(m_B) * m_V) *
                            ((m_B2 + 3.0 * m_V2 - s) * ff_T2 - lambda / (m_B2 - m_V2) * ff_T3);

        const complex<double>
                calH_perp = nonlocal_formfactor->H_perp(s),
                calH_para = nonlocal_formfactor->H_para(s),
                calH_long = nonlocal_formfactor->H_long(s);

        // Wilson coefficients
        const complex<double>
                c910_m_r = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime()),
                c910_m_l = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime()),
                c910_p_r = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime()),
                c910_p_l = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime()),
                c7_m = (wc.c7() - wc.c7prime()),
                c7_p = (wc.c7() + wc.c7prime());

        // quark masses
        const double
                m_b_msbar = model->m_b_msbar(mu()),
                m_s_msbar = model->m_s_msbar(mu());

        // normalization constant, cf. KM2005A (3.7)
        const double calN = g_fermi() * alpha_e * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * sqrt(s * beta_l(s) * sqrt_lambda / (3.0 * 1024 * power_of<5>(M_PI) * m_B));

        // vector amplitudes, cf. KM2005A (3.2) - (3.4)
        result.a_long_right = -calN * m_B2 / sqrt_s * (c910_m_r * calF_long + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_long - 16.0 * power_of<2>(M_PI) * m_B * calH_long));
        result.a_long_left  = -calN * m_B2 / sqrt_s * (c910_m_l * calF_long + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_long - 16.0 * power_of<2>(M_PI) * m_B * calH_long));

        result.a_para_right = -calN * m_B * (c910_m_r * calF_para + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_para - 16.0 * power_of<2>(M_PI) * m_B * calH_para));
        result.a_para_left  = -calN * m_B * (c910_m_l * calF_para + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_para - 16.0 * power_of<2>(M_PI) * m_B * calH_para));

        result.a_perp_right = +calN * m_B * (c910_p_r * calF_perp + 2.0 * m_B / s * ((m_b_msbar + m_s_msbar) * c7_p * calF_T_perp - 16.0 * power_of<2>(M_PI) * m_B * calH_perp));
        result.a_perp_left  = +calN * m_B * (c910_p_l * calF_perp + 2.0 * m_B / s * ((m_b_msbar + m_s_msbar) * c7_p * calF_T_perp - 16.0 * power_of<2>(M_PI) * m_B * calH_perp));

        // scalar amplitude, cf. KM2005A (3.5)
        result.a_time = calN * 2.0 * sqrt_lambda / sqrt_s * (wc.c10() - wc.c10prime()) * calF_time;

        return result;
    }
}
