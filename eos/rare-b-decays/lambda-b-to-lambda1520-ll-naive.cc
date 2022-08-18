/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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

#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll-naive.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;

    LambdaBToLambda1520DileptonAmplitudes<tag::Naive>::LambdaBToLambda1520DileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        m_c(p["mass::c"], *this)
    {
    }

    LambdaBToLambda1520DileptonAmplitudes<tag::Naive>::~LambdaBToLambda1520DileptonAmplitudes()
    {
    }

    double
    LambdaBToLambda1520DileptonAmplitudes<tag::Naive>::norm(const double & s) const
    {
        // cf. [DN:2019A], eqs. (3.18 - 3.20)
        double lambda_t2 = std::norm(model->ckm_tb() * conj(model->ckm_ts()));

        return g_fermi() * alpha_e() * std::sqrt(
                  1.0 / 3.0 / 2048.0 / power_of<5>(M_PI) / power_of<3>(m_Lb)
                  * lambda_t2 * s * std::sqrt(lambda(s))
               );
    }

    double
    LambdaBToLambda1520DileptonAmplitudes<tag::Naive>::mu_f() const
    {
        return 1.5;
    }

    double
    LambdaBToLambda1520DileptonAmplitudes<tag::Naive>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    /* Amplitudes */
    // cf. [DN:2019A], eqs. (3.18 - 3.20)
    LambdaBToLambda1520Dilepton::Amplitudes
    LambdaBToLambda1520DileptonAmplitudes<tag::Naive>::amplitudes(const double & s) const
    {
        LambdaBToLambda1520Dilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const double
            norm_s = norm(s),
            sqrt_s = std::sqrt(s),
            s_minus = power_of<2>(m_Lb - m_Lstar) - s,
            s_plus = power_of<2>(m_Lb + m_Lstar) - s;

        const double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
        const double m_b_msbar = model->m_b_msbar(mu());

        const complex<double>
            c9eff = ShortDistanceLowRecoil::c9eff(s, mu_f(), alpha_s_mu, m_b_PS(), m_c, false, false, 0.0, wc),
            c7eff = ShortDistanceLowRecoil::c7eff(s, mu_f(), alpha_s_mu, m_b_PS(), false, wc);

        const complex<double>
            wilson910_minus_right = (c9eff - wc.c9prime()) + (wc.c10() - wc.c10prime()),
            wilson910_minus_left  = (c9eff - wc.c9prime()) - (wc.c10() - wc.c10prime()),
            wilson910_plus_right  = (c9eff + wc.c9prime()) + (wc.c10() + wc.c10prime()),
            wilson910_plus_left   = (c9eff + wc.c9prime()) - (wc.c10() + wc.c10prime()),
            wilson7_plus          = (c7eff + wc.c7prime()),
            wilson7_minus         = (c7eff - wc.c7prime());

        const double
            H0Vp12 = - (m_Lb + m_Lstar) / sqrt_s * sqrt(s_plus / 6.0) * form_factors->f_long12_v(s),
            HplusVm12 = - sqrt(s_plus / 3.0) * form_factors->f_perp12_v(s),
            HplusVm32 = sqrt(s_plus) * form_factors->f_perp32_v(s),
            H0Ap12 = - (m_Lb - m_Lstar) / sqrt_s * sqrt(s_minus / 6.0) * form_factors->f_long12_a(s),
            HplusAm12 = sqrt(s_minus / 3.0) * form_factors->f_perp12_a(s),
            HplusAm32 = sqrt(s_minus) * form_factors->f_perp32_a(s);

        const double
            H0Tp12 = sqrt_s * sqrt(s_plus / 6.0) * form_factors->f_long12_t(s),
            HplusTm12 = (m_Lb + m_Lstar) * sqrt(s_plus / 3.0) * form_factors->f_perp12_t(s),
            HplusTm32 = - (m_Lb + m_Lstar) * sqrt(s_plus) * form_factors->f_perp32_t(s),
            H0T5p12 = - sqrt_s * sqrt(s_minus / 6.0) * form_factors->f_long12_t5(s),
            HplusT5m12 = (m_Lb - m_Lstar) * sqrt(s_minus / 3.0) * form_factors->f_perp12_t5(s),
            HplusT5m32 = (m_Lb - m_Lstar) * sqrt(s_minus) * form_factors->f_perp32_t5(s);

        result.b_perp1_right =   sqrt(2.0) * norm_s * (wilson910_plus_right  * HplusVm32 - 2.0 * m_b_msbar / s * wilson7_plus  * HplusTm32);
        result.b_perp1_left =    sqrt(2.0) * norm_s * (wilson910_plus_left   * HplusVm32 - 2.0 * m_b_msbar / s * wilson7_plus  * HplusTm32);
        result.b_para1_right = - sqrt(2.0) * norm_s * (wilson910_minus_right * HplusAm32 + 2.0 * m_b_msbar / s * wilson7_minus * HplusT5m32);
        result.b_para1_left =  - sqrt(2.0) * norm_s * (wilson910_minus_left  * HplusAm32 + 2.0 * m_b_msbar / s * wilson7_minus * HplusT5m32);
        result.a_perp1_right =   sqrt(2.0) * norm_s * (wilson910_plus_right  * HplusVm12 - 2.0 * m_b_msbar / s * wilson7_plus  * HplusTm12);
        result.a_perp1_left =    sqrt(2.0) * norm_s * (wilson910_plus_left   * HplusVm12 - 2.0 * m_b_msbar / s * wilson7_plus  * HplusTm12);
        result.a_para1_right = - sqrt(2.0) * norm_s * (wilson910_minus_right * HplusAm12 + 2.0 * m_b_msbar / s * wilson7_minus * HplusT5m12);
        result.a_para1_left =  - sqrt(2.0) * norm_s * (wilson910_minus_left  * HplusAm12 + 2.0 * m_b_msbar / s * wilson7_minus * HplusT5m12);
        result.a_perp0_right =   sqrt(2.0) * norm_s * (wilson910_plus_right  * H0Vp12    - 2.0 * m_b_msbar / s * wilson7_plus  * H0Tp12);
        result.a_perp0_left =    sqrt(2.0) * norm_s * (wilson910_plus_left   * H0Vp12    - 2.0 * m_b_msbar / s * wilson7_plus  * H0Tp12);
        result.a_para0_right = - sqrt(2.0) * norm_s * (wilson910_minus_right * H0Ap12    + 2.0 * m_b_msbar / s * wilson7_minus * H0T5p12);
        result.a_para0_left =  - sqrt(2.0) * norm_s * (wilson910_minus_left  * H0Ap12    + 2.0 * m_b_msbar / s * wilson7_minus * H0T5p12);

        return result;
    }
}
