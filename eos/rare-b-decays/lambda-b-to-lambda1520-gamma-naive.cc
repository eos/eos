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

#include <eos/maths/power-of.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma-naive.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;

    LambdaBToLambda1520GammaAmplitudes<tag::Naive>::LambdaBToLambda1520GammaAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o)
    {
    }

    inline double
    LambdaBToLambda1520GammaAmplitudes<tag::Naive>::mu_f() const
    {
        return 1.5;
    }

    inline double
    LambdaBToLambda1520GammaAmplitudes<tag::Naive>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    LambdaBToLambda1520Gamma::Amplitudes
    LambdaBToLambda1520GammaAmplitudes<tag::Naive>::amplitudes() const
    {
        LambdaBToLambda1520Gamma::Amplitudes result;

        // Import WC with a fake lepton flavor
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), LeptonFlavor::muon, cp_conjugate);

        const double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
        const double m_b_msbar = model->m_b_msbar(mu());

        const complex<double> c7eff = ShortDistanceLowRecoil::c7eff(0, mu_f(), alpha_s_mu, m_b_PS(), false, wc);

        const complex<double>
            wilson7_plus          = (c7eff + wc.c7prime()),
            wilson7_minus         = (c7eff - wc.c7prime());

        double lambda_t2 = std::norm(model->ckm_tb() * conj(model->ckm_ts()));

        const double norm = g_fermi() * sqrt(
                1.0 / 3.0 / 128.0 / power_of<4>(M_PI) / power_of<3>(m_Lb)
                * lambda_t2 * alpha_e() * (power_of<2>(m_Lb) - power_of<2>(m_Lstar))
            );

        result.a_perp12 = norm * m_b_msbar * wilson7_plus  * power_of<2>(m_Lb + m_Lstar) * form_factors->f_perp12_t(0);
        result.a_para12 = norm * m_b_msbar * wilson7_minus * power_of<2>(m_Lb - m_Lstar) * form_factors->f_perp12_t5(0);
        result.a_perp32 = norm * m_b_msbar * wilson7_plus  * power_of<2>(m_Lb + m_Lstar) * form_factors->f_perp32_t(0);
        result.a_para32 = norm * m_b_msbar * wilson7_minus * power_of<2>(m_Lb + m_Lstar) * form_factors->f_perp32_t5(0);

        return result;
    }
}
