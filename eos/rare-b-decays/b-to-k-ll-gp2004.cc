/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2010-2011 Christian Wacker
 * Copyright (c) 2014      Frederik Beaujean
 * Copyright (c) 2014      Christoph Bobeth
 * Copyright (c) 2021      Méril Reboud
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
#include <eos/rare-b-decays/b-to-k-ll-gp2004.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/utils/destringify.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::literals::string_literals;

    BToKDileptonAmplitudes<tag::GP2004>::BToKDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        hbar(p["QM::hbar"], *this),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c_MSbar(p["mass::c"], *this),
        m_s(p["mass::s(2GeV)"], *this),
        lambda_psd(p["B->Pll::Lambda_pseudo@LargeRecoil"], *this),
        sl_phase_psd(p["B->Pll::sl_phase_pseudo@LargeRecoil"], *this),
        opt_ccbar_resonance(o, options, "ccbar-resonance"_ok),
        opt_use_nlo(o, options, "nlo"_ok),
        ccbar_resonance(opt_ccbar_resonance.value()),
        use_nlo(opt_use_nlo.value())
    {
        Context ctx("When constructing B->Kll GP2004 amplitudes");
    }

    BToKDileptonAmplitudes<tag::GP2004>::~BToKDileptonAmplitudes()
    {
    }

    const std::vector<OptionSpecification>
    BToKDileptonAmplitudes<tag::GP2004>::options
    {
        { "ccbar-resonance"_ok, { "true"s, "false"s },  "false"s },
        { "nlo"_ok, { "true"s, "false"s },  "true"s },
    };

    // We use the PS mass except for kappa
    double
    BToKDileptonAmplitudes<tag::GP2004>::m_b_PS() const
    {
        // Actually use m_b_PS at mu_PS = 2.0 GeV
        return model->m_b_ps(2.0);
    }

    // cf. [GP:2004A], Eq. (56)
    complex<double>
    BToKDileptonAmplitudes<tag::GP2004>::c7eff(const WilsonCoefficients<BToS> & wc, const double & s) const
    {
        return ShortDistanceLowRecoil::c7eff(s, mu(), model->alpha_s(mu), m_b_PS(), use_nlo, wc);
    }

    // cf. [GP:2004A], Eq. (55), p. 10
    complex<double>
    BToKDileptonAmplitudes<tag::GP2004>::c9eff(const WilsonCoefficients<BToS> & wc, const double & s) const
    {
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
        {
            lambda_hat_u = conj(lambda_hat_u);
        }

        return ShortDistanceLowRecoil::c9eff(s, mu(), model->alpha_s(mu), m_b_PS(), model->m_c_msbar(mu), use_nlo, ccbar_resonance, lambda_hat_u, wc);
    }

    double
    BToKDileptonAmplitudes<tag::GP2004>::kappa() const
    {
        // cf. [BHvD:2010A], Eq. (3.8), p. 8
        // Use m_b_MSbar(m_b_MSbar) instead m_b_MSbar(mu), as we want kappa up to NLO only.
        return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / m_b_MSbar));
    }

    // // this is rho_1^+
    // double
    // BToKDileptonAmplitudes<tag::GP2004>::rho_1(const double & s) const
    // {
    //     WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

    //     double alpha_s = model->alpha_s(mu());

    //     return std::norm(kappa() * (2.0 * (m_b_MSbar + lambda_psd()) * m_B() / s) * (c7eff(wc, s) + wc.c7prime())
    //             + 0.5 * alpha_s / m_B * std::polar(lambda_psd(), sl_phase_psd()) + (c9eff(wc, s) + wc.c9prime()))
    //             + std::norm(wc.c10() + wc.c10prime());
    // }

    /* Amplitudes */
    BToKDilepton::Amplitudes
    BToKDileptonAmplitudes<tag::GP2004>::amplitudes(const double & s) const
    {
        BToKDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        // cf. [BF:2001A] Eq. (22 + TODO: 31)
        // cf. [BF:2001A] Eq. (22 + TODO: 30)
        double f_t_over_f_p = form_factors->f_t(s) / form_factors->f_p(s);
        double f_0_over_f_p = form_factors->f_0(s) / form_factors->f_p(s);

        double F_Tkin = f_t_over_f_p * 2.0 * std::sqrt(lambda(s)) * beta_l(s) / (m_B() + m_K());
        double F_Skin = f_0_over_f_p * 0.5 * (power_of<2>(m_B()) - power_of<2>(m_K())) / (m_b_MSbar - m_s);

        // cf. [BHP:2007A], Eq. (3.2), p. 3 and 4
        result.F_A  = wc.c10() + wc.c10prime();
        result.F_T  = F_Tkin * wc.cT();
        result.F_T5 = F_Tkin * wc.cT5();
        result.F_S  = F_Skin * (wc.cS() + wc.cSprime());
        result.F_P  = F_Skin * (wc.cP() + wc.cPprime()) + m_l() * (wc.c10() + wc.c10prime()) *
                      ((m_B() * m_B() - m_K() * m_K()) / s * (f_0_over_f_p - 1.0) - 1.0);
        result.F_V  = c9eff(wc, s) + wc.c9prime()
                      + kappa() * (2.0 * (m_b_MSbar + lambda_psd()) * m_B() / s) * (c7eff(wc, s) + wc.c7prime())
                      + 0.5 * model->alpha_s(mu) / m_B * std::polar(lambda_psd(), sl_phase_psd())
                      + 8.0 * m_l / (m_B() + m_K()) * f_t_over_f_p * wc.cT();

        return result;
    }
}
