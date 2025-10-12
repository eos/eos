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
#include <eos/rare-b-decays/b-to-kstar-ll-gp2004.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/nonlocal-form-factors/hard-scattering.hh>
#include <eos/nonlocal-form-factors/long-distance.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{
    using namespace std::literals::string_literals;

    BToKstarDileptonAmplitudes<tag::GP2004>::BToKstarDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        hbar(p["QM::hbar"], *this),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c_MSbar(p["mass::c"], *this),
        m_s(p["mass::s(2GeV)"], *this),
        opt_use_simple_sl(o, options, "simple-sl"_ok),
        use_simple_sl(opt_use_simple_sl.value()),
        lambda_long(p["B->Vll::Lambda" + std::string(use_simple_sl ? "" : "_0") + "@LowRecoil"], *this),
        lambda_par(p["B->Vll::Lambda" + std::string(use_simple_sl ? "" : "_pa") + "@LowRecoil"], *this),
        lambda_perp(p["B->Vll::Lambda" + std::string(use_simple_sl ? "" : "_pp") + "@LowRecoil"], *this),
        sl_phase_long(p["B->Vll::sl_phase" + std::string(use_simple_sl ? "" : "_0") + "@LowRecoil"], *this),
        sl_phase_par(p["B->Vll::sl_phase" + std::string(use_simple_sl ? "" : "_pa") + "@LowRecoil"], *this),
        sl_phase_perp(p["B->Vll::sl_phase" + std::string(use_simple_sl ? "" : "_pp") + "@LowRecoil"], *this),
        opt_ccbar_resonance(o, options, "ccbar-resonance"_ok),
        opt_use_nlo(o, options, "nlo"_ok),
        ccbar_resonance(opt_ccbar_resonance.value()),
        use_nlo(opt_use_nlo.value())
    {
        Context ctx("When constructing B->K^*ll GP2004 amplitudes");
    }

    BToKstarDileptonAmplitudes<tag::GP2004>::~BToKstarDileptonAmplitudes()
    {
    }

    const std::vector<OptionSpecification>
    BToKstarDileptonAmplitudes<tag::GP2004>::options
    {
        { "ccbar-resonance"_ok, { "true"s, "false"s },  "false"s },
        { "nlo"_ok,             { "true"s, "false"s },  "true"s  },
        { "simple-sl"_ok,       { "true"s, "false"s },  "false"s }
    };

    // cf. [GP2004], Eq. (56)
    complex<double>
    BToKstarDileptonAmplitudes<tag::GP2004>::c7eff(const WilsonCoefficients<BToS> & wc, const double & s) const
    {
        return ShortDistanceLowRecoil::c7eff(s, mu(), model->alpha_s(mu), m_b_PS(), use_nlo, wc);
    }

    // cf. [GP2004], Eq. (55), p. 10
    complex<double>
    BToKstarDileptonAmplitudes<tag::GP2004>::c9eff(const WilsonCoefficients<BToS> & wc, const double & s) const
    {
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
        {
            lambda_hat_u = conj(lambda_hat_u);
        }

        return ShortDistanceLowRecoil::c9eff(s, mu(), model->alpha_s(mu), m_b_PS(), model->m_c_msbar(mu), use_nlo, ccbar_resonance, lambda_hat_u, wc);
    }

    // We use the PS mass except for kappa
    double
    BToKstarDileptonAmplitudes<tag::GP2004>::m_b_PS() const
    {
        // Actually use m_b_PS at mu_PS = 2.0 GeV
        return model->m_b_ps(2.0);
    }

    double
    BToKstarDileptonAmplitudes<tag::GP2004>::kappa() const
    {
        // cf. [BHvD2010], Eq. (3.8), p. 8
        // Use m_b_MSbar(m_b_MSbar) instead m_b_MSbar(mu), as we want kappa up to NLO only.
        return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / m_b_MSbar));
    }

    double
    BToKstarDileptonAmplitudes<tag::GP2004>::norm(const double & s) const
    {
        double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

        return std::sqrt(power_of<2>(g_fermi() * alpha_e()) / 3.0 / 1024 / power_of<5>(M_PI) / m_B
                * lambda_t * lambda_t * s_hat(s) * beta_l(s)
                * std::sqrt(eos::lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
    }

    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::GP2004>::amplitudes(const double & s) const
    {
        // compute J_i, [BHvD2010], p. 26, Eqs. (A1)-(A11)
        // TODO: possibly optimize the calculation
        BToKstarDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const double m_B2 = m_B * m_B, m_Kstar2 = m_Kstar * m_Kstar, m2_diff = m_B2 - m_Kstar2;
        const double m_Kstarhat = m_Kstar / m_B;
        const double m_Kstarhat2 = power_of<2>(m_Kstarhat);
        const double s_hat = s / m_B / m_B;
        const double a_1 = form_factors->a_1(s), a_2 = form_factors->a_2(s);
        const double alpha_s = model->alpha_s(mu());
        const double norm_s = this->norm(s);
        const double lam = this->lambda(s);
        const double sqrt_lam = std::sqrt(lam);
        const double sqrt_s = std::sqrt(s);

        const complex<double> subleading_perp = 0.5 / m_B * alpha_s * std::polar(lambda_perp(), sl_phase_perp());
        const complex<double> subleading_par  = 0.5 / m_B * alpha_s * std::polar(lambda_par(), sl_phase_par());
        const complex<double> subleading_long = 0.5 / m_B * alpha_s * std::polar(lambda_long(), sl_phase_long());

        const complex<double> c_9eff = c9eff(wc, s);
        const complex<double> c_7eff = c7eff(wc, s);
        const complex<double> c910_plus_left   = (c_9eff + wc.c9prime()) - (wc.c10() + wc.c10prime());
        const complex<double> c910_plus_right  = (c_9eff + wc.c9prime()) + (wc.c10() + wc.c10prime());
        const complex<double> c910_minus_left  = (c_9eff - wc.c9prime()) - (wc.c10() - wc.c10prime());
        const complex<double> c910_minus_right = (c_9eff - wc.c9prime()) + (wc.c10() - wc.c10prime());
        const complex<double> c7_plus  = kappa() * (c_7eff + wc.c7prime()) * (2.0 * m_B / s);
        const complex<double> c7_minus = kappa() * (c_7eff - wc.c7prime()) * (2.0 * m_B / s);

        // longitudinal
        complex<double> prefactor_long = complex<double>(-1.0, 0.0) * m_B()
            / (2.0 * m_Kstarhat * (1.0 + m_Kstarhat) * std::sqrt(s_hat));
        complex<double> wilson_long1_right = c910_minus_right + c7_minus * (m_b_MSbar() - m_s() - lambda_par()) + subleading_par;
        complex<double> wilson_long1_left  = c910_minus_left  + c7_minus * (m_b_MSbar() - m_s() - lambda_par()) + subleading_par;
        complex<double> wilson_long2_right = c910_minus_right + c7_minus * (m_b_MSbar() - m_s() - lambda_long()) - subleading_long;
        complex<double> wilson_long2_left  = c910_minus_left  + c7_minus * (m_b_MSbar() - m_s() - lambda_long()) - subleading_long;

        double formfactor_long1 = (1.0 - m_Kstarhat2 - s_hat) * power_of<2>(1.0 + m_Kstarhat) * a_1;
        double formfactor_long2 = -eos::lambda(1.0, m_Kstarhat2, s_hat) * a_2;
        // cf. [BHvD2010], Eq. (3.15), p. 10
        result.a_long_right = norm_s * prefactor_long * (wilson_long1_right * formfactor_long1 + wilson_long2_right * formfactor_long2);
        result.a_long_left  = norm_s * prefactor_long * (wilson_long1_left  * formfactor_long1 + wilson_long2_left  * formfactor_long2);

        // perpendicular
        complex<double> prefactor_perp = complex<double>(1.0, 0.0) * m_B();
        complex<double> wilson_perp_right = c910_plus_right + c7_plus * (m_b_MSbar() + m_s() + lambda_perp()) - subleading_perp;
        complex<double> wilson_perp_left  = c910_plus_left  + c7_plus * (m_b_MSbar() + m_s() + lambda_perp()) - subleading_perp;

        double formfactor_perp = std::sqrt(2.0 * eos::lambda(1.0, m_Kstarhat2, s_hat)) / (1.0 + m_Kstarhat) * form_factors->v(s);
        // cf. [BHvD2010], Eq. (3.13), p. 10
        result.a_perp_right = norm_s * prefactor_perp * wilson_perp_right * formfactor_perp;
        result.a_perp_left  = norm_s * prefactor_perp * wilson_perp_left  * formfactor_perp;

        // parallel
        complex<double> prefactor_par = complex<double>(-1.0, 0.0) * m_B();
        complex<double> wilson_par_right = c910_minus_right + c7_minus * (m_b_MSbar() - m_s() - lambda_par()) + subleading_par;
        complex<double> wilson_par_left  = c910_minus_left  + c7_minus * (m_b_MSbar() - m_s() - lambda_par()) + subleading_par;
        double formfactor_par = std::sqrt(2) * (1.0 + m_Kstarhat) * a_1;
        // cf. [BHvD2010], Eq. (3.14), p. 10
        result.a_para_right = norm_s * prefactor_par * wilson_par_right * formfactor_par;
        result.a_para_left  = norm_s * prefactor_par * wilson_par_left  * formfactor_par;

        // timelike
        result.a_time = norm_s * sqrt_lam / sqrt_s
            * (2.0 * (wc.c10() - wc.c10prime()) + s / m_l / (m_b_MSbar + m_s()) * (wc.cP() - wc.cPprime()))
            * form_factors->a_0(s);

        // scalar amplitude
        result.a_scal = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s()) * form_factors->a_0(s);

        // tensor amplitudes [BHvD2012]  eqs. (B18 - B20)
        // no form factor relations used
        const double ff_T1  = form_factors->t_1(s);
        const double ff_T2  = form_factors->t_2(s);
        const double ff_T3  = form_factors->t_3(s);

        const double kin_tensor_1 = norm_s / m_Kstar * ((m_B2 + 3.0 * m_Kstar2 - s) * ff_T2 - lam / m2_diff * ff_T3);
        const double kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1;
        const double kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

        // correct the sign of C_T5 from [BHvD2012] (arXiv v4) because of inconsistent use of
        // gamma5 <-> Levi-Civita
        static const double sign = -1;

        result.a_para_perp = kin_tensor_1 * wc.cT();
        result.a_time_long  = kin_tensor_1 * sign * wc.cT5();

        result.a_time_perp = kin_tensor_2 * wc.cT();
        result.a_long_perp = kin_tensor_2 * sign * wc.cT5();

        result.a_time_para = kin_tensor_3 * sign * wc.cT5();
        result.a_long_para = kin_tensor_3 * wc.cT();

        return result;
    }
}
