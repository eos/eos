/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2025 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/rare-b-decays/b-to-k-ll-gvdv2020.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::literals::string_literals;
    using namespace std::placeholders;

    BToKDileptonAmplitudes<tag::GvDV2020>::BToKDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        f_B(p["decay-constant::B_" + o.get("q"_ok, "d")], *this),
        f_K(p["decay-constant::K_" + o.get("q"_ok, "d")], *this),
        lambda_B_p_inv(p["B::1/lambda_B_p"], *this),
        q(o, options, "q"_ok),
        opt_nonlocal_formfactor(o, options, "nonlocal-formfactor"_ok),
        nonlocal_formfactor(NonlocalFormFactor<PToP>::make("B->K::" + opt_nonlocal_formfactor.value(), p, o))
    {
        Context ctx("When constructing B->Kll GvDV2020 amplitudes");
    }

    BToKDileptonAmplitudes<tag::GvDV2020>::~BToKDileptonAmplitudes()
    {
    }

    const std::vector<OptionSpecification>
    BToKDileptonAmplitudes<tag::GvDV2020>::options
    {
        { "q"_ok, { "d"s, "u"s },  "d"s },
        { "nonlocal-formfactor"_ok, { "GvDV2020"s, "GRvDV2022order5"s, "GRvDV2022order6"s }, "GvDV2020"s }
    };

    BToKDilepton::DipoleFormFactors
    BToKDileptonAmplitudes<tag::GvDV2020>::dipole_form_factors(const double & s, const WilsonCoefficients<BToS> & wc) const
    {
        // charges of down- and up-type quarks
        static const double e_d = -1.0 / 3.0;
        static const double e_u = +2.0 / 3.0;

        // spectator contributions
        double delta_qu = 0.0, e_q = e_d;
        if (q.value() == QuarkFlavor::up)
        {
            delta_qu = 1.0;
            e_q = e_u;
        }

        // kinematics
        double m_b_PS = this->m_b_PS();

        // couplings
        double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
        double a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI;
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
            lambda_hat_u = std::conj(lambda_hat_u);


        // inverse of the "negative" moment of the B meson LCDA
        // cf. [BFS:2001], Eq. (54), p. 15
        double omega_0 = 1.0 / this->lambda_B_p_inv;
        complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

        /* Y(s) for the up and the top sector */
        // cf. [BFS:2001], Eq. (10), p. 4
//        complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
        complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

        // Use b pole mass according to [BFS:2001], Sec. 3.1, paragraph Quark Masses,
        // then replace b pole mass by the PS mass.
        complex<double> Y_top = //Y_top_c * CharmLoops::h(mu, s, m_c_pole);
                + Y_top_b * CharmLoops::h(mu, s, m_b_PS);
                + Y_top_0 * CharmLoops::h(mu, s);
                + Y_top_;
        // cf. [BFS:2004], Eq. (43), p. 24
        complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (//CharmLoops::h(mu, s, m_c_pole)
                - CharmLoops::h(mu, s));

        /* Effective wilson coefficients */
        complex<double> c8eff = ShortDistanceLowRecoil::c8eff(wc); // LO C8eff

        /* top sector */
        // cf. [BHP:2007], Eq. (B.2) and [BFS:2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
        complex<double> C0_top_psd = 1.0 * (m_B / (2.0 * m_b_PS) * Y_top); // c7eff + wc.c7prime() +
        // cf. [BHP:2007], Eq. (B.2) and [BFS:2004], Eq. (45), p. 24
        // the correct sign in front of C_7^eff is plus, as one can see by
        // comparison with [BF:2001], Eq. (63)
//        complex<double> C1f_top_psd = 1.0 * (c7eff + wc.c7prime()) * (8.0 * std::log(m_b_PS / mu) + 2.0 * L - 4.0 * (1.0 - mu_f() / m_b_PS));
        // cf. [BHP:2007], Eq. (B.2) and [BFS:2001], Eqs. (38), p. 9
        complex<double> C1nf_top_psd = -(+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive_Qsb, s)
                + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * memoise(CharmLoops::F19_massive_Qsb, s)
                    + wc.c2() * memoise(CharmLoops::F29_massive_Qsb, s)
                    + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

        /* parallel, up sector */
        // cf. [BHP:2007], Eq. (B.2) and [BFS:2004], comment before Eq. (43), p. 24
        complex<double> C0_up_psd = 1.0 * m_B / (2.0 * m_b_PS) * Y_up;
        // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
        // cf. [BFS:2004], last paragraph in Sec A.1, p. 24
        // [BFS:2004], [S:2004] have a different sign convention for F{12}{79}_massless than we!
        // Use here FF_massive - FF_massless because FF_massless is defined with an extra '-'
        // compared to [S:2004]
        complex<double> C1nf_up_psd = -(+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive_Qsb, s) - CharmLoops::F27_massless(mu, s, m_b_PS))
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * (memoise(CharmLoops::F19_massive_Qsb, s) - CharmLoops::F19_massless(mu, s, m_b_PS))
                    + wc.c2() * (memoise(CharmLoops::F29_massive_Qsb, s) - CharmLoops::F29_massless(mu, s, m_b_PS))));

        // compute the factorizing contributions
        complex<double> C_psd = C0_top_psd + lambda_hat_u * C0_up_psd
                + a_mu * (C1nf_top_psd + lambda_hat_u * C1nf_up_psd);

        /* parallel, top sector */
        // T0_top_par_p = 0, cf. [BFS:2001], Eq. (17), p. 6
        // cf. [BFS:2004], Eqs. (46)-(47), p. 25 without the \omega term.
        complex<double> T0_top_psd_m = +e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;

        /* parallel, up sector */
        // all T1f_up vanish, cf. [BFS:2004], sentence below Eq. (49), p. 25
        // cf. [BFS:2004], Eqs. (46),(48), p. 25 without the \omega term
        complex<double> T0_up_psd_m = -e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;


        // Compute the nonfactorizing contributions
        complex<double> T_psd = T0_top_psd_m + lambda_hat_u * T0_up_psd_m;

        // Subleading weak annihilation and hard spectator interaction contributions have only been
        // computed for calT_perp, not for calT_par ~ calT_psd.

        // cf. [BFS:2001], Eq. (15), and [BHP:2008], Eq. (C.4)
        BToKDilepton::DipoleFormFactors result;
        result.calT = xi_pseudo(s) * C_psd + power_of<2>(M_PI) / 3.0 * (f_B * f_K) / m_B * T_psd;

        return result;
    }

    double
    BToKDileptonAmplitudes<tag::GvDV2020>::xi_pseudo(const double & s) const
    {
        // cf. [BF:2001], Eq. (22)
        return form_factors->f_p(s);
    }

    double
    BToKDileptonAmplitudes<tag::GvDV2020>::mu_f() const
    {
        return 1.5;
    }

    double
    BToKDileptonAmplitudes<tag::GvDV2020>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    /* Amplitudes */
    BToKDilepton::Amplitudes
    BToKDileptonAmplitudes<tag::GvDV2020>::amplitudes(const double & s) const
    {
        BToKDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        auto dff = dipole_form_factors(s, wc);

        const double m_B2 = m_B * m_B, m_K2 = m_K * m_K;

        // cf. [GvDV:2020] Eq. (A.5)
        const double
            calF_plus   = form_factors->f_p(s),
            calF_time   = form_factors->f_0(s),
            calF_T_plus = s / m_B / (m_B + m_K) * form_factors->f_t(s);

        const complex<double> calH_plus = nonlocal_formfactor->H_plus(s);

        double F_Tkin = calF_T_plus / calF_plus * 2.0 * std::sqrt(lambda(s)) * beta_l(s) * m_B / s;
        double F_Skin = calF_time / calF_plus * 0.5 * (m_B2 - m_K2) / (m_b_MSbar - m_s_MSbar);

        // Wilson coefficients
        const complex<double>
            c7eff = ShortDistanceLowRecoil::c7eff(s, 0.0, 0.0, 0.0, false, wc); // LO C7eff
        const complex<double>
            c9_p  = wc.c9() + wc.c9prime(),
            c10_p = wc.c10() + wc.c10prime(),
            c7_p  = c7eff + wc.c7prime();

        // cf. [BHP:2007], Eq. (3.2), p. 3 and 4 or [BKMS:2012] (1205.5811)
        result.F_A  = c10_p;
        result.F_T  = F_Tkin * wc.cT();
        result.F_T5 = F_Tkin * wc.cT5();
        result.F_S  = F_Skin * (wc.cS() + wc.cSprime());
        result.F_P  = F_Skin * (wc.cP() + wc.cPprime()) + m_l() * c10_p *
                      ((m_B2 - m_K2) / s * (calF_time / calF_plus - 1.0) - 1.0);
        result.F_V  = c9_p
                      + 2.0 * m_b_MSbar() * m_B / s * c7_p * calF_T_plus / calF_plus
                      + 2.0 * m_b_PS() / m_B / xi_pseudo(s) * (dff.calT - 16.0 * power_of<2>(M_PI) * power_of<3>(m_B()) / m_b_PS() / s * calH_plus)
                      + 8.0 * m_l * m_B / s * calF_T_plus / calF_plus * wc.cT();

        return result;
    }
}
