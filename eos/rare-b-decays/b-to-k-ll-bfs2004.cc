/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2011      Christian Wacker
 * Copyright (c) 2014      Frederik Beaujean
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
#include <eos/rare-b-decays/b-to-k-ll-bfs2004.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::literals::string_literals;
    using namespace std::placeholders;

    BToKDileptonAmplitudes<tag::BFS2004>::BToKDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c(p["mass::c"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        f_B(p["decay-constant::B_" + o.get("q"_ok, "d")], *this),
        f_K(p["decay-constant::K_" + o.get("q"_ok, "d")], *this),
        lambda_B_p_inv(p["B::1/lambda_B_p"], *this),
        a_1(p["K::a_1@1GeV"], *this),
        a_2(p["K::a_2@1GeV"], *this),
        lambda_psd(p["B->Pll::Lambda_pseudo@LargeRecoil"], *this),
        sl_phase_psd(p["B->Pll::sl_phase_pseudo@LargeRecoil"], *this),
        q(o, options, "q"_ok)
    {
        Context ctx("When constructing B->Kll BFS2004 amplitudes");

        switch (q.value())
        {
            case QuarkFlavor::down:
                e_q = -1.0 / 3.0;
                break;

            case QuarkFlavor::up:
                e_q = 2.0 / 3.0;
                break;

            default:
                throw InternalError("Unexpected quark flavor: '" + q.str() + "'");
        }

        // Select the appropriate calculator for the QCDF integrals
        std::string qcdf_integrals(o.get("qcdf-integrals"_ok, "mixed"));
        if ("mixed" == qcdf_integrals)
        {
            qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_massless_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
            qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
        }
        else if ("numerical" == qcdf_integrals)
        {
            qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_massless_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
            qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
        }
        else if ("analytical" == qcdf_integrals)
        {
            qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_massless_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
            qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
        }
        else
        {
            throw InvalidOptionValueError("qcdf-integrals"_ok, qcdf_integrals, "mixed, numerical, analytical");
        }
    }

    BToKDileptonAmplitudes<tag::BFS2004>::~BToKDileptonAmplitudes()
    {
    }

    const std::vector<OptionSpecification>
    BToKDileptonAmplitudes<tag::BFS2004>::options
    {
        { "q"_ok, { "d"s, "u"s }, "d"s },
    };

    BToKDilepton::DipoleFormFactors
    BToKDileptonAmplitudes<tag::BFS2004>::dipole_form_factors(const double & s, const WilsonCoefficients<BToS> & wc) const
    {
        // charges of down- and up-type quarks
        static const double e_d = -1.0 / 3.0;
        static const double e_u = +2.0 / 3.0;

        // spectator contributions
        double delta_qu = (q.value() == QuarkFlavor::up ? 1.0 : 0.0);

        // kinematics
        double m_c_pole = model->m_c_pole();
        double m_b_PS = this->m_b_PS(), m_b_PS2 = m_b_PS * m_b_PS;
        double energy = this->energy(s);
        double L = -1.0 * (m_b_PS2 - s) / s * std::log(1.0 - s / m_b_PS2);

        // couplings
        double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
        double a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI;
        double alpha_s_mu_f = model->alpha_s(std::sqrt(mu() * 0.5)); // alpha_s at the factorization scale
        double a_mu_f = alpha_s_mu_f * QCD::casimir_f / 4.0 / M_PI;
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
            lambda_hat_u = std::conj(lambda_hat_u);

        // Compute the QCDF Integrals
        double invm1_psd = 3.0 * (1.0 + a_1 + a_2); // <ubar^-1>
        QCDFIntegrals<BToKstarDilepton> qcdf_0 = this->qcdf_dilepton_massless_case(s, m_B, m_K, mu, 0.0, 0.0, a_1, a_2);
        QCDFIntegrals<BToKstarDilepton> qcdf_c = this->qcdf_dilepton_charm_case(s, m_c_pole, m_B, m_K, mu, 0.0, 0.0, a_1, a_2);
        QCDFIntegrals<BToKstarDilepton> qcdf_b = this->qcdf_dilepton_bottom_case(s, m_b_PS, m_B, m_K, mu, 0.0, 0.0, a_1, a_2);

        // inverse of the "negative" moment of the B meson LCDA
        // cf. [BFS2001], Eq. (54), p. 15
        double lambda_B_p_inv = this->lambda_B_p_inv, omega_0 = 1.0 / this->lambda_B_p_inv;
        complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

        /* Y(s) for the up and the top sector */
        // cf. [BFS2001], Eq. (10), p. 4
        complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
        complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

        // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
        // then replace b pole mass by the PS mass.
        complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole);
        Y_top += Y_top_b * CharmLoops::h(mu, s, m_b_PS);
        Y_top += Y_top_0 * CharmLoops::h(mu, s);
        Y_top += Y_top_;
        // cf. [BFS2004], Eq. (43), p. 24
        complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
        // cf. [BFS2001], below Eq. (26), p. 8
        complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

        /* top sector */
        // cf. [BHP2007], Eq. (B.2) and [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
        complex<double> C0_top_psd = 1.0 * (c7eff + wc.c7prime() + m_B / (2.0 * m_b_PS) * Y_top);
        // cf. [BHP2007], Eq. (B.2) and [BFS2004], Eq. (45), p. 24
        // the correct sign in front of C_7^eff is plus, as one can see by
        // comparison with [BF2001], Eq. (63)
        complex<double> C1f_top_psd = 1.0 * (c7eff + wc.c7prime()) * (8.0 * std::log(m_b_PS / mu) + 2.0 * L - 4.0 * (1.0 - mu_f() / m_b_PS));
        // cf. [BHP2007], Eq. (B.2) and [BFS2001], Eqs. (38), p. 9
        complex<double> C1nf_top_psd = -(+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                    + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

        /* parallel, up sector */
        // cf. [BHP2007], Eq. (B.2) and [BFS2004], comment before Eq. (43), p. 24
        complex<double> C0_up_psd = 1.0 * m_B / (2.0 * m_b_PS) * Y_up;
        // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
        // cf. [BFS2004], last paragraph in Sec A.1, p. 24
        // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
        // Use here FF_massive - FF_massless because FF_massless is defined with an extra '-'
        // compared to [S2004]
        complex<double> C1nf_up_psd = -(+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                    + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

        // compute the factorizing contributions
        complex<double> C_psd = C0_top_psd + lambda_hat_u * C0_up_psd
            + a_mu * (C1f_top_psd + C1nf_top_psd + lambda_hat_u * C1nf_up_psd);

        /* parallel, top sector */
        // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
        // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
        complex<double> T0_top_psd_m = +e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
        // cf. [BHP2007], Eq. (B.2)
        complex<double> T1f_top_psd_p  = -(c7eff + wc.c7prime()) * (4.0 * m_B / energy) * invm1_psd * lambda_B_p_inv;
        // T1f_top_par_m = 0, cf. [BFS2001], Eq. (22), p. 7
        // cf. [BFS2001], Eq. (25), p. 7
        complex<double> T1nf_top_psd_p = -m_B / m_b_PS * (
                e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
        // cf. [BFS2001], Eq. (26), pp. 7-8
        complex<double> T1nf_top_psd_m = -e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                + 6.0 * m_B / m_b_PS * (
                    (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                    -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

        /* parallel, up sector */
        // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
        // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
        complex<double> T0_up_psd_m = -e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;
        // cf. [BFS2004], Eq. (50), p. 25
        complex<double> T1nf_up_psd_p = -e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
        // cf. [BFS2004], Eq. (50), p. 25 without the \omega term
        complex<double> T1nf_up_psd_m = -e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


        // Compute the nonfactorizing contributions
        complex<double> T_psd = a_mu_f * (T1f_top_psd_p + T1nf_top_psd_p + lambda_hat_u * T1nf_up_psd_p)
            + (T0_top_psd_m + lambda_hat_u * T0_up_psd_m + a_mu_f * (T1nf_top_psd_m + lambda_hat_u * T1nf_up_psd_m));

        // Subleading weak annihilation and hard spectator interaction contributions have only been
        // computed for calT_perp, not for calT_par ~ calT_psd.

        // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
        BToKDilepton::DipoleFormFactors result;
        result.calT = xi_pseudo(s) * C_psd + power_of<2>(M_PI) / 3.0 * (f_B * f_K) / m_B  * T_psd;

        return result;
    }

    double
    BToKDileptonAmplitudes<tag::BFS2004>::xi_pseudo(const double & s) const
    {
        // cf. [BF2001], Eq. (22)
        return form_factors->f_p(s);
    }

    double
    BToKDileptonAmplitudes<tag::BFS2004>::mu_f() const
    {
        return 1.5;
    }

    double
    BToKDileptonAmplitudes<tag::BFS2004>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    /* Amplitudes */
    BToKDilepton::Amplitudes
    BToKDileptonAmplitudes<tag::BFS2004>::amplitudes(const double & s) const
    {
        BToKDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        auto dff = dipole_form_factors(s, wc);

        // cf. [BF2001] Eq. (22 + TODO: 31)
        // cf. [BF2001] Eq. (22 + TODO: 30)
        double f_t_over_f_p = form_factors->f_t(s) / form_factors->f_p(s);
        double f_0_over_f_p = form_factors->f_0(s) / form_factors->f_p(s);

        double F_Tkin = f_t_over_f_p * 2.0 * std::sqrt(lambda(s)) * beta_l(s) / (m_B() + m_K());
        double F_Skin = f_0_over_f_p * 0.5 * (power_of<2>(m_B()) - power_of<2>(m_K())) / (m_b_MSbar - m_s_MSbar);

        // cf. [BHP2007], Eq. (3.2), p. 3 and 4
        result.F_A  = wc.c10() + wc.c10prime();
        result.F_T  = F_Tkin * wc.cT();
        result.F_T5 = F_Tkin * wc.cT5();
        result.F_S  = F_Skin * (wc.cS() + wc.cSprime());
        result.F_P  = F_Skin * (wc.cP() + wc.cPprime()) + m_l() * (wc.c10() + wc.c10prime()) *
                      ((m_B() * m_B() - m_K() * m_K()) / s * (f_0_over_f_p - 1.0) - 1.0);
        result.F_V  = wc.c9() + wc.c9prime()
                      + 2.0 * m_b_PS() / m_B() / xi_pseudo(s) * (dff.calT + lambda_psd / m_B * std::polar(1.0, sl_phase_psd()))
                      + 8.0 * m_l / (m_B() + m_K()) * f_t_over_f_p * wc.cT();

        return result;
    }
}
