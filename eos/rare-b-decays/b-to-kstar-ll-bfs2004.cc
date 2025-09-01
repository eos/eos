/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011      Christian Wacker
 * Copyright (c) 2014      Christoph Bobeth
 * Copyright (c) 2016-2025 Danny van Dyk
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

#include <eos/rare-b-decays/b-to-kstar-ll-bfs2004.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/nonlocal-form-factors/hard-scattering.hh>
#include <eos/nonlocal-form-factors/long-distance.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>

#include <functional>

#include <gsl/gsl_sf.h>
namespace eos
{
    using namespace std::literals::string_literals;
    using namespace std::placeholders;

    BToKstarDileptonAmplitudes<tag::BFS2004>::BToKstarDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c(p["mass::c"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        f_B(p["decay-constant::B_" + o.get("q"_ok, "d")], *this),
        f_Kstar_par(p["B->K^*::f_Kstar_par"], *this),
        f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"], *this),
        lambda_B_p_inv(p["B::1/lambda_B_p"], *this),
        a_1_par(p["K^*::a_1_para@1GeV"], *this),
        a_2_par(p["K^*::a_2_para@1GeV"], *this),
        a_1_perp(p["K^*::a_1_perp@1GeV"], *this),
        a_2_perp(p["K^*::a_2_perp@1GeV"], *this),
        uncertainty_para(p["B->K^*ll::A_para_uncertainty@LargeRecoil"], *this),
        uncertainty_perp(p["B->K^*ll::A_perp_uncertainty@LargeRecoil"], *this),
        uncertainty_long(p["B->K^*ll::A_long_uncertainty@LargeRecoil"], *this),
        uncertainty_xi_perp(p["formfactors::xi_perp_uncertainty"], *this),
        uncertainty_xi_par(p["formfactors::xi_par_uncertainty"], *this),
        q(o, options, "q"_ok),
        opt_ccbar_resonance(o, options, "ccbar-resonance"_ok),
        opt_use_nlo(o, options, "nlo"_ok),
        ccbar_resonance(opt_ccbar_resonance.value()),
        use_nlo(opt_use_nlo.value())
    {
        Context ctx("When constructing B->K^*ll BFS2004 amplitudes");

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

    BToKstarDileptonAmplitudes<tag::BFS2004>::~BToKstarDileptonAmplitudes()
    {
    }

    const std::vector<OptionSpecification>
    BToKstarDileptonAmplitudes<tag::BFS2004>::options
    {
        { "q"_ok, { "d"s, "u"s }, "d"s },
        { "ccbar-resonance"_ok, { "true"s, "false"s },  "false"s },
        { "nlo"_ok, { "true"s, "false"s },  "true"s },
    };


    BToKstarDilepton::DipoleFormFactors
    BToKstarDileptonAmplitudes<tag::BFS2004>::dipole_form_factors(const double & s, const WilsonCoefficients<BToS> & wc) const
    {
        // charges of down- and up-type quarks
        static const double e_d = -1.0/3.0;
        static const double e_u = +2.0/3.0;

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
        double invm1_par = 3.0 * (1.0 + a_1_par + a_2_par); // <ubar^-1>_par
        double invm1_perp = 3.0 * (1.0 + a_1_perp + a_2_perp); // <ubar^-1>_perp
        QCDFIntegrals<BToKstarDilepton> qcdf_0 = this->qcdf_dilepton_massless_case(s, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);
        QCDFIntegrals<BToKstarDilepton> qcdf_c = this->qcdf_dilepton_charm_case(s, m_c_pole, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);
        QCDFIntegrals<BToKstarDilepton> qcdf_b = this->qcdf_dilepton_bottom_case(s, m_b_PS, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);

        // inverse of the "negative" moment of the B meson LCDA
        // cf. [BFS:2001A], Eq. (54), p. 15
        double lambda_B_p_inv = this->lambda_B_p_inv;
        double omega_0 = 1.0 / lambda_B_p_inv;
        complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

        /* Y(s) for the up and the top sector */
        // cf. [BFS:2001A], Eq. (10), p. 4
        complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
        complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6());
        complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

        // Use b pole mass according to [BFS:2001A], Sec. 3.1, paragraph Quark Masses,
        // then replace b pole mass by the PS mass.
        complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
             + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
             + Y_top_0 * CharmLoops::h(mu, s)
             + Y_top_;
        // cf. [BFS:2004A], Eq. (43), p. 24
        complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

        /* Effective wilson coefficients */
        // cf. [BFS:2001A], below Eq. (9), p. 4
        complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
        // cf. [BFS:2001A], below Eq. (26), p. 8
        complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

        /* perpendicular, top sector */
        // cf. [BFS:2001A], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
        complex<double> C0_top_perp_left  = (c7eff - wc.c7prime()) + s / (2.0 * m_b_PS * m_B) * Y_top;
        complex<double> C0_top_perp_right = (c7eff + wc.c7prime()) + s / (2.0 * m_b_PS * m_B) * Y_top;
        // cf. [BFS:2004A], Eq. (44), p. 24
        complex<double> C1f_top_perp_left  = (c7eff - wc.c7prime()) * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
        complex<double> C1f_top_perp_right = (c7eff + wc.c7prime()) * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
        // cf. [BFS:2001A], Eqs. (34), (37), p. 9
        complex<double> C1nf_top_perp = (-1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (s / (2.0 * m_b_PS * m_B)) * (
                    wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                    + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

        /* perpendicular, up sector */
        // cf. [BFS:2004A], comment before Eq. (43), p. 24
        complex<double> C0_up_perp = s / (2.0 * m_b_PS * m_B) * Y_up;
        // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
        // cf. [BFS:2001A], Eqs. (34), (37), p. 9
        // [BFS:2004A], [S:2004A] have a different sign convention for F{12}{79}_massless than we!
        complex<double> C1nf_up_perp = (-1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                + (s / (2.0 * m_b_PS * m_B)) * (
                    wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                    + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

        /* parallel, top sector */
        // cf. [BFS:2001A], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
        complex<double> C0_top_par = -1.0 * (c7eff - wc.c7prime() + m_B / (2.0 * m_b_PS) * Y_top);
        // cf. [BFS:2004A], Eq. (45), p. 24
        complex<double> C1f_top_par = -1.0 * (c7eff - wc.c7prime()) * (8.0 * std::log(m_b_PS / mu) + 2.0 * L - 4.0 * (1.0 - mu_f() / m_b_PS));
        // cf. [BFS:2001A], Eqs. (38), p. 9
        complex<double> C1nf_top_par = (+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                    + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

        /* parallel, up sector */
        // cf. [BFS:2004A], comment before Eq. (43), p. 24
        complex<double> C0_up_par = -1.0 * m_B / (2.0 * m_b_PS) * Y_up;
        // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
        // cf. [BFS:2004A], last paragraph in Sec A.1, p. 24
        // [BFS:2004A], [S:2004A] have a different sign convention for F{12}{79}_massless than we!
        complex<double> C1nf_up_par = (+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                    + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

        // compute the factorizing contributions
        complex<double> C_perp_left  = C0_top_perp_left  + lambda_hat_u * C0_up_perp
            + a_mu * (C1f_top_perp_left  + C1nf_top_perp + lambda_hat_u * C1nf_up_perp);
        complex<double> C_perp_right = C0_top_perp_right + lambda_hat_u * C0_up_perp
            + a_mu * (C1f_top_perp_right + C1nf_top_perp + lambda_hat_u * C1nf_up_perp);
        complex<double> C_par = C0_top_par + lambda_hat_u * C0_up_par
            + a_mu * (C1f_top_par + C1nf_top_par + lambda_hat_u * C1nf_up_par);


        /* perpendicular, top sector */
        // T0_top_perp_{p,m} = 0, cf. [BFS:2001A], Eq. (17), p. 6
        // cf. [BFS:2004A], Eq. (49)
        complex<double> T1f_top_perp_p_left  = (c7eff - wc.c7prime()) * (2.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
        complex<double> T1f_top_perp_p_right = (c7eff + wc.c7prime()) * (2.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
        // T1f_top_perp_m = 0, cf. [BFS:2001A], Eq. (22), p. 7
        // cf. [BFS:2001A], Eq. (23), p. 7
        // [Christoph] Use c8 instead of c8eff
        complex<double> T1nf_top_perp_p = (-4.0 * e_d * c8eff * qcdf_0.j0_perp
            + m_B / (2.0 * m_b_PS) * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde1_perp
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0/3.0 * wc.c6() - (4.0 * m_b_PS / m_B) * (wc.c3() - wc.c4()/6.0 + 4.0 * wc.c5() - 2.0/3.0 * wc.c6())) * qcdf_b.jtilde1_perp
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() - 8.0/3.0 * wc.c6()) * qcdf_0.jtilde1_perp)) * lambda_B_p_inv;
        // T1nf_top_perp_m = 0, cf. [BFS:2001A], Eq. (17), p. 6

        /* perpendicular, up sector */
        // all T1f_up vanish, cf. [BFS:2004A], sentence below Eq. (49), p. 25
        // cf. [BFS:2004A], Eq. (50), p. 25
        complex<double> T1nf_up_perp_p = +e_u * m_B / (2.0 * m_b_PS) * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde1_perp - qcdf_0.jtilde1_perp) * lambda_B_p_inv;

        /* parallel, top sector */
        // T0_top_par_p = 0, cf. [BFS:2001A], Eq. (17), p. 6
        // cf. [BFS:2004A], Eqs. (46)-(47), p. 25 without the \omega term.
        complex<double> T0_top_par_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
        // cf. [BFS:2004A], Eq. (49), p. 25
        complex<double> T1f_top_par_p  = (c7eff - wc.c7prime()) * (4.0 * m_B / energy) * invm1_par * lambda_B_p_inv;
        // T1f_top_par_m = 0, cf. [BFS:2001A], Eq. (22), p. 7
        // cf. [BFS:2001A], Eq. (25), p. 7
        complex<double> T1nf_top_par_p = m_B / m_b_PS * (
                e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
        // cf. [BFS:2001A], Eq. (26), pp. 7-8
        complex<double> T1nf_top_par_m = e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                + 6.0 * m_B / m_b_PS * (
                    (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                    -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

        /* parallel, up sector */
        // all T1f_up vanish, cf. [BFS:2004A], sentence below Eq. (49), p. 25
        // cf. [BFS:2004A], Eqs. (46),(48), p. 25 without the \omega term
        complex<double> T0_up_par_m = +e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;
        // cf. [BFS:2004A], Eq. (50), p. 25
        complex<double> T1nf_up_par_p = +e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
        // cf. [BFS:2004A], Eq. (50), p. 25 without the \omega term
        complex<double> T1nf_up_par_m = +e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


        // Compute the nonfactorizing contributions
        complex<double> T_perp_left  = a_mu_f * (T1f_top_perp_p_left  + T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p);
        complex<double> T_perp_right = a_mu_f * (T1f_top_perp_p_right + T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p);
        complex<double> T_par = a_mu_f * (T1f_top_par_p + T1nf_top_par_p + lambda_hat_u * T1nf_up_par_p)
            + (T0_top_par_m + lambda_hat_u * T0_up_par_m + a_mu_f * (T1nf_top_par_m + lambda_hat_u * T1nf_up_par_m));

        // Compute the numerically leading power-suppressed weak annihilation contributions to order alpha_s^0
        // cf. [BFS:2004A], Eq. (51)
        complex<double> Delta_T_ann_top_perp = e_q * M_PI * M_PI * f_B / 3.0 / m_b_PS / m_B * (
                -4.0 * f_Kstar_perp * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 3.0 * wc.c5() + 4.0 * wc.c6())) * qcdf_0.j0_perp
                + 2.0 * f_Kstar_par * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 12.0 * wc.c5() + 16.0 * wc.c6())) *
                    (m_Kstar / (1.0 - s / (m_B * m_B)) * lambda_B_p_inv));
        complex<double> Delta_T_ann_up_perp = -e_q * 2.0 * M_PI * M_PI * f_B * f_Kstar_par / 3.0 / m_b_PS / m_B *
            (m_Kstar / (1.0 - s / (m_B * m_B)) * lambda_B_p_inv) * 3.0 * delta_qu * wc.c2();
        // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
        // cf. [BFS:2004A], Eqs. (52), (53)
        complex<double> Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                12.0 * c8eff * (m_b_PS / m_B) * f_Kstar_perp() * 1.0 / 3.0 * (qcdf_0.j0_perp + qcdf_0.j7_perp)
                + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (
                    (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j5_perp
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j5_perp
                    - (8.0 / 27.0) * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()) * qcdf_0.j0_perp)
                - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) * lambda_B_p_inv) * (3.0 / 4.0) * (
                    (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j6_perp
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j6_perp
                    - 8.0 / 27.0 * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6())));
        complex<double> Delta_T_hsa_up_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0) * (qcdf_c.j5_perp - qcdf_0.j5_perp)
                - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) * lambda_B_p_inv) * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0)
                    * (qcdf_c.j6_perp - qcdf_0.j6_perp));

        // Compute the sum of the numerically leading power-suppressed contributions
        complex<double> Delta_T_top_perp = Delta_T_ann_top_perp + Delta_T_hsa_top_perp;
        complex<double> Delta_T_up_perp = Delta_T_ann_up_perp + Delta_T_hsa_up_perp;
        complex<double> Delta_T_perp = Delta_T_top_perp + lambda_hat_u * Delta_T_up_perp;


        // cf. [BFS:2001A], Eq. (15), and [BHP:2008A], Eq. (C.4)
        BToKstarDilepton::DipoleFormFactors result;
        result.calT_perp_left  = xi_perp(s) * C_perp_left
            + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_left
            + Delta_T_perp;
        result.calT_perp_right = xi_perp(s) * C_perp_right
            + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_right
            + Delta_T_perp;
        result.calT_parallel = xi_par(s) * C_par
            + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_par * m_Kstar) / (m_B * energy) * T_par;

        return result;
    }

    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::xi_perp(const double & s) const
    {
        const double factor = m_B() / (m_B() + m_Kstar());
        double result = uncertainty_xi_perp * factor * form_factors->v(s);

        return result;
    }

    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::xi_par(const double & s) const
    {
        const double factor1 = (m_B() + m_Kstar()) / (2.0 * energy(s));
        const double factor2 = (1.0 - m_Kstar() / m_B());
        double result = uncertainty_xi_par * (factor1 * form_factors->a_1(s) - factor2 * form_factors->a_2(s));

        return result;
    }

    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::norm(const double & s) const
    {
        double lambda_t2 = std::norm(model->ckm_tb() * conj(model->ckm_ts()));

        return g_fermi() * alpha_e() * std::sqrt(
                  1.0 / 3.0 / 1024 / power_of<5>(M_PI) / m_B()
                  * lambda_t2 * s_hat(s) * std::sqrt(lambda(s)) * beta_l(s)
               ); // cf. [BHP:2008A], Eq. (C.6), p. 21
    }

    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::mu_f() const
    {
        return 1.5;
    }

    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    /* Amplitudes */
    // cf. [BHP:2008A], p. 20
    // cf. [BHvD:2012A], app B, eqs. (B13 - B19)
    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::BFS2004>::amplitudes(const double & s) const
    {
        BToKstarDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const double
            shat = s_hat(s),
            mbhat = m_b_PS() / m_B,
            mKhat2 = power_of<2>(m_Kstar() / m_B()),
            m_K2 = power_of<2>(m_Kstar()),
            m_B2 = power_of<2>(m_B()),
            m2_diff = m_B2 - m_K2,
            norm_s = this->norm(s),
            sqrt_lam = std::sqrt(lambda(s)),
            sqrt_s = std::sqrt(s);

        auto dff = dipole_form_factors(s, wc);

        const complex<double>
            wilson_minus_right = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime()),
            wilson_minus_left  = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime()),
            wilson_plus_right  = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime()),
            wilson_plus_left   = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime());

        // longitudinal amplitude
        const double prefactor_long = -norm_s / (2.0 * m_Kstar() * std::sqrt(s));

        const complex<double>
            a = (m2_diff - s) * 2.0 * energy(s) * xi_perp(s) - lambda(s) * m_B() / m2_diff * (xi_perp(s) - xi_par(s)),
            b = 2.0 * m_b_PS() * (
                    ((m_B2 + 3.0 * m_K2 - s) * 2.0 * energy(s) / m_B() - lambda(s) / m2_diff) * dff.calT_perp_left
                    - lambda(s) / m2_diff * dff.calT_parallel
                );

        result.a_long_right = prefactor_long * (wilson_minus_right * a + uncertainty_long() * b);
        result.a_long_left  = prefactor_long * (wilson_minus_left  * a + uncertainty_long() * b);

        // perpendicular amplitude
        const double prefactor_perp = +std::sqrt(2.0) * norm_s * m_B() * std::sqrt(eos::lambda(1.0, mKhat2, shat));

        result.a_perp_right = prefactor_perp * (wilson_plus_right * xi_perp(s) + uncertainty_perp() * (2.0 * mbhat / shat) * dff.calT_perp_right);
        result.a_perp_left  = prefactor_perp * (wilson_plus_left  * xi_perp(s) + uncertainty_perp() * (2.0 * mbhat / shat) * dff.calT_perp_right);

        // parallel amplitude
        const double prefactor_par = -std::sqrt(2.0) * norm_s * m2_diff;

        result.a_para_right = prefactor_par * (
                                wilson_minus_right * xi_perp(s) * 2.0 * energy(s) / m2_diff
                                + uncertainty_para() * 4.0 * m_b_PS() * energy(s) / s / m_B() * dff.calT_perp_left
                             );
        result.a_para_left  = prefactor_par * (
                                wilson_minus_left  * xi_perp(s) * 2.0 * energy(s) / m2_diff
                                + uncertainty_para() * 4.0 * m_b_PS() * energy(s) / s / m_B() * dff.calT_perp_left
                             );

        // timelike amplitude
        result.a_time = norm_s * sqrt_lam / sqrt_s
            * (2.0 * (wc.c10() - wc.c10prime()) + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime()))
            * form_factors->a_0(s);

        // scalar amplitude
        result.a_scal = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar) * form_factors->a_0(s);

        // tensor amplitudes [BHvD:2012A]  eqs. (B18 - B20)
        // no form factor relations used
        const double
            ff_T1  = form_factors->t_1(s),
            ff_T2  = form_factors->t_2(s),
            ff_T3  = form_factors->t_3(s),

            kin_tensor_1 = norm_s / m_Kstar() * ((m_B2 + 3.0 * m_K2 - s) * ff_T2 - lambda(s) / m2_diff * ff_T3),
            kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1,
            kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

        // correct the sign of C_T5 from [BHvD2012v4] because of inconsistent use of gamma5 <-> Levi-Civita
        static const double sign = -1;

        result.a_para_perp = kin_tensor_1 * wc.cT();
        result.a_time_long = kin_tensor_1 * sign * wc.cT5();

        result.a_time_perp = kin_tensor_2 * wc.cT();
        result.a_long_perp = kin_tensor_2 * sign * wc.cT5();

        result.a_time_para = kin_tensor_3 * sign * wc.cT5();
        result.a_long_para = kin_tensor_3 * wc.cT();

        return result;
    }

    // C9 and its corrections [BFS:2001A] eqs. (40-41)
    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::real_C9_perp(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const double
            shat = s_hat(s),
            mbhat = m_b_PS() / m_B;

        auto dff = dipole_form_factors(s, wc);

        return real(wc.c9() + uncertainty_perp() * (2.0 * mbhat / shat) * dff.calT_perp_right / xi_perp(s));
    }
    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::imag_C9_perp(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const double
            shat = s_hat(s),
            mbhat = m_b_PS() / m_B;

        auto dff = dipole_form_factors(s, wc);

        return imag(wc.c9() + uncertainty_perp() * (2.0 * mbhat / shat) * dff.calT_perp_right / xi_perp(s));
    }
    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::real_C9_para(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        auto dff = dipole_form_factors(s, wc);

        return real(wc.c9() - uncertainty_para() * 2.0 * m_b_PS() / m_B() * dff.calT_parallel / xi_par(s));
    }
    double
    BToKstarDileptonAmplitudes<tag::BFS2004>::imag_C9_para(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        auto dff = dipole_form_factors(s, wc);

        return imag(wc.c9() - uncertainty_para() * 2.0 * m_b_PS() / m_B() * dff.calT_parallel / xi_par(s));
    }
}
