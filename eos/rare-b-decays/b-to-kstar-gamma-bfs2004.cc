/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2015, 2016, 2017 Danny van Dyk
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

#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma-bfs2004.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/nonlocal-form-factors/hard-scattering.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>

#include <cmath>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;

    BToKstarGammaAmplitudes<tag::BFS2004>::BToKstarGammaAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        hbar(p["QM::hbar"], *this),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c(p["mass::c"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        f_B(p["decay-constant::B_" + o.get("q", "d")], *this),
        f_Kstar_par(p["B->K^*::f_Kstar_par"], *this),
        f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"], *this),
        lambda_B_p_inv(p["B::1/lambda_B_p"], *this),
        a_1_para(p["K^*::a_1_para@1GeV"], *this),
        a_2_para(p["K^*::a_2_para@1GeV"], *this),
        a_1_perp(p["K^*::a_1_perp@1GeV"], *this),
        a_2_perp(p["K^*::a_2_perp@1GeV"], *this),
        uncertainty_para(p["B->K^*ll::A_para_uncertainty@LargeRecoil"], *this),
        uncertainty_perp(p["B->K^*ll::A_perp_uncertainty@LargeRecoil"], *this),
        form_factors(FormFactorFactory<PToV>::create("B->K^*::" + o.get("form-factors", "BSZ2015"), p)),
        mu(p["sb::mu"], *this)
    {
        Context ctx("When constructing B->K^*gamma BFS2004 amplitudes");

        this->uses(*model);
        this->uses(*form_factors);

        // Select the appropriate calculator for the QCDF integrals
        std::string qcdf_integrals(o.get("qcdf-integrals", "mixed"));
        if ("mixed" == qcdf_integrals)
        {
            qcdf_photon_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::photon_massless_case,
                        _1, _2, _3, _4, _5, _6, _7);
            qcdf_photon_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::photon_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_photon_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::photon_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
        }
        else if ("analytical" == qcdf_integrals)
        {
            qcdf_photon_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_massless_case,
                        _1, _2, _3, _4, _5, _6, _7);
            qcdf_photon_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_photon_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
        }
        else
        {
            throw InvalidOptionValueError("qcdf-integrals", qcdf_integrals, "mixed, analytical");
        }
    }

    inline double
    BToKstarGammaAmplitudes<tag::BFS2004>::xi_perp() const
    {
        const double factor = m_B / (m_B + m_Kstar);
        double result = factor * form_factors->v(0.0);

        return result;
    }

    inline double
    BToKstarGammaAmplitudes<tag::BFS2004>::mu_f() const
    {
        return 1.5;
    }

    inline double
    BToKstarGammaAmplitudes<tag::BFS2004>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    BToKstarGamma::Amplitudes
    BToKstarGammaAmplitudes<tag::BFS2004>::amplitudes() const
    {
        // charges of down- and up-type quarks
        static const double e_d = -1.0/3.0;
        static const double e_u = +2.0/3.0;

        // spectator contributions
        double delta_qu = (q.value() == QuarkFlavor::up ? 1.0 : 0.0);

        // kinematics
        double m_c_pole = model->m_c_pole();
        double m_b_PS = this->m_b_PS();
        double energy = (m_B * m_B + m_Kstar * m_Kstar) / 2.0 / m_B;
        // L from B->K^*ll for s -> 0
        double L = +1.0;

        // couplings
        double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
        double a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI;
        double alpha_s_mu_f = model->alpha_s(std::sqrt(mu() * 0.5)); // alpha_s at the factorization scale
        double a_mu_f = alpha_s_mu_f * QCD::casimir_f / 4.0 / M_PI;
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
            lambda_hat_u = std::conj(lambda_hat_u);
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), LeptonFlavor::muon /*fake lepton flavor*/, cp_conjugate);

        // Compute the QCDF Integrals
        double invm1_perp = 3.0 * (1.0 + a_1_perp + a_2_perp); // <ubar^-1>_perp
        QCDFIntegrals<BToKstarDilepton> qcdf_0 = this->qcdf_photon_massless_case(m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
        QCDFIntegrals<BToKstarDilepton> qcdf_c = this->qcdf_photon_charm_case(m_c_pole, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);
        QCDFIntegrals<BToKstarDilepton> qcdf_b = this->qcdf_photon_bottom_case(m_b_PS, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_para, a_2_para);

        // inverse of the "negative" moment of the B meson LCDA
        // cf. [BFS2001], Eq. (54), p. 15
        double lambda_B_p_inv = this->lambda_B_p_inv;

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
        // cf. [BFS2001], below Eq. (26), p. 8
        complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

        /* perpendicular, top sector */
        // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1, s -> 0, +/- -> left/right handed
        complex<double> C0_top_perp_left  = c7eff;
        complex<double> C0_top_perp_right = wc.c7prime();
        // cf. [BFS2004], Eq. (44), p. 24, s -> 0
        complex<double> C1f_top_perp_left  = c7eff * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
        complex<double> C1f_top_perp_right = wc.c7prime() * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
        // cf. [BFS2001], Eqs. (34), (37), p. 9, s -> 0
        complex<double> C1nf_top_perp_left = (-1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), 0.0, m_b_PS, m_c_pole) + c8eff * CharmLoops::F87_massless(mu, 0.0, m_b_PS));
        const complex<double> C1nf_top_perp_right = 0.0;

        /* perpendicular, up sector */
        // cf. [BFS2004], comment before Eq. (43), p. 24, s -> 0
        const complex<double> C0_up_perp_left = 0.0;
        const complex<double> C0_up_perp_right = 0.0;
        // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
        // cf. [BFS2001], Eqs. (34), (37), p. 9
        // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
        complex<double> C1nf_up_perp_left = (-1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), 0.0, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, 0.0, m_b_PS)));
        const complex<double> C1nf_up_perp_right = 0.0;

        // compute the factorizing contributions
        complex<double> C_perp_left  = C0_top_perp_left  + lambda_hat_u * C0_up_perp_left
            + a_mu * (C1f_top_perp_left  + C1nf_top_perp_left + lambda_hat_u * C1nf_up_perp_left);
        complex<double> C_perp_right = C0_top_perp_right + lambda_hat_u * C0_up_perp_right
            + a_mu * (C1f_top_perp_right + C1nf_top_perp_right + lambda_hat_u * C1nf_up_perp_right);

        /* perpendicular, top sector */
        // T0_top_perp_{p,m} = 0, cf. [BFS2001], Eq. (17), p. 6
        // cf. [BFS2004], Eq. (49)
        complex<double> T1f_top_perp_p_left  = c7eff * (4.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
        complex<double> T1f_top_perp_p_right = wc.c7prime() * (4.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
        // T1f_top_perp_m = 0, cf. [BFS2001], Eq. (22), p. 7
        // cf. [BFS2001], Eq. (23), p. 7
        // [Christoph] Use c8 instead of c8eff
        complex<double> T1nf_top_perp_p_left = (-4.0 * e_d * c8eff * qcdf_0.j0bar_perp
            + m_B / (2.0 * m_b_PS) * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde1_perp
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0/3.0 * wc.c6() - (4.0 * m_b_PS / m_B) * (wc.c3() - wc.c4()/6.0 + 4.0 * wc.c5() - 2.0/3.0 * wc.c6())) * qcdf_b.jtilde1_perp
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() - 8.0/3.0 * wc.c6()) * qcdf_0.jtilde1_perp)) * lambda_B_p_inv;
        const complex<double> T1nf_top_perp_p_right = 0.0;
        // T1nf_top_perp_m = 0, cf. [BFS2001], Eq. (17), p. 6

        /* perpendicular, up sector */
        // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
        // cf. [BFS2004], Eq. (50), p. 25
        complex<double> T1nf_up_perp_p_left = +e_u * m_B / (2.0 * m_b_PS) * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde1_perp - qcdf_0.jtilde1_perp) * lambda_B_p_inv;
        const complex<double> T1nf_up_perp_p_right = 0.0;


        // Compute the nonfactorizing contributions
        complex<double> T_perp_left  = a_mu_f * (T1f_top_perp_p_left + T1nf_top_perp_p_left + lambda_hat_u * T1nf_up_perp_p_left);
        complex<double> T_perp_right = a_mu_f * (T1f_top_perp_p_right + T1nf_top_perp_p_right + lambda_hat_u * T1nf_up_perp_p_right);

        // Compute the numerically leading power-suppressed weak annihilation contributions to order alpha_s^0
        // cf. [BFS2004], Eq. (51)
        complex<double> Delta_T_ann_top_perp = e_q * M_PI * M_PI * f_B / 3.0 / m_b_PS / m_B * (
                -4.0 * f_Kstar_perp * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 3.0 * wc.c5() + 4.0 * wc.c6())) * qcdf_0.j0_perp
                + 2.0 * f_Kstar_par * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 12.0 * wc.c5() + 16.0 * wc.c6())) *
                    (m_Kstar * lambda_B_p_inv));
        complex<double> Delta_T_ann_up_perp = -e_q * 2.0 * M_PI * M_PI * f_B * f_Kstar_par / 3.0 / m_b_PS / m_B *
            (m_Kstar * lambda_B_p_inv) * 3.0 * delta_qu * wc.c2();
        // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
        // cf. [BFS2004], Eqs. (52), (53)
        complex<double> Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                12.0 * c8eff * (m_b_PS / m_B) * f_Kstar_perp() * 1.0 / 3.0 * (qcdf_0.j0_perp + qcdf_0.j7_perp)
                + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (
                    (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j5_perp
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j5_perp
                    - (8.0 / 27.0) * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()) * qcdf_0.j0_perp)
                - (4.0 * m_Kstar * f_Kstar_par * lambda_B_p_inv) * (3.0 / 4.0) * (
                    (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j6_perp
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j6_perp
                    - 8.0 / 27.0 * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6())));
        complex<double> Delta_T_hsa_up_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0) * (qcdf_c.j5_perp - qcdf_0.j5_perp)
                - (4.0 * m_Kstar * f_Kstar_par * lambda_B_p_inv) * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0)
                    * (qcdf_c.j6_perp - qcdf_0.j6_perp));

        // Compute the sum of the numerically leading power-suppressed contributions
        complex<double> Delta_T_top_perp = Delta_T_ann_top_perp + Delta_T_hsa_top_perp;
        complex<double> Delta_T_up_perp = Delta_T_ann_up_perp + Delta_T_hsa_up_perp;
        complex<double> Delta_T_perp = Delta_T_top_perp + lambda_hat_u * Delta_T_up_perp;

        // Form factor
        double xi_perp_zero = xi_perp();

        double calN = sqrt(alpha_e * power_of<3>(m_B) * power_of<3>(1.0 -  m_Kstar * m_Kstar / (m_B * m_B)) / (32.0 * power_of<4>(M_PI))) * g_fermi
                * model->m_b_msbar(mu())
                * abs(model->ckm_tb() * conj(model->ckm_ts()));

        // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
        complex<double> a_left  = complex<double>(0.0, +1.0) * calN * (
                xi_perp_zero * C_perp_left
                + (uncertainty_perp() + uncertainty_para()) / 2.0 * (
                    power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_left
                    + Delta_T_perp
                )
                + (uncertainty_perp() - uncertainty_para()) / 2.0 * (
                    power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_right
                    + Delta_T_perp
                )
                );
        complex<double> a_right = complex<double>(0.0, -1.0) * calN * (
                xi_perp_zero * C_perp_right
                + (uncertainty_perp() - uncertainty_para()) / 2.0 * (
                    power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_left
                    + Delta_T_perp
                )
                + (uncertainty_perp() + uncertainty_para()) / 2.0 * (
                    power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_right
                    + Delta_T_perp
                )
                );

        return BToKstarGamma::Amplitudes{
                (a_left - a_right) / sqrt(2.0),
                (a_left + a_right) / sqrt(2.0),
            };
    }
}
