/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/hard-scattering.hh>
#include <eos/rare-b-decays/long-distance.hh>
#include <eos/rare-b-decays/qcdf_integrals.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/save.hh>

#include <cmath>
#include <functional>

#include <gsl/gsl_sf.h>

#include <iostream>

namespace eos
{
    using namespace eos::btovll;
    using std::norm;

    /*
     * Decay: B -> K^* l lbar at Large Recoil, cf. [BHP2008]
     */
    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        UsedParameter hbar;

        UsedParameter m_b_MSbar;

        UsedParameter m_c;

        UsedParameter m_s_MSbar;

        UsedParameter m_B;

        UsedParameter m_Kstar;

        UsedParameter m_l;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter f_B;

        UsedParameter f_Kstar_par;

        UsedParameter f_Kstar_perp;

        UsedParameter lambda_B_p;

        UsedParameter a_1_par;

        UsedParameter a_2_par;

        UsedParameter a_1_perp;

        UsedParameter a_2_perp;

        UsedParameter uncertainty_para;

        UsedParameter uncertainty_perp;

        UsedParameter uncertainty_long;

        UsedParameter uncertainty_xi_perp;

        UsedParameter uncertainty_xi_par;

        UsedParameter tau;

        double e_q;

        char q;

        std::string lepton_flavour;

        bool cp_conjugate;

        std::string ff_relation;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "WilsonScan"), p, o)),
            parameters(p),
            hbar(p["hbar"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            m_c(p["mass::c"], u),
            m_s_MSbar(p["mass::s(2GeV)"], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            m_Kstar(p["mass::K^*_d"], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            f_B(p["decay-constant::B_" + o.get("q", "d")], u),
            f_Kstar_par(p["B->K^*::f_Kstar_par"], u),
            f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"], u),
            lambda_B_p(p["lambda_B_p"], u),
            a_1_par(p["B->K^*::a_1_par"], u),
            a_2_par(p["B->K^*::a_2_par"], u),
            a_1_perp(p["B->K^*::a_1_perp"], u),
            a_2_perp(p["B->K^*::a_2_perp"], u),
            uncertainty_para(p["B->K^*ll::A_para_uncertainty@LargeRecoil"], u),
            uncertainty_perp(p["B->K^*ll::A_perp_uncertainty@LargeRecoil"], u),
            uncertainty_long(p["B->K^*ll::A_long_uncertainty@LargeRecoil"], u),
            uncertainty_xi_perp(p["formfactors::xi_perp_uncertainty"], u),
            uncertainty_xi_par(p["formfactors::xi_par_uncertainty"], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            e_q(-1.0/3.0),
            lepton_flavour(o.get("l", "mu")),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            if (0.0 == m_l())
            {
                throw InternalError("Zero lepton mass leads to NaNs in timelike amplitudes. Use tiny lepton mass > 0!");
            }

            form_factors = FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "KMPW2010"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

            std::string spectator_quark = o.get("q", "d");
            if (spectator_quark.size() != 1)
                throw InternalError("Option q should only be one character!");

            q = spectator_quark[0];
            if (q == 'd')
            {
                e_q = -1.0 / 3.0;
            }
            else if (q == 'u')
            {
                e_q = 2.0 / 3.0;
            }
            else
            {
                throw InvalidOptionValueError("q", spectator_quark, "u, d");
            }

            ff_relation = o.get("large-recoil-ff", "BFS2004");

#if 0
            p["Abs{c7}"]  =  0.33670; // c7eff = 0.30726
            p["Abs{c9}"]  =  4.27305;
            p["Abs{c10}"] =  4.17166;
            p["c1"]       = -0.32196;
            p["c2"]       =  1.00925;
            p["c3"]       = -0.00519;
            p["c4"]       = -0.08787;
            p["c5"]       =  0.00036;
            p["c6"]       =  0.00101;
            p["c8"]       = -0.18262; // c8eff = -0.16925
            p["decay-constant::B_d"] = 0.200;
            p["decay-constant::B_u"] = 0.200;
            p["B->K^*::f_Kstar_perp@2GeV"] = 0.165138862;

            Amplitudes a = amplitudes(3.0);
            std::cout << "A_pp^L(q^2 = 3) = " << a.a_perp_left << std::endl;
            std::cout << "A_pp^R(q^2 = 3) = " << a.a_perp_right << std::endl;
            std::cout << "A_pa^L(q^2 = 3) = " << a.a_par_left << std::endl;
            std::cout << "A_pa^R(q^2 = 3) = " << a.a_par_right << std::endl;
            std::cout << "A_00^L(q^2 = 3) = " << a.a_long_left << std::endl;
            std::cout << "A_00^R(q^2 = 3) = " << a.a_long_right << std::endl;

            std::cout << "lam_up    = " << (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts())) << std::endl;
            std::cout << "xi_par(3) = " << xi_par(3.0) << std::endl;
            std::cout << "xi_perp(3)= " << xi_perp(3.0) << std::endl;
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s();
            std::cout << "C9       = " << wc.c9() << std::endl;
            std::cout << "C10      = " << wc.c10() << std::endl;
            std::cout << "C9 - C10 = " << wc.c9() - wc.c10() << std::endl;
            std::cout << "C9 + c10 = " << wc.c9() + wc.c10() << std::endl;
            DipoleFormFactors dff = calT(3.0);
            std::cout << "s = 3" << std::endl;
            std::cout << "calT_perp_left  = " << dff.calT_perp_left << std::endl;
            std::cout << "calT_perp_right = " << dff.calT_perp_right << std::endl;
            std::cout << "calT_parallel   = " << dff.calT_parallel << std::endl;
            dff = calT(0.0);
            std::cout << "s = 0" << std::endl;
            std::cout << "calT_perp_left  = " << dff.calT_perp_left << std::endl;
            std::cout << "calT_perp_right = " << dff.calT_perp_right << std::endl;
            std::cout << "calT_parallel   = " << dff.calT_parallel << std::endl;
            throw std::string("foo");
#endif
        }

        WilsonCoefficients<BToS> wilson_coefficients() const
        {
            return model->wilson_coefficients_b_to_s(lepton_flavour, cp_conjugate);
        }

        struct DipoleFormFactors
        {
            complex<double> calT_perp_left;
            complex<double> calT_perp_right;
            complex<double> calT_parallel;
        };

        DipoleFormFactors calT_BFS2004(const double & s, const WilsonCoefficients<BToS> & wc) const
        {
            // charges of down- and up-type quarks
            static const double e_d = -1.0/3.0;
            static const double e_u = +2.0/3.0;

            // spectator contributions
            double delta_qu = (q == 'u' ? 1.0 : 0.0);

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
            QCDFIntegrals::Results qcdf_0 = QCDFIntegrals::dilepton_massless_case(s, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);
            QCDFIntegrals::Results qcdf_c = QCDFIntegrals::dilepton_charm_case(s, m_c_pole, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);
            QCDFIntegrals::Results qcdf_b = QCDFIntegrals::dilepton_bottom_case(s, m_b_PS, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);

            // inverse of the "negative" moment of the B meson LCDA
            // cf. [BFS2001], Eq. (54), p. 15
            double omega_0 = lambda_B_p, lambda_B_p_inv = 1.0 / lambda_B_p;
            complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

            /* Y(s) for the up and the top sector */
            // cf. [BFS2001], Eq. (10), p. 4
            complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
            complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
            complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6());
            complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

            // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
            // then replace b pole mass by the PS mass.
            complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
                 + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
                 + Y_top_0 * CharmLoops::h(mu, s)
                 + Y_top_;
            // cf. [BFS2004], Eq. (43), p. 24
            complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

            /* Effective wilson coefficients */
            // cf. [BFS2001], below Eq. (9), p. 4
            complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
            // cf. [BFS2001], below Eq. (26), p. 8
            complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

            /* perpendicular, top sector */
            // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
            complex<double> C0_top_perp_left  = (c7eff - wc.c7prime()) + s / (2.0 * m_b_PS * m_B) * Y_top;
            complex<double> C0_top_perp_right = (c7eff + wc.c7prime()) + s / (2.0 * m_b_PS * m_B) * Y_top;
            // cf. [BFS2004], Eq. (44), p. 24
            complex<double> C1f_top_perp_left  = (c7eff - wc.c7prime()) * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
            complex<double> C1f_top_perp_right = (c7eff + wc.c7prime()) * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
            // cf. [BFS2001], Eqs. (34), (37), p. 9
            complex<double> C1nf_top_perp = (-1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (s / (2.0 * m_b_PS * m_B)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

            /* perpendicular, up sector */
            // cf. [BFS2004], comment before Eq. (43), p. 24
            complex<double> C0_up_perp = s / (2.0 * m_b_PS * m_B) * Y_up;
            // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
            // cf. [BFS2001], Eqs. (34), (37), p. 9
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
            complex<double> C1nf_up_perp = (-1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                    + (s / (2.0 * m_b_PS * m_B)) * (
                        wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                        + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

            /* parallel, top sector */
            // cf. [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
            complex<double> C0_top_par = -1.0 * (c7eff - wc.c7prime() + m_B / (2.0 * m_b_PS) * Y_top);
            // cf. [BFS2004], Eq. (45), p. 24
            complex<double> C1f_top_par = -1.0 * (c7eff - wc.c7prime()) * (8.0 * std::log(m_b_PS / mu) + 2.0 * L - 4.0 * (1.0 - mu_f() / m_b_PS));
            // cf. [BFS2001], Eqs. (38), p. 9
            complex<double> C1nf_top_par = (+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

            /* parallel, up sector */
            // cf. [BFS2004], comment before Eq. (43), p. 24
            complex<double> C0_up_par = -1.0 * m_B / (2.0 * m_b_PS) * Y_up;
            // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
            // cf. [BFS2004], last paragraph in Sec A.1, p. 24
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
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
            // T0_top_perp_{p,m} = 0, cf. [BFS2001], Eq. (17), p. 6
            // cf. [BFS2004], Eq. (49)
            complex<double> T1f_top_perp_p_left  = (c7eff - wc.c7prime()) * (2.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
            complex<double> T1f_top_perp_p_right = (c7eff + wc.c7prime()) * (2.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
            // T1f_top_perp_m = 0, cf. [BFS2001], Eq. (22), p. 7
            // cf. [BFS2001], Eq. (23), p. 7
            // [Christoph] Use c8 instead of c8eff
            complex<double> T1nf_top_perp_p = (-4.0 * e_d * c8eff * qcdf_0.j0bar_perp
                + m_B / (2.0 * m_b_PS) * (
                        e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde1_perp
                        + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0/3.0 * wc.c6() - (4.0 * m_b_PS / m_B) * (wc.c3() - wc.c4()/6.0 + 4.0 * wc.c5() - 2.0/3.0 * wc.c6())) * qcdf_b.jtilde1_perp
                        + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() - 8.0/3.0 * wc.c6()) * qcdf_0.jtilde1_perp)) * lambda_B_p_inv;
            // T1nf_top_perp_m = 0, cf. [BFS2001], Eq. (17), p. 6

            /* perpendicular, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eq. (50), p. 25
            complex<double> T1nf_up_perp_p = +e_u * m_B / (2.0 * m_b_PS) * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde1_perp - qcdf_0.jtilde1_perp) * lambda_B_p_inv;

            /* parallel, top sector */
            // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
            // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
            complex<double> T0_top_par_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
            // cf. [BFS2004], Eq. (49), p. 25
            complex<double> T1f_top_par_p  = (c7eff - wc.c7prime()) * (4.0 * m_B / energy) * invm1_par * lambda_B_p_inv;
            // T1f_top_par_m = 0, cf. [BFS2001], Eq. (22), p. 7
            // cf. [BFS2001], Eq. (25), p. 7
            complex<double> T1nf_top_par_p = m_B / m_b_PS * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
            // cf. [BFS2001], Eq. (26), pp. 7-8
            complex<double> T1nf_top_par_m = e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                    + 6.0 * m_B / m_b_PS * (
                        (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                        -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

            /* parallel, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
            complex<double> T0_up_par_m = +e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;
            // cf. [BFS2004], Eq. (50), p. 25
            complex<double> T1nf_up_par_p = +e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
            // cf. [BFS2004], Eq. (50), p. 25 without the \omega term
            complex<double> T1nf_up_par_m = +e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


            // Compute the nonfactorizing contributions
            complex<double> T_perp_left  = a_mu_f * (T1f_top_perp_p_left  + T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p);
            complex<double> T_perp_right = a_mu_f * (T1f_top_perp_p_right + T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p);
            complex<double> T_par = a_mu_f * (T1f_top_par_p + T1nf_top_par_p + lambda_hat_u * T1nf_up_par_p)
                + (T0_top_par_m + lambda_hat_u * T0_up_par_m + a_mu_f * (T1nf_top_par_m + lambda_hat_u * T1nf_up_par_m));

            // Compute the numerically leading power-suppressed weak annihilation contributions to order alpha_s^0
            // cf. [BFS2004], Eq. (51)
            complex<double> Delta_T_ann_top_perp = e_q * M_PI * M_PI * f_B / 3.0 / m_b_PS / m_B * (
                    -4.0 * f_Kstar_perp * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 3.0 * wc.c5() + 4.0 * wc.c6())) * qcdf_0.j0_perp
                    + 2.0 * f_Kstar_par * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 12.0 * wc.c5() + 16.0 * wc.c6())) *
                        (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p));
            complex<double> Delta_T_ann_up_perp = -e_q * 2.0 * M_PI * M_PI * f_B * f_Kstar_par / 3.0 / m_b_PS / m_B *
                (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p) * 3.0 * delta_qu * wc.c2();
            // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
            // cf. [BFS2004], Eqs. (52), (53)
            complex<double> Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                    12.0 * c8eff * (m_b_PS / m_B) * f_Kstar_perp() * 1.0 / 3.0 * (qcdf_0.j0_perp + qcdf_0.j7_perp)
                    + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (
                        (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j5_perp
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j5_perp
                        - (8.0 / 27.0) * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()) * qcdf_0.j0_perp)
                    - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (
                        (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j6_perp
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j6_perp
                        - 8.0 / 27.0 * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6())));
            complex<double> Delta_T_hsa_up_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                    + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0) * (qcdf_c.j5_perp - qcdf_0.j5_perp)
                    - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0)
                        * (qcdf_c.j6_perp - qcdf_0.j6_perp));

            // Compute the sum of the numerically leading power-suppressed contributions
            complex<double> Delta_T_top_perp = Delta_T_ann_top_perp + Delta_T_hsa_top_perp;
            complex<double> Delta_T_up_perp = Delta_T_ann_up_perp + Delta_T_hsa_up_perp;
            complex<double> Delta_T_perp = Delta_T_top_perp + lambda_hat_u * Delta_T_up_perp;


            // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
            DipoleFormFactors result;
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

        DipoleFormFactors calT_ABBBSW2008(const double & s, const WilsonCoefficients<BToS> & wc) const
        {
            // charges of down- and up-type quarks
            static const double
                e_d = -1.0 / 3.0,
                e_u = +2.0 / 3.0;

            // spectator contributions
            const double delta_qu = (q == 'u' ? 1.0 : 0.0);

            // kinematics
            const double
                m_c_pole = model->m_c_pole(),
                m_b_PS = this->m_b_PS(),
                energy = this->energy(s);

            // couplings
            const double
                alpha_s_mu = model->alpha_s(mu()), // alpha_s at the hard scale
                a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI,
                alpha_s_mu_f = model->alpha_s(std::sqrt(mu() * 0.5)), // alpha_s at the factorization scale
                a_mu_f = alpha_s_mu_f * QCD::casimir_f / 4.0 / M_PI;

            complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
            if (cp_conjugate)
                lambda_hat_u = std::conj(lambda_hat_u);

            QCDFIntegrals::Results
                qcdf_0 = QCDFIntegrals::dilepton_massless_case(s, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par),
                qcdf_c = QCDFIntegrals::dilepton_charm_case(s, m_c_pole, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par),
                qcdf_b = QCDFIntegrals::dilepton_bottom_case(s, m_b_PS, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);

            // inverse of the "negative" moment of the B meson LCDA
            // cf. [BFS2001], Eq. (54), p. 15
            const double
                omega_0 = lambda_B_p,
                lambda_B_p_inv = 1.0 / lambda_B_p;

            const complex<double>
                lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

            /* Effective wilson coefficients */
            // cf. [BFS2001], below Eq. (26), p. 8
            const complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

            /* perpendicular, top sector */
            // cf. [BFS2001], Eqs. (34), (37), p. 9
            const complex<double>
                C1nf_top_perp = (-1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (s / (2.0 * m_b_PS * m_B)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS))),

            /* perpendicular, up sector */
            // cf. [BFS2001], Eqs. (34), (37), p. 9
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
                C1nf_up_perp = (-1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                    + (s / (2.0 * m_b_PS * m_B)) * (
                        wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                        + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS)))),

            /* parallel, top sector */
            // cf. [BFS2001], Eqs. (38), p. 9
                C1nf_top_par = (+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS))),

            /* parallel, up sector */
            // cf. [BFS2004], last paragraph in Sec A.1, p. 24
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
                C1nf_up_par = (+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                        + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

            // compute the factorizing contributions
            // in ABBBSW2008: C0 is included in naively factorizing part and C1f = 0
            const complex<double>
                C_perp = a_mu * (C1nf_top_perp + lambda_hat_u * C1nf_up_perp),
                C_par  = a_mu * (C1nf_top_par  + lambda_hat_u * C1nf_up_par);


            /* perpendicular, top sector */
            // cf. [BFS2001], Eq. (23), p. 7
            const complex<double>
                T1nf_top_perp_p = (-4.0 * e_d * c8eff * qcdf_0.j0bar_perp
                    + m_B / (2.0 * m_b_PS) * (
                        e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde1_perp
                        + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0/3.0 * wc.c6() - (4.0 * m_b_PS / m_B) * (wc.c3() - wc.c4()/6.0 + 4.0 * wc.c5() - 2.0/3.0 * wc.c6())) * qcdf_b.jtilde1_perp
                        + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() - 8.0/3.0 * wc.c6()) * qcdf_0.jtilde1_perp)) * lambda_B_p_inv,
            // T1nf_top_perp_m = 0, cf. [BFS2001], Eq. (17), p. 6

            /* perpendicular, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eq. (50), p. 25
                T1nf_up_perp_p = +e_u * m_B / (2.0 * m_b_PS) * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde1_perp - qcdf_0.jtilde1_perp) * lambda_B_p_inv,

            /* parallel, top sector */
            // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
                T0_top_par_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv,
            // cf. [BFS2001], Eq. (25), p. 7
                T1nf_top_par_p = m_B / m_b_PS * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv,
            // cf. [BFS2001], Eq. (26), pp. 7-8
                T1nf_top_par_m = e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                    + 6.0 * m_B / m_b_PS * (
                        (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                        -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

            /* parallel, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
            const complex<double>
                T0_up_par_m = +e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv,
            // cf. [BFS2004], Eq. (50), p. 25
                T1nf_up_par_p = +e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv,
            // cf. [BFS2004], Eq. (50), p. 25 without the \omega term
                T1nf_up_par_m = +e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


            // Compute the nonfactorizing contributions
            const complex<double>
                T_perp = a_mu_f * (T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p),
                T_par  = a_mu_f * (T1nf_top_par_p  + lambda_hat_u * T1nf_up_par_p)
                         + T0_top_par_m + lambda_hat_u * T0_up_par_m
                         + a_mu_f * (T1nf_top_par_m + lambda_hat_u * T1nf_up_par_m);

            // Compute the numerically leading power-suppressed weak annihilation contributions to order alpha_s^0
            // cf. [BFS2004], Eq. (51)
            const complex<double>
                Delta_T_ann_top_perp = e_q * M_PI * M_PI * f_B / 3.0 / m_b_PS / m_B * (
                    -4.0 * f_Kstar_perp * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 3.0 * wc.c5() + 4.0 * wc.c6())) * qcdf_0.j0_perp
                    + 2.0 * f_Kstar_par * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 12.0 * wc.c5() + 16.0 * wc.c6())) *
                        (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p)),
                Delta_T_ann_up_perp = -e_q * 2.0 * M_PI * M_PI * f_B * f_Kstar_par / 3.0 / m_b_PS / m_B
                    * (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p) * 3.0 * delta_qu * wc.c2(),
            // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
            // cf. [BFS2004], Eqs. (52), (53)
                Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                    12.0 * c8eff * (m_b_PS / m_B) * f_Kstar_perp() * 1.0 / 3.0 * (qcdf_0.j0_perp + qcdf_0.j7_perp)
                    + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (
                          (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                        + (wc.c3() +  5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j5_perp
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j5_perp
                        - (8.0 / 27.0) * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()) * qcdf_0.j0_perp)
                    - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (
                          (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                        + (wc.c3() +  5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j6_perp
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j6_perp
                        - 8.0 / 27.0 * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))),
                Delta_T_hsa_up_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                    + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0) * (qcdf_c.j5_perp - qcdf_0.j5_perp)
                    - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0)
                        * (qcdf_c.j6_perp - qcdf_0.j6_perp));

            // Compute the sum of the numerically leading power-suppressed contributions
            const complex<double>
                Delta_T_top_perp = Delta_T_ann_top_perp + Delta_T_hsa_top_perp,
                Delta_T_up_perp  = Delta_T_ann_up_perp + Delta_T_hsa_up_perp,
                Delta_T_perp     = Delta_T_top_perp + lambda_hat_u * Delta_T_up_perp;


            // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
            DipoleFormFactors result;
            result.calT_perp_left  = xi_perp(s) * C_perp + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp + Delta_T_perp;
            result.calT_perp_right = result.calT_perp_left;
            result.calT_parallel   = xi_par(s) * C_par   + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_par * m_Kstar) / (m_B * energy) * T_par;

            return result;
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B() / (m_B() + m_Kstar());
            double result = uncertainty_xi_perp * factor * form_factors->v(s);

            return result;
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B() + m_Kstar()) / (2.0 * energy(s));
            const double factor2 = (1.0 - m_Kstar() / m_B());
            double result = uncertainty_xi_par * (factor1 * form_factors->a_1(s) - factor2 * form_factors->a_2(s));

            return result;
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double lam(const double & s) const
        {
            return lambda(m_B() * m_B(), m_Kstar() * m_Kstar(), s);
        }

        double norm(const double & s) const
        {
            double lambda_t2 = std::norm(model->ckm_tb() * conj(model->ckm_ts()));

            return g_fermi() * alpha_e() * std::sqrt(
                      1.0 / 3.0 / 1024 / power_of<5>(M_PI) / m_B()
                      * lambda_t2 * s_hat(s) * std::sqrt(lam(s)) * beta_l(s)
                   ); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(double s) const
        {
            return s / m_B() / m_B();
        }

        inline double mu_f() const
        {
            return 1.5;
        }

        inline double m_b_PS() const
        {
            // Actually use the PS mass at mu_f = 1.5 GeV
            return model->m_b_ps(mu_f());
        }

        inline double energy(const double & s) const
        {
            return (m_B() * m_B() + m_Kstar() * m_Kstar() - s) / (2.0 * m_B());
        }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        // cf. [BHvD2012], app B, eqs. (B13 - B19)
        Amplitudes amp_BFS2004(const double & s) const
        {
            Amplitudes result;

            WilsonCoefficients<BToS> wc = wilson_coefficients();

            const double
                shat = s_hat(s),
                mbhat = m_b_PS() / m_B,
                mKhat2 = power_of<2>(m_Kstar() / m_B()),
                m_K2 = power_of<2>(m_Kstar()),
                m_B2 = power_of<2>(m_B()),
                m2_diff = m_B2 - m_K2,
                norm_s = this->norm(s),
                sqrt_lam = std::sqrt(lam(s)),
                sqrt_s = std::sqrt(s);

            DipoleFormFactors dff = calT_BFS2004(s, wc);

            const complex<double>
                wilson_minus_right = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime()),
                wilson_minus_left  = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime()),
                wilson_plus_right  = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime()),
                wilson_plus_left   = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime());

            // longitudinal amplitude
            const double prefactor_long = -norm_s / (2.0 * m_Kstar() * std::sqrt(s));

            const complex<double>
                a = (m2_diff - s) * 2.0 * energy(s) * xi_perp(s) - lam(s) * m_B() / m2_diff * (xi_perp(s) - xi_par(s)),
                b = 2.0 * m_b_PS() * (
                        ((m_B2 + 3.0 * m_K2 - s) * 2.0 * energy(s) / m_B() - lam(s) / m2_diff) * dff.calT_perp_left
                        - lam(s) / m2_diff * dff.calT_parallel
                    );

            result.a_long_right = prefactor_long * (wilson_minus_right * a + uncertainty_long() * b);
            result.a_long_left  = prefactor_long * (wilson_minus_left  * a + uncertainty_long() * b);

            // perpendicular amplitude
            const double prefactor_perp = +std::sqrt(2.0) * norm_s * m_B() * std::sqrt(lambda(1.0, mKhat2, shat));

            result.a_perp_right = prefactor_perp * (wilson_plus_right * xi_perp(s) + uncertainty_perp() * (2.0 * mbhat / shat) * dff.calT_perp_right);
            result.a_perp_left  = prefactor_perp * (wilson_plus_left  * xi_perp(s) + uncertainty_perp() * (2.0 * mbhat / shat) * dff.calT_perp_right);

            // parallel amplitude
            const double prefactor_par = -std::sqrt(2.0) * norm_s * m2_diff;

            result.a_par_right = prefactor_par * (
                                    wilson_minus_right * xi_perp(s) * 2.0 * energy(s) / m2_diff
                                    + uncertainty_para() * 4.0 * m_b_PS() * energy(s) / s / m_B() * dff.calT_perp_left
                                 );
            result.a_par_left  = prefactor_par * (
                                    wilson_minus_left  * xi_perp(s) * 2.0 * energy(s) / m2_diff
                                    + uncertainty_para() * 4.0 * m_b_PS() * energy(s) / s / m_B() * dff.calT_perp_left
                                 );

            // timelike amplitude
            result.a_timelike = norm_s * sqrt_lam / sqrt_s
                * (2.0 * (wc.c10() - wc.c10prime()) + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime()))
                * form_factors->a_0(s);

            // scalar amplitude
            result.a_scalar = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar) * form_factors->a_0(s);

            // tensor amplitudes [BHvD2012]  eqs. (B18 - B20)
            // no form factor relations used
            const double
                ff_T1  = form_factors->t_1(s),
                ff_T2  = form_factors->t_2(s),
                ff_T3  = form_factors->t_3(s),

                kin_tensor_1 = norm_s / m_Kstar() * ((m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam(s) / m2_diff * ff_T3),
                kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1,
                kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

            // correct the sign of C_T5 from [BHvD2012v4] because of inconsistent use of gamma5 <-> Levi-Civita
            static const double sign = -1;

            result.a_par_perp = kin_tensor_1 * wc.cT();
            result.a_t_long   = kin_tensor_1 * sign * wc.cT5();

            result.a_t_perp    = kin_tensor_2 * wc.cT();
            result.a_long_perp = kin_tensor_2 * sign * wc.cT5();

            result.a_t_par     = kin_tensor_3 * sign * wc.cT5();
            result.a_long_par  = kin_tensor_3 * wc.cT();

            return result;
        }

        // cf. [ABBBSW2008]
        // cf. [BHvD2012] for tensor amplitudes
        // use full QCD form factors in leading QCDF (naively factorizing) amplitudes
        // use soft form factors in non-factorizable contributions (~ alpha_s)
        Amplitudes amp_ABBBSW2008(const double & s) const
        {
            Amplitudes result;

            WilsonCoefficients<BToS> wc = wilson_coefficients();

            const double
                shat = s_hat(s),
                sqrt_s = std::sqrt(s),
                mKhat2 = power_of<2>(m_Kstar() / m_B()),
                m_K2 = power_of<2>(m_Kstar()),
                m_B2 = power_of<2>(m_B()),
                m_sum = m_B() + m_Kstar(),
                m_diff = m_B() - m_Kstar(),
                m2_diff = m_B2 - m_K2,
                norm_s = this->norm(s),
                sqrt_lam = std::sqrt(lam(s));

            const double
                ff_V   = form_factors->v(s),
                ff_A0  = form_factors->a_0(s),
                ff_A1  = form_factors->a_1(s),
                ff_A2  = form_factors->a_2(s),
                ff_T1  = form_factors->t_1(s),
                ff_T2  = form_factors->t_2(s),
                ff_T3  = form_factors->t_3(s);

            const double
                m_c_pole = model->m_c_pole(),
                m_b_PS = this->m_b_PS();

            complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
            if (cp_conjugate)
                lambda_hat_u = std::conj(lambda_hat_u);

            /* Y(s) for the up and the top sector for effective Wilson coefficients */
            // cf. [BFS2001], Eq. (10), p. 4
            const complex<double>
                Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5(),
                Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6()),
                Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6()),
                Y_top_  = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

            // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
            // then replace b pole mass by the PS mass.
            complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
                 + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
                 + Y_top_0 * CharmLoops::h(mu, s)
                 + Y_top_;
            // cf. [BFS2004], Eq. (43), p. 24
            complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

            const complex<double>
                // cf. [BFS2001], below Eq. (9), p. 4
                c7eff = - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6(),
                c9eff = Y_top + lambda_hat_u * Y_up,
                c910_mi_r = (wc.c9() + c9eff - wc.c9prime()) + (wc.c10() - wc.c10prime()),
                c910_mi_l = (wc.c9() + c9eff - wc.c9prime()) - (wc.c10() - wc.c10prime()),
                c910_pl_r = (wc.c9() + c9eff + wc.c9prime()) + (wc.c10() + wc.c10prime()),
                c910_pl_l = (wc.c9() + c9eff + wc.c9prime()) - (wc.c10() + wc.c10prime()),
                // here b-quark mass in O_7 is MS-bar contrary to PS-scheme in [BFS2001]
                c7_mi = 2.0 * m_b_MSbar * (wc.c7() + c7eff - wc.c7prime()),
                c7_pl = 2.0 * m_b_MSbar * (wc.c7() + c7eff + wc.c7prime());

            //
            // Naive factorization part - including 1-loop effective wilson coefficients c_{7,9}^eff
            //

            // longitudinal amplitude
            const double pre_long = -norm_s / (2.0 * m_Kstar() * std::sqrt(s));

            const complex<double>
                a = (m2_diff - s) * m_sum * ff_A1 - lam(s) / m_sum * ff_A2,
                b = (m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam(s) / m2_diff * ff_T3;

            result.a_long_right = pre_long * (c910_mi_r * a + c7_mi * b);
            result.a_long_left  = pre_long * (c910_mi_l * a + c7_mi * b);

            // perpendicular amplitude
            const double pre_perp = +std::sqrt(2.0) * norm_s* m_B2 * std::sqrt(lambda(1.0, mKhat2, shat));

            result.a_perp_right = pre_perp * (c910_pl_r * ff_V / m_sum + c7_pl / s * ff_T1);
            result.a_perp_left  = pre_perp * (c910_pl_l * ff_V / m_sum + c7_pl / s * ff_T1);

            // parallel amplitude
            const double pre_par = -std::sqrt(2.0) * norm_s * m2_diff;

            result.a_par_right = pre_par * (c910_mi_r * ff_A1 / m_diff + c7_mi / s * ff_T2);
            result.a_par_left  = pre_par * (c910_mi_l * ff_A1 / m_diff + c7_mi / s * ff_T2);

            // timelike amplitude
            result.a_timelike = norm_s * sqrt_lam / sqrt_s
                * (2.0 * (wc.c10() - wc.c10prime())
                   + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime()))
                * ff_A0;

            // scalar amplitude
            result.a_scalar = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar) * ff_A0;

            // tensor amplitudes
            const double
                kin_tensor_1 = norm_s / m_Kstar() * ((m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam(s) / m2_diff * ff_T3),
                kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1,
                kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

            // correct the sign of C_T5 from [BHvD2012v4] because of inconsistent use of gamma5 <-> Levi-Civita
            static const double sign = -1;

            result.a_par_perp = kin_tensor_1 * wc.cT();
            result.a_t_long   = kin_tensor_1 * sign * wc.cT5();

            result.a_t_perp    = kin_tensor_2 * wc.cT();
            result.a_long_perp = kin_tensor_2 * sign * wc.cT5();

            result.a_t_par     = kin_tensor_3 * sign * wc.cT5();
            result.a_long_par  = kin_tensor_3 * wc.cT();

            //
            // Beyond Naive factorization part - from QCDF
            //

            DipoleFormFactors dff = calT_ABBBSW2008(s, wc);

            // these kinematical factors reduce for mKstar = 0 to [ABBBSW2008] eq. (3.46)
#if 0
            const complex<double>
                c_long = 2.0 * m_b_MSbar * (
                           ((m_B2 + 3.0 * m_K2 - s) * 2.0 * energy(s) / m_B() - lam(s) / m2_diff) * dff.calT_perp_left
                            - lam(s) / m2_diff * dff.calT_parallel
                         );

            result.a_long_right += uncertainty_long_right * pre_long * c_long;
            result.a_long_left  += uncertainty_long_left  * pre_long * c_long;

            result.a_perp_right += uncertainty_perp_right * pre_perp * 2.0 * m_b_MSbar * m_B() / s * dff.calT_perp_right;
            result.a_perp_left  += uncertainty_perp_left  * pre_perp * 2.0 * m_b_MSbar * m_B() / s * dff.calT_perp_right;

            result.a_par_right  += uncertainty_par_right  * pre_par  * 4.0 * m_b_MSbar * energy(s) / s / m_B() * dff.calT_perp_left;
            result.a_par_left   += uncertainty_par_left   * pre_par  * 4.0 * m_b_MSbar * energy(s) / s / m_B() * dff.calT_perp_left;
#endif

            // using approximation mKstar = 0 as in [ABBBSW2008] eq. (3.46)
            // [christoph] !!! setting mKstar = 0 causes quite large shifts in Br (+17%), F_L (-17%) (integrated from 1.1 to 6 GeV^2, KMPW-form factors)
            // did not test other observables
            const double
                pre_p = std::sqrt(2.0) * norm_s * 2.0 * m_b_MSbar / s * (m_B2 - s),
                pre_0 = norm_s / m_Kstar() / m_B2 / sqrt_s * m_b_MSbar * power_of<2>(m_B2 - s);

            result.a_long_right += uncertainty_long() * pre_0 * dff.calT_parallel;
            result.a_long_left  += uncertainty_long() * pre_0 * dff.calT_parallel;

            result.a_perp_right += uncertainty_perp() * pre_p * dff.calT_perp_right;
            result.a_perp_left  += uncertainty_perp() * pre_p * dff.calT_perp_right;

            result.a_par_right  -= uncertainty_para() * pre_p * dff.calT_perp_left;
            result.a_par_left   -= uncertainty_para() * pre_p * dff.calT_perp_left;

            return result;
        }

        Amplitudes amplitudes(const double & s) const
        {
            Amplitudes amp;

            if (ff_relation == "BFS2004")
                amp = amp_BFS2004(s);
            else if (ff_relation == "ABBBSW2008")
                amp = amp_ABBBSW2008(s);
            else
                throw InvalidOptionValueError("large-recoil-ff", ff_relation, "BFS2004, ABBBSW2008");
            return amp;
        }

        std::array<double, 12> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitudes(s), s, m_l());
        }

        AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return array_to_angular_coefficients(angular_coefficients_array(amplitudes(s), s, m_l()));
        }

        AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand =
                    std::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 12> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return array_to_angular_coefficients(integrated_angular_coefficients_array);
        }

        double a_fb_zero_crossing() const
        {
            // We trust QCDF results in a validity range from 0.5 GeV^2 < s < 6.0 GeV^2
            static const double min_result = 0.5;
            static const double max_result = 7.0;

            // use calT_perp / xi_perp = C_7 as start point
            WilsonCoefficients<BToS> wc = wilson_coefficients();
            const double start = -2.0 * model->m_b_msbar(mu()) * m_B() * real(wc.c7() / wc.c9());

            double result = start;
            // clamp result to QCDF validity region
            result = std::max(min_result, result);
            result = std::min(max_result, result);

            // perform a couple of Newton-Raphson steps
            for (unsigned i = 0 ; i < 100 ; ++i)
            {
                double xplus = result * 1.03;
                double xminus = result * 0.97;

                AngularCoefficients a_c_central = differential_angular_coefficients(result);
                double f = a_c_central.j6s + 0.5 * a_c_central.j6c;
                AngularCoefficients a_c_minus   = differential_angular_coefficients(xminus);
                double f_xminus = a_c_minus.j6s + 0.5 * a_c_minus.j6c;
                AngularCoefficients a_c_plus    = differential_angular_coefficients(xplus);
                double f_xplus = a_c_plus.j6s + 0.5 * a_c_plus.j6c;

                double fprime = (f_xplus - f_xminus) / (xplus - xminus);

                if (std::abs(f / fprime) < 1e-8)
                    break;

                result = result - f / fprime;
                // clamp result to QCDF validity region
                result = std::max(min_result, result);
                result = std::min(max_result, result);
            }

            return result;
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options, *this))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    double
    BToKstarDilepton<LargeRecoil>::xi_perp(const double & s) const
    {
        return _imp->xi_perp(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::xi_para(const double & s) const
    {
        return _imp->xi_par(s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        if (h == -1)
            return amp.a_long_left;
        else
            return amp.a_long_right;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        if (h == -1)
            return amp.a_perp_left;
        else
            return amp.a_perp_right;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_par(const Helicity & h, const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        if (h == -1)
            return amp.a_par_left;
        else
            return amp.a_par_right;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_timelike(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_timelike;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_scalar(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_scalar;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_par_perp(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_par_perp;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_t_long(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_t_long;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_t_perp(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_t_perp;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_t_par(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_t_par;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long_par(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_long_par;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long_perp(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return amp.a_long_perp;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_decay_width(const double & s) const
    {
        return decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_isospin_asymmetry(const double & s) const
    {
        Save<char> save_q(_imp->q, 'd');
        Save<double> save_e_q(_imp->e_q, -1.0/3.0);

        double gamma_zero = differential_decay_width(s);
        _imp->q = 'u';
        _imp->e_q = +2.0/3.0;
        double gamma_minus = differential_decay_width(s);

        return (gamma_zero - gamma_minus) / (gamma_zero + gamma_minus);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.8)
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j6s + 0.5 * a_c.j6c) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.11)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.12)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((power_of<2>(_imp->beta_l(s) * a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_re(const double & s) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.25 * _imp->beta_l(s) * a_c.j6s / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_im(const double & s) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_p_prime_4(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (15)
        return (a_c.j4 + a_c_bar.j4) / std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_p_prime_5(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (16)
        return (a_c.j5 + a_c_bar.j5) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_p_prime_6(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (17)
        return -1.0 * (a_c.j7 + a_c_bar.j7) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], eq. (A9)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j1c - a_c.j2c / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transversal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], eq. (A10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_1(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_2(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_3(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    // differential angular coefficients
    double
    BToKstarDilepton<LargeRecoil>::differential_j_1c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1c;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_1s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_2c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2c;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_2s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_3(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j3;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_3_normalized(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j3 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_3_normalized_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return (a_c.j3 + a_c_bar.j3) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j4;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j5;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_6c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6c;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_6s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_6c_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return 0.5 * (a_c.j6c + a_c_bar.j6c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_7(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j7;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_8(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j8;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_9(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j9;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_9_normalized(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j9 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_9_normalized_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return (a_c.j9 + a_c_bar.j9) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_1c_plus_j_2c_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return 0.5 * (a_c.j1c + a_c_bar.j1c + a_c.j2c + a_c_bar.j2c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_1s_minus_3j_2s_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return 0.5 * (a_c.j1s + a_c_bar.j1s - 3.0 * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_d_4(const double & s) const
    {
        double J4_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            J4_electrons = differential_j_4(s);
        }

        double J4_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            J4_muons = differential_j_4(s);
        }

        return 4.0 / 3.0 * (J4_electrons - J4_muons) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_d_5(const double & s) const
    {
        double J5_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            J5_electrons = differential_j_5(s);
        }

        double J5_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            J5_muons = differential_j_5(s);
        }

        return 3.0 / 4.0 * (J5_electrons - J5_muons) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_d_6s(const double & s) const
    {
        double J6s_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            J6s_electrons = differential_j_6s(s);
        }

        double J6s_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            J6s_muons = differential_j_6s(s);
        }

        return 3.0 / 4.0 * (J6s_electrons - J6s_muons) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_ratio_muons_electrons(const double & s) const
    {
        double gamma_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            gamma_electrons = differential_decay_width(s);
        }

        double gamma_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            gamma_muons = differential_decay_width(s);
        }

        return gamma_muons / gamma_electrons;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrated_decay_width(s_min, s_max) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double br = integrated_branching_ratio(s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrated_branching_ratio(s_min, s_max);

        return 0.5 * (br + br_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double gamma = integrated_decay_width(s_min, s_max);
        _imp->cp_conjugate = true;
        double gamma_bar = integrated_decay_width(s_min, s_max);

        return (gamma - gamma_bar) / (gamma + gamma_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_isospin_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<char> save_q(_imp->q, 'd');
        Save<double> save_e_q(_imp->e_q, -1.0/3.0);

        double gamma_zero = integrated_decay_width(s_min, s_max);
        _imp->q = 'u';
        _imp->e_q = +2.0/3.0;
        double gamma_minus = integrated_decay_width(s_min, s_max);

        return (gamma_zero - gamma_minus) / (gamma_zero + gamma_minus);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.8), p. 6
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j6s + 0.5 * a_c.j6c) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double a_fb = integrated_forward_backward_asymmetry(s_min, s_max);
        _imp->cp_conjugate = true;
        double a_fb_bar = integrated_forward_backward_asymmetry(s_min, s_max);

        return 0.5 * (a_fb + a_fb_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], eq. (A9)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j1c - a_c.j2c / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double f_l = integrated_longitudinal_polarisation(s_min, s_max);
        _imp->cp_conjugate = true;
        double f_l_bar = integrated_longitudinal_polarisation(s_min, s_max);

        return 0.5 * (f_l + f_l_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transversal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], eq. (A10)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transversal_polarisation_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double f_t = integrated_transversal_polarisation(s_min, s_max);
        _imp->cp_conjugate = true;
        double f_t_bar = integrated_transversal_polarisation(s_min, s_max);

        return 0.5 * (f_t + f_t_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2_cp_averaged(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        Save<bool> save(_imp->cp_conjugate, false);

        double a_t_2 = integrated_transverse_asymmetry_2(s_min, s_max);
        _imp->cp_conjugate = true;
        double a_t_2_bar = integrated_transverse_asymmetry_2(s_min, s_max);

        return 0.5 * (a_t_2 + a_t_2_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.11), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.12), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((power_of<2>(a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.25 * a_c.j6s / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_p_prime_4(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (15)
        return (a_c.j4 + a_c_bar.j4) / std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_p_prime_5(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (16)
        return (a_c.j5 + a_c_bar.j5) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_p_prime_6(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (17)
        return -1.0 * (a_c.j7 + a_c_bar.j7) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_1(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return  a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::a_fb_zero_crossing() const
    {
        return _imp->a_fb_zero_crossing();
    }

    // integrated angular coefficients
    double
    BToKstarDilepton<LargeRecoil>::integrated_j_1c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1c;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_1s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_2c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2c;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_2s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_3(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_3_normalized(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_3_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j3 + a_c_bar.j3) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j4;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_4_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j4 + a_c_bar.j4) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j5;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_5_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j5 + a_c_bar.j5) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_6c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6c;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_6s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_7(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j7;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_7_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j7 + a_c_bar.j7) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_8(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j8;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_8_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j8 + a_c_bar.j8) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_9(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j9;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_9_normalized(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j9 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_9_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j9 + a_c_bar.j9) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_a_9(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j9 - a_c_bar.j9) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_d_4(const double & s_min, const double & s_max) const
    {
        double J4_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            J4_electrons = integrated_j_4(s_min, s_max);
        }

        double J4_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            J4_muons = integrated_j_4(s_min, s_max);
        }

        return 3.0 / 4.0 * (J4_electrons - J4_muons) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_d_5(const double & s_min, const double & s_max) const
    {
        double J5_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            J5_electrons = integrated_j_5(s_min, s_max);
        }

        double J5_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            J5_muons = integrated_j_5(s_min, s_max);
        }

        return 3.0 / 4.0 * (J5_electrons - J5_muons) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_d_6s(const double & s_min, const double & s_max) const
    {
        double J6s_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");

            J6s_electrons = integrated_j_6s(s_min, s_max);
        }

        double J6s_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");

            J6s_muons = integrated_j_6s(s_min, s_max);
        }

        return 3.0 / 4.0 * (J6s_electrons - J6s_muons) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");
            br_electrons = integrate<GSL::QNG>(integrand, s_min, s_max);
        }

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");
            br_muons = integrate<GSL::QNG>(integrand, s_min, s_max);
        }

        return br_muons / br_electrons;
    }

    double
    BToKstarDilepton<LargeRecoil>::four_differential_decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const
    {
        // compute d^4 Gamma, cf. [BHvD2010], p. 5, Eq. (2.6)
        // Cosine squared of the angles
        double c_theta_k_2 = c_theta_k * c_theta_k;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = cos(phi);
        // Sine squared of the angles
        double s_theta_k_2 = 1.0 - c_theta_k_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_k = sqrt(s_theta_k_2);
        double s_theta_l = sqrt(s_theta_l_2);
        double s_phi = sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_k = 2.0 * c_theta_k_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_k = 2.0 * s_theta_k * c_theta_k;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = sin(2.0 * phi);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        double Gamma = decay_width(_imp->integrated_angular_coefficients(1.00, 6.00));

        double result = 3.0 / 8.0 / M_PI * (
                 a_c.j1s + (a_c.j1c - a_c.j1s) * c_theta_k_2
                +  (a_c.j2s + (a_c.j2c - a_c.j2s) * c_theta_k_2) * c_2_theta_l
                +  a_c.j3 * s_theta_k_2 * s_theta_l_2 * c_2_phi
                +  a_c.j4 * s_2_theta_k * s_2_theta_l * c_phi
                +  a_c.j5 * s_2_theta_k * s_theta_l * c_phi
                +  (a_c.j6s * s_theta_k_2 + a_c.j6c * c_theta_k_2) * c_theta_l
                +  a_c.j7 * s_2_theta_k * s_theta_l * s_phi
                +  a_c.j8 * s_2_theta_k * s_2_theta_l * s_phi
                +  a_c.j9 * s_theta_k_2 * s_theta_l_2 * s_2_phi
                ) / Gamma;

        return result;
    }

    const std::string
    BToKstarDilepton<LargeRecoil>::description = "\
The decay Bbar->Kbar^*(-> Kbar pi) l^+ l^- in the region q^2 <= 6-8 GeV^2, with l=e,mu,tau\
a charged lepton.";

    const std::string
    BToKstarDilepton<LargeRecoil>::kinematics_description_s = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    BToKstarDilepton<LargeRecoil>::kinematics_description_c_theta_l = "\
The cosine of the negatively-charged lepton l^-'s helicity angle theta_l in the l^+l^- rest frame.";

    const std::string
    BToKstarDilepton<LargeRecoil>::kinematics_description_c_theta_k = "\
The cosine of the Kbar's helicity angle theta_k in the Kbar-pi rest frame.";

    const std::string
    BToKstarDilepton<LargeRecoil>::kinematics_description_phi = "\
The azimuthal angle between the Kbar-pi plane and the l^+l^- plane.";

    /*
     * Decay: B -> K l lbar at Large Recoil
     */
    template <>
    struct Implementation<BToKDilepton<LargeRecoil>>
    {
        Parameters parameters;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter m_b_MSbar;

        UsedParameter m_c;

        UsedParameter m_s_MSbar;

        UsedParameter m_B;

        UsedParameter m_K;

        UsedParameter m_l;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter f_B;

        UsedParameter f_K;

        UsedParameter lambda_B_p;

        UsedParameter a_1;

        UsedParameter a_2;

        // Mean life times
        UsedParameter tau;

        // Estimation of subleading contributions
        UsedParameter lambda_psd;

        UsedParameter sl_phase_psd;

        // spectator quark charge
        double e_q;

        // spectator quark flavor
        char q;

        std::string lepton_flavour;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            parameters(p),
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            m_c(p["mass::c"], u),
            m_s_MSbar(p["mass::s(2GeV)"], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            m_K(p["mass::K_" + o.get("q", "d")], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            f_B(p["decay-constant::B_" + o.get("q", "d")], u),
            f_K(p["decay-constant::K_" + o.get("q", "d")], u),
            lambda_B_p(p["lambda_B_p"], u),
            a_1(p["B->K::a_1@1GeV"], u),
            a_2(p["B->K::a_2@1GeV"], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            lambda_psd(p["B->Pll::Lambda_pseudo@LargeRecoil"], u),
            sl_phase_psd(p["B->Pll::sl_phase_pseudo@LargeRecoil"], u),
            e_q(-1.0/3.0),
            lepton_flavour(o.get("l", "mu")),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<PToP>::create("B->K@" + o.get("form-factors", "KMPW2010"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

            std::string spectator_quark = o.get("q", "d");
            if (spectator_quark.size() != 1)
                throw InternalError("Option q should only be one character!");

            q = spectator_quark[0];
            if (q == 'd')
            {
                e_q = -1.0 / 3.0;
            }
            else if (q == 'u')
            {
                e_q = 2.0 / 3.0;
            }
            else
                throw InternalError("Unsupported spectator quark");
        }

        WilsonCoefficients<BToS> wilson_coefficients() const
        {
            return model->wilson_coefficients_b_to_s(lepton_flavour, cp_conjugate);
        }

        complex<double> calT(const double & s) const
        {
            // charges of down- and up-type quarks
            static const double e_d = -1.0 / 3.0;
            static const double e_u = +2.0 / 3.0;

            // spectator contributions
            double delta_qu = (q == 'u' ? 1.0 : 0.0);

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
            WilsonCoefficients<BToS> wc = wilson_coefficients();

            // Compute the QCDF Integrals
            double invm1_psd = 3.0 * (1.0 + a_1 + a_2); // <ubar^-1>
            QCDFIntegrals::Results qcdf_0 = QCDFIntegrals::dilepton_massless_case(s, m_B, m_K, mu, 0.0, 0.0, a_1, a_2);
            QCDFIntegrals::Results qcdf_c = QCDFIntegrals::dilepton_charm_case(s, m_c_pole, m_B, m_K, mu, 0.0, 0.0, a_1, a_2);
            QCDFIntegrals::Results qcdf_b = QCDFIntegrals::dilepton_bottom_case(s, m_b_PS, m_B, m_K, mu, 0.0, 0.0, a_1, a_2);

            // inverse of the "negative" moment of the B meson LCDA
            // cf. [BFS2001], Eq. (54), p. 15
            double omega_0 = lambda_B_p, lambda_B_p_inv = 1.0 / lambda_B_p;
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
            complex<double> result = xi_pseudo(s) * C_psd + power_of<2>(M_PI) / 3.0 * (f_B * f_K) / m_B  * T_psd;

            return result;
        }

        /* Form factors */
        // cf. [BF2001], Eq. (22)
        double xi_pseudo(const double & s) const
        {
            return form_factors->f_p(s);
        }

        inline double mu_f() const
        {
            return 1.5;
        }

        inline double m_b_PS() const
        {
            // Actually use the PS mass at mu_f = 1.5 GeV
            return model->m_b_ps(mu_f());
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double lam(const double & s) const
        {
            return lambda(m_B() * m_B(), m_K() * m_K(), s);
        }

        double energy(const double & s) const
        {
            return (m_B() * m_B() + m_K() * m_K() - s) / (2.0 * m_B());
        }

        // cf. [BF2001] Eq. (22 + TODO: 31)
        double f_t_over_f_p(const double & s) const
        {
            return form_factors->f_t(s) / form_factors->f_p(s);
        }

        // cf. [BF2001] Eq. (22 + TODO: 30)
        double f_0_over_f_p(const double & s) const
        {
            return form_factors->f_0(s) / form_factors->f_p(s);
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_A(const WilsonCoefficients<BToS> & wc, const double &) const
        {
            return wc.c10() + wc.c10prime();
        }

        double F_Tkin(const double & s) const
        {
            double result = 2.0 * std::sqrt(lam(s)) * beta_l(s) / (m_B() + m_K());
            result *= f_t_over_f_p(s);
            return result;
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_T(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Tkin(s) * wc.cT();
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_T5(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Tkin(s) * wc.cT5();
        }

        double F_Skin(const double & s) const
        {
            double result = 0.5 * (power_of<2>(m_B()) - power_of<2>(m_K())) / (m_b_MSbar - m_s_MSbar);
            result *= f_0_over_f_p(s);
            return result;
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_S(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Skin(s) * (wc.cS() + wc.cSprime());
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_P(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Skin(s) * (wc.cP() + wc.cPprime()) + m_l() * (wc.c10() + wc.c10prime()) *
                    ((m_B() * m_B() - m_K() * m_K()) / s * (f_0_over_f_p(s) - 1.0) - 1.0);
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_V(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            std::complex<double> result = wc.c9() + wc.c9prime();
            result += 2.0 * m_b_PS() / m_B() / xi_pseudo(s) *
                      (calT(s) + lambda_psd / m_B * std::polar(1.0, sl_phase_psd()));
            result += 8.0 * m_l / (m_B() + m_K()) * f_t_over_f_p(s) * wc.cT();
            return result;
        }

        // cf. [BHP2007], Eqs. (4.2), (4.4), (4.5), p. 5
        double N(const double & s) const
        {
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

            return power_of<2>(g_fermi * alpha_e() * lambda_t) * std::sqrt(lam(s)) * beta_l(s) * xi_pseudo(s) * xi_pseudo(s) /
                    (512.0 * power_of<5>(M_PI) * power_of<3>(m_B()));
        }

        // cf. [BHP2007], Eq. (4.2)
        double a_l(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double result = s * (power_of<2>(beta_l(s)) * std::norm(F_S(wc, s)) + std::norm(F_P(wc, s)));
            result += 0.25 * lam(s) * (std::norm(F_A(wc, s)) + std::norm(F_V(wc, s)));
            result += 2.0 * m_l * (m_B() * m_B() - m_K() * m_K() + s) * std::real(F_P(wc, s) * std::conj(F_A(wc, s)));
            result += 4.0 * m_l * m_l * m_B() * m_B() * std::norm(F_A(wc, s));

            return N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.3)
        double b_l(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double result = s * (power_of<2>(beta_l(s)) * std::real(F_S(wc, s) * std::conj(F_T(wc, s)))
                                 + std::real(F_P(wc, s) * std::conj(F_T5(wc, s))));
            result += m_l * (std::sqrt(lam(s)) * beta_l(s) * std::real(F_S(wc, s) * std::conj(F_V(wc, s)))
                             + (m_B() * m_B() - m_K() * m_K() + s) * std::real(F_T5(wc, s) * std::conj(F_A(wc, s))));

            return 2.0 * N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.4)
        double c_l(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double result = s * (power_of<2>(beta_l(s)) * std::norm(F_T(wc, s)) + std::norm(F_T5(wc, s)));
            result -= 0.25 * lam(s) * power_of<2>(beta_l(s)) * (std::norm(F_A(wc, s)) + std::norm(F_V(wc, s)));
            result += 2.0 * m_l * std::sqrt(lam(s)) * beta_l(s) * std::real(F_T(wc, s) * std::conj(F_V(wc, s)));
            return N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.8)
        double unnormalized_decay_width(const double & s) const
        {
            WilsonCoefficients<BToS> wc = wilson_coefficients();

            return 2.0 * (a_l(wc, s) + c_l(wc, s) / 3.0);
        }

        double differential_branching_ratio(const double & s) const
        {
            return unnormalized_decay_width(s) * tau() / hbar();
        }

        // cf. [BHP2007], Eq. (4.9)
        double differential_flat_term_numerator(const double & s) const
        {
            WilsonCoefficients<BToS> wc = wilson_coefficients();

            return 2.0 * (a_l(wc, s) + c_l(wc, s));
        }

        double differential_forward_backward_asymmetry_numerator(const double & s) const
        {
            WilsonCoefficients<BToS> wc = wilson_coefficients();

            return b_l(wc, s);
        }
    };

    BToKDilepton<LargeRecoil>::BToKDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKDilepton<LargeRecoil>>(new Implementation<BToKDilepton<LargeRecoil>>(parameters, options, *this))
    {
    }

    BToKDilepton<LargeRecoil>::~BToKDilepton()
    {
    }

    std::complex<double>
    BToKDilepton<LargeRecoil>::F_A(const double & s) const
    {
        return _imp->F_A(_imp->wilson_coefficients(), s);
    }

    std::complex<double>
    BToKDilepton<LargeRecoil>::F_V(const double & s) const
    {
        return _imp->F_V(_imp->wilson_coefficients(), s);
    }

    std::complex<double>
    BToKDilepton<LargeRecoil>::F_S(const double & s) const
    {
        return _imp->F_S(_imp->wilson_coefficients(), s);
    }

    std::complex<double>
    BToKDilepton<LargeRecoil>::F_P(const double & s) const
    {
        return _imp->F_P(_imp->wilson_coefficients(), s);
    }

    std::complex<double>
    BToKDilepton<LargeRecoil>::F_T(const double & s) const
    {
        return _imp->F_T(_imp->wilson_coefficients(), s);
    }

    std::complex<double>
    BToKDilepton<LargeRecoil>::F_T5(const double & s) const
    {
        return _imp->F_T5(_imp->wilson_coefficients(), s);
    }

    double
    BToKDilepton<LargeRecoil>::a_l(const double & s) const
    {
        return _imp->a_l(_imp->wilson_coefficients(), s);
    }

    double
    BToKDilepton<LargeRecoil>::b_l(const double & s) const
    {
        return _imp->b_l(_imp->wilson_coefficients(), s);
    }

    double
    BToKDilepton<LargeRecoil>::c_l(const double & s) const
    {
        return _imp->c_l(_imp->wilson_coefficients(), s);
    }

    double
    BToKDilepton<LargeRecoil>::two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        auto wc = _imp->wilson_coefficients();

        // cf. [BHP2007], Eq. (4.1)
        return _imp->a_l(wc, s) + _imp->b_l(wc, s) * c_theta_l + _imp->c_l(wc, s) * c_theta_l * c_theta_l;
    }

    double
    BToKDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToKDilepton<LargeRecoil>::differential_flat_term(const double & s) const
    {
        return _imp->differential_flat_term_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToKDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return _imp->differential_forward_backward_asymmetry_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToKDilepton<LargeRecoil>::differential_ratio_muons_electrons(const double & s) const
    {
        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");
            br_electrons = BToKDilepton<LargeRecoil>::differential_branching_ratio(s);
        }

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");
            br_muons = BToKDilepton<LargeRecoil>::differential_branching_ratio(s);
        }

        return br_muons / br_electrons;
    }

    // Integrated Observables
    double
    BToKDilepton<LargeRecoil>::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        return integrate<GSL::QNG>(f, s_min, s_max);
    }

    double
    BToKDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(std::mem_fn(&BToKDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        return integrate<GSL::QNG>(f, s_min, s_max);
    }

    double
    BToKDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> f = std::bind(&BToKDilepton<LargeRecoil>::differential_branching_ratio,
                this, std::placeholders::_1);

        double br = integrate<GSL::QNG>(f, s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrate<GSL::QNG>(f, s_min, s_max);

        return (br + br_bar) / 2.0;
    }

    double
    BToKDilepton<LargeRecoil>::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> f = std::bind(&BToKDilepton<LargeRecoil>::differential_branching_ratio,
                this, std::placeholders::_1);

        double br = integrate<GSL::QNG>(f, s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrate<GSL::QNG>(f, s_min, s_max);

        return (br - br_bar) / (br + br_bar);
    }

    double
    BToKDilepton<LargeRecoil>::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate<GSL::QNG>(num, s_min, s_max);
        double denom_integrated = integrate<GSL::QNG>(denom, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToKDilepton<LargeRecoil>::integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate<GSL::QNG>(num, s_min, s_max);
        double denom_integrated = integrate<GSL::QNG>(denom, s_min, s_max);

        _imp->cp_conjugate = true;

        num_integrated += integrate<GSL::QNG>(num, s_min, s_max);
        denom_integrated += integrate<GSL::QNG>(denom, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    // todo caching of denominator?
    double
    BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::differential_forward_backward_asymmetry_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate<GSL::QNG>(num, s_min, s_max);
        double denom_integrated = integrate<GSL::QNG>(denom, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::differential_forward_backward_asymmetry_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate<GSL::QNG>(num, s_min, s_max);
        double denom_integrated = integrate<GSL::QNG>(denom, s_min, s_max);

        _imp->cp_conjugate = true;

        num_integrated += integrate<GSL::QNG>(num, s_min, s_max);
        denom_integrated += integrate<GSL::QNG>(denom, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToKDilepton<LargeRecoil>::integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToKDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "e");
            // br_electrons = integrate<GSL::QNG>(integrand, s_min, s_max);
            br_electrons = integrate<GSL::QNG>(integrand, s_min, s_max);
        }

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_lepton_flavour(_imp->lepton_flavour, "mu");
            br_muons = integrate<GSL::QNG>(integrand, s_min, s_max);
        }

        // cf. [BHP2007], Eq. (4.10), p. 6
        return br_muons / br_electrons;
    }

    const std::string
    BToKDilepton<LargeRecoil>::description = "\
The decay B->K l^+ l^- in the region q^2 <= 6-8 GeV^2, with l=e,mu,tau\
a charged lepton.";

    const std::string
    BToKDilepton<LargeRecoil>::kinematics_description_s = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    BToKDilepton<LargeRecoil>::kinematics_description_c_theta_l = "\
The cosine of the negatively-charged lepton l^-'s helicity angle theta_l in the l^+l^- rest frame.";
}
