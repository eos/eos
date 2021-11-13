/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2020 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Frederik Beaujean
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

#include <eos/rare-b-decays/b-to-k-ll-gvdv2020.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/power_of.hh>

#include <gsl/gsl_sf.h>

namespace eos
{
    using namespace std::placeholders;

    BToKDileptonAmplitudes<tag::GvDV2020>::BToKDileptonAmplitudes(const Parameters & p,
            const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c(p["mass::c"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        f_B(p["decay-constant::B_" + o.get("q", "d")], *this),
        f_K(p["decay-constant::K_" + o.get("q", "d")], *this),
        lambda_B_p_inv(p["B::1/lambda_B_p"], *this),
        a_1(p["K::a_1@1GeV"], *this),
        a_2(p["K::a_2@1GeV"], *this),
        lambda_psd(p["B->Pll::Lambda_pseudo@LargeRecoil"], *this),
        sl_phase_psd(p["B->Pll::sl_phase_pseudo@LargeRecoil"], *this),
        opt_nonlocal_formfactor(o, "nonlocal-formfactor", { "GvDV2020", "GRvDV2021", "naive" }, "GvDV2020"),
        nonlocal_formfactor(NonlocalFormFactor<nff::PToP>::make("B->K::" + opt_nonlocal_formfactor.value(), p, o))
    {
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
            throw InvalidOptionValueError("q", spectator_quark, "u, d");

        // Select the appropriate calculator for the QCDF integrals
        std::string qcdf_integrals(o.get("qcdf-integrals", "mixed"));
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
            throw InvalidOptionValueError("qcdf-integrals", qcdf_integrals, "mixed, numerical, analytical");
        }
    }

    BToKDileptonAmplitudes<tag::GvDV2020>::~BToKDileptonAmplitudes()
    {
    }

    BToKDilepton::DipoleFormFactors
    BToKDileptonAmplitudes<tag::GvDV2020>::dipole_form_factors(const double & s, const WilsonCoefficients<BToS> & wc) const
    {
        //Implementation of calT follows BFS2004 but C1 and C2 contributions are set to zero as they are accounted for by calH

        // charges of down- and up-type quarks
        static const double e_d = -1.0 / 3.0;
        static const double e_u = +2.0 / 3.0;

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
        complex<double> Y_top_c =  6.0 * wc.c3() + 60.0 * wc.c5();
        complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

        // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
        // then replace b pole mass by the PS mass.
        complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole);
        Y_top += Y_top_b * CharmLoops::h(mu, s, m_b_PS);
        Y_top += Y_top_0 * CharmLoops::h(mu, s);
        Y_top += Y_top_;

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
                + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (m_B / (2.0 * m_b_PS)) * c8eff * CharmLoops::F89_massless(s, m_b_PS));


        // compute the factorizing contributions
        complex<double> C_psd = C0_top_psd + a_mu * (C1f_top_psd + C1nf_top_psd);

        /* parallel, top sector */
        // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
        // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
        complex<double> T0_top_psd_m = +e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
        // cf. [BHP2007], Eq. (B.2)
        complex<double> T1f_top_psd_p  = -(c7eff + wc.c7prime()) * (4.0 * m_B / energy) * invm1_psd * lambda_B_p_inv;
        // T1f_top_par_m = 0, cf. [BFS2001], Eq. (22), p. 7
        // cf. [BFS2001], Eq. (25), p. 7
        complex<double> T1nf_top_psd_p = -m_B / m_b_PS * (
                e_u * wc.c6() * qcdf_c.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
        // cf. [BFS2001], Eq. (26), pp. 7-8
        complex<double> T1nf_top_psd_m = -e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                + 6.0 * m_B / m_b_PS * (
                    (wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                    -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

        // Compute the nonfactorizing contributions
        complex<double> T_psd = a_mu_f * (T1f_top_psd_p + T1nf_top_psd_p)
            + (T0_top_psd_m + a_mu_f * T1nf_top_psd_m);

        // Subleading weak annihilation and hard spectator interaction contributions have only been
        // computed for calT_perp, not for calT_par ~ calT_psd.

        // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
        BToKDilepton::DipoleFormFactors result;
        result.calT = xi_pseudo(s) * C_psd + power_of<2>(M_PI) / 3.0 * (f_B * f_K) / m_B  * T_psd;

        return result;
    }

    double
    BToKDileptonAmplitudes<tag::GvDV2020>::xi_pseudo(const double & s) const
    {
        // cf. [BF2001], Eq. (22)
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

        // cf. [GvDV2020] Eq. (A.5)
        const double
            calF_plus   = form_factors->f_p(s),
            calF_time   = form_factors->f_0(s),
            calF_T_plus = s / m_B / (m_B + m_K) * form_factors->f_t(s);

        const complex<double> calH_plus = nonlocal_formfactor->H_plus(s);

        double F_Tkin = calF_T_plus / calF_plus * 2.0 * std::sqrt(lambda(s)) * beta_l(s) * m_B / s;
        double F_Skin = calF_time / calF_plus * 0.5 * (m_B2 - m_K2) / (m_b_MSbar - m_s_MSbar);

        // cf. [BHP2007], Eq. (3.2), p. 3 and 4
        result.F_A  = wc.c10() + wc.c10prime();
        result.F_T  = F_Tkin * wc.cT();
        result.F_T5 = F_Tkin * wc.cT5();
        result.F_S  = F_Skin * (wc.cS() + wc.cSprime());
        result.F_P  = F_Skin * (wc.cP() + wc.cPprime()) + m_l() * (wc.c10() + wc.c10prime()) *
                      ((m_B2 - m_K2) / s * (calF_time / calF_plus - 1.0) - 1.0);
        result.F_V  = wc.c9() + wc.c9prime()
                      + 2.0 * m_b_PS() / m_B / xi_pseudo(s) * (dff.calT - 16.0 * pow(M_PI, 2) * m_B / m_b_PS() * calH_plus)
                      + 8.0 * m_l * m_B / s * calF_T_plus / calF_plus * wc.cT();

        return result;
    }
}
