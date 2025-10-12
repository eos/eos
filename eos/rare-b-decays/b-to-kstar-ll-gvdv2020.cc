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
#include <eos/rare-b-decays/b-to-kstar-ll-gvdv2020.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

using namespace std;

namespace eos
{
    using namespace std::literals::string_literals;
    using namespace std::placeholders;
    using std::norm;
    using std::sqrt;

    BToKstarDileptonAmplitudes<tag::GvDV2020>::BToKstarDileptonAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        f_B(p["decay-constant::B_" + o.get("q"_ok, "d")], *this),
        f_Kstar_par(p["B->K^*::f_Kstar_par"], *this),
        lambda_B_p_inv(p["B::1/lambda_B_p"], *this),
        q(o, options, "q"_ok),
        opt_nonlocal_formfactor(o, options, "nonlocal-formfactor"_ok),
        nonlocal_formfactor(NonlocalFormFactor<PToV>::make("B->K^*::" + opt_nonlocal_formfactor.value(), p, o))
    {
        Context ctx("When constructing B->K^*ll GVdV2020 amplitudes");
    }

    const std::vector<OptionSpecification>
    BToKstarDileptonAmplitudes<tag::GvDV2020>::options
    {
        { "q"_ok, { "d"s, "u"s },  "d"s },
        { "nonlocal-formfactor"_ok, { "GvDV2020"s, "GRvDV2022order5"s }, "GvDV2020"s }
    };

    BToKstarDilepton::FormFactorCorrections
    BToKstarDileptonAmplitudes<tag::GvDV2020>::sb_contributions(const double & s, const WilsonCoefficients<BToS> & wc) const
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
        double m_B2 = power_of<2>(m_B());
        double m_V2 = power_of<2>(m_Kstar());
        double energy = this->energy(s);

        // Coupling
        double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
            lambda_hat_u = std::conj(lambda_hat_u);

        /* Effective wilson coefficients */
        complex<double> c8eff = ShortDistanceLowRecoil::c8eff(wc); // LO C8eff

        /* Y(s) for the up and the top sector */
        // cf. [BFS2001], Eq. (10), p. 4
        complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
        complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6());
        complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

        // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
        // then replace b pole mass by the PS mass.
        // CharmLoops::h(mu, s, m_c_pole) contributions are set to zero
        complex<double> Y_top = Y_top_b * CharmLoops::h(mu, s, m_b_PS)
                              + Y_top_0 * CharmLoops::h(mu, s)
                              + Y_top_;
        // cf. [BFS2004], Eq. (43), p. 24
        complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (- CharmLoops::h(mu, s));

        complex<double> Y_contribution = Y_top + lambda_hat_u * Y_up;

        ////////////////////////////////////////
        // compute the factorizing contributions
        ////////////////////////////////////////

        complex<double> vector_contribution = alpha_s_mu / 4.0 / M_PI * (
            wc.c1() * memoise(CharmLoops::F19_massive_Qsb, s) + wc.c2() * memoise(CharmLoops::F29_massive_Qsb, s)
            + c8eff * CharmLoops::F89_massless(s, m_b_PS) + lambda_hat_u * (
                wc.c1() * (memoise(CharmLoops::F19_massive_Qsb, s) - CharmLoops::F19_massless(mu, s, m_b_PS))
                + wc.c2() * (memoise(CharmLoops::F29_massive_Qsb, s) - CharmLoops::F29_massless(mu, s, m_b_PS))
            )
        );

        complex<double> tensor_contribution = alpha_s_mu / 4.0 / M_PI * (
            (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive_Qsb, s) + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
            + lambda_hat_u * (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive_Qsb, s) - CharmLoops::F27_massless(mu, s, m_b_PS))
        );

        BToKstarDilepton::FormFactorCorrections result;

        result.t   = s / 2.0 / m_B2 * (Y_contribution - vector_contribution);
        result.t_T = - m_b_PS / m_B * tensor_contribution;


        ///////////////////////////////////////////
        // Compute the nonfactorizing contributions
        ///////////////////////////////////////////

        // inverse of the "negative" moment of the B meson LCDA
        // cf. [BFS2001], Eq. (54), p. 15
        double lambda_B_p_inv = this->lambda_B_p_inv;
        double omega_0 = 1.0 / lambda_B_p_inv;
        complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

        /* parallel, top sector */
        // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
        // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
        complex<double> T0_top_par_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;

        /* parallel, up sector */
        // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
        complex<double> T0_up_par_m = +e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;
        complex<double> T_par = T0_top_par_m + lambda_hat_u * T0_up_par_m;

        result.t_wa = - m_b_PS * s * eos::lambda(m_B2, m_V2, s) / 96.0 / power_of<5>(m_B()) / (m_B2 - m_V2) * (f_B * f_Kstar_par) / energy * T_par;

        return result;
    }

    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::mu_f() const
    {
        return 1.5;
    }

    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    // For testing purposes, compute the ratio between Qc and non-Qc non-local contributions
    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::H_perp_corrections(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);
        // Contributions not probortional to Qc
        auto sb_c = sb_contributions(s, wc);

        const double
            ff_V  = form_factors->v(s),
            ff_T1 = form_factors->t_1(s);

        const double
            m_B         = this->m_B(), m_B2 = power_of<2>(m_B),
            m_V         = this->m_Kstar(), m_V2 = power_of<2>(m_V),
            lambda      = eos::lambda(m_B2, m_V2, s),
            sqrt_lambda = std::sqrt(lambda);

        const double
            calF_perp = sqrt(2.0) * sqrt_lambda / (m_B * (m_B + m_V)) * ff_V,
            calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1;

        const double abs_Hsb_perp = std::abs(1.0 / 16.0 / power_of<2>(M_PI) * (calF_perp * sb_c.t + calF_T_perp * sb_c.t_T));
        const double abs_Hc_perp = std::abs(nonlocal_formfactor->H_perp(s));

        return  abs_Hsb_perp / abs_Hc_perp;
    }

    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::H_para_corrections(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);
        // Contributions not probortional to Qc
        auto sb_c = sb_contributions(s, wc);

        const double
            ff_A1 = form_factors->a_1(s),
            ff_T2 = form_factors->t_2(s);

        const double
            m_B         = this->m_B(), m_B2 = power_of<2>(m_B),
            m_V         = this->m_Kstar(), m_V2 = power_of<2>(m_V);

        const double
            calF_para = sqrt(2.0) * (m_B + m_V) / m_B * ff_A1,
            calF_T_para = sqrt(2.0) * (m_B2 - m_V2) / m_B2 * ff_T2;

        const double abs_Hsb_para = std::abs(1.0 / 16.0 / power_of<2>(M_PI) * (calF_para * sb_c.t + calF_T_para * sb_c.t_T));
        const double abs_Hc_para = std::abs(nonlocal_formfactor->H_para(s));

        return  abs_Hsb_para / abs_Hc_para;
    }

    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::H_long_corrections(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);
        // Contributions not probortional to Qc
        auto sb_c = sb_contributions(s, wc);

        const double
            ff_A1 = form_factors->a_1(s),
            ff_A2 = form_factors->a_2(s),
            ff_T2 = form_factors->t_2(s),
            ff_T3 = form_factors->t_3(s);

        const double
            m_B         = this->m_B(), m_B2 = power_of<2>(m_B),
            m_V         = this->m_Kstar(), m_V2 = power_of<2>(m_V),
            lambda      = eos::lambda(m_B2, m_V2, s);

        const double
            calF_long = ((m_B2 - m_V2 - s) * power_of<2>(m_B + m_V) * ff_A1 - lambda * ff_A2)
                      / (2.0 * m_V * m_B2 * (m_B + m_V)),
            calF_T_long = s / (2.0 * power_of<3>(m_B) * m_V) *
                        ((m_B2 + 3.0 * m_V2 - s) * ff_T2 - lambda / (m_B2 - m_V2) * ff_T3);

        const double abs_Hsb_long = std::abs(1.0 / 16.0 / power_of<2>(M_PI) * (calF_long * sb_c.t + calF_T_long * sb_c.t_T) - sb_c.t_wa);
        const double abs_Hc_long = std::abs(nonlocal_formfactor->H_long(s));

        return  abs_Hsb_long / abs_Hc_long;
    }

    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::GvDV2020>::amplitudes(const double & s) const
    {
        BToKstarDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        // local form factors
        const double
            ff_V  = form_factors->v(s),
            ff_A0 = form_factors->a_0(s),
            ff_A1 = form_factors->a_1(s),
            ff_A2 = form_factors->a_2(s),
            ff_T1 = form_factors->t_1(s),
            ff_T2 = form_factors->t_2(s),
            ff_T3 = form_factors->t_3(s);

        // kinematics
        const double
            sqrt_s      = std::sqrt(s),
            m_B         = this->m_B(), m_B2 = power_of<2>(m_B),
            m_V         = this->m_Kstar(), m_V2 = power_of<2>(m_V),
            lambda      = eos::lambda(m_B2, m_V2, s),
            sqrt_lambda = std::sqrt(lambda);

        // vectorial form factors, cf. [GvDV2020], eq. (A.11)
        const double
            calF_perp = sqrt(2.0) * sqrt_lambda / (m_B * (m_B + m_V)) * ff_V,
            calF_para = sqrt(2.0) * (m_B + m_V) / m_B * ff_A1,
            calF_long = ((m_B2 - m_V2 - s) * power_of<2>(m_B + m_V) * ff_A1 - lambda * ff_A2)
                        / (2.0 * m_V * m_B2 * (m_B + m_V)),
            calF_time = ff_A0;

        // tensorial form factors, cf. [GvDV2020], eq. (A.11)
        const double
            calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1,
            calF_T_para = sqrt(2.0) * (m_B2 - m_V2) / m_B2 * ff_T2,
            calF_T_long = s / (2.0 * power_of<3>(m_B) * m_V) *
                        ((m_B2 + 3.0 * m_V2 - s) * ff_T2 - lambda / (m_B2 - m_V2) * ff_T3);

        // Contributions not probortional to Qc
        auto sb_c = sb_contributions(s, wc);

        const complex<double>
            calH_perp = nonlocal_formfactor->H_perp(s) - 1.0 / 16.0 / power_of<2>(M_PI) * (calF_perp * sb_c.t + calF_T_perp * sb_c.t_T),
            calH_para = nonlocal_formfactor->H_para(s) - 1.0 / 16.0 / power_of<2>(M_PI) * (calF_para * sb_c.t + calF_T_para * sb_c.t_T),
            calH_long = nonlocal_formfactor->H_long(s) - 1.0 / 16.0 / power_of<2>(M_PI) * (calF_long * sb_c.t + calF_T_long * sb_c.t_T) - sb_c.t_wa;

        // Wilson coefficients
        const complex<double>
            c7eff = ShortDistanceLowRecoil::c7eff(s, 0.0, 0.0, 0.0, false, wc); // LO C7eff
        const complex<double>
            c910_m_r = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime()),
            c910_m_l = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime()),
            c910_p_r = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime()),
            c910_p_l = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime()),
            c7_m = (c7eff - wc.c7prime()),
            c7_p = (c7eff + wc.c7prime());

        // quark masses
        const double
            m_b_msbar = model->m_b_msbar(mu()),
            m_s_msbar = model->m_s_msbar(mu());

        // normalization constant, cf. KM2005A (3.7)
        const double calN = g_fermi() * alpha_e * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * sqrt(s * beta_l(s) * sqrt_lambda / (3.0 * 1024 * power_of<5>(M_PI) * m_B));

        // vector amplitudes, cf. KM2005A (3.2) - (3.4)
        result.a_long_right = -calN * m_B / sqrt_s * (
                c910_m_r * calF_long
                + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_long - 16.0 * power_of<2>(M_PI) * m_B * calH_long)
        );
        result.a_long_left  = -calN * m_B / sqrt_s * (
                c910_m_l * calF_long
                + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_long - 16.0 * power_of<2>(M_PI) * m_B * calH_long)
        );

        result.a_para_right = -calN * (
                c910_m_r * calF_para
                + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_para - 16.0 * power_of<2>(M_PI) * m_B * calH_para)
        );
        result.a_para_left  = -calN * (
                c910_m_l * calF_para
                + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_para - 16.0 * power_of<2>(M_PI) * m_B * calH_para)
        );

        result.a_perp_right = +calN * (
                c910_p_r * calF_perp
                + 2.0 * m_B / s * ((m_b_msbar + m_s_msbar) * c7_p * calF_T_perp - 16.0 * power_of<2>(M_PI) * m_B * calH_perp)
        );
        result.a_perp_left  = +calN * (
                c910_p_l * calF_perp
                + 2.0 * m_B / s * ((m_b_msbar + m_s_msbar) * c7_p * calF_T_perp - 16.0 * power_of<2>(M_PI) * m_B * calH_perp)
        );

        // scalar amplitude, cf. KM2005A (3.5)
        result.a_time = calN / m_B * sqrt_lambda / sqrt_s * calF_time * (
            2.0 * (wc.c10() - wc.c10prime()) + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime())
        );

        // Tensor amplitudes, cf BHvD2012 (B.17)-(B.20) and GVdV2020 (A.11)
        result.a_scal = -2.0 * calN / m_B * sqrt_lambda * calF_time * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar);

        result.a_para_perp = 2.0 * calN * m_B2 / s * calF_T_long * wc.cT();
        result.a_time_long = - 2.0 * calN * m_B2 / s * calF_T_long * wc.cT5();

        result.a_time_perp = sqrt(2) * calN * m_B / sqrt_s * calF_T_perp * wc.cT();
        result.a_long_perp = - sqrt(2) * calN * m_B / sqrt_s * calF_T_perp * wc.cT5();

        result.a_long_para = sqrt(2) * calN * m_B / sqrt_s * calF_T_para * wc.cT();
        result.a_time_para = - sqrt(2) * calN * m_B / sqrt_s * calF_T_para* wc.cT5();

        return result;
    }

    // C9 and its corrections [BFS2001] eqs. (40-41)
    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::real_C9_perp(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const complex<double>
            c7eff = ShortDistanceLowRecoil::c7eff(s, 0.0, 0.0, 0.0, false, wc); // LO C7eff

        const double
            m_B = this->m_B(), m_B2 = power_of<2>(m_B),
            m_V         = this->m_Kstar(), m_V2 = power_of<2>(m_V),
            lambda      = eos::lambda(m_B2, m_V2, s),
            sqrt_lambda = std::sqrt(lambda);

        const double
            m_b_msbar = model->m_b_msbar(mu()),
            m_s_msbar = model->m_s_msbar(mu());

        const double
            ff_V  = form_factors->v(s),
            ff_T1 = form_factors->t_1(s),
            calF_perp = sqrt(2.0) * sqrt_lambda / (m_B * (m_B + m_V)) * ff_V,
            calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1;

        auto sb_c = sb_contributions(s, wc);

        const complex<double>
            calH_perp = nonlocal_formfactor->H_perp(s) - 1.0 / 16.0 / power_of<2>(M_PI) * (calF_perp * sb_c.t + calF_T_perp * sb_c.t_T);

        return real(wc.c9() + 2.0 * m_B / s * (
            (m_b_msbar + m_s_msbar) * c7eff * calF_T_perp/calF_perp
            - 16.0 * power_of<2>(M_PI) * m_B * calH_perp/calF_perp
            ));
    }
    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::imag_C9_perp(const double & s) const
    {
        WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), lepton_flavor, cp_conjugate);

        const complex<double>
            c7eff = ShortDistanceLowRecoil::c7eff(s, 0.0, 0.0, 0.0, false, wc); // LO C7eff

        const double
            m_B = this->m_B(), m_B2 = power_of<2>(m_B),
            m_V         = this->m_Kstar(), m_V2 = power_of<2>(m_V),
            lambda      = eos::lambda(m_B2, m_V2, s),
            sqrt_lambda = std::sqrt(lambda);

        const double
            m_b_msbar = model->m_b_msbar(mu()),
            m_s_msbar = model->m_s_msbar(mu());

        const double
            ff_V  = form_factors->v(s),
            ff_T1 = form_factors->t_1(s),
            calF_perp = sqrt(2.0) * sqrt_lambda / (m_B * (m_B + m_V)) * ff_V,
            calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1;

        auto sb_c = sb_contributions(s, wc);

        const complex<double>
            calH_perp = nonlocal_formfactor->H_perp(s) - 1.0 / 16.0 / power_of<2>(M_PI) * (calF_perp * sb_c.t + calF_T_perp * sb_c.t_T);

        return imag(wc.c9() + 2.0 * m_B / s * (
            (m_b_msbar + m_s_msbar) * c7eff * calF_T_perp/calF_perp
            - 16.0 * power_of<2>(M_PI) * m_B * calH_perp/calF_perp
            ));
    }
    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::real_C9_para(const double & ) const
    {
        return 0.0;
    }
    double
    BToKstarDileptonAmplitudes<tag::GvDV2020>::imag_C9_para(const double & ) const
    {
        return 0.0;
    }
}
