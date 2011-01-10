/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#include <src/rare-b-decays/charm-loops.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/form-factors.hh>
#include <src/rare-b-decays/hard-scattering.hh>
#include <src/rare-b-decays/long-distance.hh>
#include <src/utils/destringify.hh>
#include <src/utils/integrate.hh>
#include <src/utils/memoise.hh>
#include <src/utils/model.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>
#include <src/utils/save.hh>

#include <cmath>
#include <functional>
#include <utility>
#include <map>
#include <vector>

#include <gsl/gsl_sf.h>

namespace eos
{
    using std::norm;

    // Large Recoil

    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter abs_c7;

        Parameter arg_c7;

        Parameter c7prime;

        Parameter c8;

        Parameter abs_c9;

        Parameter arg_c9;

        Parameter c9prime;

        Parameter abs_c10;

        Parameter arg_c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        double m_l;

        Parameter mu;

        Parameter f_B;

        Parameter f_Kstar_par;

        Parameter f_Kstar_perp;

        Parameter lambda_B_p;

        Parameter a_1_par;

        Parameter a_2_par;

        Parameter a_1_perp;

        Parameter a_2_perp;

        Parameter uncertainty_par_left;

        Parameter uncertainty_par_right;

        Parameter uncertainty_perp_left;

        Parameter uncertainty_perp_right;

        Parameter uncertainty_long_left;

        Parameter uncertainty_long_right;

        Parameter uncertainty_xi_perp;

        Parameter uncertainty_xi_par;

        double e_q;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
            model(Model::make("SM", p)),
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            abs_c7(p["Abs{c7}"]),
            arg_c7(p["Arg{c7}"]),
            c7prime(p["c7prime"]),
            c8(p["c8"]),
            abs_c9(p["Abs{c9}"]),
            arg_c9(p["Arg{c9}"]),
            c9prime(p["c9prime"]),
            abs_c10(p["Abs{c10}"]),
            arg_c10(p["Arg{c10}"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            f_B(p["f_B"]),
            f_Kstar_par(p["B->K^*::f_Kstar_par"]),
            f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"]),
            lambda_B_p(p["lambda_B_p"]),
            a_1_par(p["B->K^*::a_1_par"]),
            a_2_par(p["B->K^*::a_1_par"]),
            a_1_perp(p["B->K^*::a_1_perp"]),
            a_2_perp(p["B->K^*::a_2_perp"]),
            uncertainty_par_left(p["B->K^*ll::A_par^L_uncertainty@LargeRecoil"]),
            uncertainty_par_right(p["B->K^*ll::A_par^R_uncertainty@LargeRecoil"]),
            uncertainty_perp_left(p["B->K^*ll::A_perp^L_uncertainty@LargeRecoil"]),
            uncertainty_perp_right(p["B->K^*ll::A_perp^R_uncertainty@LargeRecoil"]),
            uncertainty_long_left(p["B->K^*ll::A_0^L_uncertainty@LargeRecoil"]),
            uncertainty_long_right(p["B->K^*ll::A_0^R_uncertainty@LargeRecoil"]),
            uncertainty_xi_perp(p["formfactors::xi_perp_uncertainty"]),
            uncertainty_xi_par(p["formfactors::xi_par_uncertainty"]),
            e_q(-1.0/3.0)
        {
            form_factors = FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "BZ2004"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.0;//m_l = 0.10565836; // (GeV), cf. [PDG2008], p. 13
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5

            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(double s) const
        {
            return s / m_B / m_B;
        }

        inline double energy(const double & s) const
        {
            return (m_B * m_B + m_Kstar * m_Kstar - s) / (2.0 * m_B);
        }

        inline double mu_f() const
        {
            return 1.5;
        }

        inline double m_b_PS() const
        {
            // Actually use the PS mass at mu_f = 1.5 GeV
            return model->m_b_ps(1.5);
        }

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        complex<double> c7eff() const
        {
            return abs_c7() * complex<double>(cos(arg_c7()), sin(arg_c7())) - 1.0/3.0 * c3() - 4.0/9.0 * c4() - 20.0/3.0 * c5() - 80.0/9.0 * c6();
        }

        // cf. [BFS2001], below Eq. (26), p. 8
        double c8eff() const
        {
            return c8() + c3() - 1.0/6.0 * c4() + 20.0 * c5() - 10.0/3.0 * c6();
        }

        complex<double> c9() const
        {
            return abs_c9() * complex<double>(cos(arg_c9()), sin(arg_c9()));
        }

        complex<double> c10() const
        {
            return abs_c10() * complex<double>(cos(arg_c10()), sin(arg_c10()));
        }

        // cf. [BFS2001], Eq. (10), p. 4
        complex<double> Y0(const double & s) const
        {
            double Y_c = 4.0 / 3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5();
            double Y_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double Y_0 = -0.5 * (c3() + 4.0 / 3.0 * c4() + 16.0 * c5() + 64 / 3.0 * c6());
            double Y = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * CharmLoops::h(mu, s, m_c)
                + Y_b * CharmLoops::h(mu, s, m_b_PS())
                + Y_0 * CharmLoops::h(mu, s)
                + Y;
        }

        // cf. [BFS2004], ?
        complex<double> Y0u(const double & s) const
        {
            double a = 4.0 / 3.0 * c1() + c2();

            return a * (CharmLoops::h(mu, s, m_c) - CharmLoops::h(mu, s));
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);
            double result = uncertainty_xi_perp * factor * form_factors->v(s);

            return result;
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B + m_Kstar) / (2.0 * energy(s));
            const double factor2 = (1.0 - m_Kstar / m_B);
            double result = uncertainty_xi_par * (factor1 * form_factors->a_1(s) - factor2 * form_factors->a_2(s));

            return result;
        }

        /* NLO functions */
        // cf. [BFS2001], Eq. (48)
        double phi_K(const double & u, const double & a_1, const double & a_2) const
        {
            double xi = 2.0 * u - 1.0;

            return 6.0 * u * (1 - u) * (1.0 + a_1 * 3.0 * xi + a_2 * (7.5 * xi * xi - 1.5));
        }

        // cf. [BFS2001], Eq. (27), p. 8
        complex<double> t_perp(const double & s, const double & u, const double & m_q) const
        {
            if (0.0 == s)
                return t_perp_0(u, m_q);

            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            complex<double> result = (2.0 * m_B / ubar / E) * memoise(HardScattering::I1, s, u, m_q, m_B());
            if (m_q > 0.0)
                result = result + (s / ubar / ubar / E / E) * (CharmLoops::B0(x, m_q) - CharmLoops::B0(s, m_q));

            return result;
        }

        complex<double> t_perp_0(const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q, m_B2 = m_B * m_B;
            double a, a2, sign;
            complex<double> dilogArg, dilog1, dilog2;
            complex<double> LxpLxm;
            gsl_sf_result res_re, res_im;

            if (m_q > 0)
            { // m != 0
                if (1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2) > 0)
                {
                    a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2)))
                        / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * m_B2)));
                    LxpLxm = -M_PI * M_PI / 3.0 + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                        + gsl_sf_dilog(-a) + gsl_sf_dilog(-1.0/a);
                }
                else
                {
                    a = std::sqrt(4.0 * m_q2 / (m_B2 - u * m_B2) - 1.0);
                    a2 = a * a;

                    if (a2 - 1 > 0)
                        sign = +1.0;
                    else
                        sign = -1.0;

                    dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
                    gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                    dilog1 = complex<double>(res_re.val, res_im.val);
                    dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
                    gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                    dilog2 = complex<double>(res_re.val, res_im.val);

                    LxpLxm = -M_PI * M_PI / 3.0 - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                        + dilog1 + dilog2;
                }

                return 4.0 / ubar * (1.0 + 2.0 * m_q2 / ubar / m_B2 * LxpLxm);
            }
            else
            {
                return complex<double>(4.0 / ubar, 0.0);
            }
        }

        // cf. [BFS2001], Eq. (28), p. 8
        complex<double> t_par(const double & s, const double & u, const double & m_q) const
        {
            double ubar = 1.0 - u;
            double E = energy(s);
            double x = ubar * m_B * m_B + u * s;

            complex<double> result = (2.0 * m_B / ubar / E) * memoise(HardScattering::I1, s, u, m_q, m_B());
            if (m_q > 0.0)
                result = result + (x / ubar / ubar / E / E) * (CharmLoops::B0(x, m_q) - CharmLoops::B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (23), p. 7, multiplied by phi_K^*,perp
        complex<double> Tnf_perp_p(const double & s, const double & u) const
        {
            double m_b = m_b_PS();
            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;

            double a = (4.0 / 3.0 / (u + ubar * s_hat)) * c8eff();
            complex<double> ba = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_perp(s, u, m_c);
            complex<double> bb = (-1.0 / 3.0)
                * (c3 - c4 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6 - (4.0 * m_b / m_B) * (c3 - c4/6.0 + 4.0 * c5 - 2.0/3.0 * c6))
                * t_perp(s, u, m_b);
            complex<double> bc = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_perp(s, u, 0.0);

            return (a + (m_B / 2.0 / m_b) * (ba + bb + bc)) * phi_K(u, a_1_perp, a_2_perp);
        }

        // cf. [BFS2001, Eq. (25), p. 7, multiplied by phi_K^*,par
        complex<double> Tnf_par_p(const double & s, const double & u) const
        {
            double m_b = m_b_PS();

            complex<double> a = (+2.0 / 3.0) * (-c1 / 6.0 + c2 + 6.0 * c6) * t_par(s, u, m_c);
            complex<double> b = (-1.0 / 3.0) * (c3 - c1 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6) * t_par(s, u, m_b);
            complex<double> c = (-1.0 / 3.0) * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_par(s, u, 0.0);

            return (m_B / m_b) * (a + b + c) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (26), pp. 7-8, multiplied by phi_K^*,par
        complex<double> Tnf_par_m(const double & s, const double & u) const
        {
            double m_b = m_b_PS();

            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;
            double x = ubar * m_B * m_B + u * s;

            double a = (e_q * 8.0 / (ubar + u * s_hat)) * c8eff();
            complex<double> ba = (-c1 / 6.0 + c2 + c4 + 10 * c6) * CharmLoops::h(mu, x, m_c);
            complex<double> bb = (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * CharmLoops::h(mu, x, m_b);
            complex<double> bc = (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * CharmLoops::h(mu, x);
            double bd = -8.0 / 27.0 * (-7.5 * c4 + 12.0 * c5 - 32.0 * c6);

            return (a + (6.0 * m_B / m_b) * (ba + bb + bc + bd)) * phi_K(u, a_1_par, a_2_par);
        }

        // cf. [BFS2001], Eq. (36), p. 9
        double L(const double & s) const
        {
            double m_b = m_b_PS();
            double m_b2 = m_b * m_b;

            return -1.0 * (m_b2 - s) / s * std::log(1.0 - s / m_b2);
        }

        // cf. [BFS2001], Eq. (54), p. 15
        complex<double> lambda_B_m_inv(const double & s) const
        {
            if (0.0 == s)
                return complex<double>(0.0, 0.0);

            double omega_0 = lambda_B_p;
            double x = s / m_B / omega_0;
            double ei = gsl_sf_expint_Ei(x);

            complex<double> result = complex<double>(-ei, M_PI) * (std::exp(-x) / omega_0);

            return result;
        }

        // cf. [BFS2004], Eq. (51) the integrand of the first term only,
        // or [FM2002], Eq. (15) times a factor of 3.0, respectively
        double Twa_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;

            return phi_K(u, a_1_perp, a_2_perp) / (ubar + u * s_hat);
        }

        // cf. [FM2002], Eq. (17), the integrand only
        double Xperp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double denom = ubar + u * s_hat;

            return phi_K(u, a_1_perp, a_2_perp) * (1 / denom + 1 / denom / denom) / 3.0;
        }

        // cf. [FM2002], Eq. (22), p. 9
        complex<double> FV(const double & s) const
        {
            return 3.0 / 4.0 * (
                    (-c1 / 6.0 + c2 + c4 + 10.0 * c6) * CharmLoops::h(mu, s, m_c)
                    + (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * CharmLoops::h(mu, s, m_b_PS())
                    + (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * CharmLoops::h(mu, s)
                    - 8.0/27.0 * (-7.5 * c4 + 12 * c5 - 32 * c6));
        }

        // cf. [BFS2004], Eq. (52), the integrand of the second term only
        complex<double> Thsa_1_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double x = (ubar + u * s_hat) * m_B * m_B;

            return phi_K(u, a_1_perp, a_2_perp) / (ubar + u * s_hat) * FV(x);
        }

        // cf. [BFS2004], Eq. (52), the integrand of the third term only.
        // the v integration has been executed analytically
        complex<double> Thsa_2_perp(const double & s_hat, const double & u) const
        {
            double ubar = 1.0 - u;
            double x = (ubar + u * s_hat) * m_B * m_B;

            return 3.0 * u * u * (1.0 + a_1_par * (4.0 * u - 3.0) + a_2_par * (15.0 * u * u - 20.0 * u + 6.0))
                * FV(x);
        }

        complex<double> tensor_perp_hsa(const double & shat, const double & u) const
        {
            return 12.0 * c8eff() * m_b_PS() / m_B * f_Kstar_perp * Xperp(shat, u)
                    + 8.0 * f_Kstar_perp * Thsa_1_perp(shat, u)
                    - 4.0 * m_Kstar * f_Kstar_par / ((1.0 - shat) * lambda_B_p) * Thsa_2_perp(shat, u);
        }

        // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
        complex<double> C0_perp(const double & /*h*/, const double & s) const
        {
            return c7eff() + s / (2.0 * m_b_PS() * m_B) * Y0(s);
        }

        // cf. [BFS2001], Eqs. (34), (37), p. 9
        complex<double> C1_perp(const double & /*h*/, const double & s) const
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = mu * model->m_b_pole() / m_b;

            // cf. [BFS2004], Eq. (44), p. 24
            // [Christoph] Use c7 instead of c7eff
            complex<double> C_perp_f = (abs_c7() * complex<double>(cos(arg_c7()), sin(arg_c7())) - c7prime()) * (8.0 * std::log(m_b / mu) - L(s) - 4.0 * (1.0 - mu_f() / m_b));

            // cf. [BFS2001], Eq. (37), p. 9
            // [Christoph] Use c8 instead of c8eff
            complex<double> C_perp_nf = (-1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * memoise(CharmLoops::F27_massive, mu_pole, s, m_b, m_c()) + c8() * CharmLoops::F87_massless(mu_pole, s, m_b)
                    + (s / (2.0 * m_b * m_B)) * (
                        c1() * memoise(CharmLoops::F19_massive, mu_pole, s, m_b, m_c())
                        + c2() * memoise(CharmLoops::F29_massive, mu_pole, s, m_b, m_c())
                        + c8() * CharmLoops::F89_massless(s, m_b)));

            return C_perp_f + C_perp_nf;
        }

        // cf. [BFS2001], Eqs. (16), (21), (25), pp. 5-7
        complex<double> T1_perp_p(const double & s, const double & u) const
        {
            static const double e_d = (-1.0/3.0);
            static const double e_u = (+2.0/3.0);

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();
            double s_hat = s / m_B / m_B;

            // cf. [BFS2001], Eq. (20)
            // [Christoph] Use c7 instead of c7eff
            complex<double> Tf_perp_p = abs_c7() * complex<double>(cos(arg_c7()), sin(arg_c7())) * (2.0 * m_B / (1.0 - u) / energy(s));

            // cf. [BFS2001], Eq. (23)
            // [Christoph] Use c8 instead of c8eff
            complex<double> Tnf_perp_p = -4.0 * e_d * c8 / (u + (1.0 - u) * s_hat)
                + m_B / (2.0 * m_b) * (
                        e_u * (-c1 / 6.0 + c2 + 6.0 * c6) * t_perp(s, u, m_c)
                        + e_d * (c3 - c4 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6 - (4.0 * m_b / m_B) * (c3 - c4/6.0 + 4.0 * c5 - 2.0/3.0 * c6))
                        + e_d * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_perp(s, u, 0.0));

            return (Tf_perp_p + Tnf_perp_p) / lambda_B_p();
        }

        // cf. [BFS2001], Eq. (16) times phi_K^*_perp
        complex<double> T_perp(const double & s, const double & u) const
        {
            double a = model->alpha_s(sqrt(mu * 0.5)) * QCD::casimir_f / 4.0 / M_PI;

            return phi_K(u, a_1_perp, a_2_perp) * a * T1_perp_p(s, u);

            // TODO: Hard scattering + Weak annihilation from [BFS2004], Eqs. (51), (52)
        }

        // cf. [BFS2001], Eq. (15) with a = perp
        complex<double> calT_perp(const double & h, const double & s) const
        {
            return xi_perp(s) * (C0_perp(h, s) + model->alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI * C1_perp(h, s))
                + (pow(M_PI, 2) / 3.0) * (f_B * f_Kstar_perp / m_B)
                * integrate(std::function<complex<double> (const double &)>(
                                std::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::T_perp, this, s, std::placeholders::_1)),
                        64, 0.001, 0.999);
        }

        // cf. [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
        complex<double> C0_par(const double & s) const
        {
            return -1.0 * (c7eff() - c7prime() + m_B / (2.0 * m_b_PS()) * (Y0(s)));
        }

        // cf. [BFS2001], Eqs. (35), (38), p. 9
        complex<double> C1_par(const double & s) const
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = mu * model->m_b_pole() / m_b;

            // cf. [BFS2004], Eq. (45), p. 24
            // [Christoph] Use c7 instead of c7eff.
            complex<double> C_par_f = -1.0 * (abs_c7() * complex<double>(cos(arg_c7()), sin(arg_c7())) - c7prime()) * (8.0 * std::log(m_b / mu) + 2.0 * L(s) - 4.0 * (1.0 - mu_f() / m_b));
            /* for [BFS2001] version of xi_par we also needed: */
            // C_par_f += (m_B / (2.0 * m_b)) * Y0(s) * (2.0 - 2.0 * L(s));

            // cf. [BFS2001], Eq. (38), p. 9
            // [Christoph] Use c8 instead of c8eff.
            complex<double> C_par_nf = (+1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * memoise(CharmLoops::F27_massive, mu_pole, s, m_b, m_c()) + c8eff() * CharmLoops::F87_massless(mu_pole, s, m_b)
                    + (m_B / (2.0 * m_b)) * (
                        c1() * memoise(CharmLoops::F19_massive, mu_pole, s, m_b, m_c())
                        + c2() * memoise(CharmLoops::F29_massive, mu_pole, s, m_b, m_c())
                        + c8eff() * CharmLoops::F89_massless(s, m_b)));

            return C_par_f + C_par_nf;
        }

        // cf. [BFS2001], Eq. (18), p. 6 with \omega integrated out.
        complex<double> T0_par_m(const double & s) const
        {
            return -e_q * 4.0 * m_B / m_b_PS() * (c3 + 4.0/3.0 * c4 + 16.0 * c5 + 64.0/3.0 * c6) * lambda_B_m_inv(s);
        }

        // cf. [BFS2001], Eqs. (16), (21), (25), pp. 5-7
        complex<double> T1_par_p(const double & s, const double & u) const
        {
            static const double e_d = (-1.0/3.0);
            static const double e_u = (+2.0/3.0);

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();

            // cf. [BFS2001], Eq. (21)
            //complex<double> Tf_par_p = (c7() - c7prime + (s / (2.0 * m_b * m_B)) * Y0(s)) * (2.0 * pow(m_B / energy(s), 2));
            // cf. [BFS2004], Eq. (49)
            // [Christoph] Use c7 instead of c7eff.
            complex<double> Tf_par_p = abs_c7() * complex<double>(cos(arg_c7()), sin(arg_c7())) * (4.0 * m_B / (1.0 - u) / energy(s));

            // cf. [BFS2001], Eq. (25)
            complex<double> Tnf_par_p = m_B / m_b * (
                    e_u * (-c1 / 6.0 + c2 + 6.0 * c6) * t_par(s, u, m_c)
                    + e_d * (c3 - c1 / 6.0 + 16.0 * c5 + 10.0/3.0 * c6) * t_par(s, u, m_b)
                    + e_d * (c3 - c4 / 6.0 + 16.0 * c5 - 8.0/3.0 * c6) * t_par(s, u, 0.0));

            return (Tf_par_p + Tnf_par_p) / lambda_B_p();
        }

        // cf. [BFS2001], Eq. (16), (22), (26), pp. 5-8
        complex<double> T1_par_m(const double & s, const double & u) const
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = mu * model->m_b_pole() / m_b;

            double s_hat = s / m_B / m_B;
            double ubar = 1.0 - u;
            double x = ubar * m_B * m_B + u * s;

            // [Christoph] Use c8 instead of c8eff.
            return e_q * lambda_B_m_inv(s) * (
                    8.0 / (ubar + u * s_hat) * c8()
                    + 6.0 * m_B / m_b * (
                        (-c1 / 6.0 + c2 + c4 + 10 * c6) * CharmLoops::h(mu_pole, x, m_c)
                        + (c3 + 5.0/6.0 * c4 + 16.0 * c5 + 22.0/3.0 * c6) * CharmLoops::h(mu_pole, x, m_b)
                        + (c3 + 17.0/6.0 * c4 + 16.0 * c5 + 82.0/3.0 * c6) * CharmLoops::h(mu_pole, x)
                        -8.0 / 27.0 * (-7.5 * c4 + 12.0 * c5 - 32.0 * c6)));
        }

        // cf. [BFS2001], Eq. (16) times phi_K^*_par
        complex<double> T_par(const double & s, const double & u) const
        {
            double a = model->alpha_s(sqrt(mu * 0.5)) * QCD::casimir_f / 4.0 / M_PI;

            return phi_K(u, a_1_par, a_2_par) * (T0_par_m(s) + a * (T1_par_m(s, u) + T1_par_p(s, u)));
        }

        // cf. [BFS2001], Eq. (15) with a = par, and [BHP2008], Eq. (C.4)
        complex<double> calT_par(const double & s) const
        {
            return xi_par(s) * (C0_par(s) + model->alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI * C1_par(s))
                + (pow(M_PI, 2) / 3.0) * (f_B * f_Kstar_par / m_B) * (m_Kstar / energy(s))
                * integrate(std::function<complex<double> (const double &)>(
                                std::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::T_par, this, s, std::placeholders::_1)),
                        32, 0.001, 0.999);
        }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_long_left + (1.0 + h) / 2.0 * uncertainty_long_right;
            complex<double> wilson = (c9() - c9prime()) + h * (c10() - c10prime());
            double prefactor = -1.0 / (2.0 * m_Kstar * std::sqrt(s));

            complex<double> a = wilson * ((m_B * m_B - m_Kstar * m_Kstar - s) * 2.0 * energy(s) * xi_perp(s)
                -lambda(m_B * m_B, m_Kstar * m_Kstar, s) * m_B / (m_B * m_B - m_Kstar * m_Kstar) * (xi_perp(s) - xi_par(s)));
            complex<double> b = 2 * m_b_PS() * (((m_B * m_B + 3 * m_Kstar * m_Kstar - s) * 2.0 * energy(s) / m_B
                        - lambda(m_B * m_B, m_Kstar * m_Kstar, s) / (m_B * m_B - m_Kstar * m_Kstar)) * calT_perp(-1.0, s)
                    - lambda(m_B * m_B, m_Kstar * m_Kstar, s) / (m_B * m_B - m_Kstar * m_Kstar) * calT_par(s));

            return this->norm(s) * uncertainty * prefactor * (a + b);
        }

        // cf. [BHP2008], p. 20
        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double mKhat = m_Kstar / m_B;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_perp_left + (1.0 + h) / 2.0 * uncertainty_perp_right;
            double prefactor = +std::sqrt(2.0) * m_B * std::sqrt(lambda(1.0, mKhat * mKhat, shat));
            complex<double> wilson = (c9() + c9prime()) + h * (c10() + c10prime());

            return this->norm(s) * uncertainty * prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * calT_perp(+1.0, s));
        }

        // cf. [BHP2008], p. 20
        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double mKhat = m_Kstar / m_B;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_par_left + (1.0 + h) / 2.0 * uncertainty_par_right;
            double prefactor = -std::sqrt(2.0) * m_B * (1.0 - shat);
            complex<double> wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return this->norm(s) * uncertainty * prefactor * (wilson * xi_perp(s) + (2.0 * mbhat / shat) * (1.0 - mKhat * mKhat) * calT_perp(-1.0, s));
        }

        // Unormalized combinations of transversity amplitudes
        double u_1(const double & s) const
        {
            return std::norm(a_long(left_handed, s)) + std::norm(a_long(right_handed, s));
        }

        double u_2(const double & s) const
        {
            return std::norm(a_perp(left_handed, s)) + std::norm(a_perp(right_handed, s));
        }

        double u_3(const double & s) const
        {
            return std::norm(a_par(left_handed, s)) + std::norm(a_par(right_handed, s));
        }

        double u_4(const double & s) const
        {
            return real(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        }

        double u_5(const double & s) const
        {
            return real(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        }

        double u_6(const double & s) const
        {
            return real(a_par(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_par(right_handed, s)) * a_perp(right_handed, s));
        }

        double u_7(const double & s) const
        {
            return imag(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        }

        double u_8(const double & s) const
        {
            return imag(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s));
        }

        double u_9(const double & s) const
        {
            return imag(a_par(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_par(right_handed, s)) * a_perp(right_handed, s));
        }

        // Components of observables
        double unnormalized_decay_width(const double & s) const
        {
            return (u_1(s) + u_2(s) + u_3(s));
        }

        double a_fb_numerator(const double & s) const
        {
            return 1.5 * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_par(right_handed, s)) * a_perp(right_handed, s)));
        }

        double f_l_numerator(const double & s) const
        {
            return u_1(s);
        }

        double a_t_2_numerator(const double & s) const
        {
            return u_2(s) - u_3(s);
        }

        double a_t_2_denominator(const double & s) const
        {
            return u_2(s) + u_3(s);
        }

        double a_t_3_numerator(const double & s) const
        {
            return sqrt(pow(u_4(s), 2) + pow(u_7(s), 2));
        }

        double a_t_3_denominator(const double & s) const
        {
            return sqrt(u_1(s) * u_2(s));
        }

        double a_t_4_numerator(const double & s) const
        {
            return sqrt(pow(u_5(s), 2) + pow(u_8(s), 2));
        }

        double a_t_4_denominator(const double & s) const
        {
            return sqrt(pow(u_4(s), 2) + pow(u_7(s), 2));
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) / Gamma;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_decay_width(const double & s) const
    {
        return _imp->unnormalized_decay_width(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return _imp->a_fb_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        return _imp->a_t_2_numerator(s) / _imp->a_t_2_denominator(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        return _imp->a_t_3_numerator(s) / _imp->a_t_3_denominator(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        return _imp->a_t_4_numerator(s) / _imp->a_t_4_denominator(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_5(const double & s) const
    {
        return _imp->a_t_4_numerator(s) / _imp->a_t_3_denominator(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return _imp->f_l_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(std::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_fb_numerator), _imp, std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::unnormalized_decay_width), _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_unnormalized_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // Convert from asymmetry in the decay width to asymmetry in the BR
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        std::function<double (const double &)> f = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_fb_numerator), _imp, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max) / Gamma;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::f_l_numerator), _imp, std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::unnormalized_decay_width), _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_unnormalized_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // Convert from polarisation in the decay width to polarisation in the BR
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        std::function<double (const double &)> f = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::f_l_numerator), _imp, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max) / Gamma;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_2_numerator), _imp, std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_2_denominator), _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_3_numerator), _imp, std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_3_denominator), _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_4_numerator), _imp, std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_4_denominator), _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_4_numerator), _imp, std::placeholders::_1);

        std::function<double (const double &)> denom = std::bind(
                std::mem_fn(&Implementation<BToKstarDilepton<LargeRecoil>>::a_t_3_denominator), _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    // Low Recoil

    template <>
    struct Implementation<BToKstarDilepton<LowRecoil>>
    {
        std::shared_ptr<Model> model;

        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter abs_c7;

        Parameter arg_c7;

        Parameter c7prime;

        Parameter c8;

        Parameter abs_c9;

        Parameter arg_c9;

        Parameter c9prime;

        Parameter abs_c10;

        Parameter arg_c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c_MSbar;

        Parameter m_B;

        Parameter m_Kstar;

        Parameter mu;

        Parameter uncertainty_par_left;

        Parameter uncertainty_par_right;

        Parameter uncertainty_perp_left;

        Parameter uncertainty_perp_right;

        Parameter uncertainty_long_left;

        Parameter uncertainty_long_right;

        Parameter uncertainty_isgur_wise_long;

        Parameter uncertainty_isgur_wise_par;

        Parameter uncertainty_isgur_wise_perp;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        bool cp_conjugate;

        bool ccbar_resonance;

        Implementation(const Parameters & p, const ObservableOptions & o) :
            model(Model::make("SM", p)),
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            abs_c7(p["Abs{c7}"]),
            arg_c7(p["Arg{c7}"]),
            c7prime(p["c7prime"]),
            c8(p["c8"]),
            abs_c9(p["Abs{c9}"]),
            arg_c9(p["Arg{c9}"]),
            c9prime(p["c9prime"]),
            abs_c10(p["Abs{c10}"]),
            arg_c10(p["Arg{c10}"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c_MSbar(p["mass::c"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            uncertainty_par_left(p["B->K^*ll::A_par^L_uncertainty@LowRecoil"]),
            uncertainty_par_right(p["B->K^*ll::A_par^R_uncertainty@LowRecoil"]),
            uncertainty_perp_left(p["B->K^*ll::A_perp^L_uncertainty@LowRecoil"]),
            uncertainty_perp_right(p["B->K^*ll::A_perp^R_uncertainty@LowRecoil"]),
            uncertainty_long_left(p["B->K^*ll::A_0^L_uncertainty@LowRecoil"]),
            uncertainty_long_right(p["B->K^*ll::A_0^R_uncertainty@LowRecoil"]),
            uncertainty_isgur_wise_long(p["B->K^*ll::IW_long_uncertainty"]),
            uncertainty_isgur_wise_par(p["B->K^*ll::IW_par_uncertainty"]),
            uncertainty_isgur_wise_perp(p["B->K^*ll::IW_perp_uncertainty"]),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            ccbar_resonance(destringify<bool>(o.get("ccbar-resonance", "false")))
        {
            form_factors = FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "BZ2004"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");
        }

        // We use the PS mass except for kappa_1
        double m_b_PS() const
        {
            // Actually use m_b_PS at mu_PS = 2.0 GeV
            return model->m_b_ps(2.0);
        }

        inline complex<double> c7() const
        {
            double sign = (cp_conjugate ? -1.0 : +1.0);

            return abs_c7() * complex<double>(cos(arg_c7()), sign * sin(arg_c7()));
        }

        static complex<double> c7eff_nlo(const double & mu, const double & s, const double & m_b,
                const double & c1, const double & c2, const double & c8)
        {
            return -1.0 * (
                    c1 * CharmLoops::F17_massless(mu, s, m_b)
                    + c2 * CharmLoops::F27_massless(mu, s, m_b)
                    + c8 * CharmLoops::F87_massless(mu, s, m_b));
        }

        // cf. [GP2004], Eq. (56)
        complex<double> c7eff(double s) const
        {
            double m_b = m_b_PS();

            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            double lo = - 1.0/3.0 * c3 - 4.0/9.0 * c4 - 20.0/3.0 * c5 - 80.0/9.0 * c6;
            complex<double> nlo = memoise(c7eff_nlo, mu(), s, m_b, c1(), c2(), c8());

            return c7() + lo + (model->alpha_s(mu) / (4.0 * M_PI)) * nlo;
        }

        inline complex<double> c9() const
        {
            double sign = (cp_conjugate ? -1.0 : +1.0);

            return abs_c9() * complex<double>(cos(arg_c9()), sign * sin(arg_c9()));
        }

        static complex<double> c9eff_nlo_alpha_s(const double & mu, const double & s, const double & m_b,
                const double & c1, const double & c2, const double & c8)
        {
            return -1.0 * (
                    c1 * CharmLoops::F19_massless(mu, s, m_b)
                    + c2 * CharmLoops::F29_massless(mu, s, m_b)
                    + c8 * CharmLoops::F89_massless(s, m_b));
        }

        // cf. [GP2004], Eq. (55), p. 10
        complex<double> c9eff(const double & s) const
        {
            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            // Substitute pole mass by PS mass
            double m_b = m_b_PS();
            double m_c = model->m_c_msbar(mu);
            double c = -2.0 / 27.0 * (8.0 * c1() + 6.0 * c2() - 6.0 * c3() - 8.0 * c4() - 12.0 * c5() - 160.0 * c6());
            double c_0 = -2.0 / 27.0 * (48.0 * c1() + 36.0 * c2() + 198.0 * c3() - 24.0 * c4() + 1872.0 * c5() - 384.0 * c6());
            double c_b = +2.0 / 27.0 * (126.0 * c3() + 24.0 * c4 + 1368.0 * c5() + 384.0 * c6());
            complex<double> G0 = -3.0 / 8.0 * ((ccbar_resonance ? LongDistance::g_had_ccbar(s, m_c) : CharmLoops::h(mu, s)) + 4.0 / 9.0);
            complex<double> Gb = -3.0 / 8.0 * (CharmLoops::h(mu, s, m_b) + 4.0 / 9.0);

            complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
            if (cp_conjugate)
            {
                lambda_hat_u = conj(lambda_hat_u);
            }

            complex<double> lo = c_b * Gb + c_0 * G0 + c;
            complex<double> nlo_alpha_s = memoise(c9eff_nlo_alpha_s, mu(), s, m_b, c1(), c2(), c8());
            complex<double> nlo_mc = m_c * m_c / s * 8 * ((4.0/9.0 * c1() + 1.0/3.0 * c2()) * (1.0 + lambda_hat_u) + 2.0 * c3() + 20.0 * c5());

            complex<double> result = c9() + lo;
            if (! ccbar_resonance)
                result += (model->alpha_s(mu) / (4.0 * M_PI)) * nlo_alpha_s + nlo_mc;

            return result;
        }

        complex<double> c10() const
        {
            double sign = (cp_conjugate ? -1.0 : +1.0);

            return abs_c10() * complex<double>(cos(arg_c10()), sign * sin(arg_c10()));
        }

        double rho_1(const double & s) const
        {
            return std::norm(c9eff(s) + kappa() * (2.0 * m_b_MSbar * m_B / s) * c7eff(s)) + std::norm(c10());
        }

        double rho_2(const double & s) const
        {
            return real((c9eff(s) + kappa() * (2.0 * m_b_MSbar * m_B / s) * c7eff(s)) * conj(c10()));
        }

        complex<double> rho_L(const double & s) const
        {
            return c9eff(s) + kappa() * (2.0 * m_b_MSbar * m_B / s) * c7eff(s) - c10();
        }

        complex<double> rho_R(const double & s) const
        {
            return c9eff(s) + kappa() * (2.0 * m_b_MSbar * m_B / s) * c7eff(s) + c10();
        }

        double kappa() const
        {
            // cf. [BHvD2010], Eq. (?), p. ?
            // Use m_b_MSbar(m_b_MSbar) instead m_b_MSbar(mu), as we want kappa_1 up to NLO only.
            return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / m_b_MSbar));
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(const double & s) const
        {
            return s / m_B / m_B;
        }

        // Amplitudes
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstarhat = m_Kstar / m_B;
            double m_Kstarhat2 = std::pow(m_Kstarhat, 2);
            double s_hat = s / m_B / m_B;
            double a_1 = form_factors->a_1(s), a_2 = form_factors->a_2(s);

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_long_left + (1.0 + h) / 2.0 * uncertainty_long_right;
            complex<double> prefactor = complex<double>(0.0, -1.0) * m_B();
            complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime())
                + uncertainty_isgur_wise_long * kappa() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            double formfactor = ((1.0 - m_Kstarhat2 - s_hat) * std::pow(1.0 + m_Kstarhat, 2) * a_1 - lambda(1.0, m_Kstarhat2, s_hat) * a_2)
                / (2.0 * m_Kstarhat * (1.0 + m_Kstarhat) * std::sqrt(s_hat));

            return this->norm(s) * uncertainty * prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstarhat = m_Kstar / m_B;
            double m_Kstarhat2 = std::pow(m_Kstarhat, 2);

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_perp_left + (1.0 + h) / 2.0 * uncertainty_perp_right;
            complex<double> prefactor = complex<double>(0.0, 1.0) * m_B();
            complex<double> wilson = ((c9eff(s) + c9prime()) + h * (c10() + c10prime()))
                + uncertainty_isgur_wise_perp * kappa() * (c7eff(s) + c7prime()) * (2 * m_b_MSbar * m_B / s);
            double formfactor = std::sqrt(2 * lambda(1.0, m_Kstarhat2, s_hat(s))) / (1.0 + m_Kstarhat)
                * form_factors->v(s);

            return this->norm(s) * uncertainty * prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstarhat = m_Kstar / m_B;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_par_left + (1.0 + h) / 2.0 * uncertainty_par_right;
            complex<double> prefactor = complex<double>(0.0, -1.0) * m_B();
            complex<double> wilson = ((c9eff(s) - c9prime()) + h * (c10() - c10prime()))
                + uncertainty_isgur_wise_par * kappa() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            complex<double> formfactor = std::sqrt(2) * (1.0 + m_Kstarhat) * form_factors->a_1(s);

            return this->norm(s) * uncertainty * prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        // Quantity Y = Y_9 + lambda_u_hat Y_9^u + kappa_hat Y_7, the strong phase contributor of the amplitudes
        complex<double> Y(const double & s) const
        {
            return (c9eff(s) - c9()) + kappa() * (c7eff(s) - c7()) * (2.0 * m_b_MSbar * m_B / s);
        }

        // Unormalized combinations of transversity amplitudes
        double u_1(const double & s) const
        {
            return std::norm(a_long(left_handed, s)) + std::norm(a_long(right_handed, s));
        }

        double u_2(const double & s) const
        {
            return std::norm(a_perp(left_handed, s)) + std::norm(a_perp(right_handed, s));
        }

        double u_3(const double & s) const
        {
            return std::norm(a_par(left_handed, s)) + std::norm(a_par(right_handed, s));
        }

        double u_4(const double & s) const
        {
            return real(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        }

        double u_5(const double & s) const
        {
            return (real(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s)));
        }

        double u_6(const double & s) const
        {
            return (real(a_par(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_par(right_handed, s)) * a_perp(right_handed, s)));
        }

        double u_7(const double & s) const
        {
            return imag(a_long(left_handed, s) * conj(a_par(left_handed, s)) + conj(a_long(right_handed, s)) * a_par(right_handed, s));
        }

        double u_8(const double & s) const
        {
            return (imag(a_long(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_long(right_handed, s)) * a_perp(right_handed, s)));
        }

        double u_9(const double & s) const
        {
            return (imag(a_par(left_handed, s) * conj(a_perp(left_handed, s)) - conj(a_par(right_handed, s)) * a_perp(right_handed, s)));
        }

        // Components of observables
        double decay_width(const double & s) const
        {
            return (u_1(s) + u_2(s) + u_3(s));
        }

        double a_fb_numerator(const double & s) const
        {
            return 1.5 * (real(a_par(left_handed, s) * conj(a_perp(left_handed, s))) - real(a_par(right_handed, s) * conj(a_perp(right_handed, s))));
        }

        double f_l_numerator(const double & s) const
        {
            return u_1(s);
        }

        double a_t_2_numerator(const double & s) const
        {
            return (u_2(s) - u_3(s));
        }

        double a_t_2_denominator(const double & s) const
        {
            return (u_2(s) + u_3(s));
        }

        double a_t_3_numerator(const double & s) const
        {
            return sqrt(pow(u_4(s), 2) + pow(u_7(s), 2));
        }

        double a_t_3_denominator(const double & s) const
        {
            return sqrt(u_1(s) * u_2(s));
        }

        double a_t_4_numerator(const double & s) const
        {
            return sqrt(pow(u_5(s), 2) + pow(u_8(s), 2));
        }

        double a_t_4_denominator(const double & s) const
        {
            return u_4(s);
        }

        double h_1_numerator(const double & s) const
        {
            return u_4(s);
        }

        double h_1_denominator(const double & s) const
        {
            return sqrt(u_1(s) * u_3(s));
        }

        double h_2_numerator(const double & s) const
        {
            return u_5(s);
        }

        double h_2_denominator(const double & s) const
        {
            return sqrt(u_1(s) * u_2(s));
        }

        double h_3_numerator(const double & s) const
        {
            return u_6(s);
        }

        double h_3_denominator(const double & s) const
        {
            return sqrt(u_2(s) * u_3(s));
        }
    };

    BToKstarDilepton<LowRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>(new Implementation<BToKstarDilepton<LowRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LowRecoil>::~BToKstarDilepton()
    {
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    complex<double>
    BToKstarDilepton<LowRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LowRecoil>::real_y(const double & s) const
    {
        return real(_imp->Y(s));
    }

    double
    BToKstarDilepton<LowRecoil>::imag_y(const double & s) const
    {
        return imag(_imp->Y(s));
    }

    double
    BToKstarDilepton<LowRecoil>::real_c9eff(const double & s) const
    {
        return real(_imp->c9eff(s));
    }

    double
    BToKstarDilepton<LowRecoil>::imag_c9eff(const double & s) const
    {
        return imag(_imp->c9eff(s));
    }

    double
    BToKstarDilepton<LowRecoil>::rho_1(const double & s) const
    {
        return _imp->rho_1(s);
    }

    double
    BToKstarDilepton<LowRecoil>::rho_2(const double & s) const
    {
        return _imp->rho_2(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) / Gamma;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_decay_width(const double & s) const
    {
        return _imp->decay_width(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return _imp->a_fb_numerator(s) / _imp->decay_width(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        return _imp->a_t_2_numerator(s) / _imp->a_t_2_denominator(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        return _imp->a_t_3_numerator(s) / _imp->a_t_3_denominator(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        return _imp->a_t_4_numerator(s) / _imp->a_t_4_denominator(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_h_1(const double & s) const
    {
        return _imp->h_1_numerator(s) / _imp->h_1_denominator(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_h_2(const double & s) const
    {
        return _imp->h_2_numerator(s) / _imp->h_2_denominator(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_h_3(const double & s) const
    {
        return _imp->h_3_numerator(s) / _imp->h_3_denominator(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return _imp->f_l_numerator(s) / _imp->decay_width(s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_1(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double rho_1 = _imp->rho_1(s);
        _imp->cp_conjugate = true;
        double rho_1_bar = _imp->rho_1(s);

        return (rho_1 - rho_1_bar) / (rho_1 + rho_1_bar);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_2(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double rho_1 = _imp->rho_1(s), rho_2 = _imp->rho_2(s);
        _imp->cp_conjugate = true;
        double rho_1_bar = _imp->rho_1(s), rho_2_bar = _imp->rho_2(s);

        return (rho_2 / rho_1 - rho_2_bar / rho_1_bar) / (rho_2 / rho_1 + rho_2_bar / rho_1_bar);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_3(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double rho_1 = _imp->rho_1(s), rho_2 = _imp->rho_2(s);
        _imp->cp_conjugate = true;
        double rho_1_bar = _imp->rho_1(s), rho_2_bar = _imp->rho_2(s);

        return 2.0 * (rho_2 - rho_2_bar) / (rho_1 + rho_1_bar);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_mix(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double rho_1 = _imp->rho_1(s), rho_2 = _imp->rho_2(s);

        complex<double> rho_L = _imp->rho_L(s), rho_R = _imp->rho_R(s);
        _imp->cp_conjugate = true;
        complex<double> rho_L_bar = _imp->rho_L(s), rho_R_bar = _imp->rho_R(s);

        double abs2_xi_L = norm(rho_L / rho_L_bar), abs2_xi_R = norm(rho_R / rho_R_bar);

        return (2 * rho_2 * (abs2_xi_L + abs2_xi_R - 2.0) + rho_1 * (abs2_xi_R - abs2_xi_L))
            / (rho_1 * (abs2_xi_L + abs2_xi_R + 2.0) + 2 * rho_2 * (abs2_xi_R - abs2_xi_L));
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&BToKstarDilepton<LowRecoil>::differential_branching_ratio,
                this, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::a_fb_numerator, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::decay_width, _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_unnormalized_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // Convert from asymmetry in the decay width to asymmetry in the BR
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);
        std::function<double (const double &)> integrand = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::a_fb_numerator, _imp, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / Gamma;
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::f_l_numerator, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::decay_width, _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::a_t_2_numerator, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::a_t_2_denominator, _imp, std::placeholders::_1);

        return integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_4, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_1, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom2 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_2, _imp, std::placeholders::_1);

        return integrate(num1, 64, s_min, s_max)
            / sqrt(integrate(denom1, 64, s_min, s_max) * integrate(denom2, 64, s_min, s_max));
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_3_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num1 = std::bind(
                &Implementation<BToKstarDilepton<LowRecoil>>::u_5, _imp, std::placeholders::_1);

        std::function<double (const double &)> denom1 = std::bind(
                &Implementation<BToKstarDilepton<LowRecoil>>::u_4, _imp, std::placeholders::_1);

        std::function<double (const double &)> denom2 = std::bind(
                &Implementation<BToKstarDilepton<LowRecoil>>::u_7, _imp, std::placeholders::_1);

        return integrate(num1, 64, s_min, s_max)
            / sqrt(pow(integrate(denom1, 64, s_min, s_max), 2) + pow(integrate(denom2, 64, s_min, s_max), 2));
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_4_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_h_1(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_4, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_1, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom2 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_3, _imp, std::placeholders::_1);

        return integrate(num1, 64, s_min, s_max)
            / sqrt(integrate(denom1, 64, s_min, s_max) * integrate(denom2, 64, s_min, s_max));
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_h_1_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_h_1, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_h_2(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_5, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_1, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom2 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_2, _imp, std::placeholders::_1);

        return integrate(num1, 64, s_min, s_max)
            / sqrt(integrate(denom1, 64, s_min, s_max) * integrate(denom2, 64, s_min, s_max));
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_h_2_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_h_2, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_h_3(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_6, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_2, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom2 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_3, _imp, std::placeholders::_1);

        return integrate(num1, 64, s_min, s_max)
            / sqrt(integrate(denom1, 64, s_min, s_max) * integrate(denom2, 64, s_min, s_max));
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_h_3_naive(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(&BToKstarDilepton<LowRecoil>::differential_h_3, this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max) / (s_max - s_min);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_1(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        std::function<double (const double &)> integrand = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::decay_width, _imp, std::placeholders::_1);

        double gamma = integrate(integrand, 64, s_min, s_max);
        _imp->cp_conjugate = true;
        double gamma_bar = integrate(integrand, 64, s_min, s_max);

        return (gamma - gamma_bar) / (gamma + gamma_bar);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_2(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        std::function<double (const double &)> num = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::a_fb_numerator, _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::decay_width, _imp, std::placeholders::_1);

        double a_fb = integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);
        _imp->cp_conjugate = true;
        double a_fb_bar = integrate(num, 64, s_min, s_max) / integrate(denom, 64, s_min, s_max);

        // Note that in the code A_FB does not flip its sign under CP. Therefore a_fb_bar -> -a_fb_bar here.
        return (a_fb - a_fb_bar) / (a_fb + a_fb_bar);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_3(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        std::function<double (const double &)> u_1 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_1, _imp, std::placeholders::_1);
        std::function<double (const double &)> u_2 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_2, _imp, std::placeholders::_1);
        std::function<double (const double &)> u_5 = std::bind(&Implementation<BToKstarDilepton<LowRecoil>>::u_5, _imp, std::placeholders::_1);

        double iu_1 = integrate(u_1, 16, s_min, s_max), iu_2 = integrate(u_2, 16, s_min, s_max), iu_5 = integrate(u_5, 16, s_min, s_max);
        _imp->cp_conjugate = true;
        double iu_1_bar = integrate(u_1, 16, s_min, s_max), iu_2_bar = integrate(u_2, 16, s_min, s_max), iu_5_bar = integrate(u_5, 16, s_min, s_max);

        // Note that in the code U_5 does not flip its sign under CP. Therefore iu_5_bar -> -iu_5_bar here.
        return (iu_5 - iu_5_bar) / std::sqrt((iu_1 + iu_1_bar) * (iu_2 + iu_2_bar));
    }
}
