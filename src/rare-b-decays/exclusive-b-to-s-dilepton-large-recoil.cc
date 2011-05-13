/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
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
#include <src/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <src/rare-b-decays/form-factors.hh>
#include <src/rare-b-decays/hard-scattering.hh>
#include <src/rare-b-decays/long-distance.hh>
#include <src/utils/destringify.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/memoise.hh>
#include <src/utils/model.hh>
#include <src/utils/options.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <functional>

#include <gsl/gsl_sf.h>

namespace eos
{
    struct ShortDistanceLargeRecoil
    {
        struct ParameterSet
        {
            double                    m_b_PS, m_b_pole, m_c;
            double                    m_B, m_K;

            double                    mu;
            double                    mu_f;

            double                    alpha_s_mu, alpha_s_sqrt05mu;

            double                    f_B, f_K;

            double                    c1, c2, c3, c4, c5, c6, c8;
            std::complex<double>      c7, c7prime;

            double                    e_q;

            double                    a_1, a_2;

            double                    lambda_B_p;

            ParameterSet(const double & m_b_PS_, const double & m_b_pole_, const double & m_c_,
                    const double & m_B_, const double & m_K_,
                    const double & mu_, const double & mu_f_,
                    const double & alpha_s_mu_, const double & alpha_s_sqrt05mu_,
                    const double & f_B_, const double & f_K_,
                    const double & c1_, const double & c2_, const double & c3_, const double & c4_, const double & c5_, const double & c6_, const double & c8_,
                    const std::complex<double> & c7_, const std::complex<double> & c7prime_,
                    const double & e_q_,
                    const double & a_1_, const double & a_2_,
                    const double & lambda_B_p_) :
                 m_b_PS(m_b_PS_), m_b_pole(m_b_pole_), m_c(m_c_),
                 m_B(m_B_), m_K(m_K_),
                 mu(mu_), mu_f(mu_f_),
                 alpha_s_mu(alpha_s_mu_), alpha_s_sqrt05mu(alpha_s_sqrt05mu_),
                 f_B(f_B_), f_K(f_K_),
                 c1(c1_), c2(c2_), c3(c3_), c4(c4_), c5(c5_), c6(c6_), c8(c8_),
                 c7(c7_), c7prime(c7prime_),
                 e_q(e_q_),
                 a_1(a_1_), a_2(a_2_),
                 lambda_B_p(lambda_B_p_)
            { }
        };

        static double energy(const double & s, const ParameterSet & p)
        {
            return (p.m_B * p.m_B + p.m_K * p.m_K - s) / (2.0 * p.m_B);
        }

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        static complex<double> c7eff(const ParameterSet & p)
        {
            return p.c7 - 1.0/3.0 * p.c3 - 4.0/9.0 * p.c4 - 20.0/3.0 * p.c5 - 80.0/9.0 * p.c6;
        }

        // cf. [BFS2001], below Eq. (26), p. 8
        static double c8eff(const ParameterSet & p)
        {
             return p.c8 + p.c3 - 1.0/6.0 * p.c4 + 20.0 * p.c5 - 10.0/3.0 * p.c6;
        }

        // cf. [BFS2001], Eq. (10), p. 4
        static complex<double> Y0(const double & s, const ParameterSet & p)
        {
            double Y_c = 4.0 / 3.0 * p.c1 + p.c2 + 6.0 * p.c3 + 60.0 * p.c5;
            double Y_b = -0.5 * (7.0 * p.c3 + 4.0 / 3.0 * p.c4 + 76.0 * p.c5 + 64.0 / 3.0 * p.c6);
            double Y_0 = -0.5 * (p.c3 + 4.0 / 3.0 * p.c4 + 16.0 * p.c5 + 64 / 3.0 * p.c6);
            double Y = 2.0 / 9.0 * (6.0 * p.c3 + 32.0 * p.c5 + 32.0 / 3.0 * p.c6);

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * CharmLoops::h(p.mu, s, p.m_c)
                 + Y_b * CharmLoops::h(p.mu, s, p.m_b_PS)
                 + Y_0 * CharmLoops::h(p.mu, s)
                 + Y;
        }

        /* NLO functions */
        // cf. [BFS2001], Eq. (48)
        static double phi_K(const double & u, const ParameterSet & p)
        {
            double xi = 2.0 * u - 1.0;

            return 6.0 * u * (1 - u) * (1.0 + p.a_1 * 3.0 * xi + p.a_2 * (7.5 * xi * xi - 1.5));
        }

        static complex<double> t_perp_0(const double & u, const ParameterSet & p, const double & m_q)
        {
            double ubar = 1.0 - u;
            double m_q2 = m_q * m_q, m_B2 = p.m_B * p.m_B;
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

        // cf. [BFS2001], Eq. (27), p. 8
        static complex<double> t_perp(const double & s, const double & u, const ParameterSet & p, const double & m_q)
        {
            if (0.0 == s)
                return t_perp_0(u, p, m_q);

            double ubar = 1.0 - u;
            double x = ubar * p.m_B * p.m_B + u * s;
            double E = energy(s, p);

            complex<double> result = (2.0 * p.m_B / ubar / E) * memoise(HardScattering::I1, s, u, m_q, p.m_B);
            if (m_q > 0.0)
                result = result + (s / ubar / ubar / E / E) * (CharmLoops::B0(x, m_q) - CharmLoops::B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (28), p. 8
        static complex<double> t_par(const double & s, const double & u, const ParameterSet & p, const double & m_q)
        {
            double ubar = 1.0 - u;
            double x = ubar * p.m_B * p.m_B + u * s;
            double E = energy(s, p);

            complex<double> result = (2.0 * p.m_B / ubar / E) * memoise(HardScattering::I1, s, u, m_q, p.m_B);
            if (m_q > 0.0)
                result = result + (x / ubar / ubar / E / E) * (CharmLoops::B0(x, m_q) - CharmLoops::B0(s, m_q));

            return result;
        }

        // cf. [BFS2001], Eq. (36), p. 9
        static double L(const double & s, const ParameterSet & p)
        {
            double m_b_PS2 = p.m_b_PS * p.m_b_PS;

            return -1.0 * (m_b_PS2 - s) / s * std::log(1.0 - s / m_b_PS2);
        }

        // cf. [BFS2001], Eq. (54), p. 15
        static complex<double> lambda_B_m_inv(const double & s, const ParameterSet & p)
        {
            if (0.0 == s)
                return complex<double>(0.0, 0.0);

            double omega_0 = p.lambda_B_p;
            double x = s / p.m_B / omega_0;
            double ei = gsl_sf_expint_Ei(x);

            complex<double> result = complex<double>(-ei, M_PI) * (std::exp(-x) / omega_0);

            return result;
        }

        // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
        static complex<double> C0_perp(const double & h, const double & s, const ParameterSet & p)
        {
            return (c7eff(p) + h * p.c7prime) + s / (2.0 * p.m_b_PS * p.m_B) * Y0(s, p);
        }

        // cf. [BFS2001], Eqs. (34), (37), p. 9
        static complex<double> C1f_perp(const double & h, const double & s, const ParameterSet & p)
        {
            // cf. [BFS2004], Eq. (44), p. 24
            // [Christoph] Use c7 instead of c7eff
            return (p.c7 + h * p.c7prime) * (8.0 * std::log(p.m_b_PS / p.mu) - L(s, p) - 4.0 * (1.0 - p.mu_f / p.m_b_PS));
        }

        // cf. [BFS2001], Eqs. (34), (37), p. 9
        static complex<double> C1nf_perp(const double &, const double & s, const ParameterSet & p)
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = p.mu * p.m_b_pole / p.m_b_PS;

            // cf. [BFS2001], Eq. (37), p. 9
            // [Christoph] Use c8 instead of c8eff
            return (-1.0 / QCD::casimir_f) * (
                    (p.c2 - p.c1 / 6.0) * memoise(CharmLoops::F27_massive, mu_pole, s, p.m_b_PS, p.m_c) + p.c8 * CharmLoops::F87_massless(mu_pole, s, p.m_b_PS)
                    + (s / (2.0 * p.m_b_PS * p.m_B)) * (
                        p.c1 * memoise(CharmLoops::F19_massive, mu_pole, s, p.m_b_PS, p.m_c)
                        + p.c2 * memoise(CharmLoops::F29_massive, mu_pole, s, p.m_b_PS, p.m_c)
                        + p.c8 * CharmLoops::F89_massless(s, p.m_b_PS)));
        }

        // cf. [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
        static complex<double> C0_par(const double & s, const ParameterSet & p)
        {
            return -1.0 * (c7eff(p) - p.c7prime + p.m_B / (2.0 * p.m_b_PS) * Y0(s, p));
        }

        // cf. [BFS2001], Eqs. (38), p. 9
        static complex<double> C1nf_par(const double & s, const ParameterSet & p)
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = p.mu * p.m_b_pole / p.m_b_PS;

            // cf. [BFS2001], Eq. (38), p. 9
            // [Christoph] Use c8 instead of c8eff.
            return (+1.0 / QCD::casimir_f) * (
                    (p.c2 - p.c1 / 6.0) * memoise(CharmLoops::F27_massive, mu_pole, s, p.m_b_PS, p.m_c) + c8eff(p) * CharmLoops::F87_massless(mu_pole, s, p.m_b_PS)
                    + (p.m_B / (2.0 * p.m_b_PS)) * (
                        p.c1 * memoise(CharmLoops::F19_massive, mu_pole, s, p.m_b_PS, p.m_c)
                        + p.c2 * memoise(CharmLoops::F29_massive, mu_pole, s, p.m_b_PS, p.m_c)
                        + c8eff(p) * CharmLoops::F89_massless(s, p.m_b_PS)));
        }

        // cf. [BFS2001], Eqs. (35), (38), p. 9
        static complex<double> C1f_par(const double & s, const ParameterSet & p)
        {
            // cf. [BFS2004], Eq. (45), p. 24
            // [Christoph] Use c7 instead of c7eff.
            return -1.0 * (p.c7 - p.c7prime) * (8.0 * std::log(p.m_b_PS / p.mu) + 2.0 * L(s, p) - 4.0 * (1.0 - p.mu_f / p.m_b_PS));
            /* for [BFS2001] version of xi_par we also needed: */
            // C_par_f += (m_B / (2.0 * m_b)) * Y0(s) * (2.0 - 2.0 * L(s));
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> C0_pseudo(const double & s, const ParameterSet & p)
        {
            return -C0_par(s, p);
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> C1nf_pseudo(const double & s, const ParameterSet & p)
        {
            return -C1nf_par(s, p);
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> C1f_pseudo(const double & s, const ParameterSet & p)
        {
            // the correct sign in front of C_7^eff is plus, as one can see by
            // comparison with [BF2001], Eq. (63)
            return (p.c7 - p.c7prime) * (8.0 * log(p.m_b_PS / p.mu_f) + 2.0 * L(s, p) - 4.0 + 4.0 * p.mu_f / p.m_b_PS);
        }

        // cf. [BFS2001], Eq. (17), p. 6
        static complex<double> T0_par_p(const double &, const ParameterSet &)
        {
            return 0.0;
        }

        // cf. [BFS2001], Eq. (18), p. 6 with \omega integrated out.
        static complex<double> T0_par_m(const double &, const ParameterSet & p)
        {
            return -p.e_q * 4.0 * p.m_B / p.m_b_PS * (p.c3 + 4.0/3.0 * p.c4 + 16.0 * p.c5 + 64.0/3.0 * p.c6);
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> T0_pseudo_p(const double & s, const ParameterSet & p)
        {
            return -T0_par_p(s, p);
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> T0_pseudo_m(const double & s, const ParameterSet & p)
        {
            return -T0_par_m(s, p);
        }

        // cf. [BFS2001], Eqs. (16), (21), (25), pp. 5-7
        static complex<double> T1f_perp_p(const double & h, const double & s, const double & u, const ParameterSet & p)
        {
            // cf. [BFS2001], Eq. (20)
            // [Christoph] Use c7 instead of c7eff
            return (p.c7 + h * p.c7prime) * (2.0 * p.m_B / (1.0 - u) / energy(s, p));
        }

        static complex<double> T1nf_perp_p(const double &, const double & s, const double & u, const ParameterSet & p)
        {
            static const double e_d = -1.0/3.0;
            static const double e_u = +2.0/3.0;

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double s_hat = s / p.m_B / p.m_B;

            // cf. [BFS2001], Eq. (23)
            // [Christoph] Use c8 instead of c8eff
            return -4.0 * e_d * p.c8 / (u + (1.0 - u) * s_hat)
                + p.m_B / (2.0 * p.m_b_PS) * (
                        e_u * (-p.c1 / 6.0 + p.c2 + 6.0 * p.c6) * t_perp(s, u, p, p.m_c)
                        + e_d * (p.c3 - p.c4 / 6.0 + 16.0 * p.c5 + 10.0/3.0 * p.c6 - (4.0 * p.m_b_PS / p.m_B) * (p.c3 - p.c4/6.0 + 4.0 * p.c5 - 2.0/3.0 * p.c6))
                        + e_d * (p.c3 - p.c4 / 6.0 + 16.0 * p.c5 - 8.0/3.0 * p.c6) * t_perp(s, u, p, 0.0));
        }

        // cf. [BFS2001], Eqs. (16), (21), (25), pp. 5-7
        static complex<double> T1f_par_p(const double & s, const double & u, const ParameterSet & p)
        {
            // cf. [BFS2001], Eq. (21)
            //complex<double> Tf_par_p = (c7() - c7prime + (s / (2.0 * m_b * m_B)) * Y0(s)) * (2.0 * pow(m_B / energy(s), 2));
            // cf. [BFS2004], Eq. (49)
            // [Christoph] Use c7 instead of c7eff.
            return (p.c7 - p.c7prime) * (4.0 * p.m_B / (1.0 - u) / energy(s, p));
        }

        // cf. [BFS2001], Eq. (16), (22), (26), pp. 5-8
        static complex<double> T1f_par_m(const double &, const double &, const ParameterSet &)
        {
            return 0.0;
        }

        // cf. [BFS2001], Eq. (25), p. 7
        static complex<double> T1nf_par_p(const double & s, const double & u, const ParameterSet & p)
        {
            static const double e_d = -1.0/3.0;
            static const double e_u = +2.0/3.0;

            // cf. [BFS2001], Eq. (25)
            return p.m_B / p.m_b_PS * (
                    e_u * (-p.c1 / 6.0 + p.c2 + 6.0 * p.c6) * t_par(s, u, p, p.m_c)
                    + e_d * (p.c3 - p.c1 / 6.0 + 16.0 * p.c5 + 10.0/3.0 * p.c6) * t_par(s, u, p, p.m_b_PS)
                    + e_d * (p.c3 - p.c4 / 6.0 + 16.0 * p.c5 -  8.0/3.0 * p.c6) * t_par(s, u, p, 0.0));
        }

        // cf. [BFS2001], Eq. (26), pp. 7-8 with \omega integrated out.
        static complex<double> T1nf_par_m(const double & s, const double & u, const ParameterSet & p)
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = p.mu * p.m_b_pole / p.m_b_PS;

            double s_hat = s / p.m_B / p.m_B;
            double ubar = 1.0 - u;
            double x = ubar * p.m_B * p.m_B + u * s;

            // [Christoph] Use c8 instead of c8eff.
            return p.e_q * (8.0 / (ubar + u * s_hat) * p.c8
                    + 6.0 * p.m_B / p.m_b_PS * (
                        (-p.c1 / 6.0 + p.c2 + p.c4 + 10 * p.c6) * CharmLoops::h(mu_pole, x, p.m_c)
                        + (p.c3 + 5.0/6.0 * p.c4 + 16.0 * p.c5 + 22.0/3.0 * p.c6) * CharmLoops::h(mu_pole, x, p.m_b_PS)
                        + (p.c3 + 17.0/6.0 * p.c4 + 16.0 * p.c5 + 82.0/3.0 * p.c6) * CharmLoops::h(mu_pole, x)
                        -8.0 / 27.0 * (-7.5 * p.c4 + 12.0 * p.c5 - 32.0 * p.c6)));
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> T1f_pseudo_p(const double & s, const double & u, const ParameterSet & p)
        {
            return -T1f_par_p(s, u, p);
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> T1f_pseudo_m(const double &, const double &, const ParameterSet &)
        {
            return 0.0;
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> T1nf_pseudo_p(const double & s, const double & u, const ParameterSet & p)
        {
            return -T1nf_par_p(s, u, p);
        }

        // cf. [BHP2007], Eq. (B.2)
        static complex<double> T1nf_pseudo_m(const double & s, const double & u, const ParameterSet & p)
        {
            return -T1nf_par_m(s, u, p);
        }

        // cf. [BFS2001], Eq. (16) times phi_K^*_perp
        static complex<double> T_perp_sum(const double & h, const double & s, const double & u, const ParameterSet & p)
        {
             double a = p.alpha_s_sqrt05mu * QCD::casimir_f / 4.0 / M_PI;

             complex<double> result = 1.0 / p.lambda_B_p * a * (T1f_perp_p(h, s, u, p) + T1nf_perp_p(h, s, u, p));

             return phi_K(u, p) * result;

             // TODO: Hard scattering + Weak annihilation from [BFS2004], Eqs. (51), (52)
        }

        // cf. [BFS2001], Eq. (16) times phi_K^*_par
        static complex<double> T_par_sum(const double & s, const double & u, const ParameterSet & p)
        {
            double a = p.alpha_s_sqrt05mu * QCD::casimir_f / 4.0 / M_PI;

            complex<double> result = 1.0 / p.lambda_B_p * (T0_par_p(s, p) + a * (T1f_par_p(s, u, p) + T1nf_par_p(s, u, p))) +
                                   lambda_B_m_inv(s, p) * (T0_par_m(s, p) + a * (T1f_par_m(s, u, p) + T1nf_par_m(s, u, p)));

            return phi_K(u, p) * result;
        }

        // cf. [BHP2007], Eq. (B.1), p. 25
        static complex<double> T_pseudo_sum(const double & s, const double & u, const ParameterSet & p)
        {
            double a = p.alpha_s_sqrt05mu * QCD::casimir_f / 4.0 / M_PI;

            complex<double> result = 1.0 / p.lambda_B_p * (T0_pseudo_p(s, p) + a * (T1f_pseudo_p(s, u, p) + T1nf_pseudo_p(s, u, p))) +
                                   lambda_B_m_inv(s, p) * (T0_pseudo_m(s, p) + a * (T1f_pseudo_m(s, u, p) + T1nf_pseudo_m(s, u, p)));

            return phi_K(u, p) * result;
        }

        // cf. [BFS2001], Eq. (15) with a = perp
        static complex<double> calT_perp(const double & h, const double & s, const ParameterSet & p, const double & xi_perp)
        {
            std::complex<double> result = xi_perp * (C0_perp(h, s, p) + p.alpha_s_mu * QCD::casimir_f / 4.0 / M_PI * (C1f_perp(h, s, p) + C1nf_perp(h, s, p)));

           // to do: bind umschreiben auf static
            result += power_of<2>(M_PI) / 3.0 * (p.f_B * p.f_K) / p.m_B *
                    integrate(std::function<complex<double> (const double &)>(
                                std::bind(&ShortDistanceLargeRecoil::T_perp_sum, h, s, std::placeholders::_1, p)),
                        64, 0.001, 0.999);

            return result;
        }

        // cf. [BFS2001], Eq. (15) with a = par, and [BHP2008], Eq. (C.4)
        static complex<double> calT_par(const double & s, const ParameterSet & p, const double & xi_par)
        {
            std::complex<double> result = xi_par * (C0_par(s, p) + p.alpha_s_mu * QCD::casimir_f / 4.0 / M_PI * (C1f_par(s, p) + C1nf_par(s, p)));

            result += power_of<2>(M_PI) / 3.0 * (p.f_B * p.f_K) / p.m_B * (p.m_K / energy(s, p)) *
                    integrate(std::function<complex<double> (const double &)>(
                                std::bind(&ShortDistanceLargeRecoil::T_par_sum, s, std::placeholders::_1, p)),
                        32, 0.001, 0.999);

            return result;
        }

        // cf. [BHP2007], Eq. (B.1), p. 25
        static complex<double> calT_pseudo(const double & s, const ParameterSet & p, const double & xi_pseudo)
        {
            std::complex<double> result = xi_pseudo * (C0_pseudo(s, p) + p.alpha_s_mu * QCD::casimir_f / 4.0 / M_PI * (C1f_pseudo(s, p) + C1nf_pseudo(s, p)));

            // integration of omega is included in T_pseudo_sum as lambda_B_*_inv
            result += power_of<2>(M_PI) / 3.0 * (p.f_B * p.f_K) / p.m_B *
                integrate(std::function<complex<double> (const double &)>(
                               std::bind(&ShortDistanceLargeRecoil::T_pseudo_sum, s, std::placeholders::_1, p)),
                        32, 0.001, 0.999);

            return result;
        }
    };

    using std::norm;

    /*
     * Decay: B -> K^* l lbar at Large Recoil, cf. [BHP2008]
     */
    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        Parameter hbar;

        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter abs_c7;

        Parameter arg_c7;

        Parameter abs_c7prime;

        Parameter arg_c7prime;

        Parameter c8;

        Parameter abs_c9;

        Parameter arg_c9;

        Parameter abs_c9prime;

        Parameter arg_c9prime;

        Parameter abs_c10;

        Parameter arg_c10;

        Parameter abs_c10prime;

        Parameter arg_c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        double m_l;

        Parameter mu;

        Parameter alpha_e;

        Parameter g_fermi;

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

        Parameter tau;

        double e_q;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Implementation(const Parameters & p, const Options & o) :
            model(Model::make("SM", p)),
            hbar(p["hbar"]),
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            abs_c7(p["Abs{c7}"]),
            arg_c7(p["Arg{c7}"]),
            abs_c7prime(p["Abs{c7'}"]),
            arg_c7prime(p["Arg{c7'}"]),
            c8(p["c8"]),
            abs_c9(p["Abs{c9}"]),
            arg_c9(p["Arg{c9}"]),
            abs_c9prime(p["Abs{c9'}"]),
            arg_c9prime(p["Arg{c9'}"]),
            abs_c10(p["Abs{c10}"]),
            arg_c10(p["Arg{c10}"]),
            abs_c10prime(p["Abs{c10'}"]),
            arg_c10prime(p["Arg{c10'}"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B_" + o.get("q", "d")]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            alpha_e(p["QED::alpha_e(m_b)"]),
            g_fermi(p["G_Fermi"]),
            f_B(p["decay-constant::B_" + o.get("q", "d")]),
            f_Kstar_par(p["B->K^*::f_Kstar_par"]),
            f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"]),
            lambda_B_p(p["lambda_B_p"]),
            a_1_par(p["B->K^*::a_1_par"]),
            a_2_par(p["B->K^*::a_2_par"]),
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
            tau(p["life_time::B_" + o.get("q", "d")]),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "BZ2004"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            std::string spectator_quark = o.get("q", "d");
            if (spectator_quark == "d")
            {
                e_q = -1.0 / 3.0;
            }
            else if (spectator_quark == "u")
            {
                e_q = 2.0 / 3.0;
            }
            else
                throw InternalError("Unsupported spectator quark");

            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.0;//m_l = 0.10565836; // (GeV), cf. [PDG2008], p. 13
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
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

            return std::sqrt(power_of<2>(g_fermi() * alpha_e()) / 3.0 / 1024 / power_of<5>(M_PI) / m_B()
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lam(s))); // cf. [BHP2008], Eq. (C.6), p. 21
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

        double energy(const double & s) const
        {
            return (m_B() * m_B() + m_Kstar() * m_Kstar() - s) / (2.0 * m_B());
        }

        // alias c7('),c9(') and c10(')
        inline complex<double> c7() const { return std::polar(abs_c7(), (cp_conjugate ? -1.0 : +1.0) * arg_c7()); }
        inline complex<double> c7prime() const { return std::polar(abs_c7prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c7prime()); }
        inline complex<double> c9() const { return std::polar(abs_c9(), (cp_conjugate ? -1.0 : +1.0) * arg_c9()); }
        inline complex<double> c9prime() const { return std::polar(abs_c9prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c9prime()); }
        inline complex<double> c10() const { return std::polar(abs_c10(), (cp_conjugate ? -1.0 : +1.0) * arg_c10()); }
        inline complex<double> c10prime() const { return std::polar(abs_c10prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c10prime()); }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double m_b = m_b_PS();

            double h = helicity;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_long_left + (1.0 + h) / 2.0 * uncertainty_long_right;
            complex<double> wilson = (c9() - c9prime()) + h * (c10() - c10prime());
            double prefactor = -1.0 / (2.0 * m_Kstar() * std::sqrt(s));

            ShortDistanceLargeRecoil::ParameterSet p(m_b, model->m_b_pole(), m_c(),
                                                     m_B(), m_Kstar(),
                                                     mu(), mu_f(),
                                                     model->alpha_s(mu()), model->alpha_s(sqrt(mu() * 0.5)),
                                                     f_B(), f_Kstar_perp(),
                                                     c1(), c2(), c3(), c4(), c5(), c6(), c8(),
                                                     c7(), c7prime(),
                                                     e_q,
                                                     a_1_perp(), a_2_perp(),
                                                     lambda_B_p());

            std::complex<double> calTperp = ShortDistanceLargeRecoil::calT_perp(-1.0, s, p, xi_perp(s));

            p.f_K        = f_Kstar_par;
            p.a_1        = a_1_par();
            p.a_2        = a_2_par();
            std::complex<double> calTpar = ShortDistanceLargeRecoil::calT_par(s, p, xi_par(s));

            complex<double> a = wilson * ((m_B() * m_B() - m_Kstar() * m_Kstar() - s) * 2.0 * energy(s) * xi_perp(s)
                                - lam(s) * m_B() / (m_B() * m_B() - m_Kstar() * m_Kstar()) * (xi_perp(s) - xi_par(s)));
            complex<double> b = 2.0 * m_b_PS() * (((m_B() * m_B() + 3 * m_Kstar() * m_Kstar() - s) * 2.0 * energy(s) / m_B()
                                - lam(s) / (m_B() * m_B() - m_Kstar() * m_Kstar())) * calTperp
                                - lam(s) / (m_B() * m_B() - m_Kstar() * m_Kstar()) * calTpar);

            return this->norm(s) * uncertainty * prefactor * (a + b);
        }

        // cf. [BHP2008], p. 20
        complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double m_b = m_b_PS();

            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double mKhat = m_Kstar / m_B;

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_perp_left + (1.0 + h) / 2.0 * uncertainty_perp_right;
            double prefactor = +std::sqrt(2.0) * m_B() * std::sqrt(lambda(1.0, mKhat * mKhat, shat));
            complex<double> wilson = (c9() + c9prime()) + h * (c10() + c10prime());

            ShortDistanceLargeRecoil::ParameterSet p(m_b, model->m_b_pole(), m_c(),
                                                     m_B(), m_Kstar(),
                                                     mu(), mu_f(),
                                                     model->alpha_s(mu()), model->alpha_s(sqrt(mu() * 0.5)),
                                                     f_B(), f_Kstar_perp(),
                                                     c1(), c2(), c3(), c4(), c5(), c6(), c8(),
                                                     c7(), c7prime(),
                                                     e_q,
                                                     a_1_perp(), a_2_perp(),
                                                     lambda_B_p());

            return this->norm(s) * uncertainty * prefactor *
                    (wilson * xi_perp(s) + (2.0 * mbhat / shat) * ShortDistanceLargeRecoil::calT_perp(+1.0, s, p, xi_perp(s)));
        }

        // cf. [BHP2008], p. 20
        complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double m_b = m_b_PS();

            double h = helicity;
            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B();
            double mKhat = m_Kstar() / m_B();

            double uncertainty = (1.0 - h) / 2.0 * uncertainty_par_left + (1.0 + h) / 2.0 * uncertainty_par_right;
            double prefactor = -std::sqrt(2.0) * m_B() * (1.0 - shat);
            complex<double> wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            ShortDistanceLargeRecoil::ParameterSet p(m_b, model->m_b_pole(), m_c(),
                                                     m_B(), m_Kstar(),
                                                     mu(), mu_f(),
                                                     model->alpha_s(mu()), model->alpha_s(sqrt(mu() * 0.5)),
                                                     f_B(), f_Kstar_perp(),
                                                     c1(), c2(), c3(), c4(), c5(), c6(), c8(),
                                                     c7(), c7prime(),
                                                     e_q,
                                                     a_1_perp(), a_2_perp(),
                                                     lambda_B_p());

            return this->norm(s) * uncertainty * prefactor *
                    (wilson * xi_perp(s) + (2.0 * mbhat / shat) * (1.0 - mKhat * mKhat) * ShortDistanceLargeRecoil::calT_perp(-1.0, s, p, xi_perp(s)));
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

        double differential_branching_ratio(const double & s) const
        {
            return unnormalized_decay_width(s) * tau / hbar();
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

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const Options & options) :
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
        return _imp->differential_branching_ratio(s);
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

    /*
     * Decay: B -> K l lbar at Large Recoil
     */
    template <>
    struct Implementation<BToKDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        Parameter hbar;

        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter abs_c7;

        Parameter arg_c7;

        Parameter abs_c7prime;

        Parameter arg_c7prime;

        Parameter c8;

        Parameter abs_c9;

        Parameter arg_c9;

        Parameter abs_c9prime;

        Parameter arg_c9prime;

        Parameter abs_c10;

        Parameter arg_c10;

        Parameter abs_c10prime;

        Parameter arg_c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_K;

        Parameter m_e;

        Parameter m_mu;

        Parameter m_tau;

        double m_l;

        Parameter mu;

        Parameter alpha_e;

        Parameter g_fermi;

        Parameter f_B;

        Parameter f_K;

        Parameter lambda_B_p;

        Parameter a_1;

        Parameter a_2;

        // Mean life times
        Parameter tau;

        // spectator quark charge
        double e_q;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Implementation(const Parameters & p, const Options & o) :
            model(Model::make("SM", p)),
            hbar(p["hbar"]),
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            abs_c7(p["Abs{c7}"]),
            arg_c7(p["Arg{c7}"]),
            abs_c7prime(p["Abs{c7'}"]),
            arg_c7prime(p["Arg{c7'}"]),
            c8(p["c8"]),
            abs_c9(p["Abs{c9}"]),
            arg_c9(p["Arg{c9}"]),
            abs_c9prime(p["Abs{c9'}"]),
            arg_c9prime(p["Arg{c9'}"]),
            abs_c10(p["Abs{c10}"]),
            arg_c10(p["Arg{c10}"]),
            abs_c10prime(p["Abs{c10'}"]),
            arg_c10prime(p["Arg{c10'}"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_B(p["mass::B_" + o.get("q", "d")]),
            m_K(p["mass::K0"]),
            m_e(p["mass::e"]),
            m_mu(p["mass::mu"]),
            m_tau(p["mass::tau"]),
            mu(p["mu"]),
            alpha_e(p["QED::alpha_e(m_b)"]),
            g_fermi(p["G_Fermi"]),
            f_B(p["decay-constant::B_" + o.get("q", "d")]),
            f_K(p["decay-constant::K_" + o.get("q", "d")]),
            lambda_B_p(p["lambda_B_p"]),
            a_1(p["B->K::a_1@1GeV"]),
            a_2(p["B->K::a_2@1GeV"]),
            tau(p["life_time::B_" + o.get("q", "d")]),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<PToP>::create("B->K@" + o.get("form-factors", "BZ2004v2"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            std::string spectator_quark = o.get("q", "d");
            if (spectator_quark == "d")
            {
                e_q = -1.0 / 3.0;
            }
            else if (spectator_quark == "u")
            {
                e_q = 2.0 / 3.0;
            }
            else
                throw InternalError("Unsupported spectator quark");

            std::string lepton = o.get("l", "mu");
            if ("e" == lepton)
            {
                m_l = m_e;
            }
            else if ("mu" == lepton)
            {
                m_l = m_mu;
            }
            else if ("tau" == lepton)
            {
                m_l = m_tau;
            }
            else
                throw InternalError("Unknown fourth lepton generation: " + lepton);
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
            return model->m_b_ps(1.5);
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

        // alias c7('), c9(') and c10(')
        inline complex<double> c7() const { return std::polar(abs_c7(), (cp_conjugate ? -1.0 : +1.0) * arg_c7()); }
        inline complex<double> c7prime() const { return std::polar(abs_c7prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c7prime()); }
        inline complex<double> c9() const { return std::polar(abs_c9(), (cp_conjugate ? -1.0 : +1.0) * arg_c9()); }
        inline complex<double> c9prime() const { return std::polar(abs_c9prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c9prime()); }
        inline complex<double> c10() const { return std::polar(abs_c10(), (cp_conjugate ? -1.0 : +1.0) * arg_c10()); }
        inline complex<double> c10prime() const { return std::polar(abs_c10prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c10prime()); }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_A(const double &) const
        {
            return c10();
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_P(const double & s) const
        {
            return m_l * c10() *
                    ((m_B() * m_B() - m_K() * m_K()) / s * (form_factors->f_0(s) / form_factors->f_p(s) - 1.0) - 1.0);
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_V(const double & s) const
        {
            double m_b = m_b_PS();

            ShortDistanceLargeRecoil::ParameterSet p(m_b, model->m_b_pole(), m_c(),
                                                     m_B(), m_K(),
                                                     mu(), mu_f(),
                                                     model->alpha_s(mu()), model->alpha_s(sqrt(mu() * 0.5)),
                                                     f_B(), f_K(),
                                                     c1(), c2(), c3(), c4(), c5(), c6(), c8(),
                                                     c7(), c7prime(),
                                                     e_q,
                                                     a_1(), a_2(),
                                                     lambda_B_p());

            return c9() + 2.0 * m_b / m_B() * ShortDistanceLargeRecoil::calT_pseudo(s, p, xi_pseudo(s)) / xi_pseudo(s);
        }

        // cf. [BHP2007], Eqs. (4.2), (4.4), (4.5), p. 5
        double N(const double & s) const
        {
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

            return power_of<2>(g_fermi * alpha_e() * lambda_t) * std::sqrt(lam(s)) * beta_l(s) * xi_pseudo(s) * xi_pseudo(s) /
                    (512.0 * power_of<5>(M_PI) * power_of<3>(m_B()));
        }

        // cf. [BHP2007], Eq. (4.2)
        double a_l(const double & s) const
        {
            double result = s * std::norm(F_P(s));
            result += 0.25 * lam(s) * (std::norm(F_A(s)) + std::norm(F_V(s)));
            result += 2.0 * m_l * (m_B() * m_B() - m_K() * m_K() + s) * std::real(F_P(s) * std::conj(F_A(s)));
            result += 4.0 * m_l * m_l * m_B() * m_B() * std::norm(F_A(s));

            return N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.4)
        double c_l(const double & s) const
        {
            double result = N(s) * -0.25 * lam(s) * beta_l(s) * beta_l(s) * (std::norm(F_A(s)) + std::norm(F_V(s)));
            return result;
        }

        // cf. [BHP2007], Eq. (4.1)
        double unnormalized_decay_width(const double & s) const
        {
            return 2.0 * (a_l(s) + c_l(s) / 3.0);
        }

        double differential_branching_ratio(const double & s) const
        {
            return unnormalized_decay_width(s) * tau() / hbar();
        }

        double differential_flat_term_numerator(const double & s) const
        {
             return 2.0 * (a_l(s) + c_l(s));
        }
    };

    BToKDilepton<LargeRecoil>::BToKDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKDilepton<LargeRecoil>>(new Implementation<BToKDilepton<LargeRecoil>>(parameters, options))
    {
    }

    BToKDilepton<LargeRecoil>::~BToKDilepton()
    {
    }


    double
    BToKDilepton<LargeRecoil>::a_l(const double & s) const
    {
        return _imp->a_l(s);
    }

    double
    BToKDilepton<LargeRecoil>::c_l(const double & s) const
    {
        return _imp->c_l(s);
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

    // Integrated Observables
    double
    BToKDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(std::mem_fn(&BToKDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }

    double
    BToKDilepton<LargeRecoil>::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToKDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }
}
