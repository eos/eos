/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <src/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <src/rare-b-decays/form-factors.hh>
#include <src/rare-b-decays/hard-scattering.hh>
#include <src/utils/destringify.hh>
#include <src/utils/integrate.hh>
#include <src/utils/memoise.hh>
#include <src/utils/model.hh>
#include <src/utils/options.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>
#include <src/utils/save.hh>

#include <cmath>

#include <gsl/gsl_sf.h>

namespace eos
{
    template <>
    struct Implementation<BToKstarGamma>
    {
        std::shared_ptr<Model> model;

        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c8;

        Parameter abs_c7;

        Parameter arg_c7;

        Parameter abs_c7prime;

        Parameter arg_c7prime;

        Parameter a_1_perp;

        Parameter a_2_perp;

        Parameter f_B;

        Parameter f_Kstar_perp;

        Parameter lambda_B_p;

        Parameter uncertainty_xi_perp;

        Parameter m_B;

        Parameter m_Kstar;

        Parameter mu;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Implementation(const Parameters & p, const Options & o) :
            model(Model::make("SM", p)),
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c8(p["c8"]),
            abs_c7(p["Abs{c7}"]),
            arg_c7(p["Arg{c7}"]),
            abs_c7prime(p["Abs{c7'}"]),
            arg_c7prime(p["Arg{c7'}"]),
            a_1_perp(p["B->K^*::a_1_perp"]),
            a_2_perp(p["B->K^*::a_2_perp"]),
            f_B(p["f_B"]),
            f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"]),
            lambda_B_p(p["lambda_B_p"]),
            uncertainty_xi_perp(p["formfactors::xi_perp_uncertainty"]),
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"]),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            form_factors(FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "BZ2004"), p))
        {
        }

        inline double m_b_PS() const { return model->m_b_ps(std::sqrt(mu * 0.5)); }

        // alias c7('),c9(') and c10(')
        inline complex<double> c7() const { return std::polar(abs_c7(), (cp_conjugate ? -1.0 : +1.0) * arg_c7()); }
        inline complex<double> c7prime() const { return std::polar(abs_c7prime(), (cp_conjugate ? -1.0 : +1.0) * arg_c7prime()); }

        // cf. [BFS2001], Eq. (10), p. 4
        complex<double> Y0(const double & s) const
        {
            static const double lambda_QCD = 0.5; // GeV

            double Y_c = 4.0 / 3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5();
            double Y_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * c4() + 76.0 * c5() + 64.0 / 3.0 * c6());
            double Y_0 = -0.5 * (c3() + 4.0 / 3.0 * c4() + 16.0 * c5() + 64 / 3.0 * c6());
            double Y = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * CharmLoops::h(mu, s, model->m_c_pole())
                + Y_b * CharmLoops::h(mu, s, m_b_PS())
                + Y_0 * CharmLoops::h(mu, s, lambda_QCD)
                + Y;
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);
            double result = uncertainty_xi_perp * factor * form_factors->v(s);

            return result;
        }

        /* Effective wilson coefficient */
        // cf. [BFS2001], below Eq. (9), p. 4
        inline complex<double> c7eff() const
        {
            return c7() - 1.0/3.0 * c3() - 4.0/9.0 * c4() - 20.0/3.0 * c5() - 80.0/9.0 * c6();
        }

        inline double energy(const double & s) const
        {
            return (m_B * m_B + m_Kstar * m_Kstar - s) / (2.0 * m_B);
        }

        // cf. [BFS2001], Eq. (36), p. 9
        double L(const double & s) const
        {
            if (std::abs(s) < 1e-4)
                return -1.0;

            double m_b = m_b_PS();
            double m_b2 = m_b * m_b;

            return -1.0 * (m_b2 - s) / s * std::log(1.0 - s / m_b2);
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

        // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
        complex<double> C0_perp(const double & h, const double & s) const
        {
            return (c7eff() + h * c7prime()) + s / (2.0 * m_b_PS() * m_B) * Y0(s);
        }

        // cf. [BFS2001], Eqs. (34), (37), p. 9 in the limit q2 -> 0
        complex<double> C1_perp(const double & h, const double & s) const
        {
            static const double L = -1.0;

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();
            double m_c = model->m_c_pole();

            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = mu * model->m_b_pole() / m_b;

            // cf. [BFS2004], Eq. (44), p. 24
            // [Christoph] Use c7 instead of c7eff
            complex<double> C_perp_f = (c7() + h * c7prime()) * (8.0 * std::log(m_b / mu) - L - 4.0 * (1.0 - std::sqrt(mu * 0.5) / m_b));

            // cf. [BFS2001], Eq. (37), p. 9
            // [Christoph] Use c8 instead of c8eff
            complex<double> C_perp_nf = (-1.0 / QCD::casimir_f) * (
                    (c2 - c1 / 6.0) * memoise(CharmLoops::F27_massive, mu_pole, s, m_b, m_c) + c8() * CharmLoops::F87_massless(mu_pole, s, m_b));

            return C_perp_f + C_perp_nf;
        }

        // cf. [BFS2001], Eqs. (16), (21), (25), pp. 5-7
        complex<double> T1_perp_p(const double & h, const double & s, const double & u) const
        {
            static const double e_d = (-1.0/3.0);
            static const double e_u = (+2.0/3.0);

            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            double m_b = m_b_PS();
            double m_c = model->m_c_pole();
            double s_hat = s / m_B / m_B;

            // cf. [BFS2001], Eq. (20)
            // [Christoph] Use c7 instead of c7eff
            complex<double> Tf_perp_p = (c7() + h * c7prime()) * (2.0 * m_B / (1.0 - u) / energy(s));

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
        complex<double> T_perp(const double & h, const double & s, const double & u) const
        {
            double a = model->alpha_s(sqrt(mu * 0.5)) * QCD::casimir_f / 4.0 / M_PI;

            return phi_K(u, a_1_perp, a_2_perp) * a * T1_perp_p(h, s, u);

            // TODO: Hard scattering + Weak annihilation from [BFS2004], Eqs. (51), (52)
        }

        // cf. [BFS2001], Eq. (15) with a = perp
        complex<double> calT_perp(const double & h, const double & s) const
        {
            return xi_perp(s) * (C0_perp(h, s) + model->alpha_s(mu) * QCD::casimir_f / 4.0 / M_PI * C1_perp(h, s))
                + (power_of<2>(M_PI) / 3.0) * (f_B * f_Kstar_perp / m_B)
                * integrate(std::function<complex<double> (const double &)>(
                                std::bind(&Implementation<BToKstarGamma>::T_perp, this, h, s, std::placeholders::_1)),
                        64, 0.001, 0.999);
        }

        struct Amplitudes
        {
            complex<double> left;
            complex<double> right;
        };

        Amplitudes amplitudes() const
        {
            complex<double> calT_perp_plus = calT_perp(+1.0, 0.0);
            complex<double> calT_perp_minus = calT_perp(-1.0, 0.0);

            complex<double> a_left  = complex<double>(0.0, +0.5) * (calT_perp_plus + calT_perp_minus);
            complex<double> a_right = complex<double>(0.0, -0.5) * (calT_perp_plus - calT_perp_minus);

            return Amplitudes{ a_left, a_right };
        }

        double branching_ratio()
        {
            // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
            static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5

            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));
            Amplitudes a = amplitudes();

            return alpha_e * power_of<2>(g_fermi * model->m_b_msbar(mu())) * power_of<3>(m_B()) / (32.0 * power_of<4>(M_PI)) *
                power_of<3>(1.0 - power_of<2>(m_Kstar / m_B)) * lambda_t * lambda_t * (std::norm(a.left) + std::norm(a.right)) / Gamma;
        }

        double s_kstar_gamma()
        {
            Save<bool> save(this->cp_conjugate, false);

            Amplitudes abar = amplitudes();
            cp_conjugate = true;
            Amplitudes a = amplitudes();

            double phi_d = arg(power_of<2>(conj(model->ckm_td()) * model->ckm_tb()));
            complex<double> q_over_p = std::polar(1.0, -phi_d);

            double numerator = -2.0 * imag(q_over_p * (conj(a.left) * abar.right + conj(a.right) * abar.left));
            double denominator = std::norm(a.left) + std::norm(a.right) + std::norm(abar.left) + std::norm(abar.right);

            return numerator / denominator;
        }
    };

    BToKstarGamma::BToKstarGamma(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKstarGamma>(new Implementation<BToKstarGamma>(parameters, options))
    {
    }

    BToKstarGamma::~BToKstarGamma()
    {
    }

    double
    BToKstarGamma::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BToKstarGamma::s_kstar_gamma() const
    {
        return _imp->s_kstar_gamma();
    }
}
