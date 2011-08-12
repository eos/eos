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

#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/form-factors.hh>
#include <eos/rare-b-decays/hard-scattering.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/save.hh>

#include <cmath>

#include <gsl/gsl_sf.h>

namespace eos
{
    struct ShortDistanceQCDF
    {
        struct ParameterSet
        {
            double                    m_b_PS, m_b_pole, m_c;
            double                    m_B, m_K;

            double                    mu;
            double                    mu_f;

            double                    alpha_s_mu, alpha_s_sqrt05mu;

            double                    f_B, f_K;

            WilsonCoefficients<BToS>  wc;

            double                    e_q;

            double                    a_1, a_2;

            double                    lambda_B_p;

            ParameterSet(const double & m_b_PS_, const double & m_b_pole_, const double & m_c_,
                    const double & m_B_, const double & m_K_,
                    const double & mu_, const double & mu_f_,
                    const double & alpha_s_mu_, const double & alpha_s_sqrt05mu_,
                    const double & f_B_, const double & f_K_,
                    const WilsonCoefficients<BToS> & wc_,
                    const double & e_q_,
                    const double & a_1_, const double & a_2_,
                    const double & lambda_B_p_) :
                 m_b_PS(m_b_PS_), m_b_pole(m_b_pole_), m_c(m_c_),
                 m_B(m_B_), m_K(m_K_),
                 mu(mu_), mu_f(mu_f_),
                 alpha_s_mu(alpha_s_mu_), alpha_s_sqrt05mu(alpha_s_sqrt05mu_),
                 f_B(f_B_), f_K(f_K_),
                 wc(wc_),
                 e_q(e_q_),
                 a_1(a_1_), a_2(a_2_),
                 lambda_B_p(lambda_B_p_)
            { }
        };

        // for s = 0
        static double energy(const ParameterSet & p)
        {
            return (p.m_B * p.m_B + p.m_K * p.m_K) / (2.0 * p.m_B);
        }

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (9), p. 4
        static complex<double> c7eff(const ParameterSet & p)
        {
            return p.wc.c7() - 1.0/3.0 * p.wc.c3() - 4.0/9.0 * p.wc.c4() - 20.0/3.0 * p.wc.c5() - 80.0/9.0 * p.wc.c6();
        }

        // cf. [BFS2001], below Eq. (26), p. 8
        static complex<double> c8eff(const ParameterSet & p)
        {
             return p.wc.c8() + p.wc.c3() - 1.0/6.0 * p.wc.c4() + 20.0 * p.wc.c5() - 10.0/3.0 * p.wc.c6();
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

        // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
        static complex<double> C0_perp(const double & h, const ParameterSet & p)
        {
            const double p_l = (1.0 + h) / 2.0, p_r = (1.0 - h) / 2.0;

            return (p_l * c7eff(p) + p_r * p.wc.c7prime());
        }

        // cf. [BFS2001], Eqs. (34), (37), p. 9
        static complex<double> C1f_perp(const double & h, const ParameterSet & p)
        {
            const double p_l = (1.0 + h) / 2.0, p_r = (1.0 - h) / 2.0;

            // TODO
            static const double L = -1.0;

            // cf. [BFS2004], Eq. (44), p. 24
            // [Christoph] Use c7 instead of c7eff
            return (p_l * p.wc.c7() + p_r * p.wc.c7prime()) * (8.0 * std::log(p.m_b_PS / p.mu) - L - 4.0 * (1.0 - p.mu_f / p.m_b_PS));
        }

        // cf. [BFS2001], Eqs. (34), (37), p. 9
        static complex<double> C1nf_perp(const double & h, const ParameterSet & p)
        {
            // Here m_b_PS is used instead of m_b_pole, cf. [BFS2001], comment below Eq. (36), p. 9
            // Two-Loop Function are calculated for the pole mass! Use mu_pole instead
            double mu_pole = p.mu * p.m_b_pole / p.m_b_PS;

            const double p_l = (1.0 + h) / 2.0;

            // cf. [BFS2001], Eq. (37), p. 9
            // [Christoph] Use c8 instead of c8eff
            return p_l * (-1.0 / QCD::casimir_f) * (
                    (p.wc.c2() - p.wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu_pole, 0.0, p.m_b_PS, p.m_c)
                    + p.wc.c8() * CharmLoops::F87_massless(mu_pole, 0.0, p.m_b_PS));
        }

        // cf. [BFS2001], Eqs. (16), (21), (25), pp. 5-7
        static complex<double> T1f_perp_p(const double & h, const double & u, const ParameterSet & p)
        {
            const double p_l = (1.0 + h) / 2.0, p_r = (1.0 - h) / 2.0;

            // cf. [BFS2001], Eq. (20)
            // [Christoph] Use c7 instead of c7eff
            return (p_l * p.wc.c7() + p_r * p.wc.c7prime()) * (2.0 * p.m_B / (1.0 - u) / energy(p));
        }

        static complex<double> T1nf_perp_p(const double & h, const double & u, const ParameterSet & p)
        {
            static const double e_d = -1.0/3.0;
            static const double e_u = +2.0/3.0;

            const double p_l = (1.0 + h) / 2.0;

            // cf. [BFS2001], Eq. (23)
            // [Christoph] Use c8 instead of c8eff
            return p_l * (-4.0 * e_d * p.wc.c8() / u
                + p.m_B / (2.0 * p.m_b_PS) * (
                        e_u * (-p.wc.c1() / 6.0 + p.wc.c2() + 6.0 * p.wc.c6()) * t_perp_0(u, p, p.m_c)
                        + e_d * (p.wc.c3() - p.wc.c4() / 6.0 + 16.0 * p.wc.c5() + 10.0/3.0 * p.wc.c6() - (4.0 * p.m_b_PS / p.m_B) * (p.wc.c3() - p.wc.c4()/6.0 + 4.0 * p.wc.c5() - 2.0/3.0 * p.wc.c6()))
                        + e_d * (p.wc.c3() - p.wc.c4() / 6.0 + 16.0 * p.wc.c5() - 8.0/3.0 * p.wc.c6()) * t_perp_0(u, p, 0.0)));
        }

        // cf. [BFS2001], Eq. (16) times phi_K^*_perp
        static complex<double> T_perp_sum(const double & h, const double & u, const ParameterSet & p)
        {
             double a = p.alpha_s_sqrt05mu * QCD::casimir_f / 4.0 / M_PI;

             complex<double> result = 1.0 / p.lambda_B_p * a * (T1f_perp_p(h, u, p) + T1nf_perp_p(h, u, p));

             return phi_K(u, p) * result;

             // TODO: Hard scattering + Weak annihilation from [BFS2004], Eqs. (51), (52)
        }

        // cf. [BFS2001], Eq. (15) with a = perp
        static complex<double> calT_perp(const double & h, const ParameterSet & p, const double & xi_perp)
        {
            std::complex<double> result = xi_perp * (C0_perp(h, p) + p.alpha_s_mu * QCD::casimir_f / 4.0 / M_PI * (C1f_perp(h, p) + C1nf_perp(h, p)));

            result += power_of<2>(M_PI) / 3.0 * (p.f_B * p.f_K) / p.m_B *
                    integrate(std::function<complex<double> (const double &)>(
                                std::bind(&ShortDistanceQCDF::T_perp_sum, h, std::placeholders::_1, p)),
                        64, 0.001, 0.999);

            return result;
        }
    };

    template <>
    struct Implementation<BToKstarGamma>
    {
        std::shared_ptr<Model> model;

        UsedParameter a_1_perp;

        UsedParameter a_2_perp;

        UsedParameter f_B;

        UsedParameter f_Kstar_perp;

        UsedParameter lambda_B_p;

        UsedParameter m_B;

        UsedParameter m_Kstar;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            a_1_perp(p["B->K^*::a_1_perp"], u),
            a_2_perp(p["B->K^*::a_2_perp"], u),
            f_B(p["decay-constant::B_" + o.get("q", "d")], u),
            f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"], u),
            lambda_B_p(p["lambda_B_p"], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            m_Kstar(p["mass::K^*0"], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            form_factors(FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "KMPW2010"), p))
        {
            u.uses(*model);
            u.uses(*form_factors);
        }

        inline double m_b_PS() const { return model->m_b_ps(std::sqrt(mu * 0.5)); }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);
            double result = factor * form_factors->v(s);

            return result;
        }

        struct Amplitudes
        {
            complex<double> left;
            complex<double> right;
        };

        Amplitudes amplitudes() const
        {
            ShortDistanceQCDF::ParameterSet p(m_b_PS(), model->m_b_pole(), model->m_c_pole(),
                    m_B(), m_Kstar(),
                    mu(), std::sqrt(0.5 * mu()),
                    model->alpha_s(mu()), model->alpha_s(sqrt(mu() * 0.5)),
                    f_B(), f_Kstar_perp(),
                    model->wilson_coefficients_b_to_s(cp_conjugate),
                    -1.0 / 3.0,
                    a_1_perp(), a_2_perp(),
                    lambda_B_p());

            complex<double> a_left = complex<double>(0.0, +1.0) * ShortDistanceQCDF::calT_perp(+1.0, p, xi_perp(0.0));
            complex<double> a_right = complex<double>(0.0, -1.0) * ShortDistanceQCDF::calT_perp(-1.0, p, xi_perp(0.0));

            return Amplitudes{ a_left, a_right };
        }

        double branching_ratio()
        {
            // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
            static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));
            Amplitudes a = amplitudes();

            return alpha_e() * power_of<2>(g_fermi() * model->m_b_msbar(mu())) * power_of<3>(m_B()) / (32.0 * power_of<4>(M_PI)) *
                power_of<3>(1.0 - power_of<2>(m_Kstar / m_B)) * lambda_t * lambda_t * (std::norm(a.left) + std::norm(a.right)) / Gamma;
        }

        double branching_ratio_cp_averaged()
        {
            // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
            static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

            Save<bool> save(this->cp_conjugate, false);
            Amplitudes a = amplitudes();
            cp_conjugate = true;
            Amplitudes abar = amplitudes();

            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));
            return alpha_e() * power_of<2>(g_fermi() * model->m_b_msbar(mu())) * power_of<3>(m_B()) / (32.0 * power_of<4>(M_PI)) *
                power_of<3>(1.0 - power_of<2>(m_Kstar / m_B)) * lambda_t * lambda_t * (std::norm(a.left) + std::norm(a.right) + std::norm(abar.left) + std::norm(abar.right)) / 2.0 / Gamma;
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

        double c_kstar_gamma()
        {
            Save<bool> save(this->cp_conjugate, false);

            Amplitudes abar = amplitudes();
            cp_conjugate = true;
            Amplitudes a = amplitudes();

            double numerator = std::norm(a.left) + std::norm(a.right) - std::norm(abar.left) - std::norm(abar.right);
            double denominator = std::norm(a.left) + std::norm(a.right) + std::norm(abar.left) + std::norm(abar.right);

            return numerator / denominator;
        }
    };

    BToKstarGamma::BToKstarGamma(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKstarGamma>(new Implementation<BToKstarGamma>(parameters, options, *this))
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
    BToKstarGamma::branching_ratio_cp_averaged() const
    {
        return _imp->branching_ratio_cp_averaged();
    }

    double
    BToKstarGamma::s_kstar_gamma() const
    {
        return _imp->s_kstar_gamma();
    }

    double
    BToKstarGamma::c_kstar_gamma() const
    {
        return _imp->c_kstar_gamma();
    }
}
