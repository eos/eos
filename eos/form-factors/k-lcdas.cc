/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2024 Danny van Dyk
 * Copyright (c) 2022 Carolina Bolognani
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

#include <eos/form-factors/k-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<AntiKaonLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1K_0;
        UsedParameter a2K_0;

        // twist 3 parameters
        UsedParameter f3K_0;
        UsedParameter lambda3K_0;
        UsedParameter omega3K_0;

        // twist 4 parameters
        UsedParameter delta4K_0;
        UsedParameter kappa4K_0;
        UsedParameter omega4K_0;

        // mass and decay constant of the pion
        UsedParameter m_K;
        UsedParameter f_K;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1K_0(p["K::a1@1GeV"], u),
            a2K_0(p["K::a2@1GeV"], u),
            f3K_0(p["K::f3@1GeV"], u),
            lambda3K_0(p["K::lambda3@1GeV"], u),
            omega3K_0(p["K::omega3@1GeV"], u),
            delta4K_0(p["K::delta4@1GeV"], u),
            kappa4K_0(p["K::kappa4@1GeV"], u),
            omega4K_0(p["K::omega4@1GeV"], u),
            m_K(p["mass::K_u"], u),
            f_K(p["decay-constant::K_u"], u),
            _mu_c(p["QCD::mu_c"], u),
            _mu_b(p["QCD::mu_b"], u),
            _mu_t(p["QCD::mu_t"], u)
        {
        }

        inline double c_rge(const double & _mu) const
        {
            /*
             * RGE coefficient, basically
             *
             *     (alpha_s/alpha_s_0)^(1/beta_0),
             *
             * with matching between the individual n-flavor QCDs.
             */

            double mu = _mu, alpha_s_mu = model->alpha_s(mu);
            double mu_0 = 1.0, alpha_s_0 = model->alpha_s(mu_0);

            if (mu < _mu_c)
                return std::pow(alpha_s_mu / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            double alpha_s_c = model->alpha_s(_mu_c);
            double result = std::pow(alpha_s_c / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            if (mu < _mu_b)
                return result * std::pow(alpha_s_mu / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            double alpha_s_b = model->alpha_s(_mu_b);
            result *= std::pow(alpha_s_b / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            if (mu < _mu_t)
                return result * std::pow(alpha_s_mu / alpha_s_b, 1.0 / QCD::beta_function_nf_5[0]);

            throw InternalError("Implementation<AntiKaonLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a1K(const double & mu) const
        {
            return a1K_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a2K(const double & mu) const
        {
            return a2K_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double muK(const double & mu) const
        {
            return m_K * m_K / (model->m_s_msbar(mu) + model->m_ud_msbar(mu) / 2.0);
        }

        double f3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_s_msbar(mu_0);
            const double m_q_0 = model->m_ud_msbar(mu_0) / 2.0;

            return f3K_0 * std::pow(c_rge, 55.0 / 9.0)
                + 2.0 / 19.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 55.0 / 9.0)) * f_K * (m_s_0 + m_q_0)
                + 6.0 / 65.0 * (std::pow(c_rge, 55.0 / 9.0) - std::pow(c_rge, 68.0 / 9.0)) * f_K * (m_s_0 - m_q_0) * a1K_0;
        }

        double omega3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_s_msbar(mu_0);
            const double m_q_0 = model->m_ud_msbar(mu_0) / 2.0;

            return (f3K_0 * omega3K_0 * std::pow(c_rge, 104.0 / 9.0)
                + 1.0 / 170.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 104.0 / 9.0)) * f_K * (m_s_0 + m_q_0)
                + 1.0 /  10.0 * (std::pow(c_rge, 68.0 / 9.0) - std::pow(c_rge, 104.0 / 9.0)) * f_K * (m_s_0 - m_q_0) * a1K_0
                + 2.0 /  15.0 * (std::pow(c_rge, 86.0 / 9.0) - std::pow(c_rge, 104.0 / 9.0)) * f_K * (m_s_0 + m_q_0) * a2K_0) / f3K(mu);
        }

        double lambda3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_s_msbar(mu_0);
            const double m_q_0 = model->m_ud_msbar(mu_0) / 2.0;

            return (f3K_0 * lambda3K_0 * std::pow(c_rge, 139.0 / 18.0)
                - 14.0 / 67.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 139.0 / 18.0)) * f_K * (m_s_0 - m_q_0)
                + 14.0 /  5.0 * (std::pow(c_rge, 68.0 / 9.0) - std::pow(c_rge, 139.0 / 18.0)) * f_K * (m_s_0 + m_q_0) * a1K_0
                - 4.0  / 11.0 * (std::pow(c_rge, 86.0 / 9.0) - std::pow(c_rge, 139.0 / 18.0)) * f_K * (m_s_0 - m_q_0) * a2K_0) / f3K(mu);
        }

        inline double eta3K(const double & mu) const
        {
            return f3K(mu) / (f_K() * muK(mu));
        }

        double delta4K(const double & mu) const
        {
            const double c_rge  = this->c_rge(mu);

            return delta4K_0 * std::pow(c_rge, 32.0 / 9.0) + 1.0 / 8.0 * m_K * m_K * (1.0 - std::pow(c_rge, 32.0 / 9.0));
        }

        double kappa4K(const double & mu) const
        {
            const double c_rge  = this->c_rge(mu);
            const double mu_0   = 1.0;
            const double m_s_0 = model->m_s_msbar(mu_0);
            const double m_q_0 = model->m_ud_msbar(mu_0) / 2.0;

            return kappa4K_0
            - 9.0 / 40.0 * a1K_0 * (std::pow(c_rge, 32.0 / 9.0) - 1.0)
            + (m_s_0 * m_s_0 - m_q_0 * m_q_0) / (2.0 * m_K * m_K) * (std::pow(c_rge, 8.0) - 1.0);
        }

        double omega4K(const double & mu) const
        {
            const double c_rge  = this->c_rge(mu);

            return 1.0 / delta4K(mu) * omega4K_0 * delta4K_0 * std::pow(c_rge, 10.0);
        }

    };

    AntiKaonLCDAs::AntiKaonLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AntiKaonLCDAs>(new Implementation<AntiKaonLCDAs>(p, o, *this))
    {
    }

    AntiKaonLCDAs::~AntiKaonLCDAs()
    {
    }

    PseudoscalarLCDAs *
    AntiKaonLCDAs::make(const Parameters & p, const Options & o)
    {
        return new AntiKaonLCDAs(p, o);
    }

    double
    AntiKaonLCDAs::a1(const double & mu) const
    {
        return _imp->a1K(mu);
    }

    double
    AntiKaonLCDAs::a2(const double & mu) const
    {
        return _imp->a2K(mu);
    }

    double
    AntiKaonLCDAs::mu3(const double & mu) const
    {
        return _imp->muK(mu);
    }

    double
    AntiKaonLCDAs::f3(const double & mu) const
    {
        return _imp->f3K(mu);
    }

    double
    AntiKaonLCDAs::eta3(const double & mu) const
    {
        return _imp->eta3K(mu);
    }

    double
    AntiKaonLCDAs::lambda3(const double & mu) const
    {
        return _imp->lambda3K(mu);
    }

    double
    AntiKaonLCDAs::omega3(const double & mu) const
    {
        return _imp->omega3K(mu);
    }

    double
    AntiKaonLCDAs::delta4(const double & mu) const
    {
        return _imp->delta4K(mu);
    }

    double
    AntiKaonLCDAs::kappa4(const double & mu) const
    {
        return _imp->kappa4K(mu);
    }

    double
    AntiKaonLCDAs::omega4(const double & mu) const
    {
        return _imp->omega4K(mu);
    }

    double
    AntiKaonLCDAs::phi(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1_3o2(1.0, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2_3o2(2.0, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1_3o2.evaluate(x);
        const double c2 = gp_2_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1K(mu) * c1 + _imp->a2K(mu) * c2);
    }

    double
    AntiKaonLCDAs::phi3p(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(1/2)
        static const GegenbauerPolynomial gp_1_1o2(1.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_2_1o2(2.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_3_1o2(3.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_4_1o2(4.0, 1.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1_1o2.evaluate(x);
        const double c2 = gp_2_1o2.evaluate(x);
        const double c3 = gp_3_1o2.evaluate(x);
        const double c4 = gp_4_1o2.evaluate(x);

        return 1.0 + 3.0 * rhopK * (1.0 + 6.0 * a2K) - 9.0 * rhomK * a1K
            + c1 * (27.0 / 2.0 * rhopK * a1K - rhomK * (3.0 / 2.0 + 27.0 * a2K))
            + c2 * (30.0 * eta3K + 15.0 * rhopK * a2K - 3.0 * rhomK * a1K)
            + c3 * (10.0 * eta3K * lambda3K - 9.0 / 2.0 * rhomK * a2K)
            + c4 * (-3.0 * eta3K * omega3K)
            + 3.0 / 2.0 * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 3.0 / 2.0 * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u);
    }

    double
    AntiKaonLCDAs::phi3s(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1_3o2(1.0, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2_3o2(2.0, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3_3o2(3.0, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1_3o2.evaluate(x);
        const double c2 = gp_2_3o2.evaluate(x);
        const double c3 = gp_3_3o2.evaluate(x);

        const double ubar = 1.0 - u;

        return 6.0 * u * ubar * (
            1.0 + 3.0 / 2.0 * rhopK + 15.0 * rhopK * a2K - 15.0 / 2.0 * rhomK * a1K
            + c1 * (3.0 * rhopK * a1K - 15.0 / 2.0 * rhomK * a2K)
            + c2 * (5.0 * eta3K - 1.0 / 2.0 * eta3K * omega3K + 3.0 / 2.0 * rhopK * a2K)
            + c3 * (eta3K * lambda3K))
            + 9.0 * u * ubar * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 9.0 * u * ubar * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u)
        ;
    }

    double
    AntiKaonLCDAs::phi3s_d1(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        const double ubar = 1.0 - u, x = 2.0 * u - 1.0;
        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u;

        return -3.0 * (60.0 * eta3K * (-1.0 + 12.0 * u - 30.0 * u2 + 20.0 * u3)
            -3.0 * rhomK * (1.0 + x * (-1.0 * std::log(u) +  std::log(ubar) + a1K * (8.0 + 3.0 * (std::log(u) + std::log(ubar))))
            + 3.0 * a2K  * (7.0 - 30.0 * u + 30.0 * u2 + 2.0 * x * (std::log(ubar) - std::log(u))))
            + x * (2.0 + 3.0 * rhopK * (2.0 + std::log(u) + std::log(ubar) + 3.0 * a1K * (-3.0 + 6.0 * u - std::log(u) + std::log(ubar))
            + a2K * (22.0 - 60.0 * u + 60.0 * u2 + 6.0 * (std::log(u) + std::log(ubar))))))
            - 6.0 * eta3K * (10.0 * lambda3K * (1.0 - 20.0 * u + 90.0 * u2 - 140.0 * u3 + 70.0 * u4) + 3.0 * omega3K * (1.0 - 12.0 * u + 30.0 * u2 - 20.0 * u3))
        ;
    }

    double
    AntiKaonLCDAs::phi4(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double omega4K  = _imp->omega4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;
        const double phi2K    = -7.0 / 20.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, ubar3 = ubar2 * ubar, lnubar = std::log(ubar);
        const double x = 2.0 * u - 1.0;

        // Twist 4 contributions
        const double phi4T4 = 200.0 / 3.0 * delta4K * u2 * ubar2 + 20.0 * u2 * ubar2 * x * (4.0 * theta1K - 5.0 * theta2K)
            + 21.0 * delta4K * omega4K * (u * ubar * (2.0 + 13.0 * u * ubar) + (2.0 * u3 * (6.0 * u2 - 15.0 * u + 10.0) * lnu)
            + (2.0 * ubar3 * (6.0 * ubar2 - 15.0 * ubar + 10.0) * lnubar))
            + 40.0 * phi2K * (u * ubar * x * (2.0 - 3.0 * u * ubar) - (2.0 * u3 * (u - 2.0) * lnu)
            + (2.0 * ubar3 * (ubar - 2.0) * lnubar));
        const double phi4WW = 16.0 / 3.0 * m_K * m_K * kappa4K * (u * ubar * x * (1.0 - 2.0 * u * ubar)
            + (5.0 * (u - 2.0) * u3 * lnu) - (5.0 * (ubar - 2.0) * ubar3 * lnubar))
            + 4.0 * f3K / f_K * (m_s + m_ud) * u * ubar * (30.0 * (1.0 - x * (m_s - m_ud) / (m_s + m_ud))
            + 10.0 * lambda3K * (x * (1.0 - u * ubar) - (m_s - m_ud) / (m_s + m_ud) * (1.0 - 5.0 * u * ubar))
            - omega3K * (3.0 - 21.0 * u * ubar + 28.0 * u2 * ubar2 + 3.0 * x * (m_s - m_ud) / (m_s + m_ud) * (1.0 - 7.0 * u * ubar)))
            - 36.0 / 5.0 * m_K * m_K * a2K * (1.0 / 4.0 * u * ubar * (4.0 - 9.0 * u * ubar + 110.0 * u2 * ubar2)
            + (u3 * (10.0 - 15.0 * u + 6.0 * u2) * lnu) + (ubar3 * (10.0 - 15.0 * ubar + 6.0 * ubar2) * lnubar))
            + 4.0 * m_K * m_K * u * ubar * (1.0 + 3.0 * u * ubar) * (1.0 + 9.0 / 5.0 * a1K * x);

        return phi4T4 + phi4WW;
    }

    double
    AntiKaonLCDAs::phi4_d1(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double omega4K  = _imp->omega4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;
        const double phi2K    = -7.0 / 20.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, u5 = u4 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, lnubar = std::log(ubar);
        const double x = 2.0 * u - 1.0;

        // Twist 4 derivatives contributions
        const double phi4T4_d1 = 20.0 * ubar2 * lnubar * (8.0 * (1.0 + 2.0 * u) * phi2K - 63.0 * u2 * omega4K * delta4K)
            + 20.0 * u2 * lnu * (-8.0 * (2.0 * u - 3.0) * phi2K + 63.0 * ubar2 * omega4K * delta4K)
            - 20.0 / 3.0 * u * (-20.0 * (1.0 - 3.0 * u + 2.0 * u2) * delta4K + 12.0 * (-8.0 + 23.0 * u - 30.0 * u2 + 15.0 * u3) * phi2K
            + 3.0 * ubar * ((1.0 - 5.0 * u + 5.0 * u2) * (8.0 * theta1K - 10.0 * theta2K) + 21.0 * x * omega4K * delta4K));
        const double phi4WW_d1 = 36.0 / 5.0 * m_K * m_K * a1K * (-1.0 + 30.0 * u2 - 60.0 * u3 + 30.0 * u4)
            - 54.0 * m_K * m_K * a2K * u * ubar * (-1.0 + 13.0 * u - 33.0 * u2 + 22.0 * u3)
            - 16.0 / 3.0 * m_K * m_K * kappa4K * (6.0 - 15.0 * u + 35.0 * u2 - 40.0 * u3 + 20.0 * u4)
            + 8.0 / 3.0 * m_K * m_K * (u2 * lnu * (-81.0 * ubar2 * a2K + 20.0 * (-3.0 + 2.0 * u) * kappa4K)
            + ubar2 * lnubar * (81.0 * u2 * a2K - 20.0 * (2.0 * u + 1.0) * kappa4K))
            + 4.0 / f_K * (f_K * m_K * m_K * (1.0 + 4.0 * u - 18.0 * u2 + 12.0 * u3)
            + f3K * (60.0 * (m_s * (1.0 - 4.0 * u + 3.0 * u2) + m_ud * u * (2.0 - 3.0 * u))
            - 20.0 * lambda3K * (m_s * (1.0 - 10.0 * u + 24.0 * u2 - 20.0 * u3 + 5.0 * u4) + m_ud * u * (2.0 - 6.0 * u + 5.0 * u3))
            + omega3K * m_s * (-12.0 * u + 60.0 * u2 - 210.0 * u4 + 168.0 * u5)
            + omega3K * m_ud * (-6.0 + 108.0 * u - 480.0 * u2 + 840.0 * u3 - 630.0 * u4 + 168.0 * u5)));

        return phi4T4_d1 + phi4WW_d1;
    }

    double
    AntiKaonLCDAs::phi4_d2(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double omega4K  = _imp->omega4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;
        const double phi2K    = -7.0 / 20.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        // Twist 4 derivatives contributions
        const double phi4T4_d2 = 400.0 / 3.0 * (1.0 - 6.0 * u + 6.0 * u2) * delta4K
            - 20.0 * (24.0 * phi2K * (-1.0 + 7.0 * u - 15.0 * u2 + 10.0 * u3)
            + (-1.0 + 12.0 * u - 30.0 * u2 + 20.0 * u3) * (- 8.0 * theta1K + 10.0 * theta2K)
            - 21.0 * omega4K * delta4K * (1.0 - 3.0 * u + 3.0 * u2))
            + 120.0 * u * (8.0 * ubar * phi2K + 21.0 * (1.0 - 3.0 * u + 2.0 * u2) * omega4K * delta4K) * (lnu - lnubar);
        const double phi4WW_d2 = m_K * m_K * (432.0 * u * (1.0 - 3.0 * u + 2.0 * u2) * a1K
            + 54.0 * (1.0 - 32.0 * u + 142.0 * u2 - 220.0 * u3 + 110.0 * u4) * a2K
            - 80.0 / 3.0 * (-5.0 + 18.0 * u - 24.0 * u2 + 16.0 * u3) * kappa4K
            + 16.0 * u * (27.0 * (1.0 - 3.0 * u + 2.0 * u2) * a2K + 20.0 * ubar * kappa4K) * (-lnu + lnubar))
            + 16.0 / f_K * (f_K * m_K * m_K * (1.0 - 9.0 * u + 9.0 * u2)
            + f3K * (30.0 * (m_s * (-2.0 + 3.0 * u) + m_ud * (1.0 - 3.0 * u))
            - 10.0 * lambda3K * (m_s * (-5.0 + 24.0 * u - 30.0 * u2 + 10.0 * u3) + m_ud * (1.0 - 6.0 * u + 10.0 * u3))
            + omega3K * (m_s * (-3.0 + 30.0 * u - 210.0 * u3 + 210.0 * u4)
            + m_ud * (27.0 - 240.0 * u + 630.0 * u2 - 630.0 * u3 + 210.0 * u4))));

        return phi4T4_d2 + phi4WW_d2;
    }

    double
    AntiKaonLCDAs::psi4(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;

        const double lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        // Gegenbauer polynomials C_n^(1/2)
        static const GegenbauerPolynomial gp_1_1o2(1.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_2_1o2(2.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_3_1o2(3.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_4_1o2(4.0, 1.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c0 = 1.0;
        const double c1 = gp_1_1o2.evaluate(x);
        const double c2 = gp_2_1o2.evaluate(x);
        const double c3 = gp_3_1o2.evaluate(x);
        const double c4 = gp_4_1o2.evaluate(x);

        // Twist 4 contributions
        const double psi4T4 = 20.0 / 3.0 * delta4K * c2 + 5.0 * (5.0 * theta1K - theta2K) * c3;
        const double psi4WW = c0 * m_K * m_K * (1.0 + 6.0 * rhopK * (1.0 + 6.0 * a2K) - 18.0 * rhomK * a1K)
            + c1 * m_K * m_K * (-12.0 * kappa4K - 9.0 / 5.0 * a1K + 27.0 * rhopK * a1K - 3.0 * rhomK * (1.0 + 18.0 * a2K))
            + c2 * (m_K * m_K * (1.0 + 18.0 / 7.0 * a2K + 30.0 * rhopK * a2K - 6.0 * rhomK * a1K) + 60.0 * f3K / f_K * (m_s + m_ud))
            + c3 * (m_K * m_K * (9.0 / 5.0 * a1K + 16.0 / 3.0 * kappa4K - 9.0 * rhomK * a2K) + 20.0 * f3K * lambda3K / f_K * (m_s + m_ud))
            + c4 * (-9.0 / 28.0 * m_K * m_K * a2K - 6.0 * f3K * omega3K / f_K * (m_s + m_ud))
            + 6.0 * m_ud * (m_s + m_ud) * (1.0 + 3.0 * a1K + 6.0 * a2K) * lnubar
            + 6.0 * m_s * (m_s + m_ud) * (1.0 - 3.0 * a1K + 6.0 * a2K) * lnu;

        return psi4T4 + psi4WW;
    }

    double
    AntiKaonLCDAs::psi4_i(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_s_msbar(mu);
        const double m_ud = _imp->model->m_ud_msbar(mu) / 2.0;

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        // Twist 4 contributions
        const double psi4T4_i = -5.0 / 3.0 * u * ubar * (delta4K * (8.0 * u - 4.0)
            + 3.0 * (1.0 - 5.0 * u + 5.0 * u2) * (5.0 * theta1K - theta2K));
        const double psi4WW_i = 20.0 / 3.0 * m_K * m_K * kappa4K * u * (1.0 + 3.0 * u - 8.0 * u2 + 4.0 * u3)
            - 6.0 * m_s * (m_s + m_ud) * u * (-1.0 + 3.0 * a1K - 6.0 * a2K) * lnu
            + 6.0 * m_ud * (m_s + m_ud) * ubar * (-1.0 - 3.0 * a1K - 6.0 * a2K) * lnubar
            - 3.0 * u * a1K * (-6.0 * m_s * m_s + 6.0 * m_ud * m_ud
            + m_K * m_K * (rhomK * (8.0 - 6.0 * u + 4.0 * u2) - 3.0 * ubar * (u * ubar - 3.0 * rhopK)))
            - 3.0 / 4.0 * u * a2K * (48.0 * m_s * m_s + 96.0 * m_s * m_ud + 48.0 * m_ud * m_ud
            + m_K * m_K * (-3.0 + 6.0 * u + 6.0 * u2 - 15.0 * u3 + 6.0 * u4
            + 12.0 * (-7.0 + 12.0 * u - 10.0 * u2 + 5.0 * u3) * rhomK - 8.0 * (11.0 - 15.0 * u + 10.0 * u2) * rhopK))
            + u / f_K * (f_K * m_K * m_K * (2.0 - 3.0 * u + 2.0 * u2 + 3.0 * rhomK - 3.0 * u * rhomK + 6.0 * rhopK)
            - 6.0 * f_K * (m_s + m_ud) * (m_s + m_ud)
            + f3K * (m_s + m_ud) * (60.0 * (1.0 - 3.0 * u + 2.0 * u2) + 20.0 * lambda3K * (-1.0 + 6.0 * u - 10.0 * u2 + 5.0 * u3)
            + omega3K * (-6.0 + 60.0 * u - 180.0 * u2 + 210.0 * u3 - 84.0 * u4)));

        return psi4T4_i + psi4WW_i;
    }

    Diagnostics
    AntiKaonLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });

        return results;
    }

    template <>
    struct Implementation<KaonLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1K_0;
        UsedParameter a2K_0;

        // twist 3 parameters
        UsedParameter f3K_0;
        UsedParameter lambda3K_0;
        UsedParameter omega3K_0;

        // twist 4 parameters
        UsedParameter delta4K_0;
        UsedParameter kappa4K_0;
        UsedParameter omega4K_0;

        // mass and decay constant of the pion
        UsedParameter m_K;
        UsedParameter f_K;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1K_0(p["K::a1@1GeV"], u),
            a2K_0(p["K::a2@1GeV"], u),
            f3K_0(p["K::f3@1GeV"], u),
            lambda3K_0(p["K::lambda3@1GeV"], u),
            omega3K_0(p["K::omega3@1GeV"], u),
            delta4K_0(p["K::delta4@1GeV"], u),
            kappa4K_0(p["K::kappa4@1GeV"], u),
            omega4K_0(p["K::omega4@1GeV"], u),
            m_K(p["mass::K_u"], u),
            f_K(p["decay-constant::K_u"], u),
            _mu_c(p["QCD::mu_c"], u),
            _mu_b(p["QCD::mu_b"], u),
            _mu_t(p["QCD::mu_t"], u)
        {
        }

        inline double c_rge(const double & _mu) const
        {
            /*
             * RGE coefficient, basically
             *
             *     (alpha_s/alpha_s_0)^(1/beta_0),
             *
             * with matching between the individual n-flavor QCDs.
             */

            double mu = _mu, alpha_s_mu = model->alpha_s(mu);
            double mu_0 = 1.0, alpha_s_0 = model->alpha_s(mu_0);

            if (mu < _mu_c)
                return std::pow(alpha_s_mu / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            double alpha_s_c = model->alpha_s(_mu_c);
            double result = std::pow(alpha_s_c / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            if (mu < _mu_b)
                return result * std::pow(alpha_s_mu / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            double alpha_s_b = model->alpha_s(_mu_b);
            result *= std::pow(alpha_s_b / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            if (mu < _mu_t)
                return result * std::pow(alpha_s_mu / alpha_s_b, 1.0 / QCD::beta_function_nf_5[0]);

            throw InternalError("Implementation<KaonLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a1K(const double & mu) const
        {
            return -1.0 * a1K_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a2K(const double & mu) const
        {
            return a2K_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double muK(const double & mu) const
        {
            return m_K * m_K / (model->m_s_msbar(mu) + model->m_ud_msbar(mu) / 2.0);
        }

        double f3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_ud_msbar(mu_0) / 2.0; // swapped m_s with m_q
            const double m_q_0 = model->m_s_msbar(mu_0);

            return f3K_0 * std::pow(c_rge, 55.0 / 9.0)
                + 2.0 / 19.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 55.0 / 9.0)) * f_K * (m_s_0 + m_q_0)
                - 6.0 / 65.0 * (std::pow(c_rge, 55.0 / 9.0) - std::pow(c_rge, 68.0 / 9.0)) * f_K * (m_s_0 - m_q_0) * a1K_0
                ;
        }

        double omega3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_ud_msbar(mu_0) / 2.0; // swapped m_s with m_q
            const double m_q_0 = model->m_s_msbar(mu_0);

            return (f3K_0 * omega3K_0 * std::pow(c_rge, 104.0 / 9.0)
                + 1.0 / 170.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 104.0 / 9.0)) * f_K * (m_s_0 + m_q_0)
                - 1.0 /  10.0 * (std::pow(c_rge, 68.0 / 9.0) - std::pow(c_rge, 104.0 / 9.0)) * f_K * (m_s_0 - m_q_0) * a1K_0
                + 2.0 /  15.0 * (std::pow(c_rge, 86.0 / 9.0) - std::pow(c_rge, 104.0 / 9.0)) * f_K * (m_s_0 + m_q_0) * a2K_0
                ) / f3K(mu);
        }

        double lambda3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_ud_msbar(mu_0) / 2.0; // swapped m_s with m_q
            const double m_q_0 = model->m_s_msbar(mu_0);

            return (-f3K_0 * lambda3K_0 * std::pow(c_rge, 139.0 / 18.0)
                - 14.0 / 67.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 139.0 / 18.0)) * f_K * (m_s_0 - m_q_0)
                - 14.0 /  5.0 * (std::pow(c_rge, 68.0 / 9.0) - std::pow(c_rge, 139.0 / 18.0)) * f_K * (m_s_0 + m_q_0) * a1K_0
                - 4.0  / 11.0 * (std::pow(c_rge, 86.0 / 9.0) - std::pow(c_rge, 139.0 / 18.0)) * f_K * (m_s_0 - m_q_0) * a2K_0) / f3K(mu);
        }

        inline double eta3K(const double & mu) const
        {
            return f3K(mu) / (f_K() * muK(mu));
        }

        double delta4K(const double & mu) const
        {
            const double c_rge  = this->c_rge(mu);

            return delta4K_0 * std::pow(c_rge, 32.0 / 9.0) + 1.0 / 8.0 * m_K * m_K * (1.0 - std::pow(c_rge, 32.0 / 9.0));
        }

        double kappa4K(const double & mu) const
        {
            const double c_rge  = this->c_rge(mu);
            const double mu_0   = 1.0;
            const double m_s_0 = model->m_ud_msbar(mu_0) / 2.0; // swapped m_s with m_q
            const double m_q_0 = model->m_s_msbar(mu_0);

            return -kappa4K_0
            + 9.0 / 40.0 * a1K_0 * (std::pow(c_rge, 32.0 / 9.0) - 1.0)
            + (m_s_0 * m_s_0 - m_q_0 * m_q_0) / (2.0 * m_K * m_K) * (std::pow(c_rge, 8.0) - 1.0);
        }

        double omega4K(const double & mu) const
        {
            const double c_rge  = this->c_rge(mu);

            return 1.0 / delta4K(mu) * omega4K_0 * delta4K_0 * std::pow(c_rge, 10.0);
        }
    };

    KaonLCDAs::KaonLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<KaonLCDAs>(new Implementation<KaonLCDAs>(p, o, *this))
    {
    }

    KaonLCDAs::~KaonLCDAs()
    {
    }

    PseudoscalarLCDAs *
    KaonLCDAs::make(const Parameters & p, const Options & o)
    {
        return new KaonLCDAs(p, o);
    }

    double
    KaonLCDAs::a1(const double & mu) const
    {
        return _imp->a1K(mu);
    }

    double
    KaonLCDAs::a2(const double & mu) const
    {
        return _imp->a2K(mu);
    }

    double
    KaonLCDAs::mu3(const double & mu) const
    {
        return _imp->muK(mu);
    }

    double
    KaonLCDAs::f3(const double & mu) const
    {
        return _imp->f3K(mu);
    }

    double
    KaonLCDAs::eta3(const double & mu) const
    {
        return _imp->eta3K(mu);
    }

    double
    KaonLCDAs::lambda3(const double & mu) const
    {
        return _imp->lambda3K(mu);
    }

    double
    KaonLCDAs::omega3(const double & mu) const
    {
        return _imp->omega3K(mu);
    }

    double
    KaonLCDAs::delta4(const double & mu) const
    {
        return _imp->delta4K(mu);
    }

    double
    KaonLCDAs::kappa4(const double & mu) const
    {
        return _imp->kappa4K(mu);
    }

    double
    KaonLCDAs::omega4(const double & mu) const
    {
        return _imp->omega4K(mu);
    }

    double
    KaonLCDAs::phi(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1_3o2(1.0, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2_3o2(2.0, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1_3o2.evaluate(x);
        const double c2 = gp_2_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1K(mu) * c1 + _imp->a2K(mu) * c2);
    }

    double
    KaonLCDAs::phi3p(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(1/2)
        static const GegenbauerPolynomial gp_1_1o2(1.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_2_1o2(2.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_3_1o2(3.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_4_1o2(4.0, 1.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1_1o2.evaluate(x);
        const double c2 = gp_2_1o2.evaluate(x);
        const double c3 = gp_3_1o2.evaluate(x);
        const double c4 = gp_4_1o2.evaluate(x);

        return 1.0 + 3.0 * rhopK * (1.0 + 6.0 * a2K) - 9.0 * rhomK * a1K
            + c1 * (27.0 / 2.0 * rhopK * a1K - rhomK * (3.0 / 2.0 + 27.0 * a2K))
            + c2 * (30.0 * eta3K + 15.0 * rhopK * a2K - 3.0 * rhomK * a1K)
            + c3 * (10.0 * eta3K * lambda3K - 9.0 / 2.0 * rhomK * a2K)
            + c4 * (-3.0 * eta3K * omega3K)
            + 3.0 / 2.0 * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 3.0 / 2.0 * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u);
    }

    double
    KaonLCDAs::phi3s(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1_3o2(1.0, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2_3o2(2.0, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3_3o2(3.0, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1_3o2.evaluate(x);
        const double c2 = gp_2_3o2.evaluate(x);
        const double c3 = gp_3_3o2.evaluate(x);

        const double ubar = 1.0 - u;

        return 6.0 * u * ubar * (
            1.0 + 3.0 / 2.0 * rhopK + 15.0 * rhopK * a2K - 15.0 / 2.0 * rhomK * a1K
            + c1 * (3.0 * rhopK * a1K - 15.0 / 2.0 * rhomK * a2K)
            + c2 * (5.0 * eta3K - 1.0 / 2.0 * eta3K * omega3K + 3.0 / 2.0 * rhopK * a2K)
            + c3 * (eta3K * lambda3K))
            + 9.0 * u * ubar * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 9.0 * u * ubar * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u)
        ;
    }

    double
    KaonLCDAs::phi3s_d1(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        const double ubar = 1.0 - u, x = 2.0 * u - 1.0;
        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u;

        return -3.0 * (60.0 * eta3K * (-1.0 + 12.0 * u - 30.0 * u2 + 20.0 * u3)
            -3.0 * rhomK * (1.0 + x * (-1.0 * std::log(u) +  std::log(ubar) + a1K * (8.0 + 3.0 * (std::log(u) + std::log(ubar))))
            + 3.0 * a2K  * (7.0 - 30.0 * u + 30.0 * u2 + 2.0 * x * (std::log(ubar) - std::log(u))))
            + x * (2.0 + 3.0 * rhopK * (2.0 + std::log(u) + std::log(ubar) + 3.0 * a1K * (-3.0 + 6.0 * u - std::log(u) + std::log(ubar))
            + a2K * (22.0 - 60.0 * u + 60.0 * u2 + 6.0 * (std::log(u) + std::log(ubar))))))
            - 6.0 * eta3K * (10.0 * lambda3K * (1.0 - 20.0 * u + 90.0 * u2 - 140.0 * u3 + 70.0 * u4) + 3.0 * omega3K * (1.0 - 12.0 * u + 30.0 * u2 - 20.0 * u3))
        ;
    }

    double
    KaonLCDAs::phi4(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double omega4K  = _imp->omega4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;
        const double phi2K    = -7.0 / 20.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, ubar3 = ubar2 * ubar, lnubar = std::log(ubar);
        const double x = 2.0 * u - 1.0;

        // Twist 4 contributions
        const double phi4T4 = 200.0 / 3.0 * delta4K * u2 * ubar2 + 20.0 * u2 * ubar2 * x * (4.0 * theta1K - 5.0 * theta2K)
            + 21.0 * delta4K * omega4K * (u * ubar * (2.0 + 13.0 * u * ubar) + (2.0 * u3 * (6.0 * u2 - 15.0 * u + 10.0) * lnu)
            + (2.0 * ubar3 * (6.0 * ubar2 - 15.0 * ubar + 10.0) * lnubar))
            + 40.0 * phi2K * (u * ubar * x * (2.0 - 3.0 * u * ubar) - (2.0 * u3 * (u - 2.0) * lnu)
            + (2.0 * ubar3 * (ubar - 2.0) * lnubar));
        const double phi4WW = 16.0 / 3.0 * m_K * m_K * kappa4K * (u * ubar * x * (1.0 - 2.0 * u * ubar)
            + (5.0 * (u - 2.0) * u3 * lnu) - (5.0 * (ubar - 2.0) * ubar3 * lnubar))
            + 4.0 * f3K / f_K * (m_s + m_ud) * u * ubar * (30.0 * (1.0 - x * (m_s - m_ud) / (m_s + m_ud))
            + 10.0 * lambda3K * (x * (1.0 - u * ubar) - (m_s - m_ud) / (m_s + m_ud) * (1.0 - 5.0 * u * ubar))
            - omega3K * (3.0 - 21.0 * u * ubar + 28.0 * u2 * ubar2 + 3.0 * x * (m_s - m_ud) / (m_s + m_ud) * (1.0 - 7.0 * u * ubar)))
            - 36.0 / 5.0 * m_K * m_K * a2K * (1.0 / 4.0 * u * ubar * (4.0 - 9.0 * u * ubar + 110.0 * u2 * ubar2)
            + (u3 * (10.0 - 15.0 * u + 6.0 * u2) * lnu) + (ubar3 * (10.0 - 15.0 * ubar + 6.0 * ubar2) * lnubar))
            + 4.0 * m_K * m_K * u * ubar * (1.0 + 3.0 * u * ubar) * (1.0 + 9.0 / 5.0 * a1K * x);

        return phi4T4 + phi4WW;
    }

    double
    KaonLCDAs::phi4_d1(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double omega4K  = _imp->omega4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;
        const double phi2K    = -7.0 / 20.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, u5 = u4 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, lnubar = std::log(ubar);
        const double x = 2.0 * u - 1.0;

        // Twist 4 derivatives contributions
        const double phi4T4_d1 = 20.0 * ubar2 * lnubar * (8.0 * (1.0 + 2.0 * u) * phi2K - 63.0 * u2 * omega4K * delta4K)
            + 20.0 * u2 * lnu * (-8.0 * (2.0 * u - 3.0) * phi2K + 63.0 * ubar2 * omega4K * delta4K)
            - 20.0 / 3.0 * u * (-20.0 * (1.0 - 3.0 * u + 2.0 * u2) * delta4K + 12.0 * (-8.0 + 23.0 * u - 30.0 * u2 + 15.0 * u3) * phi2K
            + 3.0 * ubar * ((1.0 - 5.0 * u + 5.0 * u2) * (8.0 * theta1K - 10.0 * theta2K) + 21.0 * x * omega4K * delta4K));
        const double phi4WW_d1 = 36.0 / 5.0 * m_K * m_K * a1K * (-1.0 + 30.0 * u2 - 60.0 * u3 + 30.0 * u4)
            - 54.0 * m_K * m_K * a2K * u * ubar * (-1.0 + 13.0 * u - 33.0 * u2 + 22.0 * u3)
            - 16.0 / 3.0 * m_K * m_K * kappa4K * (6.0 - 15.0 * u + 35.0 * u2 - 40.0 * u3 + 20.0 * u4)
            + 8.0 / 3.0 * m_K * m_K * (u2 * lnu * (-81.0 * ubar2 * a2K + 20.0 * (-3.0 + 2.0 * u) * kappa4K)
            + ubar2 * lnubar * (81.0 * u2 * a2K - 20.0 * (2.0 * u + 1.0) * kappa4K))
            + 4.0 / f_K * (f_K * m_K * m_K * (1.0 + 4.0 * u - 18.0 * u2 + 12.0 * u3)
            + f3K * (60.0 * (m_s * (1.0 - 4.0 * u + 3.0 * u2) + m_ud * u * (2.0 - 3.0 * u))
            - 20.0 * lambda3K * (m_s * (1.0 - 10.0 * u + 24.0 * u2 - 20.0 * u3 + 5.0 * u4) + m_ud * u * (2.0 - 6.0 * u + 5.0 * u3))
            + omega3K * m_s * (-12.0 * u + 60.0 * u2 - 210.0 * u4 + 168.0 * u5)
            + omega3K * m_ud * (-6.0 + 108.0 * u - 480.0 * u2 + 840.0 * u3 - 630.0 * u4 + 168.0 * u5)));

        return phi4T4_d1 + phi4WW_d1;
    }

    double
    KaonLCDAs::phi4_d2(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double omega4K  = _imp->omega4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;
        const double phi2K    = -7.0 / 20.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        // Twist 4 derivatives contributions
        const double phi4T4_d2 = 400.0 / 3.0 * (1.0 - 6.0 * u + 6.0 * u2) * delta4K
            - 20.0 * (24.0 * phi2K * (-1.0 + 7.0 * u - 15.0 * u2 + 10.0 * u3)
            + (-1.0 + 12.0 * u - 30.0 * u2 + 20.0 * u3) * (- 8.0 * theta1K + 10.0 * theta2K)
            - 21.0 * omega4K * delta4K * (1.0 - 3.0 * u + 3.0 * u2))
            + 120.0 * u * (8.0 * ubar * phi2K + 21.0 * (1.0 - 3.0 * u + 2.0 * u2) * omega4K * delta4K) * (lnu - lnubar);
        const double phi4WW_d2 = m_K * m_K * (432.0 * u * (1.0 - 3.0 * u + 2.0 * u2) * a1K
            + 54.0 * (1.0 - 32.0 * u + 142.0 * u2 - 220.0 * u3 + 110.0 * u4) * a2K
            - 80.0 / 3.0 * (-5.0 + 18.0 * u - 24.0 * u2 + 16.0 * u3) * kappa4K
            + 16.0 * u * (27.0 * (1.0 - 3.0 * u + 2.0 * u2) * a2K + 20.0 * ubar * kappa4K) * (-lnu + lnubar))
            + 16.0 / f_K * (f_K * m_K * m_K * (1.0 - 9.0 * u + 9.0 * u2)
            + f3K * (30.0 * (m_s * (-2.0 + 3.0 * u) + m_ud * (1.0 - 3.0 * u))
            - 10.0 * lambda3K * (m_s * (-5.0 + 24.0 * u - 30.0 * u2 + 10.0 * u3) + m_ud * (1.0 - 6.0 * u + 10.0 * u3))
            + omega3K * (m_s * (-3.0 + 30.0 * u - 210.0 * u3 + 210.0 * u4)
            + m_ud * (27.0 - 240.0 * u + 630.0 * u2 - 630.0 * u3 + 210.0 * u4))));

        return phi4T4_d2 + phi4WW_d2;
    }

    double
    KaonLCDAs::psi4(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;

        const double lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        // Gegenbauer polynomials C_n^(1/2)
        static const GegenbauerPolynomial gp_1_1o2(1.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_2_1o2(2.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_3_1o2(3.0, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_4_1o2(4.0, 1.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c0 = 1.0;
        const double c1 = gp_1_1o2.evaluate(x);
        const double c2 = gp_2_1o2.evaluate(x);
        const double c3 = gp_3_1o2.evaluate(x);
        const double c4 = gp_4_1o2.evaluate(x);

        // Twist 4 contributions
        const double psi4T4 = 20.0 / 3.0 * delta4K * c2 + 5.0 * (5.0 * theta1K - theta2K) * c3;
        const double psi4WW = c0 * m_K * m_K * (1.0 + 6.0 * rhopK * (1.0 + 6.0 * a2K) - 18.0 * rhomK * a1K)
            + c1 * m_K * m_K * (-12.0 * kappa4K - 9.0 / 5.0 * a1K + 27.0 * rhopK * a1K - 3.0 * rhomK * (1.0 + 18.0 * a2K))
            + c2 * (m_K * m_K * (1.0 + 18.0 / 7.0 * a2K + 30.0 * rhopK * a2K - 6.0 * rhomK * a1K) + 60.0 * f3K / f_K * (m_s + m_ud))
            + c3 * (m_K * m_K * (9.0 / 5.0 * a1K + 16.0 / 3.0 * kappa4K - 9.0 * rhomK * a2K) + 20.0 * f3K * lambda3K / f_K * (m_s + m_ud))
            + c4 * (-9.0 / 28.0 * m_K * m_K * a2K - 6.0 * f3K * omega3K / f_K * (m_s + m_ud))
            + 6.0 * m_ud * (m_s + m_ud) * (1.0 + 3.0 * a1K + 6.0 * a2K) * lnubar
            + 6.0 * m_s * (m_s + m_ud) * (1.0 - 3.0 * a1K + 6.0 * a2K) * lnu;

        return psi4T4 + psi4WW;
    }

    double
    KaonLCDAs::psi4_i(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s  = _imp->model->m_ud_msbar(mu) / 2.0;
        const double m_ud = _imp->model->m_s_msbar(mu); // swapped m_s with m_ud

        const double m_K  = _imp->m_K;
        const double f_K  = _imp->f_K;

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>((m_s + m_ud) / _imp->m_K); // EOM constraints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = (m_s * m_s - m_ud * m_ud) / power_of<2>(_imp->m_K); // identical in the limit m_q -> 0
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);
        const double f3K      = _imp->f3K(mu);

        // Twist 4 coefficients
        const double delta4K  = _imp->delta4K(mu);
        const double kappa4K  = _imp->kappa4K(mu);
        const double theta1K  = 7.0 / 10.0 * a1K * delta4K;
        const double theta2K  = -7.0 / 5.0 * a1K * delta4K;

        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        // Twist 4 contributions
        const double psi4T4_i = -5.0 / 3.0 * u * ubar * (delta4K * (8.0 * u - 4.0)
            + 3.0 * (1.0 - 5.0 * u + 5.0 * u2) * (5.0 * theta1K - theta2K));
        const double psi4WW_i = 20.0 / 3.0 * m_K * m_K * kappa4K * u * (1.0 + 3.0 * u - 8.0 * u2 + 4.0 * u3)
            - 6.0 * m_s * (m_s + m_ud) * u * (-1.0 + 3.0 * a1K - 6.0 * a2K) * lnu
            + 6.0 * m_ud * (m_s + m_ud) * ubar * (-1.0 - 3.0 * a1K - 6.0 * a2K) * lnubar
            - 3.0 * u * a1K * (-6.0 * m_s * m_s + 6.0 * m_ud * m_ud
            + m_K * m_K * (rhomK * (8.0 - 6.0 * u + 4.0 * u2) - 3.0 * ubar * (u * ubar - 3.0 * rhopK)))
            - 3.0 / 4.0 * u * a2K * (48.0 * m_s * m_s + 96.0 * m_s * m_ud + 48.0 * m_ud * m_ud
            + m_K * m_K * (-3.0 + 6.0 * u + 6.0 * u2 - 15.0 * u3 + 6.0 * u4
            + 12.0 * (-7.0 + 12.0 * u - 10.0 * u2 + 5.0 * u3) * rhomK - 8.0 * (11.0 - 15.0 * u + 10.0 * u2) * rhopK))
            + u / f_K * (f_K * m_K * m_K * (2.0 - 3.0 * u + 2.0 * u2 + 3.0 * rhomK - 3.0 * u * rhomK + 6.0 * rhopK)
            - 6.0 * f_K * (m_s + m_ud) * (m_s + m_ud)
            + f3K * (m_s + m_ud) * (60.0 * (1.0 - 3.0 * u + 2.0 * u2) + 20.0 * lambda3K * (-1.0 + 6.0 * u - 10.0 * u2 + 5.0 * u3)
            + omega3K * (-6.0 + 60.0 * u - 180.0 * u2 + 210.0 * u3 - 84.0 * u4)));

        return psi4T4_i + psi4WW_i;
    }

    Diagnostics
    KaonLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });

        return results;
    }
}
