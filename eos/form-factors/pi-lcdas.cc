/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2024 Danny van Dyk
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

#include <eos/form-factors/pi-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<PionLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a2pi_0;
        UsedParameter a4pi_0;

        // twist 3 parameters
        UsedParameter f3pi_0;
        UsedParameter omega3_0;

        // twist 4 parameters
        UsedParameter delta4_0;
        UsedParameter omega4_0;

        // mass and decay constant of the pion
        UsedParameter m_pi;
        UsedParameter f_pi;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a2pi_0(p["pi::a2@1GeV"], u),
            a4pi_0(p["pi::a4@1GeV"], u),
            f3pi_0(p["pi::f3@1GeV"], u),
            omega3_0(p["pi::omega3@1GeV"], u),
            delta4_0(p["pi::delta4@1GeV"], u),
            omega4_0(p["pi::omega4@1GeV"], u),
            m_pi(p["mass::pi^+"], u),
            f_pi(p["decay-constant::pi"], u),
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

            throw InternalError("Implementation<PionLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a2pi(const double & mu) const
        {
            return a2pi_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a4pi(const double & mu) const
        {
            return a4pi_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double m_ud_msbar(const double & mu) const
        {
            return this->model->m_ud_msbar(mu);
        }

        inline double mu3(const double & mu) const
        {
            return m_pi * m_pi / this->m_ud_msbar(mu);
        }

        double f3(const double & mu) const
        {
            return f3pi_0 * std::pow(c_rge(mu), 55.0 / 9.0);
        }

        inline double eta3(const double & mu) const
        {
            return f3(mu) / (f_pi() * this->mu3(mu));
        }

        double omega3(const double & mu) const
        {
            return omega3_0 * std::pow(c_rge(mu), 49.0 / 9.0);
        }

        double delta4(const double & mu) const
        {
            return delta4_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        double omega4(const double & mu) const
        {
            return omega4_0 * std::pow(c_rge(mu), 58.0 / 9.0);
        }
    };

    PionLCDAs::PionLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<PionLCDAs>(new Implementation<PionLCDAs>(p, o, *this))
    {
    }

    PionLCDAs::~PionLCDAs()
    {
    }

    PseudoscalarLCDAs *
    PionLCDAs::make(const Parameters & p, const Options & o)
    {
        return new PionLCDAs(p, o);
    }

    double
    PionLCDAs::a2(const double & mu) const
    {
        return _imp->a2pi(mu);
    }

    double
    PionLCDAs::a4(const double & mu) const
    {
        return _imp->a4pi(mu);
    }

    double
    PionLCDAs::mu3(const double & mu) const
    {
        return _imp->mu3(mu);
    }

    double
    PionLCDAs::f3(const double & mu) const
    {
        return _imp->f3(mu);
    }

    double
    PionLCDAs::eta3(const double & mu) const
    {
        return _imp->eta3(mu);
    }

    double
    PionLCDAs::omega3(const double & mu) const
    {
        return _imp->omega3(mu);
    }

    double
    PionLCDAs::delta4(const double & mu) const
    {
        return _imp->delta4(mu);
    }

    double
    PionLCDAs::omega4(const double & mu) const
    {
        return _imp->omega4(mu);
    }

    double
    PionLCDAs::phi(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_2_3o2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4_3o2(4, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);
        const double c4 = gp_4_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a2pi(mu) * c2 + _imp->a4pi(mu) * c4);
    }

    double
    PionLCDAs::phi3p(const double & u, const double & mu) const
    {
        // Setting lambda3pi and rhopi to zero.
        const double eta3 = _imp->eta3(mu);
        const double omega3 = _imp->omega3(mu);

        // Gegenbauer polynomials C_n^(1/2)
        static const GegenbauerPolynomial gp_2_1o2(2, 1.0 / 2.0);
        static const GegenbauerPolynomial gp_4_1o2(4, 1.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_1o2.evaluate(x);
        const double c4 = gp_4_1o2.evaluate(x);

        return 1.0 + 30.0 * eta3 * c2 - 3.0 * eta3 * omega3 * c4;
    }

    double
    PionLCDAs::phi3s(const double & u, const double & mu) const
    {
        // Setting lambda3pi and rhopi to zero.
        const double eta3 = _imp->eta3(mu);
        const double omega3 = _imp->omega3(mu);

        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_2_3o2(2, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + 0.5 * eta3 * (10.0 - omega3) * c2);
    }

    double
    PionLCDAs::phi3s_d1(const double & u, const double & mu) const
    {
        // Setting lambda3pi and rhopi to zero.
        const double eta3 = _imp->eta3(mu);
        const double omega3 = _imp->omega3(mu);

        // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_2_3o2(2, 3.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);

        return -6.0 * x * (1.0 + 0.5 * eta3 * (10.0 - omega3) * c2)
            + 180.0 * u * (1.0 - u) * 0.5 * eta3 * (10.0 - omega3) * x;
    }

    double
    PionLCDAs::phi4(const double & u, const double & mu) const
    {
        const double u2 = u * u, u3 = u2 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, ubar3 = ubar2 * ubar, lnubar = std::log(ubar);

        return _imp->delta4(mu) * (200.0 / 3.0 * u2 * ubar2 + 21.0 * _imp->omega4(mu) * (
                u * ubar * (2.0 + 13.0 * u * ubar)
                + 2.0 * u3    * (6.0 * u2    - 15.0 * u    + 10.0) * lnu
                + 2.0 * ubar3 * (6.0 * ubar2 - 15.0 * ubar + 10.0) * lnubar
            ));
    }

    double
    PionLCDAs::phi4_d1(const double & u, const double & mu) const
    {
        const double u2 = u * u, u3 = u2 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, lnubar = std::log(ubar);

        return _imp->delta4(mu) * (400.0 / 3.0 * u * (1.0 - 3.0 * u + 2.0 * u2) + 21.0 * _imp->omega4(mu) * (
                2.0 + 22.0 * u - 78.0 * u2 + 52.0 * u3
                + 2.0 * u2    * (6.0 * u2 - 15.0 * u + 10.0 + 30.0 * ubar2 * lnu)
                - 2.0 * ubar2 * (6.0 * u2 +  3.0 * u +  1.0 + 30.0 * u2    * lnubar)
            ));
    }

    double
    PionLCDAs::phi4_d2(const double & u, const double & mu) const
    {
        const double u2 = u * u, lnu = std::log(u);
        const double ubar = 1.0 - u, lnubar = std::log(ubar);

        return 20.0 / 3.0 * _imp->delta4(mu) * (
                20.0 * (1.0 - 6.0 * u + 6.0 * u2)
                - 63.0 * (
                    -1.0 + 3.0 * u - 3.0 * u2
                    + 6.0 * u * (1.0 - 3.0 * u + 2.0 * u2) * (lnubar - lnu)
                ) * _imp->omega4(mu)
            );
    }

    double
    PionLCDAs::psi4(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(1/2)
        static const GegenbauerPolynomial gp_2_1o2(2, 1.0 / 2.0);
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_1o2.evaluate(x);

        return _imp->delta4(mu) * 20.0 / 3.0 * c2;
    }

    double
    PionLCDAs::psi4_i(const double & u, const double & mu) const
    {
        const double u2 = u * u;

        return _imp->delta4(mu) * 20.0 / 3.0 * u * (1.0 - 3.0 * u + 2.0 * u2);
    }

    Diagnostics
    PionLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(5.0), "RGE coefficient C(mu = 5.0 GeV)" });

        return results;
    }
}
