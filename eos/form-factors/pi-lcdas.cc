/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/utils/model.hh>
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
        UsedParameter omega3pi_0;

        // twist 4 parameters
        UsedParameter deltapipi_0;
        UsedParameter omega4pi_0;

        // mass and decay constant of the pion
        UsedParameter m_pi;
        UsedParameter f_pi;

        // matching scales for the individual n-flavour effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a2pi_0(p["pi::a2@1GeV"], u),
            a4pi_0(p["pi::a4@1GeV"], u),
            f3pi_0(p["pi::f3@1GeV"], u),
            omega3pi_0(p["pi::omega3@1GeV"], u),
            deltapipi_0(p["pi::delta^2@1GeV"], u),
            omega4pi_0(p["pi::omega4@1GeV"], u),
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
             * with matching between the individual n-flavour QCDs.
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

        inline double mupi(const double & mu) const
        {
            return m_pi * m_pi / this->m_ud_msbar(mu);
        }

        double f3pi(const double & mu) const
        {
            return f3pi_0 * std::pow(c_rge(mu), 55.0 / 9.0);
        }

        inline double eta3pi(const double & mu) const
        {
            return f3pi(mu) / (f_pi() * mupi(mu));
        }

        double omega3pi(const double & mu) const
        {
            return omega3pi_0 * std::pow(c_rge(mu), 49.0 / 9.0);
        }

        double deltapipi(const double & mu) const
        {
            return deltapipi_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        double omega4pi(const double & mu) const
        {
            return omega4pi_0 * std::pow(c_rge(mu), 58.0 / 9.0);
        }
    };

    PionLCDAs::PionLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<PionLCDAs>(new Implementation<PionLCDAs>(p, o, *this))
    {
    }

    PionLCDAs::~PionLCDAs()
    {
    }

    double
    PionLCDAs::a2pi(const double & mu) const
    {
        return _imp->a2pi(mu);
    }

    double
    PionLCDAs::a4pi(const double & mu) const
    {
        return _imp->a4pi(mu);
    }

    double
    PionLCDAs::mupi(const double & mu) const
    {
        return _imp->mupi(mu);
    }

    double
    PionLCDAs::f3pi(const double & mu) const
    {
        return _imp->f3pi(mu);
    }

    double
    PionLCDAs::eta3pi(const double & mu) const
    {
        return _imp->eta3pi(mu);
    }

    double
    PionLCDAs::omega3pi(const double & mu) const
    {
        return _imp->omega3pi(mu);
    }

    double
    PionLCDAs::deltapipi(const double & mu) const
    {
        return _imp->deltapipi(mu);
    }

    double
    PionLCDAs::omega4pi(const double & mu) const
    {
        return _imp->omega4pi(mu);
    }

    double
    PionLCDAs::phi(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        const double x = 2.0 * u - 1.0, x2 = x * x, x4 = x2 * x2;
        const double c2 = (15.0 * x2 - 3.0) / 2.0;
        const double c4 = (15.0 - 210.0 * x2 + 315.0 * x4) / 8.0;

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a2pi(mu) * c2 + _imp->a4pi(mu) * c4);
    }

    double
    PionLCDAs::phi3p(const double & u, const double & mu) const
    {
        // Setting lambda3pi and rhopi to zero.
        const double eta3pi = _imp->eta3pi(mu);
        const double omega3pi = _imp->omega3pi(mu);

        // Gegenbauer polynomials C_n^(1/2)
        const double x = 2.0 * u - 1.0, x2 = x * x, x4 = x2 * x2;
        const double c2 = (3.0 * x2 - 1.0) / 2.0;
        const double c4 = (35.0 * x4 - 30.0 * x2 + 3.0) / 8.0;

        return 1.0 + 30.0 * eta3pi * c2 - 3.0 * eta3pi * _imp->omega3pi(mu) * c4;
    }

    double
    PionLCDAs::phi3s(const double & u, const double & mu) const
    {
        // Setting lambda3pi and rhopi to zero.

        const double mupi = _imp->m_pi * _imp->m_pi / _imp->m_ud_msbar(mu);
        const double f3pi = _imp->f3pi(mu);
        const double eta3pi = _imp->eta3pi(mu);

        // Gegenbauer polynomials C_n^(3/2)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c2 = (15.0 * x2 - 3.0) / 2.0;

        return 6.0 * u * (1.0 - u) * (1.0 + 0.5 * eta3pi * (10.0 - _imp->omega3pi(mu)) * c2);
    }

    double
    PionLCDAs::phi3s_d1(const double & u, const double & mu) const
    {
        // Setting lambda3pi and rhopi to zero.
        //
        const double mupi = _imp->m_pi * _imp->m_pi / _imp->m_ud_msbar(mu);
        const double f3pi = _imp->f3pi(mu);
        const double eta3pi = _imp->eta3pi(mu);
        const double omega3pi = _imp->omega3pi(mu);

        // Gegenbauer polynomials C_n^(3/2)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c2 = (15.0 * x2 - 3.0) / 2.0;

        return -6.0 * x * (1.0 + 0.5 * eta3pi * (10.0 - omega3pi) * c2)
            + 180.0 * u * (1.0 - u) * 0.5 * eta3pi * (10.0 - omega3pi) * x;
    }

    double
    PionLCDAs::phi4(const double & u, const double & mu) const
    {
        const double u2 = u * u, u3 = u2 * u, lnu = std::log(u);
        const double ubar = 1.0 - u, ubar2 = ubar * ubar, ubar3 = ubar2 * ubar, lnubar = std::log(ubar);

        return _imp->deltapipi(mu) * (200.0 / 3.0 * u2 * ubar2 + 21.0 * _imp->omega4pi(mu) * (
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

        return _imp->deltapipi(mu) * (400.0 / 3.0 * u * (1.0 - 3.0 * u + 2.0 * u2) + 21.0 * _imp->omega4pi(mu) * (
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

        return 20.0 / 3.0 * _imp->deltapipi(mu) * (
                20.0 * (1.0 - 6.0 * u + 6.0 * u2)
                - 63.0 * (
                    -1.0 + 3.0 * u - 3.0 * u2
                    + 6.0 * u * (1.0 - 3.0 * u + 2.0 * u2) * (lnubar - lnu)
                ) * _imp->omega4pi(mu)
            );
    }

    double
    PionLCDAs::psi4(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(1/2)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c2 = (3.0 * x2 - 1.0) / 2.0;

        return _imp->deltapipi(mu) * 20.0 / 3.0 * c2;
    }

    double
    PionLCDAs::psi4_i(const double & u, const double & mu) const
    {
        const double u2 = u * u;

        return _imp->deltapipi(mu) * 20.0 / 3.0 * u * (1.0 - 3.0 * u + 2.0 * u2);
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
