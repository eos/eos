/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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

#include <eos/form-factors/rho-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<RhoLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) para Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a2para_0;
        UsedParameter a4para_0;
        UsedParameter fpara;

        // twist 2 (even) perp Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a2perp_0;
        UsedParameter a4perp_0;
        UsedParameter fperp_0;

        // mass and decay constant of the rho
        UsedParameter m_rho;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a2para_0(p["rho::a2para@1GeV"], u),
            a4para_0(p["rho::a4para@1GeV"], u),
            fpara(p["rho::fpara"], u),
            a2perp_0(p["rho::a2perp@1GeV"], u),
            a4perp_0(p["rho::a4perp@1GeV"], u),
            fperp_0(p["rho::fperp@1GeV"], u),
            m_rho(p["mass::rho^+"], u),
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

            throw InternalError("Implementation<RhoLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a2para(const double & mu) const
        {
            return a2para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a4para(const double & mu) const
        {
            return a4para_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double a2perp(const double & mu) const
        {
            return a2perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

        inline double a4perp(const double & mu) const
        {
            return a4perp_0 * std::pow(c_rge(mu), 368.0 / 45.0);
        }

        inline double fperp(const double & mu) const
        {
            // [BBKT1998A], p. 23, eq. (3.59)
            return fperp_0 * std::pow(c_rge(mu), +4.0 / 3.0);
        }
    };

    RhoLCDAs::RhoLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<RhoLCDAs>(new Implementation<RhoLCDAs>(p, o, *this)),
        gp_2_3o2(2, 3.0 / 2.0),
        gp_4_3o2(4, 3.0 / 2.0)
    {
    }

    RhoLCDAs::~RhoLCDAs()
    {
    }

    VectorLCDAs *
    RhoLCDAs::make(const Parameters & p, const Options & o)
    {
        return new RhoLCDAs(p, o);
    }

    double
    RhoLCDAs::a2para(const double & mu) const
    {
        return _imp->a2para(mu);
    }

    double
    RhoLCDAs::a4para(const double & mu) const
    {
        return _imp->a4para(mu);
    }

    double
    RhoLCDAs::fpara() const
    {
        return _imp->fpara();
    }
    double
    RhoLCDAs::a2perp(const double & mu) const
    {
        return _imp->a2perp(mu);
    }

    double
    RhoLCDAs::a4perp(const double & mu) const
    {
        return _imp->a4perp(mu);
    }

    double
    RhoLCDAs::fperp(const double &mu) const
    {
        return _imp->fperp(mu);
    }

    double
    RhoLCDAs::phipara(const double & u, const double & mu) const
    {
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);
        const double c4 = gp_4_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a2para(mu) * c2 + _imp->a4para(mu) * c4);
    }

    double
    RhoLCDAs::phiperp(const double & u, const double & mu) const
    {
        const double x = 2.0 * u - 1.0;
        const double c2 = gp_2_3o2.evaluate(x);
        const double c4 = gp_4_3o2.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a2perp(mu) * c2 + _imp->a4perp(mu) * c4);
    }

    Diagnostics
    RhoLCDAs::diagnostics() const
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
