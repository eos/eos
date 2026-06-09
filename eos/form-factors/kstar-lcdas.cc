/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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

#include <eos/form-factors/kstar-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<KstarLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (vector) Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a_1_para_0;
        UsedParameter a_2_para_0;
        UsedParameter f_para;

        // twist 2 (tensor) Gegenbauer coefficients and normalization at mu = 1 GeV
        UsedParameter a_1_perp_0;
        UsedParameter a_2_perp_0;
        UsedParameter f_perp_0;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a_1_para_0(p["K^*::a_1_para@1GeV"], u),
            a_2_para_0(p["K^*::a_2_para@1GeV"], u),
            f_para(p["K^*::f_para"], u),
            a_1_perp_0(p["K^*::a_1_perp@1GeV"], u),
            a_2_perp_0(p["K^*::a_2_perp@1GeV"], u),
            f_perp_0(p["K^*::f_perp@1GeV"], u),
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

            throw InternalError("Implementation<KstarLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a_1_para(const double & mu) const
        {
            return a_1_para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a_2_para(const double & mu) const
        {
            return a_2_para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a_1_perp(const double & mu) const
        {
            return a_1_perp_0 * std::pow(c_rge(mu), 36.0 / 9.0);
        }

        inline double a_2_perp(const double & mu) const
        {
            return a_2_perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

        inline double f_perp(const double & mu) const
        {
            // gamma_0 / (beta_0^Nf=3) = 4 / 23, see [BFS:2001A], p. 14, below eq. (48)
            return f_perp_0 * std::pow(c_rge(mu), +4.0 / 23.0 * QCD::beta_function_nf_3[0]);
        }
    };

    KstarLCDAs::KstarLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<KstarLCDAs>(new Implementation<KstarLCDAs>(p, o, *this))
    {
    }

    KstarLCDAs::~KstarLCDAs()
    {
    }

    double
    KstarLCDAs::a_1_para(const double & mu) const
    {
        return _imp->a_1_para(mu);
    }

    double
    KstarLCDAs::a_2_para(const double & mu) const
    {
        return _imp->a_2_para(mu);
    }

    double
    KstarLCDAs::f_para() const
    {
        return _imp->f_para();
    }

    double
    KstarLCDAs::a_1_perp(const double & mu) const
    {
        return _imp->a_1_perp(mu);
    }

    double
    KstarLCDAs::a_2_perp(const double & mu) const
    {
        return _imp->a_2_perp(mu);
    }

    double
    KstarLCDAs::f_perp(const double & mu) const
    {
        return _imp->f_perp(mu);
    }

    double
    KstarLCDAs::phi_2_para(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c1 = 3.0 * x;
        const double c2 = (15.0 * x2 - 3.0) / 2.0;

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a_1_para(mu) * c1 + _imp->a_2_para(mu) * c2);
    }

    double
    KstarLCDAs::phi_2_perp(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c1 = 3.0 * x;
        const double c2 = (15.0 * x2 - 3.0) / 2.0;

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a_1_perp(mu) * c1 + _imp->a_2_perp(mu) * c2);
    }

    double
    KstarLCDAs::part_inv_mom_phi_2_para(const double & u, const double & mu) const
    {
        // \int_0^u phi_2_para(v) / (1 - v)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c0 = 3.0;
        const double c1 = 6.0 * x - 3.0;
        const double c2 = (15.0 * x2 - 10.0 * x - 1.0) * 3.0 / 4.0;

        return (1.0 + x) * (1.0 + x) / 4.0 * (c0 + _imp->a_1_para(mu) * c1 + _imp->a_2_para(mu) * c2);
    }

    Diagnostics
    KstarLCDAs::diagnostics() const
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
