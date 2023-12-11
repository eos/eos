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

#include <eos/form-factors/k-star-lcdas.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<AntiKStarLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) para Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1para_0;
        UsedParameter a2para_0;
        UsedParameter a3para_0;
        UsedParameter a4para_0;
        UsedParameter fpara;

        // twist 2 (tensor) Gegenbauer coefficients and normalization at mu = 1 GeV
        UsedParameter a1perp_0;
        UsedParameter a2perp_0;
        UsedParameter a3perp_0;
        UsedParameter a4perp_0;
        UsedParameter fperp_0;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1para_0(p["K^*::a1para@1GeV"], u),
            a2para_0(p["K^*::a2para@1GeV"], u),
            a3para_0(p["K^*::a3para@1GeV"], u),
            a4para_0(p["K^*::a4para@1GeV"], u),
            fpara(p["K^*::fpara"], u),
            a1perp_0(p["K^*::a1perp@1GeV"], u),
            a2perp_0(p["K^*::a2perp@1GeV"], u),
            a3perp_0(p["K^*::a3perp@1GeV"], u),
            a4perp_0(p["K^*::a4perp@1GeV"], u),
            fperp_0(p["K^*::fperp@1GeV"], u),
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

            throw InternalError("Implementation<AntiKStarLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a1para(const double & mu) const
        {
            return a1para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a2para(const double & mu) const
        {
            return a2para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a3para(const double & mu) const
        {
            return a3para_0 * std::pow(c_rge(mu), 314.0 / 45.0);
        }

        inline double a4para(const double & mu) const
        {
            return a4para_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double a1perp(const double & mu) const
        {
            return a1perp_0 * std::pow(c_rge(mu), 36.0 / 9.0);
        }

        inline double a2perp(const double & mu) const
        {
            return a2perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

         inline double a3perp(const double & mu) const
        {
            return a3perp_0 * std::pow(c_rge(mu), 64.0 / 9.0);
        }

        inline double a4perp(const double & mu) const
        {
            return a4perp_0 * std::pow(c_rge(mu), 368.0 / 45.0);
        }

        inline double fperp(const double & mu) const
        {
            // gamma_0 / (beta_0^Nf=3) = 4 / 23, see [BFS2001], p. 14, below eq. (48)
            return fperp_0 * std::pow(c_rge(mu), +4.0 / 23.0 * QCD::beta_function_nf_3[0]);
        }
    };

    AntiKStarLCDAs::AntiKStarLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AntiKStarLCDAs>(new Implementation<AntiKStarLCDAs>(p, o, *this))
    {
    }

    AntiKStarLCDAs::~AntiKStarLCDAs()
    {
    }

    VectorLCDAs *
    AntiKStarLCDAs::make(const Parameters & p, const Options & o)
    {
        return new AntiKStarLCDAs(p, o);
    }

    double
    AntiKStarLCDAs::a1para(const double & mu) const
    {
        return _imp->a1para(mu);
    }

    double
    AntiKStarLCDAs::a2para(const double & mu) const
    {
        return _imp->a2para(mu);
    }

    double
    AntiKStarLCDAs::a3para(const double & mu) const
    {
        return _imp->a3para(mu);
    }

    double
    AntiKStarLCDAs::a4para(const double & mu) const
    {
        return _imp->a4para(mu);
    }

    double
    AntiKStarLCDAs::fpara() const
    {
        return _imp->fpara();
    }

    double
    AntiKStarLCDAs::a1perp(const double & mu) const
    {
        return _imp->a1perp(mu);
    }

    double
    AntiKStarLCDAs::a2perp(const double & mu) const
    {
        return _imp->a2perp(mu);
    }

    double
    AntiKStarLCDAs::a3perp(const double & mu) const
    {
        return _imp->a3perp(mu);
    }

    double
    AntiKStarLCDAs::a4perp(const double & mu) const
    {
        return _imp->a4perp(mu);
    }

    double
    AntiKStarLCDAs::fperp(const double & mu) const
    {
        return _imp->fperp(mu);
    }

    double
    AntiKStarLCDAs::phipara(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1para(mu) * c1 + _imp->a2para(mu) * c2 + _imp->a3para(mu) * c3 + _imp->a4para(mu) * c4);
    }

    double
    AntiKStarLCDAs::phiperp(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1perp(mu) * c1 + _imp->a2perp(mu) * c2 + _imp->a3perp(mu) * c3 + _imp->a4perp(mu) * c4);
    }

    Diagnostics
    AntiKStarLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(5.0), "RGE coefficient C(mu = 5.0 GeV)" });

        return results;
    }

    template <>
    struct Implementation<KStarLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 (even) para Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1para_0;
        UsedParameter a2para_0;
        UsedParameter a3para_0;
        UsedParameter a4para_0;
        UsedParameter fpara;

        // twist 2 (tensor) Gegenbauer coefficients and normalization at mu = 1 GeV
        UsedParameter a1perp_0;
        UsedParameter a2perp_0;
        UsedParameter a3perp_0;
        UsedParameter a4perp_0;
        UsedParameter fperp_0;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1para_0(p["K^*::a1para@1GeV"], u),
            a2para_0(p["K^*::a2para@1GeV"], u),
            a3para_0(p["K^*::a3para@1GeV"], u),
            a4para_0(p["K^*::a4para@1GeV"], u),
            fpara(p["K^*::fpara"], u),
            a1perp_0(p["K^*::a1perp@1GeV"], u),
            a2perp_0(p["K^*::a2perp@1GeV"], u),
            a3perp_0(p["K^*::a3perp@1GeV"], u),
            a4perp_0(p["K^*::a4perp@1GeV"], u),
            fperp_0(p["K^*::fperp@1GeV"], u),
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

            throw InternalError("Implementation<KStarLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a1para(const double & mu) const
        {
            return -1.0 * a1para_0 * std::pow(c_rge(mu), 32.0 / 9.0);
        }

        inline double a2para(const double & mu) const
        {
            return a2para_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a3para(const double & mu) const
        {
            return -1.0 * a3para_0 * std::pow(c_rge(mu), 314.0 / 45.0);
        }

        inline double a4para(const double & mu) const
        {
            return a4para_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double a1perp(const double & mu) const
        {
            return a1perp_0 * std::pow(c_rge(mu), 36.0 / 9.0);
        }

        inline double a2perp(const double & mu) const
        {
            return a2perp_0 * std::pow(c_rge(mu), 52.0 / 9.0);
        }

         inline double a3perp(const double & mu) const
        {
            return a3perp_0 * std::pow(c_rge(mu), 64.0 / 9.0);
        }

        inline double a4perp(const double & mu) const
        {
            return a4perp_0 * std::pow(c_rge(mu), 368.0 / 45.0);
        }

        inline double fperp(const double & mu) const
        {
            // gamma_0 / (beta_0^Nf=3) = 4 / 23, see [BFS2001], p. 14, below eq. (48)
            return fperp_0 * std::pow(c_rge(mu), +4.0 / 23.0 * QCD::beta_function_nf_3[0]);
        }
    };

    KStarLCDAs::KStarLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<KStarLCDAs>(new Implementation<KStarLCDAs>(p, o, *this))
    {
    }

    KStarLCDAs::~KStarLCDAs()
    {
    }

    VectorLCDAs *
    KStarLCDAs::make(const Parameters & p, const Options & o)
    {
        return new KStarLCDAs(p, o);
    }

    double
    KStarLCDAs::a1para(const double & mu) const
    {
        return _imp->a1para(mu);
    }

    double
    KStarLCDAs::a2para(const double & mu) const
    {
        return _imp->a2para(mu);
    }

    double
    KStarLCDAs::a3para(const double & mu) const
    {
        return _imp->a3para(mu);
    }

    double
    KStarLCDAs::a4para(const double & mu) const
    {
        return _imp->a4para(mu);
    }

    double
    KStarLCDAs::fpara() const
    {
        return _imp->fpara();
    }

    double
    KStarLCDAs::a1perp(const double & mu) const
    {
        return _imp->a1perp(mu);
    }

    double
    KStarLCDAs::a2perp(const double & mu) const
    {
        return _imp->a2perp(mu);
    }

    double
    KStarLCDAs::a3perp(const double & mu) const
    {
        return _imp->a3perp(mu);
    }

    double
    KStarLCDAs::a4perp(const double & mu) const
    {
        return _imp->a4perp(mu);
    }

    double
    KStarLCDAs::fperp(const double & mu) const
    {
        return _imp->fperp(mu);
    }

    double
    KStarLCDAs::phipara(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1para(mu) * c1 + _imp->a2para(mu) * c2 + _imp->a3para(mu) * c3 + _imp->a4para(mu) * c4);
    }

    double
    KStarLCDAs::phiperp(const double & u, const double & mu) const
    {
         // Gegenbauer polynomials C_n^(3/2)
        static const GegenbauerPolynomial gp_1(1, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_2(2, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_3(3, 3.0 / 2.0);
        static const GegenbauerPolynomial gp_4(4, 3.0 / 2.0);

        const double x = 2.0 * u - 1.0;
        const double c1 = gp_1.evaluate(x);
        const double c2 = gp_2.evaluate(x);
        const double c3 = gp_3.evaluate(x);
        const double c4 = gp_4.evaluate(x);

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1perp(mu) * c1 + _imp->a2perp(mu) * c2 + _imp->a3perp(mu) * c3 + _imp->a4perp(mu) * c4);
    }

    Diagnostics
    KStarLCDAs::diagnostics() const
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
