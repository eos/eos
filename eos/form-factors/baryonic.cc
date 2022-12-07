/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2017 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
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

#include <eos/form-factors/baryonic-impl.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/form-factors/parametric-abr2022.hh>
#include <eos/form-factors/parametric-bmrvd2022.hh>
#include <eos/form-factors/parametric-dkmr2017.hh>
#include <eos/utils/destringify.hh>

#include <map>

namespace eos
{
    /* J=1/2^+ -> J=1/2^+ Processes */

    /* Lambda_b -> Lambda */

    const constexpr char * LambdaBToLambda::label;
    const constexpr double LambdaBToLambda::tm;
    const constexpr double LambdaBToLambda::tp;
    const constexpr double LambdaBToLambda::mR2_0m;
    const constexpr double LambdaBToLambda::mR2_0p;
    const constexpr double LambdaBToLambda::mR2_1m;
    const constexpr double LambdaBToLambda::mR2_1p;

    /* Lambda_b -> Lambda_c */

    /* Form Factors according to [BFvD2014] */
    class BFvD2014FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            UsedParameter _f_long_v, _b_1_long_v;
            UsedParameter _f_long_a, _b_1_long_a;
            UsedParameter _f_perp_v, _b_1_perp_v;
            UsedParameter _f_perp_a, _b_1_perp_a;

            UsedParameter _m_lambda_b, _m_lambda;

            // Squares of the masses for the vector and axialvector Bbar_s resonances
            static constexpr double mv2 = 5.415 * 5.415;
            static constexpr double ma2 = 5.829 * 5.829;

            static double _z(const double & t, const double & tp, const double & t0)
            {
                return (std::sqrt(tp - t) - std::sqrt(tp - t0)) / (std::sqrt(tp - t) + std::sqrt(tp - t0));
            }

        public:
            BFvD2014FormFactors(const Parameters & p, const Options &) :
                _f_long_v(p["Lambda_b->Lambda::f_0^V(0)@BFvD2014"], *this),
                _b_1_long_v(p["Lambda_b->Lambda::b_1_0^V@BFvD2014"], *this),
                _f_long_a(p["Lambda_b->Lambda::f_0^A(0)@BFvD2014"], *this),
                _b_1_long_a(p["Lambda_b->Lambda::b_1_0^A@BFvD2014"], *this),
                _f_perp_v(p["Lambda_b->Lambda::f_perp^V(0)@BFvD2014"], *this),
                _b_1_perp_v(p["Lambda_b->Lambda::b_1_perp^V@BFvD2014"], *this),
                _f_perp_a(p["Lambda_b->Lambda::f_perp^A(0)@BFvD2014"], *this),
                _b_1_perp_a(p["Lambda_b->Lambda::b_1_perp^A@BFvD2014"], *this),
                _m_lambda_b(p["mass::Lambda_b"], *this),
                _m_lambda(p["mass::Lambda"], *this)
            {
            }

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, const Options & options)
            {
                return new BFvD2014FormFactors(parameters, options);
            }

            virtual double f_long_v(const double & s) const
            {
                const double tp = power_of<2>(_m_lambda_b + _m_lambda);
                const double zt = _z(s, tp, 12.0), z0 = _z(0.0, tp, 12.0);

                return _f_long_v() / (1.0 - s / mv2) * (1.0 + _b_1_long_v() * (zt - z0));
            }

            virtual double f_perp_v(const double & s) const
            {
                const double tp = power_of<2>(_m_lambda_b + _m_lambda);
                const double zt = _z(s, tp, 12.0), z0 = _z(0, tp, 12.0);

                return _f_perp_v() / (1.0 - s / mv2) * (1.0 + _b_1_perp_v() * (zt - z0));
            }

            virtual double f_long_a(const double & s) const
            {
                const double tp = power_of<2>(_m_lambda_b + _m_lambda);
                const double zt = _z(s, tp, 12.0), z0 = _z(0, tp, 12.0);

                return _f_long_a() / (1.0 - s / ma2) * (1.0 + _b_1_long_a() * (zt - z0));
            }

            virtual double f_perp_a(const double & s) const
            {
                const double tp = power_of<2>(_m_lambda_b + _m_lambda);
                const double zt = _z(s, tp, 12.0), z0 = _z(0, tp, 12.0);

                return _f_perp_a() / (1.0 - s / ma2) * (1.0 + _b_1_perp_a() * (zt - z0));
            }

            // Not yet implemented:
            virtual double f_time_v(const double &) const { throw InternalError("BFvD2014FormFactors::f_time_v(): not implemented"); }
            virtual double f_time_a(const double &) const { throw InternalError("BFvD2014FormFactors::f_time_a(): not implemented"); }
            virtual double f_perp_t(const double &) const { throw InternalError("BFvD2014FormFactors::f_perp_t(): not implemented"); }
            virtual double f_perp_t5(const double &) const { throw InternalError("BFvD2014FormFactors::f_perp_t5(): not implemented"); }
            virtual double f_long_t(const double &) const { throw InternalError("BFvD2014FormFactors::f_long_t(): not implemented"); }
            virtual double f_long_t5(const double &) const { throw InternalError("BFvD2014FormFactors::f_long_t5(): not implemented"); }
    };

    FormFactors<OneHalfPlusToOneHalfPlus>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<OneHalfPlusToOneHalfPlus>::KeyType, FormFactorFactory<OneHalfPlusToOneHalfPlus>::ValueType>
    FormFactorFactory<OneHalfPlusToOneHalfPlus>::form_factors
    {
        { "Lambda_b->Lambda::BFvD2014",   &BFvD2014FormFactors::make                      },
        { "Lambda_b->Lambda::DM2016",     &DM2016FormFactors<LambdaBToLambda>::make       },
        { "Lambda_b->Lambda::BMRvD2022",  &BMRvD2022FormFactors<LambdaBToLambda>::make    },
        { "Lambda_b->Lambda_c::DKMR2017", &DKMR2017FormFactors<LambdaBToLambdaC>::make    }
    };

    std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>>
    FormFactorFactory<OneHalfPlusToOneHalfPlus>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> result;

        auto i = FormFactorFactory<OneHalfPlusToOneHalfPlus>::form_factors.find(name);
        if (FormFactorFactory<OneHalfPlusToOneHalfPlus>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<OneHalfPlusToOneHalfPlus>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<OneHalfPlusToOneHalfPlus>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    /* J=1/2^+ -> J=1/2^- Processes */

    /* Lambda_b -> Lambda_c(2595) */

    FormFactors<OneHalfPlusToOneHalfMinus>::~FormFactors()
    {
    }

    Diagnostics
    FormFactors<OneHalfPlusToOneHalfMinus>::diagnostics() const
    {
        return { };
    }

    const std::map<FormFactorFactory<OneHalfPlusToOneHalfMinus>::KeyType, FormFactorFactory<OneHalfPlusToOneHalfMinus>::ValueType>
    FormFactorFactory<OneHalfPlusToOneHalfMinus>::form_factors
    {
        { "Lambda_b->Lambda_c(2595)::HQET",        &HQETFormFactors<OneHalfPlusToOneHalfMinus, LambdaBToLambdaC2595>::make },
    };

    std::shared_ptr<FormFactors<OneHalfPlusToOneHalfMinus>>
    FormFactorFactory<OneHalfPlusToOneHalfMinus>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfMinus>> result;

        auto i = FormFactorFactory<OneHalfPlusToOneHalfMinus>::form_factors.find(name);
        if (FormFactorFactory<OneHalfPlusToOneHalfMinus>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<OneHalfPlusToOneHalfMinus>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<OneHalfPlusToOneHalfMinus>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    /* J=1/2^+ -> J=3/2^- Processes */

    /* Lambda_b -> Lambda_c(2625) */

    /* Lambda_b -> Lambda(1520) */

    const SzegoPolynomial<5> LambdaBToLambda1520::orthonormal_polynomials(SzegoPolynomial<5>::FlatMeasure(3.42519));


    FormFactors<OneHalfPlusToThreeHalfMinus>::~FormFactors()
    {
    }

    Diagnostics
    FormFactors<OneHalfPlusToThreeHalfMinus>::diagnostics() const
    {
        return { };
    }

    const std::map<FormFactorFactory<OneHalfPlusToThreeHalfMinus>::KeyType, FormFactorFactory<OneHalfPlusToThreeHalfMinus>::ValueType>
    FormFactorFactory<OneHalfPlusToThreeHalfMinus>::form_factors
    {
        { "Lambda_b->Lambda_c(2625)::HQET",          &HQETFormFactors<OneHalfPlusToThreeHalfMinus, LambdaBToLambdaC2625>::make },
        { "Lambda_b->Lambda(1520)::ABR2022",         &ABR2022FormFactors<LambdaBToLambda1520>::make },
    };

    std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>>
    FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>> result;

        auto i = FormFactorFactory<OneHalfPlusToThreeHalfMinus>::form_factors.find(name);
        if (FormFactorFactory<OneHalfPlusToThreeHalfMinus>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<OneHalfPlusToThreeHalfMinus>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<OneHalfPlusToThreeHalfMinus>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }
}
