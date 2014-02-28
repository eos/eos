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

#include <eos/form-factors/baryonic-impl.hh>
#include <eos/utils/destringify.hh>

#include <map>

namespace eos
{
    /* J=1/2^+ -> J=1/2^+ Processes */

    /* Lambda_b -> Lambda */

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

            static constexpr double _z(const double & t, const double & tp, const double & t0)
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

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, unsigned)
            {
                return new BFvD2014FormFactors(parameters, Options());
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
    };

    FormFactors<OneHalfPlusToOneHalfPlus>::~FormFactors()
    {
    }

    std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>>
    FormFactorFactory<OneHalfPlusToOneHalfPlus>::create(const std::string & label, const Parameters & parameters)
    {
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> result;

        typedef std::tuple<std::string, std::string> KeyType;
        typedef std::function<FormFactors<OneHalfPlusToOneHalfPlus> * (const Parameters &, unsigned)> ValueType;
        static const std::map<KeyType, ValueType> form_factors
        {
            { KeyType("Lambda_b->Lambda",     "BFvD2014"),   &BFvD2014FormFactors::make     },
        };

        /*
         * Labels have the form
         *
         *   PROCESS@NAME[:SET]
         *
         * The brackets indicate the latter part to be optional.
         */

        std::string process, name, input(label);
        unsigned set(0);

        std::string::size_type sep_at(input.find('@')), sep_colon(input.find(':'));
        if (std::string::npos == sep_at)
            return result;

        if (std::string::npos != sep_colon)
        {
            set = destringify<unsigned>(input.substr(sep_colon + 1));
            input.erase(sep_colon + 1);
        }

        name = input.substr(sep_at + 1);
        process = input.substr(0, sep_at);

        auto i = form_factors.find(KeyType(process, name));
        if (form_factors.cend() == i)
            return result;

        result = std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>>(i->second(parameters, set));

        return result;
    }
}
