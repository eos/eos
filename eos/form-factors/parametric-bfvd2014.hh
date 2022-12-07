/*
 * Copyright (c) 2014-2017 Danny van Dyk
 * Copyright (c) 2022 Méril Reboud
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BFvD2014_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BFvD2014_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/models/model.hh>
#include <eos/utils/reference-name.hh>

#include <set>
#include <vector>

namespace eos
{
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
            BFvD2014FormFactors(const Parameters & parameters, const Options & options);
            virtual ~BFvD2014FormFactors() = default;

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, const Options & options);

            virtual double f_long_v(const double & s) const;
            virtual double f_perp_v(const double & s) const;
            virtual double f_long_a(const double & s) const;
            virtual double f_perp_a(const double & s) const;

            // Not yet implemented:
            virtual double f_time_v(const double &) const { throw InternalError("BFvD2014FormFactors::f_time_v(): not implemented"); }
            virtual double f_time_a(const double &) const { throw InternalError("BFvD2014FormFactors::f_time_a(): not implemented"); }
            virtual double f_perp_t(const double &) const { throw InternalError("BFvD2014FormFactors::f_perp_t(): not implemented"); }
            virtual double f_perp_t5(const double &) const { throw InternalError("BFvD2014FormFactors::f_perp_t5(): not implemented"); }
            virtual double f_long_t(const double &) const { throw InternalError("BFvD2014FormFactors::f_long_t(): not implemented"); }
            virtual double f_long_t5(const double &) const { throw InternalError("BFvD2014FormFactors::f_long_t5(): not implemented"); }

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> options;
    };
}

#endif
