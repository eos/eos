/*
 * Copyright (c) 2014, 2015, 2016, 2017 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DM2016_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DM2016_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/models/model.hh>
#include <eos/utils/reference-name.hh>

#include <set>
#include <vector>

namespace eos
{
    template <typename Process_>
    class DM2016FormFactorTraits;

    template <typename Process_>
    class DM2016FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            UsedParameter _alpha_0_time_v, _alpha_1_time_v, _alpha_2_time_v;
            UsedParameter _alpha_0_time_a, _alpha_1_time_a, _alpha_2_time_a;

            UsedParameter _alpha_0_long_v, _alpha_1_long_v, _alpha_2_long_v;
            UsedParameter _alpha_0_long_a, _alpha_1_long_a, _alpha_2_long_a;
            UsedParameter _alpha_0_perp_v, _alpha_1_perp_v, _alpha_2_perp_v;
            UsedParameter                  _alpha_1_perp_a, _alpha_2_perp_a;

            UsedParameter _alpha_0_long_t,  _alpha_1_long_t,  _alpha_2_long_t;
            UsedParameter _alpha_0_long_t5, _alpha_1_long_t5, _alpha_2_long_t5;
            UsedParameter _alpha_0_perp_t,  _alpha_1_perp_t,  _alpha_2_perp_t;
            UsedParameter                   _alpha_1_perp_t5, _alpha_2_perp_t5;

            std::unique_ptr<const DM2016FormFactorTraits<Process_>> _traits;

        public:
            DM2016FormFactors(const Parameters & parameters, const Options & options);
            virtual ~DM2016FormFactors() = default;

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, const Options & options);

            // vector current
            virtual double f_time_v(const double & s) const;
            virtual double f_long_v(const double & s) const;
            virtual double f_perp_v(const double & s) const;
            // axial vector current
            virtual double f_time_a(const double & s) const;
            virtual double f_long_a(const double & s) const;
            virtual double f_perp_a(const double & s) const;
            // tensor current
            virtual double f_long_t(const double & s) const;
            virtual double f_perp_t(const double & s) const;
            // axial tensor current
            virtual double f_long_t5(const double & s) const;
            virtual double f_perp_t5(const double & s) const;

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

    extern template class DM2016FormFactors<LambdaBToLambda>;
}

#endif
