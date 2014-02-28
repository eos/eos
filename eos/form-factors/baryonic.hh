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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_HH
#define EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/utils/parameters.hh>

#include <memory>
#include <string>

namespace eos
{
    /* Baryonic Tags */

    /*
     * J=1/2^+ -> J=1/2^+ transitions
     */
    struct OneHalfPlusToOneHalfPlus { };

    template <>
    class FormFactors<OneHalfPlusToOneHalfPlus> :
        public ParameterUser
    {
        public:
            virtual ~FormFactors();

            //virtual double f_time_v(const double & s) const = 0;
            virtual double f_long_v(const double & s) const = 0;
            virtual double f_perp_v(const double & s) const = 0;

            //virtual double f_time_a(const double & s) const = 0;
            virtual double f_long_a(const double & s) const = 0;
            virtual double f_perp_a(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<OneHalfPlusToOneHalfPlus>
    {
        public:
            static std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> create(const std::string & label, const Parameters & parameters);
    };
}

#endif
