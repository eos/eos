/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_HH 1

#include <src/utils/parameters.hh>

#include <memory>
#include <string>

namespace eos
{
    template <typename Transition_>
    class FormFactors;

    template <typename Transition_>
    class FormFactorFactory;

    /* Tags */

    /*
     * P -> V transitions
     *
     * P: (Heavy) pseudoscalar meson.
     * V: Light vector meson.
     */
    struct PToV { };

    struct PToP { };

    template <>
    class FormFactors<PToV>
    {
        public:
            virtual ~FormFactors();

            virtual double v(const double & s) const = 0;

            virtual double a_0(const double & s) const = 0;
            virtual double a_1(const double & s) const = 0;
            virtual double a_2(const double & s) const = 0;

            // TODO: dipole form factors
    };

    template <>
    class FormFactorFactory<PToV>
    {
        public:
            static std::shared_ptr<FormFactors<PToV>> create(const std::string & label, const Parameters & parameters);
    };

    template <>
    class FormFactors<PToP>
    {
        public:
            virtual ~FormFactors();

            virtual double f_p(const double & s) const = 0;
            virtual double f_0(const double & s) const = 0;
            virtual double f_t(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<PToP>
    {
        public:
            static std::shared_ptr<FormFactors<PToP>> create(const std::string & label, const Parameters & parameters);
    };
}


#endif
