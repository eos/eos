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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_VEC_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_VEC_LCDAS_HH 1

#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class VectorLCDAs :
        public ParameterUser
    {
        public:
            virtual ~VectorLCDAs() = 0;

            /* Twist 2 LCDA parameters: para Gegenbauer coefficients */
            virtual double a1para(const double & mu) const = 0;
            virtual double a2para(const double & mu) const = 0;
            virtual double a3para(const double & mu) const = 0;
            virtual double a4para(const double & mu) const = 0;
            virtual double fpara() const = 0;

            /* Twist 2 LCDA parameters: perp Gegenbauer coefficients */
            virtual double a1perp(const double & mu) const = 0;
            virtual double a2perp(const double & mu) const = 0;
            virtual double a3perp(const double & mu) const = 0;
            virtual double a4perp(const double & mu) const = 0;
            virtual double fperp(const double & mu) const = 0;

            /* Twist 2 LCDAs */
            virtual double phipara(const double & u, const double & mu) const = 0;
            virtual double phiperp(const double & u, const double & mu) const = 0;

            /* Twist 3 parameters */
            virtual double kappa3para(const double & mu) const = 0;
            virtual double omega3para(const double & mu) const = 0;
            virtual double lambda3para(const double & mu) const = 0;
            virtual double zeta3para(const double & mu) const = 0;
            virtual double lambda3paratilde(const double & mu) const = 0;
            virtual double omega3paratilde(const double & mu) const = 0;
            virtual double kappa3perp(const double & mu) const = 0;
            virtual double omega3perp(const double & mu) const = 0;
            virtual double lambda3perp(const double & mu) const = 0;


            static std::shared_ptr<VectorLCDAs> make(const std::string & name, const Parameters & parameters, const Options & options);
    };
}

#endif
