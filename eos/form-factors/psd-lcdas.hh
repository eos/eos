/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2022 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PSD_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PSD_LCDAS_HH 1

#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class PseudoscalarLCDAs
    {
        public:
            virtual ~PseudoscalarLCDAs() = 0;

            /* Twist 2 LCDA parameters: Gegenbauer coefficients */
            virtual double a1(const double & mu) const = 0;
            virtual double a2(const double & mu) const = 0;
            virtual double a3(const double & mu) const = 0;
            virtual double a4(const double & mu) const = 0;

            /* Twist 3 LCDA parameters */
            virtual double mu3(const double & mu) const = 0;
            virtual double f3(const double & mu) const = 0;
            virtual double eta3(const double & mu) const = 0;
            virtual double lambda3(const double & mu) const = 0;
            virtual double omega3(const double & mu) const = 0;

            /* Twist 4 LCDA parameter */
            virtual double delta4(const double & mu) const = 0;
            virtual double kappa4(const double & mu) const = 0;
            virtual double omega4(const double & mu) const = 0;

            /* Twist 2 LCDA */
            virtual double phi(const double & u, const double & mu) const = 0;

            /* Twist 3 LCDAs and their derivatives */
            virtual double phi3p(const double & u, const double & mu) const = 0;
            virtual double phi3s(const double & u, const double & mu) const = 0;
            virtual double phi3s_d1(const double & u, const double & mu) const = 0;

            /* Twist 4 LCDAs, their derivatives and integrals */
            virtual double phi4(const double & u, const double & mu) const = 0;
            virtual double phi4_d1(const double & u, const double & mu) const = 0;
            virtual double phi4_d2(const double & u, const double & mu) const = 0;
            virtual double psi4(const double & u, const double & mu) const = 0;
            virtual double psi4_i(const double & u, const double & mu) const = 0;

            static std::shared_ptr<PseudoscalarLCDAs> make(const std::string & name, const Parameters & parameters, const Options & options);
    };
}

#endif
