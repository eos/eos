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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PI_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PI_LCDAS_HH 1

#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class PionLCDAs :
        public ParameterUser,
        public PrivateImplementationPattern<PionLCDAs>
    {
        public:
            PionLCDAs(const Parameters &, const Options &);
            ~PionLCDAs();

            /* Twist 2 LCDA (even) Gegenbauer coefficients */
            double a2pi(const double & mu) const;
            double a4pi(const double & mu) const;

            /* Twist 3 LCDA parameters */
            double mupi(const double & mu) const;
            double f3pi(const double & mu) const;
            double eta3pi(const double & mu) const;
            double omega3pi(const double & mu) const;

            /* Twist 4 LCDA parameter */
            double deltapipi(const double & mu) const;
            double omega4pi(const double & mu) const;

            /* Twist 2 LCDA */
            double phi(const double & u, const double & mu) const;

            /* Twist 3 LCDAs and their derivatives */
            double phi3p(const double & u, const double & mu) const;
            double phi3s(const double & u, const double & mu) const;
            double phi3s_d1(const double & u, const double & mu) const;

            /* Twist 4 LCDAs, their derivatives and integrals */
            double phi4(const double & u, const double & mu) const;
            double phi4_d1(const double & u, const double & mu) const;
            double phi4_d2(const double & u, const double & mu) const;
            double psi4(const double & u, const double & mu) const;
            double psi4_i(const double & u, const double & mu) const;

            /* Internal diagnostics */
            Diagnostics diagnostics() const;
    };
}

#endif
