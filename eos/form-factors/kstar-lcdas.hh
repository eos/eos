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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_KSTAR_LCDAS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_KSTAR_LCDAS_HH 1

#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    class KstarLCDAs :
        public ParameterUser,
        public PrivateImplementationPattern<KstarLCDAs>
    {
        public:
            KstarLCDAs(const Parameters &, const Options &);
            ~KstarLCDAs();

            /* Twist 2 LCDA for the vector current: Gegenbauer coefficients */
            double a_1_para(const double & mu) const;
            double a_2_para(const double & mu) const;

            /* Twist 2 LCDA for the vector current: Normalization */
            double f_para() const;

            /* Twist 2 LCDA for the tensor current: Gegenbauer coefficients */
            double a_1_perp(const double & mu) const;
            double a_2_perp(const double & mu) const;

            /* Twist 2 LCDA for the tensor current: Normalization */
            double f_perp(const double & mu) const;

            /* Twist 2 LCDAs */
            double phi_2_perp(const double & u, const double & mu) const;
            double phi_2_para(const double & u, const double & mu) const;

            /* Twist 2 LCDAs : partial inverse moment */
            double part_inv_mom_phi_2_para(const double & u, const double & mu) const;

            /* Internal diagnostics */
            Diagnostics diagnostics() const;
    };
}

#endif
