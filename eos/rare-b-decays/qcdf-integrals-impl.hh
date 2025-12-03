/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_QCDF_INTEGRALS_IMPL_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_QCDF_INTEGRALS_IMPL_HH 1

#include <eos/maths/power-of.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    namespace impl
    {
        /* simple inverse moments of the twist-2 LCDAs */
        inline double moment_inverse_ubar(const double & a_1, const double & a_2)
        {
            return 3.0 * (1.0 + a_1 + a_2);
        }

        inline double moment_inverse_ubar2(const double & a1, const double & a2, const double & x)
        {
            return -6.0 * ((1.0 + 3.0 * a1 + 6.0 * a2) * (x + std::log(1.0 - x)) + x * x * (3.0 * a1 + 4.0 * a2 * x));
        }

        /* s > 0, cases for B->V(P)l^+l^- */

        // cf. [vD:2011A], Eq. (26), p. 3
        inline double j0(const double & sh, const double & a1, const double & a2)
        {
            double lnsh = std::log(sh), sh2 = sh * sh, sh3 = sh2 * sh, sh4 = sh2 * sh2;

            // asymptotic part
            double asymp = 3.0 * (1.0 + 2.0 * sh * lnsh - sh2) / power_of<3>(1.0 - sh);
            double gb1 = 3.0 * (1.0 + 9.0 * sh - 9.0 * sh2 - sh3 + 6.0 * sh * (1.0 + sh) * lnsh) / power_of<4>(1.0 - sh);
            double gb2 = 3.0 * (1.0 + 28.0 * sh - 28.0 * sh3 - sh4 + 12.0 * sh * (1.0 + 3.0 * sh + sh2) * lnsh) / power_of<5>(1.0 - sh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline double j0bar(const double & sh, const double & a1, const double & a2)
        {
            return j0(sh, -a1, a2);
        }
    }

}

#endif
