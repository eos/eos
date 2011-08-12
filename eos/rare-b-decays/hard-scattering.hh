/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_HARD_SCATTERING_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_HARD_SCATTERING_HH 1

#include <complex>

namespace eos
{
    using std::complex;

    struct HardScattering
    {
        // cf. [BFS2001], Eqs. (30)-(32), p. 8
        // @s    : Dilepton invariant mass
        // @u    : Relative contribution of the quark (versus ubar = 1-u for the
        //         antiquark) to the light meson's energy.
        // @m_q  : Mass of the internal loop quark
        // @m_B  : Mass of the B meson
        static complex<double> I1(const double & s, const double & u, const double & m_q, const double & m_B);
    };
}

#endif
