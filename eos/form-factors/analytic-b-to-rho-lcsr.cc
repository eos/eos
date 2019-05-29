/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
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

#include <eos/form-factors/analytic-b-to-v-lcsr-impl.hh>

namespace eos
{
    namespace lcsr
    {
        struct BToRho
        {
            constexpr static const char * B    = "B";
            constexpr static const char * m_B  = "mass::B_d";
            constexpr static const char * f_B  = "decay-constant::B_d";
            constexpr static const char * V    = "rho";
            constexpr static const char * m_V  = "mass::rho^+";
            constexpr static const char * f_V  = "decay-constant::rho";
            constexpr static const char   q_v  = 'u';
            constexpr static const char   q_s  = 'd';
            constexpr static const double chi2 = 1.0;
        };

        // B -> rho
        constexpr const char * BToRho::B;
        constexpr const char * BToRho::m_B;
        constexpr const char * BToRho::f_B;
        constexpr const char * BToRho::V;
        constexpr const char * BToRho::m_V;
        constexpr const char * BToRho::f_V;
        constexpr const char   BToRho::q_s;
    }

    template class AnalyticFormFactorBToVLCSR<lcsr::BToRho>;
}
