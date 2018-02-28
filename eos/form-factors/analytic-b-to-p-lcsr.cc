/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
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

#include <eos/form-factors/analytic-b-to-p-lcsr-impl.hh>

namespace eos
{
    namespace lcsr
    {
        // B -> pi
        constexpr const char * BToPi::P;
        constexpr const char * BToPi::m_P;
        constexpr const char * BToPi::f_P;

        // B -> K
        constexpr const char * BToK::P;
        constexpr const char * BToK::m_P;
        constexpr const char * BToK::f_P;

        // B -> D
        constexpr const char * BToD::P;
        constexpr const char * BToD::m_P;
        constexpr const char * BToD::f_P;
    }

    template class AnalyticFormFactorBToPLCSR<lcsr::BToPi>;
    template class AnalyticFormFactorBToPLCSR<lcsr::BToK>;
    template class AnalyticFormFactorBToPLCSR<lcsr::BToD>;
}
