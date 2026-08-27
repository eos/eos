/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2025 Danny van Dyk
 * Copyright (c) 2018      Ahmet Kokulu
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
    template <> struct AnalyticFormFactorBToVLCSRTraits<BToKstar>
    {
            static const constexpr char *                               label               = "B->K^*";
            static const constexpr char *                               name_B              = "mass::B_d";
            static const constexpr char *                               f_B                 = "decay-constant::B_d";
            static const constexpr char *                               name_V              = "mass::K_d^*";
            static const constexpr char *                               f_V                 = "B->K^*::f_Kstar_par";
            static const constexpr std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange);
            static const constexpr QuarkFlavor                          spectator_flavor    = QuarkFlavor::down;
            static const constexpr double                               chi2                = 1.0;
    };

    template class AnalyticFormFactorBToVLCSR<BToKstar>;
} // namespace eos
