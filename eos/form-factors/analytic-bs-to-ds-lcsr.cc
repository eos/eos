/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2025 Danny van Dyk
 * Copyright (c) 2019      Nico Gubernari
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
    template <>
    struct AnalyticFormFactorBToPLCSRProcessTraits<BsToDs>
    {
        static constexpr const char * label  = "B_s->D_s";
        static constexpr const char * name_B = "mass::B_s";
        static constexpr const char * f_B    = "decay-constant::B_s";
        static constexpr const char * name_P = "mass::D_s";
        static constexpr const char * f_P    = "decay-constant::D_s";
        static constexpr const std::tuple<QuarkFlavor, QuarkFlavor> partonic_transition = std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm);
        static constexpr const QuarkFlavor spectator_flavor                             = QuarkFlavor::strange;
        static constexpr const double chi2 = 1.0;
    };

    template class AnalyticFormFactorBToPLCSR<BsToDs>;
}
