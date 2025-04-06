/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2025 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_QUANTUM_NUMBERS_HH
#define EOS_GUARD_EOS_UTILS_QUANTUM_NUMBERS_HH 1

#include <eos/utils/stringify.hh>

#include <cstdint>
#include <ostream>
#include <string>
#include <type_traits>

namespace eos
{
    enum class LeptonFlavor : uint8_t
    {
        electron = 0,
        muon     = 1,
        tauon    = 2
    };

    std::ostream& operator<< (std::ostream & os, LeptonFlavor lf);

    enum class QuarkFlavor : uint8_t
    {
        up      = 0,
        down    = 1,
        strange = 2,
        charm   = 3,
        bottom  = 4,
        top     = 5
    };
    std::ostream& operator<< (std::ostream & os, QuarkFlavor qf);

    enum class Isospin : uint16_t
    {
        none        = 0,
        zero        = 1 << 0,
        onehalf     = 1 << 1,
        one         = 1 << 2,
        threehalves = 1 << 3,
        two         = 1 << 4
    };

    inline Isospin operator | (Isospin lhs, Isospin rhs)
    {
        using T = std::underlying_type_t<Isospin>;
        return static_cast<Isospin>(static_cast<T>(lhs) | static_cast<T>(rhs));
    }

    inline Isospin & operator |= (Isospin & lhs, Isospin rhs)
    {
        lhs = lhs | rhs;
        return lhs;
    }

    inline Isospin operator & (Isospin lhs, Isospin rhs)
    {
        using T = std::underlying_type_t<Isospin>;
        return static_cast<Isospin>(static_cast<T>(lhs) & static_cast<T>(rhs));
    }

    inline Isospin & operator &= (Isospin & lhs, Isospin rhs)
    {
        lhs = lhs & rhs;
        return lhs;
    }

    inline bool operator && (Isospin lhs, Isospin rhs)
    {
        return (lhs & rhs) != Isospin::none;
    }

    std::ostream & operator<< (std::ostream & os, Isospin i);

    enum class IsospinRepresentation : uint8_t
    {
        zero        = 0,
        one         = 1,
        two         = 2,
        onehalf     = 3,
        threehalves = 4
    };

    std::ostream& operator<< (std::ostream & os, IsospinRepresentation ir);

    enum class LightMeson : uint8_t
    {
        pi0      = 0,
        piplus   = 1,
        piminus  = 2,
        K0       = 3,
        K0bar    = 4,
        KS       = 5,
        Kplus    = 6,
        Kminus   = 7,
        eta      = 8,
        etap     = 9,
    };

    std::ostream& operator<< (std::ostream & os, LightMeson qf);
}

#endif
