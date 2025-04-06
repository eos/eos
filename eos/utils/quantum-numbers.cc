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

#include <eos/utils/join.hh>
#include <eos/utils/quantum-numbers.hh>

#include <array>
#include <vector>

namespace eos
{
    std::ostream &
    operator<< (std::ostream & os, LeptonFlavor lf)
    {
        static const std::array<std::string, 3u> names
        {
            "e", "mu", "tau"
        };

        os << names[static_cast<unsigned>(lf)];

        return os;
    }

    std::ostream &
    operator<< (std::ostream & os, QuarkFlavor qf)
    {
        static const std::array<std::string, 6u> names
        {
            "u", "d", "s", "c", "b", "t"
        };

        os << names[static_cast<unsigned>(qf)];

        return os;
    }

    std::ostream &
    operator<< (std::ostream & os, Isospin i)
    {
        static const std::array<std::string, 5u> names
        {
            "0", "1/2", "1", "3/2", "2"
        };

        std::vector<std::string> tmp;
        for (unsigned shift = 0 ; shift < 5 ; ++shift)
        {
            if (i && static_cast<Isospin>(1 << shift))
                tmp.push_back(names[shift]);
        }

        os << join(tmp.begin(), tmp.end(), "|");

        return os;
    }

    std::ostream &
    operator<< (std::ostream & os, IsospinRepresentation ir)
    {
        static const std::array<std::string, 5u> names
        {
            "0", "1", "2", "1/2", "3/2"
        };

        os << names[static_cast<unsigned>(ir)];

        return os;
    }

    std::ostream &
    operator<< (std::ostream & os, LightMeson qf)
    {
        static const std::array<std::string, 10u> names
        {
            "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime"
        };

        os << names[static_cast<unsigned>(qf)];

        return os;
    }

    std::ostream &
    operator<< (std::ostream & os, PartialWave pw)
    {
        static const std::array<std::string, 4u> names
        {
            "S", "P", "D", "F"
        };

        std::vector<std::string> tmp;
        for (unsigned shift = 0 ; shift < 4 ; ++shift)
        {
            if (pw && static_cast<PartialWave>(1 << shift))
                tmp.push_back(names[shift]);
        }

        os << join(tmp.begin(), tmp.end(), "|");

        return os;
    }

}
