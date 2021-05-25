
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/units.hh>

#include <string>

namespace eos
{
    Unit::Unit(const std::string & latex) :
        _latex(latex)
    {
    }

    const Unit &
    Unit::Undefined()
    {
        static const Unit undefined(R"(\textrm{undefined})");

        return undefined;
    }

    const Unit &
    Unit::None()
    {
        static const Unit none("1");

        return none;
    }

    const Unit &
    Unit::GeV()
    {
        static const Unit gev(R"(\textrm{GeV})");

        return gev;
    }

    const Unit &
    Unit::GeV2()
    {
        static const Unit gev_2(R"(\textrm{GeV}^2)");

        return gev_2;
    }

    const Unit &
    Unit::InverseGeV2()
    {
        static const Unit inverse_gev_2(R"(\textrm{GeV}^{-2})");

        return inverse_gev_2;
    }

    const Unit &
    Unit::InverseGeV4()
    {
        static const Unit inverse_gev_4(R"(\textrm{GeV}^{-4})");

        return inverse_gev_4;
    }

    const Unit &
    Unit::InversePicoSecond()
    {
        static const Unit inverse_pico_second(R"(\textrm{ps}^{-1})");

        return inverse_pico_second;
    }
}