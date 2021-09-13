
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
#include <vector>

namespace eos
{
    const std::string &
    Unit::latex() const
    {
        static const std::vector<std::string> representations
        {
            R"(\textrm{undefined})",
            "1",
            R"(\textrm{GeV})",
            R"(\textrm{GeV}^2)",
            R"(\textrm{GeV}^{-2})",
            R"(\textrm{GeV}^{-4})",
            R"(\textrm{ps}^{-1})",
            R"(\textrm{undefined})"
        };

        return representations[static_cast<int>(_id)];
    }

    Unit
    Unit::Undefined()
    {
        return Unit(Id::undefined);
    }

    Unit
    Unit::None()
    {
        return Unit(Id::none);
    }

    Unit
    Unit::GeV()
    {
        return Unit(Id::gev);
    }

    Unit
    Unit::GeV2()
    {
        return Unit(Id::gev2);
    }

    Unit
    Unit::InverseGeV2()
    {
        return Unit(Id::inverse_gev2);
    }

    Unit
    Unit::InverseGeV4()
    {
        return Unit(Id::inverse_gev4);
    }

    Unit
    Unit::InversePicoSecond()
    {
        return Unit(Id::inverse_ps);
    }

    bool
    Unit::operator== (const Unit & rhs) const
    {
        return this->_id == rhs._id;
    }
}