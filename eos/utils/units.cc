
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

#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/units.hh>

#include <map>
#include <string>
#include <vector>

namespace eos
{
    const std::vector<std::string> Unit::_latex_representations{
        R"(\textrm{undefined})", "1",
        R"(\textrm{GeV})",       R"(\textrm{GeV}^2)",
        R"(\textrm{GeV}^3)",     R"(\textrm{GeV}^{-1})",
        R"(\textrm{GeV}^{-2})",  R"(\textrm{GeV}^{-4})",
        R"(\textrm{s})",         R"(\textrm{s}^{-1})",
        R"(\textrm{ps}^{-1})",   R"(\textrm{GeV}\,\textrm{s})",
        R"(\textrm{fm}^2)",
    };

    const std::vector<std::string> Unit::_internal_representations{
        "undefined", "1", "GeV", "GeV^2", "GeV^3", "GeV^-1", "GeV^-2", "GeV^-4", "s", "s^-1", "ps^-1", "GeV s", "fm^2",
    };

    const std::string &
    Unit::latex() const
    {
        return _latex_representations[static_cast<int>(_id)];
    }

    const std::string &
    Unit::string() const
    {
        return _internal_representations[static_cast<int>(_id)];
    }

    Unit::Unit(const std::string & s) :
        Unit(Id::undefined)
    {
        static const std::map<std::string, Id> map{
            {    _internal_representations[static_cast<unsigned>(Id::undefined)],    Id::undefined },
            {         _internal_representations[static_cast<unsigned>(Id::none)],         Id::none },
            {          _internal_representations[static_cast<unsigned>(Id::gev)],          Id::gev },
            {         _internal_representations[static_cast<unsigned>(Id::gev2)],         Id::gev2 },
            {         _internal_representations[static_cast<unsigned>(Id::gev3)],         Id::gev3 },
            {  _internal_representations[static_cast<unsigned>(Id::inverse_gev)],  Id::inverse_gev },
            { _internal_representations[static_cast<unsigned>(Id::inverse_gev2)], Id::inverse_gev2 },
            { _internal_representations[static_cast<unsigned>(Id::inverse_gev4)], Id::inverse_gev4 },
            {            _internal_representations[static_cast<unsigned>(Id::s)],            Id::s },
            {    _internal_representations[static_cast<unsigned>(Id::inverse_s)],    Id::inverse_s },
            {   _internal_representations[static_cast<unsigned>(Id::inverse_ps)],   Id::inverse_ps },
            {        _internal_representations[static_cast<unsigned>(Id::gev_s)],        Id::gev_s },
            {          _internal_representations[static_cast<unsigned>(Id::fm2)],          Id::fm2 },
        };

        const auto i = map.find(s);
        if (map.cend() == i)
        {
            Log::instance()->message("Unit", ll_error) << "Unrecognized unit '" << s << "' encountered";
        }
        else
        {
            _id = i->second;
        }
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
    Unit::GeV3()
    {
        return Unit(Id::gev3);
    }

    Unit
    Unit::InverseGeV()
    {
        return Unit(Id::inverse_gev);
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
    Unit::Second()
    {
        return Unit(Id::s);
    }

    Unit
    Unit::InverseSecond()
    {
        return Unit(Id::inverse_s);
    }

    Unit
    Unit::InversePicoSecond()
    {
        return Unit(Id::inverse_ps);
    }

    Unit
    Unit::GeVSecond()
    {
        return Unit(Id::gev_s);
    }

    Unit
    Unit::Femtometer2()
    {
        return Unit(Id::fm2);
    }

    bool
    Unit::operator== (const Unit & rhs) const
    {
        return this->_id == rhs._id;
    }
} // namespace eos
