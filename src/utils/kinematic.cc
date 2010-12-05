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

#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/stringify.hh>

#include <map>

namespace eos
{
    template <>
    struct Implementation<Kinematics>
    {
        std::map<std::string, double> variables;
    };

    Kinematics::Kinematics() :
        PrivateImplementationPattern<Kinematics>(new Implementation<Kinematics>)
    {
    }

    Kinematics::~Kinematics()
    {
    }

    double
    Kinematics::operator[] (const std::string & variable) const
    {
        auto i(_imp->variables.find(variable));

        if (_imp->variables.end() == i)
            throw UnknownKinematicVariableError(variable);

        return i->second;
    }

    std::string
    Kinematics::as_string () const
    {
        std::string result;

        auto i(_imp->variables.cbegin()), i_end(_imp->variables.cend());
        if (i != i_end)
        {
            result = i->first + '=' + stringify(i->second);
            ++i;
        }


        for ( ; i != i_end ; ++i)
        {
            result += ", " + i->first + '=' + stringify(i->second);
        }

        return result;
    }

    void
    Kinematics::declare(const std::string & variable)
    {
        auto i(_imp->variables.find(variable));

        if (_imp->variables.end() == i)
            _imp->variables[variable] = 0.0;
    }

    void
    Kinematics::set(const std::string & variable, const double & value)
    {
        auto i(_imp->variables.find(variable));

        if (_imp->variables.end() == i)
            throw UnknownKinematicVariableError(variable);

        i->second = value;
    }

    UnknownKinematicVariableError::UnknownKinematicVariableError(const std::string & variable) throw () :
        Exception("Unknown kinematic variable: '" + variable + "'")
    {
    }
}
