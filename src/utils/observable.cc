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

#include <src/utils/observable.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <map>

namespace eos
{
    Observable::Observable()
    {
    }

    Observable::~Observable()
    {
    }

    template <>
    struct Implementation<ObservableOptions>
    {
        std::map<std::string, std::string> options;
    };

    ObservableOptions::ObservableOptions() :
        PrivateImplementationPattern<ObservableOptions>(new Implementation<ObservableOptions>)
    {
    }

    ObservableOptions::~ObservableOptions()
    {
    }

    const std::string &
    ObservableOptions::operator[] (const std::string & key) const
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() == i)
            throw UnknownOptionError(key);

        return i->second;
    }

    bool
    ObservableOptions::has(const std::string & key) const
    {
        return _imp->options.end() != _imp->options.find(key);
    }

    void
    ObservableOptions::set(const std::string & key, const std::string & value)
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() != i)
        {
            i->second = value;
        }
        else
        {
            _imp->options[key] = value;
        }
    }

    std::string
    ObservableOptions::get(const std::string & key, const std::string & default_value) const
    {
        auto i(_imp->options.find(key));

        if (_imp->options.end() == i)
            return default_value;

        return i->second;
    }

    ObservableFactory::ObservableFactory()
    {
    }

    ObservableFactory::~ObservableFactory()
    {
    }

    UnknownOptionError::UnknownOptionError(const std::string & key) throw () :
        Exception("Unknown option: '" + key + "'")
    {
    }
}
