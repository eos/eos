/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <map>

namespace eos
{
    template <>
    struct Implementation<Options>
    {
        std::map<std::string, std::string> options;

        Implementation()
        {
        }

        Implementation(const std::initializer_list<std::pair<std::string, std::string>> & _options)
        {
            for (auto o = _options.begin(), o_end = _options.end() ; o != o_end ; ++o)
            {
                options.insert(*o);
            }
        }
    };

    Options::Options() :
        PrivateImplementationPattern<Options>(new Implementation<Options>)
    {
    }

    Options::Options(const std::initializer_list<std::pair<std::string, std::string>> & options) :
        PrivateImplementationPattern<Options>(new Implementation<Options>(options))
    {
    }

    Options::~Options()
    {
    }

    bool
    Options::operator== (const Options & rhs) const
    {
        if (_imp->options.size() != rhs._imp->options.size())
            return false;

        for (auto l = _imp->options.cbegin(), l_end = _imp->options.cend(), r = rhs._imp->options.cbegin() ; l != l_end ; ++l, ++r)
        {
            if (*l != *r)
                return false;
        }

        return true;
    }

    bool
    Options::operator!= (const Options & rhs) const
    {
        return ! (*this == rhs);
    }

    const std::string &
    Options::operator[] (const std::string & key) const
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() == i)
            throw UnknownOptionError(key);

        return i->second;
    }

    bool
    Options::has(const std::string & key) const
    {
        return _imp->options.end() != _imp->options.find(key);
    }

    void
    Options::set(const std::string & key, const std::string & value)
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
    Options::get(const std::string & key, const std::string & default_value) const
    {
        auto i(_imp->options.find(key));

        if (_imp->options.end() == i)
            return default_value;

        return i->second;
    }

    std::string
    Options::as_string() const
    {
        std::string result;

        auto i(_imp->options.cbegin()), i_end(_imp->options.cend());

        if (i != i_end)
        {
            result += i->first + '=' + i->second;
            ++i;
        }

        for ( ; i != i_end ; ++i)
        {
            result += ',' + i->first + '=' + i->second;
        }

        return result;
    }

    UnknownOptionError::UnknownOptionError(const std::string & key) throw () :
        Exception("Unknown option: '" + key + "'")
    {
    }

    Options
    operator+ (const Options & lhs, const Options & rhs)
    {
        Options result;

        result._imp->options.insert(lhs._imp->options.cbegin(), lhs._imp->options.cend());
        result._imp->options.insert(rhs._imp->options.cbegin(), rhs._imp->options.cend());

        return result;
    }
}
