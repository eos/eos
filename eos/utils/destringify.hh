/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_DESTRINGIFY_HH
#define EOS_GUARD_EOS_UTILS_DESTRINGIFY_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/quantum-numbers.hh>

#include <map>
#include <ranges>
#include <sstream>
#include <string_view>

namespace eos
{
    /**
     * Exception thrown when destringification fails.
     */
    class DestringifyError :
        public Exception
    {
        public:
            DestringifyError(const std::string & str);
    };

    namespace impl
    {
        template <typename T_>
        struct SimpleDestringify
        {
            static T_ destringify(const std::string & input)
            {
                std::stringstream ss(input);
                T_ value;

                ss >> value;
                if (! ss.eof() || ss.fail())
                {
                    throw DestringifyError(input);
                }

                return value;
            }
        };

        template <typename T_> struct DoDestringify : public SimpleDestringify<T_> {};

        template <> struct DoDestringify<bool>
        {
            static bool destringify(const std::string & input)
            {
                if ("true" == input)
                    return true;

                return false;
            }
        };

        template <> struct DoDestringify<Isospin>
        {
            static Isospin destringify(const std::string & input)
            {
                static const std::map<std::string, Isospin> isospins
                {
                    { "0",   Isospin::zero        },
                    { "1",   Isospin::one         },
                    { "1/2", Isospin::onehalf     },
                    { "2",   Isospin::two         },
                    { "3/2", Isospin::threehalves }
                };

                Isospin result = Isospin::none;

                std::string::size_type i = 0, j = input.find('|');
                do
                {
                    const auto value = input.substr(i, j);

                    const auto k = isospins.find(value);
                    if (isospins.cend() == k)
                        throw DestringifyError(std::string("'") + value + "' is not a valid Isospin value");

                    result |= k->second;

                    if (std::string::npos != j)
                    {
                        i = j + 1;
                        j = input.find('|', i);
                    }
                    else
                    {
                        i = std::string::npos;
                    }
                }
                while (std::string::npos != i);

                return result;
            }
        };
    }

    template <typename T_>
    T_
    destringify(const std::string & input)
    {
        return impl::DoDestringify<T_>::destringify(input);
    }
}

#endif
