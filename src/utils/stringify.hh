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

#ifndef EOS_GUARD_SRC_UTILS_STRINGIFY_HH
#define EOS_GUARD_SRC_UTILS_STRINGIFY_HH 1

#include <string>
#include <sstream>

namespace eos
{
    namespace implementation
    {
        template <typename T_>
        struct DoStringify
        {
            static std::string stringify(const T_ & x)
            {
                std::stringstream ss;
                ss.precision(10);
                ss << x;

                return ss.str();
            }
        };

        template <>
        struct DoStringify<std::string>
        {
            static std::string stringify(const std::string & x)
            {
                return x;
            }
        };
    }

    template <typename T_>
    std::string stringify(const T_ & x)
    {
        return implementation::DoStringify<T_>::stringify(x);
    }
}

#endif
