/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_JOIN_HH
#define EOS_GUARD_EOS_UTILS_JOIN_HH 1

#include <eos/utils/stringify.hh>

namespace eos
{
    /**
     * Join a range of iterators [it, end) with the separator sep.
     */
    template <typename Iter_>
    std::string
    join(Iter_ it, const Iter_ & end, const std::string & sep = ", ")
    {
        if (it == end)
        {
            return {};
        }

        std::stringstream ss;

        ss << *(it++);

        for (; it != end; ++it)
        {
            ss << sep << *it;
        }

        return ss.str();
    }
} // namespace eos

#endif
