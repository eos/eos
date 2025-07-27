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

#ifndef EOS_GUARD_EOS_UTILS_PRIVATE_IMPLEMENTATION_PATTERN_IMPL_HH
#define EOS_GUARD_EOS_UTILS_PRIVATE_IMPLEMENTATION_PATTERN_IMPL_HH 1

#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    template <typename T_>
    PrivateImplementationPattern<T_>::PrivateImplementationPattern(Implementation<T_> * imp) :
        _imp(imp)
    {
    }

    template <typename T_> PrivateImplementationPattern<T_>::~PrivateImplementationPattern() {}
} // namespace eos

#endif
