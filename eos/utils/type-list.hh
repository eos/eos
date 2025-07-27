/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2008-2010 Ciaran McCreesh
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

#ifndef EOS_GUARD_EOS_UTILS_TYPE_LIST_HH
#define EOS_GUARD_EOS_UTILS_TYPE_LIST_HH 1

#include <eos/utils/type-list-fwd.hh>

namespace eos
{
    struct TypeListTail
    {};

    template <typename Item_, typename Tail_> struct TypeListEntry
    {
            using Item = Item_;
            using Tail = Tail_;
    };

    template <> struct MakeTypeList<>
    {
            using Type = TypeListTail;
    };

    template <typename H_, typename... T_> struct MakeTypeList<H_, T_...>
    {
            using Type = TypeListEntry<H_, typename MakeTypeList<T_...>::Type>;
    };

    template <> struct MakeTypeListConstEntry<TypeListTail>
    {
            using Type = TypeListTail;
    };

    template <typename Item_, typename Tail_> struct MakeTypeListConstEntry<TypeListEntry<Item_, Tail_>>
    {
            using Type = TypeListEntry<const Item_, typename MakeTypeListConstEntry<Tail_>::Type>;
    };

    template <typename TypeList_> struct MakeTypeListConst
    {
            using Type = typename MakeTypeListConstEntry<TypeList_>::Type;
    };

    template <typename Item_> struct TypeListContains<TypeListTail, Item_>
    {
            enum
            {
                value = 0
            };
    };

    template <typename Item_, typename Tail_> struct TypeListContains<TypeListEntry<Item_, Tail_>, Item_>
    {
            enum
            {
                value = 1
            };
    };

    template <typename NotItem_, typename Item_, typename Tail_> struct TypeListContains<TypeListEntry<NotItem_, Tail_>, Item_>
    {
            enum
            {
                value = TypeListContains<Tail_, Item_>::value
            };
    };
} // namespace eos

#endif
