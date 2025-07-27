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

#ifndef EOS_GUARD_EOS_UTILS_TYPE_LIST_FWD_HH
#define EOS_GUARD_EOS_UTILS_TYPE_LIST_FWD_HH 1

namespace eos
{
    struct TypeListTail;

    template <typename Item_, typename Tail_> struct TypeListEntry;

    template <typename...> struct MakeTypeList;

    template <> struct MakeTypeList<>;

    template <typename H_, typename... T_> struct MakeTypeList<H_, T_...>;

    template <typename TypeList_> struct MakeTypeListConst;

    template <typename Item_> struct MakeTypeListConstEntry;

    template <typename TypeList_, typename Item_> struct TypeListContains;
} // namespace eos
#endif
