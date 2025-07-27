/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2007-2011 Ciaran McCreesh
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

#ifndef EOS_GUARD_EOS_UTILS_INDIRECT_ITERATOR_IMPL_HH
#define EOS_GUARD_EOS_UTILS_INDIRECT_ITERATOR_IMPL_HH 1

#include <eos/utils/indirect-iterator.hh>

namespace eos
{
    template <typename Iter_, typename Value_>
    IndirectIterator<Iter_, Value_>::IndirectIterator() :
        _iter()
    {
    }

    template <typename Iter_, typename Value_>
    IndirectIterator<Iter_, Value_>::IndirectIterator(const IndirectIterator & i) :
        _iter(i._iter)
    {
    }

    template <typename Iter_, typename Value_>
    IndirectIterator<Iter_, Value_>::IndirectIterator(const Iter_ & i) :
        _iter(i)
    {
    }

    template <typename Iter_, typename Value_>
    IndirectIterator<Iter_, Value_> &
    IndirectIterator<Iter_, Value_>::operator= (const IndirectIterator & i)
    {
        _iter = i._iter;
        return *this;
    }

    template <typename Iter_, typename Value_>
    IndirectIterator<Iter_, Value_> &
    IndirectIterator<Iter_, Value_>::operator++ ()
    {
        ++_iter;
        return *this;
    }

    template <typename Iter_, typename Value_>
    IndirectIterator<Iter_, Value_>
    IndirectIterator<Iter_, Value_>::operator++ (int)
    {
        IndirectIterator result(*this);
        ++_iter;
        return result;
    }

    template <typename Iter_, typename Value_>
    typename IndirectIterator<Iter_, Value_>::pointer
    IndirectIterator<Iter_, Value_>::operator->() const
    {
        return &**_iter;
    }

    template <typename Iter_, typename Value_>
    typename IndirectIterator<Iter_, Value_>::reference
    IndirectIterator<Iter_, Value_>::operator* () const
    {
        return **_iter;
    }

    template <typename Iter_, typename Value_>
    typename IndirectIterator<Iter_, Value_>::underlying_iterator_type
    IndirectIterator<Iter_, Value_>::underlying_iterator()
    {
        return _iter;
    }

    template <typename Iter_, typename Value_>
    bool
    operator== (const IndirectIterator<Iter_, Value_> & lhs, const IndirectIterator<Iter_, Value_> & rhs)
    {
        return lhs._iter == rhs._iter;
    }

    template <typename Iter_, typename Value_>
    bool
    operator!= (const IndirectIterator<Iter_, Value_> & lhs, const IndirectIterator<Iter_, Value_> & rhs)
    {
        return lhs._iter != rhs._iter;
    }

    template <typename Iter_, typename Value_>
    bool
    operator< (const IndirectIterator<Iter_, Value_> & lhs, const IndirectIterator<Iter_, Value_> & rhs)
    {
        return lhs._iter < rhs._iter;
    }

    template <typename Iter_, typename Value_>
    bool
    operator> (const IndirectIterator<Iter_, Value_> & lhs, const IndirectIterator<Iter_, Value_> & rhs)
    {
        return rhs._iter < lhs._iter;
    }

    template <typename Iter_>
    IndirectIterator<Iter_>
    indirect_iterator(const Iter_ & i)
    {
        return IndirectIterator<Iter_>(i);
    }
} // namespace eos

#endif
