/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 *
 * Copied from the Paludis package manager, which is
 * Copyright (c) 2007-2009 Ciaran McCreesh
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

#ifndef EOS_GUARD_EOS_UTILS_WRAPPED_FORWARD_ITERATOR_HH
#define EOS_GUARD_EOS_UTILS_WRAPPED_FORWARD_ITERATOR_HH 1

#include <eos/utils/wrapped_forward_iterator-fwd.hh>

#include <functional>
#include <iterator>
#include <type_traits>

namespace eos
{
    /**
     * A WrappedForwardIterator is a generic wrapper around a forward iterator,
     * hiding the underlying base iterator.
     */
    template <typename Tag_, typename Value_> class WrappedForwardIterator
    {
        private:
            WrappedForwardIteratorUnderlyingIteratorHolder * _iter;

        public:
            using Tag = Tag_;

            ///@name Basic operations
            ///@{

            WrappedForwardIterator();
            ~WrappedForwardIterator();
            WrappedForwardIterator(const WrappedForwardIterator &);

            template <typename T_> WrappedForwardIterator(const T_ &);

            WrappedForwardIterator & operator= (const WrappedForwardIterator &);

            ///@}

            ///@name Standard library typedefs
            ///@{

            using value_type        = typename std::remove_reference<Value_>::type &;
            using reference         = typename std::remove_reference<Value_>::type &;
            using pointer           = typename std::remove_reference<Value_>::type *;
            using difference_type   = std::ptrdiff_t;
            using iterator_category = std::forward_iterator_tag;

            ///@}

            ///@name Increment
            ///@{

            WrappedForwardIterator & operator++ ();
            WrappedForwardIterator   operator++ (int);

            ///@}

            ///@name Dereference
            ///@{

            pointer   operator->() const;
            reference operator* () const;

            ///@}

            ///@name Equality
            ///@{

            bool operator== (const WrappedForwardIterator &) const;
            bool operator!= (const WrappedForwardIterator &) const;

            ///@}

            ///@name Underlying iterator
            ///@{

            template <typename T_> T_ &       underlying_iterator();
            template <typename T_> const T_ & underlying_iterator() const;

            ///@}
    };
} // namespace eos

#endif
