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

#ifndef EOS_GUARD_EOS_UTILS_INDIRECT_ITERATOR_HH
#define EOS_GUARD_EOS_UTILS_INDIRECT_ITERATOR_HH 1

#include <eos/utils/indirect-iterator-fwd.hh>

#include <functional>
#include <memory>
#include <type_traits>

namespace eos
{
    template <typename T_> struct IndirectIteratorValueType
    {
            using Type = typename std::iterator_traits<T_>::value_type;
    };

    template <typename T_> struct IndirectIteratorValueType<T_ *>
    {
            using Type = T_;
    };

    template <typename T_> struct IndirectIteratorValueType<std::shared_ptr<T_>>
    {
            using Type = T_;
    };

    template <typename T_> struct IndirectIteratorValueType<std::shared_ptr<const T_>>
    {
            using Type = const T_;
    };

    template <typename T_> struct IndirectIteratorValueType<const T_>
    {
            using Type = typename IndirectIteratorValueType<T_>::Type;
    };

    template <typename T_> struct IndirectIteratorValueType<T_ &>
    {
            using Type = typename IndirectIteratorValueType<T_>::Type;
    };

    /**
     * An IndirectIterator turns an iterator over T_ * or std::shared_ptr<T_> into an iterator
     * over T_.
     */
    template <typename Iter_, typename Value_> class IndirectIterator
    {
            friend bool operator== <>(const IndirectIterator &, const IndirectIterator &);
            friend bool operator!= <>(const IndirectIterator &, const IndirectIterator &);
            friend bool operator< <>(const IndirectIterator &, const IndirectIterator &);
            friend bool operator><> (const IndirectIterator &, const IndirectIterator &);

        private:
            Iter_ _iter;

        public:
            ///@name Basic operations
            ///@{

            IndirectIterator();
            IndirectIterator(const IndirectIterator &);
            IndirectIterator(const Iter_ &);

            IndirectIterator & operator= (const IndirectIterator &);

            ///@}

            ///@name Standard library typedefs
            ///@{

            typedef typename std::remove_reference<Value_>::type & value_type;
            typedef typename std::remove_reference<Value_>::type & reference;
            typedef typename std::remove_reference<Value_>::type * pointer;
            using difference_type   = std::ptrdiff_t;
            using iterator_category = std::forward_iterator_tag;

            ///@}

            ///@name Additional typedefs
            ///@{

            using underlying_iterator_type = Iter_;

            ///@}

            ///@name Increment
            ///@{

            IndirectIterator & operator++ ();
            IndirectIterator   operator++ (int);

            ///@}

            ///@name Dereference
            ///@{

            pointer   operator->() const;
            reference operator* () const;

            underlying_iterator_type underlying_iterator();

            ///@}
    };

    /**
     * Construct an IndirectIterator from another iterator.
     */
    template <typename Iter_> IndirectIterator<Iter_> indirect_iterator(const Iter_ & t);
} // namespace eos

#endif
