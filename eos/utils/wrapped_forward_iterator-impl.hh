/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Ciaran McCreesh
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

#ifndef EOS_GUARD_UTIL_WRAPPED_FORWARD_ITERATOR_IMPL_HH
#define EOS_GUARD_UTIL_WRAPPED_FORWARD_ITERATOR_IMPL_HH 1

#include <eos/utils/wrapped_forward_iterator.hh>

namespace eos
{
    template <typename Tag_, typename Value_>
    struct WrappedForwardIterator<Tag_, Value_>::Base
    {
        virtual Base * clone() const = 0;
        virtual void increment() = 0;
        virtual typename WrappedForwardIterator<Tag_, Value_>::pointer pointer() const = 0;
        virtual typename WrappedForwardIterator<Tag_, Value_>::reference reference() const = 0;
        virtual bool equal(const Base *) const = 0;
        virtual void * underlying_iterator_ptr() = 0;
        virtual const void * underlying_iterator_ptr() const = 0;

        virtual ~Base()
        {
        }
    };

    template <typename Tag_, typename Value_>
    template <typename Iter_>
    struct WrappedForwardIterator<Tag_, Value_>::BaseImpl :
        WrappedForwardIterator<Tag_, Value_>::Base
    {
        Iter_ i;

        BaseImpl(const Iter_ & ii) :
            i(ii)
        {
        }

        Base * clone() const
        {
            return new BaseImpl(i);
        }

        void increment()
        {
            ++i;
        }

        typename WrappedForwardIterator<Tag_, Value_>::reference reference() const
        {
            return *i;
        }

        typename WrappedForwardIterator<Tag_, Value_>::pointer pointer() const
        {
            return i.operator-> ();
        }

        bool equal(const Base * other) const
        {
            return i == static_cast<const BaseImpl *>(other)->i;
        }

        void * underlying_iterator_ptr()
        {
            return &i;
        }

        const void * underlying_iterator_ptr() const
        {
            return &i;
        }
    };

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>::WrappedForwardIterator() :
        _base(0)
    {
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>::~WrappedForwardIterator()
    {
        delete _base;
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>::WrappedForwardIterator(const WrappedForwardIterator & other) :
        _base(other._base ? other._base->clone() : other._base)
    {
    }

    template <typename Tag_, typename Value_>
    template <typename T_>
    WrappedForwardIterator<Tag_, Value_>::WrappedForwardIterator(const T_ & base) :
        _base(new BaseImpl<T_>(base))
    {
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_> &
    WrappedForwardIterator<Tag_, Value_>::operator= (const WrappedForwardIterator<Tag_, Value_> & other)
    {
        if (this != &other)
            _base = other._base ? other._base->clone() : other._base;
        return *this;
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_> &
    WrappedForwardIterator<Tag_, Value_>::operator++ ()
    {
        _base->increment();
        return *this;
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>
    WrappedForwardIterator<Tag_, Value_>::operator++ (int)
    {
        WrappedForwardIterator result(*this);
        _base->increment();
        return result;
    }

    template <typename Tag_, typename Value_>
    typename WrappedForwardIterator<Tag_, Value_>::pointer
    WrappedForwardIterator<Tag_, Value_>::operator-> () const
    {
        return _base->pointer();
    }

    template <typename Tag_, typename Value_>
    typename WrappedForwardIterator<Tag_, Value_>::reference
    WrappedForwardIterator<Tag_, Value_>::operator* () const
    {
        return _base->reference();
    }

    template <typename Tag_, typename Value_>
    bool
    WrappedForwardIterator<Tag_, Value_>::operator== (const WrappedForwardIterator & other) const
    {
        if (! _base)
            return ! other._base;

        return _base->equal(other._base);
    }

    template <typename Tag_, typename Value_>
    bool
    WrappedForwardIterator<Tag_, Value_>::operator!= (const WrappedForwardIterator & other) const
    {
        if (! _base)
            return 0 != other._base;

        return ! _base->equal(other._base);
    }

    template <typename Tag_, typename Value_>
    template <typename Iter_>
    Iter_ &
    WrappedForwardIterator<Tag_, Value_>::underlying_iterator()
    {
        return *static_cast<Iter_ *>(_base->underlying_iterator_ptr());
    }

    template <typename Tag_, typename Value_>
    template <typename Iter_>
    const Iter_ &
    WrappedForwardIterator<Tag_, Value_>::underlying_iterator() const
    {
        return *static_cast<const Iter_ *>(_base->underlying_iterator_ptr());
    }
}

#endif
