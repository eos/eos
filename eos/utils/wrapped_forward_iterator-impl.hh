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

#ifndef EOS_GUARD_EOS_UTILS_WRAPPED_FORWARD_ITERATOR_IMPL_HH
#define EOS_GUARD_EOS_UTILS_WRAPPED_FORWARD_ITERATOR_IMPL_HH 1

#include <eos/utils/wrapped_forward_iterator.hh>

namespace eos
{
    template <typename T_>
    void
    checked_delete(T_ * const t)
    {
        typedef char     make_sure_type_is_defined[sizeof(T_) ? 1 : -1];
        static const int make_sure_type_is_defined_again __attribute__((unused)) = sizeof(make_sure_type_is_defined);
        delete t;
    }

    template <typename WrappedIter_>
    typename WrappedForwardIteratorTraits<typename WrappedIter_::Tag>::UnderlyingIterator *
    wrapped_underlying_iterator_real_type(const WrappedIter_ &, void * i)
    {
        return reinterpret_cast<typename WrappedForwardIteratorTraits<typename WrappedIter_::Tag>::UnderlyingIterator *>(i);
    }

    template <typename WrappedIter_>
    WrappedForwardIteratorUnderlyingIteratorHolder *
    wrapped_underlying_iterator_hide_real_type(const WrappedIter_ &, typename WrappedForwardIteratorTraits<typename WrappedIter_::Tag>::UnderlyingIterator * const i)
    {
        return reinterpret_cast<WrappedForwardIteratorUnderlyingIteratorHolder *>(i);
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>::WrappedForwardIterator() :
        _iter(wrapped_underlying_iterator_hide_real_type(*this, new typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator))
    {
    }

    template <typename Tag_, typename Value_> WrappedForwardIterator<Tag_, Value_>::~WrappedForwardIterator()
    {
        checked_delete(wrapped_underlying_iterator_real_type(*this, _iter));
    }

    template <typename Tag_, typename Value_, typename Iter_> struct WrappedForwardIteratorGetBase
    {
            template <typename Traits_>
            static typename Traits_::UnderlyingIterator
            real_get_iter(const typename Traits_::UnderlyingIterator & i)
            {
                return i;
            }

            template <typename Traits_>
            static typename Traits_::UnderlyingIterator
            real_get_iter(const typename Traits_::EquivalentNonConstIterator & i)
            {
                return typename Traits_::UnderlyingIterator(
                        i.template underlying_iterator<
                                typename WrappedForwardIteratorTraits<typename WrappedForwardIteratorTraits<Tag_>::EquivalentNonConstIterator::Tag>::UnderlyingIterator>());
            }

            static typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator
            get_iter(const Iter_ & i)
            {
                return real_get_iter<WrappedForwardIteratorTraits<Tag_>>(i);
            }
    };

    template <typename Tag_, typename Value_> struct WrappedForwardIteratorGetBase<Tag_, Value_, WrappedForwardIterator<Tag_, Value_>>
    {
            static WrappedForwardIterator<Tag_, Value_>
            get_iter(const WrappedForwardIterator<Tag_, Value_> & i)
            {
                return i.template underlying_iterator<typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator>();
            }
    };

    template <typename LHS_, typename RHS_>
    void
    set_wrapped_forward_iterator_iterator(LHS_ & lhs, RHS_ rhs)
    {
        lhs = rhs;
    }

    template <typename Tag_, typename Value_>
    template <typename T_>
    WrappedForwardIterator<Tag_, Value_>::WrappedForwardIterator(const T_ & iter) :
        _iter(wrapped_underlying_iterator_hide_real_type(*this, new typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator))
    {
        set_wrapped_forward_iterator_iterator(*wrapped_underlying_iterator_real_type(*this, _iter), WrappedForwardIteratorGetBase<Tag_, Value_, T_>::get_iter(iter));
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>::WrappedForwardIterator(const WrappedForwardIterator & other) :
        _iter(wrapped_underlying_iterator_hide_real_type(
                *this, new typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator(*wrapped_underlying_iterator_real_type(other, other._iter))))
    {
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_> &
    WrappedForwardIterator<Tag_, Value_>::operator= (const WrappedForwardIterator & other)
    {
        *wrapped_underlying_iterator_real_type(*this, _iter) = *wrapped_underlying_iterator_real_type(other, other._iter);
        return *this;
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_> &
    WrappedForwardIterator<Tag_, Value_>::operator++ ()
    {
        ++*wrapped_underlying_iterator_real_type(*this, _iter);
        return *this;
    }

    template <typename Tag_, typename Value_>
    WrappedForwardIterator<Tag_, Value_>
    WrappedForwardIterator<Tag_, Value_>::operator++ (int)
    {
        WrappedForwardIterator result(*this);
        operator++ ();
        return result;
    }

    namespace impl
    {
        template <typename R_, typename T_, bool is_pointer> struct GetPointer;

        template <typename R_, typename T_> struct GetPointer<R_, T_, true>
        {
                static R_
                get_pointer(T_ * t)
                {
                    return *t;
                }
        };

        template <typename R_, typename T_> struct GetPointer<R_, T_, false>
        {
                static R_
                get_pointer(T_ * t)
                {
                    return t->operator->();
                }
        };

        template <typename R_, typename T_>
        R_
        get_pointer(T_ * t)
        {
            return GetPointer<R_, T_, std::is_pointer<T_>::value>::get_pointer(t);
        }
    } // namespace impl

    template <typename Tag_, typename Value_>
    typename WrappedForwardIterator<Tag_, Value_>::pointer
    WrappedForwardIterator<Tag_, Value_>::operator->() const
    {
        return impl::get_pointer<WrappedForwardIterator<Tag_, Value_>::pointer, typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator>(
                wrapped_underlying_iterator_real_type(*this, _iter));
    }

    namespace impl
    {
        template <typename R_, typename T_, bool is_pointer> struct GetReference;

        template <typename R_, typename T_> struct GetReference<R_, T_, true>
        {
                static R_
                get_reference(T_ * t)
                {
                    return **t;
                }
        };

        template <typename R_, typename T_> struct GetReference<R_, T_, false>
        {
                static R_
                get_reference(T_ * t)
                {
                    return t->operator* ();
                }
        };

        template <typename R_, typename T_>
        R_
        get_reference(T_ * t)
        {
            return GetReference<R_, T_, std::is_pointer<T_>::value>::get_reference(t);
        }
    } // namespace impl

    template <typename Tag_, typename Value_>
    typename WrappedForwardIterator<Tag_, Value_>::reference
    WrappedForwardIterator<Tag_, Value_>::operator* () const
    {
        return impl::get_reference<WrappedForwardIterator<Tag_, Value_>::reference, typename WrappedForwardIteratorTraits<Tag_>::UnderlyingIterator>(
                wrapped_underlying_iterator_real_type(*this, _iter));
    }

    template <typename Tag_, typename Value_>
    bool
    WrappedForwardIterator<Tag_, Value_>::operator== (const WrappedForwardIterator & other) const
    {
        return *wrapped_underlying_iterator_real_type(*this, _iter) == *wrapped_underlying_iterator_real_type(other, other._iter);
    }

    template <typename Tag_, typename Value_>
    bool
    WrappedForwardIterator<Tag_, Value_>::operator!= (const WrappedForwardIterator & other) const
    {
        return *wrapped_underlying_iterator_real_type(*this, _iter) != *wrapped_underlying_iterator_real_type(other, other._iter);
    }

    template <typename Tag_, typename Value_>
    template <typename T_>
    T_ &
    WrappedForwardIterator<Tag_, Value_>::underlying_iterator()
    {
        return *wrapped_underlying_iterator_real_type(*this, _iter);
    }

    template <typename Tag_, typename Value_>
    template <typename T_>
    const T_ &
    WrappedForwardIterator<Tag_, Value_>::underlying_iterator() const
    {
        return *wrapped_underlying_iterator_real_type(*this, _iter);
    }
} // namespace eos

#endif
