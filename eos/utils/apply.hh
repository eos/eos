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

#ifndef EOS_GUARD_SRC_UTILS_APPLY_HH
#define EOS_GUARD_SRC_UTILS_APPLY_HH 1

#include <array>
#include <functional>
#include <tuple>

namespace eos
{
    /* apply(function, std::tuple<...>) */

    namespace impl
    {
        template <typename Result_, unsigned n_> struct TupleApplicator
        {
            template <typename Function_, typename Tuple_, typename ... CallArgs_>
            static Result_ apply(const Function_ & f, const Tuple_ & t, CallArgs_ ... a)
            {
                return TupleApplicator<Result_, n_ - 1>::apply(f, t, std::get<n_ - 1>(t), a ...);
            }
        };

        template <typename Result_> struct TupleApplicator<Result_, 0>
        {
            template <typename Function_, typename Tuple_, typename ... CallArgs_>
            static Result_ apply(const Function_ & f, const Tuple_ &, CallArgs_ ... a)
            {
                return f(a ...);
            }
        };

        template <> struct TupleApplicator<void, 0>
        {
            template <typename Function_, typename Tuple_, typename ... CallArgs_>
            static void apply(const Function_ & f, const Tuple_ &, CallArgs_ ... a)
            {
                f(a ...);
            }
        };
    }

    /* freestanding function */
    template <typename Result_, typename ... FunctionArgs_, typename ... TupleElements_>
    Result_ apply(Result_ (* f)(FunctionArgs_ ...), const std::tuple<TupleElements_ ...> & t)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(TupleElements_), "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::TupleApplicator<Result_, sizeof...(FunctionArgs_)>::apply(f, t);
    }

    /* function wrapped in std::function<> */
    template <typename Result_, typename ... FunctionArgs_, typename ... TupleElements_>
    Result_ apply(const std::function<Result_ (FunctionArgs_ ...)> & f, const std::tuple<TupleElements_ ...> & t)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(TupleElements_), "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::TupleApplicator<Result_, sizeof...(FunctionArgs_)>::apply(f, t);
    }

    /* pointer to member function */
    template <typename Result_, typename Class_, typename ... FunctionArgs_, typename ... TupleElements_>
    Result_ apply(Result_ (Class_::* f)(FunctionArgs_ ...), const std::tuple<Class_ *, TupleElements_ ...> & t)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(TupleElements_), "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::TupleApplicator<Result_, 1 + sizeof...(FunctionArgs_)>::apply(std::function<Result_ (Class_ *, FunctionArgs_ & ...)>(std::mem_fn(f)), t);
    }

    /* apply(function, std::array<T_, dim_>) */

    namespace impl
    {
        template <typename Result_, unsigned long m_, unsigned long n_> struct ArrayApplicator
        {
            template <typename Function_, typename Array_, typename ... CallArgs_>
            static Result_ apply(const Function_ & f, const Array_ & a, CallArgs_ ... args)
            {
                return ArrayApplicator<Result_, m_, n_ - 1>::apply(f, a, args ..., a[m_ - n_]);
            }
        };

        template <typename Result_, unsigned long m_> struct ArrayApplicator<Result_, m_, 0>
        {
            template <typename Function_, typename Array_, typename ... CallArgs_>
            static Result_ apply(const Function_ & f, const Array_ &, CallArgs_ ... args)
            {
                return f(args ...);
            }
        };

        template <unsigned long m_> struct ArrayApplicator<void, m_, 0>
        {
            template <typename Function_, typename Array_, typename ... CallArgs_>
            static void apply(const Function_ & f, const Array_ &, CallArgs_ ... args)
            {
                f(args ...);
            }
        };
    }

    /* freestanding function */
    template <typename Result_, typename ArrayBaseType_, unsigned long n_, typename ... FunctionArgs_>
    Result_ apply(Result_ (* f)(FunctionArgs_ ...), const std::array<ArrayBaseType_, n_> & a)
    {
        static_assert(sizeof...(FunctionArgs_) == n_, "Cannot apply function of N parameters to array of M elements, N != M");

        return impl::ArrayApplicator<Result_, n_, n_>::apply(f, a);
    }

    /* function wrapped in std::function<> */
    template <typename Result_, typename ArrayBaseType_, unsigned long n_, typename ... FunctionArgs_>
    Result_ apply(const std::function<Result_ (FunctionArgs_ ...)> & f, const std::array<ArrayBaseType_, n_> & a)
    {
        static_assert(sizeof...(FunctionArgs_) == n_, "Cannot apply function of N parameters to array of M elements, N != M");

        return impl::ArrayApplicator<Result_, n_, n_>::apply(f, a);
    }

    /* pointer to member function */
    template <typename Result_, typename ArrayBaseType_, unsigned long n_, typename Class_, typename ... FunctionArgs_>
    Result_ apply(Result_ (Class_::* f)(FunctionArgs_ ...), Class_ * c, const std::array<ArrayBaseType_, n_> & a)
    {
        static_assert(sizeof...(FunctionArgs_) == n_, "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::ArrayApplicator<Result_, n_, n_>::apply(std::function<Result_ (Class_ *, FunctionArgs_ & ...)>(std::mem_fn(f)), a, c);
    }

    /* pointer to member function */
    template <typename Result_, typename ArrayBaseType_, typename Class_>
    Result_ apply(Result_ (Class_::* f)() const, const Class_ * c, const std::array<ArrayBaseType_, 0ul> &)
    {
        return std::mem_fn(f)(c);
    }

    /* pointer to member function wrapped in std::function<> */
    template <typename Result_, typename ArrayBaseType_, unsigned long n_, typename Class_, typename ... FunctionArgs_>
    Result_ apply(const std::function<Result_ (Class_ *, FunctionArgs_ ...)> & f, Class_ * c, const std::array<ArrayBaseType_, n_> & a)
    {
        static_assert(sizeof...(FunctionArgs_) == n_, "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::ArrayApplicator<Result_, n_, n_>::apply(f, a, c);
    }
}

#endif
