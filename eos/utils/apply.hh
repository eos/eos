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

#include <functional>
#include <tuple>

namespace eos
{
    namespace impl
    {
        template <typename Result_, unsigned n_> struct Applicator
        {
            template <typename Function_, typename Tuple_, typename ... CallArgs_>
            static Result_ apply(const Function_ & f, const Tuple_ & t, CallArgs_ ... a)
            {
                return Applicator<Result_, n_ - 1>::apply(f, t, std::get<n_ - 1>(t), a ...);
            }
        };

        template <typename Result_> struct Applicator<Result_, 0>
        {
            template <typename Function_, typename Tuple_, typename ... CallArgs_>
            static Result_ apply(const Function_ & f, const Tuple_ &, CallArgs_ ... a)
            {
                return f(a ...);
            }
        };

        template <> struct Applicator<void, 0>
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

        return impl::Applicator<Result_, sizeof...(FunctionArgs_)>::apply(f, t);
    }

    /* function wrapped in std::function<> */
    template <typename Result_, typename ... FunctionArgs_, typename ... TupleElements_>
    Result_ apply(const std::function<Result_ (FunctionArgs_ ...)> & f, const std::tuple<TupleElements_ ...> & t)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(TupleElements_), "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::Applicator<Result_, sizeof...(FunctionArgs_)>::apply(f, t);
    }

    /* pointer to member function */
    template <typename Result_, typename Class_, typename ... FunctionArgs_, typename ... TupleElements_>
    Result_ apply(Result_ (Class_::* f)(FunctionArgs_ ...), const std::tuple<Class_ *, TupleElements_ ...> & t)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(TupleElements_), "Cannot apply function of N parameters to tuple of M elements, N != M");

        return impl::Applicator<Result_, 1 + sizeof...(FunctionArgs_)>::apply(std::function<Result_ (Class_ *, FunctionArgs_ & ...)>(std::mem_fn(f)), t);
    }
}

#endif
