/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2015, 2018 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_TUPLE_MAKER_HH
#define EOS_GUARD_EOS_UTILS_TUPLE_MAKER_HH 1

#include <eos/utils/kinematic.hh>

#include <array>
#include <tuple>
#include <utility>

namespace eos
{
    namespace impl
    {
        template <typename T_, typename U_> struct ConvertTo
        {
                using Type = U_;
        };

        template <unsigned n_> struct TupleMaker
        {
                template <typename Decay_, typename... TupleElements_, typename... ResultElements_>
                static auto
                make(const Kinematics & k, const std::tuple<TupleElements_...> & t, const Decay_ * d, ResultElements_... r)
                        -> std::tuple<const Decay_ *, typename ConvertTo<TupleElements_, KinematicVariable>::Type...>
                {
                    return TupleMaker<n_ - 1>::make(k, t, d, k[static_cast<const char *>(std::get<n_ - 1>(t))], r...);
                }
        };

        template <> struct TupleMaker<0>
        {
                template <typename Decay_, typename... TupleElements_, typename... ResultElements_>
                static auto
                make(const Kinematics &, const std::tuple<TupleElements_...> &, const Decay_ * d, ResultElements_... r)
                        -> std::tuple<const Decay_ *, typename ConvertTo<TupleElements_, KinematicVariable>::Type...>
                {
                    return std::make_tuple(d, r...);
                }
        };

        template <typename T_> struct TupleSize;

        template <typename... TupleElements_> struct TupleSize<std::tuple<TupleElements_...>>
        {
                static const unsigned long size = sizeof...(TupleElements_);
        };

        template <typename T_, typename Tuple_, std::size_t... Indices_>
        auto
        make_array(const Tuple_ & tuple, std::index_sequence<Indices_...>) -> std::array<T_, std::tuple_size<Tuple_>::value>
        {
            return { { std::get<Indices_>(tuple)... } };
        }

        template <typename T_, typename Tuple_>
        auto
        make_array(const Tuple_ & tuple) -> std::array<T_, std::tuple_size<Tuple_>::value>
        {
            using Indices_ = std::make_index_sequence<std::tuple_size<Tuple_>::value>;

            return make_array<T_>(tuple, Indices_{});
        }

        template <typename T_>
        std::array<T_, 0>
        make_array(const std::tuple<> &)
        {
            return {};
        }
    } // namespace impl
} // namespace eos

#endif
