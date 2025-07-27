/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2022 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#ifndef EOS_GUARD_EOS_UTILS_MEMOISE_HH
#define EOS_GUARD_EOS_UTILS_MEMOISE_HH 1

#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/instantiation_policy.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/mutex.hh>

#include <cstdint>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <vector>

/* We are using a custom std::hash for tuple<FunctionPtr, double, ..., double> */
namespace std
{
    template <typename U_>
    static uint64_t
    __hash_one(const U_ & u)
    {
        static_assert(sizeof(U_) == sizeof(uint32_t) || sizeof(U_) == sizeof(uint64_t), "Need to specialize __hash_one for non 32- and 64-bit data types");

        if (sizeof(U_) == sizeof(int64_t))
        {
            return *reinterpret_cast<const uint64_t *>(&u);
        }
        else if (sizeof(U_) == sizeof(int32_t))
        {
            return static_cast<uint64_t>(*reinterpret_cast<const uint32_t *>(&u));
        }
    }

    template <unsigned n_, typename... T_> struct __TupleHasher
    {
            static uint64_t
            hash(const std::tuple<T_...> & t)
            {
                return __TupleHasher<n_ - 1, T_...>::hash(t) ^ __hash_one<decltype(std::get<n_>(t))>(std::get<n_>(t));
            }
    };

    template <typename... T_> struct __TupleHasher<0, T_...>
    {
            static uint64_t
            hash(const std::tuple<T_...> & t)
            {
                return __hash_one<decltype(std::get<0>(t))>(std::get<0>(t));
            }
    };

    template <typename... T_> struct hash<std::tuple<T_...>>
    {
            size_t
            operator() (const std::tuple<T_...> & t) const
            {
                return __TupleHasher<sizeof...(T_) - 1, T_...>::hash(t);
            }
    };
} // namespace std

namespace eos
{
    namespace implementation
    {
        template <typename T_> struct ResultOf;

        template <typename Result_, typename Class_, typename... Args_> struct ResultOf<Result_ (Class_::*)(Args_...)>
        {
                using Type = Result_;
        };

        template <typename Result_, typename... Args_> struct ResultOf<Result_ (*)(Args_...)>
        {
                using Type = Result_;
        };
    } // namespace implementation

    class MemoisationControl : public InstantiationPolicy<MemoisationControl, Singleton>
    {
        private:
            Mutex * const _mutex;

            std::vector<std::function<void()>> _clear_functions;

        public:
            MemoisationControl();

            ~MemoisationControl();

            void register_clear_function(const std::function<void()> & clear_function);

            void clear();
    };

    template <typename Result_, typename... Params_> class Memoiser : public InstantiationPolicy<Memoiser<Result_, Params_...>, Singleton>
    {
        public:
            using FunctionType = Result_ (*)(const Params_ &...);
            using KeyType      = std::tuple<FunctionType, Params_...>;

        private:
            Mutex * const _mutex;

            std::unordered_map<KeyType, Result_> _memoisations;

        public:
            Memoiser() :
                _mutex(new Mutex)
            {
                MemoisationControl::instance()->register_clear_function(std::bind(&Memoiser<Result_, Params_...>::clear, this));
            }

            ~Memoiser() { delete _mutex; }

            Result_
            operator() (const FunctionType & f, const Params_ &... p)
            {
                Lock l(*_mutex);

                KeyType key(f, p...);
                auto    i = _memoisations.find(key);

                if (_memoisations.end() != i)
                {
                    return i->second;
                }

                Result_ result = f(p...);

                if (_memoisations.size() > 100000u)
                {
                    _memoisations.clear();
                }

                _memoisations.insert(std::pair<KeyType, Result_>(key, result));

                return result;
            }

            void
            clear()
            {
                Lock l(*_mutex);

                _memoisations.clear();
            }

            unsigned
            number_of_memoisations() const
            {
                Lock l(*_mutex);

                return _memoisations.size();
            }
    };

    template <typename FunctionType_, typename... Params>
    typename implementation::ResultOf<FunctionType_>::Type
    memoise(FunctionType_ f, const Params &... p)
    {
        return (*Memoiser<typename implementation::ResultOf<FunctionType_>::Type, Params...>::instance())(f, p...);
    }

    template <typename FunctionType_, typename... Params>
    unsigned
    number_of_memoisations(FunctionType_, const Params &...)
    {
        return Memoiser<typename implementation::ResultOf<FunctionType_>::Type, Params...>::instance()->number_of_memoisations();
    }
} // namespace eos

#endif
