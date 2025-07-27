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

#ifndef EOS_GUARD_EOS_UTILS_NAMED_VALUE_HH
#define EOS_GUARD_EOS_UTILS_NAMED_VALUE_HH 1

#include <eos/utils/named-value-fwd.hh>

#include <string>
#include <type_traits>
#include <utility>

namespace eos
{
    /**
     * A NamedValue is used to hold a member of type V_ for a class.
     *
     * NamedValue is used to simplify 'plain old data' style classes, and to
     * provide compiler-time-checked named parameters for functions. Use
     * thestruct.themember() and thestruct.themember() = to access the real
     * underlying values.
     *
     * Usually a struct containing NamedValue objects will be constructed using
     * the make_named_values() function. For each NamedValue object,
     * make_named_values() takes a parameter in the form
     * <code>n::whatever_K_is() = the_value</code>.
     *
     * In all cases, NamedValue members are listed in name-sorted order, and
     * the same name is used for K_ and the member name.
     */
    template <typename K_, typename V_> class NamedValue
    {
            static_assert(! std::is_reference<V_>::value, "Tried to make a NamedValue hold a reference");

        private:
            V_ _value;

        public:
            using KeyType   = K_;
            using ValueType = V_;

            template <typename T_>
            NamedValue(const NamedValue<K_, T_> & v) :
                _value(v())
            {
            }

            template <typename T_>
            NamedValue(NamedValue<K_, T_> && v) :
                _value(std::move(v()))
            {
            }

            explicit NamedValue(const V_ & v) :
                _value(v)
            {
            }

            explicit NamedValue(V_ && v) :
                _value(v)
            {
            }

            NamedValue(const NamedValue & v) :
                _value(v._value)
            {
            }

            NamedValue(NamedValue && v) :
                _value(std::move(v._value))
            {
            }

            NamedValue &
            operator= (const NamedValue & v)
            {
                _value = v._value;
                return *this;
            }

            V_ &
            operator() ()
            {
                return _value;
            }

            const V_ &
            operator() () const
            {
                return _value;
            }
    };

    /**
     * A Name is used to make the assignment for NamedValue keys work.
     */
    template <typename T_> class Name
    {
        public:
            template <typename V_>
            NamedValue<Name<T_>, V_>
            operator= (const V_ & v) const
            {
                return NamedValue<Name<T_>, V_>(v);
            }

            template <typename V_>
            NamedValue<Name<T_>, typename std::remove_reference<V_>::type>
            operator= (V_ && v) const
            {
                return NamedValue<Name<T_>, typename std::remove_reference<V_>::type>(v);
            }

            NamedValue<Name<T_>, std::string>
            operator= (const char * const v) const
            {
                return NamedValue<Name<T_>, std::string>(std::string(v));
            }
    };

    /**
     * Bizarre oddity in C++98: you can only use an initialiser list when using
     * equals to initialise a newly constructed object. C++0x fixes this, but GCC
     * 4.4 is buggy, so for now we can't use braces directly...
     */
    template <typename R_, typename... T_>
    R_
    make_named_values(T_ &&... a)
    {
        R_ result = { std::forward<T_>(a)... };
        return result;
    }
} // namespace eos

#endif
