/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_STRONG_TYPEDEF_HH
#define EOS_GUARD_EOS_UTILS_STRONG_TYPEDEF_HH 1

namespace eos
{
    /*! \brief A strong typedef for creating type-safe identifiers. */
    template <typename T_, typename Tag_> class StrongTypedef
    {
        private:
            T_ _value;

        public:
            StrongTypedef() = default;

            explicit StrongTypedef(const T_ & value) :
                _value(value)
            {
            }

            StrongTypedef(const StrongTypedef & other) = default;

            StrongTypedef(StrongTypedef && other) = default;

            StrongTypedef & operator= (const StrongTypedef & other) = default;

            explicit
            operator T_ & () noexcept
            {
                return _value;
            }

            explicit
            operator const T_ & () const noexcept
            {
                return _value;
            }

            const T_ &
            value() const
            {
                return _value;
            }

            friend bool
            operator== (const StrongTypedef & lhs, const StrongTypedef & rhs)
            {
                return lhs._value == rhs._value;
            }
    };
} // namespace eos

#endif /* EOS_GUARD_EOS_UTILS_STRONG_TYPEDEF_HH */
