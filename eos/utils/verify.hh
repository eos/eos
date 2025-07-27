/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_UTILS_VERIFY_HH
#define EOS_GUARD_EOS_UTILS_VERIFY_HH 1

#include <eos/utils/exception.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    /*!
     * VerifiedRangeError is thrown whenever a VerifiedRange is assigned a value
     * that exceeds its allowed range.
     */
    class VerifiedRangeError : public Exception
    {
        protected:
            VerifiedRangeError(const std::string & message);
    };

    /*!
     * VerifiedRangeOverflow is thrown whenever a VerifiedRange is assigned a value
     * that is larger than its allowed maximum.
     */
    struct VerifiedRangeOverflow : public VerifiedRangeError
    {
            VerifiedRangeOverflow(const std::string & value, const std::string & maximum);
    };

    /*!
     * VerifiedRangeUnderflow is thrown whenever a VerifiedRange is assigned a value
     * that is smaller than its allowed minimum.
     */
    struct VerifiedRangeUnderflow : public VerifiedRangeError
    {
            VerifiedRangeUnderflow(const std::string & value, const std::string & minimum);
    };

    /*!
     * VerifiedRange is a wrapper around a variable that allows assignment
     * only within a given range.
     *
     * @code
     * VerifiedRange<double> value_between_zero_and_one(0.0, 1.0, 0.5);
     * test_range = +10; // not allowed, throws a VerifiedRangeOverflow
     * test_range = -10; // not allowed, throws a VerifiedRangeUnderflow
     * @endcode
     */
    template <typename T_> class VerifiedRange
    {
        private:
            T_ _min;

            T_ _max;

            T_ _value;

            const T_ &
            _verify(const T_ & t)
            {
                if (t < _min)
                {
                    throw VerifiedRangeUnderflow(stringify(t), stringify(_min));
                }

                if (t > _max)
                {
                    throw VerifiedRangeOverflow(stringify(t), stringify(_max));
                }

                return t;
            }

        public:
            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param min   Minimal value that is allowed for any value of this object.
             * @param max   Maximal value that is allowed for any value of this object.
             * @param value Initial value of this object.
             */
            VerifiedRange(const T_ & min, const T_ & max, const T_ & value) :
                _min(min),
                _max(max),
                _value(_verify(value))
            {
            }

            ///@}

            ///@name Operations
            ///@{
            /// Conversion operator its native type T_
            operator T_ () const { return _value; }

            /*!
             * Assignment operator.
             *
             * @param rhs Right-hand-side of an assignment to this object.
             */
            const T_ &
            operator= (const T_ & rhs)
            {
                _value = _verify(rhs);
                return _value;
            }

            /// Retrieve the minimal value that is allowed for this object.
            T_
            min() const
            {
                return _min;
            }

            /// Retrieve the maximal value that is allowed for this object.
            T_
            max() const
            {
                return _max;
            }

            ///@}
    };
} // namespace eos

#endif
