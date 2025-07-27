/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_STRINGIFY_HH
#define EOS_GUARD_EOS_UTILS_STRINGIFY_HH 1

#include <sstream>
#include <string>

namespace eos
{
    namespace implementation
    {
        template <typename T_> struct DoStringify
        {
                static std::string
                stringify(const T_ & x, unsigned precision)
                {
                    std::stringstream ss;
                    ss.precision(precision);
                    ss << x;

                    return ss.str();
                }
        };

        template <> struct DoStringify<std::string>
        {
                static std::string
                stringify(const std::string & x, unsigned)
                {
                    return x;
                }
        };
    } // namespace implementation

    /*!
     * Stringify an arbritrary (scalar) data type.
     *
     * @param x         Object that shall be stringified.
     * @param precision (Optional) floating point precision for the stringification.
     */
    template <typename T_>
    std::string
    stringify(const T_ & x, const unsigned & precision = 10)
    {
        return implementation::DoStringify<T_>::stringify(x, precision);
    }

    /*!
     * Stringify a range of iterators.
     *
     * @param begin     Iterator pointing to the first element of the range.
     * @param end       Iterator pointing beyond the last element of the range.
     * @param precision (Optional) floating point precision for the stringification.
     */
    template <typename Iterator_>
    std::string
    stringify(const Iterator_ & begin, const Iterator_ & end, const unsigned & precision = 10)
    {
        std::stringstream ss;
        ss.precision(precision);
        ss << '(';

        for (Iterator_ i = begin; i != end; ++i)
        {
            ss << ' ' << *i;
        }

        ss << " )";

        return ss.str();
    }

    /*!
     * Stringify a C-style square matrix
     *
     * @param m         Pointer to the first matrix element.
     * @param dim       Dimension of the matrix.
     * @param precision (Optional) floating point precision for the stringification.
     */
    template <typename T_>
    std::string
    stringify(const T_ * m, const unsigned & dim, const unsigned & precision = 10)
    {
        std::stringstream ss;
        ss.precision(precision);
        ss << "\n(";

        for (unsigned i = 0; i < dim; ++i)
        {
            ss << '(';

            for (unsigned j = 0; j < dim; ++j)
            {
                ss << m[i * dim + j];

                if (j != dim - 1)
                {
                    ss << ", ";
                }
            }

            ss << ")\n";
        }

        ss << " )";

        return ss.str();
    }

    template <typename Container_>
    std::string
    stringify_container(const Container_ & container, unsigned precision = 10)
    {
        return stringify(container.begin(), container.end(), precision);
    }
} // namespace eos

#endif
