/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_MATRIX_HH
#define EOS_GUARD_SRC_UTILS_MATRIX_HH 1

#include <array>

namespace eos
{
    /* Addition */

    /* matrix plus matrix */
    template <typename T_, std::size_t m_, std::size_t n_>
    std::array<std::array<T_, n_>, m_> operator+ (const std::array<std::array<T_, n_>, m_> & x,
            const std::array<std::array<T_, n_>, m_> & y)
    {
        std::array<std::array<T_, n_>, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] = x[i][j] + y[i][j];
            }
        }

        return result;
    }

    /* matrix minus matrix */
    template <typename T_, std::size_t m_, std::size_t n_>
    std::array<std::array<T_, n_>, m_> operator- (const std::array<std::array<T_, n_>, m_> & x,
            const std::array<std::array<T_, n_>, m_> & y)
    {
        std::array<std::array<T_, n_>, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] = x[i][j] - y[i][j];
            }
        }

        return result;
    }

    /* vector plus vector */
    template <typename T_, std::size_t m_>
    std::array<T_, m_> operator+ (const std::array<T_, m_> & x,
            const std::array<T_, m_> & y)
    {
        std::array<T_, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            result[i] = x[i] + y[i];
        }

        return result;
    }

    /* vector minus vector */
    template <typename T_, std::size_t m_>
    std::array<T_, m_> operator- (const std::array<T_, m_> & x,
            const std::array<T_, m_> & y)
    {
        std::array<T_, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            result[i] = x[i] - y[i];
        }

        return result;
    }

    /* Multiplication */

    /* matrix times matrix */
    template <typename T_, std::size_t m_, std::size_t n_, std::size_t o_>
    std::array<std::array<T_, n_>, m_> operator* (const std::array<std::array<T_, o_>, m_> & x,
            const std::array<std::array<T_, n_>, o_> & y)
    {
        std::array<std::array<T_, n_>, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] = 0.0;

                for (std::size_t k(0) ; k < o_ ; ++k)
                {
                    result[i][j] += x[i][k] * y[k][j];
                }
            }
        }

        return result;
    }

    /* matrix times vector, right-multiplication */
    template <typename T_, std::size_t m_, std::size_t n_>
    std::array<T_, m_> operator* (const std::array<std::array<T_, n_>, m_> & x, const std::array<T_, n_> & y)
    {
        std::array<T_, m_> result;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            result[i] = 0.0;

            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i] += x[i][j] * y[j];
            }
        }

        return result;
    }

    /* matrix times vector, left-multiplication */
    template <typename T_, std::size_t m_, std::size_t n_>
    std::array<T_, n_> operator^ (const std::array<T_, m_> & x, const std::array<std::array<T_, n_>, m_> & y)
    {
        std::array<T_, n_> result;
        for (std::size_t i(0) ; i < n_ ; ++i)
        {
            result[i] = 0.0;
            for (std::size_t j(0) ; j < m_ ; ++j)
            {
                result[i] += x[j] * y[j][i];
            }
        }

        return result;
    }

    /* scalar times matrix */
    template <typename T_, std::size_t m_, std::size_t n_>
    std::array<std::array<T_, n_>, m_> operator* (const T_ & x, const std::array<std::array<T_, n_>, m_> & y)
    {
        std::array<std::array<T_, n_>, m_> result = y;
        for (std::size_t i(0) ; i < m_ ; ++i)
        {
            for (std::size_t j(0) ; j < n_ ; ++j)
            {
                result[i][j] *= x;
            }
        }

        return result;
    }

    /* scalar times vector */
    template <typename T_, std::size_t n_>
    std::array<T_, n_> operator* (const T_ & x, const std::array<T_, n_> & y)
    {
        std::array<T_, n_> result = y;
        for (std::size_t i(0) ; i < n_ ; ++i)
        {
            result[i] *= x;
        }

        return result;
    }

    /* vector times vector */
    template <typename T_, std::size_t n_>
    T_ dot(const std::array<T_, n_> & x, const std::array<T_, n_> & y)
    {
        T_ result = 0.0;
        for (std::size_t i(0) ; i < n_ ; ++i)
        {
            result += x[i] * y[i];
        }

        return result;
    }

    /* component-wise multiplication of two vectors */
    template <typename T_, std::size_t n_>
    std::array<T_, n_> mult(const std::array<T_, n_> & x, const std::array<T_, n_> & y)
    {
        std::array<T_, n_> result;
        for (std::size_t i(0) ; i < n_ ; ++i)
        {
            result[i] = x[i] * y[i];
        }

        return result;
    }

    /* component-wise division of two vectors */
    template <typename T_, std::size_t n_>
    std::array<T_, n_> divide(const std::array<T_, n_> & x, const std::array<T_, n_> & y)
    {
        std::array<T_, n_> result;
        for (std::size_t i(0) ; i < n_ ; ++i)
        {
            result[i] = x[i] / y[i];
        }

        return result;
    }
}

#endif
