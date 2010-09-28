/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_MATRIX_HH
#define WFITTER_GUARD_SRC_UTILS_MATRIX_HH 1

#include <array>

namespace wf
{
    /* Addition */

    /* matrix plus matrix */
    template <unsigned long m_, unsigned long n_>
    std::array<std::array<double, n_>, m_> operator+ (const std::array<std::array<double, n_>, m_> & x,
            const std::array<std::array<double, n_>, m_> & y)
    {
        std::array<std::array<double, n_>, m_> result;
        for (unsigned i(0) ; i < m_ ; ++i)
        {
            for (unsigned j(0) ; j < n_ ; ++j)
            {
                result[i][j] = x[i][j] + y[i][j];
            }
        }

        return result;
    }

    /* matrix minus matrix */
    template <unsigned long m_, unsigned long n_>
    std::array<std::array<double, n_>, m_> operator- (const std::array<std::array<double, n_>, m_> & x,
            const std::array<std::array<double, n_>, m_> & y)
    {
        std::array<std::array<double, n_>, m_> result;
        for (unsigned i(0) ; i < m_ ; ++i)
        {
            for (unsigned j(0) ; j < n_ ; ++j)
            {
                result[i][j] = x[i][j] - y[i][j];
            }
        }

        return result;
    }

    /* vector plus vector */
    template <unsigned long m_>
    std::array<double, m_> operator+ (const std::array<double, m_> & x,
            const std::array<double, m_> & y)
    {
        std::array<double, m_> result;
        for (unsigned i(0) ; i < m_ ; ++i)
        {
            result[i] = x[i] + y[i];
        }

        return result;
    }

    /* Multiplication */

    /* matrix times matrix */
    template <unsigned long m_, unsigned long n_, unsigned long o_>
    std::array<std::array<double, n_>, m_> operator* (const std::array<std::array<double, o_>, m_> & x,
            const std::array<std::array<double, n_>, o_> & y)
    {
        std::array<std::array<double, n_>, m_> result;
        for (unsigned i(0) ; i < m_ ; ++i)
        {
            for (unsigned j(0) ; j < n_ ; ++j)
            {
                result[i][j] = 0.0;

                for (unsigned k(0) ; k < o_ ; ++k)
                {
                    result[i][j] += x[i][k] * y[k][j];
                }
            }
        }

        return result;
    }

    /* matrix times vector */
    template <unsigned long m_, unsigned long n_>
    std::array<double, m_> operator* (const std::array<std::array<double, n_>, m_> & x, const std::array<double, n_> & y)
    {
        std::array<double, m_> result;
        for (unsigned i(0) ; i < m_ ; ++i)
        {
            result[i] = 0.0;

            for (unsigned j(0) ; j < n_ ; ++j)
            {
                result[i] += x[i][j] * y[j];
            }
        }

        return result;
    }

    /* scalar times matrix */
    template <unsigned long m_, unsigned long n_>
    std::array<std::array<double, n_>, m_> operator* (const double & x, const std::array<std::array<double, n_>, m_> & y)
    {
        std::array<std::array<double, n_>, m_> result = y;
        for (unsigned i(0) ; i < m_ ; ++i)
        {
            for (unsigned j(0) ; j < n_ ; ++j)
            {
                result[i][j] *= x;
            }
        }

        return result;
    }

    /* scalar times vector */
    template <unsigned long n_>
    std::array<double, n_> operator* (const double & x, const std::array<double, n_> & y)
    {
        std::array<double, n_> result = y;
        for (unsigned i(0) ; i < n_ ; ++i)
        {
            result[i] *= x;
        }

        return result;
    }
}

#endif
