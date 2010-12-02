/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_INTEGRATE_HH
#define EOS_GUARD_SRC_UTILS_INTEGRATE_HH 1

#include <src/utils/complex.hh>

#include <functional>

namespace eos
{
    /// @{
    /*!
     * Numerically integrate functions of one real-valued parameter.
     *
     * Uses the Delta^2-Rule by Aitkin to refine the result.
     *
     * @param f      Integrand.
     * @param n      Number of evaluations, must be a power of 2.
     * @param a      Lower limit of the domain of integration.
     * @param b      Upper limit of the domain of integration.
     */
    double integrate(const std::function<double (const double &)> & f, unsigned n, const double & a, const double & b);
    complex<double> integrate(const std::function<complex<double> (const double &)> & f, unsigned n, const double & a, const double & b);
    /// @}
}

#endif
