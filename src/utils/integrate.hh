/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_INTEGRATE_HH
#define EOS_GUARD_SRC_UTILS_INTEGRATE_HH 1

#include <src/utils/complex.hh>

#include <functional>

namespace eos
{
    // TODO: Do a proper MonteCarlo integration
    double integrate(const std::function<double (const double &)> & f, unsigned n, const double & a, const double & b);

    complex<double> integrate(const std::function<complex<double> (const double &)> & f, unsigned n, const double & a, const double & b);
}

#endif
