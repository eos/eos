/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_INTEGRATE_HH
#define WFITTER_GUARD_SRC_UTILS_INTEGRATE_HH 1

#include <src/utils/complex.hh>

#include <tr1/functional>

namespace wf
{
    // TODO: Do a proper MonteCarlo integration
    double integrate(const std::tr1::function<double (const double &)> & f, unsigned n, const double & a, const double & b);

    Complex<double> integrate(const std::tr1::function<Complex<double> (const double &)> & f, unsigned n, const double & a, const double & b);
}

#endif
