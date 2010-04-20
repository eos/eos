/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_WILSON_COEFFICIENTS_HH
#define WFITTER_GUARD_SRC_UTILS_WILSON_COEFFICIENTS_HH 1

#include <src/utils/parameters.hh>

namespace wf
{
    void calculate_wilson_coefficients(const double & mu, Parameters & parameters);
}

#endif
