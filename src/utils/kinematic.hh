/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_KINEMATIC_HH
#define WFITTER_GUARD_SRC_UTILS_KINEMATIC_HH 1

namespace wf
{
    inline double lambda(const double & a, const double & b, const double & c)
    {
        return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c);
    }
}

#endif
