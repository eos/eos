/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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

#ifndef BTOPI_FF_GUARD_EOS_UTILS_DERIVATIVE_HH
#define BTOPI_FF_GUARD_EOS_UTILS_DERIVATIVE_HH 1

#include <functional>

namespace eos
{
    /*
     * Compute the derivative of a given function f with respect
     * to the only variable x, at point x0. The method of derivation
     * is governed by the template argument Method_.
     */
    template <unsigned order_, typename Method_>
    double derivative(const std::function<double (const double &)> & f, const double & x0);

    namespace deriv
    {
        class TwoSided{ };
    }

    template <>
    double derivative<1u, deriv::TwoSided>(const std::function<double (const double &)> & f, const double & x0);
    template <>
    double derivative<2u, deriv::TwoSided>(const std::function<double (const double &)> & f, const double & x0);
}

#endif
