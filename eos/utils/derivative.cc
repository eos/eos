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

#include <eos/utils/derivative.hh>

#include <limits>
#include <cmath>

namespace eos
{
    template <>
    double derivative<1u, deriv::TwoSided>(const std::function<double (const double &)> & f, const double & x0)
    {
        static const double eps = std::numeric_limits<double>::epsilon(), sqrteps = std::sqrt(eps);

        double h;
        if (x0 > 4.0 * sqrteps)
        {
            h = sqrteps * x0;
        }
        else
        {
            h = 4.0 * sqrteps;
        }

        // use a 5-point stencil for the central difference quotient
        double numerator = f(x0 - 2.0 * h) - 8.0 * f(x0 - h) + 8.0 * f(x0 + h) - f(x0 + 2.0 * h);
        double denominator = 12.0 * h;

        return numerator / denominator;
    }

    template <>
    double derivative<2u, deriv::TwoSided>(const std::function<double (const double &)> & f, const double & x0)
    {
        static const double eps = std::numeric_limits<double>::epsilon(), sqrteps = std::sqrt(eps), sqrt2eps = std::sqrt(sqrteps);

        double h;
        if (x0 > 4.0 * sqrt2eps)
        {
            h = sqrt2eps * x0;
        }
        else
        {
            h = 4.0 * sqrt2eps;
        }

        // use a 5-point stencil for the central difference quotient
        double numerator = -f(x0 - 2.0 * h) + 16.0 * f(x0 - h) - 30.0 * f(x0) + 16.0 * f(x0 + h) - f(x0 + 2.0 * h);
        double denominator = 12.0 * h * h;

        return numerator / denominator;
    }
}
