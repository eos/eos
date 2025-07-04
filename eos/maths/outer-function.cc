/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024-2025 Florian Herren
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

#include <eos/maths/integrate.hh>
#include <eos/maths/outer-function.hh>
#include <eos/utils/exception.hh>

namespace eos
{

    complex<double> outer(const std::function<complex<double> (const complex<double> &)> & f, complex<double> z, double relative_precision)
    {
        if (std::abs(z) >= 1.0)
            throw InternalError("Trying to evaluate outer function outside of unit disk. This is not yet supported.");

        std::function<complex<double>(const double &)> integrand = [&] (const double & t)
        {
            double x = std::abs(f(std::exp(complex<double>(0, t))));
            if (std::abs(x) < 1e-16)
            {
                throw InternalError("Trying to compute outer function of function with zero on unit circle");
            }
            else if (!isfinite(x))
            {
                throw InternalError("Trying to compute outer function of function with pole on unit circle");
            }
            return (std::exp(complex<double>(0, t)) + z) / (std::exp(complex<double>(0, t)) - z) * std::log(x);
        };

        return std::exp(integrate(integrand, 0.0, 2 * M_PI, cubature::Config().epsrel(relative_precision)) / 2.0 / M_PI );
    }

}
