/*
 * Copyright (c) 2023 Viktor Kuschke
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

#include <eos/maths/power-of.hh>
#include <eos/maths/gegenbauer-polynomial.hh>

#include <cmath>
#include <vector>

namespace eos
{
    std::vector<double> GegenbauerPolynomial::_calculate_coefficients() const
    {
        const double epsilon = 1e-8;
        std::vector<double> result(_half_order + 1);

        for (unsigned int i = 0 ; i <= _half_order ; i++)
        {
            const double entry = pow(2.0, 2.0 * i + _r) * pow(-1.0, 1.0 * (_half_order - i)) *
                tgamma(_alpha + 1.0 * (_half_order + i) + _r) / (tgamma(1.0 * (_half_order - i + 1)) * tgamma(2.0 * i + _r + 1.0)); // denominator: (half_order - i)! * (2i + r)!;

            result.at(i) = entry;
        }

        double gamma_alpha = (std::abs(_alpha) < epsilon) ? 1.0 : tgamma(_alpha);
        for (unsigned int i = 0 ; i <= _half_order ; i++)
        {
            result[i] /= gamma_alpha;
        }

        return result;
    }

    double GegenbauerPolynomial::evaluate(const double & z) const
    {
        double result = 0.0;
        double x = (1 - _r) * 1.0 + _r * z;

        for (double entry : _coefficients)
        {
            result += entry * x;
            x *= z * z;
        }

        return result;
    }
}
