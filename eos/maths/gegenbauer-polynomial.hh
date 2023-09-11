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

#ifndef EOS_GUARD_EOS_MATHS_GEGENBAUER_POLYNOMIALS_HH
#define EOS_GUARD_EOS_MATHS_GEGENBAUER_POLYNOMIALS_HH 1

#include <eos/maths/power-of.hh>

#include <vector>

namespace eos
{
    class GegenbauerPolynomial
    {
        private:
            const unsigned int _order; // order of the Gegenbauer polynomial

            const double _alpha; // parameter of the Gegenbauer polynomial with alpha > -0.5

            const double _r; // remainder of (order / 2)

            const unsigned int _half_order; // rounded down half of the order: floor(order / 2)

            std::vector<double> _calculate_coefficients() const;

            std::vector<double> _coefficients;

        public:
            GegenbauerPolynomial(const unsigned int & order, const double & alpha) :
                _order(order),
                _alpha(alpha),
                _r((_order % 2 == 0) ? 0.0 : 1.0),
                _half_order(floor(static_cast<double>(order) / 2.0)),
                _coefficients(_calculate_coefficients())
            {
            }

            double evaluate(const double & z) const;
    };
}

#endif
