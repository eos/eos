/*
 * Copyright (c) 2022 Méril Reboud
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

#ifndef EOS_GUARD_EOS_MATHS_LAGRANGE_POLYNOMIAL_HH
#define EOS_GUARD_EOS_MATHS_LAGRANGE_POLYNOMIAL_HH 1

#include <eos/maths/complex.hh>

#include <array>

namespace eos
{
    /*
     * Representation of a complex Lagrange polynomial.
     * The polynomial L satisfies L(x_values) = y_values.
     */
    template <unsigned order_>
    class LagrangePolynomial
    {
        private:
            const std::array<complex<double>, order_ + 1> _x_values;

        public:
            LagrangePolynomial(std::array<complex<double>, order_ + 1 > && x_values) :
                _x_values(x_values)
            {
            }

            LagrangePolynomial(const LagrangePolynomial &) = default;
            LagrangePolynomial(LagrangePolynomial &&) = default;
            ~LagrangePolynomial() = default;
            LagrangePolynomial & operator= (LagrangePolynomial &&) = default;
            LagrangePolynomial & operator= (const LagrangePolynomial &) = default;

            // Evaluate the Lagrange polynomial
            complex<double> operator() (const std::array<complex<double>, order_ + 1 > & y_values, const complex<double> & z) const
            {
                complex<double> result = 0;
                complex<double> partial_result;

                for (unsigned i = 0; i <= order_; i++)
                {
                    partial_result = 1;
                    for (unsigned j = 0; j <= order_; j++)
                    {
                        if (i != j)
                        {
                            partial_result *= (z - _x_values[j])/(_x_values[i] - _x_values[j]);
                        }
                    }
                    result += partial_result * y_values[i];
                }

                return result;
            }
    };
}

#endif
