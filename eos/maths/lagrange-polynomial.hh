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
#include <eos/maths/szego-polynomial.hh>

#include <array>
#include <numeric>

namespace eos
{
    template <unsigned order_, unsigned n_>
    class LagrangePolynomialDerivative;

    template <unsigned order_, unsigned n_>
    class LagrangePolynomialCoefficients;

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
                _x_values(std::move(x_values))
            {
            }

            LagrangePolynomial(const std::array<complex<double>, order_ + 1 > & x_values) :
                _x_values(x_values)
            {
            }

            LagrangePolynomial(const LagrangePolynomial &) = default;
            LagrangePolynomial(LagrangePolynomial &&) = default;
            ~LagrangePolynomial() = default;
            LagrangePolynomial & operator= (LagrangePolynomial &&) = default;
            LagrangePolynomial & operator= (const LagrangePolynomial &) = default;

            // Evaluate the Lagrange polynomial
            complex<double> constexpr operator() (const std::array<complex<double>, order_ + 1> & y_values, const complex<double> & z) const
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

            // Returns the coefficients of the polynomial in the monomial basis (z^n)_n
            std::array<complex<double>, order_ + 1> get_coefficients(const std::array<complex<double>, order_ + 1> & y_values) const
            {
                LagrangePolynomialCoefficients<order_, order_> lpc;

                return lpc.fill_coefficients(_x_values, y_values);
            }
    };


    template <unsigned order_>
    class LagrangePolynomialDerivative<order_, 0>
    {
        public:
            complex<double> constexpr evaluate(const std::array<complex<double>, order_ + 1> & x_values,
                                               const std::array<complex<double>, order_ + 1> & y_values,
                                               const complex<double> & z) const
            {
                LagrangePolynomial<order_> polynomial(x_values);

                return polynomial(y_values, z);
            }
    };

    template <unsigned n_>
    class LagrangePolynomialDerivative<0, n_>
    {
        public:
            complex<double> constexpr evaluate(const std::array<complex<double>, 1> &,
                                               const std::array<complex<double>, 1> &,
                                               const complex<double> &) const
            {
                return 0.;
            }
    };

    template <>
    class LagrangePolynomialDerivative<0, 0>
    {
        public:
            complex<double> constexpr evaluate(const std::array<complex<double>, 1> &,
                                               const std::array<complex<double>, 1> & y_values,
                                               const complex<double> &) const
            {
                return y_values[0];
            }
    };

    template <unsigned order_, unsigned n_>
    class LagrangePolynomialDerivative
    {
        private:
            std::array<complex<double>, order_ + 1> _make_derivative(const std::array<complex<double>, order_ + 1> & x_values,
                                                                     const std::array<complex<double>, order_ + 1> & y_values,
                                                                     const complex<double> & z) const
            {
                std::array<complex<double>, order_ + 1> derivatives;

                for (unsigned i = 0 ; i < order_ + 1 ; ++i)
                {
                    std::array<complex<double>, order_> new_y_values;

                    for (unsigned k = 0 ; k < order_ + 1 ; ++k)
                    {
                        if (k < i)
                        {
                            new_y_values[k] = y_values[k] / (x_values[k] - x_values[i]);
                        }
                        else if (k > i)
                        {
                            new_y_values[k - 1] = y_values[k] / (x_values[k] - x_values[i]);
                        }
                    }

                    LagrangePolynomialDerivative<order_ - 1, n_ - 1> d;

                    std::array<complex<double>, order_> new_x_values;
                    std::copy(x_values.begin(), x_values.begin() + i, new_x_values.begin());
                    std::copy(x_values.begin() + i + 1, x_values.end(), new_x_values.begin() + i);

                    derivatives[i] = d.evaluate(new_x_values, new_y_values, z);
                }

                return derivatives;
            }


        public:
            complex<double> constexpr evaluate(const std::array<complex<double>, order_ + 1> & x_values,
                                               const std::array<complex<double>, order_ + 1> & y_values,
                                               const complex<double> & z) const
            {
                std::array<complex<double>, order_ + 1> derivatives = _make_derivative(x_values, y_values, z);

                return std::accumulate(derivatives.begin(), derivatives.end(), complex<double>(0., 0.));
            }
    };

    template <unsigned order_>
    class LagrangePolynomialCoefficients<order_, 0>
    {
        public:
            unsigned factorial() const
            {
                return 1;
            }

            // Coefficients of the Lagrange polynomial in the monomial basis
            std::array<complex<double>, order_ + 1> fill_coefficients(const std::array<complex<double>, order_ + 1> & x_values,
                                                                      const std::array<complex<double>, order_ + 1> & y_values) const
            {
                std::array<complex<double>, order_ + 1> coefficients;

                LagrangePolynomialDerivative<order_, 0> dL;

                coefficients[0] = dL.evaluate(x_values, y_values, 0.);

                return coefficients;
            }
    };

    template <unsigned order_, unsigned n_>
    class LagrangePolynomialCoefficients
    {
        public:
            unsigned factorial() const
            {
                LagrangePolynomialCoefficients<order_, n_ - 1> lpc;

                return n_ * lpc.factorial();
            }

            // Coefficients of the Lagrange polynomial in the monomial basis
            std::array<complex<double>, order_ + 1> fill_coefficients(const std::array<complex<double>, order_ + 1> & x_values,
                                                                      const std::array<complex<double>, order_ + 1> & y_values) const
            {
                LagrangePolynomialCoefficients<order_, n_ - 1> lpc;

                std::array<complex<double>, order_ + 1> coefficients = lpc.fill_coefficients(x_values, y_values);

                LagrangePolynomialDerivative<order_, n_> dL;

                coefficients[n_] = 1.0 / this->factorial() * dL.evaluate(x_values, y_values, 0.);

                return coefficients;
            }
    };
}

#endif
