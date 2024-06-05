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

#include <test/test.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/lagrange-polynomial.hh>

#include <iostream>

using namespace test;
using namespace eos;

class LagrangePolynomialTest :
    public TestCase
{
    public:
        LagrangePolynomialTest() :
            TestCase("lagrange_polynomial_test")
        {
        }

        virtual void run() const
        {
            // test case
            {

                // Evaluation
                {
                    LagrangePolynomial<3u> L{
                        { -1.0,  0.0, 2.0, complex<double>(0.0, 1.0) },
                    };

                    const std::array<complex<double>, 4u> y_values{ 1.0, -2.0, complex<double>(1.0, 1.0), 7.0 };

                    TEST_CHECK_RELATIVE_ERROR(real(L(y_values, -1.0)),                         1.0,                       1.0e-5);
                    TEST_CHECK_NEARLY_EQUAL(  imag(L(y_values, -1.0)),                         0.,                        1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(L(y_values,  2.0),                             complex<double>(1.0, 1.0), 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(real(L(y_values, complex<double>(0.0,  1.0))),   7.0,                       1.0e-5);
                    TEST_CHECK_NEARLY_EQUAL(  imag(L(y_values, complex<double>(0.0,  1.0))),   0.,                        1.0e-5);

                    TEST_CHECK_RELATIVE_ERROR_C(L(y_values,  1.0),                       complex<double>(-3.13333333, -6.4), 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(L(y_values, complex<double>(0.0, -1.0)), complex<double>(-12.8666666,  6.4), 1.0e-5);
                }

                // 0th order derivative
                {
                    LagrangePolynomialDerivative<3u, 0u> L;

                    const std::array<complex<double>, 4u> x_values{ -1.0,  0.0, 2.0, complex<double>(0.0, 1.0) };
                    const std::array<complex<double>, 4u> y_values{  1.0, -2.0, complex<double>(1.0, 1.0), 7.0 };

                    complex<double> Lat1 = L.evaluate(x_values, y_values, 1.0);
                    TEST_CHECK_RELATIVE_ERROR_C(Lat1, complex<double>(-3.13333333, -6.4), 1.0e-5);
                }

                // 1st order derivative
                {
                    LagrangePolynomialDerivative<3u, 1u> dL;

                    const std::array<complex<double>, 4u> x_values{ -1.0,  0.0, 2.0, complex<double>(0.0, 1.0) };
                    const std::array<complex<double>, 4u> y_values{  1.0, -2.0, complex<double>(1.0, 1.0), 7.0 };

                    complex<double> dLat0 = dL.evaluate(x_values, y_values, 0.0);
                    TEST_CHECK_RELATIVE_ERROR_C(dLat0, complex<double>(-2.63333333, -6.5666666), 1.0e-5);
                    complex<double> dLat1pi = dL.evaluate(x_values, y_values, complex<double>(1.0, 1.0));
                    TEST_CHECK_RELATIVE_ERROR_C(dLat1pi, complex<double>(-14.5666666, -7.7), 1.0e-5);
                }

                // 2nd order derivative
                {
                    LagrangePolynomialDerivative<3u, 2u> dL;

                    const std::array<complex<double>, 4u> x_values{ -1.0,  0.0, 2.0, complex<double>(0.0, 1.0) };
                    const std::array<complex<double>, 4u> y_values{  1.0, -2.0, complex<double>(1.0, 1.0), 7.0 };

                    complex<double> d2Lat0 = dL.evaluate(x_values, y_values, 0.0);
                    TEST_CHECK_RELATIVE_ERROR_C(d2Lat0, complex<double>(1.8666667, -6.4), 1.0e-5);
                    complex<double> d2Lat1pi = dL.evaluate(x_values, y_values, complex<double>(1.0, 1.0));
                    TEST_CHECK_RELATIVE_ERROR_C(d2Lat1pi, complex<double>(-14.9333333, 17.2), 1.0e-5);
                }

                // get coefficients
                {
                    LagrangePolynomial<3u> L{
                        { -1.0,  0.0, 2.0, complex<double>(0.0, 1.0) },
                    };
                    const std::array<complex<double>, 4u> y_values{  1.0, -2.0, complex<double>(1.0, 1.0), 7.0 };

                    auto coefficients = L.get_coefficients(y_values);

                    TEST_CHECK_RELATIVE_ERROR(real(coefficients[0]), -2.,                                 1.0e-5);
                    TEST_CHECK_NEARLY_EQUAL(  imag(coefficients[0]),  0.,                                 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(coefficients[1], complex<double>(-2.6333333, -6.5666667), 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(coefficients[2], complex<double>( 0.9333333, -3.2),       1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(coefficients[3], complex<double>( 0.5666667,  3.3666667), 1.0e-5);
                }
            }
        }
} lagrange_polynomial_test;
