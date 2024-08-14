/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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
#include <eos/maths/szego-polynomial.hh>

#include <cmath>
#include <gsl/gsl_matrix.h>

using namespace test;
using namespace eos;

class SzegoPolynomialTest :
    public TestCase
{
    public:
        SzegoPolynomialTest() :
            TestCase("szego_polynomial_test")
        {
        }

        virtual void run() const
        {
            // Test the evaluation
            {
                const auto p = SzegoPolynomial<5u>::FlatMeasure(2.47895); // norm of the measure

                {
                    const auto [p0, p1, p2, p3, p4, p5] = p(-0.1);

                    TEST_CHECK_RELATIVE_ERROR(p0, +0.6351351032391984, 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p1, -0.847745,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p2, +1.54489,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p3, -2.82388,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p4, +5.166081,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p5, -9.464102,           1.0e-5);
                }

                {
                    const auto [p0, p1, p2, p3, p4, p5] = p(+0.0);

                    TEST_CHECK_RELATIVE_ERROR(p0, +0.6351351032391984, 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p1, -0.749503,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p2, +1.304458,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p3, -2.237009,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p4, +3.834085,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p5, -6.577731,           1.0e-5);
                }

                {
                    const auto [p0, p1, p2, p3, p4, p5] = p(+0.1);

                    TEST_CHECK_RELATIVE_ERROR(p0, +0.6351351032391984, 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p1, -0.651261,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p2, +1.09668,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p3, -1.76188,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p4, +2.82970,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p5, -4.54691,            1.0e-5);
                }

                {
                    const auto [p0, p1, p2, p3, p4, p5] = p(complex<double>(0.4, 0.916515139));

                    TEST_CHECK_RELATIVE_ERROR(real(p0), +0.6351351032391984,                    1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p1, complex<double>(-0.35653430,  0.90040363),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p2, complex<double>(-0.70241504, -0.85660005),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p3, complex<double>( 1.06106731, -0.16955667),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p4, complex<double>(-0.31896585,  0.83869890),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p5, complex<double>(-0.40550250, -0.45614671),  1.0e-5);
                }
            }

            // Test the coefficients
            {
                const auto p = SzegoPolynomial<5u>::FlatMeasure(2.47895); // norm of the measure

                gsl_matrix * coefficient_matrix = p.coefficient_matrix();

                TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(coefficient_matrix, 0, 0),  0.6351351032391984,  1.0e-5);
                TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(coefficient_matrix, 1, 2), -2.24105,  1.0e-5);
                TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(coefficient_matrix, 1, 5),  24.1447,  1.0e-5);
                TEST_CHECK_RELATIVE_ERROR(gsl_matrix_get(coefficient_matrix, 3, 4), -12.6392,  1.0e-5);
            }

            // Test the derivatives
            {
                const auto p = SzegoPolynomial<5u>::FlatMeasure(2.47895); // norm of the measure

                {
                    const auto [p0, p1, p2, p3, p4, p5] = p.derivatives(-0.1);

                    TEST_CHECK_NEARLY_EQUAL(std::real(p0),    0.0,                1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(std::real(p1), +0.982422,           1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(std::real(p2), -2.56765,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(std::real(p3), +6.48297,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(std::real(p4), -15.2203,            1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(std::real(p5), +34.0797,            1.0e-5);
                }

                {
                    const auto [p0, p1, p2, p3, p4, p5] = p.derivatives(complex<double>(0.4, 0.916515139));

                    TEST_CHECK_NEARLY_EQUAL(p0,                      0.0,                  1.0e-5);
                    TEST_CHECK_NEARLY_EQUAL(p1.imag(),               0.0,                  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(p1.real(),             0.982422,             1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p2, complex<double>(-0.934628,  2.99338),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p3, complex<double>(-4.83802,  -4.15038),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p4, complex<double>( 9.43519,  -4.10475),  1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR_C(p5, complex<double>(-0.881431,  14.1859),  1.0e-5);
                }
            }

        }
} szego_polynomial_test;
