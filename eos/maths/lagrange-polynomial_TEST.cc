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
                LagrangePolynomial<3u> L{
                    { -1.0,  0.0, 2.0, complex<double>(0.0, 1.0) },
                };

                {
                    const std::array<complex<double>, 4u> y_values{ 1.0, -2.0, complex<double>(1.0, 1.0), 7.0 };

                    TEST_CHECK_RELATIVE_ERROR(L(y_values, -1.0),                      1.0,                       1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(L(y_values,  2.0),                      complex<double>(1.0, 1.0), 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(L(y_values, complex<double>(0.0, 1.0)), 7.0,                       1.0e-5);

                    TEST_CHECK_RELATIVE_ERROR(L(y_values,  1.0),                      complex<double>(-3.13333333,-6.4), 1.0e-5);
                    TEST_CHECK_RELATIVE_ERROR(L(y_values, complex<double>(0.0,-1.0)), complex<double>(-12.8666666, 6.4), 1.0e-5);
                }
            }
        }
} lagrange_polynomial_test;
