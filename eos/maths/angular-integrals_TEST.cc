/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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
#include <eos/maths/angular-integrals.hh>

using namespace test;
using namespace eos;

class AngularIntegralTest :
    public TestCase
{
    public:
        AngularIntegralTest() :
            TestCase("angular_integral_test")
        {
        }

        virtual void run() const
        {

            constexpr double eps = 1e-5;

            {
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(1,0,1,0,1,0),       0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(1,0,1,0,2,0),   0.36515, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(1,0,2,0,3,0),  -0.29277, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(2,0,1,0,3,0),  -0.29277, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(1,1,2,1,3,-2), -0.30861, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(2,1,1,1,3,-2), -0.30861, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(2,2,1,1,3,-3),  0.37796, eps);
                TEST_CHECK_NEARLY_EQUAL(wigner_3j(2,2,4,1,6,-3), -0.13993, eps);

                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(0,0,1,0),    0.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(2,1,1,1),    0.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(0,0,0,0),    2.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(2,0,2,0),    0.40000, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(1,1,1,1),    1.33333, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(4,3,4,3), 1120.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(1,1,2,0),    0.19635, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(5,1,2,2),    0.13806, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(1,-1,2,0),  -0.09818, eps);
                TEST_CHECK_NEARLY_EQUAL(two_legendre_integral(3,1,2,-2),   0.03682, eps);

                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(0,0,0,0,0,0),    2.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(2,0,2,0,0,0),    0.40000, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(1,1,1,1,0,0),    1.33333, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(4,3,4,3,0,0), 1120.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(2,0,2,0,4,0),    0.11429, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(2,1,2,1,4,0),   -0.45714, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(1,0,2,0,3,0),    0.17143, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(1,1,2,1,3,0),   -0.34286, eps);
                TEST_CHECK_NEARLY_EQUAL(three_legendre_integral(1,0,2,1,2,0),   -0.29452, eps);
            }
        }
} angular_integral_test;
