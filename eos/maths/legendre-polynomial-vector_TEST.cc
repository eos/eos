/*
 * Copyright (c) 2024 Florian Herren
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
#include <eos/maths/legendre-polynomial-vector.hh>

using namespace test;
using namespace eos;

class LegendrePolynomialVectorTest :
    public TestCase
{
    public:
        LegendrePolynomialVectorTest() :
            TestCase("legendre_polynomial_vector_test")
        {
        }

        virtual void run() const
        {
            // test case
            {
                // Evaluate P
                {
                    LegendrePVector<6u> P;

                    auto resP = P(0.5);

                    TEST_CHECK_RELATIVE_ERROR(resP[0],  1.0,            1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resP[1],  0.5,            1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resP[2],  -0.125,         1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resP[3],  -0.4375,        1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resP[4],  -0.2890625,     1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resP[5],  0.08984375,     1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resP[6],  0.3232421875,   1.0e-10);
                }

                // Evaluate P
                {
                    LegendreReQVector<6u> Q;

                    auto resQ15 = Q(1.5);

                    TEST_CHECK_RELATIVE_ERROR(resQ15[0],    0.8047189562170503,     1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resQ15[1],    0.20707843432557507,    1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resQ15[2],    0.06356699912401897,    1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resQ15[3],    0.02086520825966291,    1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resQ15[4],    0.007095922338601302,   1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resQ15[5],    0.0024668237064868373,  1.0e-10);
                    TEST_CHECK_RELATIVE_ERROR(resQ15[6],    0.0008704965773399678,  1.0e-10);
                }
            }
        }
} legendre_polynomial_vector_test;
