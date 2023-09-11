/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <test/test.hh>
#include <eos/maths/gegenbauer-polynomial.hh>

#include <cmath>
#include <array>

using namespace test;
using namespace eos;

class GegenbauerPolynomialTest :
    public TestCase
{
    public:
        GegenbauerPolynomialTest() :
            TestCase("gegenbauer_polynomial_test")
        {
        }

        virtual void run() const
        {
            {
                const GegenbauerPolynomial gp0(2, 0.0);
                const GegenbauerPolynomial gp1(3, 0.5);
                const GegenbauerPolynomial gp2(4, 1.5);

                {
                    const double z[5] = {-1.0, -0.3, 0.0, 0.7, 1.0};
                    double g[3];
                    const double v[3][5] = {{  1.0, -0.82,      -1.0,   -0.02,       1.0 },
                                            { -1.0,  0.3825,     0.0,   -0.1925,     1.0 },
                                            { 15.0, -0.1685625,  1.875, -1.5335625, 15.0 }};

                    for (int i = 0 ; i < 5 ; i++)
                    {
                        g[0] = gp0.evaluate(z[i]);
                        g[1] = gp1.evaluate(z[i]);
                        g[2] = gp2.evaluate(z[i]);

                        TEST_CHECK_NEARLY_EQUAL(g[0], v[0][i], 1.0e-5);
                        TEST_CHECK_NEARLY_EQUAL(g[1], v[1][i], 1.0e-5);
                        TEST_CHECK_NEARLY_EQUAL(g[2], v[2][i], 1.0e-5);
                    }
                }
            }
        }
} gegenbauer_polynomial_test;
