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
#include <eos/maths/complex.hh>
#include <eos/maths/outer-function.hh>

using namespace test;
using namespace eos;

// Trivial case: Phi(z) = 1
complex<double> test_func_1(const complex<double> & z)
{
    return z;
}

// No zeros in the unit disk: Phi(z) = 2 + z
complex<double> test_func_2(const complex<double> & z)
{
    return 2.0 + z;
}

// Pole and zero inside disk: Phi(z) = (2 - z) / (2 + z)
complex<double> test_func_3(const complex<double> & z)
{
    return (1.0 - 2.0 * z) / (1.0 + 2.0 * z);
}

// Pole at 0 and zero on circle: Phi(z) = 1 - z
complex<double> test_func_4(const complex<double> & z)
{
    return (1.0 - z) / z;
}

// Zero at 0 and pole on circle: Phi(z) = 1 / (1 - z)
complex<double> test_func_5(const complex<double> & z)
{
    return z / (1.0 - z);
}

/////////
/// TESTS
/////////
class OuterFunctionTest :
    public TestCase
{
    public:
        OuterFunctionTest() :
            TestCase("OuterFunction tests")
        {
        }

        virtual void run() const
        {

            constexpr double eps = 1e-5;

            {
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_1, 0.0, 1024).real(),                       1.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_1, 0.0, 1024).imag(),                       0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_1, 0.5, 1024).real(),                       1.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_1, 0.5, 1024).imag(),                       0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_1, complex<double>(0.1,0.2), 1024).real(),  1.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_1, complex<double>(0.1,0.2), 1024).imag(),  0.0,        eps);

                TEST_CHECK_NEARLY_EQUAL(outer(test_func_2, 0.0, 1024).real(),                       2.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_2, 0.0, 1024).imag(),                       0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_2, 0.5, 1024).real(),                       2.5,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_2, 0.5, 1024).imag(),                       0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_2, complex<double>(0.1,0.2), 1024).real(),  2.1,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_2, complex<double>(0.1,0.2), 1024).imag(),  0.2,        eps);

                TEST_CHECK_NEARLY_EQUAL(outer(test_func_3, 0.0, 1024).real(),                       1.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_3, 0.0, 1024).imag(),                       0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_3, 0.5, 1024).real(),                       0.6,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_3, 0.5, 1024).imag(),                       0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_3, complex<double>(0.1,0.2), 1024).real(),  0.887640,   eps);
                TEST_CHECK_NEARLY_EQUAL(outer(test_func_3, complex<double>(0.1,0.2), 1024).imag(), -0.179775,   eps);

                TEST_CHECK_THROWS(InternalError, outer(test_func_1, 2.0, 1024));
                TEST_CHECK_THROWS(InternalError, outer(test_func_4, 0.0, 1024));
                TEST_CHECK_THROWS(InternalError, outer(test_func_5, 0.0, 1024));
            }
    }
} outer_function_test;
