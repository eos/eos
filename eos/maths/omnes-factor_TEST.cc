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
#include <eos/maths/omnes-factor.hh>
#include <eos/maths/omnes-factor-impl.hh>

#include <cmath>
#include <array>

using namespace test;
using namespace eos;

double test_phase(const double & s)
{
    if (s < 5.0)
    {
        return -std::atan(std::sqrt(s - 4.0) / (s - 16.0));
    }
    else
    {
        return M_PI / 2.0 + std::atan((s - 16.0) / std::sqrt(s - 4.0));
    }
}

/////////
/// TESTS
/////////
class OmnesFactorTest :
    public TestCase
{
    public:
        OmnesFactorTest() :
            TestCase("OmnesFactor tests")
        {
        }

        virtual void run() const
        {

            constexpr double eps = 1e-5;

            {
                std::array<double, 4> intervals = {4.0, 10.0, 25.0, 50.0};
                OmnesFactor<50, 4> O(intervals, test_phase, 1.0);
                OmnesFactor<50, 4> O2(intervals, test_phase, O.get_weights());

                TEST_CHECK_NEARLY_EQUAL(O(-25.0),       0.36072,    eps);
                TEST_CHECK_NEARLY_EQUAL(O(-12.5),       0.51385,    eps);
                TEST_CHECK_NEARLY_EQUAL(O(-1.5),        0.84313,    eps);
                TEST_CHECK_NEARLY_EQUAL(O(1.0),         1.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(O(2.5),         1.13632,    eps);
                TEST_CHECK_NEARLY_EQUAL(O(3.9),         1.34760,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(4.1)),    1.40556,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(8.0)),    2.02906,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(12.0)),   3.41542,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(16.0)),   4.83013,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(16.1)),   4.80814,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(25.01)),  1.65518,    eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O(30.0)),   1.12298,    eps);
                TEST_CHECK_NEARLY_EQUAL(O2(1.0),        1.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(abs(O2(16.1)),  4.80814,    eps);
            }
    }
} omnes_factor_test;
