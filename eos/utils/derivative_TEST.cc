/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/utils/derivative.hh>

#include <cmath>
#include <limits>

using namespace test;
using namespace eos;

class DerivativeTest :
    public TestCase
{
    public:
        DerivativeTest() :
            TestCase("derivative_test")
        {
        }

        static double f1(const double & x)
        {
            return 6.0 * x * (1.0 - x);
        }

        static double f2(const double & x)
        {
            return f1(x) / (1.0 - x);
        }

        static double f3(const double & x)
        {
            return std::exp(-x);
        }

        static double f4(const double & x)
        {
            return std::log(x);
        }

        static double f5(const double & x)
        {
            return std::cos(x);
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            // two sided, at x = 0
            {
                double f1_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f1), 0.0);
                double f1_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f1), 0.0);
                TEST_CHECK_NEARLY_EQUAL(f1_d1,   6.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f1_d2, -12.0, eps);

                double f2_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f2), 0.0);
                double f2_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f2), 0.0);
                TEST_CHECK_NEARLY_EQUAL(f2_d1,   6.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f2_d2,   0.0, eps);

                double f3_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f3), 0.0);
                double f3_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f3), 0.0);
                TEST_CHECK_NEARLY_EQUAL(f3_d1,  -1.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f3_d2,  +1.0, eps);

                // skip f4, since log(x) diverges for x -> 0

                double f5_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f5), 0.0);
                double f5_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f5), 0.0);
                TEST_CHECK_NEARLY_EQUAL(f5_d1,   0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f5_d2,  -1.0, eps);
            }

            // two sided, at x != 0
            {
                double f1_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f1), 0.5);
                double f1_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f1), 0.5);
                TEST_CHECK_NEARLY_EQUAL(f1_d1,   0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f1_d2, -12.0, eps);

                double f2_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f2), 0.5);
                double f2_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f2), 0.5);
                TEST_CHECK_NEARLY_EQUAL(f2_d1, 6.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f2_d2, 0.0, eps);

                double f3_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f3), 0.5);
                double f3_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f3), 0.5);
                TEST_CHECK_NEARLY_EQUAL(f3_d1, -std::exp(-0.5), eps);
                TEST_CHECK_NEARLY_EQUAL(f3_d2, +std::exp(-0.5), eps);

                double f4_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f4), 0.5);
                double f4_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f4), 0.5);
                TEST_CHECK_NEARLY_EQUAL(f4_d1, +2.0, eps);
                TEST_CHECK_NEARLY_EQUAL(f4_d2, -4.0, eps);

                double f5_d1 = derivative<1u, deriv::TwoSided>(std::function<double (const double &)>(&f5), 0.5);
                double f5_d2 = derivative<2u, deriv::TwoSided>(std::function<double (const double &)>(&f5), 0.5);
                TEST_CHECK_NEARLY_EQUAL(f5_d1, -std::sin(0.5), eps);
                TEST_CHECK_NEARLY_EQUAL(f5_d2, -std::cos(0.5), eps);
            }
        }
} derivative_test;
