/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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
#include <eos/utils/integrate.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace eos;

class IntegrateTest :
    public TestCase
{
    public:
        IntegrateTest() :
            TestCase("integrate_test")
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

        virtual void run() const
        {
            double q1 = integrate(std::function<double (const double &)>(&f1), 16, 0.0, 1.0), i1 = 1.0;
            std::cout << "\\int_0.0^1.0 f1(x) dx = " << q1 << ", eps = " << std::abs(i1 - q1) / q1 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i1 - q1) / i1 < 0.01);

            double q2 = integrate(std::function<double (const double &)>(&f2), 16, 0.00, 0.999999), i2 = 3.0;
            std::cout << "\\int_0.0^1.0 f2(x) dx = " << q2 << ", eps = " << std::abs(i2 - q2) / q2 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i2 - q2) / i2 < 0.01);

            double q3 = integrate(std::function<double (const double &)>(&f3), 16, 0.00, 10.0), i3 = 1.0 - std::exp(-10.0);
            std::cout << "\\int_0.0^10.0 f3(x) dx = " << q3 << ", eps = " << std::abs(i3 - q3) / q3 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i3 - q3) / i3 < 0.01);

            double q4 = integrate(std::function<double (const double &)>(&f4), 16, 1.0, std::exp(1)), i4 = 1.0;
            std::cout << "\\int_0.0^10.0 f4(x) dx = " << q4 << ", eps = " << std::abs(i4 - q4) / q4 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i4 - q4) / i4 < 0.01);
        }
} model_test;
