/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk, 2014 Stephan Jahn
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
#include <eos/utils/integrate-impl.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace eos;
using namespace eos::integration;
using std::array;

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
            double q1 = integrate(std::function<double (const double &)>(&f1), 0.0, 1.0, 16), i1 = 1.0;
            std::cout << "\\int_0.0^1.0 f1(x) dx = " << q1 << ", eps = " << std::abs(i1 - q1) / q1 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i1 - q1) / i1 < 0.01);

            double q2 = integrate(std::function<double (const double &)>(&f2), 0.00, 0.999999, 16), i2 = 3.0;
            std::cout << "\\int_0.0^1.0 f2(x) dx = " << q2 << ", eps = " << std::abs(i2 - q2) / q2 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i2 - q2) / i2 < 0.01);

            double q3 = integrate(std::function<double (const double &)>(&f3), 0.00, 10.0, 16), i3 = 1.0 - std::exp(-10.0);
            std::cout << "\\int_0.0^10.0 f3(x) dx = " << q3 << ", eps = " << std::abs(i3 - q3) / q3 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i3 - q3) / i3 < 0.01);

            double q4 = integrate(std::function<double (const double &)>(&f4), 1.0, std::exp(1), 16), i4 = 1.0;
            std::cout << "\\int_0.0^10.0 f4(x) dx = " << q4 << ", eps = " << std::abs(i4 - q4) / q4 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i4 - q4) / i4 < 0.01);
        }
} model_test;

class IntegrateDefaultNTest :
    public TestCase
{
    public:
        IntegrateDefaultNTest() :
            TestCase("integrate_default_n_test")
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

        static complex<double> fc(const double & x)
        {
            return complex<double>(1.,2.) + x;
        }

        static array<double, 5> farray(const double & x)
        {
            return array<double, 5>{1.,2.*x,3.,4.,5.};
        }

        virtual void run() const
        {
            TEST_CHECK_EQUAL(get_n(), 64);
            set_n(16);
            TEST_CHECK_EQUAL(get_n(), 16);

            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, set_n(3));
            TEST_CHECK_EQUAL(get_n(), 16);

            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, set_n(35));
            TEST_CHECK_EQUAL(get_n(), 16);

            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, integrate(std::function<double (const double &)>(&f1), 0.0, 1.0, 8));
            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, integrate(std::function<double (const double &)>(&f1), 0.0, 1.0, 33));

            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, integrate(std::function<complex<double> (const double &)>(&fc), 0.0, 1.0, 2));
            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, integrate(std::function<complex<double> (const double &)>(&fc), 0.0, 1.0, 73));

            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, integrate(std::function<std::array<double, 5> (const double &)>(&farray), 0.0, 1.0, 3));
            TEST_CHECK_THROWS(InvalidNumberOfEvaluations, integrate(std::function<std::array<double, 5> (const double &)>(&farray), 0.0, 1.0, 93));

            double q1 = integrate(std::function<double (const double &)>(&f1), 0.0, 1.0), i1 = 1.0;
            std::cout << "\\int_0.0^1.0 f1(x) dx = " << q1 << ", eps = " << std::abs(i1 - q1) / q1 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i1 - q1) / i1 < 0.01);

            double q2 = integrate(std::function<double (const double &)>(&f2), 0.00, 0.999999), i2 = 3.0;
            std::cout << "\\int_0.0^1.0 f2(x) dx = " << q2 << ", eps = " << std::abs(i2 - q2) / q2 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i2 - q2) / i2 < 0.01);

            double q3 = integrate(std::function<double (const double &)>(&f3), 0.00, 10.0), i3 = 1.0 - std::exp(-10.0);
            std::cout << "\\int_0.0^10.0 f3(x) dx = " << q3 << ", eps = " << std::abs(i3 - q3) / q3 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i3 - q3) / i3 < 0.01);

            double q4 = integrate(std::function<double (const double &)>(&f4), 1.0, std::exp(1)), i4 = 1.0;
            std::cout << "\\int_0.0^10.0 f4(x) dx = " << q4 << ", eps = " << std::abs(i4 - q4) / q4 << " over 16 points" << std::endl;
            TEST_CHECK(std::abs(i4 - q4) / i4 < 0.01);
        }
} integrate_default_n_test;
