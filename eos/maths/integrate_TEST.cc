/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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
#include <eos/maths/integrate-impl.hh>

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

        std::function<std::array<double, 4> (const std::array<double, 1> &)> f6 = [&] (const std::array<double, 1> & x)
        {
            return std::array<double, 4>{f4(x[0]), f4(x[0]), f4(x[0]), f4(x[0])};
        };

        std::function<std::array<double, 4> (const double &)> f7 = [&] (const double & x)
        {
            return std::array<double, 4>{std::log(x), 2 * std::log(x), 3 * std::log(x), 4 * std::log(x)};
        };

        virtual void run() const
        {
            constexpr double eps = 0.01;

            double q1 = integrate1D(std::function<double (const double &)>(&f1), 16, 0.0, 1.0), i1 = 1.0;
            std::cout << "\\int_0.0^1.0 f1(x) dx = " << q1 << ", eps = " << std::abs(i1 - q1) / q1 << " over 16 points" << std::endl;
            TEST_CHECK_RELATIVE_ERROR(i1, q1, eps);

            double q2 = integrate1D(std::function<double (const double &)>(&f2), 16, 0.00, 0.999999), i2 = 3.0;
            std::cout << "\\int_0.0^1.0 f2(x) dx = " << q2 << ", eps = " << std::abs(i2 - q2) / q2 << " over 16 points" << std::endl;
            TEST_CHECK_RELATIVE_ERROR(i2, q2, eps);

            double q3 = integrate1D(std::function<double (const double &)>(&f3), 16, 0.00, 10.0), i3 = 1.0 - std::exp(-10.0);
            std::cout << "\\int_0.0^10.0 f3(x) dx = " << q3 << ", eps = " << std::abs(i3 - q3) / q3 << " over 16 points" << std::endl;
            TEST_CHECK_RELATIVE_ERROR(i3, q3, eps);

            auto f4obj = std::function<double (const double &)>(&f4);
            double q4 = integrate1D(f4obj, 16, 1.0, std::exp(1)), i4 = 1.0;
            std::cout << "\\int_0.0^exp(1) f4(x) dx = " << q4 << ", eps = " << std::abs(i4 - q4) / q4 << " over 16 points" << std::endl;
            TEST_CHECK_RELATIVE_ERROR(i4, q4, eps);

            auto config_QNG = GSL::QNG::Config().epsrel(eps);
            q4 = integrate<GSL::QNG>(f4obj, 1.0, std::exp(1), config_QNG);
            std::cout << "\\int_0.0^exp(1) f4(x) dx = " << q4 << ", eps = " << std::abs(i4 - q4) / q4 << " with QNG" << std::endl;
            TEST_CHECK_RELATIVE_ERROR(i4, q4, eps);

            auto config_QAGS = GSL::QAGS::Config().epsrel(1e-12);
            q4 = integrate<GSL::QAGS>(f4obj, 1.0, std::exp(1), config_QAGS);
            std::cout << "\\int_0.0^exp(1) f4(x) dx = " << q4 << ", eps = " << std::abs(i4 - q4) / q4 << " with QAGS" << std::endl;
            TEST_CHECK_RELATIVE_ERROR(i4, q4, eps);

            auto config_cubature = cubature::Config().epsrel(eps);
            auto f4lam = [](const std::array<double, 1> &args) -> double {
                return f4(args[0]);
            };
            constexpr std::array<double, 1> a_4 { 1.0 };
            static const std::array<double, 1> b_4 { std::exp(1) }; // std::exp isn't constexpr (yet)
            q4 = integrate(cubature::fdd_s_v<1>(f4lam), a_4, b_4, config_cubature);
            TEST_CHECK_RELATIVE_ERROR(i4, q4, eps);

            // Morokoff test function
            constexpr size_t dim = 4;
            constexpr std::array<double, dim> a_5 { 0, 0, 0, 0 };
            constexpr std::array<double, dim> b_5 { 1, 1, 1, 1 };
            auto f5lam = [&dim](const std::array<double, dim> &args) -> double {
                double p = 1.0 / dim;
                double prod = pow(1 + p, dim);
                unsigned int i;
                for (i = 0; i < dim; i++)
                    prod *= pow(args[i], p);
                return prod;
            };
            auto q5 = integrate(cubature::fdd_s_v<dim>(f5lam), a_5, b_5, config_cubature);
            TEST_CHECK_RELATIVE_ERROR(q5, 1.0, eps);

            std::array<double, 4> q6 = integrate(f6, std::array<double, 1>{1.0}, std::array<double, 1>{std::exp(1)}, config_cubature);
            TEST_CHECK_RELATIVE_ERROR(i4, q6[0], eps);
            TEST_CHECK_RELATIVE_ERROR(i4, q6[1], eps);

            std::array<double, 4> q7 = integrate(f7, 1.0, std::exp(1), config_cubature);
            TEST_CHECK_RELATIVE_ERROR(    i4, q7[0], eps);
            TEST_CHECK_RELATIVE_ERROR(2 * i4, q7[1], eps);
            TEST_CHECK_RELATIVE_ERROR(3 * i4, q7[2], eps);
            TEST_CHECK_RELATIVE_ERROR(4 * i4, q7[3], eps);

            double q8 = integrate(cubature::fdd_s_s(f4), 1.0, std::exp(1), config_cubature);
            TEST_CHECK_RELATIVE_ERROR(i4, q8, eps);
        }
} model_test;
