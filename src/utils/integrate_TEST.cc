/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/integrate.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace wf;

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

        virtual void run() const
        {
            double q1 = integrate(std::tr1::function<double (const double &)>(&f1), 250, 0.0, 1.0), i1 = 1.0;
            std::cout << "\\int_0.0^1.0 f1(x) dx = " << q1 << " over 250 points" << std::endl;
            TEST_CHECK(std::abs(i1 - q1) / i1 < 0.04);

            double q2 = integrate(std::tr1::function<double (const double &)>(&f2), 1000, 0.00, 0.999), i2 = 3.0;
            std::cout << "\\int_0.0^1.0 f2(x) dx = " << q2 << " over 250 points" << std::endl;
            TEST_CHECK(std::abs(i2 - q2) / i2 < 0.04);
        }
} model_test;
