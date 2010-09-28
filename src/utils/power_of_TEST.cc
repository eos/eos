/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/power_of.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace wf;

class PowerOfTest :
    public TestCase
{
    public:
        PowerOfTest() :
            TestCase("power_of_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-14;

            TEST_CHECK_NEARLY_EQUAL(1.0,      power_of<0>(1.2), eps);
            TEST_CHECK_NEARLY_EQUAL(1.2,      power_of<1>(1.2), eps);
            TEST_CHECK_NEARLY_EQUAL(1.44,     power_of<2>(1.2), eps);
            TEST_CHECK_NEARLY_EQUAL(1.728,    power_of<3>(1.2), eps);
            TEST_CHECK_NEARLY_EQUAL(2.0736,   power_of<4>(1.2), eps);
            TEST_CHECK_NEARLY_EQUAL(2.48832,  power_of<5>(1.2), eps);
            TEST_CHECK_NEARLY_EQUAL(2.985984, power_of<6>(1.2), eps);

            TEST_CHECK_NEARLY_EQUAL(1.0,      power_of<0>(0.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.4,      power_of<1>(0.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.16,     power_of<2>(0.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.064,    power_of<3>(0.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0256,   power_of<4>(0.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.01024,  power_of<5>(0.4), eps);
            TEST_CHECK_NEARLY_EQUAL(0.004096, power_of<6>(0.4), eps);
        }
} power_of_test;
