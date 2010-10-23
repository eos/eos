/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/apply.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

struct TestPointerToMemberFunction
{
    double _x, _y;
    void f(const double & x, const double & y)
    {
        _x = x; _y = y;
    }
};

class ApplyTest :
    public TestCase
{
    public:
        ApplyTest() :
            TestCase("apply_test")
        {
        }

        static double test_2(const double & x)
        {
            return x;
        }

        static double test_3_function()
        {
            return std::exp(1.0);
        }

        virtual void run() const
        {
            TestPointerToMemberFunction test_1;

            apply(&TestPointerToMemberFunction::f, std::make_tuple(&test_1, 1.0, 2.0));
            TEST_CHECK_EQUAL(test_1._x, 1.0);
            TEST_CHECK_EQUAL(test_1._y, 2.0);

            double result;

            result = apply(&ApplyTest::test_2, std::make_tuple(M_PI));
            TEST_CHECK_EQUAL(result, M_PI);

            std::function<double ()> test_3(&ApplyTest::test_3_function);
            result = apply(test_3, std::make_tuple());
            TEST_CHECK_EQUAL(result, std::exp(1.0));
        }
} apply_test;
