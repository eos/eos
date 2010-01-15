/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/complex.hh>

#include <cmath>
#include <limits>

using namespace test;
using namespace wf;

class ComplexTest :
    public TestCase
{
    public:
        ComplexTest() :
            TestCase("complex_test")
        {
        }

        virtual void run() const
        {
            Complex<double> z1 = Complex<double>::Cartesian(3.0, 4.0);
            Complex<double> z2 = Complex<double>::Cartesian(3.0, -4.0);

            TEST_CHECK(std::numeric_limits<double>::epsilon() >= std::abs(z1.absolute() - 5.0));
            TEST_CHECK(std::numeric_limits<double>::epsilon() >= std::abs(z2.absolute() - 5.0));

            TEST_CHECK(std::numeric_limits<double>::epsilon() >= std::abs((z1 + z2).absolute() - 6.0));

            TEST_CHECK(std::numeric_limits<double>::epsilon() >= std::abs((z1 * z2).absolute() - 25.0));
        }
} model_test;
