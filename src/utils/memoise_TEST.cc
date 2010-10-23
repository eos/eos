/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/memoise.hh>

#include <complex>

using namespace test;
using namespace eos;

class MemoiseTest :
    public TestCase
{
    public:
        MemoiseTest() :
            TestCase("memoise_test")
        {
        }

        static double f1(const double & x, const double & y)
        {
            return x / y;
        }

        static std::complex<double> f2(const double & x, const double & y)
        {
            return std::complex<double>(x, y);
        }

        virtual void run() const
        {
            /* f1 */
            {
                TEST_CHECK_EQUAL(0.5, f1(1.0, 2.0));
                TEST_CHECK_EQUAL(2.0, f1(2.0, 1.0));

                // First round of memoisation
                TEST_CHECK_EQUAL(0.5, memoise(f1, 1.0, 2.0));
                TEST_CHECK_EQUAL(2.0, memoise(f1, 2.0, 1.0));

                // Second round of memoisation
                TEST_CHECK_EQUAL(0.5, memoise(f1, 1.0, 2.0));
                TEST_CHECK_EQUAL(2.0, memoise(f1, 2.0, 1.0));
            }

            /* f2 */
            {
                TEST_CHECK_EQUAL(1.0, real(f2(1.0, 2.0)));
                TEST_CHECK_EQUAL(2.0, imag(f2(1.0, 2.0)));

                // First round of memoisation
                TEST_CHECK_EQUAL(1.0, real(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(2.0, imag(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(2.0, real(memoise(f2, 2.0, 1.0)));
                TEST_CHECK_EQUAL(1.0, imag(memoise(f2, 2.0, 1.0)));

                // Second round of memoisation
                TEST_CHECK_EQUAL(1.0, real(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(2.0, imag(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(2.0, real(memoise(f2, 2.0, 1.0)));
                TEST_CHECK_EQUAL(1.0, imag(memoise(f2, 2.0, 1.0)));
            }
        }
} memoise_test;
