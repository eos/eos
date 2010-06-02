/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/charm-loops.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace wf;

class OneLoopTest :
    public TestCase
{
    public:
        OneLoopTest() :
            TestCase("one_loop_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christoph Bobeth's result from May 2010 */

            /* One-Loop */
            {
                static const double mu = 4.2, s = 1.0, m_c = 1.4, m_b = 4.8, eps = 0.00001;
                TEST_CHECK_NEARLY_EQUAL(+1.57192, real(CharmLoops::h(mu, s)), eps);
                TEST_CHECK_NEARLY_EQUAL(+1.39626, imag(CharmLoops::h(mu, s)), eps);

                TEST_CHECK_NEARLY_EQUAL(+0.58013, CharmLoops::h(mu, s, m_c), eps);

                TEST_CHECK_NEARLY_EQUAL(-0.55926, CharmLoops::h(mu, s, m_b), eps);
            }
        }
} one_loop_test;

class TwoLoopTest :
    public TestCase
{
    public:
        TwoLoopTest() :
            TestCase("two_loop_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christoph Bobeth's result from May 2010 */

            /* Two-Loop, massless */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, eps = 0.0000001;
                TEST_CHECK_NEARLY_EQUAL(- 0.8832611, real(CharmLoops::F17(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.6937322, imag(CharmLoops::F17(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 5.2995666, real(CharmLoops::F27(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 4.1623936, imag(CharmLoops::F27(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 3.3632062, real(CharmLoops::F19(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 6.9078480, imag(CharmLoops::F19(mu, s, m_b)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 3.4455298, real(CharmLoops::F29(mu, s, m_b)), eps);
                TEST_CHECK_NEARLY_EQUAL(+24.6919276, imag(CharmLoops::F29(mu, s, m_b)), eps);

                // F87, F89?
            }

            /* Two-Loop, massive */
            {
                static const double mu = 4.2, s = 6.0, m_b = 4.6, m_c = 1.2, eps = 0.0000001;
                TEST_CHECK_NEARLY_EQUAL(- 0.73093991, real(CharmLoops::F17(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.17771334, imag(CharmLoops::F17(mu, s, m_b, m_c)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 4.38563254, real(CharmLoops::F27(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 1.06627403, imag(CharmLoops::F27(mu, s, m_b, m_c)), eps);

                TEST_CHECK_NEARLY_EQUAL(-34.40870331, real(CharmLoops::F19(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(- 0.25864665, imag(CharmLoops::F19(mu, s, m_b, m_c)), eps);

                TEST_CHECK_NEARLY_EQUAL(+ 6.27364439, real(CharmLoops::F29(mu, s, m_b, m_c)), eps);
                TEST_CHECK_NEARLY_EQUAL(+ 1.55195807, imag(CharmLoops::F29(mu, s, m_b, m_c)), eps);

                // F87, F89?
            }
        }
} two_loop_test;
