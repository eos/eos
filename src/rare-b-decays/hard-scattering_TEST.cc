/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/hard-scattering.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace wf;

class HardScatteringTest :
    public TestCase
{
    public:
        HardScatteringTest() :
            TestCase("hard_scattering_test")
        {
        }

        virtual void run() const
        {
            /* One-Loop */
            {
                static const double s = 1.0, m_c = 1.4, m_B = 5.279, eps = 0.0000001;

                TEST_CHECK_NEARLY_EQUAL(+0.7005873, real(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.2097160, imag(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0317353, real(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.5061933, imag(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.2769282, real(HardScattering::I1(s, 0.9, m_c, m_B)), eps);
                TEST_CHECK_NEARLY_EQUAL(+0.0000000, imag(HardScattering::I1(s, 0.9, m_c, m_B)), eps);
            }
        }
} hard_scattering_test;
