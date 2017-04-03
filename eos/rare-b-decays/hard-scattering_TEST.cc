/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2017 Danny van Dyk
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
#include <eos/rare-b-decays/hard-scattering.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

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
                static const double m_c = 1.4, m_B = 5.279;

                /* s = -4.0 */
                {
                    static const double s = -4.0, eps = 1.0e-5;

                    TEST_CHECK_NEARLY_EQUAL(+0.582643,  real(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.209591,  real(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(+0.161063,  real(HardScattering::I1(s, 0.9, m_c, m_B)), eps);

                    TEST_CHECK_NEARLY_EQUAL(-1.00951,   imag(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(-1.03736,   imag(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(+0.0000000, imag(HardScattering::I1(s, 0.9, m_c, m_B)), eps);
                }

                /* s = +1.0 */
                {
                    static const double s = +1.0, eps = 1.0e-7;

                    TEST_CHECK_NEARLY_EQUAL(+0.7005873, real(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(+0.0317353, real(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(-0.2769282, real(HardScattering::I1(s, 0.9, m_c, m_B)), eps);

                    TEST_CHECK_NEARLY_EQUAL(-1.2097160, imag(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(-1.5061933, imag(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(+0.0000000, imag(HardScattering::I1(s, 0.9, m_c, m_B)), eps);
                }

                /* s = +8.0 */
                {
                    static const double s = +8.0, eps = 1.0e-4;

                    TEST_CHECK_NEARLY_EQUAL(+1.62683,   real(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(+1.73189,   real(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(+1.91415,   real(HardScattering::I1(s, 0.9, m_c, m_B)), eps);

                    TEST_CHECK_NEARLY_EQUAL(-1.46251,   imag(HardScattering::I1(s, 0.1, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(-2.06062,   imag(HardScattering::I1(s, 0.5, m_c, m_B)), eps);
                    TEST_CHECK_NEARLY_EQUAL(-4.45749,   imag(HardScattering::I1(s, 0.9, m_c, m_B)), eps);
                }
            }
        }
} hard_scattering_test;
