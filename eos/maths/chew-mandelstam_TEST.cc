/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2026 Simon Mutke
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/maths/chew-mandelstam.hh>
#include <eos/maths/power-of.hh>

#include <test/test.hh>

#include <cmath>
#include <functional>
#include <string>
#include <utility>
#include <vector>

using namespace test;
using namespace eos;

class ChewMandelstamSWaveTest : public TestCase
{
    public:
        ChewMandelstamSWaveTest() :
            TestCase("chew_mandelstam_s_wave_test")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-5;

            // electron mass, cf. eos/scattering/ee-to-ccbar, EEChannel
            constexpr double m  = 5.10999e-4;
            constexpr double m1 = 1.86;
            constexpr double m2 = 2.01;

            // real s, above threshold
            {
                const complex<double> result = chew_mandelstam::s_wave(20.0, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11496166, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01989437, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(2.0, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.10038034, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01989436, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(0.5, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.09160146, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01989435, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(20.0, m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11496166, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01989437, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(2.0, m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.10038034, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01989436, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(0.5, m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.09160146, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01989435, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(40.0, m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.01074883, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01573072, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(20.0, m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00349196, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00996453, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(15.0, m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0000194850, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0007801256, eps);
            }

            // real s, below threshold: the function is real-valued
            {
                const complex<double> result = chew_mandelstam::s_wave(-20.0, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11496167, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(-20.0, m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11496167, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(-2.0, m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.01319441, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // at threshold
            {
                const complex<double> result = chew_mandelstam::s_wave(4 * m * m, m);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(4 * m * m, m, m);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave((m1 + m2) * (m1 + m2), m1, m2);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // genuinely complex s, e.g. as used to move onto another Riemann sheet
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(20.0, 5.0), m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11515361, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01834302, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(2.0, -1.0), m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.10108689, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.01695827, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(-5.0, 0.3), m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.10619424, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00037950, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(20.0, 5.0), m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11515361, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.01834302, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(2.0, -1.0), m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.10108689, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.01695827, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(-5.0, 0.3), m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.10619424, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00037950, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(20.0, 5.0), m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0070441584, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0084444293, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(2.0, -1.0), m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0120721769, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.0003153895, eps);
            }
            {
                const complex<double> result = chew_mandelstam::s_wave(complex<double>(-5.0, 0.3), m1, m2);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0139095167, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000669477, eps);
            }
        }
} chew_mandelstam_s_wave_test;

class ChewMandelstamPWaveTest : public TestCase
{
    public:
        ChewMandelstamPWaveTest() :
            TestCase("chew_mandelstam_p_wave_test")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-5;

            // D^+ mass and effective momentum scale, cf. eos/scattering/ee-to-ccbar, PWavePPChannel
            constexpr double mD = 1.86965;
            constexpr double m1 = 1.86;
            constexpr double m2 = 2.01;
            constexpr double q0 = 0.5;

            // real s, below threshold: the function is real-valued
            {
                const complex<double> result = chew_mandelstam::p_wave(-20.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0125291, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(2.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0081324, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(-20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0125291, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(2.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0081324, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(-20.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.012448605, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(2.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.008264793, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // at threshold
            {
                const complex<double> result = chew_mandelstam::p_wave(4 * mD * mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(4 * mD * mD, mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave((m1 + m2) * (m1 + m2), m1, m2, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // real s, above threshold
            {
                const complex<double> result = chew_mandelstam::p_wave(20.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00024738, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0093576, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00024738, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0093576, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(20.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.0007870605, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0083085871, eps);
            }

            // effective channel mass, cf. eos/scattering/ee-to-ccbar, EffChannel
            {
                constexpr double      meff   = 0.1349768;
                const complex<double> result = chew_mandelstam::p_wave(20.0, meff, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0172797, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0189092, eps);
            }
            {
                constexpr double      meff   = 0.1349768;
                const complex<double> result = chew_mandelstam::p_wave(20.0, meff, meff, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0172797, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0189092, eps);
            }

            // right at threshold, where the Taylor-expanded branch is taken
            {
                const double          mp2    = power_of<2>(2.0 * mD);
                const complex<double> result = chew_mandelstam::p_wave(mp2, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, 1e-6);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, 1e-6);
            }
            {
                const double          mp2    = power_of<2>(2.0 * mD);
                const complex<double> result = chew_mandelstam::p_wave(mp2, mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, 1e-6);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, 1e-6);
            }

            // at the critical points delta and s1plus, s1minus
            {
                const double          delta  = 4.0 * (mD * mD - q0 * q0);
                const complex<double> result = chew_mandelstam::p_wave(delta, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.001973137240, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const double a = m1 * m1 + m2 * m2 - 2.0 * q0 * q0;
                const double b = std::sqrt(a * a - power_of<2>((m1 + m2) * (m1 - m2)));

                const complex<double> result_plus  = chew_mandelstam::p_wave(a + b, m1, m2, q0);
                const complex<double> result_minus = chew_mandelstam::p_wave(a - b, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result_plus), -0.001923828789, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result_plus), 0.0, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_minus), -0.008806210696, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result_minus), 0.0, eps);
            }

            // genuinely complex s, e.g. as used to move onto another Riemann sheet
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(20.0, 5.0), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0032251, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0081452, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(2.0, -1.0), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0081419, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.0003102, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(-5.0, 0.3), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0099493, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000658867, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(20.0, 5.0), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0032251, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0081452, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(2.0, -1.0), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0081419, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.0003102, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(-5.0, 0.3), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0099493, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000658867, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(20.0, 5.0), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00304153, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00740049, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(2.0, -1.0), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00827295, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.000288635, eps);
            }
            {
                const complex<double> result = chew_mandelstam::p_wave(complex<double>(-5.0, 0.3), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0099739, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000625791, eps);
            }
        }
} chew_mandelstam_p_wave_test;

class ChewMandelstamDWaveTest : public TestCase
{
    public:
        ChewMandelstamDWaveTest() :
            TestCase("chew_mandelstam_d_wave_test")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-5;

            // D^+ mass and effective momentum scale, cf. eos/scattering/ee-to-ccbar, PWavePPChannel
            constexpr double mD = 1.86965;
            constexpr double m1 = 1.86;
            constexpr double m2 = 2.01;
            constexpr double q0 = 0.5;

            // real s, below threshold: the function is real-valued
            {
                const complex<double> result = chew_mandelstam::d_wave(-20.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.009340026, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(2.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0054647, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(-20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.009340026, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(2.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0054647, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(-20.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.009333018, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(2.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0056257, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // at the critical points s1plus, s1minus, s2plus, s2minus, which are complex for l = 2
            {
                const complex<double> u1 = 3.0 / 2.0 * complex<double>(-1.0, std::sqrt(3.0));
                const complex<double> u2 = -3.0 / 2.0 * complex<double>(1.0, std::sqrt(3.0));
                const complex<double> a1 = m1 * m1 + m2 * m2 + 2.0 * q0 * q0 * u1;
                const complex<double> b1 = std::sqrt(a1 * a1 - power_of<2>((m1 + m2) * (m1 - m2)));
                const complex<double> a2 = m1 * m1 + m2 * m2 + 2.0 * q0 * q0 * u2;
                const complex<double> b2 = std::sqrt(a2 * a2 - power_of<2>((m1 + m2) * (m1 - m2)));

                const complex<double> result_1plus  = chew_mandelstam::d_wave(a1 + b1, m1, m2, q0);
                const complex<double> result_1minus = chew_mandelstam::d_wave(a1 - b1, m1, m2, q0);
                const complex<double> result_2plus  = chew_mandelstam::d_wave(a2 + b2, m1, m2, q0);
                const complex<double> result_2minus = chew_mandelstam::d_wave(a2 - b2, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result_1plus), -0.001641809996, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_1plus), 0.001614695243, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_1minus), -0.006088914216, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_1minus), -1.042989707e-6, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_2plus), -0.001641809996, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_2plus), -0.001614695243, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_2minus), -0.006088914216, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_2minus), 1.042989707e-6, eps);
            }

            // at threshold
            {
                const complex<double> result = chew_mandelstam::d_wave(4 * mD * mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(4 * mD * mD, mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave((m1 + m2) * (m1 + m2), m1, m2, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // real s, above threshold
            {
                const complex<double> result = chew_mandelstam::d_wave(20.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00336492, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00624622, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00336492, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00624622, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(20.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00344569, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00509587, eps);
            }

            // genuinely complex s, e.g. as used to move onto another Riemann sheet
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(20.0, 5.0), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0003150954, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0062365592, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(2.0, -1.0), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0054716, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.000261292, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(-5.0, 0.3), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00702802, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000577958, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(20.0, 5.0), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0003150954, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0062365592, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(2.0, -1.0), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0054716, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.000261292, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(-5.0, 0.3), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00702802, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000577958, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(20.0, 5.0), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0003740028, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0055794287, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(2.0, -1.0), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00563175, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.000245385, eps);
            }
            {
                const complex<double> result = chew_mandelstam::d_wave(complex<double>(-5.0, 0.3), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0071068, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000551838, eps);
            }
        }
} chew_mandelstam_d_wave_test;

class ChewMandelstamFWaveTest : public TestCase
{
    public:
        ChewMandelstamFWaveTest() :
            TestCase("chew_mandelstam_f_wave_test")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-5;

            // D^+ mass and effective momentum scale, cf. eos/scattering/ee-to-ccbar, PWavePPChannel
            constexpr double mD = 1.86965;
            constexpr double m1 = 1.86;
            constexpr double m2 = 2.01;
            constexpr double q0 = 0.5;

            // real s, below threshold: the function is real-valued
            {
                const complex<double> result = chew_mandelstam::f_wave(-20.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00698812, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(2.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00371294, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(-20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00698812, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(2.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00371294, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(-20.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00701956, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(2.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00386655, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // at the critical points s1plus, s1minus, which are real for l = 3, as well as s2plus, s2minus, s3plus, s3minus, which are complex for l = 3
            {
                const double          acbr = std::cbrt((75.0 * std::sqrt(5.0) - 151.0) / 2.0);
                const double          u1   = -2.0 + acbr - 11.0 / acbr;
                const complex<double> u2   = -2.0 + acbr * complex<double>(-1, std::sqrt(3.0)) / 2.0 + 5.5 / acbr * complex<double>(1, std::sqrt(3.0));
                const complex<double> u3   = -2.0 - acbr * complex<double>(1, std::sqrt(3.0)) / 2.0 - 5.5 / acbr * complex<double>(-1, std::sqrt(3.0));
                const double          a1   = m1 * m1 + m2 * m2 + 2.0 * q0 * q0 * u1;
                const double          b1   = std::sqrt(a1 * a1 - power_of<2>((m1 + m2) * (m1 - m2)));
                const complex<double> a2   = m1 * m1 + m2 * m2 + 2.0 * q0 * q0 * u2;
                const complex<double> b2   = std::sqrt(a2 * a2 - power_of<2>((m1 + m2) * (m1 - m2)));
                const complex<double> a3   = m1 * m1 + m2 * m2 + 2.0 * q0 * q0 * u3;
                const complex<double> b3   = std::sqrt(a3 * a3 - power_of<2>((m1 + m2) * (m1 - m2)));

                const complex<double> result_1plus  = chew_mandelstam::f_wave(a1 + b1, m1, m2, q0);
                const complex<double> result_1minus = chew_mandelstam::f_wave(a1 - b1, m1, m2, q0);
                const complex<double> result_2plus  = chew_mandelstam::f_wave(a2 + b2, m1, m2, q0);
                const complex<double> result_2minus = chew_mandelstam::f_wave(a2 - b2, m1, m2, q0);
                const complex<double> result_3plus  = chew_mandelstam::f_wave(a3 + b3, m1, m2, q0);
                const complex<double> result_3minus = chew_mandelstam::f_wave(a3 - b3, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result_1plus), -0.002053710801, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result_1plus), 0.0, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_1minus), -0.004241819684, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result_1minus), 0.0, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_2plus), -0.00107689, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_2plus), 0.00250824, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_2minus), -0.004244756, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_2minus), -1.5560258e-6, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_3plus), -0.00107689, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_3plus), -0.00250824, eps);

                TEST_CHECK_RELATIVE_ERROR(real(result_3minus), -0.004244756, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result_3minus), 1.5560258e-6, eps);
            }

            // at threshold
            {
                const complex<double> result = chew_mandelstam::f_wave(4 * mD * mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(4 * mD * mD, mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave((m1 + m2) * (m1 + m2), m1, m2, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // real s, above threshold
            {
                const complex<double> result = chew_mandelstam::f_wave(20.0, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00394589, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00255427, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00394589, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00255427, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(20.0, m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00343814, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00172857, eps);
            }

            // genuinely complex s, e.g. as used to move onto another Riemann sheet
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(20.0, 5.0), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00110454, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00406132, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(2.0, -1.0), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00371747, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.000209013, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(-5.0, 0.3), mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00499524, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000485, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(20.0, 5.0), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00110454, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00406132, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(2.0, -1.0), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00371747, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.000209013, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(-5.0, 0.3), mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00499524, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000485, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(20.0, 5.0), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.000861406, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.00359757, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(2.0, -1.0), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00387061, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), -0.0001984, eps);
            }
            {
                const complex<double> result = chew_mandelstam::f_wave(complex<double>(-5.0, 0.3), m1, m2, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.00509165, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0000466051, eps);
            }
        }
} chew_mandelstam_f_wave_test;

class ChewMandelstamMasslessTest : public TestCase
{
    public:
        ChewMandelstamMasslessTest() :
            TestCase("chew_mandelstam_massless_test")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-6;

            static constexpr double m  = 1.86;
            static constexpr double q0 = 0.5;

            // reference values from the dispersion integral of rho(s) = (s - m^2) / (16 pi s), subtracted at threshold
            const std::vector<std::pair<complex<double>, complex<double>>> references{
                {          -2.0,                        -0.00788656433562 },
                {           1.0,                        -0.00531368162968 },
                { { 20.0, 5.0 },   { -0.00926101976069, 0.0146846040815 } },
                { { 2.0, -1.0 }, { -0.00426185662505, -0.00147813852605 } },
            };
            for (const auto & [s, reference] : references)
            {
                for (const auto & result : { chew_mandelstam::s_wave(s, m, 0.0), chew_mandelstam::s_wave(s, 0.0, m) })
                {
                    TEST_CHECK_NEARLY_EQUAL(real(result), real(reference), eps * std::abs(reference));
                    TEST_CHECK_NEARLY_EQUAL(imag(result), imag(reference), eps * std::abs(reference));
                }
            }

            // Im s_wave = rho above threshold
            for (const double s : { 5.0, 20.0 })
            {
                TEST_CHECK_RELATIVE_ERROR(imag(chew_mandelstam::s_wave(s, m, 0.0)), (s - m * m) / (16.0 * M_PI * s), eps);
            }

            // the limit m2 -> 0 is approached as O(m2 log m2)
            static constexpr double small     = 1e-5;
            static constexpr double eps_limit = 1e-3;

            const auto check_limit = [&](const complex<double> & a, const complex<double> & b)
            {
                TEST_CHECK_NEARLY_EQUAL(real(a), real(b), eps_limit * std::abs(b));
                TEST_CHECK_NEARLY_EQUAL(imag(a), imag(b), eps_limit * std::abs(b));
            };

            for (const complex<double> s : { complex<double>(-2.0), complex<double>(5.0), complex<double>(20.0, 5.0), complex<double>(2.0, -1.0) })
            {
                check_limit(chew_mandelstam::p_wave(s, m, 0.0, q0), chew_mandelstam::p_wave(s, m, small, q0));
                check_limit(chew_mandelstam::d_wave(s, m, 0.0, q0), chew_mandelstam::d_wave(s, m, small, q0));
                check_limit(chew_mandelstam::f_wave(s, m, 0.0, q0), chew_mandelstam::f_wave(s, m, small, q0));
            }

            // at a pole of the Blatt-Weisskopf form factor, which uses the derivative of s_wave
            {
                const double          a = m * m - 2.0 * q0 * q0;
                const complex<double> s = a + std::sqrt(complex<double>(a * a - power_of<4>(m)));

                check_limit(chew_mandelstam::p_wave(s, m, 0.0, q0), chew_mandelstam::p_wave(s, m, small, q0));
            }
        }
} chew_mandelstam_massless_test;

class ChewMandelstamCoincidentPolesTest : public TestCase
{
    public:
        ChewMandelstamCoincidentPolesTest() :
            TestCase("chew_mandelstam_coincident_poles_test")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps   = 1e-6;
            constexpr double delta = 1e-6;

            // real root of the F-wave Blatt-Weisskopf denominator u^3 + 6 u^2 + 45 u + 225
            const double acbr = std::cbrt((75.0 * std::sqrt(5.0) - 151.0) / 2.0);
            const double u1   = -2.0 + acbr - 11.0 / acbr;

            // the value at coincident poles must match the average of its neighbours
            const auto check = [&](const std::function<complex<double>(const double &)> & f, const double & q0)
            {
                const complex<double> result    = f(q0);
                const complex<double> reference = 0.5 * (f(q0 * (1.0 - delta)) + f(q0 * (1.0 + delta)));

                TEST_CHECK_NEARLY_EQUAL(real(result), real(reference), eps * std::abs(reference));
                TEST_CHECK_NEARLY_EQUAL(imag(result), imag(reference), eps * std::abs(reference));
            };

            for (const complex<double> s : { complex<double>(20.0), complex<double>(5.0), complex<double>(-2.0), complex<double>(20.0, 5.0) })
            {
                for (const auto & [m1, m2] : {
                         std::pair{  1.0,  2.0 },
                         std::pair{ 1.86, 2.01 }
                })
                {
                    const auto p_wave = [&](const double & q0) { return chew_mandelstam::p_wave(s, m1, m2, q0); };
                    const auto f_wave = [&](const double & q0) { return chew_mandelstam::f_wave(s, m1, m2, q0); };

                    // poles coincide for q0^2 u = -m1^2 or -m2^2, at u = -1 for the P wave and u = u1 for the F wave
                    check(p_wave, m1);
                    check(p_wave, m2);
                    check(f_wave, m1 / std::sqrt(-u1));
                    check(f_wave, m2 / std::sqrt(-u1));
                }

                // equal masses, where the coincident poles sit at s = 0
                const double m      = 1.86;
                const auto   f_wave = [&](const double & q0) { return chew_mandelstam::f_wave(s, m, q0); };
                check(f_wave, m / std::sqrt(-u1));
            }
        }
} chew_mandelstam_coincident_poles_test;

class ChewMandelstamDiscontinuityTest : public TestCase
{
    public:
        ChewMandelstamDiscontinuityTest() :
            TestCase("chew_mandelstam_discontinuity_test")
        {
        }

        virtual void
        run() const
        {
            using Function = std::function<complex<double>(const complex<double> &)>;

            struct Channel
            {
                    std::string name;
                    Function    cm;
                    Function    disc;
                    double      m1, m2;
            };

            static constexpr double me = 5.10999e-4;
            static constexpr double mD = 1.86965;
            static constexpr double m1 = 1.86;
            static constexpr double m2 = 2.01;
            static constexpr double q0 = 0.5;

            using namespace chew_mandelstam;

            // clang-format off
            const std::vector<Channel> channels{
                { "S, m",      [](const complex<double> & s) { return s_wave(s, me); },          [](const complex<double> & s) { return s_wave_disc(s, me); },          me, me },
                { "S, m, m",   [](const complex<double> & s) { return s_wave(s, me, me); },      [](const complex<double> & s) { return s_wave_disc(s, me, me); },      me, me },
                { "S, m1, m2", [](const complex<double> & s) { return s_wave(s, m1, m2); },      [](const complex<double> & s) { return s_wave_disc(s, m1, m2); },      m1, m2 },
                { "P, m",      [](const complex<double> & s) { return p_wave(s, mD, q0); },      [](const complex<double> & s) { return p_wave_disc(s, mD, q0); },      mD, mD },
                { "P, m, m",   [](const complex<double> & s) { return p_wave(s, mD, mD, q0); },  [](const complex<double> & s) { return p_wave_disc(s, mD, mD, q0); },  mD, mD },
                { "P, m1, m2", [](const complex<double> & s) { return p_wave(s, m1, m2, q0); },  [](const complex<double> & s) { return p_wave_disc(s, m1, m2, q0); },  m1, m2 },
                { "D, m",      [](const complex<double> & s) { return d_wave(s, mD, q0); },      [](const complex<double> & s) { return d_wave_disc(s, mD, q0); },      mD, mD },
                { "D, m, m",   [](const complex<double> & s) { return d_wave(s, mD, mD, q0); },  [](const complex<double> & s) { return d_wave_disc(s, mD, mD, q0); },  mD, mD },
                { "D, m1, m2", [](const complex<double> & s) { return d_wave(s, m1, m2, q0); },  [](const complex<double> & s) { return d_wave_disc(s, m1, m2, q0); },  m1, m2 },
                { "F, m",      [](const complex<double> & s) { return f_wave(s, mD, q0); },      [](const complex<double> & s) { return f_wave_disc(s, mD, q0); },      mD, mD },
                { "F, m, m",   [](const complex<double> & s) { return f_wave(s, mD, mD, q0); },  [](const complex<double> & s) { return f_wave_disc(s, mD, mD, q0); },  mD, mD },
                { "F, m1, m2", [](const complex<double> & s) { return f_wave(s, m1, m2, q0); },  [](const complex<double> & s) { return f_wave_disc(s, m1, m2, q0); },  m1, m2 },
            };
            // clang-format on

            constexpr double eps = 1e-6;

            // small imaginary part relative to the real part of s
            const auto ie = [](const double & x) { return complex<double>(0.0, 1e-11 * std::abs(x)); };

            const auto check_close = [](const complex<double> & a, const complex<double> & b, const double & eps)
            {
                TEST_CHECK_NEARLY_EQUAL(real(a), real(b), eps * std::abs(b));
                TEST_CHECK_NEARLY_EQUAL(imag(a), imag(b), eps * std::abs(b));
            };

            for (const auto & c : channels)
            {
                const double s_plus  = power_of<2>(c.m1 + c.m2);
                const double s_minus = power_of<2>(c.m1 - c.m2);

                // above threshold, the discontinuity of the Chew-Mandelstam function across the real axis
                for (const double x : { s_plus + 0.5, s_plus + 2.0, s_plus + 10.0 })
                {
                    const complex<double> disc = c.disc(x);

                    TEST_CHECK_MSG(imag(disc) > 0.0, c.name);
                    TEST_CHECK_RELATIVE_ERROR(imag(c.cm(x)), 0.5 * imag(disc), eps);
                    check_close(c.disc(x + ie(x)), disc, eps);
                    check_close(c.cm(x + ie(x)) - c.cm(x - ie(x)), disc, eps);

                    // the second sheet, approached from below, continues the first sheet's upper lip
                    check_close(c.cm(x - ie(x)) - c.disc(x - ie(x)), c.cm(x + ie(x)), eps);
                }

                // between pseudothreshold and threshold, real-valued and without a cut
                for (const double x : { s_minus + 0.25 * (s_plus - s_minus), s_minus + 0.75 * (s_plus - s_minus) })
                {
                    const complex<double> disc = c.disc(x);

                    TEST_CHECK_NEARLY_EQUAL(imag(disc), 0.0, eps * std::abs(disc));
                    check_close(c.disc(x + ie(x)), c.disc(x - ie(x)), eps);
                }

                // below pseudothreshold, a cut across which the function changes sign
                for (const double x : { s_minus - 1.0, s_minus - 10.0 })
                {
                    check_close(c.disc(x + ie(x)), -c.disc(x - ie(x)), eps);
                }

                // no cut off the real axis, in particular not at Re(s) = m1^2 + m2^2
                for (const double y : { -1.0, 1.0 })
                {
                    const double x = c.m1 * c.m1 + c.m2 * c.m2;

                    check_close(c.disc(complex<double>(x - 1e-9, y)), c.disc(complex<double>(x + 1e-9, y)), eps);
                }
            }
        }
} chew_mandelstam_discontinuity_test;
