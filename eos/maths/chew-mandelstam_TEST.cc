/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
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
#include <eos/utils/exception.hh>

#include <test/test.hh>

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
            constexpr double m = 5.10999e-4;

            // real s, above threshold
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

            // real s, below threshold: the function is real-valued
            {
                const complex<double> result = chew_mandelstam::s_wave(-20.0, m, m);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.11496167, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, eps);
            }

            // genuinely complex s, e.g. as used to move onto another Riemann sheet
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

            // only equal masses are implemented
            {
                TEST_CHECK_THROWS(InternalError, chew_mandelstam::s_wave(20.0, m, 2.0 * m));
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
            constexpr double q0 = 0.5;

            // real s, below threshold: the function is real-valued
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

            // real s, above threshold
            {
                const complex<double> result = chew_mandelstam::p_wave(20.0, mD, mD, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), 0.00024738, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0093576, eps);
            }

            // effective channel mass, cf. eos/scattering/ee-to-ccbar, EffChannel
            {
                constexpr double      meff   = 0.1349768;
                const complex<double> result = chew_mandelstam::p_wave(20.0, meff, meff, q0);

                TEST_CHECK_RELATIVE_ERROR(real(result), -0.0172797, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(result), 0.0189092, eps);
            }

            // right at threshold, where the Taylor-expanded branch is taken
            {
                const double          mp2    = power_of<2>(2.0 * mD);
                const complex<double> result = chew_mandelstam::p_wave(mp2, mD, mD, q0);

                TEST_CHECK_NEARLY_EQUAL(real(result), 0.0, 1e-6);
                TEST_CHECK_NEARLY_EQUAL(imag(result), 0.0, 1e-6);
            }

            // genuinely complex s, e.g. as used to move onto another Riemann sheet
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

            // only equal masses are implemented
            {
                TEST_CHECK_THROWS(InternalError, chew_mandelstam::p_wave(20.0, mD, 2.0 * mD, q0));
            }
        }
} chew_mandelstam_p_wave_test;
