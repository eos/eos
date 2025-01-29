/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/top-loops.hh>
#include <eos/models/wilson-coefficients.hh>

#include <test/test.hh>

#include <gsl/gsl_sf_clausen.h>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class WilsonCoefficientsTest : public TestCase
{
    public:
        WilsonCoefficientsTest() :
            TestCase("wilson_coefficients_test")
        {
        }

        virtual void
        run() const
        {
            /* Test for 5 active flavors, evolving from mu_0c = 80, mu_0t = 120 to mu = m_b = 4.2 */
            {
                static const double            alpha_s_0 = 0.119918, alpha_s = 0.220996781524799;
                static const double            mu_0c = 80, mu_0t = 120;
                static const double            m_t_mu_0c = 175.960604211196, m_t_mu_0t = 170.067834276559, m_W = 80.403;
                static const double            log_charm = 2.0 * std::log(mu_0c / m_W), log_top = std::log(mu_0t / m_t_mu_0t);
                static const double            x_c = power_of<2>(m_t_mu_0c / m_W), x_t = power_of<2>(m_t_mu_0t / m_W);
                static const double            sw2 = 0.23122;
                static const double            nf  = 5.0;
                static const QCD::BetaFunction beta{
                    {
                     23.0 / 3.0,
                     116.0 / 3.0,
                     9769.0 / 54.0,
                     4826.1563287908967, }
                };

                /* Charm Sector */
                std::array<complex<double>, 15> initial_charm_qcd_0{
                    { 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                };
                std::array<complex<double>, 15> initial_charm_qcd_1{
                    { -15.0 - 6.0 * log_charm,
                     0.0, 0.0,
                     7.0 / 9.0 - 2.0 / 3.0 * log_charm,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 0.0,
                     0.0, 23.0 / 36.0,
                     1.0 / 3.0,
                     -0.25 / sw2 - 38.0 / 27.0,
                     0.25 / sw2 }
                };
                std::array<complex<double>, 15> initial_charm_qcd_2{
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                };
                initial_charm_qcd_2[0] = -(16.0 * x_c + 8.0) * std::sqrt(4.0 * x_c - 1.0) * gsl_sf_clausen(2.0 * std::asin(1.0 / 2.0 / std::sqrt(x_c)))
                                         + (16.0 * x_c + 20.0 / 3.0) * std::log(x_c) + 32.0 * x_c + 112.0 / 9.0 - 7987.0 / 72.0 - 17.0 / 3.0 * M_PI * M_PI - 475.0 / 6.0 * log_charm
                                         - 17.0 * log_charm * log_charm;
                initial_charm_qcd_2[1]  = -127.0 / 18.0 - 4.0 / 3.0 * M_PI * M_PI - 46.0 / 3.0 * log_charm - 4.0 * log_charm * log_charm;
                initial_charm_qcd_2[2]  = 680.0 / 243.0 + 20.0 / 81.0 * M_PI * M_PI + 68.0 / 81.0 * log_charm + 20.0 / 27.0 * log_charm * log_charm;
                initial_charm_qcd_2[3]  = -950.0 / 243.0 - 10.0 / 81.0 * M_PI * M_PI - 124.0 / 27.0 * log_charm - 10.0 / 27.0 * log_charm * log_charm;
                initial_charm_qcd_2[4]  = -68.0 / 243.0 - 2.0 / 81.0 * M_PI * M_PI - 14.0 / 81.0 * log_charm - 2.0 / 27.0 * log_charm * log_charm;
                initial_charm_qcd_2[5]  = -85.0 / 162.0 - 5.0 / 108.0 * M_PI * M_PI - 35.0 / 108.0 * log_charm - 5.0 / 36.0 * log_charm * log_charm;
                initial_charm_qcd_2[11] = -713.0 / 243.0 - 4.0 / 81.0 * log_charm;
                initial_charm_qcd_2[12] = -91.0 / 324.0 + 4.0 / 27.0 * log_charm;
                initial_charm_qcd_2[13] = -1.0 / sw2 - 524.0 / 729.0 + 128.0 / 243.0 * M_PI * M_PI + 16.0 / 3.0 * log_charm + 128.0 / 81.0 * log_charm * log_charm;
                initial_charm_qcd_2[14] = 1.0 / sw2;

                /* Top Sector */
                std::array<complex<double>, 15> initial_top_qcd_0{
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                };
                std::array<complex<double>, 15> initial_top_qcd_1{
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                };
                initial_top_qcd_1[3]  = TopLoops::E0(x_t);
                initial_top_qcd_1[11] = -0.5 * TopLoops::A0(x_t);
                initial_top_qcd_1[12] = -0.5 * TopLoops::F0(x_t);
                initial_top_qcd_1[13] = (1.0 - 4.0 * sw2) / sw2 * TopLoops::C0(x_t) - TopLoops::B0(x_t) / sw2 - TopLoops::D0(x_t);
                initial_top_qcd_1[14] = (TopLoops::B0(x_t) - TopLoops::C0(x_t)) / sw2;

                std::array<complex<double>, 15> initial_top_qcd_2{
                    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
                };
                initial_top_qcd_2[2]  = TopLoops::G1(x_t, log_top);
                initial_top_qcd_2[3]  = TopLoops::E1(x_t, log_top);
                initial_top_qcd_2[4]  = -0.1 * TopLoops::G1(x_t, log_top) + 2.0 / 15.0 * TopLoops::E0(x_t);
                initial_top_qcd_2[5]  = -3.0 / 16.0 * TopLoops::E1(x_t, log_top) + 0.25 * TopLoops::E0(x_t);
                initial_top_qcd_2[11] = -0.5 * TopLoops::A1(x_t, log_top);
                initial_top_qcd_2[12] = -0.5 * TopLoops::F1(x_t, log_top);
                initial_top_qcd_2[13] = (1.0 - 4.0 * sw2) / sw2 * TopLoops::C1(x_t, log_top) - TopLoops::B1(x_t, log_top) / sw2 - TopLoops::D1(x_t, log_top);
                initial_top_qcd_2[14] = (TopLoops::B1(x_t, log_top) - TopLoops::C1(x_t, log_top)) / sw2;

                WilsonCoefficients<BToS> downscaled_charm = evolve(initial_charm_qcd_0, initial_charm_qcd_1, initial_charm_qcd_2, alpha_s_0, alpha_s, nf, beta);
                WilsonCoefficients<BToS> downscaled_top   = evolve(initial_top_qcd_0, initial_top_qcd_1, initial_top_qcd_2, alpha_s_0, alpha_s, nf, beta);

                WilsonCoefficients<BToS> wc = downscaled_top;
                wc._sm_like_coefficients    = wc._sm_like_coefficients + complex<double>(-1.0, 0.0) * downscaled_charm._sm_like_coefficients;

                static const double eps = 1e-10;
                TEST_CHECK_NEARLY_EQUAL(real(wc.c1()), -0.280597305725000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c2()), +1.009729183490000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c3()), -0.005821929583920, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c4()), -0.084146804435200, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c5()), +0.000405161899729, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c6()), +0.001090753111680, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c7()), -0.326845973907000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c8()), -0.176431447889000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c9()), +4.274693790550000, eps);
                TEST_CHECK_NEARLY_EQUAL(real(wc.c10()), -4.227994731400000, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c1()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c2()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c3()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c4()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c5()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c6()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c7()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c8()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c9()), +0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(wc.c10()), +0.0, eps);
            }
        }
} wilson_coefficients_test;
