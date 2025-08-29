/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2025 Danny van Dyk
 * Copyright (c) 2018 Nico Gubernari
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2022 Philip Lüghausen
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
#include <eos/form-factors/heavy-meson-lcdas-exponential.hh>

#include <eos/models/model.hh>

#include <array>
#include <cmath>
#include <limits>
#include <numeric>

using namespace test;
using namespace eos;
using namespace eos::heavy_meson_lcdas;

class ExponentialTest :
    public TestCase
{
    public:
        ExponentialTest() :
            TestCase("b_lcdas_exponential_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            /* m_s = u */
            /* test cases in the limit lambda_E2 = lambda_H2 as used in [KMPW:2010] */
            {
                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"] = 2.1739;
                p["B::lambda_E^2"]   = 0.3174;
                p["B::lambda_H^2"]   = 0.3174;

                /* Two-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "u" } });

                    // phi_plus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(1.0),    0.537484,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(2.0),    0.122258,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(3.0),    0.0208569,  eps);

                    // phi_minus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(1.0),   0.247243,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(2.0),   0.0281194,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(3.0),   0.00319806, eps);

                    // phi_bar at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(1.0),    -0.247243,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(2.0),    -0.0562387,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(3.0),    -0.00959419, eps);

                    // phi_bar' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(1.0),  0.290242,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(2.0),  0.0941387,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(3.0),  0.0176589,  eps);

                    // g_+ at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(1.0),      0.107355,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(2.0),      0.0505172,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(3.0),      0.0131656,  eps);

                    // g_+' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(1.0),  -0.0133214,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(2.0),  -0.0581194,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(3.0),  -0.0196547,  eps);

                    // g_+'' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(1.0),  -0.218476,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(2.0),   0.0409186,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(3.0),   0.0263377,  eps);

                    // g_- WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(1.0),     0.0852988,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(2.0),     0.0194024,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(3.0),     0.00330999,  eps);

                    // g_-' WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(1.0),  -0.100133,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(2.0),  -0.0324779,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(3.0),  -0.00609231, eps);

                    // g_-'' WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(1.0),  0.0322491,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(2.0),  0.0495145,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(3.0),  0.0108456,   eps);

                    // g_bar (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar(1.0),     -0.0349921,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar(2.0),     0.000610866, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar(3.0),     0.0199498,   eps);

                    // g_bar' (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d1(1.0),  0.0220565,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d1(2.0),  0.0311148,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d1(3.0),  0.00985562,  eps);

                    // g_bar'' (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d2(1.0),  0.0868119,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d2(2.0),  -0.0256415,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d2(3.0),  -0.0135624,  eps);

                    // g_bar''' (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d3(1.0),  -0.250726,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d3(2.0),  -0.00859595, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d3(3.0),  0.0154921,   eps);
                }

                /* Three-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "u" } });

                    // phi_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.1), 0.0         ,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(2.0, 0.1), 0.0         ,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.1), 0.0         ,  eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.3), 0.0         ,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(2.0, 0.3), 0.0         ,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.3), 0.0         ,  eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.5), 0.0         ,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(2.0, 0.5), 0.0         ,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.5), 0.0         ,  eps);

                    // phi_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.1), 2.1623400e-3,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(2.0, 0.1), 2.4592700e-4,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.1), 2.7969700e-5,  eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.3), 1.2599200e-2,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(2.0, 0.3), 1.4329300e-3,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.3), 1.6296900e-4,  eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.5), 2.2657700e-2,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(2.0, 0.5), 2.5769000e-3,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.5), 2.9307600e-4,  eps);

                    // phi_bar_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.1),  0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(2.0, 0.1),  0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.1),  0.0         , eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.3),  0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(2.0, 0.3),  0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.3),  0.0         , eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.5),  0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(2.0, 0.5),  0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.5),  0.0         , eps);

                    // phi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.1),  7.7511400e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(2.0, 0.1),  8.6326900e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.1),  8.7329500e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.3),  4.5163200e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(2.0, 0.3),  5.0299600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.3),  5.0883800e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.5),  8.1219100e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(2.0, 0.5),  9.0456300e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.5),  9.1506800e-2, eps);

                    // phi_bar2_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(1.0, 0.1), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(2.0, 0.1), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(3.0, 0.1), 0.0         , eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(1.0, 0.3), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(2.0, 0.3), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(3.0, 0.3), 0.0         , eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(1.0, 0.5), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(2.0, 0.5), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(3.0, 0.5), 0.0         , eps);

                    // phi_bar2_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(1.0, 0.1), 7.6172000e-5, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(2.0, 0.1), 8.6631700e-6, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(3.0, 0.1), 9.8527700e-7, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(1.0, 0.3), 1.4953400e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(2.0, 0.3), 1.7006800e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(3.0, 0.3), 1.9342100e-5, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(1.0, 0.5), 5.0731900e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(2.0, 0.5), 5.7698300e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(3.0, 0.5), 6.5621200e-5, eps);

                    // phi_bar_bar_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.1), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(2.0, 0.1), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.1), 0.0         , eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.3), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(2.0, 0.3), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.3), 0.0         , eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.5), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(2.0, 0.5), 0.0         , eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.5), 0.0         , eps);

                    // phi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.1), 2.7304657e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(2.0, 0.1), 3.0410063e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.1), 3.0763246e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.3), 5.3602167e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(2.0, 0.3), 5.9698433e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.3), 6.0391772e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.5), 1.8185399e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(2.0, 0.5), 2.0253656e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.5), 2.0488882e-2, eps);

                    // psi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.1), 2.5708500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(2.0, 0.1), 3.7447800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.1), 3.9785600e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.3), 4.9931400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(2.0, 0.3), 7.2731800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.3), 7.7272200e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.5), 5.3876400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(2.0, 0.5), 7.8478300e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.5), 8.3377400e-2, eps);

                    // psi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.1), 1.3838600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(2.0, 0.1), 2.0157800e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.1), 2.1416200e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.3), 9.4221500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(2.0, 0.3), 1.3724600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.3), 1.4581400e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.5), 2.0025200e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(2.0, 0.5), 2.9169400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.5), 3.0990400e-2, eps);

                    // chi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.1), 2.5708500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(2.0, 0.1), 3.7447800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.1), 3.9785600e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.3), 4.9931400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(2.0, 0.3), 7.2731800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.3), 7.7272200e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.5), 5.3876400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(2.0, 0.5), 7.8478300e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.5), 8.3377400e-2, eps);

                    // chi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.1), 1.3838600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(2.0, 0.1), 2.0157800e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.1), 2.1416200e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.3), 9.4221500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(2.0, 0.3), 1.3724600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.3), 1.4581400e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.5), 2.0025200e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(2.0, 0.5), 2.9169400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.5), 3.0990400e-2, eps);
                }
                {
                    Exponential B(p, Options{ { "q"_ok, "u" } });

                    // psi_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(1.0, 0.1), 1.0811700e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(2.0, 0.1), 0.1229630e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(3.0, 0.1), 0.0139848e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(1.0, 0.3), 6.2996000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(2.0, 0.3), 0.7164640e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(3.0, 0.3), 0.0814847e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(1.0, 0.5),11.3289000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(2.0, 0.5), 1.2884500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_A(3.0, 0.5), 0.1465389e-3, eps);

                    // psi_V at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(1.0, 0.1), 1.0811700e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(2.0, 0.1), 0.1229630e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(3.0, 0.1), 0.0139848e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(1.0, 0.3), 6.2996000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(2.0, 0.3), 0.7164640e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(3.0, 0.3), 0.0814847e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(1.0, 0.5),11.3289000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(2.0, 0.5), 1.2884500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_V(3.0, 0.5), 0.1465389e-3, eps);

                    // X_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(1.0, 0.1),  20.5422000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(2.0, 0.1),   4.7955700e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(3.0, 0.1),   0.8251050e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.X_A(1.0, 0.3),  35.6977000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(2.0, 0.3),   8.8363900e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(3.0, 0.3),   1.5482100e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.X_A(1.0, 0.5),  33.9866000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(2.0, 0.5),   9.0191600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.X_A(3.0, 0.5),   1.6119200e-3, eps);

                    // Y_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(1.0, 0.1),  25.6237000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(2.0, 0.1),   6.9105400e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(3.0, 0.1),   1.2404500e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(1.0, 0.3),  46.6170000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(2.0, 0.3),  13.0635000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(3.0, 0.3),   2.3684900e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(1.0, 0.5),  46.9015000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(2.0, 0.5),  13.7091000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Y_A(3.0, 0.5),   2.5116600e-3, eps);
                }

                /* Auxiliary functions to three-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "u" } });

                    // Xbar_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(1.0, 0.1), 2.1832900e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(2.0, 0.1), 3.3131500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(3.0, 0.1), 3.5419100e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(1.0, 0.3), 2.7349800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(2.0, 0.3), 4.7582000e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(3.0, 0.3), 5.1830300e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(1.0, 0.5), 1.3266900e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(2.0, 0.5), 3.3250100e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Xbar_A(3.0, 0.5), 3.7624000e-2, eps);

                    // Ybar_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(1.0, 0.1), 0.7671260e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(2.0, 0.1), 2.2868900e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(3.0, 0.1), 2.6226600e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(1.0, 0.3), 0.3608450e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(2.0, 0.3), 3.1841500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(3.0, 0.3), 3.8216800e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(1.0, 0.5),-0.8289310e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(2.0, 0.5), 2.0788800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.Ybar_A(3.0, 0.5), 2.7510200e-2, eps);
                }
            }

            /* test cases for lambda_E2 != lambda_H2 */
            {
                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"] = 2.1739;
                p["B::lambda_E^2"]   = 0.3174;
                p["B::lambda_H^2"]   = 1.2696;

                /* Two-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "u" } });

                    // phi_plus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(1.0),    0.537484,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(2.0),    0.122258,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(3.0),    0.0208569,  eps);

                    // phi_minus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(1.0),   0.125491,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(2.0),   0.0528102,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(3.0),   0.017946,   eps);

                    // phi_bar at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(1.0),    -0.236493,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(2.0),    -0.023229,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(3.0),     0.00125142, eps);

                    // phi_bar' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(1.0),  0.411994,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(2.0),  0.0694479,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(3.0),  0.0029109,  eps);

                    // g_+ at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(1.0),      0.122808,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(2.0),      0.057547,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(3.0),      0.0149645,  eps);

                    // g_+' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(1.0),  -0.0160088,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(2.0),  -0.0663718,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(3.0),  -0.0223661,  eps);

                    // g_+'' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(1.0),  -0.248914,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(2.0),   0.0470913,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(3.0),   0.0300247,  eps);

                    // g_- WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(1.0), 0.0852988,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(2.0), 0.0194024,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(3.0), 0.00330999,  eps);

                    // g_-' WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(1.0), -0.100133,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(2.0), -0.0324779,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(3.0), -0.00609231, eps);

                    // g_-'' WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(1.0), 0.0322491,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(2.0), 0.0495145,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(3.0), 0.0108456,   eps);

                    // g_bar (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar(1.0),    -0.0251981,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar(2.0),    0.0219975,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar(3.0),    0.0452796,   eps);

                    // g_bar' (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d1(1.0), 0.0375092,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d1(2.0), 0.0381447,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d1(3.0), 0.0116545,   eps);

                    // g_bar'' (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d2(1.0), 0.0841245,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d2(2.0), -0.033894,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d2(3.0), -0.0162738,  eps);

                    // g_bar''' (partial WW) at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d3(1.0), -0.281164,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d3(2.0), -0.00242324, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_bar_d3(3.0), 0.0191791,   eps);
                }

                /* Three-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "u" } });

                    // phi_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.1),-7.0511100e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(2.0, 0.1),-1.6038700e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.1),-2.7361600e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.3),-4.1084300e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(2.0, 0.3),-9.3451800e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.3),-1.5942700e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.5),-7.3888400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(2.0, 0.5),-1.6805900e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.5),-2.8670400e-3, eps);

                    // phi_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.1), 5.4058500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(2.0, 0.1), 6.1481700e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.1), 6.9924100e-5, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.3), 3.1498000e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(2.0, 0.3), 3.5823200e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.3), 4.0742300e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.5), 5.6644400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(2.0, 0.5), 6.4422600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.5), 7.3268900e-4, eps);

                    // phi_bar_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.1),-8.3832000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(2.0, 0.1),-1.2211300e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.1),-1.2973600e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.3),-4.8846000e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(2.0, 0.3),-7.1150700e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.3),-7.5592400e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.5),-8.7842000e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(2.0, 0.5),-1.2795400e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.5),-1.3594100e-1, eps);

                    // phi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.1), 1.9377800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(2.0, 0.1), 2.1581700e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.1), 2.1832400e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.3), 1.1290800e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(2.0, 0.3), 1.2574900e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.3), 1.2721000e-1, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.5), 2.0304800e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(2.0, 0.5), 2.2614100e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.5), 2.2876700e-1, eps);

                    // phi_bar2_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(1.0, 0.1),-2.4838700e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(2.0, 0.1),-5.6498900e-5, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(3.0, 0.1),-9.6385800e-6, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(1.0, 0.3),-4.8761200e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(2.0, 0.3),-1.1091400e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(3.0, 0.3),-1.8921600e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(1.0, 0.5),-1.6543000e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(2.0, 0.5),-3.7629300e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_3(3.0, 0.5),-6.4194700e-4, eps);

                    // phi_bar2_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(1.0, 0.1), 1.9043000e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(2.0, 0.1), 2.1657900e-5, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(3.0, 0.1), 2.4631900e-6, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(1.0, 0.3), 3.7383600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(2.0, 0.3), 4.2517000e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(3.0, 0.3), 4.8355300e-5, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(1.0, 0.5), 1.2683000e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(2.0, 0.5), 1.4424600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar2_4(3.0, 0.5), 1.6405300e-4, eps);

                    // phi_bar_bar_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.1),-2.9531189e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(2.0, 0.1),-4.3016143e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.1),-4.5701494e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.3),-5.7973103e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(2.0, 0.3),-8.4445612e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.3),-8.9717263e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.5),-1.9668310e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(2.0, 0.5),-2.8649535e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.5),-3.0438000e-2, eps);

                    // phi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.1), 6.8261600e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(2.0, 0.1), 7.6025200e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.1), 7.6908100e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.3), 1.3400500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(2.0, 0.3), 1.4924600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.3), 1.5097900e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.5), 4.5463500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(2.0, 0.5), 5.0634100e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.5), 5.1222200e-2, eps);

                    // psi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.1), 2.5708500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(2.0, 0.1), 3.7447800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.1), 3.9785600e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.3), 4.9931400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(2.0, 0.3), 7.2731800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.3), 7.7272200e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.5), 5.3876400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(2.0, 0.5), 7.8478300e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.5), 8.3377400e-2, eps);

                    // psi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.1), 1.3838600e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(2.0, 0.1), 2.0157800e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.1), 2.1416200e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.3), 9.4221500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(2.0, 0.3), 1.3724600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.3), 1.4581400e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.5), 2.0025200e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(2.0, 0.5), 2.9169400e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.5), 3.0990400e-2, eps);

                    // chi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.1), 1.0283400e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(2.0, 0.1), 1.4979100e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.1), 1.5914200e-1, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.3), 1.9972600e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(2.0, 0.3), 2.9092700e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.3), 3.0908900e-1, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.5), 2.1550600e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(2.0, 0.5), 3.1391300e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.5), 3.3351000e-1, eps);

                    // chi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.1), 5.5354400e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(2.0, 0.1), 8.0631200e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.1), 8.5664700e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.3), 3.7688600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(2.0, 0.3), 5.4898500e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.3), 5.8325600e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.5), 8.0100900e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(2.0, 0.5), 1.1667800e-1, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.5), 1.2396100e-1, eps);
                }
            }

            /* m_s = s */
            {
                Parameters p = Parameters::Defaults();
                p["B_s::1/lambda_B_p"] = 1.69348;
                p["B_s::lambda_E^2"]   = 0.5;
                p["B_s::lambda_H^2"]   = 0.8;

                /* Two-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "s" } });

                    // phi_plus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(1.0),   0.527341,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(2.0),   0.193933,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_plus(3.0),   0.0534902,  eps);

                    // phi_minus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(1.0),  0.283025,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(2.0),  0.0570498,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_minus(3.0),  0.0142975,  eps);

                    // phi_bar at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(1.0),   -0.315957,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(2.0),   -0.106926,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar(3.0),   -0.0269353,  eps);

                    // phi_bar' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(1.0), 0.244316,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(2.0), 0.136883,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_d1(3.0), 0.0391927,  eps);

                    // g_+ at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(1.0),     0.138237,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(2.0),     0.104952,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus(3.0),     0.0441973,  eps);

                    // g_+' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(1.0),  0.0487100,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(2.0), -0.0704267,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d1(3.0), -0.0447601,  eps);

                    // g_+'' at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(1.0), -0.265357,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(2.0), -0.00652033, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_plus_d2(3.0),  0.0352610,  eps);

                    // g_- WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(1.0),  0.137909,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(2.0),  0.050717,   eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW(3.0),  0.0139886,  eps);

                    // g_-' WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(1.0),  -0.0956371, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(2.0),  -0.0605297, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d1(3.0),  -0.0190266, eps);

                    // g_-'' WW at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(1.0), -0.0715865,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(2.0),  0.0595617,  eps);
                    TEST_CHECK_NEARLY_EQUAL(B.g_minusWW_d2(3.0),  0.0243247,  eps);
                }

                /* Three-particle LCDAs */
                {
                    Exponential B(p, Options{ { "q"_ok, "s" } });

                    // phi_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.1), -1.0810700e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.1), -1.0965700e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(1.0, 0.3), -6.9343100e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_3(3.0, 0.3), -7.0337400e-4, eps);

                    // phi_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.1),  2.7662800e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.1),  9.3531600e-5, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(1.0, 0.3),  1.7743700e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_4(3.0, 0.3),  5.9993900e-4, eps);

                    // phi_bar_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.1), -1.0347100e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.1), -1.9725500e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(1.0, 0.3), -6.6369500e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_3(3.0, 0.3), -1.2652500e-2, eps);

                    // phi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.1),  7.2500400e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.1),  8.8283000e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(1.0, 0.3),  4.6503900e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_4(3.0, 0.3),  5.6627300e-2, eps);

                    // phi_bar_bar_3 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.1), -3.6001500e-5, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.1), -6.8632300e-5, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(1.0, 0.3), -7.5733800e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_3(3.0, 0.3), -1.4437700e-3, eps);

                    // phi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.1),  2.5225600e-4, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.1),  3.0717000e-4, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(1.0, 0.3),  5.3065300e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.phi_bar_bar_4(3.0, 0.3),  6.4617000e-3, eps);

                    // psi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.1),  2.0366600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.1),  3.8826300e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(1.0, 0.3),  4.3545800e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_4(3.0, 0.3),  8.3014400e-2, eps);

                    // psi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.1),  1.0783300e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.1),  2.0557000e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(1.0, 0.3),  7.7941000e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.psi_bar_bar_4(3.0, 0.3),  1.4858400e-2, eps);

                    // chi_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.1),  3.2586600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.1),  6.2122100e-2, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(1.0, 0.3),  6.9673200e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_4(3.0, 0.3),  1.3282300e-1, eps);

                    // chi_bar_bar_4 at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.1 GeV, 0.3 GeV and 0.5 GeV
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.1),  1.7253300e-3, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.1),  3.2891200e-3, eps);

                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(1.0, 0.3),  1.2470600e-2, eps);
                    TEST_CHECK_NEARLY_EQUAL(B.chi_bar_bar_4(3.0, 0.3),  2.3773500e-2, eps);
                }
            }

            /* coefficient interface */
            {
                Parameters p = Parameters::Defaults();
                Exponential B(p, Options{ { "q"_ok, "u" } });
                auto [c, c_end] = B.coefficient_range(1.0);

                std::array<double, 9> ref = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                for (auto it = c; it != c_end; ++it)
                {
                    TEST_CHECK_NEARLY_EQUAL(*it, ref[std::distance(c, it)], 1e-15);
                }

            }
        }
} b_lcdas_exponential_test;
