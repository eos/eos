/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2025 Danny van Dyk
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
#include <eos/observable.hh>
#include <eos/b-decays/b-to-d-pi-l-nu.hh>

#include <eos/maths/integrate.hh>

using namespace test;
using namespace eos;

class BToDPiLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToDPiLeptonNeutrinoTest() :
            TestCase("b_to_d_pi_l_nu_test")
        {
        }

        virtual void run() const
        {

            // agreement with Martin Jung in 3/2/1 model
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u"].set(1.86723);
                p["mass::D_d"].set(1.86723);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);

                Options o{
                    { "l"_ok,             "tau"       },
                    { "q"_ok,             "d"         },
                    { "z-order-lp"_ok,    "3"         },
                    { "z-order-slp"_ok,   "2"         },
                    { "z-order-sslp"_ok,  "1"         },
                    { "form-factors"_ok,  "BGJvD2019" }
                };

                BToDPiLeptonNeutrino d(p, o);

                const double eps = 1e-5;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_lepton_polarization(3.157, 10.689), 0.484992,     eps);
            }

            {
                Parameters p = Parameters::Defaults();

                p[     "B(*)->D(*)::xi'(1)@HQET"] = -1.3060e+00;
                p[    "B(*)->D(*)::xi''(1)@HQET"] = 2.0 * +1.2200e+00;
                p[   "B(*)->D(*)::xi'''(1)@HQET"] = +0.0000e+00;
                p[   "B(*)->D(*)::chi_2(1)@HQET"] = -5.8000e-02;
                p[  "B(*)->D(*)::chi_2'(1)@HQET"] = -1.0000e-03;
                p[ "B(*)->D(*)::chi_2''(1)@HQET"] = +0.0000e+00;
                p[  "B(*)->D(*)::chi_3'(1)@HQET"] = +3.5000e-02;
                p[ "B(*)->D(*)::chi_3''(1)@HQET"] = +0.0000e+00;
                p[     "B(*)->D(*)::eta(1)@HQET"] = +3.5800e-01;
                p[    "B(*)->D(*)::eta'(1)@HQET"] = +4.4000e-02;
                p[   "B(*)->D(*)::eta''(1)@HQET"] = +0.0000e+00;
                p[     "B(*)->D(*)::l_1(1)@HQET"] = +4.8500e-01;
                p[     "B(*)->D(*)::l_2(1)@HQET"] = -2.2990e+00;
                p[     "B(*)->D(*)::l_3(1)@HQET"] = +0.0000e+00;
                p[     "B(*)->D(*)::l_4(1)@HQET"] = +0.0000e+00;
                p[     "B(*)->D(*)::l_5(1)@HQET"] = +0.0000e+00;
                p[     "B(*)->D(*)::l_6(1)@HQET"] = +0.0000e+00;
                p[    "B(*)->D(*)::l_1'(1)@HQET"] = +0.0000e+00;
                p[    "B(*)->D(*)::l_2'(1)@HQET"] = +0.0000e+00;
                p[    "B(*)->D(*)::l_3'(1)@HQET"] = +0.0000e+00;
                p[    "B(*)->D(*)::l_4'(1)@HQET"] = +0.0000e+00;
                p[    "B(*)->D(*)::l_5'(1)@HQET"] = +0.0000e+00;
                p[    "B(*)->D(*)::l_6'(1)@HQET"] = +0.0000e+00;

                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u"].set(1.86723);
                p["mass::D_d"].set(1.86723);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);

                Options o{
                    { "l"_ok,             "mu"        },
                    { "q"_ok,             "d"         },
                    { "z-order-lp"_ok,    "2"         },
                    { "z-order-slp"_ok,   "1"         },
                    { "z-order-sslp"_ok,  "0"         },
                    { "form-factors"_ok,  "BGJvD2019" }
                };

                BToDPiLeptonNeutrino d(p, o);

                const double eps = 1e-5;

                // distribution in cos(theta_D)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-1.00), 0.785481, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.80), 0.631321, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.60), 0.511419, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.40), 0.425775, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.20), 0.374389, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.00), 0.357260, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.20), 0.374389, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.40), 0.425775, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.60), 0.511419, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.80), 0.631321, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+1.00), 0.785481, eps);

                // distribution in cos(theta_l)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-1.00), 0.133416, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.80), 0.254246, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.60), 0.358307, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.40), 0.445598, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.20), 0.516120, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.00), 0.569873, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.20), 0.606856, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.40), 0.627069, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.60), 0.630513, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.80), 0.617188, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+1.00), 0.587094, eps);

                // distribution in chi
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 1.5708),  0.186796, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 2.35619), 0.168542, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 3.14159), 0.144789, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 3.45575), 0.149418, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 5.49779), 0.149768, eps);

                // normalization of the integrated distributions
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0,  0.0),     0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0, +1.0),     1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0,  0.0),     0.38658, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0, +1.0),     1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi( 0.0,  +M_PI), 0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi(-M_PI, +M_PI), 1.0,     eps);
            }
        }
} b_to_d_pi_l_nu_test;
