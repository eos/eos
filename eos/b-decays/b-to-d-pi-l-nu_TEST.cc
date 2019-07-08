/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
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

#include <eos/utils/integrate.hh>

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
            Parameters p = Parameters::Defaults();

            {
                p[      "B(*)->D(*)::rho^2@HQET"] = +1.3060e+00;
                p[          "B(*)->D(*)::c@HQET"] = +1.2200e+00;
                p[          "B(*)->D(*)::d@HQET"] = +0.0000e+00;
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

                Options o{ { "z3", "0" } };

                BToDPiLeptonNeutrino d(p, o);

                const double eps = 1e-5;

                // distribution in cos(theta_D)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-1.00), 0.796712, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.80), 0.636487, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.60), 0.511868, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.40), 0.422855, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(-0.20), 0.369447, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d( 0.00), 0.351644, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.20), 0.369447, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.40), 0.422855, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.60), 0.511868, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+0.80), 0.636487, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_d(+1.00), 0.796712, eps);

                // distribution in cos(theta_l)
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-1.00), 0.133742, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.80), 0.255721, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.60), 0.360502, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.40), 0.448085, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(-0.20), 0.518470, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l( 0.00), 0.571658, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.20), 0.607647, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.40), 0.626439, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.60), 0.628033, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+0.80), 0.612429, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_l(+1.00), 0.579628, eps);

                // distribution in chi
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 1.5708),  0.186515, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 2.35619), 0.168450, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 3.14159), 0.144939, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 3.45575), 0.149521, eps);
                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_chi( 5.49779), 0.149860, eps);

                // normalization of the integrated distributions
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0,  0.0),     0.49999, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_d(-1.0, +1.0),     1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0,  0.0),     0.38853, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_l(-1.0, +1.0),     1.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi( 0.0,  +M_PI), 0.50000, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_chi(-M_PI, +M_PI), 1.0,     eps);
            }
        }
} b_to_d_pi_l_nu_test;
