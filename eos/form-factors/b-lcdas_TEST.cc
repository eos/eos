/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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
#include <eos/form-factors/b-lcdas.hh>

#include <eos/utils/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BMesonLCDAsTest :
    public TestCase
{
    public:
        BMesonLCDAsTest() :
            TestCase("b_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["lambda_B_p"] = 0.460;

            /* Two-particle LCDAs */
            {
                BMesonLCDAs B(p, Options{ });

                // phi_plus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.537484,   B.phi_plus(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.122258,   B.phi_plus(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0208569,  B.phi_plus(3.0),   eps);

                // phi_minus at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.247242,   B.phi_minus(1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0281194,  B.phi_minus(2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00319806, B.phi_minus(3.0),  eps);
            }

            /* Three-particle LCDAs */
            {
                BMesonLCDAs B(p, Options{ });

                // psi_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.0811700e-3, B.psi_A(1.0, 0.1), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1229630e-3, B.psi_A(2.0, 0.1), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0139848e-3, B.psi_A(3.0, 0.1), eps);

                TEST_CHECK_NEARLY_EQUAL( 6.2996000e-3, B.psi_A(1.0, 0.3), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.7164640e-3, B.psi_A(2.0, 0.3), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0814847e-3, B.psi_A(3.0, 0.3), eps);

                TEST_CHECK_NEARLY_EQUAL(11.3289000e-3, B.psi_A(1.0, 0.5), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.2884500e-3, B.psi_A(2.0, 0.5), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1465389e-3, B.psi_A(3.0, 0.5), eps);

                // psi_V at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.0811700e-3, B.psi_V(1.0, 0.1), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1229630e-3, B.psi_V(2.0, 0.1), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0139848e-3, B.psi_V(3.0, 0.1), eps);

                TEST_CHECK_NEARLY_EQUAL( 6.2996000e-3, B.psi_V(1.0, 0.3), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.7164640e-3, B.psi_V(2.0, 0.3), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0814847e-3, B.psi_V(3.0, 0.3), eps);

                TEST_CHECK_NEARLY_EQUAL(11.3289000e-3, B.psi_V(1.0, 0.5), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.2884500e-3, B.psi_V(2.0, 0.5), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1465389e-3, B.psi_V(3.0, 0.5), eps);

                // X_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                TEST_CHECK_NEARLY_EQUAL(20.5422000e-3, B.X_A(1.0, 0.1),   eps);
                TEST_CHECK_NEARLY_EQUAL( 4.7955700e-3, B.X_A(2.0, 0.1),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.8251050e-3, B.X_A(3.0, 0.1),   eps);

                TEST_CHECK_NEARLY_EQUAL(35.6977000e-3, B.X_A(1.0, 0.3),   eps);
                TEST_CHECK_NEARLY_EQUAL( 8.8363900e-3, B.X_A(2.0, 0.3),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.5482100e-3, B.X_A(3.0, 0.3),   eps);

                TEST_CHECK_NEARLY_EQUAL(33.9866000e-3, B.X_A(1.0, 0.5),   eps);
                TEST_CHECK_NEARLY_EQUAL( 9.0191600e-3, B.X_A(2.0, 0.5),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.6119200e-3, B.X_A(3.0, 0.5),   eps);

                // Y_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                TEST_CHECK_NEARLY_EQUAL(25.6237000e-3, B.Y_A(1.0, 0.1),   eps);
                TEST_CHECK_NEARLY_EQUAL( 6.9105400e-3, B.Y_A(2.0, 0.1),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.2404500e-3, B.Y_A(3.0, 0.1),   eps);

                TEST_CHECK_NEARLY_EQUAL(46.6170000e-3, B.Y_A(1.0, 0.3),   eps);
                TEST_CHECK_NEARLY_EQUAL(13.0635000e-3, B.Y_A(2.0, 0.3),   eps);
                TEST_CHECK_NEARLY_EQUAL( 2.3684900e-3, B.Y_A(3.0, 0.3),   eps);

                TEST_CHECK_NEARLY_EQUAL(46.9015000e-3, B.Y_A(1.0, 0.5),   eps);
                TEST_CHECK_NEARLY_EQUAL(13.7091000e-3, B.Y_A(2.0, 0.5),   eps);
                TEST_CHECK_NEARLY_EQUAL( 2.5116600e-3, B.Y_A(3.0, 0.5),   eps);
            }

            /* Auxiliary functions to three-particle LCDAs */
            {
                BMesonLCDAs B(p, Options{ });

                // Xbar_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                TEST_CHECK_NEARLY_EQUAL( 2.1832900e-2, B.Xbar_A(1.0, 0.1),eps);
                TEST_CHECK_NEARLY_EQUAL( 3.3131500e-2, B.Xbar_A(2.0, 0.1),eps);
                TEST_CHECK_NEARLY_EQUAL( 3.5419100e-2, B.Xbar_A(3.0, 0.1),eps);

                TEST_CHECK_NEARLY_EQUAL( 2.7349800e-2, B.Xbar_A(1.0, 0.3),eps);
                TEST_CHECK_NEARLY_EQUAL( 4.7582000e-2, B.Xbar_A(2.0, 0.3),eps);
                TEST_CHECK_NEARLY_EQUAL( 5.1830300e-2, B.Xbar_A(3.0, 0.3),eps);

                TEST_CHECK_NEARLY_EQUAL( 1.3266900e-2, B.Xbar_A(1.0, 0.5),eps);
                TEST_CHECK_NEARLY_EQUAL( 3.3250100e-2, B.Xbar_A(2.0, 0.5),eps);
                TEST_CHECK_NEARLY_EQUAL( 3.7624000e-2, B.Xbar_A(3.0, 0.5),eps);

                // Ybar_A at omega = 1.0 GeV, 2.0 GeV and 3.0 GeV, and xi = 0.5 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.7671260e-2, B.Ybar_A(1.0, 0.1),eps);
                TEST_CHECK_NEARLY_EQUAL( 2.2868900e-2, B.Ybar_A(2.0, 0.1),eps);
                TEST_CHECK_NEARLY_EQUAL( 2.6226600e-2, B.Ybar_A(3.0, 0.1),eps);

                TEST_CHECK_NEARLY_EQUAL( 0.3608450e-2, B.Ybar_A(1.0, 0.3),eps);
                TEST_CHECK_NEARLY_EQUAL( 3.1841500e-2, B.Ybar_A(2.0, 0.3),eps);
                TEST_CHECK_NEARLY_EQUAL( 3.8216800e-2, B.Ybar_A(3.0, 0.3),eps);

                TEST_CHECK_NEARLY_EQUAL(-0.8289310e-2, B.Ybar_A(1.0, 0.5),eps);
                TEST_CHECK_NEARLY_EQUAL( 2.0788800e-2, B.Ybar_A(2.0, 0.5),eps);
                TEST_CHECK_NEARLY_EQUAL( 2.7510200e-2, B.Ybar_A(3.0, 0.5),eps);
            }
        }
} b_lcdas_test;
