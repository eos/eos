/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/form-factors/pi-lcdas.hh>

#include <eos/models/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class PionLCDAsTest :
    public TestCase
{
    public:
        PionLCDAsTest() :
            TestCase("pi_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"] = 0.1176;
            p["mass::d(2GeV)"]    = 0.0048;
            p["mass::u(2GeV)"]    = 0.0032;
            p["pi::a2@1GeV"] = 0.17;
            p["pi::a4@1GeV"] = 0.06;
            p["pi::f3@1GeV"] = 0.0045;
            p["pi::omega3@1GeV"] = -1.5;
            p["pi::delta4@1GeV"] = 0.18;
            p["pi::omega4@1GeV"] = 0.2;
            p["decay-constant::pi"] = 0.1302;
            p["mass::pi^+"] = 0.13957;
            p["decay-constant::pi"] = 0.1304;

            /* Diagnostics */
            {
                PionLCDAs pi(p, Options{ });
                Diagnostics diagnostics = pi.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94850, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92874, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91708, 1e-5), // c_rge(mu = 4.0 GeV)
                    std::make_pair(+0.90893, 1e-5), // c_rge(mu = 5.0 GeV)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* Twist 2 */
            {
                PionLCDAs pi(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.a1(1.0),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a1(2.0),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a1(3.0),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a2(1.0),  0.17,       eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a2(2.0),  0.126731,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a2(3.0),  0.112741,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a3(1.0),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a3(2.0),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a3(3.0),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a4(1.0),  0.06,       eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a4(2.0),  0.0391213,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.a4(3.0),  0.0329952,  eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.3, 1.0), 1.09617,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.5, 1.0), 1.28625,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.7, 1.0), 1.09617,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(1.0, 1.0), 0.0,      eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.3, 2.0), 1.14718,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.5, 2.0), 1.32488,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(0.7, 2.0), 1.14718,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi(1.0, 2.0), 0.0,      eps);
            }

            /* Twist 3 */
            {
                PionLCDAs pi(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.f3(1.0),       0.0045,           eps);
                TEST_CHECK_NEARLY_EQUAL(pi.f3(2.0),       0.003257533016,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.f3(3.0),       0.002864248153,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.mu3(1.0),      1.846891546,      eps);
                TEST_CHECK_NEARLY_EQUAL(pi.mu3(2.0),      2.434973113,      eps);
                TEST_CHECK_NEARLY_EQUAL(pi.mu3(3.0),      2.697044518,      eps);
                TEST_CHECK_NEARLY_EQUAL(pi.eta3(1.0),     0.01868720856,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.eta3(2.0),     0.0102592843,     eps);
                TEST_CHECK_NEARLY_EQUAL(pi.eta3(3.0),     0.008143983155,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.omega3(1.0),  -1.5,              eps);
                TEST_CHECK_NEARLY_EQUAL(pi.omega3(2.0),  -1.124800783,      eps);
                TEST_CHECK_NEARLY_EQUAL(pi.omega3(3.0),  -1.002982205,      eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi3p(0.1, 1.0),    1.23826201,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3p(0.2, 1.0),    0.9881149354, eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3p(0.3, 1.0),    0.8447555249, eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi3p(0.1, 2.0),    1.133511907,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3p(0.2, 2.0),    0.9981866083, eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3p(0.3, 2.0),    0.9160656408, eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s(0.1, 1.0),    0.7314560406, eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s(0.2, 1.0),    1.083769561,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s(0.3, 1.0),    1.219383352,  eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s(0.1, 2.0),    0.6416920521, eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s(0.2, 2.0),    1.025740317,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s(0.3, 2.0),    1.238428959,  eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s_d1(0.1, 1.0), 5.109423904,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s_d1(0.2, 1.0), 2.207592431,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s_d1(0.3, 1.0), 0.6981685276, eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s_d1(0.1, 2.0), 4.964350791,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s_d1(0.2, 2.0), 2.860421439,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi3s_d1(0.3, 2.0), 1.496070648,  eps);
            }

            /* Twist 4 */
            {
                PionLCDAs pi(p, Options{ });

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.psi4(0.1, 1.0),     0.55200,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.psi4(0.2, 1.0),     0.04800,        eps);
                TEST_CHECK_NEARLY_EQUAL(pi.psi4(0.3, 1.0),    -0.31200,        eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.psi4(0.1, 2.0),     0.457400047,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.psi4(0.2, 2.0),     0.03977391713,  eps);
                TEST_CHECK_NEARLY_EQUAL(pi.psi4(0.3, 2.0),    -0.2585304614,   eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi4(0.1, 1.0),     0.1251436378,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4(0.2, 1.0),     0.3419189663,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4(0.3, 1.0),     0.5319181541,   eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi4(0.1, 2.0),     0.0970109977,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4(0.2, 2.0),     0.2750149814,   eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4(0.3, 2.0),     0.4401094031,   eps);

                // phi4_d1, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d1(0.1, 1.0),  2.009344535,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d1(0.2, 1.0),  2.145732807,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d1(0.3, 1.0),  1.591165695,    eps);

                // phi4_d1, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d1(0.1, 2.0),  1.59767376,     eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d1(0.2, 2.0),  1.815871935,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d1(0.3, 2.0),  1.420124887,    eps);

                // phi4_d2, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d2(0.1, 1.0),  7.725680617,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d2(0.2, 1.0), -3.251003946,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d2(0.3, 1.0), -7.102416399,    eps);

                // phi4_d2, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d2(0.1, 2.0),  7.194679231,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d2(0.2, 2.0), -1.686311876,    eps);
                TEST_CHECK_NEARLY_EQUAL(pi.phi4_d2(0.3, 2.0), -5.678881509,    eps);
            }
        }
} pi_lcdas_test;
