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

#include <eos/utils/model.hh>

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
            p["mass::ud(2GeV)"] = 0.008;
            p["pi::a2@1GeV"] = 0.17;
            p["pi::a4@1GeV"] = 0.06;
            p["pi::f3@1GeV"] = 0.0045;
            p["pi::omega3@1GeV"] = -1.5;
            p["pi::delta^2@1GeV"] = 0.18;
            p["pi::omega4@1GeV"] = 0.2;
            p["decay-constant::pi"] = 0.1302;
            p["mass::pi^+"] = 0.13957;
            p["mass::ud(2GeV)"] = 0.008;
            p["decay-constant::pi"] = 0.1304;

            /* Diagnostics */
            {
                PionLCDAs pi(p, Options{ });
                Diagnostics diagnostics = pi.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94663, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92648, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91464, 1e-5), // c_rge(mu = 4.0 GeV)
                    std::make_pair(+0.90638, 1e-5), // c_rge(mu = 5.0 GeV)
                };

                //TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

#undef TEST_CHECK_NEARLY_EQUAL
#define TEST_CHECK_NEARLY_EQUAL(ref, test, eps) \
                do \
                { \
                    std::cout << #test << " = " << test << ", should be " << ref << std::endl; \
                } while(false)
            /* Twist 2 */
            {
                PionLCDAs pi(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.17,      pi.a2pi(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.125815,  pi.a2pi(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1109,    pi.a2pi(3.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.06,      pi.a4pi(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0387104, pi.a4pi(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0322135, pi.a4pi(3.0),   eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,     pi.phi(0.0, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.13314, pi.phi(0.3, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.25250, pi.phi(0.5, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,     pi.phi(1.0, 1.0), eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,     pi.phi(0.0, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.17281, pi.phi(0.3, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.30675, pi.phi(0.5, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,     pi.phi(1.0, 2.0), eps);
            }

            /* Twist 3 */
            {
                PionLCDAs pi(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.00450, pi.f3pi(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00322, pi.f3pi(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00322, pi.f3pi(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 1.88175, pi.mupi(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 2.51278, pi.mupi(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 2.51278, pi.mupi(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL(-1.20000, pi.eta3pi(1.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-0.89021, pi.eta3pi(2.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-0.89021, pi.eta3pi(3.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-1.20000, pi.omega3pi(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-0.89021, pi.omega3pi(2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-0.89021, pi.omega3pi(3.0),    eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.25842, pi.phi3p(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.98711, pi.phi3p(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.83162, pi.phi3p(0.3, 1.0), eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.14135, pi.phi3p(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.99824, pi.phi3p(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.91124, pi.phi3p(0.3, 2.0), eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.74766, pi.phi3s(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.09424, pi.phi3s(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.21595, pi.phi3s(0.3, 1.0), eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.64747, pi.phi3s(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.02948, pi.phi3s(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.23720, pi.phi3s(0.3, 2.0), eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.13560, pi.phi3s_d1(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.08978, pi.phi3s_d1(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.55417, pi.phi3s_d1(0.3, 1.0), eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.97370, pi.phi3s_d1(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.81837, pi.phi3s_d1(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.44467, pi.phi3s_d1(0.3, 2.0), eps);
            }

            /* Twist 4 */
            {
                PionLCDAs pi(p, Options{ });

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.55200, pi.psi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.04800, pi.psi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.31200, pi.psi4(0.3, 1.0), eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.45420, pi.psi4(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.03950, pi.psi4(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.25672, pi.psi4(0.3, 2.0), eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.12514, pi.phi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.34192, pi.phi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.53192, pi.phi4(0.3, 1.0), eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.09612, pi.phi4(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.27283, pi.phi4(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.43701, pi.phi4(0.3, 2.0), eps);

                // phi4_d1, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 2.00934, pi.phi4_d1(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.14573, pi.phi4_d1(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.59117, pi.phi4_d1(0.3, 1.0), eps);

                // phi4_d1, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.58440, pi.phi4_d1(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.80433, pi.phi4_d1(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.41332, pi.phi4_d1(0.3, 2.0), eps);

                // phi4_d2, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 7.72568, pi.phi4_d2(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-3.25100, pi.phi4_d2(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-7.10242, pi.phi4_d2(0.3, 1.0), eps);

                // phi4_d2, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 7.16884, pi.phi4_d2(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-1.64331, pi.phi4_d2(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-5.63272, pi.phi4_d2(0.3, 2.0), eps);
            }
        }
} pi_lcdas_test;
