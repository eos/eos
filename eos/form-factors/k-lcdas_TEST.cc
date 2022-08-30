/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
 * Copyright (c) 2022 Carolina Bolognani
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
#include <eos/form-factors/k-lcdas.hh>

#include <eos/models/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class KaonLCDAsTest :
    public TestCase
{
    public:
        KaonLCDAsTest() :
            TestCase("k_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"]    =  0.1176;
            p["mass::s(2GeV)"]       =  0.095;
            p["mass::u(2GeV)"]       =  0.0085;  // we use 8.5 MeV for twice the average u/d quark mass
            p["mass::d(2GeV)"]       =  0.0;
            p["K::a1@1GeV"]          = -0.0525;
            p["K::a2@1GeV"]          =  0.106;
            p["K::f3@1GeV"]          =  0.0045;
            p["mass::K_u"]           =  0.49368;
            p["decay-constant::K_u"] =  0.1561;

            /* Diagnostics */
            {
                KaonLCDAs k(p, Options{ });
                Diagnostics diagnostics = k.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    // Using running alpha_s from eos these values are reproduced by Carolina Bolognani
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94850, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92874, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91708, 1e-5), // c_rge(mu = 4.0 GeV)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* Twist 2 */
            {
                KaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(-0.0525,     k.a1(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL(-0.04350,    k.a1(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL(-0.04037,    k.a1(3.0),   eps);

                TEST_CHECK_NEARLY_EQUAL( 0.106,      k.a2(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.079020,   k.a2(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.070298,   k.a2(3.0),   eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(0.0, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.29931, k.phi(0.3, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.2615,  k.phi(0.5, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.14055, k.phi(0.7, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(1.0, 1.0), eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(0.0, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.29591, k.phi(0.3, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.3222,  k.phi(0.5, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.16435, k.phi(0.7, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(1.0, 2.0), eps);
            }

            /* Twist 3 */
            {
                KaonLCDAs k(p, Options{ });

                // Using alpha_s, m_s_msbar, m_ud_msbar from eos, with omega3 and lambda3 tested as f3K.omega3(mu), reproduced by Carolina Bolognani
                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.00450,   k.f3(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00343,   k.f3(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00308,   k.f3(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 1.86233,   k.mu(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 2.45562,   k.mu(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 2.71996,   k.mu(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01548,   k.eta3(1.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00894,   k.eta3(2.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00725,   k.eta3(3.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-0.00675,   k.omega3(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-0.00361,   k.omega3(2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-0.00281,   k.omega3(3.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00720,   k.lambda3(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00422,   k.lambda3(2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.003364,  k.lambda3(3.0),    eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.05106,   k.phi3p(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.06602,   k.phi3p(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.00465,   k.phi3p(0.3, 1.0),    eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.02952,   k.phi3p(0.1, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.02602,   k.phi3p(0.2, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.98740,   k.phi3p(0.3, 2.0),    10.0 * eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.45052,   k.phi3s(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.88126,   k.phi3s(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.13286,   k.phi3s(0.3, 1.0),    eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.493549,  k.phi3s(0.1, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.909772,  k.phi3s(0.2, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.17754,   k.phi3s(0.3, 2.0),    10.0 * eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.076227,  k.phi3s_d1(0.1, 1.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 3.418348,  k.phi3s_d1(0.2, 1.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.674345,  k.phi3s_d1(0.3, 1.0), 50.0 * eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.869314,  k.phi3s_d1(0.1, 2.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 3.411135,  k.phi3s_d1(0.2, 2.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.982612,  k.phi3s_d1(0.3, 2.0), 50.0 * eps);
            }

            /* Twist 4 */
            {
                KaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.18,       k.delta4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15437,    k.delta4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14544,    k.delta4(3.0),   eps);

                TEST_CHECK_NEARLY_EQUAL(-0.09,       k.kappa4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL(-0.103105,   k.kappa4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL(-0.107073,   k.kappa4(3.0),   eps);

                TEST_CHECK_NEARLY_EQUAL( 0.036,      k.omega4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.021217,   k.omega4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.017189,   k.omega4(3.0),   eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.3663714,  k.phi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.8035105,  k.phi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.1313687,  k.phi4(0.3, 1.0), eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.3254996,  k.phi4(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.7047138,  k.phi4(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.9951251,  k.phi4(0.3, 2.0), eps);

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.7499464,  k.psi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0239059,  k.psi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.4302608,  k.psi4(0.3, 1.0), eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.6285379,  k.psi4(0.1, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0001853,  k.psi4(0.2, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3769193,  k.psi4(0.3, 2.0), 5.0 * eps);
            }
        }
} k_lcdas_test;
