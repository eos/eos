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

class AntiKaonLCDAsTest :
    public TestCase
{
    public:
        AntiKaonLCDAsTest() :
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
            p["K::omega3@1GeV"]      = -1.5;
            p["K::lambda3@1GeV"]     =  1.6;
            p["K::delta4@1GeV"]     =  0.18;
            p["K::kappa4@1GeV"]      = -0.09;
            p["K::omega4@1GeV"]      =  0.2;
            p["mass::K_u"]           =  0.49368;
            p["decay-constant::K_u"] =  0.1561;

            /* Diagnostics */
            {
                AntiKaonLCDAs k(p, Options{ });
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
                AntiKaonLCDAs k(p, Options{ });

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
                AntiKaonLCDAs k(p, Options{ });

                // Using alpha_s, m_s_msbar, m_ud_msbar from eos, with omega3 and lambda3 tested as f3K.omega3(mu), reproduced by Carolina Bolognani
                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.00450,   k.f3(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00344,   k.f3(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00309,   k.f3(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 1.86254,   k.mu3(1.0),       eps);
                TEST_CHECK_NEARLY_EQUAL( 2.45562,   k.mu3(2.0),       eps);
                TEST_CHECK_NEARLY_EQUAL( 2.71990,   k.mu3(3.0),       eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01548,   k.eta3(1.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00897,   k.eta3(2.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00728,   k.eta3(3.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-1.5,       k.omega3(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-1.05546,   k.omega3(2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-0.91576,   k.omega3(3.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.6,       k.lambda3(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.23634,   k.lambda3(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.10056,   k.lambda3(3.0),   eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.04067,   k.phi3p(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.06101,   k.phi3p(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.00200,   k.phi3p(0.3, 1.0),    eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.02362,   k.phi3p(0.1, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.02319,   k.phi3p(0.2, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.98585,   k.phi3p(0.3, 2.0),    10.0 * eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.439681,   k.phi3s(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.870671,   k.phi3s(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.125266,   k.phi3s(0.3, 1.0),    eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.487512,  k.phi3s(0.1, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.903761,  k.phi3s(0.2, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.173090,   k.phi3s(0.3, 2.0),    10.0 * eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.052666,  k.phi3s_d1(0.1, 1.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 3.439031,  k.phi3s_d1(0.2, 1.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.710885,  k.phi3s_d1(0.3, 1.0), 50.0 * eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.854949,  k.phi3s_d1(0.1, 2.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 3.421594,  k.phi3s_d1(0.2, 2.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 2.002037,  k.phi3s_d1(0.3, 2.0), 50.0 * eps);
            }

            /* Twist 4 */
            {
                AntiKaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.18,       k.delta4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15437,    k.delta4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14544,    k.delta4(3.0),   eps);

                TEST_CHECK_NEARLY_EQUAL(-0.09,       k.kappa4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL(-0.103105,   k.kappa4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL(-0.107073,   k.kappa4(3.0),   eps);

                TEST_CHECK_NEARLY_EQUAL( 0.2,        k.omega4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.137439,   k.omega4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.118188,   k.omega4(3.0),   eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.3663714,  k.phi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.8035002,  k.phi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.1313548,  k.phi4(0.3, 1.0), eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.3255290,  k.phi4(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.7048003,  k.phi4(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.9952677,  k.phi4(0.3, 2.0), eps);

                // phi4_d1, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.4619481,  k.phi4_d1(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 3.9998219,  k.phi4_d1(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.4640627,  k.phi4_d1(0.3, 1.0), eps);

                // phi4_d1, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 3.8610347,  k.phi4_d1(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 3.4987248,  k.phi4_d1(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.2209078,  k.phi4_d1(0.3, 2.0), eps);

                // phi4_d2, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.1756668,  k.phi4_d2(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-11.852937,  k.phi4_d2(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-17.533827,  k.phi4_d2(0.3, 1.0), eps);

                // phi4_d2, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.0958549,  k.phi4_d2(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-9.5516860,  k.phi4_d2(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-14.990626,  k.phi4_d2(0.3, 2.0), eps);

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.7608706,  k.psi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0317497,  k.psi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.4249862,  k.psi4(0.3, 1.0), eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.6338835,  k.psi4(0.1, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0037173,  k.psi4(0.2, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3743138,  k.psi4(0.3, 2.0), 5.0 * eps);

                // psi4_i, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.11623111, k.psi4_i(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15368080, k.psi4_i(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.13175167, k.psi4_i(0.3, 1.0), eps);

                // psi4_i, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.1019395,  k.psi4_i(0.1, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1317095,  k.psi4_i(0.2, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1111354,  k.psi4_i(0.3, 2.0), 5.0 * eps);
            }
        }
} anti_k_lcdas_test;

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
            p["mass::s(2GeV)"]       =  0.5 * 0.0085;
            p["mass::u(2GeV)"]       =  2.0 * 0.095;  // we use 8.5 MeV for twice the average u/d quark mass
            p["mass::d(2GeV)"]       =  0.0;
            p["K::a1@1GeV"]          =  0.0525;
            p["K::a2@1GeV"]          =  0.106;
            p["K::f3@1GeV"]          =  0.0045;
            p["K::omega3@1GeV"]      = -1.5;
            p["K::lambda3@1GeV"]     = -1.6;
            p["K::delta4@1GeV"]     =  0.18;
            p["K::kappa4@1GeV"]      =  0.09;
            p["K::omega4@1GeV"]      =  0.2;
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
                const double eps = 1.0e-5;
                KaonLCDAs k(p, Options{ });

                // Using alpha_s, m_s_msbar, m_ud_msbar from eos, with omega3 and lambda3 tested as f3K.omega3(mu), reproduced by Carolina Bolognani
                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.00450,   k.f3(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00344,   k.f3(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00309,   k.f3(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 1.86254,   k.mu3(1.0),       eps);
                TEST_CHECK_NEARLY_EQUAL( 2.45562,   k.mu3(2.0),       eps);
                TEST_CHECK_NEARLY_EQUAL( 2.71990,   k.mu3(3.0),       eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01548,   k.eta3(1.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00897,   k.eta3(2.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00728,   k.eta3(3.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-1.5,       k.omega3(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-1.05546,   k.omega3(2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-0.91576,   k.omega3(3.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.6,       k.lambda3(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.23634,   k.lambda3(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 1.10056,   k.lambda3(3.0),   eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.04067,   k.phi3p(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.06101,   k.phi3p(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.00200,   k.phi3p(0.3, 1.0),    eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.02362,   k.phi3p(0.1, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.02319,   k.phi3p(0.2, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.98585,   k.phi3p(0.3, 2.0),    10.0 * eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.439681,   k.phi3s(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.870671,   k.phi3s(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.125266,   k.phi3s(0.3, 1.0),    eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.487512,  k.phi3s(0.1, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.903761,  k.phi3s(0.2, 2.0),    10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.173090,   k.phi3s(0.3, 2.0),    10.0 * eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.052666,  k.phi3s_d1(0.1, 1.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 3.439031,  k.phi3s_d1(0.2, 1.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 1.710885,  k.phi3s_d1(0.3, 1.0), 50.0 * eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.854949,  k.phi3s_d1(0.1, 2.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 3.421594,  k.phi3s_d1(0.2, 2.0), 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 2.002037,  k.phi3s_d1(0.3, 2.0), 50.0 * eps);
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

                TEST_CHECK_NEARLY_EQUAL( 0.2,        k.omega4(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.137439,   k.omega4(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.118188,   k.omega4(3.0),   eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.3663714,  k.phi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.8035002,  k.phi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.1313548,  k.phi4(0.3, 1.0), eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.3255290,  k.phi4(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.7048003,  k.phi4(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.9952677,  k.phi4(0.3, 2.0), eps);

                // phi4_d1, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.4619481,  k.phi4_d1(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 3.9998219,  k.phi4_d1(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.4640627,  k.phi4_d1(0.3, 1.0), eps);

                // phi4_d1, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 3.8610347,  k.phi4_d1(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 3.4987248,  k.phi4_d1(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.2209078,  k.phi4_d1(0.3, 2.0), eps);

                // phi4_d2, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.1756668,  k.phi4_d2(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-11.852937,  k.phi4_d2(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-17.533827,  k.phi4_d2(0.3, 1.0), eps);

                // phi4_d2, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.0958549,  k.phi4_d2(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-9.5516860,  k.phi4_d2(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-14.990626,  k.phi4_d2(0.3, 2.0), eps);

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.7608706,  k.psi4(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0317497,  k.psi4(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.4249862,  k.psi4(0.3, 1.0), eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.6338835,  k.psi4(0.1, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0037173,  k.psi4(0.2, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3743138,  k.psi4(0.3, 2.0), 5.0 * eps);

                // psi4_i, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.11623111, k.psi4_i(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15368080, k.psi4_i(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.13175167, k.psi4_i(0.3, 1.0), eps);

                // psi4_i, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.1019395,  k.psi4_i(0.1, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1317095,  k.psi4_i(0.2, 2.0), 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1111354,  k.psi4_i(0.3, 2.0), 5.0 * eps);
            }
        }
} k_lcdas_test;
