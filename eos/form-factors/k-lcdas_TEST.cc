/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
 * Copyright (c) 2022 Carolina Bolognani
 * Copyright (c) 2024 Stefan Meiser
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
            p["K::delta4@1GeV"]      =  0.18;
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
                TEST_CHECK_NEARLY_EQUAL(k.a1(1.0), -0.0525,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.a1(2.0), -0.04350,     eps);
                TEST_CHECK_NEARLY_EQUAL(k.a1(3.0), -0.04037,     eps);

                TEST_CHECK_NEARLY_EQUAL(k.a2(1.0),  0.106,       eps);
                TEST_CHECK_NEARLY_EQUAL(k.a2(2.0),  0.079020,    eps);
                TEST_CHECK_NEARLY_EQUAL(k.a2(3.0),  0.070298,    eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.3, 1.0), 1.29931,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.5, 1.0), 1.2615,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.7, 1.0), 1.14055,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(1.0, 1.0), 0.0,      eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.3, 2.0), 1.29591,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.5, 2.0), 1.3222,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.7, 2.0), 1.16435,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(1.0, 2.0), 0.0,      eps);
            }

            /* Twist 3 */
            {
                AntiKaonLCDAs k(p, Options{ });

                // Using alpha_s, m_s_msbar, m_ud_msbar from eos, with omega3 and lambda3 tested as f3K.omega3(mu), reproduced by Carolina Bolognani
                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.f3(1.0),       0.00450,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.f3(2.0),       0.00344,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.f3(3.0),       0.00309,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.mu3(1.0),      1.86254,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.mu3(2.0),      2.45562,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.mu3(3.0),      2.71990,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.eta3(1.0),     0.01548,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.eta3(2.0),     0.00897,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.eta3(3.0),     0.00728,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega3(1.0),  -1.5,       eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega3(2.0),  -1.05546,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega3(3.0),  -0.91576,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.lambda3(1.0),  1.6,       eps);
                TEST_CHECK_NEARLY_EQUAL(k.lambda3(2.0),  1.23634,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.lambda3(3.0),  1.10056,   eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.1, 1.0), 1.04067, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.2, 1.0), 1.06101, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.3, 1.0), 1.00200, eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.1, 2.0), 1.02362,  10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.2, 2.0), 1.02319,  10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.3, 2.0), 0.98585,  10.0 * eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.1, 1.0), 0.439681, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.2, 1.0), 0.870671, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.3, 1.0), 1.125266, eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.1, 2.0), 0.487512, 10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.2, 2.0), 0.903761, 10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.3, 2.0), 1.173090, 10.0 * eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.1, 1.0), 5.052666, 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.2, 1.0), 3.439031, 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.3, 1.0), 1.710885, 50.0 * eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.1, 2.0), 4.854949, 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.2, 2.0), 3.421594, 50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.3, 2.0), 2.002037, 50.0 * eps);
            }

            /* Twist 4 */
            {
                AntiKaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.delta4(1.0),  0.18,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.delta4(2.0),  0.15437,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.delta4(3.0),  0.14544,   eps);

                TEST_CHECK_NEARLY_EQUAL(k.kappa4(1.0), -0.070363,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.kappa4(2.0), -0.083466,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.kappa4(3.0), -0.087433,  eps);

                TEST_CHECK_NEARLY_EQUAL(k.omega4(1.0),  0.2,       eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega4(2.0),  0.137439,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega4(3.0),  0.118188,  eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.1, 1.0), 0.3546348,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.2, 1.0), 0.7872947,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.3, 1.0), 1.1168668,  eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.1, 2.0), 0.3137956,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.2, 2.0), 0.6885918,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.3, 2.0), 0.9807768,  eps);

                // phi4_d1, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.1, 1.0),  4.3817199,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.2, 1.0),  3.9886215,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.3, 1.0),  2.5065593,  eps);

                // phi4_d1, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.1, 2.0),  3.7808052,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.2, 2.0),  3.4875212,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.3, 2.0),  2.2634101,  eps);

                // phi4_d2, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.1, 1.0),  5.9155005,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.2, 1.0), -11.228389,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.3, 1.0), -17.091844,  eps);

                // phi4_d2, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.1, 2.0),  4.8356851,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.2, 2.0), -8.9271334,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.3, 2.0), -14.548637,  eps);

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.1, 1.0),  0.8047785210,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.2, 1.0),  0.0753944603,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.3, 1.0), -0.3907843998,  eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.1, 2.0),  0.6777910058,  5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.2, 2.0),  0.0473704719,  5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.3, 2.0), -0.3401055477,  5.0 * eps);

                // psi4_i, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.1, 1.0), 0.120139861, eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.2, 1.0), 0.162056387, eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.3, 1.0), 0.144083365, eps);

                // psi4_i, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.1, 2.0), 0.105845125, 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.2, 2.0), 0.140082543, 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.3, 2.0), 0.123465300, 5.0 * eps);
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
            p["K::delta4@1GeV"]      =  0.18;
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
                TEST_CHECK_NEARLY_EQUAL(k.a1(1.0), -0.0525,    eps);
                TEST_CHECK_NEARLY_EQUAL(k.a1(2.0), -0.04350,   eps);
                TEST_CHECK_NEARLY_EQUAL(k.a1(3.0), -0.04037,   eps);

                TEST_CHECK_NEARLY_EQUAL(k.a2(1.0),  0.106,     eps);
                TEST_CHECK_NEARLY_EQUAL(k.a2(2.0),  0.079020,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.a2(3.0),  0.070298,  eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.0, 1.0), 0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.3, 1.0), 1.29931, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.5, 1.0), 1.2615,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.7, 1.0), 1.14055, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(1.0, 1.0), 0.0,     eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.0, 2.0), 0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.3, 2.0), 1.29591, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.5, 2.0), 1.3222,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(0.7, 2.0), 1.16435, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi(1.0, 2.0), 0.0,     eps);
            }

            /* Twist 3 */
            {
                const double eps = 1.0e-5;
                KaonLCDAs k(p, Options{ });

                // Using alpha_s, m_s_msbar, m_ud_msbar from eos, with omega3 and lambda3 tested as f3K.omega3(mu), reproduced by Carolina Bolognani
                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.f3(1.0),       0.00450,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.f3(2.0),       0.00344,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.f3(3.0),       0.00309,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.mu3(1.0),      1.86254,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.mu3(2.0),      2.45562,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.mu3(3.0),      2.71990,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.eta3(1.0),     0.01548,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.eta3(2.0),     0.00897,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.eta3(3.0),     0.00728,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega3(1.0),  -1.5,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega3(2.0),  -1.05546,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega3(3.0),  -0.91576,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.lambda3(1.0),  1.6,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.lambda3(2.0),  1.23634,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.lambda3(3.0),  1.10056,  eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.1, 1.0), 1.04067,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.2, 1.0), 1.06101,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.3, 1.0), 1.00200,  eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.1, 2.0), 1.02362,  10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.2, 2.0), 1.02319,  10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3p(0.3, 2.0), 0.98585,  10.0 * eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.1, 1.0), 0.439681, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.2, 1.0), 0.870671, eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.3, 1.0), 1.125266, eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.1, 2.0), 0.487512, 10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.2, 2.0), 0.903761, 10.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s(0.3, 2.0), 1.173090, 10.0 * eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.1, 1.0), 5.052666,  50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.2, 1.0), 3.439031,  50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.3, 1.0), 1.710885,  50.0 * eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.1, 2.0), 4.854949,  50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.2, 2.0), 3.421594,  50.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi3s_d1(0.3, 2.0), 2.002037,  50.0 * eps);
            }

            /* Twist 4 */
            {
                KaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.delta4(1.0),  0.18,     eps);
                TEST_CHECK_NEARLY_EQUAL(k.delta4(2.0),  0.15437,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.delta4(3.0),  0.14544,  eps);

                TEST_CHECK_NEARLY_EQUAL(k.kappa4(1.0), -0.070363, eps);
                TEST_CHECK_NEARLY_EQUAL(k.kappa4(2.0), -0.083466, eps);
                TEST_CHECK_NEARLY_EQUAL(k.kappa4(3.0), -0.087433, eps);

                TEST_CHECK_NEARLY_EQUAL(k.omega4(1.0),  0.2,      eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega4(2.0),  0.137439, eps);
                TEST_CHECK_NEARLY_EQUAL(k.omega4(3.0),  0.118188, eps);

                // phi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.1, 1.0), 0.3546348,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.2, 1.0), 0.7872947,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.3, 1.0), 1.1168668,  eps);

                // phi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.1, 2.0), 0.3137956,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.2, 2.0), 0.6885918,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4(0.3, 2.0), 0.9807768,  eps);

                // phi4_d1, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.1, 1.0),  4.3817199,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.2, 1.0),  3.9886215,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.3, 1.0),  2.5065593,  eps);

                // phi4_d1, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.1, 2.0),  3.7808052,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.2, 2.0),  3.4875212,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d1(0.3, 2.0),  2.2634101,  eps);

                // phi4_d2, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.1, 1.0),  5.9155005,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.2, 1.0), -11.228389,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.3, 1.0), -17.091844,  eps);

                // phi4_d2, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.1, 2.0),  4.8356851,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.2, 2.0), -8.9271334,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.phi4_d2(0.3, 2.0), -14.548637,  eps);

                // psi4, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.1, 1.0),  0.8047785210,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.2, 1.0),  0.0753944603,  eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.3, 1.0), -0.3907843998,  eps);

                // psi4, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.1, 2.0),  0.6777910058,  5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.2, 2.0),  0.0473704719,  5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4(0.3, 2.0), -0.3401055477,  5.0 * eps);

                // psi4_i, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.1, 1.0), 0.120139861, eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.2, 1.0), 0.162056387, eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.3, 1.0), 0.144083365, eps);

                // psi4_i, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.1, 2.0), 0.105845125, 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.2, 2.0), 0.140082543, 5.0 * eps);
                TEST_CHECK_NEARLY_EQUAL(k.psi4_i(0.3, 2.0), 0.123465300, 5.0 * eps);
            }
        }
} k_lcdas_test;
