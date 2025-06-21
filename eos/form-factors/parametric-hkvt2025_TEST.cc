/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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
#include <eos/form-factors/parametric-hkvt2025-impl.hh>

using namespace test;
using namespace eos;

class BToPiPiHKVT2025FormFactorsTest :
    public TestCase
{
    public:
        BToPiPiHKVT2025FormFactorsTest() :
            TestCase("b_to_pipi_hkvt2025_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->pipi::chi_1p_a@HKvT2025"]  = 5.742e-4;
                p["B->pipi::chi_1m_v@HKvT2025"]  = 5.742e-4;
                p["B->pipi::tp_a@HKvT2025"]  = 29.8;
                p["B->pipi::tp_v@HKvT2025"]  = 29.8;
                p["B->pipi::t0@HKvT2025"]    = 0.0;
                p["B->pipi::sin0@HKvT2025"]   = 0.976144;
                p["B->pipi::sin1@HKvT2025"]   = 0.81;
                p["B->pipi::s0@HKvT2025"]    = 0.0729;

                p["mass::pi^+@HKvT2025"]     = 0.13957;
                p["mass::B_u@HKvT2025"]      = 5.279;

                p["pipi->pipi::Gamman0_K@HKvT2025"] =  5.78022130e-01;
                p["B->pipi::a^F1_1_1_0_0@HKvT2025"]   = -1.66817478e-02;
                p["B->pipi::a^F1_1_1_1_0@HKvT2025"]   = -4.85913680e-02;
                p["B->pipi::a^F1_1_1_0_1@HKvT2025"]   =  1.12230134e-04;
                p["B->pipi::a^f_1_1_0_0@HKvT2025"]    = -3.99693932e-02;
                p["B->pipi::a^f_1_1_1_0@HKvT2025"]    = -1.46214509e-02;
                p["B->pipi::a^f_1_1_0_1@HKvT2025"]    = -1.28391535e-02;
                p["B->pipi::a^g_1_1_0_0@HKvT2025"]    =  2.39237676e-02;
                p["B->pipi::a^g_1_1_1_0@HKvT2025"]    =  2.71707816e-02;
                p["B->pipi::a^g_1_1_0_1@HKvT2025"]    =  2.31612080e-04;
                p["B->pipi::a^F1_0_2_0_0@HKvT2025"]   =  6.59969573e-03;
                p["B->pipi::a^F1_0_2_1_0@HKvT2025"]   =  1.63960580e-02;
                p["B->pipi::a^F1_0_2_0_1@HKvT2025"]   = -1.15792204e-03;
                p["B->pipi::a^f_0_2_0_0@HKvT2025"]    =  8.43339186e-03;
                p["B->pipi::a^f_0_2_1_0@HKvT2025"]    =  8.56825464e-03;
                p["B->pipi::a^f_0_2_0_1@HKvT2025"]    = -3.41325494e-03;
                p["B->pipi::a^g_0_2_0_0@HKvT2025"]    =  3.33162140e-02;
                p["B->pipi::a^g_0_2_1_0@HKvT2025"]    =  4.38104301e-02;
                p["B->pipi::a^g_0_2_0_1@HKvT2025"]    = -3.15353047e-03;
                p["B->pipi::a^F1_0_0_0_0@HKvT2025"]   = -6.13947145e-03;
                p["B->pipi::a^F1_0_0_1_0@HKvT2025"]   = -1.69777024e-02;
                p["B->pipi::a^F1_0_0_0_1@HKvT2025"]   = -1.63231083e-02;

                Options oo;
                oo.declare("C"_ok, "+-");
                oo.declare("I"_ok, "0|1");
                HKVT2025FormFactors<BToPiPi, PToPP> ff(p, oo);

                Diagnostics  diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.0        ,  eps), // y(s  =  4*0.135^2)
                    std::make_pair( -0.009364   ,  eps), // y(s  =  0.1)

                    std::make_pair(  0.0        ,  eps), // z(q2 =  0) - axial
                    std::make_pair(  0.0        ,  eps), // z(q2 =  0) - vector
                    std::make_pair( -0.101852   ,  eps), // z(q2 = 10) - axial
                    std::make_pair( -0.101852   ,  eps), // z(q2 = 10) - vector

                    std::make_pair(  0.430391   ,  eps), // p_0(z = 0.0, k2 = 0.1)
                    std::make_pair( -0.069127   ,  eps), // p_1(z = 0.0, k2 = 0.1)
                    std::make_pair(  0.076481   ,  eps), // p_2(z = 0.0, k2 = 0.1)
                    std::make_pair(  0.430391   ,  eps), // p_0(z = z(q2 = 10, k2 = 0.1))
                    std::make_pair( -0.113525   ,  eps), // p_1(z = z(q2 = 10, k2 = 0.1))
                    std::make_pair(  0.089456   ,  eps), // p_2(z = z(q2 = 10, k2 = 0.1))

                    std::make_pair(  1.0        ,  eps), // p_0(y = y(k2 = 0.1), 0)
                    std::make_pair( -0.213216   ,  eps), // p_1(y = y(k2 = 0.1), 0)
                    std::make_pair(  0.045461   ,  eps), // p_2(y = y(k2 = 0.1), 0)

                    std::make_pair(  1.0        ,  eps), // p_0(y = y(k2 = 0.1), 0)
                    std::make_pair( -0.209985   ,  eps), // p_1(y = y(k2 = 0.1), 0)
                    std::make_pair(  0.038999   ,  eps), // p_2(y = y(k2 = 0.1), 0)

                    std::make_pair(  1.0        ,  eps), // p_0(y = y(k2 = 0.1), 0)
                    std::make_pair( -0.206883   ,  eps), // p_1(y = y(k2 = 0.1), 0)
                    std::make_pair(  0.035121   ,  eps), // p_2(y = y(k2 = 0.1), 0)

                    std::make_pair(  0.005054   ,  eps), // phi_g(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)
                    std::make_pair(  0.005123   ,  eps), // phi_g(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)
                    std::make_pair(  0.222888   ,  eps), // phi_g(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)

                    std::make_pair(  0.000443   ,  eps), // phi_f(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)
                    std::make_pair(  0.000460   ,  eps), // phi_f(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)
                    std::make_pair(  0.023202   ,  eps), // phi_f(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)

                    std::make_pair(  1.41e-05   ,  eps), // phi_F1(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)
                    std::make_pair(  1.50e-05   ,  eps), // phi_F1(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)
                    std::make_pair(  0.000449   ,  eps), // phi_F1(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)
                    std::make_pair(  0.001005   ,  eps), // phi_F1(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 0)

                    std::make_pair(  0.004818   ,  eps), // phi_F2(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)
                    std::make_pair(  0.004885   ,  eps), // phi_F2(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)
                    std::make_pair(  0.122691   ,  eps), // phi_F2(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)
                    std::make_pair(  0.000186   ,  eps), // phi_F2(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 0)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(ff.saturation_1p_a(),   0.077122, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.saturation_1m_v(),   0.299065, eps);
            }

        }
} b_to_pipi_hkvt2025_form_factors_test;
