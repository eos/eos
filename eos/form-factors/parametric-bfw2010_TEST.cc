/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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
#include <eos/form-factors/parametric-bfw2010-impl.hh>

using namespace test;
using namespace eos;

class BToKBFW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKBFW2010FormFactorsTest() :
            TestCase("b_to_k_bfw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K::a^f+_0@BFW2010"]     =  0.01;
                p["B->K::a^f+_1@BFW2010"]     = -0.02;
                p["B->K::a^fT_0@BFW2010"]     =  0.03;
                p["B->K::a^fT_1@BFW2010"]     = -0.04;
                p["B->K::a^f0_1@BFW2010"]     =  0.05;

                p["mass::B_d@BSZ2015"]        =  5.279;
                p["mass::K_d@BSZ2015"]        =  0.492;
                p["mass::B_s@BSZ2015"]        =  5.367;
                p["mass::B_s^*@BSZ2015"]      =  5.416;
                p["mass::B_s,0@BSZ2015"]      =  5.711;
                p["mass::B_s,1@BSZ2015"]      =  5.750;

                // Optimized t0 = (mB + mK) * (sqrt(mB) - sqrt(mK))^2
                p["B->K::t0@BFW2010"]         =  14.703305673;


                BFW2010FormFactors<BToK, PToP> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                for (const auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.144596,  eps), // z(q2 =  0)
                    std::make_pair(  0.0562957, eps), // z(q2 = 10)
                    std::make_pair(  0.398942,  eps), // p_0(z = 0.0)
                    std::make_pair(  0,         eps), // p_1(z = 0.0)
                    std::make_pair(  0.398942,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair(  0.0224588, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0669446, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0639514, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0607295, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0870644, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0874959, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0878445, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.1660948, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.1551728, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.1438773, eps), // phi_f_t(z = z(q2 =  4))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against Nico's implementation
                TEST_CHECK_NEARLY_EQUAL( ff.f_0( -1.0),  0.0916481  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(  1.0),  0.0876659  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(  4.0),  0.081212   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0( 25.0), -0.0220068  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p( -1.0),  0.0855434  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p(  1.0),  0.0941634  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p(  4.0),  0.110051   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p( 25.0),  1.12023    , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t( -1.0),  0.119334   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t(  1.0),  0.131088   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t(  4.0),  0.15294    , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t( 25.0),  1.65961    , eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),       0.00255580,   eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v_0(),     0.0005,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v_perp(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v_para(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),       0.0005 / 3.0, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a_0(),     0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a_perp(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a_para(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t_0(),     0.0025,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t_perp(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t_para(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),       0.0025 / 3.0, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5_0(),    0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5_perp(), 0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5_para(), 0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),      0,            eps);

                // Test everything for tp smaller than the scalar resonance B_s0
                p["B->K::tp@BFW2010"]         =  30.261001;

                diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference2
                {
                    std::make_pair(  0.164809,  eps), // z(q2 =  0)
                    std::make_pair(  0.0659398, eps), // z(q2 = 10)
                    std::make_pair(  0.398942,  eps), // p_0(z = 0.0)
                    std::make_pair(  0,         eps), // p_1(z = 0.0)
                    std::make_pair(  0.398942,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair(  0.0263062, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0709813, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0677758, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0643178, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0962603, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0968854, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0974312, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.1681133, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.1566374, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.1447653, eps), // phi_f_t(z = z(q2 =  4))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference2);

                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff.f_0( -1.0),  0.0570724  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(  1.0),  0.0535204  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(  4.0),  0.0477983  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0( 25.0), -0.0356356  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p( -1.0),  0.0528831  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p(  1.0),  0.0579401  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p(  4.0),  0.0671252  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p( 25.0),  0.54073    , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t( -1.0),  0.0793976  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t(  1.0),  0.0866462  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t(  4.0),  0.099963   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t( 25.0),  0.869138   , eps);
            }
        }
} b_to_k_bfw2010_form_factors_test;

class BToKstarBFW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKstarBFW2010FormFactorsTest() :
            TestCase("b_to_kstar_bfw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::a^V_0@BFW2010"]     =  0.01;
                p["B->K^*::a^V_1@BFW2010"]     = -0.02;
                p["B->K^*::a^A0_0@BFW2010"]    =  0.03;
                p["B->K^*::a^A0_1@BFW2010"]    = -0.04;
                p["B->K^*::a^A1_1@BFW2010"]    =  0.05;
                p["B->K^*::a^A12_1@BFW2010"]   = -0.06;
                p["B->K^*::a^T1_0@BFW2010"]    =  0.07;
                p["B->K^*::a^T1_1@BFW2010"]    = -0.08;
                p["B->K^*::a^T2_1@BFW2010"]    =  0.09;
                p["B->K^*::a^T23_1@BFW2010"]   = -0.10;

                p["mass::B_d@BSZ2015"]        =  5.279;
                p["mass::K_d^*@BSZ2015"]      =  0.896;
                p["mass::B_s@BSZ2015"]        =  5.367;
                p["mass::B_s^*@BSZ2015"]      =  5.416;
                p["mass::B_s,0@BSZ2015"]      =  5.711;
                p["mass::B_s,1@BSZ2015"]      =  5.750;

                // Optimized t0 = (mB + mK*) * (sqrt(mB) - sqrt(mK*))^2
                p["B->K^*::t0@BFW2010"]       =  11.271194912;

                BFW2010FormFactors<BToKstar, PToV> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                for (const auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.109126,  eps), // z_a(q2 =  0)
                    std::make_pair(  0.115965,  eps), // z_v(q2 =  0)
                    std::make_pair(  0.015044,  eps), // z_a(q2 = 10)
                    std::make_pair(  0.016197,  eps), // z_v(q2 = 10)
                    std::make_pair(  0.500293,  eps), // p_0(z = 0.0)
                    std::make_pair( -0.256101,  eps), // p_1(z = 0.0)
                    std::make_pair(  0.324605,  eps), // p_2(z = 0.0)
                    std::make_pair( -0.395358,  eps), // p_3(z = 0.0)
                    std::make_pair(  0.474303,  eps), // p_4(z = 0.0)
                    std::make_pair( -0.565749,  eps), // p_5(z = 0.0)
                    std::make_pair(  0.500293,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.246998,  eps), // p_1(z = z(q2 = 10))
                    std::make_pair(  0.317589,  eps), // p_2(z = z(q2 = 10))
                    std::make_pair( -0.385009,  eps), // p_3(z = z(q2 = 10))
                    std::make_pair(  0.459807,  eps), // p_4(z = z(q2 = 10))
                    std::make_pair( -0.545962,  eps), // p_5(z = z(q2 = 10))

                    std::make_pair(  0.220729, eps), // phi_v(z = z(q2 = -2.0))
                    std::make_pair(  0.208336, eps), // phi_v(z = z(q2 =  1.0))
                    std::make_pair(  0.195468, eps), // phi_v(z = z(q2 =  4.0))
                    std::make_pair(  0.197626, eps), // phi_a_0(z = z(q2 = -2.0))
                    std::make_pair(  0.186676, eps), // phi_a_0(z = z(q2 =  1.0))
                    std::make_pair(  0.175318, eps), // phi_a_0(z = z(q2 =  4.0))
                    std::make_pair(  0.144187, eps), // phi_a_1(z = z(q2 = -2.0))
                    std::make_pair(  0.145708, eps), // phi_a_1(z = z(q2 =  1.0))
                    std::make_pair(  0.147222, eps), // phi_a_1(z = z(q2 =  4.0))
                    std::make_pair(  0.054581, eps), // phi_a_12(z = z(q2 = -2.0))
                    std::make_pair(  0.056459, eps), // phi_a_12(z = z(q2 =  1.0))
                    std::make_pair(  0.058497, eps), // phi_a_12(z = z(q2 =  4.0))
                    std::make_pair(  0.149022, eps), // phi_t_1(z = z(q2 = -2.0))
                    std::make_pair(  0.144143, eps), // phi_t_1(z = z(q2 =  1.0))
                    std::make_pair(  0.138866, eps), // phi_t_1(z = z(q2 =  4.0))
                    std::make_pair(  0.067858, eps), // phi_t_2(z = z(q2 = -2.0))
                    std::make_pair(  0.070194, eps), // phi_t_2(z = z(q2 =  1.0))
                    std::make_pair(  0.072727, eps), // phi_t_2(z = z(q2 =  4.0))
                    std::make_pair(  0.062179, eps), // phi_t_23(z = z(q2 = -2.0))
                    std::make_pair(  0.062835, eps), // phi_t_23(z = z(q2 =  1.0))
                    std::make_pair(  0.063488, eps), // phi_t_23(z = z(q2 =  4.0))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                const double tm = (BToKstar::m_B - BToKstar::m_V) * (BToKstar::m_B - BToKstar::m_V);

                const double factora12a0 = (BToKstar::m_B * BToKstar::m_B - BToKstar::m_V * BToKstar::m_V) / 8.0 / BToKstar::m_B / BToKstar::m_V;
                const double factora12a1 = (BToKstar::m_B + BToKstar::m_V) * (BToKstar::m_B * BToKstar::m_B - BToKstar::m_V * BToKstar::m_V - tm)
                                           / 16.0 / BToKstar::m_B / BToKstar::m_V / BToKstar::m_V;
                const double factort23t2 = (BToKstar::m_B + BToKstar::m_V) * (BToKstar::m_B * BToKstar::m_B + 3.0 * BToKstar::m_V * BToKstar::m_V - tm)
                                           / 8.0 / BToKstar::m_B / BToKstar::m_V / BToKstar::m_V;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(tm) , factora12a1 * ff.a_1(tm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(tm) , factort23t2 * ff.t_2(tm) , eps);

                // Test against Nico's implementation
                TEST_CHECK_NEARLY_EQUAL( ff.v   ( -1.0),  0.057080, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.v   (  1.0),  0.061331, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.v   (  4.0),  0.068976, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.v   ( 25.0),  0.397226, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 ( -1.0),  0.197014, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 (  1.0),  0.212504, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 (  4.0),  0.2407  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 ( 25.0),  1.80208 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 ( -1.0),  0.360964, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 (  1.0),  0.355625, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 (  4.0),  0.347392, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 ( 25.0),  0.276115, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12( -1.0),  0.143508, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(  1.0),  0.149035, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(  4.0),  0.157662, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12( 25.0),  0.252678, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 ( -1.0),  0.479676, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 (  1.0),  0.503824, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 (  4.0),  0.546488, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 ( 25.0),  2.092917, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 ( -1.0),  0.502260, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 (  1.0),  0.480413, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 (  4.0),  0.447007, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 ( 25.0),  0.161138, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23( -1.0),  0.345526, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(  1.0),  0.356229, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(  4.0),  0.373785, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23( 25.0),  0.641210, eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),       0.0025,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v_0(),     0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v_perp(),  0.0005,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v_para(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),       0.0005 / 3.0, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a_0(),     0.0036050,    eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a_perp(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a_para(),  0.0025,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),       0.0020350,    eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t_0(),     0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t_perp(),  0.0113,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t_para(),  0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),       0.0113 / 3.0, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5_0(),    0.01,         eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5_perp(), 0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5_para(), 0.0178789,    eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),      0.0092929,    eps);
            }
        }
} b_to_kstar_bfw2010_form_factors_test;
