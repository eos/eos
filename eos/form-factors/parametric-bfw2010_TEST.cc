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
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.144596,  eps), // z(q2 =  0)
                    std::make_pair(  0.0562957, eps), // z(q2 = 10)
                    std::make_pair(  0.398942,  eps), // p_0(z = 0.0)
                    std::make_pair(  0,         eps), // p_1(z = 0.0)
                    std::make_pair(  0.398942,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair(  0.0224588, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0386505, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0369224, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0350622, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0870644, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0874959, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0878445, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.0958949, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.0895891, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.0830676, eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  0.018232,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against Nico's implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( -1.0), 0.157128,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  1.0), 0.153515,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  4.0), 0.147786,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.0736626, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( -1.0), 0.148165,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  1.0), 0.163096,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  4.0), 0.190614,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 1.9403,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( -1.0), 0.206692,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  1.0), 0.227051,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  4.0), 0.2649,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 2.87453,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),       0.00283240578, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),       0,             eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),       0.0005,        eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),       0,             eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),       0.0025,        eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),      0,             eps);

                // Test everything for tp smaller than the scalar resonance B_s0
                p["B->K::tp@BFW2010"]         =  30.261001;

                BFW2010FormFactors<BToK, PToP> ff2(p, Options{ });

                diagnostics = ff2.diagnostics();
                static const std::vector<std::pair<double, double>> reference2
                {
                    std::make_pair(  0.164809,  eps), // z(q2 =  0)
                    std::make_pair(  0.0659398, eps), // z(q2 = 10)
                    std::make_pair(  0.465369,  eps), // p_0(z = 0.0)
                    std::make_pair( -0.1574336, eps), // p_1(z = 0.0)
                    std::make_pair(  0.465369,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.1250388, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0409811, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0391304, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0371339, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0962603, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0968854, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0974312, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.0970603, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.0904346, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.0835803, eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  0.0542382,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference2);

                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(0.0),           ff2.f_p(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff2.f_0( -1.0), 0.224009, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(  1.0), 0.218985, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(  4.0), 0.211081, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0( 25.0), 0.120399, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p( -1.0), 0.214148, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p(  1.0), 0.229426, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p(  4.0), 0.257016, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p( 25.0), 1.59695,  eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t( -1.0), 0.251754, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t(  1.0), 0.2723,   eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t(  4.0), 0.309838, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t( 25.0), 2.36706,  eps);
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
                p["B->K^*::tp_v@BFW2010"]     =  30.261001;
                p["B->K^*::tp_a@BFW2010"]     =  31.764496;

                BFW2010FormFactors<BToKstar, PToV> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
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

                    std::make_pair(  0.127438, eps), // phi_v(z = z(q2 = -2.0))
                    std::make_pair(  0.120283, eps), // phi_v(z = z(q2 =  1.0))
                    std::make_pair(  0.112854, eps), // phi_v(z = z(q2 =  4.0))
                    std::make_pair(  0.197626, eps), // phi_a_0(z = z(q2 = -2.0))
                    std::make_pair(  0.186676, eps), // phi_a_0(z = z(q2 =  1.0))
                    std::make_pair(  0.175318, eps), // phi_a_0(z = z(q2 =  4.0))
                    std::make_pair(  0.083246, eps), // phi_a_1(z = z(q2 = -2.0))
                    std::make_pair(  0.084125, eps), // phi_a_1(z = z(q2 =  1.0))
                    std::make_pair(  0.084999, eps), // phi_a_1(z = z(q2 =  4.0))
                    std::make_pair(  0.031512, eps), // phi_a_12(z = z(q2 = -2.0))
                    std::make_pair(  0.032597, eps), // phi_a_12(z = z(q2 =  1.0))
                    std::make_pair(  0.033773, eps), // phi_a_12(z = z(q2 =  4.0))
                    std::make_pair(  0.086038, eps), // phi_t_1(z = z(q2 = -2.0))
                    std::make_pair(  0.083221, eps), // phi_t_1(z = z(q2 =  1.0))
                    std::make_pair(  0.080174, eps), // phi_t_1(z = z(q2 =  4.0))
                    std::make_pair(  0.039178, eps), // phi_t_2(z = z(q2 = -2.0))
                    std::make_pair(  0.040526, eps), // phi_t_2(z = z(q2 =  1.0))
                    std::make_pair(  0.041989, eps), // phi_t_2(z = z(q2 =  4.0))
                    std::make_pair(  0.035899, eps), // phi_t_23(z = z(q2 = -2.0))
                    std::make_pair(  0.036278, eps), // phi_t_23(z = z(q2 =  1.0))
                    std::make_pair(  0.036654, eps), // phi_t_23(z = z(q2 =  4.0))

                    std::make_pair(  0.10207 , eps), // a_A1_0
                    std::make_pair( -0.009349, eps), // a_A12_0
                    std::make_pair(  0.098888, eps), // a_T2_0
                    std::make_pair(  0.013503, eps), // a_T23_0
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
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( -1.0),  0.098866,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  1.0),  0.10623,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  4.0),  0.119471,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.688017,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( -1.0),  0.197014,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  1.0),  0.212504,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  4.0),  0.2407,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.80208,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( -1.0),  0.502886,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  1.0),  0.494491,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  4.0),  0.48148,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  0.361849,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( -1.0),  0.140285,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  1.0),  0.152292,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  4.0),  0.170921,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  0.362945,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( -1.0),  0.830824,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  1.0),  0.87265,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  4.0),  0.946545,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  3.62504,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( -1.0),  0.869941,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  1.0),  0.8321,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  4.0),  0.774239,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0),  0.2791,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( -1.0),  0.59847,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  1.0),  0.617008,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  4.0),  0.647415,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  1.11061,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),       0.0025,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),       0.0005,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),       0.0166056,    eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),       0.0113,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),      0.0280612,    eps);
            }
        }
} b_to_kstar_bfw2010_form_factors_test;
