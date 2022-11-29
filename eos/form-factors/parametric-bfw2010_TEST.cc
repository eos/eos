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

                    std::make_pair(  0.0386505, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0369224, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0350622, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0870644, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0874959, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0878445, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.0958949, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.0895891, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.0830676, eps), // phi_f_t(z = z(q2 =  4))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against Nico's implementation
                TEST_CHECK_NEARLY_EQUAL( ff.f_0( -1.0),  0.157128  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(  1.0),  0.153515  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(  4.0),  0.147786  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_0( 25.0),  0.0736626 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p( -1.0),  0.148165  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p(  1.0),  0.163096  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p(  4.0),  0.190614  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_p( 25.0),  1.9403    , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t( -1.0),  0.206692  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t(  1.0),  0.227051  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t(  4.0),  0.2649    , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_t( 25.0),  2.87453   , eps);
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
                    std::make_pair( 0.0972622, eps), // z_a(q2 =  0)
                    std::make_pair(  0.102919, eps), // z_v(q2 =  0)
                    std::make_pair( 0.0131099, eps), // z_a(q2 = 10)
                    std::make_pair( 0.0140219, eps), // z_v(q2 = 10)
                    std::make_pair(  0.469725, eps), // p_0(z = 0.0)
                    std::make_pair( -0.169168, eps), // p_1(z = 0.0)
                    std::make_pair(  0.201915, eps), // p_2(z = 0.0)
                    std::make_pair( -0.231682, eps), // p_3(z = 0.0)
                    std::make_pair(  0.260728, eps), // p_4(z = 0.0)
                    std::make_pair( -0.290629, eps), // p_5(z = 0.0)
                    std::make_pair(  0.469725, eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.162167, eps), // p_1(z = z(q2 = 10))
                    std::make_pair(  0.198502, eps), // p_2(z = z(q2 = 10))
                    std::make_pair( -0.227156, eps), // p_3(z = z(q2 = 10))
                    std::make_pair(  0.255026, eps), // p_4(z = z(q2 = 10))
                    std::make_pair( -0.283591, eps), // p_5(z = z(q2 = 10))

                    std::make_pair(  0.123474, eps), // phi_v(z = z(q2 = -2.0))
                    std::make_pair(  0.116704, eps), // phi_v(z = z(q2 =  1.0))
                    std::make_pair(  0.109688, eps), // phi_v(z = z(q2 =  4.0))
                    std::make_pair(  0.190961, eps), // phi_a_0(z = z(q2 = -2.0))
                    std::make_pair(  0.180579, eps), // phi_a_0(z = z(q2 =  1.0))
                    std::make_pair(  0.169826, eps), // phi_a_0(z = z(q2 =  4.0))
                    std::make_pair( 0.0750053, eps), // phi_a_1(z = z(q2 = -2.0))
                    std::make_pair( 0.0756391, eps), // phi_a_1(z = z(q2 =  1.0))
                    std::make_pair( 0.0762512, eps), // phi_a_1(z = z(q2 =  4.0))
                    std::make_pair( 0.0271313, eps), // phi_a_12(z = z(q2 = -2.0))
                    std::make_pair( 0.0279491, eps), // phi_a_12(z = z(q2 =  1.0))
                    std::make_pair( 0.0288233, eps), // phi_a_12(z = z(q2 =  4.0))
                    std::make_pair( 0.0795766, eps), // phi_t_1(z = z(q2 = -2.0))
                    std::make_pair( 0.0769082, eps), // phi_t_1(z = z(q2 =  1.0))
                    std::make_pair( 0.0740315, eps), // phi_t_1(z = z(q2 =  4.0))
                    std::make_pair( 0.0337313, eps), // phi_t_2(z = z(q2 = -2.0))
                    std::make_pair(  0.034748, eps), // phi_t_2(z = z(q2 =  1.0))
                    std::make_pair( 0.0358349, eps), // phi_t_2(z = z(q2 =  4.0))
                    std::make_pair( 0.0323453, eps), // phi_t_23(z = z(q2 = -2.0))
                    std::make_pair( 0.0326186, eps), // phi_t_23(z = z(q2 =  1.0))
                    std::make_pair( 0.0328826, eps), // phi_t_23(z = z(q2 =  4.0))
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
                TEST_CHECK_NEARLY_EQUAL( ff.v   ( -1.0),  0.116949 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.v   (  1.0),  0.127039 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.v   (  4.0),  0.145481 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.v   ( 25.0),  1.13972  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 ( -1.0),  0.213511 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 (  1.0),  0.231922 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 (  4.0),  0.265746 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_0 ( 25.0),  2.35769  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 ( -1.0),  0.663848 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 (  1.0),  0.662256 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 (  4.0),  0.660755 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_1 ( 25.0),  0.780837 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12( -1.0),  0.148145 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(  1.0),  0.170255 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(  4.0),  0.20602  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12( 25.0),  0.783297 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 ( -1.0),  1.0843   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 (  1.0),  1.15219  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 (  4.0),  1.27434  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1 ( 25.0),  6.98411  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 ( -1.0),  1.13435  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 (  1.0),  1.09969  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 (  4.0),  1.04678  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_2 ( 25.0),  0.595764 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23( -1.0),  0.74392  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(  1.0),  0.782935 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(  4.0),  0.849064 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23( 25.0),  2.27539  , eps);
            }
        }
} b_to_kstar_bfw2010_form_factors_test;
