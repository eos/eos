/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Nico Gubernari
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

//#define EOS_GENERATE_TESTS 1

#include <test/test.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/mesonic-hqet.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/form-factors/unitarity-bounds.hh>

#include <vector>

using namespace test;
using namespace eos;

class HQETUnitarityBoundsTest :
    public TestCase
{
    public:
        HQETUnitarityBoundsTest() :
            TestCase("unitarity_bounds_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            // q=u,d only
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.14;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +1.88;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = -3.29;
                p["B(*)->D(*)::xi''''(1)@HQET"]  =  0.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = -0.06;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  =  0.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +0.06;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = +0.04;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = -0.05;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.60;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -0.02;
                p["B(*)->D(*)::eta''(1)@HQET"]   = -0.04;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.12;
                p["B(*)->D(*)::l_1'(1)@HQET"]    = -5.78;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -1.89;
                p["B(*)->D(*)::l_2'(1)@HQET"]    = -3.14;
                p["B(*)->D(*)::l_3(1)@HQET"]     = +0.86;
                p["B(*)->D(*)::l_3'(1)@HQET"]    = +0.06;
                p["B(*)->D(*)::l_4(1)@HQET"]     = -2.02;
                p["B(*)->D(*)::l_4'(1)@HQET"]    = -0.05;
                p["B(*)->D(*)::l_5(1)@HQET"]     = +3.79;
                p["B(*)->D(*)::l_5'(1)@HQET"]    = -1.40;
                p["B(*)->D(*)::l_6(1)@HQET"]     = +3.53;
                p["B(*)->D(*)::l_6'(1)@HQET"]    = +0.04;
                p["B(*)->D(*)::a@HQET"]          =  1.0;

                BGLCoefficients bgl(p, Options{ });

                // B -> D form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL( 0.009171994425, bgl.V1_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.02087801786,  bgl.V1_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.07571211328,  bgl.V1_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.04848134891,  bgl.S1_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1309252213,   bgl.S1_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1689469063,   bgl.S1_a2(), eps);
                // }}}

                // B -> D^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL( 0.007552254359, bgl.A1_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.003247814948, bgl.A1_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.06359704257,  bgl.A1_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.005059073459, bgl.A5_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.005920314047, bgl.A5_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.162949878,    bgl.A5_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.007950422174, bgl.V4_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.01113445173,  bgl.V4_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0793764426,   bgl.V4_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.03360836283,  bgl.P1_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.07697180543,  bgl.P1_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3364350931,   bgl.P1_a2(), eps);
                // }}}

                // B^* -> D form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL( 0.04205833422,  bgl.P2_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.09669309662,  bgl.P2_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3371268552,   bgl.P2_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01057138680,  bgl.V5_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.02172979075,  bgl.V5_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.09017432407,  bgl.V5_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01307983275,  bgl.A2_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.00456493219,  bgl.A2_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1286945256,   bgl.A2_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.009476549398, bgl.A6_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.03275163135,  bgl.A6_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.01927000667,  bgl.A6_a2(), eps);
                // }}}

                // B^* -> D^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL( 0.03714987562,  bgl.S2_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1038763165,   bgl.S2_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.02671612864,  bgl.S2_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.02626892897,  bgl.S3_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.03788634407,  bgl.S3_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.5607854015,   bgl.S3_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.03441632174,  bgl.P3_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0493480833,   bgl.P3_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.353166906,    bgl.P3_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.007114668099, bgl.V2_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.01215894776,  bgl.V2_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.07368710736,  bgl.V2_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.005433355962, bgl.V3_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0141136924,   bgl.V3_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.05967391347,  bgl.V3_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00887933675,  bgl.V6_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.02332273196,  bgl.V6_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.10294312583,  bgl.V6_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0077557331,   bgl.V7_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.01100581798,  bgl.V7_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.084538259042, bgl.V7_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.006060562096, bgl.A3_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.006112265975, bgl.A3_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.05906053745,  bgl.A3_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.006060562096, bgl.A4_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01057755554,  bgl.A4_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1089654253,   bgl.A4_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.005790611152, bgl.A7_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.008695881626, bgl.A7_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.03326665319,  bgl.A7_a2(), eps);
                // }}}
            }

            // q = s
            {
                Parameters p = Parameters::Defaults();
                p["B_s(*)->D_s(*)::xi'(1)@HQET"]     = -1.25;
                p["B_s(*)->D_s(*)::xi''(1)@HQET"]    =  2.2;
                p["B_s(*)->D_s(*)::xi'''(1)@HQET"]   = -2.0;
                p["B_s(*)->D_s(*)::xi''''(1)@HQET"]  =  0.0;
                p["B_s(*)->D_s(*)::chi_2(1)@HQET"]   = -0.07;
                p["B_s(*)->D_s(*)::chi_2'(1)@HQET"]  =  0.007;
                p["B_s(*)->D_s(*)::chi_2''(1)@HQET"] = -0.1;
                p["B_s(*)->D_s(*)::chi_3'(1)@HQET"]  =  0.03;
                p["B_s(*)->D_s(*)::chi_3''(1)@HQET"] =  0.05;
                p["B_s(*)->D_s(*)::eta(1)@HQET"]     =  0.72;
                p["B_s(*)->D_s(*)::eta'(1)@HQET"]    =  0.03;
                p["B_s(*)->D_s(*)::eta''(1)@HQET"]   =  0.1;
                p["B_s(*)->D_s(*)::l_1(1)@HQET"]     =  0.2;
                p["B_s(*)->D_s(*)::l_1'(1)@HQET"]    = -3.0;
                p["B_s(*)->D_s(*)::l_2(1)@HQET"]     = -2.2;
                p["B_s(*)->D_s(*)::l_2'(1)@HQET"]    =  1.0;
                p["B_s(*)->D_s(*)::l_3(1)@HQET"]     =  0.5;
                p["B_s(*)->D_s(*)::l_3'(1)@HQET"]    = -0.2;
                p["B_s(*)->D_s(*)::l_4(1)@HQET"]     = -1.0;
                p["B_s(*)->D_s(*)::l_4'(1)@HQET"]    =  0.3;
                p["B_s(*)->D_s(*)::l_5(1)@HQET"]     =  1.4;
                p["B_s(*)->D_s(*)::l_5'(1)@HQET"]    =  2.0;
                p["B_s(*)->D_s(*)::l_6(1)@HQET"]     =  2.5;
                p["B_s(*)->D_s(*)::l_6'(1)@HQET"]    =  0.3;
                p["B_s(*)->D_s(*)::a@HQET"]          =  1.0;

                BGLCoefficients bgl(p, Options{ });
                // B_s -> D_s form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(0.004999810032, bgl.V1s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.006479848759, bgl.V1s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0235512688, bgl.V1s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.02839021048, bgl.S1s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.06800885888, bgl.S1s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.05818215908, bgl.S1s_a2(), eps);
                // }}}

                // B_s -> D_s^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(0.003376760806, bgl.A1s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.01092541431, bgl.A1s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.04900478381, bgl.A1s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.002187179324, bgl.A5s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.009908394091, bgl.A5s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.05986829205, bgl.A5s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.005143979142, bgl.V4s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.004461420957, bgl.V4s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.07510579372, bgl.V4s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.01867864367, bgl.P1s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.002616399546, bgl.P1s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.279931139, bgl.P1s_a2(), eps);
                // }}}

                // B_s^* -> D_s form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(0.03003225295, bgl.P2s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.05940031439, bgl.P2s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1588341602, bgl.P2s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.005801579333, bgl.V5s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.007579360997, bgl.V5s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.05215419135, bgl.V5s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.004906916467, bgl.A2s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.006432668511, bgl.A2s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.03297738509, bgl.A2s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.003438583196, bgl.A6s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.008792660599, bgl.A6s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.01115177087, bgl.A6s_a2(), eps);
                // }}}

                // B_s^* -> D_s^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(0.02489211964, bgl.S2s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.01569145473, bgl.S2s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.5254139564, bgl.S2s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.0176013866, bgl.S3s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.001243731446, bgl.S3s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3846074579, bgl.S3s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.02032035007, bgl.P3s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0008076552162, bgl.P3s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3569918932, bgl.P3s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.004462280149, bgl.V2s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.002222965926, bgl.V2s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.07984499414, bgl.V2s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.003324691584, bgl.V3s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.002470518091, bgl.V3s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.06407761242, bgl.V3s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.005774428399, bgl.V6s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.003482583946, bgl.V6s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1054201708, bgl.V6s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.005283927977, bgl.V7s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.004327768631, bgl.V7s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0886453271, bgl.V7s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.00282325239, bgl.A3s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.01118697484, bgl.A3s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.04938683782, bgl.A3s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.00282325239, bgl.A4s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.01221827038, bgl.A4s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.04708621016, bgl.A4s_a2(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.002610239861, bgl.A7s_a0(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.004633227622, bgl.A7s_a1(), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.06998978489, bgl.A7s_a2(), eps);
                // }}}

                //throw std::string("foo");
            }
        }
} unitarity_bounds_test;
