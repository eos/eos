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
                TEST_CHECK_NEARLY_EQUAL(bgl.V1_a0(),  0.009171994425, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1_a1(), -0.02087801786,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1_a2(), -0.07571211328,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1_a0(),  0.04848134891,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1_a1(), -0.1309252213,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1_a2(), -0.1689469063,   eps);
                // }}}

                // B -> D^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.A1_a0(),  0.007552254359, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1_a1(),  0.003247814948, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1_a2(), -0.06359704257,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5_a0(),  0.005059073459, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5_a1(),  0.005920314047, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5_a2(), -0.162949878,    eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4_a0(),  0.007950422174, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4_a1(), -0.01113445173,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4_a2(), -0.0793764426,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1_a0(),  0.03360836283,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1_a1(), -0.07697180543,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1_a2(), -0.3364350931,   eps);
                // }}}

                // B^* -> D form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.P2_a0(),  0.04205833422,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2_a1(), -0.09669309662,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2_a2(), -0.3371268552,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5_a0(),  0.01057138680,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5_a1(), -0.02172979075,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5_a2(), -0.09017432407,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2_a0(),  0.01307983275,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2_a1(), -0.00456493219,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2_a2(), -0.1286945256,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6_a0(),  0.009476549398, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6_a1(), -0.03275163135,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6_a2(), -0.01927000667,  eps);
                // }}}

                // B^* -> D^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.S2_a0(),  0.03714987562,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2_a1(), -0.1038763165,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2_a2(),  0.02671612864,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3_a0(),  0.02626892897,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3_a1(), -0.03788634407,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3_a2(), -0.5607854015,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3_a0(),  0.03441632174,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3_a1(), -0.0493480833,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3_a2(), -0.353166906,    eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2_a0(),  0.007114668099, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2_a1(), -0.01215894776,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2_a2(), -0.07368710736,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3_a0(),  0.005433355962, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3_a1(), -0.0141136924,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3_a2(), -0.05967391347,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6_a0(),  0.00887933675,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6_a1(), -0.02332273196,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6_a2(), -0.10294312583,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7_a0(),  0.0077557331,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7_a1(), -0.01100581798,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7_a2(), -0.084538259042, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3_a0(),  0.006060562096, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3_a1(),  0.006112265975, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3_a2(), -0.05906053745,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4_a0(),  0.006060562096, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4_a1(),  0.01057755554,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4_a2(), -0.1089654253,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7_a0(),  0.005790611152, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7_a1(), -0.008695881626, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7_a2(), -0.03326665319,  eps);
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
                TEST_CHECK_NEARLY_EQUAL(bgl.V1s_a0(), 0.005000217918,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1s_a1(), -0.006466750504, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1s_a2(), -0.02361042361,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1s_a0(), 0.02839282209,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1s_a1(), -0.06789224686,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1s_a2(), 0.05756204815,   eps);
                // }}}

                // B_s -> D_s^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.A1s_a0(), 0.00337612586,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1s_a1(), 0.01093213157,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1s_a2(), -0.04897929728,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5s_a0(), 0.002186707326,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5s_a1(), 0.00992311644,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5s_a2(), -0.05984468038,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4s_a0(), 0.005143733896,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4s_a1(), -0.004447453494, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4s_a2(), -0.07513693178,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1s_a0(), 0.01868222787,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1s_a1(), -0.00256410114,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1s_a2(), -0.2800710952,   eps);
                // }}}

                // B_s^* -> D_s form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.P2s_a0(), 0.03005624013,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2s_a1(), -0.05938452843,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2s_a2(), -0.1592400379,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5s_a0(), 0.005803257378,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5s_a1(), -0.007572282267, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5s_a2(), -0.05220022977,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2s_a0(), 0.004908142725,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2s_a1(), 0.006442527051,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2s_a2(), -0.03299777215,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6s_a0(), 0.003439343211,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6s_a1(), -0.008766690699, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6s_a2(), -0.01121396414,  eps);
                // }}}

                // B_s^* -> D_s^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.S2s_a0(), 0.02490042677,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2s_a1(), -0.01559912616,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2s_a2(), -0.5258204053,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3s_a0(), 0.01760726063,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3s_a1(), 0.001307787675,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3s_a2(), -0.3848882849,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3s_a0(), 0.0203329674,    eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3s_a1(), -0.0007609840494,eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3s_a2(), -0.3572204296,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2s_a0(), 0.004464032737,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2s_a1(), -0.002213376599, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2s_a2(), -0.07988386599,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3s_a0(), 0.003325891417,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3s_a1(), -0.002463495914, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3s_a2(), -0.0641042013,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6s_a0(), 0.005775324706,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6s_a1(), -0.003474061239, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6s_a2(), -0.1054484923,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7s_a0(), 0.005284912023,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7s_a1(), -0.004319366843, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7s_a2(), -0.08867635877,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3s_a0(), 0.00282394552,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3s_a1(), 0.01119342895,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3s_a2(), -0.04938129245,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4s_a0(), 0.00282394552,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4s_a1(), 0.01222442229,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4s_a2(), -0.04708050602,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7s_a0(), 0.002610754138,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7s_a1(), 0.004655520028,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7s_a2(), -0.06999610402,  eps);
                // }}}

                //throw std::string("foo");
            }
        }
} unitarity_bounds_test;

class OPEUnitarityBoundsTest :
    public TestCase
{
    public:
        OPEUnitarityBoundsTest() :
            TestCase("ope_unitarity_bounds_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-4;

            // values
            {
                Parameters p = Parameters::Defaults();

                OPEUnitarityBounds bounds(p, Options{ });

                TEST_CHECK_RELATIVE_ERROR(bounds.bound_1m(),  5.693e-04,  eps);
                TEST_CHECK_RELATIVE_ERROR(bounds.bound_1p(),  3.3094e-04, eps);
                TEST_CHECK_RELATIVE_ERROR(bounds.bound_0m(), 24.928e-03,  eps);
                TEST_CHECK_RELATIVE_ERROR(bounds.bound_0p(),  4.6391e-03, eps);
            }
        }
} ope_unitarity_bounds_test;
