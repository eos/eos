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
                TEST_CHECK_NEARLY_EQUAL(bgl.V1_a0(),  0.009176992,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1_a1(), -0.020874910,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1_a2(), -0.075783486,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1_a0(),  0.048503393,  eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1_a1(), -0.130821879,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1_a2(), -0.169762556,   eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.fT_a0(),  0.020080524,      eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.fT_a1(), -0.023524298,       eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.fT_a2(), -0.179347292,       eps);
                // }}}

                // B -> D^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.A1_a0(),   0.007556149, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1_a1(),   0.003260930, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1_a2(),  -0.063636665, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5_a0(),   0.005061682, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5_a1(),   0.005955652, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5_a2(),  -0.163047053, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4_a0(),   0.007953916, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4_a1(),  -0.011128853, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4_a2(),  -0.079433489, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1_a0(),   0.033632658, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1_a1(),  -0.076948722, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1_a2(),  -0.336758100, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T1_a0(),   0.008684675, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T1_a1(),  -0.014754705, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T1_a2(),  -0.086120077, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T2_a0(),   0.002276531, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T2_a1(),  -0.000143740, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T2_a2(),  -0.019106882, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T23_a0(),  0.006796872, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T23_a1(),  0.009192890, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T23_a2(), -0.116562690, eps);
                // }}}

                // B^* -> D form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.P2_a0(),      0.042093130, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2_a1(),     -0.096687471, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2_a2(),     -0.337554764, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5_a0(),      0.010575624, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5_a1(),     -0.021725498, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5_a2(),     -0.090241526, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2_a0(),      0.013086119, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2_a1(),     -0.004551056, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2_a2(),     -0.128767639, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6_a0(),      0.009481104, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6_a1(),     -0.032699174, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6_a2(),     -0.019488169, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T1bar_a0(),   0.011381095, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T1bar_a1(),  -0.024888540, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T1bar_a2(),  -0.099617764, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T2bar_a0(),  -0.004385960, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T2bar_a1(),   0.007940773, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T2bar_a2(),   0.032473434, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T23bar_a0(), -0.012107282, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T23bar_a1(),  0.001859758, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.T23bar_a2(),  0.127139400, eps);
                // }}}

                // B^* -> D^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.S2_a0(),  0.037167931, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2_a1(), -0.103766362, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2_a2(),  0.026104969, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3_a0(),  0.026281696, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3_a1(), -0.037808594, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3_a2(), -0.561217556, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3_a0(),  0.034444832, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3_a1(), -0.049307203, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3_a2(), -0.353556202, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2_a0(),  0.007118838, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2_a1(), -0.012150318, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2_a2(), -0.073751134, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3_a0(),  0.005436304, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3_a1(), -0.014107590, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3_a2(), -0.059719187, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6_a0(),  0.008882384, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6_a1(), -0.023316688, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6_a2(), -0.102993404, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7_a0(),  0.007758780, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7_a1(), -0.010999774, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7_a2(), -0.084588537, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3_a0(),  0.006063688, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3_a1(),  0.006124483, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3_a2(), -0.059089083, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4_a0(),  0.006063688, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4_a1(),  0.010589772, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4_a2(), -0.108993971, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7_a0(),  0.005793597, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7_a1(), -0.008654110, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7_a2(), -0.033357940, eps);
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
                TEST_CHECK_NEARLY_EQUAL(bgl.V1s_a0(),  0.005000218, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1s_a1(), -0.006466751, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V1s_a2(), -0.023610421, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1s_a0(),  0.028392821, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1s_a1(), -0.067892250, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S1s_a2(),  0.057562073, eps);
                // }}}

                // B_s -> D_s^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.A1s_a0(),  0.003376126, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1s_a1(),  0.010932131, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A1s_a2(), -0.048979296, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5s_a0(),  0.002186707, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5s_a1(),  0.009923116, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A5s_a2(), -0.059844679, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4s_a0(),  0.005143734, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4s_a1(), -0.004447454, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V4s_a2(), -0.075136930, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1s_a0(),  0.018682227, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1s_a1(), -0.002564102, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P1s_a2(), -0.280071084, eps);
                // }}}

                // B_s^* -> D_s form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.P2s_a0(),  0.030056239, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2s_a1(), -0.059384529, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P2s_a2(), -0.159240018, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5s_a0(),  0.005803257, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5s_a1(), -0.007572282, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V5s_a2(), -0.052200227, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2s_a0(),  0.004908143, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2s_a1(),  0.006442527, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A2s_a2(), -0.032997771, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6s_a0(),  0.003439343, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6s_a1(), -0.008766692, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A6s_a2(), -0.011213961, eps);
                // }}}

                // B_s^* -> D_s^* form factors
                // {{{
                TEST_CHECK_NEARLY_EQUAL(bgl.S2s_a0(),  0.024900426, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2s_a1(), -0.015599130, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S2s_a2(), -0.525820382, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3s_a0(),  0.017607260, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3s_a1(),  0.001307785, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.S3s_a2(), -0.384888269, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3s_a0(),  0.020332967, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3s_a1(), -0.000760986, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.P3s_a2(), -0.357220416, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2s_a0(),  0.004464033, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2s_a1(), -0.002213377, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V2s_a2(), -0.079883864, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3s_a0(),  0.003325891, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3s_a1(), -0.002463496, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V3s_a2(), -0.064104200, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6s_a0(),  0.005775325, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6s_a1(), -0.003474061, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V6s_a2(), -0.105448490, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7s_a0(),  0.005284912, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7s_a1(), -0.004319367, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.V7s_a2(), -0.088676357, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3s_a0(),  0.002823945, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3s_a1(),  0.011193429, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A3s_a2(), -0.049381292, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4s_a0(),  0.002823945, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4s_a1(),  0.012224422, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A4s_a2(), -0.047080505, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7s_a0(),  0.002610754, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7s_a1(),  0.004655519, eps);
                TEST_CHECK_NEARLY_EQUAL(bgl.A7s_a2(), -0.069996103, eps);
                // }}}
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
