/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
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
#include <eos/form-factors/parametric-kmpw2010.hh>

using namespace test;
using namespace eos;

class BToKstarKMPW2010FormFactorsTest:
    public TestCase
{
    public:
        BToKstarKMPW2010FormFactorsTest() :
            TestCase("b_to_kstar_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5e-6;

            // central values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.36; p["B->K^*::b^V_1@KMPW2010"]  =  -4.8;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.29; p["B->K^*::b^A0_1@KMPW2010"] = -18.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.25; p["B->K^*::b^A1_1@KMPW2010"] =  +0.34;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.23; p["B->K^*::b^A2_1@KMPW2010"] =  -0.85;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.31; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.31; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.22; p["B->K^*::b^T3_1@KMPW2010"] = -10.3;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(ff->v( 2.3),   0.421964, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 4.6),   0.497285, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 6.9),   0.590298, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 9.2),   0.707366, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(11.5),   0.858208, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(13.8),   1.058404, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(16.1),   1.334545, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(18.4),   1.735915, eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 2.3), 0.410875, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 4.6), 0.560285, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 6.9), 0.747894, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 9.2), 0.988041, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(11.5), 1.302880, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(13.8), 1.728452, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(16.1), 2.327350, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(18.4), 3.218300, eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 2.3), 0.266631, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 4.6), 0.285779, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 6.9), 0.308075, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 9.2), 0.334383, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(11.5), 0.365917, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(13.8), 0.404440, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(16.1), 0.452618, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(18.4), 0.514681, eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 2.3), 0.250198, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 4.6), 0.273759, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 6.9), 0.301563, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 9.2), 0.334822, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(11.5), 0.375252, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(13.8), 0.425358, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(16.1), 0.488950, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(18.4), 0.572108, eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 2.3), 0.362235, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 4.6), 0.425698, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 6.9), 0.504029, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 9.2), 0.602573, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(11.5), 0.729487, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(13.8), 0.897853, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(16.1), 1.12999 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(18.4), 1.46727 , eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 2.3), 0.35026 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 4.6), 0.397838, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 6.9), 0.45472 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 9.2), 0.523652, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(11.5), 0.608536, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(13.8), 0.715103, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(16.1), 0.852101, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(18.4), 1.03355 , eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 2.3), 0.276523, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 4.6), 0.344214, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 6.9), 0.426193, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 9.2), 0.526794, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(11.5), 0.6522  , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(13.8), 0.811519, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(16.1), 1.01871 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(18.4), 1.29621 , eps);
            }

            // raised values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.59; p["B->K^*::b^V_1@KMPW2010"]  =  -4.0;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.39; p["B->K^*::b^A0_1@KMPW2010"] = -16.9;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.41; p["B->K^*::b^A1_1@KMPW2010"] =  +1.2;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.42; p["B->K^*::b^A2_1@KMPW2010"] =  +2.03;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.49; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6 + 0.81;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.49; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2 + 2.1;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.39; p["B->K^*::b^T3_1@KMPW2010"] = -10.3 + 2.5;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(ff->v( 2.3),   0.683006, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 4.6),   0.795816, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 6.9),   0.934831, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 9.2),   1.109443, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(11.5),   1.333989, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(13.8),   1.631445, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(16.1),   2.041006, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(18.4),   2.635294, eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 2.3), 0.543363, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 4.6), 0.732820, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 6.9), 0.970583, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 9.2), 1.274770, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(11.5), 1.673378, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(13.8), 2.211932, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(16.1), 2.969510, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(18.4), 4.096084, eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 2.3), 0.430965, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 4.6), 0.454709, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 6.9), 0.481882, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 9.2), 0.513361, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(11.5), 0.550366, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(13.8), 0.594652, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(16.1), 0.648841, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(18.4), 0.717050, eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 2.3), 0.435238, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 4.6), 0.451990, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 6.9), 0.470539, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 9.2), 0.491252, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(11.5), 0.514616, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(13.8), 0.541305, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(16.1), 0.572274, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(18.4), 0.608954, eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 2.3), 0.56538 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 4.6), 0.656751, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 6.9), 0.769278, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 9.2), 0.910534, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(11.5), 1.09208 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(13.8), 1.33244 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(16.1), 1.66321 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(18.4), 2.14292 , eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 2.3), 0.535223, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 4.6), 0.588078, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 6.9), 0.650576, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 9.2), 0.725486, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(11.5), 0.816729, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(13.8), 0.930041, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(16.1), 1.07414 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(18.4), 1.26297 , eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 2.3), 0.472753, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 4.6), 0.571573, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 6.9), 0.690926, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 9.2), 0.837005, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(11.5), 1.01864 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(13.8), 1.24883 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(16.1), 1.54748 , eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(18.4), 1.94657 , eps);            }

            // lowered values
            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::F^V(0)@KMPW2010"]  = 0.24; p["B->K^*::b^V_1@KMPW2010"]  =  -5.2;
                p["B->K^*::F^A0(0)@KMPW2010"] = 0.22; p["B->K^*::b^A0_1@KMPW2010"] = -21.2;
                p["B->K^*::F^A1(0)@KMPW2010"] = 0.15; p["B->K^*::b^A1_1@KMPW2010"] =  -0.46;
                p["B->K^*::F^A2(0)@KMPW2010"] = 0.13; p["B->K^*::b^A2_1@KMPW2010"] =  -2.2;
                p["B->K^*::F^T1(0)@KMPW2010"] = 0.21; p["B->K^*::b^T1_1@KMPW2010"] =  -4.6 - 0.41;
                p["B->K^*::F^T2(0)@KMPW2010"] = 0.21; p["B->K^*::b^T2_1@KMPW2010"] =  -3.2 - 2.2;
                p["B->K^*::F^T3(0)@KMPW2010"] = 0.12; p["B->K^*::b^T3_1@KMPW2010"] = -10.3 - 3.1;
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::KMPW2010", p, Options{ });
                TEST_CHECK(0 != ff.get());

                // V
                TEST_CHECK_RELATIVE_ERROR(ff->v( 2.3),    0.2830474, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 4.6),    0.3354240, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 6.9),    0.4001625, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 9.2),    0.4817170, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(11.5),    0.5868880, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(13.8),    0.7265850, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(16.1),    0.9194250, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(18.4),    1.1999232, eps);

                // A_0
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 2.3),  0.3236650, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 4.6),  0.4519470, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 6.9),  0.6131980, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 9.2),  0.8198090, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(11.5),  1.0909360, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(13.8),  1.4577380, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(16.1),  1.9743500, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(18.4),  2.7434500, eps);

                // A_1
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 2.3),  0.1621258, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 4.6),  0.1762210, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 6.9),  0.1927954, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 9.2),  0.2125507, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(11.5),  0.2364774, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(13.8),  0.2660210, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(16.1),  0.3033756, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(18.4),  0.3520390, eps);

                // A_2
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 2.3),  0.1445568, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 4.6),  0.1616854, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 6.9),  0.1820760, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 9.2),  0.2066815, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(11.5),  0.2368550, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(13.8),  0.2745790, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(16.1),  0.3228774, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(18.4),  0.3865910, eps);

                // t_1
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 2.3),  0.246944,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 4.6),  0.291875,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 6.9),  0.347386,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 9.2),  0.417288,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(11.5),  0.507397,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(13.8),  0.627041,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(16.1),  0.792141,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(18.4),  1.03221 ,  eps);

                // t_2
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 2.3),  0.24554 ,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 4.6),  0.287805,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 6.9),  0.338645,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 9.2),  0.400627,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(11.5),  0.477404,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(13.8),  0.574348,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(16.1),  0.699678,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(18.4),  0.866585,  eps);

                // t_3
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 2.3),  0.157488,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 4.6),  0.202489,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 6.9),  0.257115,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 9.2),  0.324297,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(11.5),  0.40822 ,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(13.8),  0.515053,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(16.1),  0.654254,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(18.4),  0.841037,  eps);
            }
        }
} b_to_kstar_kmpw2010_form_factors_test;

class BToKKMPW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKKMPW2010FormFactorsTest() :
            TestCase("b_to_k_kmpw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-10;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K::KMPW2010", p, Options{ });
            TEST_CHECK(0 != ff.get());

            TEST_CHECK_NEARLY_EQUAL(ff->f_p( 2.3), 0.3844144474375, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p( 4.6), 0.4381960494587, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p( 6.9), 0.5043916755865, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p( 9.2), 0.5874824689246, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p(11.5), 0.6943252577904, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p(13.8), 0.8359400978914, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p(16.1), 1.0311821367700, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p(18.4), 1.3150955485260, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.7), 1.7607545970800, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_t( 2.3), 0.4417893713152, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t( 4.6), 0.5045402114062, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t( 6.9), 0.5818230155226, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t( 9.2), 0.6788902047561, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t(11.5), 0.8037805019019, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t(13.8), 0.9694161202369, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t(16.1), 1.1979110863960, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t(18.4), 1.5303740143290, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.7), 2.0525389704830, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_0( 2.3), 0.3691335300041, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0( 4.6), 0.4001522415666, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0( 6.9), 0.4333067770539, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0( 9.2), 0.4689025169045, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0(11.5), 0.5073172894003, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0(13.8), 0.5490271105427, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0(16.1), 0.5946449022903, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0(18.4), 0.6449812193912, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.7), 0.7011445924499, eps);
        }
} b_to_k_kmpw2010_form_factors_test;
