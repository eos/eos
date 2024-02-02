/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Frederik Beaujean
 * Copyright (c) 2018 Ahmet Kokulu
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
#include <eos/form-factors/parametric-bsz2015-impl.hh>

using namespace test;
using namespace eos;

class BToPiBSZ2015FormFactorsTest :
    public TestCase
{
    public:
        BToPiBSZ2015FormFactorsTest() :
            TestCase("b_to_pi_bsz2015_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                // test case created by using the known relations among the BCL2008 parameters
                // for the highest power.
                static const double eps = 1e-5;

                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::BSZ2015", p, Options{ });

                p["B->pi::alpha^f+_0@BSZ2015"] = 1.0;
                p["B->pi::alpha^f+_1@BSZ2015"] = 0.0;
                p["B->pi::alpha^f+_2@BSZ2015"] = 0.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.21408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.54479, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 2.12312, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 3.39360, eps);

                p["B->pi::alpha^f+_1@BSZ2015"] = 1.0;
                p["B->pi::alpha^f+_2@BSZ2015"] = 2.0;

                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p( 5.0), 1.16581, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(10.0), 1.42261, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(15.0), 1.88375, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_p(20.0), 2.97499, eps);

                p["B->pi::alpha^f0_1@BSZ2015"] = 1.5;
                p["B->pi::alpha^f0_2@BSZ2015"] = 1.5;

                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0( 5.0), 1.11998, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(10.0), 1.28572, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(15.0), 1.53862, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_0(20.0), 2.00499, eps);

                p["B->pi::alpha^fT_0@BSZ2015"] =  1.0;
                p["B->pi::alpha^fT_1@BSZ2015"] = -1.0;
                p["B->pi::alpha^fT_2@BSZ2015"] =  2.5;

                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 0.0), 1.00000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t( 5.0), 1.27271, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(10.0), 1.73442, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(15.0), 2.64425, eps);
                TEST_CHECK_NEARLY_EQUAL(ff->f_t(20.0), 4.99850, eps);
            }
        }
} b_to_pi_bsz2015_form_factors_test;

class BToDstarBSZ2015FormFactorsTest :
    public TestCase
{
    public:
        BToDstarBSZ2015FormFactorsTest() :
            TestCase("b_to_dstar_bsz2015_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5.1e-3;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->D^*::BSZ2015", p, Options{ });
            TEST_CHECK(ff.get() != nullptr);

            /*
            * for the B->D^* case no SSE parameters are available yet, so for the moment we set them to zero and once the GKvD2018 calculation is fitted we can then update the test case later on.
            */
            p["B->D^*::alpha^A0_0@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^A0_1@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^A0_2@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^A1_0@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^A1_1@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^A1_2@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^A12_1@BSZ2015"] = 0.0;
            p["B->D^*::alpha^A12_2@BSZ2015"] = 0.0;
            p["B->D^*::alpha^V_0@BSZ2015"  ] = 0.0;
            p["B->D^*::alpha^V_1@BSZ2015"  ] = 0.0;
            p["B->D^*::alpha^V_2@BSZ2015"  ] = 0.0;
            p["B->D^*::alpha^T1_0@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^T1_1@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^T1_2@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^T2_1@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^T2_2@BSZ2015" ] = 0.0;
            p["B->D^*::alpha^T23_0@BSZ2015"] = 0.0;
            p["B->D^*::alpha^T23_1@BSZ2015"] = 0.0;
            p["B->D^*::alpha^T23_2@BSZ2015"] = 0.0;

            TEST_CHECK_NEARLY_EQUAL(ff->a_0(0.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_0(2.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_0(4.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_0(6.1), 0.0, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->a_1(0.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_1(2.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_1(4.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_1(6.1), 0.0, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->a_2(0.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_2(2.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_2(4.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_2(6.1), 0.0, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->v(0.1),   0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->v(2.1),   0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->v(4.1),   0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->v(6.1),   0.0, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->t_1(0.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_1(2.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_1(4.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_1(6.1), 0.0, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->t_2(0.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_2(2.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_2(4.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_2(6.1), 0.0, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->t_3(0.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_3(2.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_3(4.1), 0.0, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_3(6.1), 0.0, eps);
        }
} b_to_dstar_bsz2015_form_factors_test;

class BToKstarBSZ2015FormFactorsTest :
    public TestCase
{
    public:
        BToKstarBSZ2015FormFactorsTest() :
            TestCase("b_to_kstar_bsz2015_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 5.1e-3;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::BSZ2015", p, Options{ });
            TEST_CHECK(ff.get() != nullptr);

            /* compare with values from David Straub,
             * use his values of the parameters
             */
            p["B->K^*::alpha^A0_0@BSZ2015" ] = +0.39;
            p["B->K^*::alpha^A0_1@BSZ2015" ] = -1.15;
            p["B->K^*::alpha^A0_2@BSZ2015" ] = +2.08;
            p["B->K^*::alpha^A1_0@BSZ2015" ] = +0.29;
            p["B->K^*::alpha^A1_1@BSZ2015" ] = +0.31;
            p["B->K^*::alpha^A1_2@BSZ2015" ] = +0.72;
            p["B->K^*::alpha^A12_1@BSZ2015"] = +0.57;
            p["B->K^*::alpha^A12_2@BSZ2015"] = +0.14;
            p["B->K^*::alpha^V_0@BSZ2015"  ] = +0.37;
            p["B->K^*::alpha^V_1@BSZ2015"  ] = -1.08;
            p["B->K^*::alpha^V_2@BSZ2015"  ] = +2.47;
            p["B->K^*::alpha^T1_0@BSZ2015" ] = +0.31;
            p["B->K^*::alpha^T1_1@BSZ2015" ] = -0.96;
            p["B->K^*::alpha^T1_2@BSZ2015" ] = +2.01;
            p["B->K^*::alpha^T2_1@BSZ2015" ] = +0.42;
            p["B->K^*::alpha^T2_2@BSZ2015" ] = +2.02;
            p["B->K^*::alpha^T23_0@BSZ2015"] = +0.79;
            p["B->K^*::alpha^T23_1@BSZ2015"] = +1.26;
            p["B->K^*::alpha^T23_2@BSZ2015"] = +1.96;

            TEST_CHECK_NEARLY_EQUAL(ff->a_0(0.1), 0.393136, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_0(2.1), 0.440394, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_0(4.1), 0.496878, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_0(6.1), 0.565342, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->a_1(0.1), 0.289606, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_1(2.1), 0.3039,   eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_1(4.1), 0.319847, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_1(6.1), 0.337861, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->a_2(0.1), 0.248569, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_2(2.1), 0.272718, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_2(4.1), 0.300677, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->a_2(6.1), 0.333431, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->v(0.1),   0.367312, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->v(2.1),   0.411249, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->v(4.1),   0.463812, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->v(6.1),   0.527595, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->t_1(0.1), 0.3094,   eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_1(2.1), 0.346962, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_1(4.1), 0.391946, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_1(6.1), 0.446575, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->t_2(0.1), 0.308387, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_2(2.1), 0.322844, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_2(4.1), 0.339239, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_2(6.1), 0.358183, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->t_3(0.1), 0.184952, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_3(2.1), 0.200925, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_3(4.1), 0.219004, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->t_3(6.1), 0.239587, eps);
        }
} b_to_kstar_bsz2015_form_factors_test;

class BToRhoBSZ2015FormFactorsTest :
public TestCase
{
public:
    BToRhoBSZ2015FormFactorsTest() :
    TestCase("b_to_rho_bsz2015_form_factors_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 5.1e-3;

        Parameters p = Parameters::Defaults();
        std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->rho::BSZ2015", p, Options{ });
        TEST_CHECK(ff.get() != nullptr);

        /*
         * use David Straub's values of the SSE parameters from LCSR only
         */
        p["B->rho::alpha^A0_0@BSZ2015" ] = +0.36;
        p["B->rho::alpha^A0_1@BSZ2015" ] = -0.83;
        p["B->rho::alpha^A0_2@BSZ2015" ] = +1.33;
        p["B->rho::alpha^A1_0@BSZ2015" ] = +0.26;
        p["B->rho::alpha^A1_1@BSZ2015" ] = +0.39;
        p["B->rho::alpha^A1_2@BSZ2015" ] = +0.16;
        p["B->rho::alpha^A12_1@BSZ2015"] = +0.76;
        p["B->rho::alpha^A12_2@BSZ2015"] = +0.46;
        p["B->rho::alpha^V_0@BSZ2015"  ] = +0.33;
        p["B->rho::alpha^V_1@BSZ2015"  ] = -0.86;
        p["B->rho::alpha^V_2@BSZ2015"  ] = +1.80;
        p["B->rho::alpha^T1_0@BSZ2015" ] = +0.27;
        p["B->rho::alpha^T1_1@BSZ2015" ] = -0.74;
        p["B->rho::alpha^T1_2@BSZ2015" ] = +1.45;
        p["B->rho::alpha^T2_1@BSZ2015" ] = +0.47;
        p["B->rho::alpha^T2_2@BSZ2015" ] = +0.58;
        p["B->rho::alpha^T23_0@BSZ2015"] = +0.75;
        p["B->rho::alpha^T23_1@BSZ2015"] = +1.90;
        p["B->rho::alpha^T23_2@BSZ2015"] = +2.93;

        TEST_CHECK_NEARLY_EQUAL(ff->a_0(0.1), 0.36186,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_0(2.1), 0.402772, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_0(4.1), 0.4521,   eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_0(6.1), 0.512422, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->a_1(0.1), 0.260532, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_1(2.1), 0.271749, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_1(4.1), 0.284225, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_1(6.1), 0.29821,  eps);

        TEST_CHECK_NEARLY_EQUAL(ff->a_2(0.1), 0.226525, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_2(2.1), 0.247972, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_2(4.1), 0.272771, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->a_2(6.1), 0.301645, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->v(0.1),   0.331752, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->v(2.1),   0.370391, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->v(4.1),   0.417199, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->v(6.1),   0.474694, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->t_1(0.1), 0.271458, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_1(2.1), 0.303616, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_1(4.1), 0.342573, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_1(6.1), 0.390422, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->t_2(0.1), 0.270508, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_2(2.1), 0.28128,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_2(4.1), 0.29338,  eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_2(6.1), 0.307107, eps);

        TEST_CHECK_NEARLY_EQUAL(ff->t_3(0.1), 0.179298, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_3(2.1), 0.196067, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_3(4.1), 0.215543, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->t_3(6.1), 0.238325, eps);
    }
} b_to_rho_bsz2015_form_factors_test;
