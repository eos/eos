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

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.21408, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.54479, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.12312, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(3.39360, ff->f_p(20.0), eps);

                p["B->pi::alpha^f+_1@BSZ2015"] = 1.0;
                p["B->pi::alpha^f+_2@BSZ2015"] = 2.0;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_p( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.16581, ff->f_p( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.42261, ff->f_p(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.88375, ff->f_p(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.97499, ff->f_p(20.0), eps);

                p["B->pi::alpha^f0_1@BSZ2015"] = 1.5;
                p["B->pi::alpha^f0_2@BSZ2015"] = 1.5;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_0( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.11998, ff->f_0( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.28572, ff->f_0(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.53862, ff->f_0(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.00499, ff->f_0(20.0), eps);

                p["B->pi::alpha^fT_0@BSZ2015"] =  1.0;
                p["B->pi::alpha^fT_1@BSZ2015"] = -1.0;
                p["B->pi::alpha^fT_2@BSZ2015"] =  2.5;

                TEST_CHECK_NEARLY_EQUAL(1.00000, ff->f_t( 0.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.27271, ff->f_t( 5.0), eps);
                TEST_CHECK_NEARLY_EQUAL(1.73442, ff->f_t(10.0), eps);
                TEST_CHECK_NEARLY_EQUAL(2.64425, ff->f_t(15.0), eps);
                TEST_CHECK_NEARLY_EQUAL(4.99850, ff->f_t(20.0), eps);
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

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_0(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_1(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->a_2(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->v(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_1(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_2(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.0, ff->t_3(6.1), eps);
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

            TEST_CHECK_NEARLY_EQUAL(0.393136, ff->a_0(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.440394, ff->a_0(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.496878, ff->a_0(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.565342, ff->a_0(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.289606, ff->a_1(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.3039,   ff->a_1(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.319847, ff->a_1(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.337861, ff->a_1(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.248569, ff->a_2(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.272718, ff->a_2(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.300677, ff->a_2(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.333431, ff->a_2(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.367312, ff->v(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.411249, ff->v(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.463812, ff->v(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.527595, ff->v(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.3094,   ff->t_1(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.346962, ff->t_1(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.391946, ff->t_1(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.446575, ff->t_1(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.308387, ff->t_2(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.322844, ff->t_2(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.339239, ff->t_2(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.358183, ff->t_2(6.1), eps);

            TEST_CHECK_NEARLY_EQUAL(0.184952, ff->t_3(0.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.200925, ff->t_3(2.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.219004, ff->t_3(4.1), eps);
            TEST_CHECK_NEARLY_EQUAL(0.239587, ff->t_3(6.1), eps);
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

        TEST_CHECK_NEARLY_EQUAL(0.36186,  ff->a_0(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.402772, ff->a_0(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.4521,   ff->a_0(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.512422, ff->a_0(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.260532, ff->a_1(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.271749, ff->a_1(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.284225, ff->a_1(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.29821,  ff->a_1(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.226525, ff->a_2(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.247972, ff->a_2(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.272771, ff->a_2(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.301645, ff->a_2(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.331752, ff->v(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.370391, ff->v(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.417199, ff->v(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.474694, ff->v(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.271458, ff->t_1(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.303616, ff->t_1(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.342573, ff->t_1(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.390422, ff->t_1(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.270508, ff->t_2(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.28128,  ff->t_2(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.29338,  ff->t_2(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.307107, ff->t_2(6.1), eps);

        TEST_CHECK_NEARLY_EQUAL(0.179298, ff->t_3(0.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.196067, ff->t_3(2.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.215543, ff->t_3(4.1), eps);
        TEST_CHECK_NEARLY_EQUAL(0.238325, ff->t_3(6.1), eps);
    }
} b_to_rho_bsz2015_form_factors_test;
