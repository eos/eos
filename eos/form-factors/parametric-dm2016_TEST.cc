/*
 * Copyright (c) 2014, 2015, 2016, 2017 Danny van Dyk
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
#include <eos/form-factors/parametric-dm2016.hh>
#include <eos/form-factors/parametric-dm2016-impl.hh>

using namespace test;
using namespace eos;

class DM2016FormFactorsTest :
    public TestCase
{
public:
    DM2016FormFactorsTest() :
        TestCase("dm2016_form_factors_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 1e-3;

        Parameters p = Parameters::Defaults();
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::DM2016", p);
        // fix z-expansion parameters @DM2016 - all a_2's are zero by default
        p["Lambda_b->Lambda::a_0_time^V@DM2016"]    =  0.3725;
        p["Lambda_b->Lambda::a_1_time^V@DM2016"]    = -0.9389;
        p["Lambda_b->Lambda::a_2_time^V@DM2016"]    = +0.0000;

        p["Lambda_b->Lambda::a_0_long^V@DM2016"]    =  0.4221;
        p["Lambda_b->Lambda::a_1_long^V@DM2016"]    = -1.1386;
        p["Lambda_b->Lambda::a_2_long^V@DM2016"]    =  0.0000;

        p["Lambda_b->Lambda::a_0_perp^V@DM2016"]    = +0.5182;
        p["Lambda_b->Lambda::a_1_perp^V@DM2016"]    = -1.3495;
        p["Lambda_b->Lambda::a_2_perp^V@DM2016"]    = +0.0000;

        p["Lambda_b->Lambda::a_0_time^A@DM2016"]    = +0.4028;
        p["Lambda_b->Lambda::a_1_time^A@DM2016"]    = -1.0290;
        p["Lambda_b->Lambda::a_2_time^A@DM2016"]    = -0.0000;

        p["Lambda_b->Lambda::a_0_long^A@DM2016"]    =  0.3563;
        p["Lambda_b->Lambda::a_1_long^A@DM2016"]    = -1.0612;
        p["Lambda_b->Lambda::a_2_long^A@DM2016"]    = -0.0000;

        // a_0_perp^A == a_0_long^A
        p["Lambda_b->Lambda::a_1_perp^A@DM2016"]    = -1.1357;
        p["Lambda_b->Lambda::a_2_perp^A@DM2016"]    = -0.0000;

        p["Lambda_b->Lambda::a_0_long^T@DM2016"]    =  0.4960;
        p["Lambda_b->Lambda::a_1_long^T@DM2016"]    = -1.1275;
        p["Lambda_b->Lambda::a_2_long^T@DM2016"]    =  0.0000;

        p["Lambda_b->Lambda::a_0_perp^T@DM2016"]    =  0.3876;
        p["Lambda_b->Lambda::a_1_perp^T@DM2016"]    = -0.9623;
        p["Lambda_b->Lambda::a_2_perp^T@DM2016"]    = -0.0000;

        p["Lambda_b->Lambda::a_0_long^T5@DM2016"]   =  0.3403;
        p["Lambda_b->Lambda::a_1_long^T5@DM2016"]   = -0.7697;
        p["Lambda_b->Lambda::a_2_long^T5@DM2016"]   = -0.0000;

        // a_0_perp^T5 == a_0_long^T5
        p["Lambda_b->Lambda::a_1_perp^T5@DM2016"]   = -0.8008;
        p["Lambda_b->Lambda::a_2_perp^T5@DM2016"]   =  0.0000;

        p["mass::Lambda_b"]                             = 5.61951;
        p["mass::Lambda"]                               = 1.1157;

        // test of the ten baryonic FFs
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v( 0.0),  0.1562487236, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v( 5.0),  0.2275900892, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v(10.0),  0.3417796225, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_v(15.0),  0.5422197258, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v( 0.0),  0.1598530053, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v( 5.0),  0.2459877325, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v(10.0),  0.3910960039, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_v(15.0),  0.6661662119, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v( 0.0),  0.2073776573, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v( 5.0),  0.3131482081, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v(10.0),  0.4907199746, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v(15.0),  0.8262231044, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a( 0.0),  0.1657965242, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a( 5.0),  0.2489615034, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a(10.0),  0.3895093452, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_time_a(15.0),  0.6583338379, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a( 0.0),  0.1118800889, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a( 5.0),  0.1803545293, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a(10.0),  0.2912190226, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_a(15.0),  0.4874045714, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a( 0.0),  0.0947209451, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a( 5.0),  0.1635457816, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a(10.0),  0.2758041356, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a(15.0),  0.4758361015, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t( 0.0),  0.2363096025, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t( 5.0),  0.3376352752, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t(10.0),  0.5056514080, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t(15.0),  0.8193320135, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t( 0.0),  0.1659591402, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t( 5.0),  0.2450975156, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t(10.0),  0.3773577270, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t(15.0),  0.6261686996, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5( 0.0), 0.1630195575, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5( 5.0), 0.2272722059, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5(10.0), 0.3285959782, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_long_t5(15.0), 0.5033819243, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5( 0.0), 0.1558564787, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5( 5.0), 0.2202553997, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5(10.0), 0.3221610388, eps);
        TEST_CHECK_NEARLY_EQUAL(ff->f_perp_t5(15.0), 0.4985526705, eps);
    }
} dm2016_form_factors_test;
