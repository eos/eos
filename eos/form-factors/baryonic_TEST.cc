/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/form-factors/form-factors.hh>

#include <cmath>
#include <limits>

using namespace test;
using namespace eos;

class BFvD2014FormFactorsTest :
    public TestCase
{
    public:
        BFvD2014FormFactorsTest() :
            TestCase("bfvd2014_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-3;

            Parameters p = Parameters::Defaults();
            std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::BFvD2014", p, Options{ });

            p["Lambda_b->Lambda::f_0^V(0)@BFvD2014"]    =  0.33;
            p["Lambda_b->Lambda::b_1_0^V@BFvD2014"]     = -1.75;
            p["Lambda_b->Lambda::f_0^A(0)@BFvD2014"]    =  0.31;
            p["Lambda_b->Lambda::b_1_0^A@BFvD2014"]     = -0.52;
            p["Lambda_b->Lambda::f_perp^V(0)@BFvD2014"] =  0.34;
            p["Lambda_b->Lambda::b_1_perp^V@BFvD2014"]  = -1.58;
            p["Lambda_b->Lambda::f_perp^A(0)@BFvD2014"] =  0.31;
            p["Lambda_b->Lambda::b_1_perp^A@BFvD2014"]  = -0.24;
            p["mass::Lambda_b"]                         = 5.6194;
            p["mass::Lambda"]                           = 1.1157;

            TEST_CHECK_NEARLY_EQUAL(0.330, ff->f_long_v( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.418, ff->f_long_v( 5.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.555, ff->f_long_v(10.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.794, ff->f_long_v(15.0), eps);
            TEST_CHECK_NEARLY_EQUAL(1.302, ff->f_long_v(20.0), eps);

            TEST_CHECK_NEARLY_EQUAL(0.310, ff->f_long_a( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.369, ff->f_long_a( 5.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.453, ff->f_long_a(10.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.584, ff->f_long_a(15.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.810, ff->f_long_a(20.0), eps);

            TEST_CHECK_NEARLY_EQUAL(0.340, ff->f_perp_v( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.429, ff->f_perp_v( 5.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.567, ff->f_perp_v(10.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.806, ff->f_perp_v(15.0), eps);
            TEST_CHECK_NEARLY_EQUAL(1.315, ff->f_perp_v(20.0), eps);

            TEST_CHECK_NEARLY_EQUAL(0.310, ff->f_perp_a( 0.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.366, ff->f_perp_a( 5.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.446, ff->f_perp_a(10.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.568, ff->f_perp_a(15.0), eps);
            TEST_CHECK_NEARLY_EQUAL(0.780, ff->f_perp_a(20.0), eps);
        }
} bfvd2014_form_factors_test;

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
        TEST_CHECK_NEARLY_EQUAL(0.1562487236, ff->f_time_v( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2275900892, ff->f_time_v( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.3417796225, ff->f_time_v(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.5422197258, ff->f_time_v(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1598530053, ff->f_long_v( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2459877325, ff->f_long_v( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.3910960039, ff->f_long_v(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.6661662119, ff->f_long_v(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2073776573, ff->f_perp_v( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.3131482081, ff->f_perp_v( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.4907199746, ff->f_perp_v(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.8262231044, ff->f_perp_v(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1657965242, ff->f_time_a( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2489615034, ff->f_time_a( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.3895093452, ff->f_time_a(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.6583338379, ff->f_time_a(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1118800889, ff->f_long_a( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1803545293, ff->f_long_a( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2912190226, ff->f_long_a(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.4874045714, ff->f_long_a(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.0947209451, ff->f_perp_a( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1635457816, ff->f_perp_a( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2758041356, ff->f_perp_a(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.4758361015, ff->f_perp_a(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2363096025, ff->f_long_t( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.3376352752, ff->f_long_t( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.5056514080, ff->f_long_t(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.8193320135, ff->f_long_t(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1659591402, ff->f_perp_t( 0.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.2450975156, ff->f_perp_t( 5.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.3773577270, ff->f_perp_t(10.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.6261686996, ff->f_perp_t(15.0),  eps);
        TEST_CHECK_NEARLY_EQUAL(0.1630195575, ff->f_long_t5( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.2272722059, ff->f_long_t5( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.3285959782, ff->f_long_t5(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.5033819243, ff->f_long_t5(15.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.1558564787, ff->f_perp_t5( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.2202553997, ff->f_perp_t5( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.3221610388, ff->f_perp_t5(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.4985526705, ff->f_perp_t5(15.0), eps);
    }
} dm2016_form_factors_test;

class DKMR2017FormFactorsTest :
public TestCase
{
public:
    DKMR2017FormFactorsTest() :
    TestCase("dkmr2017_form_factors_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 1e-3;

        Parameters p = Parameters::Defaults();
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda_c::DKMR2017", p);
        // fix z-expansion parameters @DKMR2017 - all a_2's are zero by default
        // Note that a0_perp_A = a0_long_A ve a0_perp_T5 = a0_long_T5
        p["Lambda_b->Lambda_c::a_0_time^V@DKMR2017"]    =  0.33;
        p["Lambda_b->Lambda_c::a_1_time^V@DKMR2017"]    = -0.43;
        p["Lambda_b->Lambda_c::a_2_time^V@DKMR2017"]    = +0.00;

        p["Lambda_b->Lambda_c::a_0_long^V@DKMR2017"]    =  0.44;
        p["Lambda_b->Lambda_c::a_1_long^V@DKMR2017"]    =  0.27;
        p["Lambda_b->Lambda_c::a_2_long^V@DKMR2017"]    =  0.00;

        p["Lambda_b->Lambda_c::a_0_perp^V@DKMR2017"]    = -0.52;
        p["Lambda_b->Lambda_c::a_1_perp^V@DKMR2017"]    = +0.42;
        p["Lambda_b->Lambda_c::a_2_perp^V@DKMR2017"]    = +0.00;

        p["Lambda_b->Lambda_c::a_0_time^A@DKMR2017"]    = +1.75;
        p["Lambda_b->Lambda_c::a_1_time^A@DKMR2017"]    = -1.1;
        p["Lambda_b->Lambda_c::a_2_time^A@DKMR2017"]    = -0.00;

        p["Lambda_b->Lambda_c::a_0_long^A@DKMR2017"]    =  2.14;
        p["Lambda_b->Lambda_c::a_1_long^A@DKMR2017"]    = -1.58;
        p["Lambda_b->Lambda_c::a_2_long^A@DKMR2017"]    = -0.00;
        //
        p["Lambda_b->Lambda_c::a_1_perp^A@DKMR2017"]    = -1.30;
        p["Lambda_b->Lambda_c::a_2_perp^A@DKMR2017"]    = -0.00;

        p["Lambda_b->Lambda_c::a_0_long^T@DKMR2017"]    =  0.71;
        p["Lambda_b->Lambda_c::a_1_long^T@DKMR2017"]    =  0.21;
        p["Lambda_b->Lambda_c::a_2_long^T@DKMR2017"]    =  0.00;

        p["Lambda_b->Lambda_c::a_0_perp^T@DKMR2017"]    = -0.24;
        p["Lambda_b->Lambda_c::a_1_perp^T@DKMR2017"]    = -0.24;
        p["Lambda_b->Lambda_c::a_2_perp^T@DKMR2017"]    = -0.00;

        p["Lambda_b->Lambda_c::a_0_long^T5@DKMR2017"]   =  0.14;
        p["Lambda_b->Lambda_c::a_1_long^T5@DKMR2017"]   = -0.14;
        p["Lambda_b->Lambda_c::a_2_long^T5@DKMR2017"]   = -0.00;
        //
        p["Lambda_b->Lambda_c::a_1_perp^T5@DKMR2017"]   =  3.14;
        p["Lambda_b->Lambda_c::a_2_perp^T5@DKMR2017"]   =  0.00;

        p["mass::Lambda_b"]                             = 5.61951;
        p["mass::Lambda_c"]                             = 2.2865;
        // test of the ten baryonic FFs
        TEST_CHECK_NEARLY_EQUAL(0.299748, ff->f_time_v( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.351122, ff->f_time_v( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.419267, ff->f_time_v(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.513241, ff->f_time_v(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(0.461852, ff->f_long_v( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.517426, ff->f_long_v( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.589584, ff->f_long_v(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.687469, ff->f_long_v(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(-0.486008, ff->f_perp_v( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(-0.571162, ff->f_perp_v( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(-0.687539, ff->f_perp_v(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(-0.855001, ff->f_perp_v(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(1.65909, ff->f_time_a( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(1.94289, ff->f_time_a( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(2.3313,  ff->f_time_a(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(2.89206, ff->f_time_a(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(2.03046, ff->f_long_a( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(2.33035, ff->f_long_a( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(2.72177, ff->f_long_a(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(3.25185, ff->f_long_a(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(2.04987, ff->f_perp_a( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(2.34308, ff->f_perp_a( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(2.72459, ff->f_perp_a(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(3.23947, ff->f_perp_a(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(0.726996, ff->f_long_t( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.822619, ff->f_long_t( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.948552, ff->f_long_t(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(1.122315, ff->f_long_t(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(-0.259424, ff->f_perp_t( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(-0.287293, ff->f_perp_t( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(-0.322751, ff->f_perp_t(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(-0.369644, ff->f_perp_t(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(0.130294, ff->f_long_t5( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.150786, ff->f_long_t5( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.177691, ff->f_long_t5(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.214357, ff->f_long_t5(15.0), eps);

        TEST_CHECK_NEARLY_EQUAL(0.357693, ff->f_perp_t5( 0.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.299983, ff->f_perp_t5( 5.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.210694, ff->f_perp_t5(10.0), eps);
        TEST_CHECK_NEARLY_EQUAL(0.069375, ff->f_perp_t5(15.0), eps);
    }
} dkmr2017_form_factors_test;

class HQETOneHalfFormFactorsTest :
    public TestCase
{
    public:
        HQETOneHalfFormFactorsTest() :
            TestCase("hqet_one_half_form_factors_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["Lambda_b->Lambda_c^*::zeta(q^2_max)@HQET"] =  1.00;
            p["Lambda_b->Lambda_c^*::delta_3b@HQET"]      = -0.14;
            p["Lambda_b->Lambda_c^*::rho@HQET"]           =  0.25;
            p["Lambda_b->Lambda_c^*::rho_3b@HQET"]        =  0.25;

            auto ff = FormFactorFactory<OneHalfPlusToOneHalfMinus>::create("Lambda_b->Lambda_c(2595)::HQET", p);

            //Diagnostics diag = ff->diagnostics();
            //for (auto d : diag)
            //{
            //    std::cout << d.description << ": " << d.value << std::endl;
            //}

            static constexpr double eps   = 1.0e-5;
            static constexpr double s_max = 9.1643031076;

            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_time_v( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0296321, ff->f_long_v( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.1021055, ff->f_perp_v( s_max), eps);

            TEST_CHECK_NEARLY_EQUAL(-0.0071049, ff->f_time_a( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_long_a( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_perp_a( s_max), eps);

            TEST_CHECK_NEARLY_EQUAL( 0.9739324, ff->f_time_v( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.2487909, ff->f_long_v( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.1707046, ff->f_perp_v( s_max - 3.0), eps);

            TEST_CHECK_NEARLY_EQUAL( 0.1809885, ff->f_time_a( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.8387119, ff->f_long_a( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.8467812, ff->f_perp_a( s_max - 3.0), eps);
        }
} hqet_one_half_form_factors_test;

class HQETThreeHalfFormFactorsTest :
    public TestCase
{
    public:
        HQETThreeHalfFormFactorsTest() :
            TestCase("hqet_three_half_form_factors_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["Lambda_b->Lambda_c^*::zeta(q^2_max)@HQET"] =  1.00;
            p["Lambda_b->Lambda_c^*::delta_3b@HQET"]      = -0.14;
            p["Lambda_b->Lambda_c^*::rho@HQET"]           =  0.25;
            p["Lambda_b->Lambda_c^*::rho_3b@HQET"]        =  0.25;

            auto ff = FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda_c(2625)::HQET", p);

            //Diagnostics diag = ff->diagnostics();
            //for (auto d : diag)
            //{
            //    std::cout << d.description << ": " << d.value << std::endl;
            //}

            static constexpr double eps   = 1.0e-5;
            static constexpr double s_max = 8.9484739600;

            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_time12_v( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.0583747, ff->f_long12_v( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.2814541, ff->f_perp12_v( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0249132, ff->f_perp32_v( s_max), eps);

            TEST_CHECK_NEARLY_EQUAL(-0.0403027, ff->f_time12_a( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_long12_a( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_perp12_a( s_max), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0,       ff->f_perp32_a( s_max), eps);

            TEST_CHECK_NEARLY_EQUAL( 0.7297828, ff->f_time12_v( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.1034446, ff->f_long12_v( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL(-0.1125738, ff->f_perp12_v( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0408266, ff->f_perp32_v( s_max - 3.0), eps);

            TEST_CHECK_NEARLY_EQUAL( 0.1251307, ff->f_time12_a( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.7419465, ff->f_long12_a( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.7770269, ff->f_perp12_a( s_max - 3.0), eps);
            TEST_CHECK_NEARLY_EQUAL( 0.0089753, ff->f_perp32_a( s_max - 3.0), eps);
        }
} hqet_three_half_form_factors_test;
