/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
            std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda@BFvD2014", p);

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

            auto ff = FormFactorFactory<OneHalfPlusToOneHalfMinus>::create("Lambda_b->Lambda_c(2595)@HQET", p);

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

            auto ff = FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda_c(2625)@HQET", p);

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
