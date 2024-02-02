/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
 * Copyright (c) 2017 Elena Graverini
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

            TEST_CHECK_NEARLY_EQUAL(ff->f_time_v( s_max),  0.0,        eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long_v( s_max),  0.0296321,  eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v( s_max), -0.1021055,  eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_time_a( s_max), -0.0071049,  eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long_a( s_max),  0.0,        eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a( s_max),  0.0,        eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_time_v( s_max - 3.0), 0.9739324, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long_v( s_max - 3.0), 0.2487909, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp_v( s_max - 3.0), 0.1707046, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_time_a( s_max - 3.0), 0.1809885, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long_a( s_max - 3.0), 0.8387119, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp_a( s_max - 3.0), 0.8467812, eps);
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

            TEST_CHECK_NEARLY_EQUAL(ff->f_time12_v( s_max),  0.0,       eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long12_v( s_max), -0.0583747, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp12_v( s_max), -0.2814541, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp32_v( s_max),  0.0249132, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_time12_a( s_max), -0.0403027, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long12_a( s_max),  0.0,       eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp12_a( s_max),  0.0,       eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp32_a( s_max),  0.0,       eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_time12_v( s_max - 3.0), 0.7297828, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long12_v( s_max - 3.0), 0.1034446, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp12_v( s_max - 3.0),-0.1125738, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp32_v( s_max - 3.0), 0.0408266, eps);

            TEST_CHECK_NEARLY_EQUAL(ff->f_time12_a( s_max - 3.0), 0.1251307, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_long12_a( s_max - 3.0), 0.7419465, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp12_a( s_max - 3.0), 0.7770269, eps);
            TEST_CHECK_NEARLY_EQUAL(ff->f_perp32_a( s_max - 3.0), 0.0089753, eps);
        }
} hqet_three_half_form_factors_test;
