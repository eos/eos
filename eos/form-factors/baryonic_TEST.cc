/*
 * Copyright (c) 2023 Danny van Dyk
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
#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/form-factors.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

using namespace test;
using namespace eos;

class OneHalfPlusToOneHalfPlusFormFactorsTest :
    public TestCase
{
    public:
        OneHalfPlusToOneHalfPlusFormFactorsTest() :
            TestCase("one_half_plus_to_one_half_plus_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Foo->Bar::DM2015",         parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::FooBar", parameter, options));
            }
        }
} one_half_plus_to_one_half_plus_form_factor_test;

class OneHalfPlusToOneHalfMinusFormFactorsTest :
    public TestCase
{
    public:
        OneHalfPlusToOneHalfMinusFormFactorsTest() :
            TestCase("one_half_plus_to_one_half_minus_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<OneHalfPlusToOneHalfMinus>::create("Foo->Bar::HQET",             parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<OneHalfPlusToOneHalfMinus>::create("Lambda_b->Lambda_c::FooBar", parameter, options));
            }
        }
} one_half_plus_to_one_half_minus_form_factor_test;

class OneHalfPlusToThreeHalfMinusFormFactorsTest :
    public TestCase
{
    public:
        OneHalfPlusToThreeHalfMinusFormFactorsTest() :
            TestCase("one_half_plus_to_three_half_minus_form_factor_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Foo->Bar::ABR2022",              parameter, options));
                TEST_CHECK_THROWS(NoSuchFormFactorError, FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda(1520)::FooBar", parameter, options));
            }
        }
} one_half_plus_to_three_half_minus_form_factor_test;
