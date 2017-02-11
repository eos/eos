/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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
#include <eos/utils/apply.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

namespace
{
    struct TestPointerToMemberFunction
    {
        double _x, _y;

        double nullary_function()
        {
            _x = -1.0;
            _y = -7.0;

            return M_PI;
        }

        void binary_function(const double & x, const double & y)
        {
            _x = x;
            _y = y;
        }
    };

    static double static_nullary_function()
    {
        return std::exp(1.0);
    }

    static double static_unary_function(const double & x)
    {
        return x;
    }

    static double static_ternary_function(const double & x, const double & y, const double & z)
    {
        return x * x + y * y + z * z - 2.0 * (x * y + y * z + z * x);
    }
}

class ApplyTupleTest :
    public TestCase
{
    public:
        ApplyTupleTest() :
            TestCase("apply_tuple_test")
        {
        }

        virtual void run() const
        {
            // Test pointer to memeber function
            {
                TestPointerToMemberFunction test;

                auto result = apply(&::TestPointerToMemberFunction::nullary_function, std::make_tuple(&test));
                TEST_CHECK_EQUAL(result, M_PI);
                TEST_CHECK_EQUAL(test._x, -1.0);
                TEST_CHECK_EQUAL(test._y, -7.0);

                apply(&::TestPointerToMemberFunction::binary_function, std::make_tuple(&test, 1.0, 2.0));
                TEST_CHECK_EQUAL(test._x, 1.0);
                TEST_CHECK_EQUAL(test._y, 2.0);
            }

            // Test static free-standing unary function
            {
                auto result = apply(&::static_unary_function, std::make_tuple(M_PI));

                TEST_CHECK_EQUAL(result, M_PI);
            }

            // Test static free-standing ternay function
            {
                auto result = apply(&::static_ternary_function, std::make_tuple(1.0, 0.5, 0.5));

                TEST_CHECK_EQUAL(result, -1.0);
            }

            // Test std::function<>-wrapped static member function
            {
                std::function<double ()> f(&::static_nullary_function);
                auto result = apply(f, std::make_tuple());

                TEST_CHECK_EQUAL(result, std::exp(1.0));
            }
        }
} apply_tuple_test;

class ApplyArrayTest :
    public TestCase
{
    public:
        ApplyArrayTest() :
            TestCase("apply_array_test")
        {
        }

        virtual void run() const
        {
            // Test pointer to memeber function
            {
                TestPointerToMemberFunction test;

                auto result = apply(&::TestPointerToMemberFunction::nullary_function, &test, std::array<double, 0ul>());

                TEST_CHECK_EQUAL(result, M_PI);
                TEST_CHECK_EQUAL(test._x, -1.0);
                TEST_CHECK_EQUAL(test._y, -7.0);

                apply(&::TestPointerToMemberFunction::binary_function, &test, std::array<double, 2ul>{ 1.0, 2.0 });

                TEST_CHECK_EQUAL(test._x, 1.0);
                TEST_CHECK_EQUAL(test._y, 2.0);
            }

            // Test std::function<>-wrapped pointer to memeber function
            {
                TestPointerToMemberFunction test;

                auto f_nullary = std::function<double (::TestPointerToMemberFunction *)>(&::TestPointerToMemberFunction::nullary_function);
                auto result = apply(f_nullary, &test, std::array<double, 0ul>());

                TEST_CHECK_EQUAL(result, M_PI);
                TEST_CHECK_EQUAL(test._x, -1.0);
                TEST_CHECK_EQUAL(test._y, -7.0);

                auto f_binary = std::function<void (::TestPointerToMemberFunction *, const double &, const double &)>(&::TestPointerToMemberFunction::binary_function);
                apply(f_binary, &test, std::array<double, 2ul>{ 1.0, 2.0 });

                TEST_CHECK_EQUAL(test._x, 1.0);
                TEST_CHECK_EQUAL(test._y, 2.0);
            }

            // Test static free-standing unary function
            {
                auto result = apply(&::static_unary_function, std::array<double, 1ul>{ M_PI });

                TEST_CHECK_EQUAL(result, M_PI);
            }

            // Test static free-standing ternay function
            {
                auto result = apply(&::static_ternary_function, std::array<double, 3ul>{ 1.0, 0.5, 0.5 });

                TEST_CHECK_EQUAL(result, -1.0);
            }

            // Test std::function<>-wrapped static free-standing function
            {
                std::function<double ()> f(&::static_nullary_function);
                auto result = apply(f, std::array<double, 0ul>());

                TEST_CHECK_EQUAL(result, std::exp(1.0));
            }
        }
} apply_array_test;
