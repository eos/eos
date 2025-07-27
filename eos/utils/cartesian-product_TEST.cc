/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
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

#include <eos/utils/cartesian-product.hh>

#include <test/test.hh>

#include <complex>
#include <iostream>

using namespace test;
using namespace eos;

class CartesianProductTest : public TestCase
{
    public:
        CartesianProductTest() :
            TestCase("cartesian_product_test")
        {
        }

        static void
        check_tuples(const std::vector<double> & a, const std::vector<double> & b)
        {
            for (auto i = a.cbegin(), i_end = a.cend(), j = b.cbegin(); i != i_end; ++i, ++j)
            {
                TEST_CHECK_EQUAL(*i, *j);
            }
        }

        virtual void
        run() const
        {
            /* No Product at all */
            {
                static const std::vector<double>              input  = { -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
                static const std::vector<std::vector<double>> result = {
                    { -0.5 }, { -0.4 }, { -0.3 }, { -0.2 }, { -0.1 }, { 0.0 }, { 0.1 }, { 0.2 }, { 0.3 }, { 0.4 }, { 0.5 },
                };

                CartesianProduct<std::vector<double>> cp;
                cp.over(input);

                auto j = result.cbegin();
                for (auto i = cp.begin(); cp.end() != i; ++i, ++j)
                {
                    check_tuples(*i, *j);
                }
            }

            {
                static const std::vector<double> input1 = { 1.0, 2.0 };
                static const std::vector<double> input2 = { 10.0, 20.0, 30.0 };
                static const std::vector<double> input3 = { 100.0, 200.0, 300.0, 400.0 };
                static const std::vector<double> input4 = { 1000.0, 2000.0 };

                CartesianProduct<std::vector<double>> cp;
                cp.over(input1);
                cp.over(input2);
                cp.over(input3);
                cp.over(input4);

                auto cp_it = cp.begin();

                TEST_CHECK_EQUAL(cp_it, cp.begin());

                static const std::vector<double> result1 = { 1.0, 10.0, 100.0, 1000.0 };
                TEST_CHECK_EQUAL(*cp_it, result1);

                ++cp_it;
                static const std::vector<double> result2 = { 1.0, 10.0, 100.0, 2000.0 };
                TEST_CHECK_EQUAL(*cp_it, result2);

                cp_it                                    += 10;
                static const std::vector<double> result3  = { 1.0, 20.0, 200.0, 2000.0 };
                TEST_CHECK_EQUAL(*cp_it, result3);

                if (cp_it == cp.begin())
                {
                    TEST_CHECK_FAILED("cp_it should not be equal cp.begin()");
                }

                if (cp_it == cp.end())
                {
                    TEST_CHECK_FAILED("cp_it should not be equal cp.end()");
                }

                cp_it += 35;

                if (cp_it == cp.end())
                {
                    TEST_CHECK_FAILED("cp_it should not be equal cp.end()");
                }

                static const std::vector<double> result4 = { 2.0, 30.0, 400.0, 1000.0 };
                TEST_CHECK_EQUAL(*cp_it, result4);

                cp_it += 2;
                TEST_CHECK_EQUAL(cp_it, cp.end());

                ++cp_it;
                TEST_CHECK_EQUAL(cp_it, cp.end());
            }
        }
} cartesian_product_test;
