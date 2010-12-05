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
#include <src/utils/cartesian-product.hh>

#include <complex>

#include <iostream>

using namespace test;
using namespace eos;

class CartesianProductTest :
    public TestCase
{
    public:
        CartesianProductTest() :
            TestCase("cartesian_product_test")
        {
        }

        static void check_tuples(const std::vector<double> & a, const std::vector<double> & b)
        {
            for (auto i = a.cbegin(), i_end = a.cend(), j = b.cbegin() ; i != i_end; ++i, ++j)
            {
                TEST_CHECK_EQUAL(*i, *j);
            }
        }

        virtual void run() const
        {
            /* No Product at all */
            {
                static const std::vector<double> input = { -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
                static const std::vector<std::vector<double>> result = {
                    { -0.5 },
                    { -0.4 },
                    { -0.3 },
                    { -0.2 },
                    { -0.1 },
                    {  0.0 },
                    {  0.1 },
                    {  0.2 },
                    {  0.3 },
                    {  0.4 },
                    {  0.5 },
                };

                CartesianProduct<std::vector<double>> cp;
                cp.over(input);

                auto j = result.cbegin();
                for (auto i = cp.begin() ; cp.end() != i ; ++i, ++j)
                {
                    check_tuples(*i, *j);
                }
            }
        }
} cartesian_product_test;
