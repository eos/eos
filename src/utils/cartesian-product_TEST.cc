/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
