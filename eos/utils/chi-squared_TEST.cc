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
#include <eos/utils/chi-squared.hh>

#include <cmath>
#include <vector>

using namespace test;
using namespace eos;

struct ChiSquaredInput
{
    double theory_min, theory_central, theory_max;
    double experiment_min, experiment_central, experiment_max;
};

class ChiSquaredTest :
    public TestCase
{
    public:
        ChiSquaredTest() :
            TestCase("chi_squared_test")
        {
        }

        virtual void run() const
        {
            static const std::vector<ChiSquaredInput> inputs
            {
                { 0.2, 0.3, 0.4, -0.1, 0.1, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.1, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.1, 0.7 },

                { 0.2, 0.3, 0.4, -0.1, 0.2, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.2, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.2, 0.7 },

                { 0.2, 0.3, 0.4, -0.1, 0.3, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.3, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.3, 0.7 },

                { 0.2, 0.3, 0.4, -0.1, 0.4, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.4, 0.4 },
                { 0.2, 0.3, 0.4, -0.4, 0.4, 0.7 },
            };

            /* with theory offset */
            {
                static const std::vector<double> results
                {
                    0.04,
                    0.015625,
                    0.0082645,

                    0.0,
                    0.0,
                    0.0,

                    0.0,
                    0.0,
                    0.0,

                    0.0,
                    0.0,
                    0.0,
                };

                auto r = results.cbegin();
                for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i, ++r)
                {
                    static const double eps = 1e-5;
                    TEST_CHECK_NEARLY_EQUAL(*r, ChiSquared::with_theory_offset(i->theory_min, i->theory_central, i->theory_max, i->experiment_min, i->experiment_central, i->experiment_max), eps);
                }
            }

            /* combined-ucnertainties */
            {
                static const std::vector<double> results
                {
                    0.153846,
                    0.0615385,
                    0.0327869,

                    0.0384615,
                    0.0153846,
                    0.00819672,

                    0.0,
                    0.0,
                    0.0,

                    0.0384615,
                    0.0153846,
                    0.00819672,
                };

                auto r = results.cbegin();
                for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i, ++r)
                {
                    static const double eps = 1e-5;
                    TEST_CHECK_NEARLY_EQUAL(*r, ChiSquared::with_combined_uncertainties(i->theory_min, i->theory_central, i->theory_max, i->experiment_min, i->experiment_central, i->experiment_max), eps);
                }
            }
        }
} chi_squared_test;
