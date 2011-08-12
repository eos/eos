/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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
#include <eos/utils/rvalue.hh>

using namespace test;
using namespace eos;

class RValueTest :
    public TestCase
{
    public:
        RValueTest() :
            TestCase("rvalue_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-14;

            // R-value calculation
            // checked against results of markov_chain_sampler.py
            {
                std::vector<double> chain_means{ 4.2, 4.25, 4.22 };
                std::vector<double> chain_variances{ 0.1, 0.15, 0.19 };

                unsigned points = 500;

                TEST_CHECK_RELATIVE_ERROR(RValue::approximation(chain_means, chain_variances, points), 1.0011584199407115, eps);
                TEST_CHECK_RELATIVE_ERROR(RValue::gelman_rubin(chain_means, chain_variances, points), 1.0176292831481546, eps);

                // for more points visited, R-value increases
                points *= 3.0;

                TEST_CHECK_RELATIVE_ERROR(RValue::approximation(chain_means, chain_variances, points), 1.0018240939164496, eps);
                TEST_CHECK_RELATIVE_ERROR(RValue::gelman_rubin(chain_means, chain_variances, points), 1.0183054631320092,eps);
            }
        }
} rvalue_test;
