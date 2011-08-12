/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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
#include <eos/utils/lcda.hh>

#include <cmath>

using namespace test;
using namespace eos;

class LCDATest :
    public TestCase
{
    public:
        LCDATest() :
            TestCase("lcda_test")
        {
        }

        virtual void run() const
        {
            // Gegenbauer evolution
            {
                // Evolve from mu_0 = 2.0 GeV to mu = 4.2 GeV.
                // We test the evolution by using Gegenbauer moments a_n = 1 for all n.

                // eta = alpha_s(mu) / alpha_s(mu_0)
                static const double eta = 0.741005;
                static const double eps = 1e-5;

                TEST_CHECK_RELATIVE_ERROR(LCDA::evolve_gegenbauer_moment(1.0, 1, eta, QCD::beta_function_nf_4), 0.879948, eps);
                TEST_CHECK_RELATIVE_ERROR(LCDA::evolve_gegenbauer_moment(1.0, 2, eta, QCD::beta_function_nf_4), 0.818869, eps);
                TEST_CHECK_RELATIVE_ERROR(LCDA::evolve_gegenbauer_moment(1.0, 3, eta, QCD::beta_function_nf_4), 0.778031, eps);
                TEST_CHECK_RELATIVE_ERROR(LCDA::evolve_gegenbauer_moment(1.0, 4, eta, QCD::beta_function_nf_4), 0.747550, eps);
            }
        }
} lcda_test;
