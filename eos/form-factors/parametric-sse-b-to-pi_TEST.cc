/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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
#include <eos/form-factors/parametric-sse-impl-p-to-p.hh>

using namespace test;
using namespace eos;

class BToPiSSEFormFactorsTest :
    public TestCase
{
    public:
        BToPiSSEFormFactorsTest() :
            TestCase("b_to_pi_sse_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::SSE", p, Options{ });
                TEST_CHECK(ff.get() != nullptr);

                // arbitrary, non-zero coefficients: the equation of motion f_0(0) = f_+(0)
                // must hold identically, regardless of their values.
                p["B->pi::alpha^f+_0@SSE"] =  0.4;
                p["B->pi::alpha^f+_1@SSE"] = -1.1;
                p["B->pi::alpha^f+_2@SSE"] =  0.9;
                p["B->pi::alpha^f0_1@SSE"] =  0.05;
                p["B->pi::alpha^f0_2@SSE"] = -0.37;
                p["B->pi::alpha^fT_0@SSE"] =  0.27;
                p["B->pi::alpha^fT_1@SSE"] = -1.4;
                p["B->pi::alpha^fT_2@SSE"] = -1.3;

                // key regression assertion: the equation of motion must be enforced exactly.
                // This FAILS on the buggy code (which reuses the raw coefficient a^f+_0) and
                // PASSES once the constant term of f_0 absorbs the z(0)-dependent difference.
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(0.0), ff->f_p(0.0), 1.0e-12);

                // pin the whole parametrisation to reference values computed with the
                // corrected f_0 (the z expansion is not shifted, so z(0) != 0).
                static const double eps = 1.0e-6;
                TEST_CHECK_RELATIVE_ERROR(ff->f_p( 0.0),  0.16058302054963150, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0( 0.0),  0.16058302054963147, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p(10.0),  0.35133196523016341, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(10.0),  0.25619406664625294, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t(10.0), -0.05125633612234586, eps);
            }
        }
} b_to_pi_sse_form_factors_test;
