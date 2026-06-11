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
#include <eos/form-factors/parametric-sse-impl-p-to-v.hh>

using namespace test;
using namespace eos;

class BToKstarSSEFormFactorsTest :
    public TestCase
{
    public:
        BToKstarSSEFormFactorsTest() :
            TestCase("b_to_kstar_sse_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> K^* */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::SSE", p, Options{ });
                TEST_CHECK(ff.get() != nullptr);

                // arbitrary, non-zero coefficients: the kinematic constraints
                //   T_2(0)  = T_1(0)   and   A_12(0) = R * A_0(0)
                // must hold identically, regardless of their values.
                p["B->K^*::alpha^A0_0@SSE" ] =  0.39;
                p["B->K^*::alpha^A0_1@SSE" ] = -1.15;
                p["B->K^*::alpha^A0_2@SSE" ] =  2.08;
                p["B->K^*::alpha^A1_0@SSE" ] =  0.29;
                p["B->K^*::alpha^A1_1@SSE" ] =  0.31;
                p["B->K^*::alpha^A1_2@SSE" ] =  0.72;
                p["B->K^*::alpha^A12_1@SSE"] =  0.57;
                p["B->K^*::alpha^A12_2@SSE"] =  0.14;
                p["B->K^*::alpha^V_0@SSE"  ] =  0.37;
                p["B->K^*::alpha^V_1@SSE"  ] = -1.08;
                p["B->K^*::alpha^V_2@SSE"  ] =  2.47;
                p["B->K^*::alpha^T1_0@SSE" ] =  0.31;
                p["B->K^*::alpha^T1_1@SSE" ] = -0.96;
                p["B->K^*::alpha^T1_2@SSE" ] =  2.01;
                p["B->K^*::alpha^T2_1@SSE" ] =  0.42;
                p["B->K^*::alpha^T2_2@SSE" ] =  2.02;
                p["B->K^*::alpha^T23_0@SSE"] =  0.79;
                p["B->K^*::alpha^T23_1@SSE"] =  1.26;
                p["B->K^*::alpha^T23_2@SSE"] =  1.96;

                // key regression assertions: the constraints must be enforced exactly.
                // These FAIL on the buggy code (which reuses the raw leading coefficients)
                // and PASS once the constant terms absorb the z(0)-dependent difference.
                const double m_B = p["mass::B_d@BSZ2015"];
                const double m_V = p["mass::K_d^*@BSZ2015"];
                const double R   = (power_of<2>(m_B) - power_of<2>(m_V)) / (8.0 * m_B * m_V);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(0.0),  ff->t_1(0.0),     1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->a_12(0.0), R * ff->a_0(0.0), 1.0e-12);

                // pin the parametrisation to reference values computed with the corrected
                // constant terms (the z expansion is not shifted, so z(0) != 0).
                static const double eps = 1.0e-6;
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 5.0), 0.40592825381243042, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_12(5.0), 0.23195540023380720, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 5.0), 0.31972506951478813, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 5.0), 0.25433272732089035, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(  5.0), 0.38600046508416613, eps);
            }
        }
} b_to_kstar_sse_form_factors_test;
