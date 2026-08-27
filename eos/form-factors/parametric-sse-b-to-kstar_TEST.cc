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

#include <eos/form-factors/parametric-sse-impl-p-to-v.hh>

#include <test/test.hh>

#include <array>
#include <cmath>

using namespace test;
using namespace eos;

class BToKstarSSEFormFactorsTest : public TestCase
{
    public:
        BToKstarSSEFormFactorsTest() :
            TestCase("b_to_kstar_sse_form_factors_test")
        {
        }

        // The SSE P -> V parametrisation has 17 free coefficients: the leading coefficients
        // of A_1, A_12, T_2 and T_23 are fixed by the four kinematic relations and are
        // therefore absent from the parameter list.
        static void
        set_coefficients(Parameters & p, const std::array<double, 17> & a)
        {
            p["B->K^*::alpha^A0_0@SSE"]  = a[0];
            p["B->K^*::alpha^A0_1@SSE"]  = a[1];
            p["B->K^*::alpha^A0_2@SSE"]  = a[2];
            p["B->K^*::alpha^A1_1@SSE"]  = a[3];
            p["B->K^*::alpha^A1_2@SSE"]  = a[4];
            p["B->K^*::alpha^A12_1@SSE"] = a[5];
            p["B->K^*::alpha^A12_2@SSE"] = a[6];
            p["B->K^*::alpha^V_0@SSE"]   = a[7];
            p["B->K^*::alpha^V_1@SSE"]   = a[8];
            p["B->K^*::alpha^V_2@SSE"]   = a[9];
            p["B->K^*::alpha^T1_0@SSE"]  = a[10];
            p["B->K^*::alpha^T1_1@SSE"]  = a[11];
            p["B->K^*::alpha^T1_2@SSE"]  = a[12];
            p["B->K^*::alpha^T2_1@SSE"]  = a[13];
            p["B->K^*::alpha^T2_2@SSE"]  = a[14];
            p["B->K^*::alpha^T23_1@SSE"] = a[15];
            p["B->K^*::alpha^T23_2@SSE"] = a[16];
        }

        // The four kinematic relations
        //   A_12(0)   = R * A_0(0) ,   T_2(0)    = T_1(0) ,
        //   A_12(t_-) = R * A_1(t_-),  T_23(t_-) = K * T_2(t_-)
        // are imposed identically, i.e. they must hold for arbitrary values of the
        // expansion coefficients, not just for the default ones.
        static void
        check_relations(Parameters & p, const FormFactors<PToV> & ff)
        {
            const double m_B = p["mass::B_d@HME"];
            const double m_V = p["mass::K_d^*@HME"];
            const double R   = (power_of<2>(m_B) - power_of<2>(m_V)) / (8.0 * m_B * m_V);
            const double K   = power_of<2>(m_B + m_V) / (4.0 * m_B * m_V);
            const double tm  = power_of<2>(m_B - m_V);

            // the two relations at q^2 = 0
            TEST_CHECK_RELATIVE_ERROR(ff.t_2(0.0), ff.t_1(0.0), 1.0e-12);
            TEST_CHECK_RELATIVE_ERROR(ff.a_12(0.0), R * ff.a_0(0.0), 1.0e-12);

            // the two relations at the endpoint q^2 = t_-, which remove the spurious poles
            // of A_2 and T_3 at t_-, where lambda(q^2) vanishes
            TEST_CHECK_RELATIVE_ERROR(ff.a_12(tm), R * ff.a_1(tm), 1.0e-12);
            TEST_CHECK_RELATIVE_ERROR(ff.t_23(tm), K * ff.t_2(tm), 1.0e-12);

            // ... so that A_2 and T_3 stay finite as q^2 -> t_-. Without the endpoint
            // relations both grow like 1 / (t_- - q^2) and exceed 10^5 already at the
            // first of the points below.
            for (const double & delta : { 1.0e-5, 1.0e-7, 1.0e-9 })
            {
                TEST_CHECK(std::abs(ff.a_2(tm - delta)) < 10.0);
                TEST_CHECK(std::abs(ff.t_3(tm - delta)) < 10.0);
            }
        }

        virtual void
        run() const
        {
            /* B -> K^* */
            {
                Parameters                         p  = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::SSE", p, Options{});
                TEST_CHECK(ff.get() != nullptr);

                // arbitrary, non-zero coefficients
                set_coefficients(p,
                                 {
                                     {
                                      0.39, -1.15,
                                      2.08, // A_0
                                         0.31, 0.72, // A_1
                                         0.57, 0.14, // A_12
                                         0.37, -1.08,
                                      2.47, // V
                                         0.31, -0.96,
                                      2.01, // T_1
                                         0.42, 2.02, // T_2
                                         1.26, 1.96 // T_23
                                     }
                });
                check_relations(p, *ff);

                // pin the parametrisation to reference values computed with all four
                // kinematic relations imposed (the z expansion is not shifted, so
                // z(0) != 0 and z(t_-) != 0).
                static const double eps = 1.0e-6;
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(5.0), 0.38368207018168554, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(5.0), 0.17997546366262482, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_12(5.0), 0.20927027174690238, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(5.0), 0.13550526035327826, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(5.0), 0.37244023813471744, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(5.0), 0.30729000563591313, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(5.0), 0.23012633798064686, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_23(5.0), 0.57505875066164758, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(5.0), 0.15287397791633492, eps);

                // the relations are imposed identically: repeat the checks for two
                // further, unrelated sets of randomised coefficients.
                set_coefficients(p,
                                 {
                                     {
                                      -0.82,
                                      1.73, -0.46, // A_0
                                         2.31, -1.09, // A_1
                                         -1.64,
                                      0.95, // A_12
                                         1.18, -2.07,
                                      0.33, // V
                                         -0.51,
                                      2.66, -1.42, // T_1
                                         0.87, -2.35, // T_2
                                         -1.91,
                                      1.04 // T_23
                                     }
                });
                check_relations(p, *ff);

                set_coefficients(p,
                                 {
                                     {
                                      1.47, -0.28,
                                      3.12, // A_0
                                         -2.64,
                                      0.19, // A_1
                                         0.73, -3.05, // A_12
                                         -1.36,
                                      2.81, -0.62, // V
                                         2.09, -1.77,
                                      0.44, // T_1
                                         -2.18,
                                      1.35, // T_2
                                         0.66, -2.93 // T_23
                                     }
                });
                check_relations(p, *ff);
            }
        }
} b_to_kstar_sse_form_factors_test;
