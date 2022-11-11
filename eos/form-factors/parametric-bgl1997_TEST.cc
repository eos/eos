/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Christoph Bobeth
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
#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/form-factors/parametric-bgl1997-impl.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BGL1997FormFactorsTest :
    public TestCase
{
    public:
        BGL1997FormFactorsTest() :
            TestCase("BGL1997_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();

            /* Outer function */
            {
                BGL1997FormFactors<BToDstar> ff(p, Options{ });
                const double m_B(BToDstar::m_B);
                const double m_V(BToDstar::m_V);
                const double t_0((m_B + m_V)* power_of<2>(std::sqrt(m_B) - std::sqrt(m_V)));

                //TEST_CHECK_NEARLY_EQUAL( 0.0535063, ff._z(1.0, 2.0), eps);

                //_phi(s, t_0, K, a, b, c, chi);
                TEST_CHECK_NEARLY_EQUAL( 0.0331832,  ff._phi(-2.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0324458,  ff._phi(+1.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0316657,  ff._phi(+4.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.488275,   ff._phi(-2.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.470779,   ff._phi(+1.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.452784,   ff._phi(+4.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.00817026, ff._phi(-2.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00822179, ff._phi(+1.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00827232, ff._phi(+4.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.0928163,  ff._phi(-2.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0934017,  ff._phi(+1.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0939757,  ff._phi(+4.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2), eps);
            }

            /* B -> D^* FFs */
            {
                const double mB(BToDstar::m_B);
                const double mV(BToDstar::m_V);
                const double t_m = (mB - mV) * (mB - mV);
                const double r   = mV / mB, wmax = (mB * mB + mV * mV) / (2.0 * mB * mV);
                const double F2factor = (1.0 + r) / ((1.0 - r) * (1.0 + wmax) * r * mB * mB);
                BGL1997FormFactors<BToDstar> ff(p, Options{ });

                p["B->D^*::a^g_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^g_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^g_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^g_3@BGL1997"] = 0.4e-02;

                p["B->D^*::a^f_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^f_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^f_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^f_3@BGL1997"] = 0.4e-02;

                /* F1_0 parameter determined by identity F1(t_-) = (mB - mV) * f(t_-) */
                p["B->D^*::a^F1_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F1_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F1_3@BGL1997"] = 0.4e-02;

                /* F2_0 parameter determined by identity between F2 and F1 at q2 = 0 */
                p["B->D^*::a^F2_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F2_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F2_3@BGL1997"] = 0.4e-02;

                p["B->D^*::a^T1_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^T1_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^T1_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^T1_3@BGL1997"] = 0.4e-02;

                /* a_T2_0 determined from T1(0) = T2(0) */
                p["B->D^*::a^T2_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^T2_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^T2_3@BGL1997"] = 0.4e-02;

                /* a_T23_0 determined from T2(t_-) = 8 mB mV^2 / ((mB + mV) * (mB^2 + 3 mV^2 - t_-)) T23(t_-) */
                p["B->D^*::a^T23_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^T23_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^T23_3@BGL1997"] = 0.4e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.0120945, ff.g(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0131032, ff.g(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0143205, ff.g(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.598215, ff.f(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.620714, ff.f(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.647914, ff.f(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 2.1225e-4, ff.a_F1_0(), 1.0e-8);
                TEST_CHECK_NEARLY_EQUAL( 3.3597400, ff.F1(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 3.1657300, ff.F1(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.9555500, ff.F1(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 5.8478e-3, ff.a_F2_0(), 1.0e-7);
                TEST_CHECK_NEARLY_EQUAL( 0.2546860, ff.F2(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2803997, ff.F2(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.3115880, ff.F2(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff.F1(t_m), (mB - mV) * ff.f(t_m),  eps);
                TEST_CHECK_NEARLY_EQUAL( ff.F2(0.0), F2factor  * ff.F1(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.0869380, ff.t_1(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0928776, ff.t_1(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1000316, ff.t_1(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 2.4837e-4, ff.a_T2_0(),  1.0e-8);
                TEST_CHECK_NEARLY_EQUAL( 0.0935896, ff.t_2(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0893296, ff.t_2(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0847745, ff.t_2(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 6.3477e-4, ff.a_T23_0(),  1.0e-8);
                TEST_CHECK_NEARLY_EQUAL( 0.0802478, ff.t_23(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0820204, ff.t_23(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0842156, ff.t_23(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0),  ff.t_2(0.0),                                                                      eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(t_m), (mB + mV) * (mB * mB + 3.0 * mV * mV - t_m) / (8.0 * mB * mV * mV) * ff.t_2(t_m), eps);

                p["B->D^*::a^g_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^g_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^g_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^g_3@BGL1997"] = 0.1e-02;

                p["B->D^*::a^f_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^f_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^f_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^f_3@BGL1997"] = 0.1e-02;

                /* F1_0 parameter determined by identity F1(t_-) = (mB - mV) * f(t_-) */
                p["B->D^*::a^F1_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^F1_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F1_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F1_3@BGL1997"] = 0.1e-02;

                /* F2_0 parameter determined by identity between F2 and F1 at q2 = 0 */
                p["B->D^*::a^F2_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F2_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F2_3@BGL1997"] = 0.1e-02;

                p["B->D^*::a^T1_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^T1_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^T1_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^T1_3@BGL1997"] = 0.1e-02;

                /* a_T2_0 determined from T1(0) = T2(0) */
                p["B->D^*::a^T2_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^T2_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^T2_3@BGL1997"] = 0.1e-02;

                /* a_T23_0 determined from T2(t_-) = 8 mB mV^2 / ((mB + mV) * (mB^2 + 3 mV^2 - t_-)) T23(t_-) */
                p["B->D^*::a^T23_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^T23_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^T23_3@BGL1997"] = 0.1e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.0461237, ff.g(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0508857, ff.g(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0566738, ff.g(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 2.281365, ff.f(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.410521, ff.f(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.564134, ff.f(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 7.3850e-4, ff.a_F1_0(), 1.0e-8);
                TEST_CHECK_NEARLY_EQUAL( 9.8438000, ff.F1(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 9.8358614, ff.F1(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 9.8422500, ff.F1(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 1.789e-2, ff.a_F2_0(), 2.0e-6);
                TEST_CHECK_NEARLY_EQUAL( 0.773754, ff.F2(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.854227, ff.F2(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.951860, ff.F2(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff.F1(t_m), (mB - mV) * ff.f(t_m),  eps);
                TEST_CHECK_NEARLY_EQUAL( ff.F2(0.0), F2factor  * ff.F1(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.331549, ff.t_1(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.360687, ff.t_1(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.395877, ff.t_1(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 1.099e-3, ff.a_T2_0(),  1.0e-7);
                TEST_CHECK_NEARLY_EQUAL( 0.347005, ff.t_2(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.352278, ff.t_2(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.358958, ff.t_2(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 3.118e-3, ff.a_T23_0(),  5.0e-7);
                TEST_CHECK_NEARLY_EQUAL( 0.363425, ff.t_23(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.382876, ff.t_23(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.406004, ff.t_23(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0),  ff.t_2(0.0),                                                                      eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(t_m), (mB + mV) * (mB * mB + 3.0 * mV * mV - t_m) / (8.0 * mB * mV * mV) * ff.t_2(t_m), eps);
            }

            /* B -> D FFs*/
            {
                BGL1997FormFactors<BToD> ff(p, Options{ });

                p["B->D::a^f+_0@BGL1997"] = 0.1e-02;
                p["B->D::a^f+_1@BGL1997"] = 0.2e-02;
                p["B->D::a^f+_2@BGL1997"] = 0.3e-02;
                p["B->D::a^f+_3@BGL1997"] = 0.4e-02;

                p["B->D::a^f0_0@BGL1997"] = 0.1e-02;
                p["B->D::a^f0_1@BGL1997"] = 0.2e-02;
                p["B->D::a^f0_2@BGL1997"] = 0.3e-02;
                p["B->D::a^f0_3@BGL1997"] = 0.4e-02;

                p["B->D::a^fT_0@BGL1997"] = 0.1e-02;
                p["B->D::a^fT_1@BGL1997"] = 0.2e-02;
                p["B->D::a^fT_2@BGL1997"] = 0.3e-02;
                p["B->D::a^fT_3@BGL1997"] = 0.4e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.0862919, ff.f_p(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0911714, ff.f_p(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0970192, ff.f_p(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.439279, ff.f_0(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.435522, ff.f_0(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.432496, ff.f_0(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.041750, ff.f_t(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.044758, ff.f_t(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.048361, ff.f_t(+4.0), eps);

                p["B->D::a^f+_0@BGL1997"] = 0.4e-02;
                p["B->D::a^f+_1@BGL1997"] = 0.3e-02;
                p["B->D::a^f+_2@BGL1997"] = 0.2e-02;
                p["B->D::a^f+_3@BGL1997"] = 0.1e-02;

                p["B->D::a^f0_0@BGL1997"] = 0.4e-02;
                p["B->D::a^f0_1@BGL1997"] = 0.3e-02;
                p["B->D::a^f0_2@BGL1997"] = 0.2e-02;
                p["B->D::a^f0_3@BGL1997"] = 0.1e-02;

                p["B->D::a^fT_0@BGL1997"] = 0.4e-02;
                p["B->D::a^fT_1@BGL1997"] = 0.3e-02;
                p["B->D::a^fT_2@BGL1997"] = 0.2e-02;
                p["B->D::a^fT_3@BGL1997"] = 0.1e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.327129, ff.f_p(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.352242, ff.f_p(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.382318, ff.f_p(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 1.66529, ff.f_0(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.68264, ff.f_0(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.70431, ff.f_0(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.158273, ff.f_t(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.172925, ff.f_t(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.190572, ff.f_t(+4.0), eps);
            }
        }
} BGL1997_form_factor_test;
