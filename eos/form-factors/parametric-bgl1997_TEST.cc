/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Christoph Bobeth
 * Copyright (c) 2023 Nico Gubernari
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

            /* Outer function: t_0 = t_0,opt */
            {
                BGL1997FormFactors<BToDstar, PToV> ff(p, Options{ });
                const double m_B(BToDstar::m_B);
                const double m_V(BToDstar::m_V);
                const double t_0((m_B + m_V)* power_of<2>(std::sqrt(m_B) - std::sqrt(m_V)));

                p["mass::B_c^*@BSZ2015"] = 6.329;
                p["mass::B_c^*[1]@BSZ2015"] = 6.910;
                p["mass::B_c^*[2]@BSZ2015"] = 7.020;
                p["mass::B_c,1@BSZ2015"] = 6.739;
                p["mass::B_c,1[1]@BSZ2015"] = 6.750;
                p["mass::B_c,1[2]@BSZ2015"] = 7.145;
                p["mass::B_c,1[3]@BSZ2015"] = 7.150;
                p["mass::B_c@BSZ2015"] = 6.275;
                p["mass::B_c[1]@BSZ2015"] = 6.871;
                p["mass::B_c[2]@BSZ2015"] = 7.250;
                p["mass::B_c,0@BSZ2015"] = 6.704;
                p["mass::B_c,0[1]@BSZ2015"] = 7.122;

                p["b->c::chiOPE[1^-_V]"] = 5.131e-04;
                p["b->c::chiOPE[0^+_V]"] = 6.204e-03;
                p["b->c::chiOPE[1^+_A]"] = 3.894e-04;
                p["b->c::chiOPE[0^-_A]"] = 19.421e-03;
                p["b->c::chiOPE[1^-_T]"] = 0.0004897959183673469;
                p["b->c::chiOPE[1^+_T5]"] = 0.0002715419501133787;

                //TEST_CHECK_NEARLY_EQUAL(ff._z(1.0, 2.0), 0.0535063, eps);

                //_phi(s, t_0, K, a, b, c, chi);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 48, 3, 3, 2, 3.1e-03),  0.0331832,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 48, 3, 3, 2, 3.1e-03),  0.0324458,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 48, 3, 3, 2, 3.1e-03),  0.0316657,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 48, 3, 3, 1, 3.1e-03),  0.488275,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 48, 3, 3, 1, 3.1e-03),  0.470779,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 48, 3, 3, 1, 3.1e-03),  0.452784,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 16, 1, 1, 1, 3.1e-03),  0.00817026, eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 16, 1, 1, 1, 3.1e-03),  0.00822179, eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 16, 1, 1, 1, 3.1e-03),  0.00827232, eps);

                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2),  0.0928163,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2),  0.0934017,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2),  0.0939757,  eps);
            }

            /* B -> D^* FFs */
            {
                const double mB(BToDstar::m_B);
                const double mV(BToDstar::m_V);
                const double t_m = (mB - mV) * (mB - mV);
                const double r   = mV / mB, wmax = (mB * mB + mV * mV) / (2.0 * mB * mV);
                const double F2factor = (1.0 + r) / ((1.0 - r) * (1.0 + wmax) * r * mB * mB);
                BGL1997FormFactors<BToDstar, PToV> ff(p, Options{ });

                p["mass::B_c^*@BSZ2015"] = 6.329;
                p["mass::B_c^*[1]@BSZ2015"] = 6.910;
                p["mass::B_c^*[2]@BSZ2015"] = 7.020;
                p["mass::B_c,1@BSZ2015"] = 6.739;
                p["mass::B_c,1[1]@BSZ2015"] = 6.750;
                p["mass::B_c,1[2]@BSZ2015"] = 7.145;
                p["mass::B_c,1[3]@BSZ2015"] = 7.150;
                p["mass::B_c@BSZ2015"] = 6.275;
                p["mass::B_c[1]@BSZ2015"] = 6.871;
                p["mass::B_c[2]@BSZ2015"] = 7.250;
                p["mass::B_c,0@BSZ2015"] = 6.704;
                p["mass::B_c,0[1]@BSZ2015"] = 7.122;

                p["b->c::chiOPE[1^-_V]"] = 5.131e-04;
                p["b->c::chiOPE[0^+_V]"] = 6.204e-03;
                p["b->c::chiOPE[1^+_A]"] = 3.894e-04;
                p["b->c::chiOPE[0^-_A]"] = 19.421e-03;
                p["b->c::chiOPE[1^-_T]"] = 0.0004897959183673469;
                p["b->c::chiOPE[1^+_T5]"] = 0.0002715419501133787;

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

                TEST_CHECK_NEARLY_EQUAL(ff.g(-2.0),  0.0120945, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+1.0),  0.0131032, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+4.0),  0.0143205, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f(-2.0),  0.598215,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+1.0),  0.620714,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+4.0),  0.647914,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F1_0(),  2.1225e-4, 1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(-2.0),  3.3597400, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+1.0),  3.1657300, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+4.0),  2.9555500, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F2_0(),  5.8478e-3, 1.0e-7);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(-2.0),  0.2546860, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+1.0),  0.2803997, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+4.0),  0.3115880, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.F1(t_m), (mB - mV) * ff.f(t_m),  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(0.0), F2factor  * ff.F1(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(-2.0),  0.0869380, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+1.0),  0.0928776, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+4.0),  0.1000316, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T2_0(),   2.4837e-4, 1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(-2.0),  0.0935896, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+1.0),  0.0893296, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+4.0),  0.0847745, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T23_0(),  6.3477e-4, 1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(-2.0), 0.0802478, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+1.0), 0.0820204, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+4.0), 0.0842156, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(0.0),   ff.t_2(0.0),                                                                     eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(t_m), (mB + mV) * (mB * mB + 3.0 * mV * mV - t_m) / (8.0 * mB * mV * mV) * ff.t_2(t_m), eps);

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

                TEST_CHECK_NEARLY_EQUAL(ff.g(-2.0),  0.0461237, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+1.0),  0.0508857, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+4.0),  0.0566738, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f(-2.0),  2.281365,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+1.0),  2.410521,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+4.0),  2.564134,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F1_0(), 7.3850e-4, 1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(-2.0), 9.8438000, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+1.0), 9.8358614, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+4.0), 9.8422500, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F2_0(), 1.789e-2,  2.0e-6);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(-2.0), 0.773754,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+1.0), 0.854227,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+4.0), 0.951860,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.F1(t_m), (mB - mV) * ff.f(t_m),  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(0.0), F2factor  * ff.F1(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(-2.0),  0.331549, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+1.0),  0.360687, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+4.0),  0.395877, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T2_0(),   1.099e-3, 1.0e-7);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(-2.0),  0.347005, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+1.0),  0.352278, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+4.0),  0.358958, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T23_0(),  3.118e-3, 5.0e-7);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(-2.0), 0.363425, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+1.0), 0.382876, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+4.0), 0.406004, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(0.0),  ff.t_2(0.0),                                                                      eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(t_m), (mB + mV) * (mB * mB + 3.0 * mV * mV - t_m) / (8.0 * mB * mV * mV) * ff.t_2(t_m), eps);
            }

            /* B -> D FFs*/
            {
                BGL1997FormFactors<BToD, PToP> ff(p, Options{ });

                p["mass::B_c^*@BSZ2015"] = 6.329;
                p["mass::B_c^*[1]@BSZ2015"] = 6.910;
                p["mass::B_c^*[2]@BSZ2015"] = 7.020;
                p["mass::B_c,0@BSZ2015"] = 6.704;
                p["mass::B_c,0[1]@BSZ2015"] = 7.122;

                p["b->c::chiOPE[1^-_V]"] = 5.131e-04;
                p["b->c::chiOPE[0^+_V]"] = 6.204e-03;
                p["b->c::chiOPE[1^-_T]"] = 0.0004897959184;

                p["mass::D_u@BSZ2015"] = 1.870;

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

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-2.0), 0.0862919, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+1.0), 0.0911714, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+4.0), 0.0970192, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0(-2.0), 0.439279,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+1.0), 0.435522,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+4.0), 0.432496,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t(-2.0), 0.041750,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+1.0), 0.044758,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+4.0), 0.048361,  eps);

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

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-2.0), 0.327129, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+1.0), 0.352242, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+4.0), 0.382318, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0(-2.0), 1.66529,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+1.0), 1.68264,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+4.0), 1.70431,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t(-2.0), 0.158273, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+1.0), 0.172925, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+4.0), 0.190572, eps);
            }

            /* Adapt the parameters */
            p["B->D::t_0@BGL1997"]  = 11.6213;
            p["B->D^*::t_0@BGL1997"] = 10.6844;

            /* Outer function: t_0 = t_- */
            {
                BGL1997FormFactors<BToDstar, PToV> ff(p, Options{ });
                const double t_0(10.684);

                p["mass::B_c^*@BSZ2015"] = 6.329;
                p["mass::B_c^*[1]@BSZ2015"] = 6.910;
                p["mass::B_c^*[2]@BSZ2015"] = 7.020;
                p["mass::B_c,1@BSZ2015"] = 6.739;
                p["mass::B_c,1[1]@BSZ2015"] = 6.750;
                p["mass::B_c,1[2]@BSZ2015"] = 7.145;
                p["mass::B_c,1[3]@BSZ2015"] = 7.150;
                p["mass::B_c@BSZ2015"] = 6.275;
                p["mass::B_c[1]@BSZ2015"] = 6.871;
                p["mass::B_c[2]@BSZ2015"] = 7.250;
                p["mass::B_c,0@BSZ2015"] = 6.704;
                p["mass::B_c,0[1]@BSZ2015"] = 7.122;

                p["b->c::chiOPE[1^-_V]"] = 5.131e-04;
                p["b->c::chiOPE[0^+_V]"] = 6.204e-03;
                p["b->c::chiOPE[1^+_A]"] = 3.894e-04;
                p["b->c::chiOPE[0^-_A]"] = 19.421e-03;
                p["b->c::chiOPE[1^-_T]"] = 0.0004897959183673469;
                p["b->c::chiOPE[1^+_T5]"] = 0.0002715419501133787;

                //_phi(s, t_0, K, a, b, c, chi);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 48, 3, 3, 2, 3.1e-03),                  0.0332310,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 48, 3, 3, 2, 3.1e-03),                  0.0324799,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 48, 3, 3, 2, 3.1e-03),                  0.0316857,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 48, 3, 3, 1, 3.1e-03),                  0.488978,    eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 48, 3, 3, 1, 3.1e-03),                  0.471272,    eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 48, 3, 3, 1, 3.1e-03),                  0.4530700,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 16, 1, 1, 1, 3.1e-03),                  0.00818202,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 16, 1, 1, 1, 3.1e-03),                  0.00823041,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 16, 1, 1, 1, 3.1e-03),                  0.00827755,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff._phi(-2.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2), 0.0929520,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+1.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2), 0.0935017,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff._phi(+4.0, t_0, 1.4153, 1, 1, 1, 4.79e-03 / 4.2 / 4.2), 0.0940372,   eps);
            }

            /* B -> D^* FFs */
            {
                const double mB(BToDstar::m_B);
                const double mV(BToDstar::m_V);
                const double t_m = (mB - mV) * (mB - mV);
                const double r   = mV / mB, wmax = (mB * mB + mV * mV) / (2.0 * mB * mV);
                const double F2factor = (1.0 + r) / ((1.0 - r) * (1.0 + wmax) * r * mB * mB);
                BGL1997FormFactors<BToDstar, PToV> ff(p, Options{ });

                p["mass::B_c^*@BSZ2015"] = 6.329;
                p["mass::B_c^*[1]@BSZ2015"] = 6.910;
                p["mass::B_c^*[2]@BSZ2015"] = 7.020;
                p["mass::B_c,1@BSZ2015"] = 6.739;
                p["mass::B_c,1[1]@BSZ2015"] = 6.750;
                p["mass::B_c,1[2]@BSZ2015"] = 7.145;
                p["mass::B_c,1[3]@BSZ2015"] = 7.150;
                p["mass::B_c@BSZ2015"] = 6.275;
                p["mass::B_c[1]@BSZ2015"] = 6.871;
                p["mass::B_c[2]@BSZ2015"] = 7.250;
                p["mass::B_c,0@BSZ2015"] = 6.704;
                p["mass::B_c,0[1]@BSZ2015"] = 7.122;

                p["b->c::chiOPE[1^-_V]"] = 5.131e-04;
                p["b->c::chiOPE[0^+_V]"] = 6.204e-03;
                p["b->c::chiOPE[1^+_A]"] = 3.894e-04;
                p["b->c::chiOPE[0^-_A]"] = 19.421e-03;
                p["b->c::chiOPE[1^-_T]"] = 0.0004897959183673469;
                p["b->c::chiOPE[1^+_T5]"] = 0.0002715419501133787;

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

                TEST_CHECK_NEARLY_EQUAL(ff.g(-2.0), 0.0128101,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+1.0), 0.0138737,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+4.0), 0.0151567,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f(-2.0), 0.633612,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+1.0), 0.657215,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+4.0), 0.685746,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F1_0(), 1.6743e-4,  1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(-2.0), 3.5930000,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+1.0), 3.3807000,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+4.0), 3.1500100,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F2_0(), 6.1959e-3,  1.0e-7);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(-2.0), 0.2720690,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+1.0), 0.2996300,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+4.0), 0.3330650,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.F1(t_m), (mB - mV) * ff.f(t_m),  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(0.0), F2factor  * ff.F1(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(-2.0), 0.0920823,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+1.0), 0.0983392,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+4.0), 0.1058720,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T2_0(),  2.0313e-4,  1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(-2.0), 0.0992379,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+1.0), 0.0945217,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+4.0), 0.0894548,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T23_0(),  6.0662e-4,  1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(-2.0), 0.0843247,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+1.0), 0.0860918,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+4.0), 0.0882836,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(0.0),  ff.t_2(0.0),                                                                      eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(t_m), (mB + mV) * (mB * mB + 3.0 * mV * mV - t_m) / (8.0 * mB * mV * mV) * ff.t_2(t_m), eps);

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

                TEST_CHECK_NEARLY_EQUAL(ff.g(-2.0), 0.0470640,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+1.0), 0.0519358,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.g(+4.0), 0.0578580,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f(-2.0), 2.327870,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+1.0), 2.460270,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f(+4.0), 2.61771,    eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F1_0(), 6.6972e-4,  1.0e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(-2.0), 10.073300,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+1.0), 10.063300,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F1(+4.0), 10.066900,  10*eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_F2_0(), 1.8241e-2,  2.0e-6);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(-2.0), 0.791514,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+1.0), 0.874152,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(+4.0), 0.974439,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff.F1(t_m), (mB - mV) * ff.f(t_m),  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F2(0.0), F2factor  * ff.F1(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(-2.0), 0.338307,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+1.0), 0.368130,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_1(+4.0), 0.404149,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T2_0(),  1.0357e-3,  1.0e-7);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(-2.0), 0.354159,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+1.0), 0.359504,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_2(+4.0), 0.366256,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff.a_T23_0(),  3.093e-3,  5.0e-7);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(-2.0), 0.370258,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+1.0), 0.390147,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(+4.0), 0.413790,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.t_1(0.0),  ff.t_2(0.0),                                                                      eps);
                TEST_CHECK_NEARLY_EQUAL(ff.t_23(t_m), (mB + mV) * (mB * mB + 3.0 * mV * mV - t_m) / (8.0 * mB * mV * mV) * ff.t_2(t_m), eps);
            }

            /* B -> D FFs*/
            {
                BGL1997FormFactors<BToD, PToP> ff(p, Options{ });

                p["mass::B_c^*@BSZ2015"] = 6.329;
                p["mass::B_c^*[1]@BSZ2015"] = 6.910;
                p["mass::B_c^*[2]@BSZ2015"] = 7.020;
                p["mass::B_c,0@BSZ2015"] = 6.704;
                p["mass::B_c,0[1]@BSZ2015"] = 7.122;

                p["b->c::chiOPE[1^-_V]"] = 5.131e-04;
                p["b->c::chiOPE[0^+_V]"] = 6.204e-03;
                p["b->c::chiOPE[1^-_T]"] = 0.0004897959183673469;

                p["mass::D_u@BSZ2015"] = 1.870;

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

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-2.0), 0.0922009,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+1.0), 0.0973758,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+4.0), 0.103574,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0(-2.0), 0.469359,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+1.0), 0.465160,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+4.0), 0.461717,   eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t(-2.0), 0.044609,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+1.0), 0.047804,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+4.0), 0.051628,   eps);

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

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-2.0), 0.334757,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+1.0), 0.360564,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(+4.0), 0.391470,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0(-2.0), 1.70412,   eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+1.0), 1.722398,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(+4.0), 1.745110,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t(-2.0), 0.161964,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+1.0), 0.177010,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(+4.0), 0.195134,  eps);
            }
        }
} BGL1997_form_factor_test;
