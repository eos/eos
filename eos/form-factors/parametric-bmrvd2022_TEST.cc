/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2022 Danny van Dyk
 * Copyright (c) 2021-2022 Muslem Rahimi
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
#include <eos/form-factors/parametric-bmrvd2022.hh>
#include <eos/form-factors/parametric-bmrvd2022-impl.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BMRvD2022FormFactorsTest :
    public TestCase
{
    public:
        BMRvD2022FormFactorsTest() :
            TestCase("BMRvD2022_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            /* Lb -> L FFs */
            {
                Parameters p = Parameters::Defaults();
                // f_time_v
                p["Lambda_b->Lambda::a^(t,V)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(t,V)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(t,V)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(t,V)_4@BMRvD2022"] = +0.5;
                // f_long_v
                p["Lambda_b->Lambda::a^(0,V)_0@BMRvD2022"] = +0.1;
                p["Lambda_b->Lambda::a^(0,V)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(0,V)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(0,V)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(0,V)_4@BMRvD2022"] = +0.5;
                // f_perp_v
                p["Lambda_b->Lambda::a^(perp,V)_0@BMRvD2022"] = +0.1;
                p["Lambda_b->Lambda::a^(perp,V)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(perp,V)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(perp,V)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(perp,V)_4@BMRvD2022"] = +0.5;
                // f_time_a
                p["Lambda_b->Lambda::a^(t,A)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(t,A)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(t,A)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(t,A)_4@BMRvD2022"] = +0.5;
                // f_long_a
                p["Lambda_b->Lambda::a^(0,A)_0@BMRvD2022"] = +0.1;
                p["Lambda_b->Lambda::a^(0,A)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(0,A)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(0,A)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(0,A)_4@BMRvD2022"] = +0.5;
                // f_perp_a
                p["Lambda_b->Lambda::a^(perp,A)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(perp,A)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(perp,A)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(perp,A)_4@BMRvD2022"] = +0.5;
                // f_long_t
                p["Lambda_b->Lambda::a^(0,T)_0@BMRvD2022"] = +0.1;
                p["Lambda_b->Lambda::a^(0,T)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(0,T)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(0,T)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(0,T)_4@BMRvD2022"] = +0.5;
                // f_perp_t
                p["Lambda_b->Lambda::a^(perp,T)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(perp,T)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(perp,T)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(perp,T)_4@BMRvD2022"] = +0.5;
                // f_long_t5
                p["Lambda_b->Lambda::a^(0,T5)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(0,T5)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(0,T5)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(0,T5)_4@BMRvD2022"] = +0.5;
                // f_perp_t5
                p["Lambda_b->Lambda::a^(perp,T5)_0@BMRvD2022"] = +0.1;
                p["Lambda_b->Lambda::a^(perp,T5)_1@BMRvD2022"] = -0.2;
                p["Lambda_b->Lambda::a^(perp,T5)_2@BMRvD2022"] = +0.3;
                p["Lambda_b->Lambda::a^(perp,T5)_3@BMRvD2022"] = -0.4;
                p["Lambda_b->Lambda::a^(perp,T5)_4@BMRvD2022"] = +0.5;
                // Resonance masses
                p["mass::B_s@BSZ2015"]   = 5.367;
                p["mass::B_s,0@BSZ2015"] = 5.711;
                p["mass::B_s^*@BSZ2015"] = 5.416;
                p["mass::B_s,1@BSZ2015"] = 5.750;
                // Fix tp_a to tp_v to match the initial publication [BMRvD:2022A]
                p["Lambda_b->Lambda::tp_a@BMRvD2022"] = p["Lambda_b->Lambda::tp_v@BMRvD2022"].evaluate();


                BMRvD2022FormFactors<LambdaBToLambda> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( 0.230324,  eps), // z(q2 =  0)
                    std::make_pair( 0.107523,  eps), // phi(q2 =  0, f_time^V)
                    std::make_pair( 0.0511848, eps), // phi(q2 =  0, f_long^V)
                    std::make_pair( 0.12409,   eps), // phi(q2 =  0, f_perp^V)
                    std::make_pair( 0.213071,  eps), // phi(q2 =  0, f_time^A)
                    std::make_pair( 0.0253141, eps), // phi(q2 =  0, f_long^A)
                    std::make_pair( 0.0917756, eps), // phi(q2 =  0, f_perp^A)
                    std::make_pair( 0.107298,  eps), // phi(q2 =  0, f_long^T)
                    std::make_pair( 0.0885165, eps), // phi(q2 =  0, f_perp^T)
                    std::make_pair( 0.0797895, eps), // phi(q2 =  0, f_long^T5)
                    std::make_pair( 0.044016,  eps), // phi(q2 =  0, f_perp^T5)
                    std::make_pair( 0.144329,  eps), // z(q2 = 10)
                    std::make_pair( 0.0978787, eps), // phi(q2 = 10, f_time^V)
                    std::make_pair( 0.0516911, eps), // phi(q2 = 10, f_long^V)
                    std::make_pair( 0.115081,  eps), // phi(q2 = 10, f_perp^V)
                    std::make_pair( 0.197602,  eps), // phi(q2 = 10, f_time^A)
                    std::make_pair( 0.0250934, eps), // phi(q2 = 10, f_long^A)
                    std::make_pair( 0.083544,  eps), // phi(q2 = 10, f_perp^A)
                    std::make_pair( 0.0995077, eps), // phi(q2 = 10, f_long^T)
                    std::make_pair( 0.089392,  eps), // phi(q2 = 10, f_perp^T)
                    std::make_pair( 0.0726330, eps), // phi(q2 = 10, f_long^T5)
                    std::make_pair( 0.0436322, eps), // phi(q2 = 10, f_perp^T5)
                    std::make_pair( 0.557107,  eps), // p_0(z = 0.0)
                    std::make_pair(-0.440501,  eps), // p_1(z = 0.0)
                    std::make_pair( 0.633591,  eps), // p_2(z = 0.0)
                    std::make_pair(-0.884426,  eps), // p_3(z = 0.0)
                    std::make_pair( 1.22601,   eps), // p_4(z = 0.0)
                    std::make_pair(-1.69717,   eps), // p_5(z = 0.0)
                    std::make_pair( 0.557107,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair(-0.337996,  eps), // p_1(z = z(q2 = 10))
                    std::make_pair( 0.511499,  eps), // p_2(z = z(q2 = 10))
                    std::make_pair(-0.664035,  eps), // p_3(z = z(q2 = 10))
                    std::make_pair( 0.863691,  eps), // p_4(z = z(q2 = 10))
                    std::make_pair(-1.12016,   eps), // p_5(z = z(q2 = 10))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // form factor values
                TEST_CHECK_NEARLY_EQUAL(ff.f_time_v ( 0.0), 33.25224570, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_v ( 0.0), 33.25224570, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_v ( 0.0), 13.71590521, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_time_a ( 0.0), 39.04787827, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_a ( 0.0), 39.04787827, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_a ( 0.0), 43.76278718, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_t ( 0.0), 15.86249033, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_t ( 0.0), 22.45684904, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_t5( 0.0), 21.61605973, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_t5( 0.0), 22.45684904, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_time_v (10.0), 40.87968764, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_v (10.0), 45.45565962, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_v (10.0), 20.41738244, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_time_a (10.0), 51.93704349, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_a (10.0), 48.08756014, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_a (10.0), 51.98391863, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_t (10.0), 23.61277120, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_t (10.0), 30.02964424, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_t5(10.0), 27.11308541, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_t5(10.0), 27.65566599, eps);

                const auto tm = power_of<2>(LambdaBToLambda::m1 - LambdaBToLambda::m2);
                TEST_CHECK_NEARLY_EQUAL(ff.f_time_v (tm),   59.87052714, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_v (tm),   88.45600970, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_v (tm),   44.89007817, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_time_a (tm),   96.34996613, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_a (tm),   70.92712078, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_a (tm),   70.92712078, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_t (tm),   51.91552581, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_t (tm),   56.60821914, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_long_t5(tm),   40.79093961, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_perp_t5(tm),   40.79093961, eps);
            }


            /* Lc -> L FFs */
            {
                Parameters p = Parameters::Defaults();
                // f_time_v
                p["Lambda_c->Lambda::a^(t,V)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(t,V)_2@BMRvD2022"] = +0.3;
                // f_long_v
                p["Lambda_c->Lambda::a^(0,V)_0@BMRvD2022"] = +0.1;
                p["Lambda_c->Lambda::a^(0,V)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(0,V)_2@BMRvD2022"] = +0.3;
                // f_perp_v
                p["Lambda_c->Lambda::a^(perp,V)_0@BMRvD2022"] = +0.1;
                p["Lambda_c->Lambda::a^(perp,V)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(perp,V)_2@BMRvD2022"] = +0.3;
                // f_time_a
                p["Lambda_c->Lambda::a^(t,A)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(t,A)_2@BMRvD2022"] = +0.3;
                // f_long_a
                p["Lambda_c->Lambda::a^(0,A)_0@BMRvD2022"] = +0.1;
                p["Lambda_c->Lambda::a^(0,A)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(0,A)_2@BMRvD2022"] = +0.3;
                // f_perp_a
                p["Lambda_c->Lambda::a^(perp,A)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(perp,A)_2@BMRvD2022"] = +0.3;
                // f_long_t
                p["Lambda_c->Lambda::a^(0,T)_0@BMRvD2022"] = +0.1;
                p["Lambda_c->Lambda::a^(0,T)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(0,T)_2@BMRvD2022"] = +0.3;
                // f_perp_t
                p["Lambda_c->Lambda::a^(perp,T)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(perp,T)_2@BMRvD2022"] = +0.3;
                // f_long_t5
                p["Lambda_c->Lambda::a^(0,T5)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(0,T5)_2@BMRvD2022"] = +0.3;
                // f_perp_t5
                p["Lambda_c->Lambda::a^(perp,T5)_0@BMRvD2022"] = +0.1;
                p["Lambda_c->Lambda::a^(perp,T5)_1@BMRvD2022"] = -0.2;
                p["Lambda_c->Lambda::a^(perp,T5)_2@BMRvD2022"] = +0.3;
                // Resonance masses
                p["mass::D_s@BSZ2015"]   = 1.968;
                p["mass::D_s,0@BSZ2015"] = 2.318;
                p["mass::D_s^*@BSZ2015"] = 2.112;
                p["mass::D_s,1@BSZ2015"] = 2.460;

                BMRvD2022FormFactors<LambdaCToLambda> ff(p, Options{ });

                // form factor values
                TEST_CHECK_RELATIVE_ERROR(ff.f_time_v (0.0), 24.45540515, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_v (0.0), 24.45540515, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_v (0.0), 12.42764544, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_time_a (0.0), 38.98340666, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_a (0.0), 38.98340666, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_a (0.0), 39.46347833, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_t (0.0), 15.08658937, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_t (0.0), 26.12119786, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_t5(0.0), 26.05081285, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_t5(0.0), 26.12119786, eps);

                TEST_CHECK_RELATIVE_ERROR(ff.f_time_v (3.0), 35.10118112, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_v (3.0), 54.25190515, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_v (3.0), 32.79010848, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_time_a (3.0), 129.3257177, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_a (3.0), 51.2073047, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_a (3.0), 49.87488124, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_t (3.0), 39.80568195, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_t (3.0), 51.72500778, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_t5(3.0), 34.29782468, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_t5(3.0), 34.31193557, eps);

                const auto tm = LambdaCToLambda::tm;
                TEST_CHECK_RELATIVE_ERROR(ff.f_time_v (tm), 27.62364063, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_v (tm), 30.97555432, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_v (tm), 16.84345485, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_time_a (tm), 53.59744051, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_a (tm), 42.88588404, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_a (tm), 42.88588404, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_t (tm), 20.44717867, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_t (tm), 31.66136446, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_long_t5(tm), 28.73608948, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.f_perp_t5(tm), 28.73608948, eps);
            }

        }
} bmrvd2022_form_factor_test;
