/*
 * Copyright (c) 2022 Méril Reboud
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
#include <eos/form-factors/parametric-abr2022.hh>
#include <eos/form-factors/parametric-abr2022-impl.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

using namespace test;
using namespace eos;

class ABR2022FormFactorsTest :
    public TestCase
{
    public:
        ABR2022FormFactorsTest() :
            TestCase("ABR2022_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();

                const auto ratio = (LambdaBToLambda1520::m1 + LambdaBToLambda1520::m2) / (LambdaBToLambda1520::m1 - LambdaBToLambda1520::m2);

                p["Lambda_b->Lambda(1520)::a^(t12,V)_1@ABR2022"]     =  0.1;
                p["Lambda_b->Lambda(1520)::a^(012,V)_1@ABR2022"]     = -0.2;
                p["Lambda_b->Lambda(1520)::a^(perp12,V)_1@ABR2022"]  =  0.3;
                p["Lambda_b->Lambda(1520)::a^(perp32,V)_0@ABR2022"]  = -0.4;
                p["Lambda_b->Lambda(1520)::a^(perp32,V)_1@ABR2022"]  =  0.5;
                p["Lambda_b->Lambda(1520)::a^(t12,A)_1@ABR2022"]     = -0.6;
                p["Lambda_b->Lambda(1520)::a^(012,A)_1@ABR2022"]     =  0.7;
                p["Lambda_b->Lambda(1520)::a^(perp12,A)_1@ABR2022"]  = -0.8;
                p["Lambda_b->Lambda(1520)::a^(perp32,A)_0@ABR2022"]  =  0.9;
                p["Lambda_b->Lambda(1520)::a^(perp32,A)_1@ABR2022"]  = -1.0;
                p["Lambda_b->Lambda(1520)::a^(012,T)_1@ABR2022"]     =  1.1;
                p["Lambda_b->Lambda(1520)::a^(perp12,T)_1@ABR2022"]  = -1.2;
                p["Lambda_b->Lambda(1520)::a^(perp32,T)_0@ABR2022"]  =  1.3;
                p["Lambda_b->Lambda(1520)::a^(perp32,T)_1@ABR2022"]  = -1.4;
                p["Lambda_b->Lambda(1520)::a^(012,T5)_1@ABR2022"]    =  1.5;
                p["Lambda_b->Lambda(1520)::a^(perp12,T5)_1@ABR2022"] = -1.6;
                p["Lambda_b->Lambda(1520)::a^(perp32,T5)_1@ABR2022"] =  1.7;

                ABR2022FormFactors<LambdaBToLambda1520> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                for (const auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.0875203,   eps), // z(q2 =  0)
                    std::make_pair( -0.00144262,  eps), // z(q2 = 10)
                    std::make_pair(  0.540328,    eps), // p_0(z = 0.0)
                    std::make_pair( -0.382764,    eps), // p_1(z = 0.0)
                    std::make_pair(  0.529468,    eps), // p_2(z = 0.0)
                    std::make_pair( -0.708331,    eps), // p_3(z = 0.0)
                    std::make_pair(  0.939069,    eps), // p_4(z = 0.0)
                    std::make_pair( -1.24208,     eps), // p_5(z = 0.0)
                    std::make_pair(  0.540328,    eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.383719,    eps), // p_1(z = z(q2 = 10))
                    std::make_pair(  0.530618,    eps), // p_2(z = z(q2 = 10))
                    std::make_pair( -0.710289,    eps), // p_3(z = z(q2 = 10))
                    std::make_pair(  0.942231,    eps), // p_4(z = z(q2 = 10))
                    std::make_pair( -1.24702,     eps), // p_5(z = z(q2 = 10))
                    std::make_pair(  0.0598709,   eps), // phi_time12_v(z = z(q2 = 1))
                    std::make_pair(  0.0174033,   eps), // phi_long12_v(z = z(q2 = 1))
                    std::make_pair(  0.0394989,   eps), // phi_perp12_v(z = z(q2 = 1))
                    std::make_pair(  0.0684142,   eps), // phi_perp32_v(z = z(q2 = 1))
                    std::make_pair(  0.0718986,   eps), // phi_time12_a(z = z(q2 = 1))
                    std::make_pair(  0.0142028,   eps), // phi_long12_a(z = z(q2 = 1))
                    std::make_pair(  0.0561359,   eps), // phi_perp12_a(z = z(q2 = 1))
                    std::make_pair(  0.0972303,   eps), // phi_perp32_a(z = z(q2 = 1))
                    std::make_pair(  0.0341431,   eps), // phi_long12_t(z = z(q2 = 1))
                    std::make_pair(  0.0300871,   eps), // phi_perp12_t(z = z(q2 = 1))
                    std::make_pair(  0.0521123,   eps), // phi_perp32_t(z = z(q2 = 1))
                    std::make_pair(  0.0487881,   eps), // phi_long12_t5(z = z(q2 = 1))
                    std::make_pair(  0.0246875,   eps), // phi_perp12_t5(z = z(q2 = 1))
                    std::make_pair(  0.04276,     eps), // phi_perp32_t5(z = z(q2 = 1))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                // Test end-point relations

                TEST_CHECK_NEARLY_EQUAL( ff.f_time12_v(0.0),    power_of<2>(ratio) * ff.f_long12_v(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_a(0.0),    power_of<2>(ratio) * ff.f_time12_a(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_t5(0.0),   power_of<2>(ratio) * ff.f_perp12_t(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp32_t5(0.0), - power_of<2>(ratio) * ff.f_perp32_t(0.0), eps);

                const auto tm = LambdaBToLambda1520::tm;

                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_v(tm),  - ff.f_perp32_v(tm),                       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_v(tm),    2.0 / ratio * ff.f_perp32_v(tm),         eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_time12_a(tm),    0.0,                                     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_a(tm),    ff.f_long12_a(tm) + ff.f_perp32_a(tm),   eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_t(tm),  - ff.f_perp32_t(tm),                       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_t(tm),    2.0 * ratio * ff.f_perp32_t(tm),         eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_t5(tm),   ff.f_perp12_t5(tm) - ff.f_perp32_t5(tm), eps);

                // Cross-check against Marzia's Mathematica notebook

                p["Lambda_b->Lambda(1520)::a^(t12,V)_1@ABR2022"]     =  0.1;
                p["Lambda_b->Lambda(1520)::a^(012,V)_1@ABR2022"]     =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp12,V)_1@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,V)_0@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,V)_1@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(t12,A)_1@ABR2022"]     =  0.1;
                p["Lambda_b->Lambda(1520)::a^(012,A)_1@ABR2022"]     =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp12,A)_1@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,A)_0@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,A)_1@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(012,T)_1@ABR2022"]     =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp12,T)_1@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,T)_0@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,T)_1@ABR2022"]  =  0.1;
                p["Lambda_b->Lambda(1520)::a^(012,T5)_1@ABR2022"]    =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp12,T5)_1@ABR2022"] =  0.1;
                p["Lambda_b->Lambda(1520)::a^(perp32,T5)_1@ABR2022"] =  0.1;

                TEST_CHECK_NEARLY_EQUAL( ff.f_time12_v(1.0),    5.42425   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_v(1.0),    1.73399   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_v(1.0),    0.281696  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp32_v(1.0),    0.641278  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_time12_a(1.0),    0.338578  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_a(1.0),    1.01076   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_a(1.0),    0.423661  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp32_a(1.0),    0.259616  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_t(1.0),    1.86808   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_t(1.0),    0.369817  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp32_t(1.0),    0.841883  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_long12_t5(1.0),   3.33155   , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp12_t5(1.0),   1.1898    , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.f_perp32_t5(1.0),   -2.60108  , eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0.266777  , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0527454 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.0166528 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.0564762 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0588784 , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.150838  , eps);
            }
        }
} abr2022_form_factor_test;
