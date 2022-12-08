/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
#include <eos/form-factors/parametric-bfw2010-impl.hh>

using namespace test;
using namespace eos;

class BToKBFW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKBFW2010FormFactorsTest() :
            TestCase("b_to_k_bfw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K::a^f+_0@BFW2010"]     =  0.01;
                p["B->K::a^f+_1@BFW2010"]     = -0.02;
                p["B->K::a^fT_0@BFW2010"]     =  0.03;
                p["B->K::a^fT_1@BFW2010"]     = -0.04;
                p["B->K::a^f0_1@BFW2010"]     =  0.05;

                BFW2010FormFactors<BToK, PToP> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                for (const auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.144553,  eps), // z(q2 =  0)
                    std::make_pair(  0.0562515, eps), // z(q2 = 10)
                    std::make_pair(  0.398942,  eps), // p_0(z = 0.0)
                    std::make_pair(  0,         eps), // p_1(z = 0.0)
                    std::make_pair(  0.398942,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair(  0.0224411, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0386502, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0369222, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0350620, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0870638, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0874954, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0878441, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.0958942, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.0895886, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.0830672, eps), // phi_f_t(z = z(q2 =  4))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations

                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);
            }
        }
} b_to_k_bfw2010_form_factors_test;

class BToKstarBFW2010FormFactorsTest :
    public TestCase
{
    public:
        BToKstarBFW2010FormFactorsTest() :
            TestCase("b_to_kstar_bfw2010_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::a^V_0@BFW2010"]     =  0.01;
                p["B->K^*::a^V_1@BFW2010"]     = -0.02;
                p["B->K^*::a^A0_0@BFW2010"]    =  0.03;
                p["B->K^*::a^A0_1@BFW2010"]    = -0.04;
                p["B->K^*::a^A1_0@BFW2010"]    =  0.05;
                p["B->K^*::a^A1_1@BFW2010"]    = -0.06;
                p["B->K^*::a^A12_1@BFW2010"]   =  0.07;
                p["B->K^*::a^T1_0@BFW2010"]    = -0.08;
                p["B->K^*::a^T1_1@BFW2010"]    =  0.09;
                p["B->K^*::a^T2_1@BFW2010"]    = -0.10;
                p["B->K^*::a^T23_0@BFW2010"]   =  0.11;
                p["B->K^*::a^T23_1@BFW2010"]   = -0.12;

                BFW2010FormFactors<BToKstar, PToV> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                for (const auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.0996732, eps), // z_a(q2 =  0)
                    std::make_pair(  0.1055,    eps), // z_v(q2 =  0)
                    std::make_pair(  0.015544,  eps), // z_a(q2 = 10)
                    std::make_pair(  0.016631,  eps), // z_v(q2 = 10)
                    std::make_pair(  0.470141,  eps), // p_0(z = 0.0)
                    std::make_pair( -0.170298,  eps), // p_1(z = 0.0)
                    std::make_pair(  0.203424,  eps), // p_2(z = 0.0)
                    std::make_pair( -0.233582,  eps), // p_3(z = 0.0)
                    std::make_pair(  0.263066,  eps), // p_4(z = 0.0)
                    std::make_pair( -0.293474,  eps), // p_5(z = 0.0)
                    std::make_pair(  0.470141,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.161982,  eps), // p_1(z = z(q2 = 10))
                    std::make_pair(  0.199363,  eps), // p_2(z = z(q2 = 10))
                    std::make_pair( -0.228172,  eps), // p_3(z = z(q2 = 10))
                    std::make_pair(  0.256243,  eps), // p_4(z = z(q2 = 10))
                    std::make_pair( -0.285044,  eps), // p_5(z = z(q2 = 10))

                    std::make_pair(  0.1235120, eps), // phi_v(z = z(q2 = -2.0))
                    std::make_pair(  0.1167340, eps), // phi_v(z = z(q2 =  1.0))
                    std::make_pair(  0.1097080, eps), // phi_v(z = z(q2 =  4.0))
                    std::make_pair(  0.1910140, eps), // phi_a_0(z = z(q2 = -2.0))
                    std::make_pair(  0.1806190, eps), // phi_a_0(z = z(q2 =  1.0))
                    std::make_pair(  0.1698540, eps), // phi_a_0(z = z(q2 =  4.0))
                    std::make_pair(  0.0750259, eps), // phi_a_1(z = z(q2 = -2.0))
                    std::make_pair(  0.0756559, eps), // phi_a_1(z = z(q2 =  1.0))
                    std::make_pair(  0.0762639, eps), // phi_a_1(z = z(q2 =  4.0))
                    std::make_pair(  0.0271387, eps), // phi_a_12(z = z(q2 = -2.0))
                    std::make_pair(  0.0279553, eps), // phi_a_12(z = z(q2 =  1.0))
                    std::make_pair(  0.0288281, eps), // phi_a_12(z = z(q2 =  4.0))
                    std::make_pair(  0.0796013, eps), // phi_t_1(z = z(q2 = -2.0))
                    std::make_pair(  0.0769276, eps), // phi_t_1(z = z(q2 =  1.0))
                    std::make_pair(  0.0740455, eps), // phi_t_1(z = z(q2 =  4.0))
                    std::make_pair(  0.0337405, eps), // phi_t_2(z = z(q2 = -2.0))
                    std::make_pair(  0.0347557, eps), // phi_t_2(z = z(q2 =  1.0))
                    std::make_pair(  0.0358409, eps), // phi_t_2(z = z(q2 =  4.0))
                    std::make_pair(  0.0323542, eps), // phi_t_23(z = z(q2 = -2.0))
                    std::make_pair(  0.0326259, eps), // phi_t_23(z = z(q2 =  1.0))
                    std::make_pair(  0.0328880, eps), // phi_t_23(z = z(q2 =  4.0))
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations

                const double factor = (BToKstar::m_B * BToKstar::m_B - BToKstar::m_V * BToKstar::m_V) / 8 / BToKstar::m_B / BToKstar::m_V;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factor * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0),           ff.t_2(0.0), eps);
            }
        }
} b_to_kstar_bfw2010_form_factors_test;
