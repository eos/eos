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
#include <eos/form-factors/parametric-se-impl-p-to-p.hh>

using namespace test;
using namespace eos;

class BToKSEFormFactorsTest :
    public TestCase
{
    public:
        BToKSEFormFactorsTest() :
            TestCase("b_to_k_se_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K::a^f+_0@SE"]     =  0.01;
                p["B->K::a^f+_1@SE"]     = -0.02;
                p["B->K::a^fT_0@SE"]     =  0.03;
                p["B->K::a^fT_1@SE"]     = -0.04;
                p["B->K::a^f0_1@SE"]     =  0.05;

                p["mass::B_d@BSZ2015"]        =  5.279;
                p["mass::K_d@BSZ2015"]        =  0.492;
                p["mass::B_s@BSZ2015"]        =  5.367;
                p["mass::B_s^*@BSZ2015"]      =  5.416;
                p["mass::B_s,0@BSZ2015"]      =  5.711;
                p["mass::B_s,1@BSZ2015"]      =  5.750;

                // Optimized t0 = (mB + mK) * (sqrt(mB) - sqrt(mK))^2
                p["B->K::t0@SE"]         =  14.703305673;


                SEFormFactors<BToK, PToP> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.144596,  eps), // z(q2 =  0)
                    std::make_pair(  0.0562957, eps), // z(q2 = 10)
                    std::make_pair(  0.398942,  eps), // p_0(z = 0.0)
                    std::make_pair(  0,         eps), // p_1(z = 0.0)
                    std::make_pair(  0.398942,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair(  0.0224588, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0386505, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0369224, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0350622, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0870644, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0874959, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0878445, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.0958949, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.0895891, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.0830676, eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  0.018232,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against Nico's implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( -1.0), 0.157128,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  1.0), 0.153515,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  4.0), 0.147786,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.0736626, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( -1.0), 0.148165,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  1.0), 0.163096,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  4.0), 0.190614,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 1.9403,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( -1.0), 0.206692,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  1.0), 0.227051,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  4.0), 0.2649,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 2.87453,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),       0.00283240578, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),       0,             eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),       0.0005,        eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),       0,             eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),       0.0025,        eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),      0,             eps);

                // Test everything for tp smaller than the scalar resonance B_s0
                p["B->K::tp@SE"]         =  30.261001;

                SEFormFactors<BToK, PToP> ff2(p, Options{ });

                diagnostics = ff2.diagnostics();
                static const std::vector<std::pair<double, double>> reference2
                {
                    std::make_pair(  0.164809,  eps), // z(q2 =  0)
                    std::make_pair(  0.0659398, eps), // z(q2 = 10)
                    std::make_pair(  0.465369,  eps), // p_0(z = 0.0)
                    std::make_pair( -0.1574336, eps), // p_1(z = 0.0)
                    std::make_pair(  0.465369,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.1250388, eps), // p_1(z = z(q2 = 10))

                    std::make_pair(  0.0409811, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.0391304, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.0371339, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.0962603, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.0968854, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.0974312, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.0970603, eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.0904346, eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.0835803, eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  0.0542382,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference2);

                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(0.0),           ff2.f_p(0.0), eps);

                TEST_CHECK_NEARLY_EQUAL( ff2.f_0( -1.0), 0.224009, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(  1.0), 0.218985, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(  4.0), 0.211081, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0( 25.0), 0.120399, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p( -1.0), 0.214148, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p(  1.0), 0.229426, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p(  4.0), 0.257016, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_p( 25.0), 1.59695,  eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t( -1.0), 0.251754, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t(  1.0), 0.2723,   eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t(  4.0), 0.309838, eps);
                TEST_CHECK_NEARLY_EQUAL( ff2.f_t( 25.0), 2.36706,  eps);
            }
        }
} b_to_k_se_form_factors_test;
