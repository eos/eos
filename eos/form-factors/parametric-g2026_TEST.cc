/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Nico Gubernari
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
#include <eos/form-factors/parametric-g2026-impl.hh>

#include <algorithm>

using namespace test;
using namespace eos;

class BToKG2026FormFactorsTest :
    public TestCase
{
    public:
        BToKG2026FormFactorsTest() :
            TestCase("b_to_k_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K::a^f+_0@G2026"]     =  0.01;
                p["B->K::a^f+_1@G2026"]     = -0.02;
                p["B->K::a^f0_1@G2026"]     =  0.05;
                p["B->K::a^fT_0@G2026"]     =  0.03;
                p["B->K::a^fT_1@G2026"]     = -0.04;

                p["mass::B_d@BSZ2015"]      =  5.279;
                p["mass::K_d@BSZ2015"]      =  0.494;

                p["mass::B_s,A^0[1]@G2026"] = 5.367;
                p["mass::B_s,V^0[1]@G2026"] = 5.711;
                p["mass::B_s,V^1[1]@G2026"] = 5.415;
                p["mass::B_s,A^1[1]@G2026"] = 5.829;

                // Optimized s0V = s_th * (1 - sqrt(1 - s_- / s_th))
                p["B->K::t0_v@G2026"]       =  14.682165;
                p["B->K::tp_v@G2026"]       =  30.272004;
                p["B->K::Q2@G2026"]         =  0.0;
                p["B->K::eta@G2026"]        =  2.0;

                p["B->K::tchi_V0@G2026"]    =  1.42e-2;
                p["B->K::tchi_V1@G2026"]    =  7.30e-4;
                p["B->K::tchi_T1@G2026"]    =  4.55e-4;

                G2026FormFactors<BToK, PToP> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B->K::t0_v@G2026"));
                TEST_CHECK(has_used_parameter("B->K::tp_v@G2026"));
                TEST_CHECK(has_used_parameter("B->K::eta@G2026"));

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.16439553,  eps), // z(q2 =  0)
                    std::make_pair(  0.052958941, eps), // z(q2 = 11)
                    std::make_pair( -0.059753422, eps), // z(q2 = 18)
                    std::make_pair(  0.18728718,  eps), // z(q2 = -3)

                    std::make_pair(  1.0,           eps), // z_0(z = 0)
                    std::make_pair(  0.16439553,    eps), // z_1(z = 0)
                    std::make_pair(  0.027025890,   eps), // z_2(z = 0)
                    std::make_pair(  0.0044429356,  eps), // z_3(z = 0)
                    std::make_pair(  0.00073039876, eps), // z_4(z = 0)
                    std::make_pair(  0.00012007429, eps), // z_5(z = 0)

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = 19))
                    std::make_pair( -0.080897132,    eps), // z_1(z = z(q2 = 19))
                    std::make_pair(  0.0065443460,   eps), // z_2(z = z(q2 = 19))
                    std::make_pair( -0.00052941882,  eps), // z_3(z = z(q2 = 19))
                    std::make_pair(  0.000042828465, eps), // z_4(z = z(q2 = 19))
                    std::make_pair( -0.0000034647,   eps), // z_5(z = z(q2 = 19))

                    std::make_pair(  1.0,           eps), // z_0(z = z(q2 = -3))
                    std::make_pair(  0.18728718,    eps), // z_1(z = z(q2 = -3))
                    std::make_pair(  0.035076488,   eps), // z_2(z = z(q2 = -3))
                    std::make_pair(  0.0065693766,  eps), // z_3(z = z(q2 = -3))
                    std::make_pair(  0.0012303600,  eps), // z_4(z = z(q2 = -3))
                    std::make_pair(  0.00023043066, eps), // z_5(z = z(q2 = -3))

                    std::make_pair(  0.09925017, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.09477547, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.08994813, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.06678005, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.06447112, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.06187060, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.2435218,  eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.2269179,  eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.2097406,  eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  -0.0017120938,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.055553910, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.11501982,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 0.97754178,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.15219810,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.083367774, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0),-0.52688301,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.084602925, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.16545090,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 1.5162935,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),     0.0025029313, eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),     0,            eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),     0,            eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),     0.0005,       eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),     0.0025,       eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),    0,            eps); // [G:2026A]: saturation_AT1

                // Test everything for different s0V and sV
                p["B->K::t0_v@G2026"]    = -4.0;
                p["B->K::tp_v@G2026"]    =  33.327529;

                G2026FormFactors<BToK, PToP> ff2(p, Options{ });

                diagnostics = ff2.diagnostics();
                static const std::vector<std::pair<double, double>> reference2
                {
                    std::make_pair( -0.028329254,  eps), // z(q2 =  0)
                    std::make_pair( -0.12777540,   eps), // z(q2 = 11)
                    std::make_pair( -0.21891875,   eps), // z(q2 = 18)
                    std::make_pair( -0.0067887129, eps), // z(q2 = -3)

                    std::make_pair(  1.0,            eps), // z_0(z = 0)
                    std::make_pair( -0.028329254,    eps), // z_1(z = 0)
                    std::make_pair(  0.00080254664,  eps), // z_2(z = 0)
                    std::make_pair( -0.000022735548, eps), // z_3(z = 0)
                    std::make_pair(  6.4408111e-7,   eps), // z_4(z = 0)
                    std::make_pair( -1.8246337e-8,   eps), // z_5(z = 0)

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = 19))
                    std::make_pair( -0.23491673,     eps), // z_1(z = z(q2 = 19))
                    std::make_pair(  0.055185872,    eps), // z_2(z = z(q2 = 19))
                    std::make_pair( -0.012964085,    eps), // z_3(z = z(q2 = 19))
                    std::make_pair(  0.0030454805,   eps), // z_4(z = z(q2 = 19))
                    std::make_pair( -0.00071543433,  eps), // z_5(z = z(q2 = 19))

                    std::make_pair(  1.0,            eps), // z_0
                    std::make_pair( -0.0067887129,   eps), // z_1
                    std::make_pair(  0.000046086623, eps), // z_2
                    std::make_pair( -3.1286885e-7,   eps), // z_3
                    std::make_pair(  2.1239768e-9,   eps), // z_4
                    std::make_pair( -1.4419069e-11,  eps), // z_5

                    std::make_pair(  0.09239401, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.08860722, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.08450451, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.2153878,  eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.2172881,  eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.2190775,  eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.2375216,  eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.2227716,  eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.2074462,  eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  0.040519167,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference2);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(0.0),           ff2.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.14861611, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.27516510, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 2.3634555,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.27005819, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.23639101, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.19320897, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.16757222, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.32309113, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 3.2127419,  eps);
            }
        }
} b_to_k_g2026_form_factors_test;

class BToPiG2026FormFactorsTest :
    public TestCase
{
    public:
        BToPiG2026FormFactorsTest() :
            TestCase("b_to_pi_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->pi::a^f+_0@G2026"]    =  0.01;
                p["B->pi::a^f+_1@G2026"]    = -0.02;
                p["B->pi::a^f0_1@G2026"]    =  0.05;
                p["B->pi::a^fT_0@G2026"]    =  0.03;
                p["B->pi::a^fT_1@G2026"]    = -0.04;

                p["mass::B_d@BSZ2015"]      =  5.279;
                p["mass::pi^0@BSZ2015"]     =  0.135;

                p["mass::B_d,V^0[1]@G2026"] = 10.0e+10; // very high value to effectively remove fake pole
                p["mass::B_d,V^1[1]@G2026"] = 5.325;

                p["B->pi::t0_v@G2026"]      =  20.170454;
                p["B->pi::tp_v@G2026"]      =  29.311396;
                p["B->pi::Q2@G2026"]        =  0.0;
                p["B->pi::eta@G2026"]       =  1.5;

                p["B->pi::tchi_V0@G2026"]   =  1.5e-2;
                p["B->pi::tchi_V1@G2026"]   =  5.7e-4;
                p["B->pi::tchi_T1@G2026"]   =  4.0e-4;

                G2026FormFactors<BToPi, PToP> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("mass::B_d,V^0[1]@G2026"));
                TEST_CHECK(has_used_parameter("mass::B_d,V^1[1]@G2026"));
                TEST_CHECK(has_used_parameter("B->pi::eta@G2026"));

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.034571698, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.10096914,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 2.3292231,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.10667532,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.076547304, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0),-0.046932887, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.078078319, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.17672171,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 3.8490245,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),     0.00250559, eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),     0,          eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),     0,          eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),     0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),     0.0025,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),    0,          eps); // [G:2026A]: saturation_AT1
            }
        }
} b_to_pi_g2026_form_factors_test;

class BsToKG2026FormFactorsTest :
    public TestCase
{
    public:
        BsToKG2026FormFactorsTest() :
            TestCase("b_s_to_k_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B_s->K::a^f+_0@G2026"]   =  0.01;
                p["B_s->K::a^f+_1@G2026"]   = -0.02;
                p["B_s->K::a^f0_1@G2026"]   =  0.05;
                p["B_s->K::a^fT_0@G2026"]   =  0.03;
                p["B_s->K::a^fT_1@G2026"]   = -0.04;

                p["mass::B_s@BSZ2015"]      =  5.367;
                p["mass::K_d@BSZ2015"]      =  0.494;

                p["mass::B_u,V^0[1]@G2026"] = 10.0e+10; // very high value to effectively remove fake pole
                p["mass::B_u,V^1[1]@G2026"] = 5.325;

                p["B_s->K::t0_v@G2026"]     =  15.264615;
                p["B_s->K::tp_v@G2026"]     =  29.311396;
                p["B_s->K::Q2@G2026"]       =  0.0;
                p["B_s->K::eta@G2026"]      =  1.0;

                p["B_s->K::tchi_V0@G2026"]  =  1.5e-2;
                p["B_s->K::tchi_V1@G2026"]  =  5.7e-4;
                p["B_s->K::tchi_T1@G2026"]  =  4.0e-4;

                G2026FormFactors<BsToK, PToP> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B_s->K::eta@G2026"));

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.063702474, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.13724289,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 1.3965556,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.15146182,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.11005836,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0),-0.014571330, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.11059611,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.22061292,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 2.4124228,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),     0.00264205, eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),     0,          eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),     0,          eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),     0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),     0.0025,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),    0,          eps); // [G:2026A]: saturation_AT1
            }
        }
} b_s_to_k_g2026_form_factors_test;

class BToDG2026FormFactorsTest :
    public TestCase
{
    public:
        BToDG2026FormFactorsTest() :
            TestCase("b_to_d_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->D::a^f+_0@G2026"]     =  0.01;
                p["B->D::a^f+_1@G2026"]     = -0.02;
                p["B->D::a^f0_1@G2026"]     =  0.05;
                p["B->D::a^fT_0@G2026"]     =  0.03;
                p["B->D::a^fT_1@G2026"]     = -0.04;

                p["mass::B_d@BSZ2015"]      =  5.279;
                p["mass::D_u@BSZ2015"]      =  1.865;

                p["mass::B_c,A^0[1]@G2026"] = 6.274;
                p["mass::B_c,A^0[2]@G2026"] = 6.871;
                p["mass::B_c,V^0[1]@G2026"] = 6.707;
                p["mass::B_c,V^1[1]@G2026"] = 6.328;
                p["mass::B_c,V^1[2]@G2026"] = 6.922;
                p["mass::B_c,A^1[1]@G2026"] = 6.739;

                p["B->D::t0_v@G2026"]       =  6.2048829;
                p["B->D::tp_v@G2026"]       =  41.075281;
                p["B->D::Q2@G2026"]         =  0.0;
                p["B->D::eta@G2026"]        =  2.0;

                p["B->D::tchi_V0@G2026"]    =  6.5e-3;
                p["B->D::tchi_V1@G2026"]    =  5.3e-4;
                p["B->D::tchi_T1@G2026"]    =  4.7e-4;

                G2026FormFactors<BToD, PToP> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B->D::eta@G2026"));

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.19686877, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.33253391, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 0.98858265, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.32941156, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.29342209, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.18895608, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.31540983, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.54169034, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 1.7384875,  eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),     0.0028941315, eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),     0,            eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),     0,            eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),     0.0005,       eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),     0.0025,       eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),    0,            eps); // [G:2026A]: saturation_AT1
            }
        }
} b_to_d_g2026_form_factors_test;

class BsToDsG2026FormFactorsTest :
    public TestCase
{
    public:
        BsToDsG2026FormFactorsTest() :
            TestCase("b_s_to_d_s_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B_s->D_s::a^f+_0@G2026"]   =  0.01;
                p["B_s->D_s::a^f+_1@G2026"]   = -0.02;
                p["B_s->D_s::a^f0_1@G2026"]   =  0.05;
                p["B_s->D_s::a^fT_0@G2026"]   =  0.03;
                p["B_s->D_s::a^fT_1@G2026"]   = -0.04;

                p["mass::B_s@BSZ2015"]        =  5.367;
                p["mass::D_s@BSZ2015"]        =  1.968;

                p["mass::B_c,A^0[1]@G2026"]   = 6.274;
                p["mass::B_c,A^0[2]@G2026"]   = 6.871;
                p["mass::B_c,V^0[1]@G2026"]   = 6.707;
                p["mass::B_c,V^1[1]@G2026"]   = 6.328;
                p["mass::B_c,V^1[2]@G2026"]   = 6.922;
                p["mass::B_c,A^1[1]@G2026"]   = 6.739;

                p["B_s->D_s::t0_v@G2026"]     =  6.1252757;
                p["B_s->D_s::tp_v@G2026"]     =  41.075281;
                p["B_s->D_s::Q2@G2026"]       =  0.0;
                p["B_s->D_s::eta@G2026"]      =  1.0;

                p["B_s->D_s::tchi_V0@G2026"]  =  6.5e-3;
                p["B_s->D_s::tchi_V1@G2026"]  =  5.3e-4;
                p["B_s->D_s::tchi_T1@G2026"]  =  4.7e-4;

                G2026FormFactors<BsToDs, PToP> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B_s->D_s::eta@G2026"));

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.27009881, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.45093416, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 1.2950788,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.44943919, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.39848725, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.25260665, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.44404500, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.75389239, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 2.3378782,  eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),     0.0028840563, eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),     0,            eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),     0,            eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),     0.0005,       eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),     0.0025,       eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(),    0,            eps); // [G:2026A]: saturation_AT1
            }
        }
} b_s_to_d_s_g2026_form_factors_test;

class BsToKstarG2026FormFactorsTest :
    public TestCase
{
    public:
        BsToKstarG2026FormFactorsTest() :
            TestCase("b_s_to_kstar_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B_s->K^*::a^V_0@G2026"]     =  0.01;
                p["B_s->K^*::a^V_1@G2026"]     = -0.02;
                p["B_s->K^*::a^A0_0@G2026"]    =  0.03;
                p["B_s->K^*::a^A0_1@G2026"]    = -0.04;
                p["B_s->K^*::a^A1_1@G2026"]    =  0.05;
                p["B_s->K^*::a^A12_1@G2026"]   = -0.06;
                p["B_s->K^*::a^T1_0@G2026"]    =  0.07;
                p["B_s->K^*::a^T1_1@G2026"]    = -0.08;
                p["B_s->K^*::a^T2_1@G2026"]    =  0.09;
                p["B_s->K^*::a^T23_1@G2026"]   = -0.10;

                p["mass::B_s@BSZ2015"]         =  5.367;
                p["mass::K_d^*@BSZ2015"]       =  0.892;

                p["mass::B_d,A^0[1]@G2026"]    = 5.279;
                p["mass::B_d,V^1[1]@G2026"]    = 5.325;
                p["mass::B_d,A^1[1]@G2026"]    = 5.726;

                p["B_s->K^*::t0_v@G2026"]      =  11.785640;
                p["B_s->K^*::t0_a@G2026"]      =  11.785640;
                p["B_s->K^*::tp_v@G2026"]      =  29.311396;
                p["B_s->K^*::tp_a@G2026"]      =  29.811600;
                p["B_s->K^*::Q2@G2026"]        =  0.0;
                p["B_s->K^*::eta@G2026"]       =  1.0;

                p["B_s->K^*::tchi_A0@G2026"]   =  1.5e-2;
                p["B_s->K^*::tchi_V1@G2026"]   =  5.7e-4;
                p["B_s->K^*::tchi_A1@G2026"]   =  5.7e-4;
                p["B_s->K^*::tchi_T1@G2026"]   =  4.0e-4;
                p["B_s->K^*::tchi_AT1@G2026"]  =  4.0e-4;

                G2026FormFactors<BsToKstar, PToV> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("mass::B_d,A^0[1]@G2026"));
                TEST_CHECK(has_used_parameter("mass::B_d,V^1[1]@G2026"));
                TEST_CHECK(has_used_parameter("mass::B_d,A^1[1]@G2026"));
                TEST_CHECK(has_used_parameter("B_s->K^*::t0_v@G2026"));
                TEST_CHECK(has_used_parameter("B_s->K^*::t0_a@G2026"));
                TEST_CHECK(has_used_parameter("B_s->K^*::tp_v@G2026"));
                TEST_CHECK(has_used_parameter("B_s->K^*::tp_a@G2026"));
                TEST_CHECK(has_used_parameter("B_s->K^*::eta@G2026"));

                // Test end-point relations
                const double mB = p["mass::B_s@BSZ2015"];
                const double mV = p["mass::K_d^*@BSZ2015"];
                const double sm = (mB - mV) * (mB - mV);

                const double factora12a0 = (mB * mB - mV * mV) / 8.0 / mB / mV;
                const double factora12a1 = (mB + mV) * (mB * mB - mV * mV - sm) / 16.0 / mB / mV / mV;
                const double factort23t2 = (mB + mV) * (mB * mB + 3.0 * mV * mV - sm) / 8.0 / mB / mV / mV;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(sm) , factora12a1 * ff.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(sm) , factort23t2 * ff.t_2(sm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.024422627, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.052891327, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.50851848,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.072078737, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.13688178,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.3770939,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  1.3058481,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  1.3587422,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  1.6975007,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.23023283,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.17281056,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  2.0631331,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  0.37526998,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  0.57930495,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  2.9529116,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  0.91725585,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  0.42910793,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0), -1.6406239,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0), -2.5396854,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0), -2.2455291,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0), -0.052814319, eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0,          eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.0025,     eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.0094409,  eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0113,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.0191505,  eps); // [G:2026A]: saturation_AT1
            }
        }
} b_s_to_kstar_g2026_form_factors_test;

class BToKstarG2026FormFactorsTest :
    public TestCase
{
    public:
        BToKstarG2026FormFactorsTest() :
            TestCase("b_to_kstar_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::a^V_0@G2026"]    =  0.01;
                p["B->K^*::a^V_1@G2026"]    = -0.02;
                p["B->K^*::a^A0_0@G2026"]   =  0.03;
                p["B->K^*::a^A0_1@G2026"]   = -0.04;
                p["B->K^*::a^A1_1@G2026"]   =  0.05;
                p["B->K^*::a^A12_1@G2026"]  = -0.06;
                p["B->K^*::a^T1_0@G2026"]   =  0.07;
                p["B->K^*::a^T1_1@G2026"]   = -0.08;
                p["B->K^*::a^T2_1@G2026"]   =  0.09;
                p["B->K^*::a^T23_1@G2026"]  = -0.10;

                p["mass::B_d@BSZ2015"]      =  5.279;
                p["mass::K_d^*@BSZ2015"]    =  0.892;

                p["mass::B_s,A^0[1]@G2026"] = 5.367;
                p["mass::B_s,V^0[1]@G2026"] = 5.711;
                p["mass::B_s,V^1[1]@G2026"] = 5.415;
                p["mass::B_s,A^1[1]@G2026"] = 5.829;

                // Optimized s0A = s0V = s_th * (1 - sqrt(1 - s_- / s_th))
                p["B->K^*::t0_v@G2026"]     =  11.299192;
                p["B->K^*::t0_a@G2026"]     =  11.299192;
                p["B->K^*::tp_v@G2026"]     =  30.272004;
                p["B->K^*::tp_a@G2026"]     =  31.775769;
                p["B->K^*::Q2@G2026"]       =  0.0;
                p["B->K^*::eta@G2026"]      =  2.0;

                p["B->K^*::tchi_A0@G2026"]  =  1.52e-2;
                p["B->K^*::tchi_V1@G2026"]  =  6.97e-4;
                p["B->K^*::tchi_A1@G2026"]  =  6.54e-4;
                p["B->K^*::tchi_T1@G2026"]  =  4.83e-4;
                p["B->K^*::tchi_AT1@G2026"] =  4.12e-4;

                G2026FormFactors<BToKstar, PToV> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B->K^*::eta@G2026"));

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.11627577,   eps), // z(q2 =  0, sV)
                    std::make_pair(  0.0039115965, eps), // z(q2 = 11, sV)
                    std::make_pair( -0.098769544,  eps), // z(q2 = 18, sA)
                    std::make_pair(  0.13164137,   eps), // z(q2 = -3, sA)

                    std::make_pair(  1.0,            eps), // z_0(z = 0, sV)
                    std::make_pair(  0.11627577,     eps), // z_1(z = 0, sV)
                    std::make_pair(  0.013520054,    eps), // z_2(z = 0, sV)
                    std::make_pair(  0.0015720547,   eps), // z_3(z = 0, sV)
                    std::make_pair(  0.00018279187,  eps), // z_4(z = 0, sV)
                    std::make_pair(  0.000021254265, eps), // z_5(z = 0, sV)

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = 19, sV))
                    std::make_pair( -0.12944094,     eps), // z_1(z = z(q2 = 19, sV))
                    std::make_pair(  0.016754956,    eps), // z_2(z = z(q2 = 19, sV))
                    std::make_pair( -0.0021687773,   eps), // z_3(z = z(q2 = 19, sV))
                    std::make_pair(  0.00028072856,  eps), // z_4(z = z(q2 = 19, sV))
                    std::make_pair( -0.000036337769, eps), // z_5(z = z(q2 = 19, sV))

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = -3, sA))
                    std::make_pair(  0.13164137,     eps), // z_1(z = z(q2 = -3, sA))
                    std::make_pair(  0.017329451,    eps), // z_2(z = z(q2 = -3, sA))
                    std::make_pair(  0.0022812726,   eps), // z_3(z = z(q2 = -3, sA))
                    std::make_pair(  0.00030030986,  eps), // z_4(z = z(q2 = -3, sA))
                    std::make_pair(  0.000039533201, eps), // z_5(z = z(q2 = -3, sA))

                    std::make_pair(  0.3152006,  eps), // phi_v(z = z(q2 = -2))
                    std::make_pair(  0.2974732,  eps), // phi_v(z = z(q2 =  1))
                    std::make_pair(  0.2790677,  eps), // phi_v(z = z(q2 =  4))

                    std::make_pair(  0.5025545,  eps), // phi_a_0(z = z(q2 = -2))
                    std::make_pair(  0.4746616,  eps), // phi_a_0(z = z(q2 =  1))
                    std::make_pair(  0.4457290,  eps), // phi_a_0(z = z(q2 =  4))
                    std::make_pair(  0.05659679, eps), // phi_a_1(z = z(q2 = -2))
                    std::make_pair(  0.05492700, eps), // phi_a_1(z = z(q2 =  1))
                    std::make_pair(  0.05304411, eps), // phi_a_1(z = z(q2 =  4))
                    std::make_pair(  0.02133889, eps), // phi_a_12(z = z(q2 = -2))
                    std::make_pair(  0.02119834, eps), // phi_a_12(z = z(q2 =  1))
                    std::make_pair(  0.02099215, eps), // phi_a_12(z = z(q2 =  4))

                    std::make_pair(  0.2089456, eps), // phi_t_1(z = z(q2 = -2))
                    std::make_pair(  0.2020819, eps), // phi_t_1(z = z(q2 =  1))
                    std::make_pair(  0.1946578, eps), // phi_t_1(z = z(q2 =  4))

                    std::make_pair(  0.02732397, eps), // phi_t_2(z = z(q2 = -2))
                    std::make_pair(  0.02714399, eps), // phi_t_2(z = z(q2 =  1))
                    std::make_pair(  0.02687998, eps), // phi_t_2(z = z(q2 =  4))
                    std::make_pair(  0.02493920, eps), // phi_t_23(z = z(q2 = -2))
                    std::make_pair(  0.02420341, eps), // phi_t_23(z = z(q2 =  1))
                    std::make_pair(  0.02337372, eps), // phi_t_23(z = z(q2 =  4))

                    std::make_pair(  0.00808554091, eps), // a_A12_0
                    std::make_pair(  0.05171365,    eps), // a_A1_0
                    std::make_pair(  0.001710846,   eps), // a_T2_0
                    std::make_pair(  -0.02595125,   eps), // a_T23_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                const double mB = p["mass::B_d@BSZ2015"];
                const double mV = p["mass::K_d^*@BSZ2015"];
                const double sm = (mB - mV) * (mB - mV);

                const double factora12a0 = (mB * mB - mV * mV) / 8.0 / mB / mV;
                const double factora12a1 = (mB + mV) * (mB * mB - mV * mV - sm) / 16.0 / mB / mV / mV;
                const double factort23t2 = (mB + mV) * (mB * mB + 3.0 * mV * mV - sm) / 8.0 / mB / mV / mV;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(sm) , factora12a1 * ff.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(sm) , factort23t2 * ff.t_2(sm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.019856051, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.041885252, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.33524673,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.059631278, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.11372013,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.0100095,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  0.99821129,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  1.0420689,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  1.3310723,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.19263127,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.14210536,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  1.6126605,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  0.30288816,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  0.46431383,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  2.0639932,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  0.72944768,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  0.34650161,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0), -1.1708205,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0), -1.6950877,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0), -1.4553959,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  0.080068039, eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0,          eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.0025,     eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.00883968, eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0113,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.0187764,  eps); // [G:2026A]: saturation_AT1


                // Test everything for different s0A = s0V and sA = sV
                p["B->K^*::t0_v@G2026"] = -4.0;
                p["B->K^*::t0_a@G2026"] = -4.0;
                p["B->K^*::tp_v@G2026"] =  38.081241;
                p["B->K^*::tp_a@G2026"] =  38.081241;

                // Test end-point relations

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(sm) , factora12a1 * ff.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(sm) , factort23t2 * ff.t_2(sm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.060354646, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.12451300,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  1.2368289,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.12014397,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.23335484,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  2.3503654,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  0.69051914,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  0.69204300,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  0.93515203,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.021451283, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.18777326,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  0.92333822,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  0.81003891,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  1.3711616,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  9.6563017,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  1.4521739,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  1.1807589,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0),  0.85018484,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0),  0.67859760,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0),  1.0230491,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  2.8515852,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0,          eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.0025,     eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.00969611, eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0113,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.0215922,  eps); // [G:2026A]: saturation_AT1
            }
        }
} b_to_kstar_g2026_form_factors_test;

class BsToPhiG2026FormFactorsTest :
    public TestCase
{
    public:
        BsToPhiG2026FormFactorsTest() :
            TestCase("b_s_to_phi_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B_s->phi::a^V_0@G2026"]     =  0.01;
                p["B_s->phi::a^V_1@G2026"]     = -0.02;
                p["B_s->phi::a^V_2@G2026"]     = -0.03;
                p["B_s->phi::a^A0_0@G2026"]    =  0.03;
                p["B_s->phi::a^A0_1@G2026"]    = -0.04;
                p["B_s->phi::a^A0_2@G2026"]    = -0.05;
                p["B_s->phi::a^A1_1@G2026"]    =  0.05;
                p["B_s->phi::a^A1_2@G2026"]    =  0.06;
                p["B_s->phi::a^A12_1@G2026"]   = -0.06;
                p["B_s->phi::a^A12_2@G2026"]   = -0.07;
                p["B_s->phi::a^T1_0@G2026"]    =  0.07;
                p["B_s->phi::a^T1_1@G2026"]    = -0.08;
                p["B_s->phi::a^T1_2@G2026"]    = -0.09;
                p["B_s->phi::a^T2_1@G2026"]    =  0.09;
                p["B_s->phi::a^T2_2@G2026"]    =  0.10;
                p["B_s->phi::a^T23_1@G2026"]   = -0.10;
                p["B_s->phi::a^T23_2@G2026"]   = -0.11;

                p["mass::B_s@BSZ2015"]         =  5.367;
                p["mass::phi@BSZ2015"]         =  1.020;

                p["mass::B_s,A^0[1]@G2026"]    = 5.367;
                p["mass::B_s,V^0[1]@G2026"]    = 5.711;
                p["mass::B_s,V^1[1]@G2026"]    = 5.415;
                p["mass::B_s,A^1[1]@G2026"]    = 5.829;

                // Optimized s0A = s0V = s_th * (1 - sqrt(1 - s_- / s_th))
                p["B_s->phi::t0_v@G2026"]      =  10.906048;
                p["B_s->phi::t0_a@G2026"]      =  10.906048;
                p["B_s->phi::tp_v@G2026"]      =  30.272004;
                p["B_s->phi::tp_a@G2026"]      =  31.775769;
                p["B_s->phi::Q2@G2026"]        =  0.0;
                p["B_s->phi::eta@G2026"]       =  1.0;

                p["B_s->phi::tchi_A0@G2026"]   =  1.52e-2;
                p["B_s->phi::tchi_V1@G2026"]   =  6.97e-4;
                p["B_s->phi::tchi_A1@G2026"]   =  6.54e-4;
                p["B_s->phi::tchi_T1@G2026"]   =  4.83e-4;
                p["B_s->phi::tchi_AT1@G2026"]  =  4.12e-4;

                G2026FormFactors<BsToPhi, PToV> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B_s->phi::eta@G2026"));

                // Test end-point relations
                const double mB = p["mass::B_s@BSZ2015"];
                const double mV = p["mass::phi@BSZ2015"];
                const double sm = (mB - mV) * (mB - mV);

                const double factora12a0 = (mB * mB - mV * mV) / 8.0 / mB / mV;
                const double factora12a1 = (mB + mV) * (mB * mB - mV * mV - sm) / 16.0 / mB / mV / mV;
                const double factort23t2 = (mB + mV) * (mB * mB + 3.0 * mV * mV - sm) / 8.0 / mB / mV / mV;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(sm) , factora12a1 * ff.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(sm) , factort23t2 * ff.t_2(sm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.021904674, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.056655476, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.34513055,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.074096030, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.15014971,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.1164464,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  1.3401280,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  1.3386879,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  1.8119204,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.34618089,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.18490812,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  1.5125224,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  0.38239900,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  0.61315297,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  2.2710191,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  1.1193157,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  0.42963852,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0), -1.1156022,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0), -2.1654880,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0), -1.6727656,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0), -0.40022028,  eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0,         eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.005,     eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.0171064, eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0014,    eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0194,    eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.0407473, eps); // [G:2026A]: saturation_AT1
            }
        }
} b_s_to_phi_g2026_form_factors_test;

class BToDstarG2026FormFactorsTest :
    public TestCase
{
    public:
        BToDstarG2026FormFactorsTest() :
            TestCase("b_to_dstar_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->D^*::a^V_0@G2026"]     =  0.01;
                p["B->D^*::a^V_1@G2026"]     = -0.02;
                p["B->D^*::a^A0_0@G2026"]    =  0.03;
                p["B->D^*::a^A0_1@G2026"]    = -0.04;
                p["B->D^*::a^A1_1@G2026"]    =  0.05;
                p["B->D^*::a^A12_1@G2026"]   = -0.06;
                p["B->D^*::a^T1_0@G2026"]    =  0.07;
                p["B->D^*::a^T1_1@G2026"]    = -0.08;
                p["B->D^*::a^T2_1@G2026"]    =  0.09;
                p["B->D^*::a^T23_1@G2026"]   = -0.10;

                p["mass::B_d@BSZ2015"]       =  5.279;
                p["mass::D_u^*@BSZ2015"]     =  2.007;

                p["mass::B_c,A^0[1]@G2026"]  = 6.274;
                p["mass::B_c,A^0[2]@G2026"]  = 6.871;
                p["mass::B_c,V^0[1]@G2026"]  = 6.707;
                p["mass::B_c,V^1[1]@G2026"]  = 6.328;
                p["mass::B_c,V^1[2]@G2026"]  = 6.922;
                p["mass::B_c,A^1[1]@G2026"]  = 6.739;

                p["B->D^*::t0_v@G2026"]      =  5.6540972;
                p["B->D^*::t0_a@G2026"]      =  5.6540972;
                p["B->D^*::tp_v@G2026"]      =  41.075281;
                p["B->D^*::tp_a@G2026"]      =  41.770369;
                p["B->D^*::Q2@G2026"]        =  0.0;
                p["B->D^*::eta@G2026"]       =  2.0;

                p["B->D^*::tchi_A0@G2026"]   =  1.7e-2;
                p["B->D^*::tchi_V1@G2026"]   =  5.3e-4;
                p["B->D^*::tchi_A1@G2026"]   =  3.9e-4;
                p["B->D^*::tchi_T1@G2026"]   =  4.7e-4;
                p["B->D^*::tchi_AT1@G2026"]  =  2.2e-4;

                G2026FormFactors<BToDstar, PToV> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B->D^*::eta@G2026"));
                // Test end-point relations
                const double mB = p["mass::B_d@BSZ2015"];
                const double mV = p["mass::D_u^*@BSZ2015"];
                const double sm = (mB - mV) * (mB - mV);

                const double factora12a0 = (mB * mB - mV * mV) / 8.0 / mB / mV;
                const double factora12a1 = (mB + mV) * (mB * mB - mV * mV - sm) / 16.0 / mB / mV / mV;
                const double factort23t2 = (mB + mV) * (mB * mB + 3.0 * mV * mV - sm) / 8.0 / mB / mV / mV;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(sm) , factora12a1 * ff.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(sm) , factort23t2 * ff.t_2(sm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.071362907, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.13160062,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.45769644,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.17937997,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.31649092,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.1115797,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  0.64762210,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  0.64128635,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  0.64672796,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.021963829, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.10456170,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  0.41609745,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  1.0112868,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  1.5218151,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  3.7954573,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  1.7051950,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  1.3357433,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0),  0.60563906,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0),  0.95933568,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0),  1.2327428,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  2.1412154,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0,          eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.0025,     eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.00804022, eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0113,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.0233911,  eps); // [G:2026A]: saturation_AT1
            }
        }
} b_to_dstar_g2026_form_factors_test;

class BsToDsstarG2026FormFactorsTest :
    public TestCase
{
    public:
        BsToDsstarG2026FormFactorsTest() :
            TestCase("b_s_to_d_sstar_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B_s->D_s^*::a^V_0@G2026"]     =  0.01;
                p["B_s->D_s^*::a^V_1@G2026"]     = -0.02;
                p["B_s->D_s^*::a^A0_0@G2026"]    =  0.03;
                p["B_s->D_s^*::a^A0_1@G2026"]    = -0.04;
                p["B_s->D_s^*::a^A1_1@G2026"]    =  0.05;
                p["B_s->D_s^*::a^A12_1@G2026"]   = -0.06;
                p["B_s->D_s^*::a^T1_0@G2026"]    =  0.07;
                p["B_s->D_s^*::a^T1_1@G2026"]    = -0.08;
                p["B_s->D_s^*::a^T2_1@G2026"]    =  0.09;
                p["B_s->D_s^*::a^T23_1@G2026"]   = -0.10;

                p["mass::B_s@BSZ2015"]           =  5.367;
                p["mass::D_s^*@BSZ2015"]         =  2.112;

                p["mass::B_c,A^0[1]@G2026"]      = 6.274;
                p["mass::B_c,A^0[2]@G2026"]      = 6.871;
                p["mass::B_c,V^0[1]@G2026"]      = 6.707;
                p["mass::B_c,V^1[1]@G2026"]      = 6.328;
                p["mass::B_c,V^1[2]@G2026"]      = 6.922;
                p["mass::B_c,A^1[1]@G2026"]      = 6.739;

                p["B_s->D_s^*::t0_v@G2026"]      =  5.5753758;
                p["B_s->D_s^*::t0_a@G2026"]      =  5.5753758;
                p["B_s->D_s^*::tp_v@G2026"]      =  41.075281;
                p["B_s->D_s^*::tp_a@G2026"]      =  41.770369;
                p["B_s->D_s^*::Q2@G2026"]        =  0.0;
                p["B_s->D_s^*::eta@G2026"]       =  1.0;

                p["B_s->D_s^*::tchi_A0@G2026"]   =  1.7e-2;
                p["B_s->D_s^*::tchi_V1@G2026"]   =  5.3e-4;
                p["B_s->D_s^*::tchi_A1@G2026"]   =  3.9e-4;
                p["B_s->D_s^*::tchi_T1@G2026"]   =  4.7e-4;
                p["B_s->D_s^*::tchi_AT1@G2026"]  =  2.2e-4;

                G2026FormFactors<BsToDsstar, PToV> ff(p, Options{ });

                const auto has_used_parameter = [&](const std::string & name)
                {
                    return std::find(ff.begin(), ff.end(), p[name].id()) != ff.end();
                };

                TEST_CHECK(has_used_parameter("B_s->D_s^*::eta@G2026"));

                // Test end-point relations
                const double mB = p["mass::B_s@BSZ2015"];
                const double mV = p["mass::D_s^*@BSZ2015"];
                const double sm = (mB - mV) * (mB - mV);

                const double factora12a0 = (mB * mB - mV * mV) / 8.0 / mB / mV;
                const double factora12a1 = (mB + mV) * (mB * mB - mV * mV - sm) / 16.0 / mB / mV / mV;
                const double factort23t2 = (mB + mV) * (mB * mB + 3.0 * mV * mV - sm) / 8.0 / mB / mV / mV;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(sm) , factora12a1 * ff.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(sm) , factort23t2 * ff.t_2(sm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.10049093, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.18325111, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.61740229, eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.24594183,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.42917085,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.4605161,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  0.87982219,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  0.86769133,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  0.86490672,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.029291323, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.13571614,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  0.53593903,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  1.3863425,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  2.0633724,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  4.9863597,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  2.3263090,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  1.8133843,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0),  0.80715380,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0),  1.2938586,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0),  1.6554173,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  2.8425204,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0p_v(),  0,          eps); // [G:2026A]: saturation_V0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_0m_a(),  0.0025,     eps); // [G:2026A]: saturation_A0
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_a(),  0.00802717, eps); // [G:2026A]: saturation_A1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_v(),  0.0005,     eps); // [G:2026A]: saturation_V1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1m_t(),  0.0113,     eps); // [G:2026A]: saturation_T1
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_1p_t5(), 0.0234089,  eps); // [G:2026A]: saturation_AT1

                // Test end-point relations for different s0A != s0V
                Parameters p_endpoint = p;
                p_endpoint["B_s->D_s^*::t0_v@G2026"] = 5.6918798;
                p_endpoint["B_s->D_s^*::t0_a@G2026"] = 5.6842827;

                G2026FormFactors<BsToDsstar, PToV> ff_endpoint(p_endpoint, Options{ });

                TEST_CHECK_NEARLY_EQUAL( ff_endpoint.a_12(0.0), factora12a0 * ff_endpoint.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff_endpoint.t_1(0.0) ,               ff_endpoint.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff_endpoint.a_12(sm) , factora12a1 * ff_endpoint.a_1(sm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff_endpoint.t_23(sm) , factort23t2 * ff_endpoint.t_2(sm) , eps);
            }
        }
} b_s_to_d_sstar_g2026_form_factors_test;
