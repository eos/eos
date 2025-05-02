/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2025 Danny van Dyk
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
#include <eos/form-factors/parametric-bgjvd2019.hh>

#include <vector>

using namespace test;
using namespace eos;

class BToDHQETFormFactorsTest :
    public TestCase
{
    public:
        BToDHQETFormFactorsTest() :
            TestCase("b_to_d_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            // using z_* with a = 1.0 and LP z-order = 2 and SLP z-order 2 and SSLP z-order 1
            // Martin's best-fit point
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -0.849472;
                p["B(*)->D(*)::xi''(1)@HQET"]    =  2.0 * 0.583711;
                p["B(*)->D(*)::xi'''(1)@HQET"]   =  0.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  =  0.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = -0.0600533;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  =  6.97061e-6;
                p["B(*)->D(*)::chi_2''(1)@HQET"] =  0.0314499;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  =  0.0400298;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = -0.039123;
                p["B(*)->D(*)::eta(1)@HQET"]     =  0.604052;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -0.00545745;
                p["B(*)->D(*)::eta''(1)@HQET"]   = -0.268764;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.111274;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.01963;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0687349;
                p["B(*)->D(*)::l_4(1)@HQET"]     = -2.02231;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  4.21978;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  4.52949;
                p["B(*)->D(*)::l_1'(1)@HQET"]    = -15.0241;
                p["B(*)->D(*)::l_2'(1)@HQET"]    = -9.43754;
                p["B(*)->D(*)::l_3'(1)@HQET"]    = -0.616533;
                p["B(*)->D(*)::l_4'(1)@HQET"]    = +0.604533;
                p["B(*)->D(*)::l_5'(1)@HQET"]    = +0.115125;
                p["B(*)->D(*)::l_6'(1)@HQET"]    = -1.4777;
                p["B(*)->D(*)::a@HQET"]          =  1.0;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp"_ok,   "2" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToD, PToP> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1, LP z^3, SLP z^1" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+0.0, eps), // LP z^3 terms enabled?
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps), // SLP z^2 terms enabled?

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+0.541418, eps), // w = 2.10
                    std::make_pair(+0.656849, eps), // w = 1.60
                    std::make_pair(+0.920648, eps), // w = 1.10
                    std::make_pair(+0.958955, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(-0.0480609, eps), // w = 2.10
                    std::make_pair(-0.0557318, eps), // w = 1.60
                    std::make_pair(-0.0599029, eps), // w = 1.10
                    std::make_pair(-0.0600146, eps), // w = 1.05
                    std::make_pair(-0.0600533, eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(+0.027665,  eps), // w = 2.10
                    std::make_pair(+0.0183516, eps), // w = 1.60
                    std::make_pair(+0.00381496,eps), // w = 1.10
                    std::make_pair(+0.00195355,eps), // w = 1.05
                    std::make_pair( 0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(+0.495827, eps), // w = 2.10
                    std::make_pair(+0.563923, eps), // w = 1.60
                    std::make_pair(+0.602227, eps), // w = 1.10
                    std::make_pair(+0.603451, eps), // w = 1.05
                    std::make_pair(+0.604052, eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(+0.669971, eps), // h_{p}
                    std::make_pair(-0.043089, eps), // h_{m}
                    std::make_pair(+0.775406, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.820365, eps), // h_{p}
                    std::make_pair(-0.049890, eps), // h_{m}
                    std::make_pair(+0.936907, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.024402, eps), // h_{p}
                    std::make_pair(-0.061313, eps), // h_{m}
                    std::make_pair(+1.159468, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);
            }

            // using z_* with a = 1.0 and LP z-order = 3 and SLP z-order 1
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp"_ok,   "3" },
                    { "z-order-slp"_ok,  "1" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToD, PToP> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1, LP z^3, SLP z^1" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+1.665540, eps), // w = 2.10
                    std::make_pair(+0.764544, eps), // w = 1.60
                    std::make_pair(+0.865908, eps), // w = 1.10
                    std::make_pair(+0.928869, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(-0.373019,  eps), // w = 2.10
                    std::make_pair(-0.0239773, eps), // w = 1.60
                    std::make_pair(+0.402425,  eps), // w = 1.10
                    std::make_pair(+0.450615,  eps), // w = 1.05
                    std::make_pair(+0.5,       eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-1.30953,   eps), // w = 2.10
                    std::make_pair(-0.785966,  eps), // w = 1.60
                    std::make_pair(-0.146363,  eps), // w = 1.10
                    std::make_pair(-0.0740769, eps), // w = 1.05
                    std::make_pair( 0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.841274, eps), // w = 2.10
                    std::make_pair(-0.404972, eps), // w = 1.60
                    std::make_pair(+0.128031, eps), // w = 1.10
                    std::make_pair(+0.188269, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(-0.371631, eps), // h_{p}
                    std::make_pair(-0.138570, eps), // h_{m}
                    std::make_pair(-0.062348, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.137702, eps), // h_{p}
                    std::make_pair(-0.112754, eps), // h_{m}
                    std::make_pair(+0.404549, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.037166, eps), // h_{p}
                    std::make_pair(-0.086163, eps), // h_{m}
                    std::make_pair(+1.271200, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.317099, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.273187, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.721643, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.311925, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.198352, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.544512, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), -0.043808, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.514150, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +0.952830, eps);
            }

            // using z_* with a = 1.0 and LP z-order = 4 and SLP z-order 2
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp"_ok,   "4" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToD, PToP> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1, LP z^4, SLP z^2" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps),

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+2.012713, eps), // w = 2.10
                    std::make_pair(+0.809594, eps), // w = 1.60
                    std::make_pair(+0.865962, eps), // w = 1.10
                    std::make_pair(+0.928873, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(+0.198603, eps), // w = 2.10
                    std::make_pair(+0.181937, eps), // w = 1.60
                    std::make_pair(+0.409565, eps), // w = 1.10
                    std::make_pair(+0.452445, eps), // w = 1.05
                    std::make_pair(+0.5,      eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-0.642637,  eps), // w = 2.10
                    std::make_pair(-0.545733,  eps), // w = 1.60
                    std::make_pair(-0.138032,  eps), // w = 1.10
                    std::make_pair(-0.0719429, eps), // w = 1.05
                    std::make_pair(+0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.412558, eps), // w = 2.10
                    std::make_pair(-0.250536, eps), // w = 1.60
                    std::make_pair(+0.133386, eps), // w = 1.10
                    std::make_pair(+0.189641, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(-0.177041, eps), // h_{p}
                    std::make_pair(-0.127333, eps), // h_{m}
                    std::make_pair(+0.112774, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.202289, eps), // h_{p}
                    std::make_pair(-0.108829, eps), // h_{m}
                    std::make_pair(+0.462087, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.037166, eps), // h_{p}
                    std::make_pair(-0.086163, eps), // h_{m}
                    std::make_pair(+1.271200, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.112387, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.335849, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.737288, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.121697, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.252192, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.557402, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), +0.145495, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.571629, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +0.967181, eps);
            }

            // using z_* with a = 1.25 and LP z-order = 4 and SLP z-order 2
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.25;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp"_ok,   "4" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToD, PToP> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1.25, LP z^4, SLP z^2" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps),

                    /* z(w) */
                    std::make_pair(-0.09904841, eps), // w = 1.10
                    std::make_pair(-0.10501000, eps), // w = 1.05
                    std::make_pair(-0.11111111, eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+2.029054, eps), // w = 2.10
                    std::make_pair(+0.810852, eps), // w = 1.60
                    std::make_pair(+0.865963, eps), // w = 1.10
                    std::make_pair(+0.928873, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(+0.212853, eps), // w = 2.10
                    std::make_pair(+0.184995, eps), // w = 1.60
                    std::make_pair(+0.409585, eps), // w = 1.10
                    std::make_pair(+0.452447, eps), // w = 1.05
                    std::make_pair(+0.5,      eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-0.6259680, eps), // w = 2.10
                    std::make_pair(-0.5421554, eps), // w = 1.60
                    std::make_pair(-0.1380090, eps), // w = 1.10
                    std::make_pair(-0.0719399, eps), // w = 1.05
                    std::make_pair(+0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.401804, eps), // w = 2.10
                    std::make_pair(-0.248228, eps), // w = 1.60
                    std::make_pair(+0.133401, eps), // w = 1.10
                    std::make_pair(+0.189643, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(-0.175024, eps), // h_{p}
                    std::make_pair(-0.127232, eps), // h_{m}
                    std::make_pair(+0.114625, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.202638, eps), // h_{p}
                    std::make_pair(-0.108808, eps), // h_{m}
                    std::make_pair(+0.462399, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.037166, eps), // h_{p}
                    std::make_pair(-0.086163, eps), // h_{m}
                    std::make_pair(+1.271200, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.110318, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.336163, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.737324, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.119776, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.252462, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.557432, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), +0.147437, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.571917, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +0.967214, eps);
            }

            // using z_* with a = 1.25 and LP z-order = 4 and SLP z-order 2, and l_3 to l_6 non-zero
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     = +1.2;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     = -2.2;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     = +2.1;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     = +3.1;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.25;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp"_ok,   "4" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToD, PToP> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1.25, LP z^4, SLP z^2" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps),

                    /* z(w) */
                    std::make_pair(-0.09904841, eps), // w = 1.10
                    std::make_pair(-0.10501000, eps), // w = 1.05
                    std::make_pair(-0.11111111, eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+2.029054, eps), // w = 2.10
                    std::make_pair(+0.810852, eps), // w = 1.60
                    std::make_pair(+0.865963, eps), // w = 1.10
                    std::make_pair(+0.928873, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(+0.212853, eps), // w = 2.10
                    std::make_pair(+0.184995, eps), // w = 1.60
                    std::make_pair(+0.409585, eps), // w = 1.10
                    std::make_pair(+0.452447, eps), // w = 1.05
                    std::make_pair(+0.5,      eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-0.6259680, eps), // w = 2.10
                    std::make_pair(-0.5421554, eps), // w = 1.60
                    std::make_pair(-0.1380090, eps), // w = 1.10
                    std::make_pair(-0.0719399, eps), // w = 1.05
                    std::make_pair(+0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.401804, eps), // w = 2.10
                    std::make_pair(-0.248228, eps), // w = 1.60
                    std::make_pair(+0.133401, eps), // w = 1.10
                    std::make_pair(+0.189643, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(-0.175024, eps), // h_{p}
                    std::make_pair(-0.177422, eps), // h_{m}
                    std::make_pair(+0.164815, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.202638, eps), // h_{p}
                    std::make_pair(-0.164244, eps), // h_{m}
                    std::make_pair(+0.517835, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.037166, eps), // h_{p}
                    std::make_pair(-0.158401, eps), // h_{m}
                    std::make_pair(+1.343437, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.083074, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.366777, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.772101, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.101892, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.262041, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.562339, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), +0.204498, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.636037, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +1.040053, eps);
            }
        }
} b_to_d_hqet_form_factors_test;

class BToDstarHQETFormFactorsTest :
    public TestCase
{
    public:
        BToDstarHQETFormFactorsTest() :
            TestCase("b_to_dstar_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            // using z_* with a = 1.0 and LP z-order = 2 and SLP z-order 2 and SSLP z-order 1
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -0.849472;
                p["B(*)->D(*)::xi''(1)@HQET"]    =  2.0 * 0.583711;
                p["B(*)->D(*)::xi'''(1)@HQET"]   =  0.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  =  0.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = -0.0600533;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  =  6.97061e-6;
                p["B(*)->D(*)::chi_2''(1)@HQET"] =  0.0314499;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  =  0.0400298;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = -0.039123;
                p["B(*)->D(*)::eta(1)@HQET"]     =  0.604052;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -0.00545745;
                p["B(*)->D(*)::eta''(1)@HQET"]   = -0.268764;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.111274;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.01963;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0687349;
                p["B(*)->D(*)::l_4(1)@HQET"]     = -2.02231;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  4.21978;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  4.52949;
                p["B(*)->D(*)::l_1'(1)@HQET"]    = -15.0241;
                p["B(*)->D(*)::l_2'(1)@HQET"]    = -9.43754;
                p["B(*)->D(*)::l_3'(1)@HQET"]    = -0.616533;
                p["B(*)->D(*)::l_4'(1)@HQET"]    = +0.604533;
                p["B(*)->D(*)::l_5'(1)@HQET"]    = +0.115125;
                p["B(*)->D(*)::l_6'(1)@HQET"]    = -1.4777;
                p["B(*)->D(*)::a@HQET"]          =  1.0;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp"_ok,   "2" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToDstar, PToV> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1, LP z^3, SLP z^1" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps),

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+0.541418, eps), // w = 2.10
                    std::make_pair(+0.656849, eps), // w = 1.60
                    std::make_pair(+0.920648, eps), // w = 1.10
                    std::make_pair(+0.958955, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(-0.0480609, eps), // w = 2.10
                    std::make_pair(-0.0557318, eps), // w = 1.60
                    std::make_pair(-0.0599029, eps), // w = 1.10
                    std::make_pair(-0.0600146, eps), // w = 1.05
                    std::make_pair(-0.0600533, eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(+0.027665,  eps), // w = 2.10
                    std::make_pair(+0.0183516, eps), // w = 1.60
                    std::make_pair(+0.00381496,eps), // w = 1.10
                    std::make_pair(+0.00195355,eps), // w = 1.05
                    std::make_pair( 0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(+0.495827, eps), // w = 2.10
                    std::make_pair(+0.563923, eps), // w = 1.60
                    std::make_pair(+0.602227, eps), // w = 1.10
                    std::make_pair(+0.603451, eps), // w = 1.05
                    std::make_pair(+0.604052, eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(+0.601021, eps), // h_{A_1}
                    std::make_pair(-0.180626, eps), // h_{A_2}
                    std::make_pair(+0.598766, eps), // h_{A_3}
                    std::make_pair(+0.693275, eps), // h_{V}
                    std::make_pair(+0.630681, eps), // h_{T_1}
                    std::make_pair(-0.070441, eps), // h_{T_2}
                    std::make_pair(-0.095175, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.727568, eps), // h_{A_1}
                    std::make_pair(-0.223398, eps), // h_{A_2}
                    std::make_pair(+0.712446, eps), // h_{A_3}
                    std::make_pair(+0.840248, eps), // h_{V}
                    std::make_pair(+0.770473, eps), // h_{T_1}
                    std::make_pair(-0.082214, eps), // h_{T_2}
                    std::make_pair(-0.128579, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.899260, eps), // h_{A_1}
                    std::make_pair(-0.282762, eps), // h_{A_2}
                    std::make_pair(+0.864723, eps), // h_{A_3}
                    std::make_pair(+1.041364, eps), // h_{V}
                    std::make_pair(+0.961350, eps), // h_{T_1}
                    std::make_pair(-0.096996, eps), // h_{T_2}
                    std::make_pair(-0.177567, eps), // h_{T_3}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);
            }

            // using z_* with a = 1.0 and LP z-order = 3 and SLP z-order 1
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0;

                auto oo = Options{
                    { "z-order-lp"_ok,   "3" },
                    { "z-order-slp"_ok,  "1" },
                    { "z-order-sslp"_ok, "1" },
                };
                HQETFormFactors<BToDstar, PToV> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1, LP z^3, SLP z^1" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+0.0, eps),

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+1.665540, eps), // w = 2.10
                    std::make_pair(+0.764544, eps), // w = 1.60
                    std::make_pair(+0.865908, eps), // w = 1.10
                    std::make_pair(+0.928869, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(-0.373019,  eps), // w = 2.10
                    std::make_pair(-0.0239773, eps), // w = 1.60
                    std::make_pair(+0.402425,  eps), // w = 1.10
                    std::make_pair(+0.450615,  eps), // w = 1.05
                    std::make_pair(+0.5,       eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-1.30953,   eps), // w = 2.10
                    std::make_pair(-0.785966,  eps), // w = 1.60
                    std::make_pair(-0.146363,  eps), // w = 1.10
                    std::make_pair(-0.0740769, eps), // w = 1.05
                    std::make_pair( 0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.841274, eps), // w = 2.10
                    std::make_pair(-0.404972, eps), // w = 1.60
                    std::make_pair(+0.128031, eps), // w = 1.10
                    std::make_pair(+0.188269, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(+0.684812, eps), // h_{A_1}
                    std::make_pair(-0.075634, eps), // h_{A_2}
                    std::make_pair(+0.702037, eps), // h_{A_3}
                    std::make_pair(+0.896998, eps), // h_{V}
                    std::make_pair(+0.686863, eps), // h_{T_1}
                    std::make_pair(-0.104330, eps), // h_{T_2}
                    std::make_pair(-0.195575, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.724634, eps), // h_{A_1}
                    std::make_pair(-0.027783, eps), // h_{A_2}
                    std::make_pair(+0.611764, eps), // h_{A_3}
                    std::make_pair(+0.965904, eps), // h_{V}
                    std::make_pair(+0.749274, eps), // h_{T_1}
                    std::make_pair(-0.133736, eps), // h_{T_2}
                    std::make_pair(-0.354839, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.899905, eps), // h_{A_1}
                    std::make_pair(+0.036348, eps), // h_{A_2}
                    std::make_pair(+0.552732, eps), // h_{A_3}
                    std::make_pair(+1.217624, eps), // h_{V}
                    std::make_pair(+0.961994, eps), // h_{T_1}
                    std::make_pair(-0.198494, eps), // h_{T_2}
                    std::make_pair(-0.665817, eps), // h_{T_3}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);
            }

            // using z_* with a = 1.0 and LP z-order = 4 and SLP z-order 2
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0;

                auto oo = Options{
                    { "z-order-lp"_ok,   "4" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToDstar, PToV> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1, LP z^4, SLP z^2" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps),

                    /* z(w) */
                    std::make_pair(0.01219690, eps), // w = 1.10
                    std::make_pair(0.00617307, eps), // w = 1.05
                    std::make_pair(0.0,        eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+2.012713, eps), // w = 2.10
                    std::make_pair(+0.809594, eps), // w = 1.60
                    std::make_pair(+0.865962, eps), // w = 1.10
                    std::make_pair(+0.928873, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(+0.198603, eps), // w = 2.10
                    std::make_pair(+0.181937, eps), // w = 1.60
                    std::make_pair(+0.409565, eps), // w = 1.10
                    std::make_pair(+0.452445, eps), // w = 1.05
                    std::make_pair(+0.5,      eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-0.642637,  eps), // w = 2.10
                    std::make_pair(-0.545733,  eps), // w = 1.60
                    std::make_pair(-0.138032,  eps), // w = 1.10
                    std::make_pair(-0.0719429, eps), // w = 1.05
                    std::make_pair(+0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.412558, eps), // w = 2.10
                    std::make_pair(-0.250536, eps), // w = 1.60
                    std::make_pair(+0.133386, eps), // w = 1.10
                    std::make_pair(+0.189641, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(+0.680851, eps), // h_{A_1}
                    std::make_pair(-0.034496, eps), // h_{A_2}
                    std::make_pair(+0.635846, eps), // h_{A_3}
                    std::make_pair(+0.891781, eps), // h_{V}
                    std::make_pair(+0.683839, eps), // h_{T_1}
                    std::make_pair(-0.111367, eps), // h_{T_2}
                    std::make_pair(-0.256559, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.721994, eps), // h_{A_1}
                    std::make_pair(-0.015242, eps), // h_{A_2}
                    std::make_pair(+0.589781, eps), // h_{A_3}
                    std::make_pair(+0.962020, eps), // h_{V}
                    std::make_pair(+0.746808, eps), // h_{T_1}
                    std::make_pair(-0.135515, eps), // h_{T_2}
                    std::make_pair(-0.372939, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.899905, eps), // h_{A_1}
                    std::make_pair(+0.036348, eps), // h_{A_2}
                    std::make_pair(+0.552732, eps), // h_{A_3}
                    std::make_pair(+1.217624, eps), // h_{V}
                    std::make_pair(+0.961994, eps), // h_{T_1}
                    std::make_pair(-0.198494, eps), // h_{T_2}
                    std::make_pair(-0.665817, eps), // h_{T_3}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);
            }

            // using z_* with a = 1.25 and LP z-order = 4 and SLP z-order 2
            {
                Parameters p = Parameters::Defaults();
                p["B(*)->D(*)::xi'(1)@HQET"]     = -1.5;
                p["B(*)->D(*)::xi''(1)@HQET"]    = +3.0;
                p["B(*)->D(*)::xi'''(1)@HQET"]   = +6.0;
                p["B(*)->D(*)::xi''''(1)@HQET"]  = -9.0;
                p["B(*)->D(*)::chi_2(1)@HQET"]   = +0.5;
                p["B(*)->D(*)::chi_2'(1)@HQET"]  = -1.0;
                p["B(*)->D(*)::chi_2''(1)@HQET"] = +2.0;
                p["B(*)->D(*)::chi_3'(1)@HQET"]  = -1.5;
                p["B(*)->D(*)::chi_3''(1)@HQET"] = +2.5;
                p["B(*)->D(*)::eta(1)@HQET"]     = +0.25;
                p["B(*)->D(*)::eta'(1)@HQET"]    = -1.25;
                p["B(*)->D(*)::eta''(1)@HQET"]   = +1.75;
                p["B(*)->D(*)::l_1(1)@HQET"]     = +0.5;
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_1'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_2'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_3'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_4'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_5'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::l_6'(1)@HQET"]    =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.25;

                auto oo = Options{
                    { "z-order-lp"_ok,   "4" },
                    { "z-order-slp"_ok,  "2" },
                    { "z-order-sslp"_ok, "1" }
                };
                HQETFormFactors<BToDstar, PToV> ff(p, oo);

                Diagnostics diag = ff.diagnostics();
                //std::cout << "a = 1.25, LP z^4, SLP z^2" << std::endl;
                //for (auto d : diag)
                //{
                //    std::cout << d.description << ": " << d.value << std::endl;
                //}
                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                    /* Options */
                    std::make_pair(+1.0, eps),
                    std::make_pair(+1.0, eps),
                    std::make_pair(+0.0, eps),
                    std::make_pair(+1.0, eps),

                    /* z(w) */
                    std::make_pair(-0.09904841, eps), // w = 1.10
                    std::make_pair(-0.10501000, eps), // w = 1.05
                    std::make_pair(-0.11111111, eps), // w = 1.00

                    /* xi(w) */
                    std::make_pair(+2.029054, eps), // w = 2.10
                    std::make_pair(+0.810852, eps), // w = 1.60
                    std::make_pair(+0.865963, eps), // w = 1.10
                    std::make_pair(+0.928873, eps), // w = 1.05
                    std::make_pair(+1.000000, eps), // w = 1.00

                    /* chi2(w) */
                    std::make_pair(+0.212853, eps), // w = 2.10
                    std::make_pair(+0.184995, eps), // w = 1.60
                    std::make_pair(+0.409585, eps), // w = 1.10
                    std::make_pair(+0.452447, eps), // w = 1.05
                    std::make_pair(+0.5,      eps), // w = 1.00

                    /* chi3(w) */
                    std::make_pair(-0.6259680, eps), // w = 2.10
                    std::make_pair(-0.5421554, eps), // w = 1.60
                    std::make_pair(-0.1380090, eps), // w = 1.10
                    std::make_pair(-0.0719399, eps), // w = 1.05
                    std::make_pair(+0.0,       eps), // w = 1.00

                    /* eta(w) */
                    std::make_pair(-0.401804, eps), // w = 2.10
                    std::make_pair(-0.248228, eps), // w = 1.60
                    std::make_pair(+0.133401, eps), // w = 1.10
                    std::make_pair(+0.189643, eps), // w = 1.05
                    std::make_pair(+0.25,     eps), // w = 1.00

                    /* r(w) */
                    std::make_pair(+0.967945, eps), // w = 1.1
                    std::make_pair(+0.999767, eps), // w = 1.0007
                    std::make_pair(+0.999967, eps), // w = 1.0001
                    std::make_pair(+0.999983, eps), // w = 1.00005
                    std::make_pair(+1.0,      eps), // w = 1.0

                    /* Omega(w, z = 0.25) */
                    std::make_pair(+1.294026, eps), // w = 1.1
                    std::make_pair(+1.310389, eps), // w = 1.0007
                    std::make_pair(+1.310476, eps), // w = 1.0001
                    std::make_pair(+1.310483, eps), // w = 1.00005
                    std::make_pair(+1.310491, eps), // w = 1.0

                    /* Omega(w, z = 0.20) */
                    std::make_pair(+1.403808, eps), // w = 1.1
                    std::make_pair(+1.414099, eps), // w = 1.0007
                    std::make_pair(+1.414149, eps), // w = 1.0001
                    std::make_pair(+1.414153, eps), // w = 1.00005
                    std::make_pair(+1.414157, eps), // w = 1.0

                    /* WCs at (w = 1.2, z = 0.20) */
                    std::make_pair(-0.591250, eps), // C_{S  }
                    std::make_pair(+0.659746, eps), // C_{P  }
                    std::make_pair(+1.123905, eps), // C_{V_1}
                    std::make_pair(-0.454499, eps), // C_{V_2}
                    std::make_pair(-0.162046, eps), // C_{V_3}
                    std::make_pair(-0.127091, eps), // C_{A_1}
                    std::make_pair(-1.247185, eps), // C_{A_2}
                    std::make_pair( 0.316106, eps), // C_{A_3}
                    std::make_pair(+0.694295, eps), // C_{T_1}
                    std::make_pair(-0.931381, eps), // C_{T_2}
                    std::make_pair( 0.319615, eps), // C_{T_3}

                    /* WCs at (w = 1.0, z = 0.25) */
                    std::make_pair(-0.666667, eps), // C_{S  }
                    std::make_pair(+0.666667, eps), // C_{P  }
                    std::make_pair(+0.977157, eps), // C_{V_1}
                    std::make_pair(-0.478135, eps), // C_{V_2}
                    std::make_pair(-0.188532, eps), // C_{V_3}
                    std::make_pair(-0.356176, eps), // C_{A_1}
                    std::make_pair(-1.250411, eps), // C_{A_2}
                    std::make_pair( 0.381601, eps), // C_{A_3}
                    std::make_pair(+0.413987, eps), // C_{T_1}
                    std::make_pair(-0.956270, eps), // C_{T_2}
                    std::make_pair( 0.377063, eps), // C_{T_3}

                    /* HQET form factors at w = 1.4 */
                    std::make_pair(+0.680901, eps), // h_{A_1}
                    std::make_pair(-0.034071, eps), // h_{A_2}
                    std::make_pair(+0.635238, eps), // h_{A_3}
                    std::make_pair(+0.891845, eps), // h_{V}
                    std::make_pair(+0.683899, eps), // h_{T_1}
                    std::make_pair(-0.111456, eps), // h_{T_2}
                    std::make_pair(-0.257232, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.721983, eps), // h_{A_1}
                    std::make_pair(-0.015175, eps), // h_{A_2}
                    std::make_pair(+0.589665, eps), // h_{A_3}
                    std::make_pair(+0.962004, eps), // h_{V}
                    std::make_pair(+0.746799, eps), // h_{T_1}
                    std::make_pair(-0.135525, eps), // h_{T_2}
                    std::make_pair(-0.373038, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.899905, eps), // h_{A_1}
                    std::make_pair(+0.036348, eps), // h_{A_2}
                    std::make_pair(+0.552732, eps), // h_{A_3}
                    std::make_pair(+1.217624, eps), // h_{V}
                    std::make_pair(+0.961994, eps), // h_{T_1}
                    std::make_pair(-0.198494, eps), // h_{T_2}
                    std::make_pair(-0.665817, eps), // h_{T_3}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);
            }
        }
} b_to_dstar_hqet_form_factors_test;

class BstarToDHQETFormFactorsTest :
    public TestCase
{
    public:
        BstarToDHQETFormFactorsTest() :
            TestCase("bstar_to_d_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            Parameters p = Parameters::Defaults();
            p["B(*)->D(*)::xi'(1)@HQET"]     = -0.849472;
            p["B(*)->D(*)::xi''(1)@HQET"]    =  2.0 * 0.583711;
            p["B(*)->D(*)::xi'''(1)@HQET"]   =  0.0;
            p["B(*)->D(*)::xi''''(1)@HQET"]  =  0.0;
            p["B(*)->D(*)::chi_2(1)@HQET"]   = -0.0600533;
            p["B(*)->D(*)::chi_2'(1)@HQET"]  =  6.97061e-6;
            p["B(*)->D(*)::chi_2''(1)@HQET"] =  0.0314499;
            p["B(*)->D(*)::chi_3'(1)@HQET"]  =  0.0400298;
            p["B(*)->D(*)::chi_3''(1)@HQET"] = -0.039123;
            p["B(*)->D(*)::eta(1)@HQET"]     =  0.604052;
            p["B(*)->D(*)::eta'(1)@HQET"]    = -0.00545745;
            p["B(*)->D(*)::eta''(1)@HQET"]   = -0.268764;
            p["B(*)->D(*)::l_1(1)@HQET"]     = +0.111274;
            p["B(*)->D(*)::l_2(1)@HQET"]     = -2.01963;
            p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0687349;
            p["B(*)->D(*)::l_4(1)@HQET"]     = -2.02231;
            p["B(*)->D(*)::l_5(1)@HQET"]     =  4.21978;
            p["B(*)->D(*)::l_6(1)@HQET"]     =  4.52949;
            p["B(*)->D(*)::l_1'(1)@HQET"]    = -15.0241;
            p["B(*)->D(*)::l_2'(1)@HQET"]    = -9.43754;
            p["B(*)->D(*)::l_3'(1)@HQET"]    = -0.616533;
            p["B(*)->D(*)::l_4'(1)@HQET"]    = +0.604533;
            p["B(*)->D(*)::l_5'(1)@HQET"]    = +0.115125;
            p["B(*)->D(*)::l_6'(1)@HQET"]    = -1.4777;
            p["B(*)->D(*)::a@HQET"]          =  1.0;
            p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
            p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

            auto oo = Options{
                { "z-order-lp"_ok,   "2" },
                { "z-order-slp"_ok,  "2" },
                { "z-order-sslp"_ok, "1" }
            };
            HQETFormFactors<BstarToD, VToP> ff(p, oo);

            Diagnostics diag = ff.diagnostics();
            //for (auto d : diag)
            //{
            //    std::cout << d.description << ": " << d.value << std::endl;
            //}
            static const std::vector<std::pair<double, double>> ref
            {
                /* Inputs */
                std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                /* Options */
                std::make_pair(+0.0, eps), // LP z^3 terms enabled?
                std::make_pair(+0.0, eps),
                std::make_pair(+0.0, eps),
                std::make_pair(+1.0, eps), // SLP z^2 terms enabled?

                /* z(w) */
                std::make_pair(0.01219690, eps), // w = 1.10
                std::make_pair(0.00617307, eps), // w = 1.05
                std::make_pair(0.0,        eps), // w = 1.00

                /* xi(w) */
                std::make_pair(+0.541418, eps), // w = 2.10
                std::make_pair(+0.656849, eps), // w = 1.60
                std::make_pair(+0.920648, eps), // w = 1.10
                std::make_pair(+0.958955, eps), // w = 1.05
                std::make_pair(+1.000000, eps), // w = 1.00

                /* chi2(w) */
                std::make_pair(-0.0480609, eps), // w = 2.10
                std::make_pair(-0.0557318, eps), // w = 1.60
                std::make_pair(-0.0599029, eps), // w = 1.10
                std::make_pair(-0.0600146, eps), // w = 1.05
                std::make_pair(-0.0600533, eps), // w = 1.00

                /* chi3(w) */
                std::make_pair(+0.027665,  eps), // w = 2.10
                std::make_pair(+0.0183516, eps), // w = 1.60
                std::make_pair(+0.00381496,eps), // w = 1.10
                std::make_pair(+0.00195355,eps), // w = 1.05
                std::make_pair( 0.0,       eps), // w = 1.00

                /* eta(w) */
                std::make_pair(+0.495827, eps), // w = 2.10
                std::make_pair(+0.563923, eps), // w = 1.60
                std::make_pair(+0.602227, eps), // w = 1.10
                std::make_pair(+0.603451, eps), // w = 1.05
                std::make_pair(+0.604052, eps), // w = 1.00

                /* r(w) */
                std::make_pair(+0.967945, eps), // w = 1.1
                std::make_pair(+0.999767, eps), // w = 1.0007
                std::make_pair(+0.999967, eps), // w = 1.0001
                std::make_pair(+0.999983, eps), // w = 1.00005
                std::make_pair(+1.0,      eps), // w = 1.0

                /* Omega(w, z = 0.25) */
                std::make_pair(+1.294026, eps), // w = 1.1
                std::make_pair(+1.310389, eps), // w = 1.0007
                std::make_pair(+1.310476, eps), // w = 1.0001
                std::make_pair(+1.310483, eps), // w = 1.00005
                std::make_pair(+1.310491, eps), // w = 1.0

                /* Omega(w, z = 0.20) */
                std::make_pair(+1.403808, eps), // w = 1.1
                std::make_pair(+1.414099, eps), // w = 1.0007
                std::make_pair(+1.414149, eps), // w = 1.0001
                std::make_pair(+1.414153, eps), // w = 1.00005
                std::make_pair(+1.414157, eps), // w = 1.0

                /* WCs at (w = 1.2, z = 0.20) */
                std::make_pair(-0.591250, eps), // C_{S  }
                std::make_pair(+0.659746, eps), // C_{P  }
                std::make_pair(+1.123905, eps), // C_{V_1}
                std::make_pair(-0.454499, eps), // C_{V_2}
                std::make_pair(-0.162046, eps), // C_{V_3}
                std::make_pair(-0.127091, eps), // C_{A_1}
                std::make_pair(-1.247185, eps), // C_{A_2}
                std::make_pair( 0.316106, eps), // C_{A_3}
                std::make_pair(+0.694295, eps), // C_{T_1}
                std::make_pair(-0.931381, eps), // C_{T_2}
                std::make_pair( 0.319615, eps), // C_{T_3}

                /* WCs at (w = 1.0, z = 0.25) */
                std::make_pair(-0.666667, eps), // C_{S  }
                std::make_pair(+0.666667, eps), // C_{P  }
                std::make_pair(+0.977157, eps), // C_{V_1}
                std::make_pair(-0.478135, eps), // C_{V_2}
                std::make_pair(-0.188532, eps), // C_{V_3}
                std::make_pair(-0.356176, eps), // C_{A_1}
                std::make_pair(-1.250411, eps), // C_{A_2}
                std::make_pair( 0.381601, eps), // C_{A_3}
                std::make_pair(+0.413987, eps), // C_{T_1}
                std::make_pair(-0.956270, eps), // C_{T_2}
                std::make_pair( 0.377063, eps), // C_{T_3}

                /* HQET form factors at w = 1.4 */
                std::make_pair(+0.638478, eps), // h_{Abar1}
                std::make_pair(-0.082947, eps), // h_{Abar2}
                std::make_pair(+0.712066, eps), // h_{Abar3}
                std::make_pair(+0.760677, eps), // h_{Vbar}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.779087, eps), // h_{Abar1}
                std::make_pair(-0.103046, eps), // h_{Abar2}
                std::make_pair(+0.866333, eps), // h_{Abar3}
                std::make_pair(+0.928731, eps), // h_{Vbar}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.969229, eps), // h_{Abar1}
                std::make_pair(-0.130890, eps), // h_{Abar2}
                std::make_pair(+1.078436, eps), // h_{Abar3}
                std::make_pair(+1.160604, eps), // h_{Vbar}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} bstar_to_d_hqet_form_factors_test;

class BstarToDstarHQETFormFactorsTest :
    public TestCase
{
    public:
        BstarToDstarHQETFormFactorsTest() :
            TestCase("bstar_to_dstar_hqet_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-6;

            Parameters p = Parameters::Defaults();
            p["B(*)->D(*)::xi'(1)@HQET"]     = -0.849472;
            p["B(*)->D(*)::xi''(1)@HQET"]    =  2.0 * 0.583711;
            p["B(*)->D(*)::xi'''(1)@HQET"]   =  0.0;
            p["B(*)->D(*)::xi''''(1)@HQET"]  =  0.0;
            p["B(*)->D(*)::chi_2(1)@HQET"]   = -0.0600533;
            p["B(*)->D(*)::chi_2'(1)@HQET"]  =  6.97061e-6;
            p["B(*)->D(*)::chi_2''(1)@HQET"] =  0.0314499;
            p["B(*)->D(*)::chi_3'(1)@HQET"]  =  0.0400298;
            p["B(*)->D(*)::chi_3''(1)@HQET"] = -0.039123;
            p["B(*)->D(*)::eta(1)@HQET"]     =  0.604052;
            p["B(*)->D(*)::eta'(1)@HQET"]    = -0.00545745;
            p["B(*)->D(*)::eta''(1)@HQET"]   = -0.268764;
            p["B(*)->D(*)::l_1(1)@HQET"]     = +0.111274;
            p["B(*)->D(*)::l_2(1)@HQET"]     = -2.01963;
            p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0687349;
            p["B(*)->D(*)::l_4(1)@HQET"]     = -2.02231;
            p["B(*)->D(*)::l_5(1)@HQET"]     =  4.21978;
            p["B(*)->D(*)::l_6(1)@HQET"]     =  4.52949;
            p["B(*)->D(*)::l_1'(1)@HQET"]    = -15.0241;
            p["B(*)->D(*)::l_2'(1)@HQET"]    = -9.43754;
            p["B(*)->D(*)::l_3'(1)@HQET"]    = -0.616533;
            p["B(*)->D(*)::l_4'(1)@HQET"]    = +0.604533;
            p["B(*)->D(*)::l_5'(1)@HQET"]    = +0.115125;
            p["B(*)->D(*)::l_6'(1)@HQET"]    = -1.4777;
            p["B(*)->D(*)::a@HQET"]          =  1.0;
            p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
            p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

            auto oo = Options{
                { "z-order-lp"_ok,   "2" },
                { "z-order-slp"_ok,  "2" },
                { "z-order-sslp"_ok, "1" }
            };
            HQETFormFactors<BstarToDstar, VToV> ff(p, oo);

            Diagnostics diag = ff.diagnostics();
            //for (auto d : diag)
            //{
            //    std::cout << d.description << ": " << d.value << std::endl;
            //}
            static const std::vector<std::pair<double, double>> ref
            {
                /* Inputs */
                std::make_pair(+0.288815, eps), // z  = m_c^1S / m_b^1S
                std::make_pair(+1.875618, eps), // wz = 1/2 (z + 1/z)

                /* Options */
                std::make_pair(+0.0, eps), // LP z^3 terms enabled?
                std::make_pair(+0.0, eps),
                std::make_pair(+0.0, eps),
                std::make_pair(+1.0, eps), // SLP z^2 terms enabled?

                /* z(w) */
                std::make_pair(0.01219690, eps), // w = 1.10
                std::make_pair(0.00617307, eps), // w = 1.05
                std::make_pair(0.0,        eps), // w = 1.00

                /* xi(w) */
                std::make_pair(+0.541418, eps), // w = 2.10
                std::make_pair(+0.656849, eps), // w = 1.60
                std::make_pair(+0.920648, eps), // w = 1.10
                std::make_pair(+0.958955, eps), // w = 1.05
                std::make_pair(+1.000000, eps), // w = 1.00

                /* chi2(w) */
                std::make_pair(-0.0480609, eps), // w = 2.10
                std::make_pair(-0.0557318, eps), // w = 1.60
                std::make_pair(-0.0599029, eps), // w = 1.10
                std::make_pair(-0.0600146, eps), // w = 1.05
                std::make_pair(-0.0600533, eps), // w = 1.00

                /* chi3(w) */
                std::make_pair(+0.027665,  eps), // w = 2.10
                std::make_pair(+0.0183516, eps), // w = 1.60
                std::make_pair(+0.00381496,eps), // w = 1.10
                std::make_pair(+0.00195355,eps), // w = 1.05
                std::make_pair( 0.0,       eps), // w = 1.00

                /* eta(w) */
                std::make_pair(+0.495827, eps), // w = 2.10
                std::make_pair(+0.563923, eps), // w = 1.60
                std::make_pair(+0.602227, eps), // w = 1.10
                std::make_pair(+0.603451, eps), // w = 1.05
                std::make_pair(+0.604052, eps), // w = 1.00

                /* r(w) */
                std::make_pair(+0.967945, eps), // w = 1.1
                std::make_pair(+0.999767, eps), // w = 1.0007
                std::make_pair(+0.999967, eps), // w = 1.0001
                std::make_pair(+0.999983, eps), // w = 1.00005
                std::make_pair(+1.0,      eps), // w = 1.0

                /* Omega(w, z = 0.25) */
                std::make_pair(+1.294026, eps), // w = 1.1
                std::make_pair(+1.310389, eps), // w = 1.0007
                std::make_pair(+1.310476, eps), // w = 1.0001
                std::make_pair(+1.310483, eps), // w = 1.00005
                std::make_pair(+1.310491, eps), // w = 1.0

                /* Omega(w, z = 0.20) */
                std::make_pair(+1.403808, eps), // w = 1.1
                std::make_pair(+1.414099, eps), // w = 1.0007
                std::make_pair(+1.414149, eps), // w = 1.0001
                std::make_pair(+1.414153, eps), // w = 1.00005
                std::make_pair(+1.414157, eps), // w = 1.0

                /* WCs at (w = 1.2, z = 0.20) */
                std::make_pair(-0.591250, eps), // C_{S  }
                std::make_pair(+0.659746, eps), // C_{P  }
                std::make_pair(+1.123905, eps), // C_{V_1}
                std::make_pair(-0.454499, eps), // C_{V_2}
                std::make_pair(-0.162046, eps), // C_{V_3}
                std::make_pair(-0.127091, eps), // C_{A_1}
                std::make_pair(-1.247185, eps), // C_{A_2}
                std::make_pair( 0.316106, eps), // C_{A_3}
                std::make_pair(+0.694295, eps), // C_{T_1}
                std::make_pair(-0.931381, eps), // C_{T_2}
                std::make_pair( 0.319615, eps), // C_{T_3}

                /* WCs at (w = 1.0, z = 0.25) */
                std::make_pair(-0.666667, eps), // C_{S  }
                std::make_pair(+0.666667, eps), // C_{P  }
                std::make_pair(+0.977157, eps), // C_{V_1}
                std::make_pair(-0.478135, eps), // C_{V_2}
                std::make_pair(-0.188532, eps), // C_{V_3}
                std::make_pair(-0.356176, eps), // C_{A_1}
                std::make_pair(-1.250411, eps), // C_{A_2}
                std::make_pair( 0.381601, eps), // C_{A_3}
                std::make_pair(+0.413987, eps), // C_{T_1}
                std::make_pair(-0.956270, eps), // C_{T_2}
                std::make_pair( 0.377063, eps), // C_{T_3}

                /* HQET form factors at w = 1.4 */
                std::make_pair(+0.614822, eps), // h_{1}
                std::make_pair(-0.007155, eps), // h_{2}
                std::make_pair(+0.843665, eps), // h_{3}
                std::make_pair(+0.768290, eps), // h_{4}
                std::make_pair(+0.069277, eps), // h_{5}
                std::make_pair(+0.053484, eps), // h_{6}
                std::make_pair(+0.567453, eps), // h_{7}
                std::make_pair(-0.043281, eps), // h_{8}
                std::make_pair(+0.110892, eps), // h_{9}
                std::make_pair(+0.064994, eps), // h_{10}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.757916, eps), // h_{1}
                std::make_pair(-0.009228, eps), // h_{2}
                std::make_pair(+1.039112, eps), // h_{3}
                std::make_pair(+0.937434, eps), // h_{4}
                std::make_pair(+0.097394, eps), // h_{5}
                std::make_pair(+0.067764, eps), // h_{6}
                std::make_pair(+0.707043, eps), // h_{7}
                std::make_pair(-0.048579, eps), // h_{8}
                std::make_pair(+0.147242, eps), // h_{9}
                std::make_pair(+0.081839, eps), // h_{10}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.954434, eps), // h_{1}
                std::make_pair(-0.012042, eps), // h_{2}
                std::make_pair(+1.303165, eps), // h_{3}
                std::make_pair(+1.167815, eps), // h_{4}
                std::make_pair(+0.139116, eps), // h_{5}
                std::make_pair(+0.088099, eps), // h_{6}
                std::make_pair(+0.899260, eps), // h_{7}
                std::make_pair(-0.055498, eps), // h_{8}
                std::make_pair(+0.200220, eps), // h_{9}
                std::make_pair(+0.105747, eps), // h_{10}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} bstar_to_dstar_hqet_form_factors_test;
