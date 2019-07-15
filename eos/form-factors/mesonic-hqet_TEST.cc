/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Danny van Dyk
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
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/mesonic-hqet.hh>
#include <eos/form-factors/mesonic-impl.hh>

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
                    { "z-order-lp",   "2" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(+0.672776, eps), // h_{p}
                    std::make_pair(-0.042011,eps), // h_{m}
                    std::make_pair(+0.777201, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.821816, eps), // h_{p}
                    std::make_pair(-0.048623, eps), // h_{m}
                    std::make_pair(+0.937161, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.023833, eps), // h_{p}
                    std::make_pair(-0.059743, eps), // h_{m}
                    std::make_pair(+1.157401, eps), // h_{T}
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
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp",   "3" },
                    { "z-order-slp",  "1" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(-0.360342, eps), // h_{p}
                    std::make_pair(-0.135928, eps), // h_{m}
                    std::make_pair(-0.053639, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.144245, eps), // h_{p}
                    std::make_pair(-0.110649, eps), // h_{m}
                    std::make_pair(+0.409050, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.036235, eps), // h_{p}
                    std::make_pair(-0.084644, eps), // h_{m}
                    std::make_pair(+1.268821, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.306008,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.279042,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.723895,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.301328,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.203886,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.546977,  eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), -0.0341907, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.5188460, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +0.9540820, eps);
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
                p["B(*)->D(*)::a@HQET"]          =  1.0;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp",   "4" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(-0.167795, eps), // h_{p}
                    std::make_pair(-0.124919, eps), // h_{m}
                    std::make_pair(+0.119669, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.208125, eps), // h_{p}
                    std::make_pair(-0.106803, eps), // h_{m}
                    std::make_pair(+0.465961, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.036235, eps), // h_{p}
                    std::make_pair(-0.084644, eps), // h_{m}
                    std::make_pair(+1.268821, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.103395, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.341034, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.739372, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.113066, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.257142, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.559725, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), +0.153143, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.575696, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +0.968275, eps);
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
                p["B(*)->D(*)::a@HQET"]          =  1.25;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp",   "4" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(-0.165798, eps), // h_{p}
                    std::make_pair(-0.124819, eps), // h_{m}
                    std::make_pair(+0.121503, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.208470, eps), // h_{p}
                    std::make_pair(-0.106783, eps), // h_{m}
                    std::make_pair(+0.466270, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.036235, eps), // h_{p}
                    std::make_pair(-0.084644, eps), // h_{m}
                    std::make_pair(+1.268821, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.101346, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.341344, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.739408, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.111165, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.257408, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.559754, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), +0.155065, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.575981, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +0.968308, eps);
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
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     = +1.2;
                p["B(*)->D(*)::l_4(1)@HQET"]     = -2.2;
                p["B(*)->D(*)::l_5(1)@HQET"]     = +2.1;
                p["B(*)->D(*)::l_6(1)@HQET"]     = +3.1;
                p["B(*)->D(*)::a@HQET"]          =  1.25;
                p["mass::B_d"]                   =  5.27942; // mixture of B0 and B+ masses
                p["mass::D_u"]                   =  1.86723; // mixture of D0 and D+ masses

                auto oo = Options{
                    { "z-order-lp",   "4" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(-0.165798, eps), // h_{p}
                    std::make_pair(-0.173589, eps), // h_{m}
                    std::make_pair(+0.170272, eps), // h_{T}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.208470, eps), // h_{p}
                    std::make_pair(-0.160648, eps), // h_{m}
                    std::make_pair(+0.520135, eps), // h_{T}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+1.036235, eps), // h_{p}
                    std::make_pair(-0.154835, eps), // h_{m}
                    std::make_pair(+1.339013, eps), // h_{T}
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 4.0), -0.0748736, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 8.0), +0.3710917, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), +0.7732000, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 4.0), -0.0937868, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0( 8.0), +0.2667160, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_0(10.0), +0.5645230, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 4.0), +0.210510,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t( 8.0), +0.638286,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_t(10.0), +1.039084,  eps);
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
                    { "z-order-lp",   "2" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(+0.604662, eps), // h_{A_1}
                    std::make_pair(-0.180031, eps), // h_{A_2}
                    std::make_pair(+0.602671, eps), // h_{A_3}
                    std::make_pair(+0.697766, eps), // h_{V}
                    std::make_pair(+0.634045, eps), // h_{T_1}
                    std::make_pair(-0.071109, eps), // h_{T_2}
                    std::make_pair(-0.095800, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.730339, eps), // h_{A_1}
                    std::make_pair(-0.222682, eps), // h_{A_2}
                    std::make_pair(+0.715769, eps), // h_{A_3}
                    std::make_pair(+0.844071, eps), // h_{V}
                    std::make_pair(+0.773011, eps), // h_{T_1}
                    std::make_pair(-0.082976, eps), // h_{T_2}
                    std::make_pair(-0.129127, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.900672, eps), // h_{A_1}
                    std::make_pair(-0.281876, eps), // h_{A_2}
                    std::make_pair(+0.867138, eps), // h_{A_3}
                    std::make_pair(+1.044116, eps), // h_{V}
                    std::make_pair(+0.962607, eps), // h_{T_1}
                    std::make_pair(-0.097881, eps), // h_{T_2}
                    std::make_pair(-0.177963, eps), // h_{T_3}
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
                p["B(*)->D(*)::l_2(1)@HQET"]     = -2.0;
                p["B(*)->D(*)::l_3(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_4(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_5(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::l_6(1)@HQET"]     =  0.0;
                p["B(*)->D(*)::a@HQET"]          =  1.0;

                auto oo = Options{
                    { "z-order-lp",   "3" },
                    { "z-order-slp",  "1" },
                    { "z-order-sslp", "1" },
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(+0.681567, eps), // h_{A_1}
                    std::make_pair(-0.075176, eps), // h_{A_2}
                    std::make_pair(+0.699622, eps), // h_{A_3}
                    std::make_pair(+0.892280, eps), // h_{V}
                    std::make_pair(+0.683813, eps), // h_{T_1}
                    std::make_pair(-0.102237, eps), // h_{T_2}
                    std::make_pair(-0.193310, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.723234, eps), // h_{A_1}
                    std::make_pair(-0.028124, eps), // h_{A_2}
                    std::make_pair(+0.613043, eps), // h_{A_3}
                    std::make_pair(+0.962703, eps), // h_{V}
                    std::make_pair(+0.747939, eps), // h_{T_1}
                    std::make_pair(-0.131398, eps), // h_{T_2}
                    std::make_pair(-0.350402, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.901299, eps), // h_{A_1}
                    std::make_pair(+0.034805, eps), // h_{A_2}
                    std::make_pair(+0.560228, eps), // h_{A_3}
                    std::make_pair(+1.216434, eps), // h_{V}
                    std::make_pair(+0.963234, eps), // h_{T_1}
                    std::make_pair(-0.195455, eps), // h_{T_2}
                    std::make_pair(-0.657192, eps), // h_{T_3}
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
                p["B(*)->D(*)::a@HQET"]          =  1.0;

                auto oo = Options{
                    { "z-order-lp",   "4" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(+0.678391, eps), // h_{A_1}
                    std::make_pair(-0.034634, eps), // h_{A_2}
                    std::make_pair(+0.635054, eps), // h_{A_3}
                    std::make_pair(+0.887824, eps), // h_{V}
                    std::make_pair(+0.681577, eps), // h_{T_1}
                    std::make_pair(-0.109242, eps), // h_{T_2}
                    std::make_pair(-0.253432, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.720845, eps), // h_{A_1}
                    std::make_pair(-0.015762, eps), // h_{A_2}
                    std::make_pair(+0.591566, eps), // h_{A_3}
                    std::make_pair(+0.959068, eps), // h_{V}
                    std::make_pair(+0.745724, eps), // h_{T_1}
                    std::make_pair(-0.133174, eps), // h_{T_2}
                    std::make_pair(-0.368245, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.901299, eps), // h_{A_1}
                    std::make_pair(+0.034805, eps), // h_{A_2}
                    std::make_pair(+0.560228, eps), // h_{A_3}
                    std::make_pair(+1.216434, eps), // h_{V}
                    std::make_pair(+0.963234, eps), // h_{T_1}
                    std::make_pair(-0.195455, eps), // h_{T_2}
                    std::make_pair(-0.657192, eps), // h_{T_3}
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
                p["B(*)->D(*)::a@HQET"]          =  1.25;

                auto oo = Options{
                    { "z-order-lp",   "4" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
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
                    std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                    std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                    std::make_pair(+0.678448, eps), // h_{A_1}
                    std::make_pair(-0.034215, eps), // h_{A_2}
                    std::make_pair(+0.634463, eps), // h_{A_3}
                    std::make_pair(+0.887896, eps), // h_{V}
                    std::make_pair(+0.681645, eps), // h_{T_1}
                    std::make_pair(-0.109330, eps), // h_{T_2}
                    std::make_pair(-0.254095, eps), // h_{T_3}

                    /* HQET form factors at w = 1.2 */
                    std::make_pair(+0.720835, eps), // h_{A_1}
                    std::make_pair(-0.015696, eps), // h_{A_2}
                    std::make_pair(+0.591454, eps), // h_{A_3}
                    std::make_pair(+0.959053, eps), // h_{V}
                    std::make_pair(+0.745716, eps), // h_{T_1}
                    std::make_pair(-0.133184, eps), // h_{T_2}
                    std::make_pair(-0.368343, eps), // h_{T_3}

                    /* HQET form factors at w = 1.0 */
                    std::make_pair(+0.901299, eps), // h_{A_1}
                    std::make_pair(+0.034805, eps), // h_{A_2}
                    std::make_pair(+0.560228, eps), // h_{A_3}
                    std::make_pair(+1.216434, eps), // h_{V}
                    std::make_pair(+0.963234, eps), // h_{T_1}
                    std::make_pair(-0.195455, eps), // h_{T_2}
                    std::make_pair(-0.657192, eps), // h_{T_3}
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
                { "z-order-lp",   "2" },
                { "z-order-slp",  "2" },
                { "z-order-sslp", "1" }
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
                std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                std::make_pair(+0.641146, eps), // h_{Abar1}
                std::make_pair(-0.083134, eps), // h_{Abar2}
                std::make_pair(+0.713701, eps), // h_{Abar3}
                std::make_pair(+0.762580, eps), // h_{Vbar}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.780449, eps), // h_{Abar1}
                std::make_pair(-0.103274, eps), // h_{Abar2}
                std::make_pair(+0.866394, eps), // h_{Abar3}
                std::make_pair(+0.929112, eps), // h_{Vbar}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.968659, eps), // h_{Abar1}
                std::make_pair(-0.131176, eps), // h_{Abar2}
                std::make_pair(+1.076130, eps), // h_{Abar3}
                std::make_pair(+1.158692, eps), // h_{Vbar}
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
                { "z-order-lp",   "2" },
                { "z-order-slp",  "2" },
                { "z-order-sslp", "1" }
            };
            HQETFormFactors<BstarToDstar, VToV> ff(p, oo);

            Diagnostics diag = ff.diagnostics();
            for (auto d : diag)
            {
                std::cout << d.description << ": " << d.value << std::endl;
            }
            static const std::vector<std::pair<double, double>> ref
            {
                /* Inputs */
                std::make_pair(+0.292994, eps), // z  = m_c^1S / m_b^1S
                std::make_pair(+1.853019, eps), // wz = 1/2 (z + 1/z)

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
                std::make_pair(+0.618277, eps), // h_{1}
                std::make_pair(-0.008016, eps), // h_{2}
                std::make_pair(+0.846868, eps), // h_{3}
                std::make_pair(+0.772781, eps), // h_{4}
                std::make_pair(+0.069989, eps), // h_{5}
                std::make_pair(+0.053556, eps), // h_{6}
                std::make_pair(+0.570940, eps), // h_{7}
                std::make_pair(-0.043757, eps), // h_{8}
                std::make_pair(+0.111398, eps), // h_{9}
                std::make_pair(+0.065181, eps), // h_{10}

                /* HQET form factors at w = 1.2 */
                std::make_pair(+0.760573, eps), // h_{1}
                std::make_pair(-0.010198, eps), // h_{2}
                std::make_pair(+1.041098, eps), // h_{3}
                std::make_pair(+0.941258, eps), // h_{4}
                std::make_pair(+0.098042, eps), // h_{5}
                std::make_pair(+0.067854, eps), // h_{6}
                std::make_pair(+0.709718, eps), // h_{7}
                std::make_pair(-0.049133, eps), // h_{8}
                std::make_pair(+0.147650, eps), // h_{9}
                std::make_pair(+0.082067, eps), // h_{10}

                /* HQET form factors at w = 1.0 */
                std::make_pair(+0.955846, eps), // h_{1}
                std::make_pair(-0.013154, eps), // h_{2}
                std::make_pair(+1.303370, eps), // h_{3}
                std::make_pair(+1.170567, eps), // h_{4}
                std::make_pair(+0.139625, eps), // h_{5}
                std::make_pair(+0.088212, eps), // h_{6}
                std::make_pair(+0.900672, eps), // h_{7}
                std::make_pair(-0.056156, eps), // h_{8}
                std::make_pair(+0.200448, eps), // h_{9}
                std::make_pair(+0.106033, eps), // h_{10}
            };

            TEST_CHECK_DIAGNOSTICS(diag, ref);
        }
} bstar_to_dstar_hqet_form_factors_test;
