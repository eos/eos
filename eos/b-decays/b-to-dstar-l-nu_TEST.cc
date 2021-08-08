/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2019, 2021 Danny van Dyk
 * Copyright (c) 2021 Christoph Bobeth
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
#include <eos/observable.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class BToVectorLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToVectorLeptonNeutrinoTest() :
            TestCase("b_to_dstar_l_nu_test")
        {
        }

        virtual void run() const
        {
            // comparison with Martin Jung in 3/2/1 model
            // l = electron
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["CKM::abs(V_cb)"].set(1.0);
                p["mass::e"].set(0.000001);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);
                p["life_time::B_d"].set(1.520e-12);

                Options o{
                    { "l",             "e"       },
                    { "model",         "CKMScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    },
                    { "integration-points",    "4096"    }
                };

                BToVectorLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                // Christoph Bobeth: Adjusted test case because increased number of integration points
                //                   in numerical integration from 256 -> 4096
                TEST_CHECK_NEARLY_EQUAL( 33.3247,     d.integrated_branching_ratio(0.001, 10.689), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.546,       d.integrated_f_L(0.001, 10.689),  eps);
                auto ir = d.prepare(0.001, 10.689);
                TEST_CHECK_NEARLY_EQUAL( 0.409302220, d.integrated_S1c(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.255523335, d.integrated_S1s(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.409302220, d.integrated_S2c(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.085174445, d.integrated_S2s(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.134468151, d.integrated_S3 (0.001, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.231808464, d.integrated_S4 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.165381861, d.integrated_S5 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,         d.integrated_S6c(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.200153929, d.integrated_S6s(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,         d.integrated_S7 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,         d.integrated_S8 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,         d.integrated_S9 (ir),  eps);
            }

            // comparison with Martin Jung in 3/2/1 model
            // l = tau
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.06919);
                p["B(*)->D(*)::xi''(1)@HQET"].set(1.66581);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(-2.91356);
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0600518);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0000101998);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.085385);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0400584);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.0522346);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.586099);
                p["B(*)->D(*)::eta'(1)@HQET"].set(-0.0233426);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.0288193);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.113962);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-8.15957);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.90706);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.16227);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-3.41861);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(5.6966);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-1.89465);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(0.220492);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.97017);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.34767);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.98608);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.42747);

                p["CKM::abs(V_cb)"].set(1.0);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);
                p["life_time::B_d"].set(1.520e-12);

                Options o{
                    { "l",             "tau"     },
                    { "model",         "CKMScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    },
                    { "integration-points",    "4096"    }
                };

                BToVectorLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                TEST_CHECK_NEARLY_EQUAL( 8.213,        d.integrated_branching_ratio(3.157, 10.689), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.475,        d.integrated_f_L(3.157, 10.689),  eps);
                auto ir = d.prepare(3.157, 10.689);
                TEST_CHECK_NEARLY_EQUAL( 0.4325856250, d.integrated_S1c(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2779590234, d.integrated_S1s(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1287773345, d.integrated_S2c(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0471441750, d.integrated_S2s(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0819412032, d.integrated_S3 (3.157, 10.689),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1057578408, d.integrated_S4 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2056068494, d.integrated_S5 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.2766922602, d.integrated_S6c(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1598442669, d.integrated_S6s(ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,          d.integrated_S7 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,          d.integrated_S8 (ir),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,          d.integrated_S9 (ir),  eps);
            }

            // New physics comparison with Martin Jung:
            {
                Parameters p = Parameters::Defaults();

                p["B(*)->D(*)::xi'(1)@HQET"].set(-1.10422);
                p["B(*)->D(*)::xi''(1)@HQET"].set(2* 0.912531);
                p["B(*)->D(*)::xi'''(1)@HQET"].set(6* (-0.565251));
                p["B(*)->D(*)::chi_2(1)@HQET"].set(-0.0648414);
                p["B(*)->D(*)::chi_2'(1)@HQET"].set(-0.0138642);
                p["B(*)->D(*)::chi_2''(1)@HQET"].set(-0.0850267);
                p["B(*)->D(*)::chi_3'(1)@HQET"].set(0.0380865);
                p["B(*)->D(*)::chi_3''(1)@HQET"].set(-0.104527);
                p["B(*)->D(*)::eta(1)@HQET"].set(0.626659);
                p["B(*)->D(*)::eta'(1)@HQET"].set(0.0402494);
                p["B(*)->D(*)::eta''(1)@HQET"].set(-0.338511);
                p["B(*)->D(*)::l_1(1)@HQET"].set(0.178431);
                p["B(*)->D(*)::l_1'(1)@HQET"].set(-6.99412);
                p["B(*)->D(*)::l_2(1)@HQET"].set(-1.58674);
                p["B(*)->D(*)::l_2'(1)@HQET"].set(-3.32927);
                p["B(*)->D(*)::l_3(1)@HQET"].set(-4.0127);
                p["B(*)->D(*)::l_3'(1)@HQET"].set(6.64143);
                p["B(*)->D(*)::l_4(1)@HQET"].set(-2.22468);
                p["B(*)->D(*)::l_4'(1)@HQET"].set(-0.607232);
                p["B(*)->D(*)::l_5(1)@HQET"].set(4.79295);
                p["B(*)->D(*)::l_5'(1)@HQET"].set(-2.06147);
                p["B(*)->D(*)::l_6(1)@HQET"].set(1.95851);
                p["B(*)->D(*)::l_6'(1)@HQET"].set(1.22043);

                p["CKM::abs(V_cb)"]            =  0.041996951916414726;
                p["mass::e"].set(0.00000000000001);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);
                p["life_time::B_d"].set(1.520e-12);
 
                p["cbmunumu::Re{cVL}"].set(+1.2);
                p["cbmunumu::Im{cVL}"].set(+0.0);
                p["cbmunumu::Re{cVR}"].set(-0.4*1.2);
                p["cbmunumu::Im{cVR}"].set(+0.2*1.2);
                p["cbmunumu::Re{cSL}"].set(+0.2*1.2);
                p["cbmunumu::Im{cSL}"].set(+0.6*1.2);
                p["cbmunumu::Re{cSR}"].set(+0.8*1.2);
                p["cbmunumu::Im{cSR}"].set(+0.3*1.2);
                p["cbmunumu::Re{cT}"].set(-0.1*1.2);
                p["cbmunumu::Im{cT}"].set(+0.2*1.2);

                p["cbmunumu::mu"].set(2.295);  // to get m_b,c in amp_P comparable to Martin

                Options o{
                    { "l",             "mu"      },
                    { "model",         "WilsonScan" },
                    { "q",             "d"       },
                    { "z-order-lp",    "3"       },
                    { "z-order-slp",   "2"       },
                    { "z-order-sslp",  "1"       },
                    { "form-factors",  "HQET"    },
                    { "integration-points",    "4096"    }
                };

                BToVectorLeptonNeutrino d(p, o);

                const double eps = 1e-3;
                const double q2min = p["mass::mu"]* p["mass::mu"];
                const double q2max = 10.689;
                // 
                auto ir = d.prepare(q2min, q2max);
                TEST_CHECK_NEARLY_EQUAL( 0.362439,  d.integrated_S1c(ir), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.268109,  d.integrated_S1s(ir), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.228862,  d.integrated_S2c(ir), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0375834, d.integrated_S2s(ir), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0600368, d.integrated_S3 (q2min, q2max), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0897816, d.integrated_S4 (ir), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0837827, d.integrated_S5 (ir), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0716409, d.integrated_S6c(ir), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0433597, d.integrated_S6s(ir), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0205058, d.integrated_A7 (ir), eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0113015, d.integrated_A8 (ir), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.013735,  d.integrated_A9 (ir), eps);
            }

            // SM tests cf. [DSD2014]
            {
                Parameters p1 = Parameters::Defaults();
                /*
                 * for the TEST case below the B->D^* SSE parameters are randomly chosen. However, the correlations together with the EOM conditions among the FFs are respected in this choice. Namely, alpha^A0_0 is correlated with alpha^A12_0, and also alpha^T1_0 should be the same as alpha^T2_0.
                 */
                p1["B->D^*::alpha^A0_0@BSZ2015" ] = +1.0;
                p1["B->D^*::alpha^A0_1@BSZ2015" ] = +0.24;
                p1["B->D^*::alpha^A0_2@BSZ2015" ] = +0.21;
                p1["B->D^*::alpha^A1_0@BSZ2015" ] = +0.5;
                p1["B->D^*::alpha^A1_1@BSZ2015" ] = +0.4;
                p1["B->D^*::alpha^A1_2@BSZ2015" ] = +0.3;
                p1["B->D^*::alpha^A12_1@BSZ2015"] = +0.72;
                p1["B->D^*::alpha^A12_2@BSZ2015"] = +1.33;
                p1["B->D^*::alpha^V_0@BSZ2015"  ] = +0.01;
                p1["B->D^*::alpha^V_1@BSZ2015"  ] = +0.02;
                p1["B->D^*::alpha^V_2@BSZ2015"  ] = +0.03;
                p1["B->D^*::alpha^T1_0@BSZ2015" ] = +0.27;
                p1["B->D^*::alpha^T1_1@BSZ2015" ] = -0.74;
                p1["B->D^*::alpha^T1_2@BSZ2015" ] = +1.45;
                p1["B->D^*::alpha^T2_1@BSZ2015" ] = +0.47;
                p1["B->D^*::alpha^T2_2@BSZ2015" ] = +0.58;
                p1["B->D^*::alpha^T23_0@BSZ2015"] = +0.75;
                p1["B->D^*::alpha^T23_1@BSZ2015"] = +1.90;
                p1["B->D^*::alpha^T23_2@BSZ2015"] = +2.93;
                p1["mass::B_d"]                   = +5.279;
                p1["mass::D_d^*"]                 = +2.0103;
                // by default, all other couplings are zero in eos
                p1["CKM::abs(V_cb)"]            =  0.041996951916414726;
                p1["cbmunumu::Re{cVL}"]         = +1.0066;  // include Sirlin correction
                p1["cbtaunutau::Re{cVL}"]       = +1.0066;  // include Sirlin correction

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BSZ2015");
                oo.set("q", "d");
                oo.set("integration-points", "4096");

                BToVectorLeptonNeutrino d(p1, oo);

                const double eps = 1e-3;

                // the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 10.68), 25.4230, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(4.0, 10.68), 0.000494949, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_f_L(4.0, 10.68), 0.737489, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_1(4.0, 10.68), -0.130926, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_2(4.0, 10.68), 0.00266046, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_3(4.0, 10.68), 0.230111, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_1(4.0, 10.68), 0.0, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_2(4.0, 10.68), 0.0, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_3(4.0, 10.68), 0.0, eps);

                Kinematics k
                {
                    { "q2_mu_min",   4.00 }, { "q2_mu_max",  10.68 },
                    { "q2_tau_min",  4.00 }, { "q2_tau_max", 10.68 },
                };
                auto obs_RDst = Observable::make("B->D^*lnu::R_D^*", p1, k, oo);
                TEST_CHECK_RELATIVE_ERROR(0.379092, obs_RDst->evaluate(), eps);
            }

            // NP tests cf. [DSD2014]
            {
                const double etaEW = 1.0066;

                Parameters p3 = Parameters::Defaults();
                /*
                 * for the TEST case below the B->D^* SSE parameters are randomly chosen.
                 * However, the correlations together with the EOM conditions among the FFs
                 * are respected in this choice. Namely, alpha^A0_0 is correlated with alpha^A12_0,
                 * and also alpha^T1_0 should be the same as alpha^T2_0.
                 */
                p3["B->D^*::alpha^A0_0@BSZ2015" ] = +1.0;
                p3["B->D^*::alpha^A0_1@BSZ2015" ] = +0.24;
                p3["B->D^*::alpha^A0_2@BSZ2015" ] = +0.21;
                p3["B->D^*::alpha^A1_0@BSZ2015" ] = +0.5;
                p3["B->D^*::alpha^A1_1@BSZ2015" ] = +0.4;
                p3["B->D^*::alpha^A1_2@BSZ2015" ] = +0.3;
                p3["B->D^*::alpha^A12_1@BSZ2015"] = +0.72;
                p3["B->D^*::alpha^A12_2@BSZ2015"] = +1.33;
                p3["B->D^*::alpha^V_0@BSZ2015"  ] = +0.01;
                p3["B->D^*::alpha^V_1@BSZ2015"  ] = +0.02;
                p3["B->D^*::alpha^V_2@BSZ2015"  ] = +0.03;
                p3["B->D^*::alpha^T1_0@BSZ2015" ] = +0.27;
                p3["B->D^*::alpha^T1_1@BSZ2015" ] = -0.74;
                p3["B->D^*::alpha^T1_2@BSZ2015" ] = +1.45;
                p3["B->D^*::alpha^T2_1@BSZ2015" ] = +0.47;
                p3["B->D^*::alpha^T2_2@BSZ2015" ] = +0.58;
                p3["B->D^*::alpha^T23_0@BSZ2015"] = +0.75;
                p3["B->D^*::alpha^T23_1@BSZ2015"] = +1.90;
                p3["B->D^*::alpha^T23_2@BSZ2015"] = +2.93;
                p3["mass::B_d"]                   = +5.279;
                p3["mass::D_d^*"]                 = +2.0103;
                // fix scale
                p3["cbmunumu::mu"]                = +4.18;
                p3["cbtaunutau::mu"]              = +4.18;
                // mb(mb)
                p3["mass::b(MSbar)"]              = +4.18;
                // mc(mc)
                p3["mass::c"]                     = +1.275;
                // CKM
                p3["CKM::abs(V_cb)"]              =  0.041996951916414726;
                // mu mode
                p3["cbmunumu::Re{cVL}"]         = +1.0 * etaEW;
                p3["cbmunumu::Im{cVL}"]         = -2.0 * etaEW;
                p3["cbmunumu::Re{cVR}"]         = +2.0 * etaEW;
                p3["cbmunumu::Im{cVR}"]         = -2.0 * etaEW;
                p3["cbmunumu::Re{cSL}"]         = +3.0 * etaEW;
                p3["cbmunumu::Im{cSL}"]         = -3.0 * etaEW;
                p3["cbmunumu::Re{cSR}"]         = +4.0 * etaEW;
                p3["cbmunumu::Im{cSR}"]         = -4.0 * etaEW;
                p3["cbmunumu::Re{cT}"]          = +5.0 * etaEW;
                p3["cbmunumu::Im{cT}"]          = -5.0 * etaEW;
                // tau mode
                p3["cbtaunutau::Re{cVL}"]       = +1.0 * etaEW;
                p3["cbtaunutau::Im{cVL}"]       = -5.0 * etaEW;
                p3["cbtaunutau::Re{cVR}"]       = +2.1 * etaEW;
                p3["cbtaunutau::Im{cVR}"]       = -6.0 * etaEW;
                p3["cbtaunutau::Re{cSL}"]       = +3.1 * etaEW;
                p3["cbtaunutau::Im{cSL}"]       = -7.0 * etaEW;
                p3["cbtaunutau::Re{cSR}"]       = +4.1 * etaEW;
                p3["cbtaunutau::Im{cSR}"]       = -8.0 * etaEW;
                p3["cbtaunutau::Re{cT}"]        = +5.1 * etaEW;
                p3["cbtaunutau::Im{cT}"]        = -9.0 * etaEW;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BSZ2015");
                oo.set("integration-points", "4096");

                BToVectorLeptonNeutrino d(p3, oo);

                const double eps = 1e-3;

                // the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 10.68), 3431.13, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(4.0, 10.68), 0.0409932, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_f_L(4.0, 10.68), 0.50729, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_1(4.0, 10.68), 0.184031, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_2(4.0, 10.68), -0.0282197, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_3(4.0, 10.68), -0.42545, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_1(4.0, 10.68), 0.0000348895, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_2(4.0, 10.68), 0.000268975, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_3(4.0, 10.68), -0.0000320251, eps);

                Kinematics k
                {
                    { "q2_mu_min",   4.00 }, { "q2_mu_max",  10.68 },
                    { "q2_tau_min",  4.00 }, { "q2_tau_max", 10.68 },
                };
                auto obs_RDst = Observable::make("B->D^*lnu::R_D^*", p3, k, oo);
                TEST_CHECK_RELATIVE_ERROR(1.20331, obs_RDst->evaluate(), eps);
            }
        }
} b_to_dstar_l_nu_test;
