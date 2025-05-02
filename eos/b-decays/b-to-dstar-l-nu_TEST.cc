/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2019 Ahmet Kokulu
 * Copyright (c) 2019-2025 Danny van Dyk
 * Copyright (c) 2021      Christoph Bobeth
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
#include <eos/maths/complex.hh>
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

                Options o
                {
                    { "V"_ok,                  "D^*"       },
                    { "q"_ok,                  "d"         },
                    { "l"_ok,                  "e"         },
                    { "model"_ok,              "CKM"       },
                    { "z-order-lp"_ok,         "3"         },
                    { "z-order-slp"_ok,        "2"         },
                    { "z-order-sslp"_ok,       "1"         },
                    { "form-factors"_ok,       "BGJvD2019" },
                    { "integration-points"_ok, "4096"      }
                };

                Kinematics k
                {
                    { "q2_min", 0.001 }, { "q2_max", 10.689 },
                };

                const double eps = 1e-3;
                // Christoph Bobeth: Adjusted test case because increased number of integration points
                //                   in numerical integration from 256 -> 4096
                BToVectorLeptonNeutrino d(p, o);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(0.001, 10.689), 33.3247, eps);
                auto ir = d.prepare(0.001, 10.689);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_f_L(ir), 0.546, eps);

                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_1c", p, k, o)->evaluate(),  0.409302220, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_1s", p, k, o)->evaluate(),  0.255523335, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_2c", p, k, o)->evaluate(), -0.409302220, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_2s", p, k, o)->evaluate(),  0.085174445, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_3",  p, k, o)->evaluate(), -0.134468151, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_4",  p, k, o)->evaluate(),  0.231808464, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_5",  p, k, o)->evaluate(),  0.165381861, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_6c", p, k, o)->evaluate(),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_6s", p, k, o)->evaluate(),  0.200153929, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_7",  p, k, o)->evaluate(),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_8",  p, k, o)->evaluate(),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_9",  p, k, o)->evaluate(),  0.0,         eps);
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

                Options o
                {
                    { "V"_ok,                  "D^*"       },
                    { "q"_ok,                  "d"         },
                    { "l"_ok,                  "tau"       },
                    { "model"_ok,              "CKM"       },
                    { "z-order-lp"_ok,         "3"         },
                    { "z-order-slp"_ok,        "2"         },
                    { "z-order-sslp"_ok,       "1"         },
                    { "form-factors"_ok,       "BGJvD2019" },
                    { "integration-points"_ok, "4096"      }
                };

                Kinematics k
                {
                    { "q2_min", 3.157 }, { "q2_max", 10.689 },
                };

                const double eps = 1e-3;

                BToVectorLeptonNeutrino d(p, o);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(3.157, 10.689), 8.213, eps);
                auto ir = d.prepare(3.157, 10.689);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_f_L(ir), 0.475, eps);

                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_1c", p, k, o)->evaluate(),  0.4325856250, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_1s", p, k, o)->evaluate(),  0.2779590234, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_2c", p, k, o)->evaluate(), -0.1287773345, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_2s", p, k, o)->evaluate(),  0.0471441750, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_3",  p, k, o)->evaluate(), -0.0819412032, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_4",  p, k, o)->evaluate(),  0.1057578408, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_5",  p, k, o)->evaluate(),  0.2056068494, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_6c", p, k, o)->evaluate(), -0.2766922602, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_6s", p, k, o)->evaluate(),  0.1598442669, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_7",  p, k, o)->evaluate(),  0.0,          eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_8",  p, k, o)->evaluate(),  0.0,          eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_9",  p, k, o)->evaluate(),  0.0,          eps);
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

                Options o
                {
                    { "V"_ok,                  "D^*"       },
                    { "q"_ok,                  "d"         },
                    { "l"_ok,                  "mu"        },
                    { "model"_ok,              "WET"       },
                    { "z-order-lp"_ok,         "3"         },
                    { "z-order-slp"_ok,        "2"         },
                    { "z-order-sslp"_ok,       "1"         },
                    { "form-factors"_ok,       "BGJvD2019" },
                    { "integration-points"_ok, "4096"      }
                };

                Kinematics k
                {
                    { "q2_min", p["mass::mu"]* p["mass::mu"] }, { "q2_max", 10.689 },
                };

                const double eps = 1e-3;

                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_1c", p, k, o)->evaluate(),  0.362439,  eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_1s", p, k, o)->evaluate(),  0.268109,  eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_2c", p, k, o)->evaluate(), -0.228862,  eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_2s", p, k, o)->evaluate(), -0.0375834, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_3",  p, k, o)->evaluate(), -0.0600368, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_4",  p, k, o)->evaluate(),  0.0897816, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_5",  p, k, o)->evaluate(),  0.0837827, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_6c", p, k, o)->evaluate(), -0.0716409, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::S_6s", p, k, o)->evaluate(),  0.0433597, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::A_7",  p, k, o)->evaluate(),  0.0205058, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::A_8",  p, k, o)->evaluate(), -0.0113015, eps);
                TEST_CHECK_NEARLY_EQUAL(Observable::make("B->D^*lnu::A_9",  p, k, o)->evaluate(),  0.013735,  eps);
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

                Options oo
                {
                    { "V"_ok,                  "D^*"     },
                    { "q"_ok,                  "d"       },
                    { "model"_ok,              "WET"     },
                    { "form-factors"_ok,       "BSZ2015" },
                    { "integration-points"_ok, "4096"    }
                };

                BToVectorLeptonNeutrino d(p1, oo);

                const double eps = 1e-3;

                auto ir = d.prepare(4.0, 10.68);

                // the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 10.68), 25.4230, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(ir), 0.000494949, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_f_L(ir),     0.737489, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_1(ir),  -0.130926, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_2(ir),   0.00266046, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_3(ir),   0.230111, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_1(ir), 0.0, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_2(ir), 0.0, eps);
                //TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_3(ir), 0.0, eps);

                Kinematics k
                {
                    { "q2_mu_min",   4.00 }, { "q2_mu_max",  10.68 },
                    { "q2_tau_min",  4.00 }, { "q2_tau_max", 10.68 },
                };
                auto obs_RDst = Observable::make("B->D^*lnu::R_D^*", p1, k, oo);
                TEST_CHECK_RELATIVE_ERROR(obs_RDst->evaluate(), 0.379092, eps);
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

                Options oo
                {
                    { "V"_ok,                  "D^*"     },
                    { "q"_ok,                  "d"       },
                    { "model"_ok,              "WET"     },
                    { "form-factors"_ok,       "BSZ2015" },
                    { "integration-points"_ok, "4096"    }
                };

                BToVectorLeptonNeutrino d(p3, oo);

                const double eps = 1e-3;

                // the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 10.68), 3431.13, eps);
                auto ir = d.prepare(4.0, 10.68);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(ir), 0.0409932, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_f_L(ir),    0.50729, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_1(ir),  0.184031, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_2(ir), -0.0282197, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_c_3(ir), -0.42545, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_1(ir),  0.0000348895, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_2(ir),  0.000268975, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_t_3(ir), -0.0000320251, eps);

                Kinematics k
                {
                    { "q2_mu_min",   4.00 }, { "q2_mu_max",  10.68 },
                    { "q2_tau_min",  4.00 }, { "q2_tau_max", 10.68 },
                };
                auto obs_RDst = Observable::make("B->D^*lnu::R_D^*", p3, k, oo);
                TEST_CHECK_RELATIVE_ERROR(obs_RDst->evaluate(), 1.20331, eps);
            }

            // Check of consistency
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_cb)"].set(1.0);
                p["B(*)->D(*)::a@HQET"].set(1.000);
                p["mass::B_d"].set(5.27942);
                p["mass::D_u^*"].set(2.01000);
                p["mass::D_d^*"].set(2.01000);
                p["life_time::B_d"].set(1.520e-12);

                Options o
                {
                    { "V"_ok,             "D^*"       },
                    { "q"_ok,             "d"         },
                    { "model"_ok,         "CKM"       },
                    { "z-order-lp"_ok,    "3"         },
                    { "z-order-slp"_ok,   "2"         },
                    { "z-order-sslp"_ok,  "1"         },
                    { "form-factors"_ok,  "BGJvD2019" },
                };
                Kinematics k
                {
                    { "q2_mu_min", 1.00 }, { "q2_mu_max", 10.68 },
                    { "q2_e_min",  1.00 }, { "q2_e_max",  10.68 },
                };

                auto obs_BRbar    = Observable::make("B->D^*lnu::BRbar",        p, k, o);
                auto obs_deltaBR  = Observable::make("B->D^*lnu::DeltaBR",      p, k, o);

                k =
                {
                    { "q2_min",  1.00 }, { "q2_max",  10.68 },
                };
                auto obs_e_BR  = Observable::make("B->D^*lnu::BR;l=e",    p, k, o);
                auto obs_mu_BR = Observable::make("B->D^*lnu::BR;l=mu",   p, k, o);

                const double eps = 1e-5;
                TEST_CHECK_RELATIVE_ERROR(
                    0.5 * (obs_e_BR->evaluate() + obs_mu_BR->evaluate()),
                    obs_BRbar->evaluate(),
                    eps
                );
                TEST_CHECK_RELATIVE_ERROR(
                    obs_mu_BR->evaluate() - obs_e_BR->evaluate(),
                    obs_deltaBR->evaluate(),
                    eps
                );
            }

        }
} b_to_dstar_l_nu_test;
