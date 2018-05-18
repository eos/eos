/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Ahmet Kokulu
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
#include <eos/b-decays/b-to-dstar-l-nu.hh>
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

class BToDstarLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToDstarLeptonNeutrinoTest() :
            TestCase("b_to_dstar_l_nu_test")
        {
        }

        virtual void run() const
        {

            // tests for normalized SM Br (|V_cb|=1)
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
                p1["B->D^*::alpha^A12_0@BSZ2015"] = +0.280646;
                p1["B->D^*::alpha^A12_1@BSZ2015"] = +0.72;
                p1["B->D^*::alpha^A12_2@BSZ2015"] = +1.33;
                p1["B->D^*::alpha^V_0@BSZ2015"  ] = +0.01;
                p1["B->D^*::alpha^V_1@BSZ2015"  ] = +0.02;
                p1["B->D^*::alpha^V_2@BSZ2015"  ] = +0.03;
                p1["B->D^*::alpha^T1_0@BSZ2015" ] = +0.27;
                p1["B->D^*::alpha^T1_1@BSZ2015" ] = -0.74;
                p1["B->D^*::alpha^T1_2@BSZ2015" ] = +1.45;
                p1["B->D^*::alpha^T2_0@BSZ2015" ] = +0.27;
                p1["B->D^*::alpha^T2_1@BSZ2015" ] = +0.47;
                p1["B->D^*::alpha^T2_2@BSZ2015" ] = +0.58;
                p1["B->D^*::alpha^T23_0@BSZ2015"] = +0.75;
                p1["B->D^*::alpha^T23_1@BSZ2015"] = +1.90;
                p1["B->D^*::alpha^T23_2@BSZ2015"] = +2.93;
                p1["mass::B_d"]                   = +5.279;
                p1["mass::D^*_d"]                 = +2.0103;
                p1["b->cmunumu::Re{cVL}"]         = +1.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BSZ2015");
                oo.set("q", "d");

                BToDstarLeptonNeutrino d(p1, oo);

                const double eps = 1e-5;

                // test for different s-bins - the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 8.0), 17.786, eps);
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(8.0, 10.68), 7.30794, eps);
                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(0.011164, 10.68), 44.3961, eps);
            }

            // tests for normalized NP Br (|V_cb|=1)
            {
                Parameters p2 = Parameters::Defaults();
                /*
                 * for the TEST case below the B->D^* SSE parameters are randomly chosen. However, the correlations together with the EOM conditions among the FFs are respected in this choice. Namely, alpha^A0_0 is correlated with alpha^A12_0, and also alpha^T1_0 should be the same as alpha^T2_0.
                 */
                p2["B->D^*::alpha^A0_0@BSZ2015" ] = +1.0;
                p2["B->D^*::alpha^A0_1@BSZ2015" ] = +0.24;
                p2["B->D^*::alpha^A0_2@BSZ2015" ] = +0.21;
                p2["B->D^*::alpha^A1_0@BSZ2015" ] = +0.5;
                p2["B->D^*::alpha^A1_1@BSZ2015" ] = +0.4;
                p2["B->D^*::alpha^A1_2@BSZ2015" ] = +0.3;
                p2["B->D^*::alpha^A12_0@BSZ2015"] = +0.280646;
                p2["B->D^*::alpha^A12_1@BSZ2015"] = +0.72;
                p2["B->D^*::alpha^A12_2@BSZ2015"] = +1.33;
                p2["B->D^*::alpha^V_0@BSZ2015"  ] = +0.01;
                p2["B->D^*::alpha^V_1@BSZ2015"  ] = +0.02;
                p2["B->D^*::alpha^V_2@BSZ2015"  ] = +0.03;
                p2["B->D^*::alpha^T1_0@BSZ2015" ] = +0.27;
                p2["B->D^*::alpha^T1_1@BSZ2015" ] = -0.74;
                p2["B->D^*::alpha^T1_2@BSZ2015" ] = +1.45;
                p2["B->D^*::alpha^T2_0@BSZ2015" ] = +0.27;
                p2["B->D^*::alpha^T2_1@BSZ2015" ] = +0.47;
                p2["B->D^*::alpha^T2_2@BSZ2015" ] = +0.58;
                p2["B->D^*::alpha^T23_0@BSZ2015"] = +0.75;
                p2["B->D^*::alpha^T23_1@BSZ2015"] = +1.90;
                p2["B->D^*::alpha^T23_2@BSZ2015"] = +2.93;
                p2["mass::B_d"]                   = +5.279;
                p2["mass::D^*_d"]                 = +2.0103;
                // fix scale
                p2["mu"]                          = +4.18;
                // mb(mb)
                p2["mass::b(MSbar)"]              = +4.18;
                // mc(mc)
                p2["mass::c"]                     = +1.275;
                p2["b->cmunumu::Re{cVL}"]         = +1.0;
                p2["b->cmunumu::Im{cVL}"]         = -1.0;
                p2["b->cmunumu::Re{cVR}"]         = +2.0;
                p2["b->cmunumu::Im{cVR}"]         = -2.0;
                p2["b->cmunumu::Re{cSL}"]         = +3.0;
                p2["b->cmunumu::Im{cSL}"]         = -3.0;
                p2["b->cmunumu::Re{cSR}"]         = +4.0;
                p2["b->cmunumu::Im{cSR}"]         = -4.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BSZ2015");

                BToDstarLeptonNeutrino d(p2, oo);

                const double eps = 1e-5;

                // test for different s-bins - the default lepton is muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(4.0, 8.0), 38.9976, eps);
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(8.0, 10.68), 15.3374, eps);
                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.normalized_integrated_branching_ratio(0.011164, 10.68), 94.2214, eps);
            }
        }
} b_to_dstar_l_nu_test;
