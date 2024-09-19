/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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
#include <eos/maths/complex.hh>
#include <eos/rare-b-decays/b-to-ll.hh>

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <iostream>

using namespace test;
using namespace eos;

class BToDileptonTest :
    public TestCase
{
    public:
        BToDileptonTest() :
            TestCase("b_to_dilepton_test")
        {
        }

        virtual void run() const
        {
            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["sbmumu::mu"] = 4.2;
                p["sbee::mu"]   = 4.2;
                p["b->smumu::Re{c10}"]  = -4.150;
                p["b->smumu::Re{c10'}"] =  0.000;
                p["b->see::Re{c10}"]    = -4.150;
                p["b->see::Re{c10'}"]   =  0.000;
                // PDG 2010 CKM parameters
                p["CKM::A"]             =  0.812;
                p["CKM::lambda"]        =  0.22543;
                p["CKM::rhobar"]        =  0.144;
                p["CKM::etabar"]        =  0.342;
                p["CKM::abs(V_ub)"]     =  0.003540950873054711;
                p["CKM::arg(V_ub)"]     = -1.1728563751359748;
                p["CKM::abs(V_cb)"]     =  0.04126451344307112;
                p["CKM::arg(V_cb)"]     =  0.0;
                p["CKM::abs(V_tb)"]     =  0.9991419776905534;
                p["CKM::arg(V_tb)"]     =  0.0;
                p["CKM::abs(V_td)"]     =  0.008576901910577167;
                p["CKM::arg(V_td)"]     = -0.37951557931964897;
                p["CKM::abs(V_us)"]     =  0.22542858674178629;
                p["CKM::arg(V_us)"]     =  0.0;
                p["CKM::abs(V_cs)"]     =  0.9734167680132911;
                p["CKM::arg(V_cs)"]     = -3.119448393424795e-05;
                p["CKM::abs(V_ts)"]     =  0.04051834255894421;
                p["CKM::arg(V_ts)"]     = -3.123445879630718;

                p["decay-constant::B_s"] = 0.2276;
                p["mass::B_s"] = 5.3663;
                p["life_time::B_s"] = 1.472e-12;
                p["life_time::Delta_B_s"] = 0.104e12;

                static const double eps = 1e-4;

                // B_s -> mu^+ mu^-
                {
                    Options oo;
                    oo.declare("model", "WET");
                    oo.declare("q", "s");
                    oo.declare("l", "mu");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           +3.03452e-09, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), +3.28604e-09, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.cp_asymmetry_del_gamma(),              +1, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.cp_asymmetry_mixing_S(),                  0, eps);
                }

                // B_s -> e^+ e^-
                {
                    Options oo;
                    oo.declare("model", "WET");
                    oo.declare("q", "s");
                    oo.declare("l", "e");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           +7.10333e-14, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), +7.69211e-14, eps);
                }
            }

            // new physics
            {
                Parameters p = Parameters::Defaults();

                // test large NP contributions
                p["sbmumu::mu"] = 4.2;
                p["b->smumu::Re{c10}"] = -4.196294696 + 3;
                p["b->smumu::Im{c10}"] = 2.5;
                p["b->smumu::Re{c10'}"] = 4;
                p["b->smumu::Im{c10'}"] = 3.5;
                p["b->smumu::Re{cS}"] = 0.5;
                p["b->smumu::Im{cS}"] = 1;
                p["b->smumu::Re{cS'}"] = 0.6;
                p["b->smumu::Im{cS'}"] = 1.1;
                p["b->smumu::Re{cP}"] = 0.7;
                p["b->smumu::Im{cP}"] = 1.2;
                p["b->smumu::Re{cP'}"] = 0.8;
                p["b->smumu::Im{cP'}"] = 1.3;

                // 2013 default values
                p["CKM::A"]          = +0.827;
                p["CKM::lambda"]     =  0.22535;
                p["CKM::rhobar"]     =  0.132;
                p["CKM::etabar"]     =  0.350;
                p["CKM::abs(V_ub)"]  =  0.003631275231633653;
                p["CKM::arg(V_ub)"]  = -1.210765774253535;
                p["CKM::abs(V_cb)"]  =  0.041996951916414726;
                p["CKM::arg(V_cb)"]  =  0.0;
                p["CKM::abs(V_tb)"]  =  0.9991111344469873;
                p["CKM::arg(V_tb)"]  =  0.0;
                p["CKM::abs(V_us)"]  =  0.22534851424944366;
                p["CKM::arg(V_us)"]  =  0.0;
                p["CKM::abs(V_cs)"]  =  0.9734061815416853;
                p["CKM::arg(V_cs)"]  = -3.304199362533668e-05;
                p["CKM::abs(V_ts)"]  =  0.04121212396309175;
                p["CKM::arg(V_ts)"]  = -3.1230250224697222;

                p["mass::B_s"] = 5.36677;
                p["life_time::B_s"] = 1.516e-12;
                p["life_time::Delta_B_s"] = 0.081e12;
                p["decay-constant::B_s"] = 0.2276;

                static const double eps = 1e-4;

                // B_s -> mu^+ mu^-
                {
                    Options oo;
                    oo.declare("model", "WET");
                    oo.declare("scan-mode", "cartesian");
                    oo.declare("q", "s");
                    oo.declare("l", "mu");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           2.030257955e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), 2.098985874e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.cp_asymmetry_del_gamma(),              0.4878740356, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.cp_asymmetry_mixing_S(),               0.4617576325, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.effective_lifetime(),                  2.387625253e+12, eps);
                }
            }
        }
} b_to_dilepton_test;
