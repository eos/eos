/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <eos/utils/complex.hh>

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <iostream>

using namespace test;
using namespace eos;

class BToKstarGammaTest :
    public TestCase
{
    public:
        BToKstarGammaTest() :
            TestCase("b_to_kstar_gamma_test")
        {
        }

        virtual void run() const
        {
            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["Abs{c10}"] = 4.150;
                p["Arg{c10}"] = M_PI;
                p["Abs{c10'}"] = 0.000;
                p["Arg{c10'}"] = M_PI;
                // PDG 2010 CKM parameters
                p["CKM::A"] = 0.812;
                p["CKM::lambda"] = 0.22543;
                p["CKM::rhobar"] = 0.144;
                p["CKM::etabar"] = 0.342;
                p["decay-constant::B_d"] = 0.212;
                p["decay-constant::B_d"] = 0.256;
                p["mass::B_d"] = 5.2795;
                p["mass::B_s"] = 5.3663;
                p["life_time::B_d"] = 1.525e-12;
                p["life_time::B_s"] = 1.472e-12;
                p["life_time::Delta_B_d"] = 0.0;
                p["life_time::Delta_B_s"] = 0.104e12;

                const double eps = 1e-4;

                // B_d -> mu^+ mu^-
                {
                    Options oo;
                    oo.set("model", "WilsonScan");
                    oo.set("q", "d");
                    oo.set("l", "mu");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           +1.75327e-10, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), +1.75327e-10, eps);
                }

                // B_d -> e^+ e^-
                {
                    Options oo;
                    oo.set("model", "WilsonScan");
                    oo.set("q", "d");
                    oo.set("l", "e");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           +4.10420e-15, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), +4.10420e-15, eps);
                }

                // B_s -> mu^+ mu^-
                {
                    Options oo;
                    oo.set("model", "WilsonScan");
                    oo.set("q", "s");
                    oo.set("l", "mu");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           +3.03452e-09, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), +3.28604e-09, eps);
                }

                // B_s -> e^+ e^-
                {
                    Options oo;
                    oo.set("model", "WilsonScan");
                    oo.set("q", "s");
                    oo.set("l", "e");

                    BToDilepton d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_time_zero(),           +7.10333e-14, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_untagged_integrated(), +7.69211e-14, eps);
                }
            }
        }
} b_to_kstar_gamma_test;
