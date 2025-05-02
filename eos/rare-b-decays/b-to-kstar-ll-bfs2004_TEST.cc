/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
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
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class BToKstarDileptonBFS2004BobethCompatibilityTest :
    public TestCase
{
    public:
    BToKstarDileptonBFS2004BobethCompatibilityTest() :
        TestCase("b_to_kstar_dilepton_BFS2004_bobeth_compatibility_test")
    {
    }

    virtual void run() const
    {
        {
            Parameters p = Parameters::Defaults();
            p["CKM::abs(V_ub)"] =  0.003631275231633653;
            p["CKM::arg(V_ub)"] = -1.210765774253535;
            p["CKM::abs(V_cb)"] =  0.041996951916414726;
            p["CKM::arg(V_cb)"] =  0.0;
            p["CKM::abs(V_tb)"] =  0.9991111344469873;
            p["CKM::arg(V_tb)"] =  0.0;
            p["CKM::abs(V_us)"] =  0.22534851424944366;
            p["CKM::arg(V_us)"] =  0.0;
            p["CKM::abs(V_cs)"] =  0.9734061815416853;
            p["CKM::arg(V_cs)"] = -3.304199362533668e-05;
            p["CKM::abs(V_ts)"] =  0.04121212396309175;
            p["CKM::arg(V_ts)"] = -3.1230250224697222;
            p["b->s::c1"]      = -0.3231323312;
            p["b->s::c2"]      = 1.009301831;
            p["b->s::c3"]      = -0.005233499106;
            p["b->s::c4"]      = -0.08829686414;
            p["b->s::c5"]      = 0.0003601965805;
            p["b->s::c6"]      = 0.001020749573;
            p["sb::mu"]        = 4.2;
            p["b->s::Re{c7}"]  = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"]  = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"]      = -0.1827530948;
            p["sbmumu::mu"]         = 4.2;
            p["b->smumu::Re{c9}"]   = 4.294489364 + 1;
            p["b->smumu::Im{c9}"]   = 0.5;
            p["b->smumu::Re{c9'}"]  = 2;
            p["b->smumu::Im{c9'}"]  = 1.5;
            p["b->smumu::Re{c10}"]  = -4.196294696 + 3;
            p["b->smumu::Im{c10}"]  = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;
            p["K^*::a_1_para@1GeV"] = 0.1;
            p["K^*::a_1_perp@1GeV"] = 0.1;
            p["K^*::a_2_para@1GeV"] = 0.1;
            p["K^*::a_2_perp@1GeV"] = 0.1;
            p["B::1/lambda_B_p"] = 1.0 / 0.485;

            Options oo
            {
                {"model"_ok, "WET"},
                {"scan-mode"_ok, "cartesian"},
                {"tag"_ok, "BFS2004"},
                {"qcdf-integrals"_ok, "mixed"},
                {"form-factors"_ok, "KMPW2010"},
                {"l"_ok, "mu"},
                {"q"_ok, "d"}
            };

            static const double eps = 0.72e-2;
            static const double q2 = 6.0;

            BToKstarDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-1.120616135e-10, +6.005404351e-12), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>(+4.337275083e-11, +3.591794269e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-4.177379962e-11, +1.649925628e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>(+5.963768892e-11, +3.601537199e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>(+4.352686602e-11, +5.276889886e-12), 2.2*eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>(+9.326590159e-11, +1.116294121e-10), eps);
       }

       // scalar and tensor
       {
            // important to agree on alpha_s, can change values by 1%
            Parameters p = Parameters::Defaults();
            p["CKM::abs(V_ub)"] =  0.003631275231633653;
            p["CKM::arg(V_ub)"] = -1.210765774253535;
            p["CKM::abs(V_cb)"] =  0.041996951916414726;
            p["CKM::arg(V_cb)"] =  0.0;
            p["CKM::abs(V_tb)"] =  0.9991111344469873;
            p["CKM::arg(V_tb)"] =  0.0;
            p["CKM::abs(V_us)"] =  0.22534851424944366;
            p["CKM::arg(V_us)"] =  0.0;
            p["CKM::abs(V_cs)"] =  0.9734061815416853;
            p["CKM::arg(V_cs)"] = -3.304199362533668e-05;
            p["CKM::abs(V_ts)"] =  0.04121212396309175;
            p["CKM::arg(V_ts)"] = -3.1230250224697222;
            p["b->s::c1"]      = -0.3231323312;
            p["b->s::c2"]      = 1.009301831;
            p["b->s::c3"]      = -0.005233499106;
            p["b->s::c4"]      = -0.08829686414;
            p["b->s::c5"]      = 0.0003601965805;
            p["b->s::c6"]      = 0.001020749573;
            p["sb::mu"]        = 4.2;
            p["b->s::Re{c7}"]  = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"]  = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"]      = -0.1827530948;
            p["sbmumu::mu"]         = 4.2;
            p["b->smumu::Re{c9}"]   = 4.294489364 + 1;
            p["b->smumu::Im{c9}"]   = 0.5;
            p["b->smumu::Re{c9'}"]  = 2;
            p["b->smumu::Im{c9'}"]  = 1.5;
            p["b->smumu::Re{c10}"]  = -4.196294696 + 3;
            p["b->smumu::Im{c10}"]  = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;
            p["b->smumu::Re{cS}"]   = 0.5;
            p["b->smumu::Im{cS}"]   = 1;
            p["b->smumu::Re{cS'}"]  = 0.6;
            p["b->smumu::Im{cS'}"]  = 1.1;
            p["b->smumu::Re{cP}"]   = 0.7;
            p["b->smumu::Im{cP}"]   = 1.2;
            p["b->smumu::Re{cP'}"]  = 0.8;
            p["b->smumu::Im{cP'}"]  = 1.3;
            p["b->smumu::Re{cT}"]   = 0.9;
            p["b->smumu::Im{cT}"]   = 1.4;
            p["b->smumu::Re{cT5}"]  = 1.0;
            p["b->smumu::Im{cT5}"]  = 1.5;

            p["mass::s(2GeV)"] = 0.12;

            Options oo
            {
                {"model"_ok, "WET"},
                {"scan-mode"_ok, "cartesian"},
                {"tag"_ok, "BFS2004"},
                {"form-factors"_ok, "KMPW2010"},
                {"l"_ok, "mu"},
                {"q"_ok, "u"}
            };

            BToKstarDilepton d(p, oo);

            static const double q2 = 6.0;
            auto amps = d.amplitudes(q2);

            double eps = 3e-2;
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-1.121022032e-10, 5.991646324e-12), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 4.333216107e-11, 3.590418466e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.353425835e-11, 5.287884397e-12), eps);

            eps = 1.2e-2;
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 9.326263994e-11, 1.117078741e-10), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-4.176304959e-11, 1.651237347e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 5.964843896e-11, 3.602848917e-11), eps);

            // nearly identically implemented, only difference from alpha_s
            eps = 2e-4;
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-2.247078271e-10, -6.370327589e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.185643583e-12,  2.185643583e-12), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 1.47786e-11,      2.2989e-11),      eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>(-1.64207e-11,     -2.46311e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 2.4322e-11,       3.78342e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>(-2.70244e-11,     -4.05366e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>(-3.24769e-11,     -4.87154e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.92292e-11,      4.54677e-11),     eps);
       }
    }
} b_to_kstar_dilepton_BFS2004_bobeth_compatibility_test;
