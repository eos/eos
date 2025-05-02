/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/rare-b-decays/bs-to-phi-charmonium.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class BsToPhiCharmoniumGvDV2020Test :
    public TestCase
{
    public:
    BsToPhiCharmoniumGvDV2020Test() :
            TestCase("bs_to_phi_charmonium_GvDV2020_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_s"]                             = 5.366;
            p["mass::phi"]                             = 1.020;
            p["mass::J/psi"]                           = 3.0969;
            p["mass::psi(2S)"]                         = 3.6860;
            p["mass::D^0"]                             = 1.86723;
            p["b->sccbar::t_0"]                        = 4.0;
            p["b->sccbar::t_s"]                        = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]            = 1.81e-4;

            p["B_s->phiccbar::Re{alpha_0^perp}@GvDV2020"]  = 0.02;
            p["B_s->phiccbar::Im{alpha_0^perp}@GvDV2020"]  = 0.03;
            p["B_s->phiccbar::Re{alpha_1^perp}@GvDV2020"]  = 0.04;
            p["B_s->phiccbar::Im{alpha_1^perp}@GvDV2020"]  = 0.05;
            p["B_s->phiccbar::Re{alpha_2^perp}@GvDV2020"]  = 0.06;
            p["B_s->phiccbar::Im{alpha_2^perp}@GvDV2020"]  = 0.07;
            p["B_s->phiccbar::Re{alpha_0^para}@GvDV2020"]  = 0.08;
            p["B_s->phiccbar::Im{alpha_0^para}@GvDV2020"]  = 0.09;
            p["B_s->phiccbar::Re{alpha_1^para}@GvDV2020"]  = 0.010;
            p["B_s->phiccbar::Im{alpha_1^para}@GvDV2020"]  = 0.011;
            p["B_s->phiccbar::Re{alpha_2^para}@GvDV2020"]  = 0.012;
            p["B_s->phiccbar::Im{alpha_2^para}@GvDV2020"]  = 0.013;
            p["B_s->phiccbar::Re{alpha_0^long}@GvDV2020"]  = 0.014;
            p["B_s->phiccbar::Im{alpha_0^long}@GvDV2020"]  = 0.015;
            p["B_s->phiccbar::Re{alpha_1^long}@GvDV2020"]  = 0.016;
            p["B_s->phiccbar::Im{alpha_1^long}@GvDV2020"]  = 0.017;
            p["B_s->phiccbar::Re{alpha_2^long}@GvDV2020"]  = 0.018;
            p["B_s->phiccbar::Im{alpha_2^long}@GvDV2020"]  = 0.019;

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

            Options oo
            {
                {"model"_ok,               "WET"},
                {"nonlocal-formfactor"_ok, "GvDV2020"},
                {"psi"_ok,                 "J/psi"},
                {"q"_ok,                   "s"}
            };

            BsToPhiCharmonium c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.branching_ratio(),  104.52058, 1e-5);

        }
} bs_to_phi_charmonium_GvDV2020_test;
