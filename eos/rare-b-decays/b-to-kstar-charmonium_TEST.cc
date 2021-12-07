/*
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/rare-b-decays/b-to-kstar-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class BToKstarCharmoniumGvDV2020Test :
    public TestCase
{
    public:
    BToKstarCharmoniumGvDV2020Test() :
            TestCase("b_to_kstar_charmonium_GvDV2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]                               = 5.279;
            p["mass::K_d^*"]                             = 0.896;
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 4.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;

            p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 0.02;
            p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 0.03;
            p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 0.04;
            p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 0.05;
            p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 0.06;
            p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 0.07;
            p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 0.08;
            p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 0.09;
            p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 0.010;
            p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 0.011;
            p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 0.012;
            p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 0.013;
            p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 0.014;
            p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 0.015;
            p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 0.016;
            p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 0.017;
            p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 0.018;
            p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 0.019;

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
                {"model",               "WET"},
                {"q",                   "d"},
                {"nonlocal-formfactor", "GvDV2020"},
                {"psi",                 "J/psi"}
            };

            BToKstarCharmonium c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.branching_ratio(),  89.546064, eps);

            TEST_CHECK_RELATIVE_ERROR(c.S_1c_LHCb(),  0.0168472445, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_1s_LHCb(),  0.7373645665, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_3_LHCb(),   0.2111739095, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_4_LHCb(),   0.0485842716, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_8_LHCb(),  -0.0040291232, eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_9_LHCb(),   0.0117985619, eps);

        }
} b_to_kstar_charmonium_GvDV2020_test;


class BToKstarCharmoniumGRvDV2021Test :
    public TestCase
{
    public:
    BToKstarCharmoniumGRvDV2021Test() :
            TestCase("b_to_kstar_charmonium_GRvDV2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]                                = 5.279;
            p["mass::K_d^*"]                              = 0.896;
            p["mass::J/psi"]                              = 3.0969;
            p["mass::psi(2S)"]                            = 3.6860;
            p["mass::D^0"]                                = 1.86723;
            p["b->sccbar::t_0"]                           = 4.0;
            p["b->sccbar::t_s"]                           = -17.4724;
            p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;

            p["B->K^*ccbar::Re{alpha_0^perp}@GRvDV2021"]  = 0.02;
            p["B->K^*ccbar::Im{alpha_0^perp}@GRvDV2021"]  = 0.03;
            p["B->K^*ccbar::Re{alpha_1^perp}@GRvDV2021"]  = 0.04;
            p["B->K^*ccbar::Im{alpha_1^perp}@GRvDV2021"]  = 0.05;
            p["B->K^*ccbar::Re{alpha_2^perp}@GRvDV2021"]  = 0.06;
            p["B->K^*ccbar::Im{alpha_2^perp}@GRvDV2021"]  = 0.07;
            p["B->K^*ccbar::Re{alpha_0^para}@GRvDV2021"]  = 0.08;
            p["B->K^*ccbar::Im{alpha_0^para}@GRvDV2021"]  = 0.09;
            p["B->K^*ccbar::Re{alpha_1^para}@GRvDV2021"]  = 0.010;
            p["B->K^*ccbar::Im{alpha_1^para}@GRvDV2021"]  = 0.011;
            p["B->K^*ccbar::Re{alpha_2^para}@GRvDV2021"]  = 0.012;
            p["B->K^*ccbar::Im{alpha_2^para}@GRvDV2021"]  = 0.013;
            p["B->K^*ccbar::Re{alpha_0^long}@GRvDV2021"]  = 0.014;
            p["B->K^*ccbar::Im{alpha_0^long}@GRvDV2021"]  = 0.015;
            p["B->K^*ccbar::Re{alpha_1^long}@GRvDV2021"]  = 0.016;
            p["B->K^*ccbar::Im{alpha_1^long}@GRvDV2021"]  = 0.017;
            p["B->K^*ccbar::Re{alpha_2^long}@GRvDV2021"]  = 0.018;
            p["B->K^*ccbar::Im{alpha_2^long}@GRvDV2021"]  = 0.019;

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
                {"model",               "WET"},
                {"q",                   "d"},
                {"nonlocal-formfactor", "GRvDV2021"},
                {"psi",                 "J/psi"}
            };

            BToKstarCharmonium c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.branching_ratio(),  49.088606, eps);

            TEST_CHECK_RELATIVE_ERROR(c.S_1c_LHCb(),  0.005445996085,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_1s_LHCb(),  0.745915502935,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_3_LHCb(),  -0.448341852948,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_4_LHCb(),   0.050729332097,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_8_LHCb(),  -0.002147694092,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.S_9_LHCb(),   0.034990413101,  eps);

        }
} b_to_kstar_charmonium_GRvDV2021_test;