/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Thomas Blake
 * Copyright (c) 2019 Danny van Dyk
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
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>

using namespace test;
using namespace eos;

class LambdaBToLambdaDileptonLowRecoilTest :
    public TestCase
{
    public:
        LambdaBToLambdaDileptonLowRecoilTest() :
            TestCase("lambda_b_to_lambda_dilepton_test")
        {
        }

        virtual void run() const
        {
            // Standard Model
            {
                static const double eps = 1e-4;

                // unpolarised SM
                {
                    Options oo;
                    oo.declare("model",                   "WET");
                    oo.declare("q",                       "d");
                    oo.declare("l",                       "mu");
                    oo.declare("production-polarisation", "unpolarised");

                    Parameters p = Parameters::Defaults();
                    p["mass::Lambda_b"] =  5.6194;
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

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), 8.250965481e-08, eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),  0.3550388404,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),  0.2899223192,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.2437315574,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2054611527,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.158558312,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),  0.1838396079,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.02081000733,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.09222727907,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),  6.268957094e-05,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0003204765254, eps);

                    TEST_CHECK_EQUAL(d.integrated_m11(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m12(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m13(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m14(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m15(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m16(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m17(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m18(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m19(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m20(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m21(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m22(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m23(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m24(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m25(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m26(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m27(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m28(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m29(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m30(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m31(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m32(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m33(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m34(15.0,19.0), 0);
                }

                // LHCb-polarised SM
                {
                    Options oo
                    {
                        {"model",                   "WET"},
                        {"q",                       "d"},
                        {"l",                       "mu"},
                        {"production-polarisation", "LHCb"}
                    };

                    Parameters p = Parameters::Defaults();
                    p["mass::Lambda_b"] =  5.6194;
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

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), 8.250965481e-08, eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),   0.3550388404,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),   0.2899223192,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.2437315574,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2054611527,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.158558312,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),   0.1838396079,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.02081000733,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.09222727907,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),   6.268957094e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0003204765254, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m11(15.0,19.0), -0.004383443052,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m12(15.0,19.0),  0.01481853383,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m13(15.0,19.0), -0.01718127177,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m14(15.0,19.0),  0.002508288396,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m15(15.0,19.0), -0.01116780774,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m16(15.0,19.0),  0.009388539593,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m17(15.0,19.0),  0.005590173438,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m18(15.0,19.0),  0.001256615897,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m19(15.0,19.0),  1.10643968e-05,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m20(15.0,19.0), -6.528475969e-06, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m21(15.0,19.0),  9.830834396e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m22(15.0,19.0), -0.0001927931588, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m23(15.0,19.0), -0.01876460111,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m24(15.0,19.0),  0.02062474436,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m25(15.0,19.0), -7.118517544e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m26(15.0,19.0),  0.0001096620272, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m27(15.0,19.0),  0.01335263112,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m28(15.0,19.0), -0.01194850199,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m29(15.0,19.0),  0,               eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m30(15.0,19.0), -2.156083411e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m31(15.0,19.0),  0,               eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m32(15.0,19.0), -0.002612205688,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m33(15.0,19.0), -0.002820345385,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m34(15.0,19.0),  1.055598677e-05, eps);
                }

                // unpolarised BMP
                {
                    Options oo
                    {
                        {"model",                   "WET"},
                        {"q",                       "d"},
                        {"l",                       "mu"},
                        {"production-polarisation", "unpolarised"}
                    };

                    Parameters p = Parameters::Defaults();
                    p["sbmumu::mu"] = 4.2;
                    p["b->smumu::Re{c9}"]  = +3.2734;
                    p["b->smumu::Re{c9'}"] = +1.0000;
                    p["mass::Lambda_b"]    =  5.6194;
                    p["CKM::abs(V_ub)"]    =  0.003631275231633653;
                    p["CKM::arg(V_ub)"]    = -1.210765774253535;
                    p["CKM::abs(V_cb)"]    =  0.041996951916414726;
                    p["CKM::arg(V_cb)"]    =  0.0;
                    p["CKM::abs(V_tb)"]    =  0.9991111344469873;
                    p["CKM::arg(V_tb)"]    =  0.0;
                    p["CKM::abs(V_us)"]    =  0.22534851424944366;
                    p["CKM::arg(V_us)"]    =  0.0;
                    p["CKM::abs(V_cs)"]    =  0.9734061815416853;
                    p["CKM::arg(V_cs)"]    = -3.304199362533668e-05;
                    p["CKM::abs(V_ts)"]    =  0.04121212396309175;
                    p["CKM::arg(V_ts)"]    = -3.1230250224697222;

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);


                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), 6.367037677e-08, eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),   0.3572380627,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),   0.2855238745,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.234809903,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2080679346,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.1618563161,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),   0.144676431,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.01842966894,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.01208460395,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),   0.000632757932,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0004266646439, eps);

                    TEST_CHECK_EQUAL(d.integrated_m11(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m12(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m13(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m14(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m15(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m16(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m17(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m18(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m19(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m20(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m21(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m22(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m23(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m24(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m25(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m26(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m27(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m28(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m29(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m30(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m31(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m32(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m33(15.0,19.0), 0);
                    TEST_CHECK_EQUAL(d.integrated_m34(15.0,19.0), 0);
                }

                // LHCb-polarised BMP
                {
                    Options oo
                    {
                        {"model",                   "WET"},
                        {"q",                       "d"},
                        {"l",                       "mu"},
                        {"production-polarisation", "LHCb"}
                    };

                    Parameters p = Parameters::Defaults();
                    p["sbmumu::mu"] = 4.2;
                    p["b->smumu::Re{c9}"]  = +3.2734;
                    p["b->smumu::Re{c9'}"] = +1.0000;
                    p["mass::Lambda_b"]    =  5.6194;
                    p["CKM::abs(V_ub)"]    =  0.003631275231633653;
                    p["CKM::arg(V_ub)"]    = -1.210765774253535;
                    p["CKM::abs(V_cb)"]    =  0.041996951916414726;
                    p["CKM::arg(V_cb)"]    =  0.0;
                    p["CKM::abs(V_tb)"]    =  0.9991111344469873;
                    p["CKM::arg(V_tb)"]    =  0.0;
                    p["CKM::abs(V_us)"]    =  0.22534851424944366;
                    p["CKM::arg(V_us)"]    =  0.0;
                    p["CKM::abs(V_cs)"]    =  0.9734061815416853;
                    p["CKM::arg(V_cs)"]    = -3.304199362533668e-05;
                    p["CKM::abs(V_ts)"]    =  0.04121212396309175;
                    p["CKM::arg(V_ts)"]    = -3.1230250224697222;

                    LambdaBToLambdaDilepton<LowRecoil> d(p, oo);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.0), 6.367037677e-08, eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m1(15.0,19.0),   0.3572380627,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m2(15.0,19.0),   0.2855238745,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m3(15.0,19.0),  -0.234809903,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m4(15.0,19.0),  -0.2080679346,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m5(15.0,19.0),  -0.1618563161,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m6(15.0,19.0),   0.144676431,     eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m7(15.0,19.0),  -0.01842966894,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m8(15.0,19.0),  -0.01208460395,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m9(15.0,19.0),   0.000632757932,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m10(15.0,19.0), -0.0004266646439, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m11(15.0,19.0), -0.004318842853,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m12(15.0,19.0),  0.01512675851,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m13(15.0,19.0), -0.01352116178,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m14(15.0,19.0),  0.00276243053,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m15(15.0,19.0), -0.01099837965,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m16(15.0,19.0),  0.009044877463,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m17(15.0,19.0),  0.00303700215,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m18(15.0,19.0),  0.001302784642,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m19(15.0,19.0), -0.0003501096568, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m20(15.0,19.0), -8.691650257e-06, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m21(15.0,19.0),  7.954450124e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m22(15.0,19.0), -0.0002566741022, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m23(15.0,19.0), -0.01899696629,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m24(15.0,19.0),  0.01711256795,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m25(15.0,19.0), -6.720120155e-06, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m26(15.0,19.0),  0.0001459979315, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m27(15.0,19.0),  0.01337143711,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m28(15.0,19.0), -0.01169255641,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m29(15.0,19.0),  0,               eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m30(15.0,19.0),  0.0002233789303, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m31(15.0,19.0),  0,               eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m32(15.0,19.0), -0.00096592011,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m33(15.0,19.0), -0.00181996628,   eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_m34(15.0,19.0),  0.0001388284105, eps);
                }
            }
        }
} lambda_b_to_lambda_dilepton_low_recoil_test;
