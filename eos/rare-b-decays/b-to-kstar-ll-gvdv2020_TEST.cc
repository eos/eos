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
#include <eos/observable.hh>
#include <eos/maths/complex.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>

#include <iostream>

using namespace test;
using namespace eos;

class BToKstarDileptonGvDV2020Test :
    public TestCase
{
    public:
    BToKstarDileptonGvDV2020Test() :
        TestCase("b_to_kstar_dilepton_GvDV2020_test")
    {
    }

    virtual void run() const
    {
        {
            Parameters p = Parameters::Defaults();

            p["mass::B_d"]                               = 5.279;
            p["mass::K_d^*"]                             = 0.896;
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 4.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;

            p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 0.0002;
            p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 0.0003;
            p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 0.0004;
            p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 0.0005;
            p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 0.0006;
            p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 0.0007;
            p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 0.0008;
            p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 0.0009;
            p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 0.0010;
            p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 0.0011;
            p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 0.0012;
            p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 0.0013;
            p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 0.0014;
            p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 0.0015;
            p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 0.0016;
            p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 0.0017;
            p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 0.0018;
            p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 0.0019;

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
            p["sb::mu"] = 4.2;
            p["b->s::Re{c7}"] = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"] = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"] = -0.1827530948;
            p["sbmumu::mu"] = 4.2;
            p["b->smumu::Re{c9}"] = 4.294489364 + 1;
            p["b->smumu::Im{c9}"] = 0.5;
            p["b->smumu::Re{c9'}"] = 2;
            p["b->smumu::Im{c9'}"] = 1.5;
            p["b->smumu::Re{c10}"] = -4.196294696 + 3;
            p["b->smumu::Im{c10}"] = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;

            Options oo
            {
                {"model",                   "WET"},
                {"tag",                     "GvDV2020"},
                {"nonlocal-formfactors",    "GvDV2020"},
                {"form-factors",            "BSZ2015"},
                {"l",                       "mu"},
                {"q",                       "d"}
            };

            static const double eps = 1e-5;
            static const double q2 = 6.0;

            BToKstarDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_left),  -1.37556e-10,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_left),  -2.97866e-11,  eps);
            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_right),  7.31512e-12,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_right), -1.90691e-12,  eps);
            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_left),   5.95688e-12,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_left),   7.3417e-11,   eps);
            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_right),  1.07997e-10,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_right),  9.3054e-11,   eps);
            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_left),   1.60422e-11,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_left),  -2.65284e-11,  eps);
            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_right),  6.31325e-11,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_right),  7.42463e-11,  eps);
            TEST_CHECK_RELATIVE_ERROR(real(amps.a_time),       -1.45112e-10,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_time),       -2.79261e-11,  eps);
       }
    }
} b_to_kstar_dilepton_GvDV2020_test;

class BToKstarDileptonJavierTest :
    public TestCase
{
    public:
    BToKstarDileptonJavierTest() :
            TestCase("b_to_kstar_dilepton_Javier_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]      = 5;
            p["mass::K_d^*"]    = 0.9;
            p["mass::mu"]       = 1e-15;
            p["mass::b(MSbar)"] = 4;
            p["mass::s(2GeV)"]  = 0.;

            p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 0.01;
            p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 0.01;
            p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 0.01;
            p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 0.;
            p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 0.;

            p["B->K^*::alpha^A0_0@BSZ2015"]  = 1.0;
            p["B->K^*::alpha^A0_1@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^A0_2@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^A1_0@BSZ2015"]  = 1.0;
            p["B->K^*::alpha^A1_1@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^A1_2@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^A12_1@BSZ2015"] = 0.0;
            p["B->K^*::alpha^A12_2@BSZ2015"] = 0.0;
            p["B->K^*::alpha^V_0@BSZ2015"]  = 1.0;
            p["B->K^*::alpha^V_1@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^V_2@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^T1_0@BSZ2015"]  = 1.0;
            p["B->K^*::alpha^T1_1@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^T1_2@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^T2_1@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^T2_2@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^T23_0@BSZ2015"]  = 1.0;
            p["B->K^*::alpha^T23_1@BSZ2015"]  = 0.0;
            p["B->K^*::alpha^T23_2@BSZ2015"]  = 0.0;

            p["b->s::Re{c7}"] = 1.;
            p["b->s::Im{c7}"] = 0.;
            p["b->s::c8"] = 0;
            p["b->smumu::Re{c9}"] = 4;
            p["b->smumu::Im{c9}"] = 0.;
            p["b->smumu::Re{c10}"] = -4;
            p["b->smumu::Im{c10}"] = 0.;
            p["b->s::c1"] = 0.;
            p["b->s::c2"] = 0.;
            p["b->s::c3"] = 0.;
            p["b->s::c4"] = 0.;
            p["b->s::c5"] = 0.;
            p["b->s::c6"] = 0.;
            p["b->smumu::Re{cS}"] = 0.;
            p["b->smumu::Im{cS}"] = 0.;
            p["b->smumu::Re{cS'}"] = 0.;
            p["b->smumu::Im{cS'}"] = 0.;
            p["b->smumu::Re{cT}"] = 0.;
            p["b->smumu::Im{cT}"] = 0.;
            p["b->smumu::Re{cT5}"] = 0.;
            p["b->smumu::Im{cT5}"] = 0.;

            Options oo
            {
                {"model",                   "WET"},
                {"tag",                     "GvDV2020"},
                {"nonlocal-formfactors",    "GvDV2020"},
                {"form-factors",            "BSZ2015"},
                {"l",                       "mu"},
                {"q",                       "d"}
            };

            static const double eps = 1e-5;
            static const double q2 = 1.0;

            auto nff = NonlocalFormFactor<PToV>::make("B->K^*::GvDV2020", p, oo);
            TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(q2)),  0.0063082395,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(q2)),  0.,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(q2)),  0.0063082395,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(q2)),  0.,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(q2)), -0.0001722913,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(q2)),  0.,            eps);

            BToKstarDilepton c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.differential_j_1s(q2),  2.15865e-20, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_1c(q2),  1.35571e-19, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_2s(q2),  7.19550e-21, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_2c(q2), -1.35571e-19, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_3(q2),  -5.80454e-21, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_4(q2),   3.69120e-20, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_5(q2),  -4.23428e-20, eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_6s(q2), -2.21891e-20, eps);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_6c(q2),  0.0,  1e-20);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_7(q2),   0.0,  1e-20);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_8(q2),   0.0,  1e-20);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_9(q2),   0.0,  1e-20);

        }
} b_to_kstar_dilepton_Javier_test;

class BToKstarDileptonJavierGRvDV2022Test :
    public TestCase
{
    public:
    BToKstarDileptonJavierGRvDV2022Test() :
            TestCase("b_to_kstar_dilepton_Javier_GRvDV2022_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]      = 5.279;
            p["life_time::B_d"] = 1.519e-12;
            p["mass::K_d^*"]    = 0.896;
            p["mass::D^0"]      = 1.86;
            p["mass::mu"]       = 0.1056583715;
            p["mass::b(MSbar)"] = 4.1884; // Ad hoc value to have mb ~ 4.18 GeV @ mu = 1.5 GeV
            p["mass::s(2GeV)"]  = 0.;
            p["mass::J/psi"]    = 3.0969;
            p["mass::psi(2S)"]  = 3.6861;

            p["B->K^*::alpha^A0_0@BSZ2015"]  = 3.41906028e-01;
            p["B->K^*::alpha^A0_1@BSZ2015"]  = -1.11803025e+00;
            p["B->K^*::alpha^A0_2@BSZ2015"]  = 2.16717478e+00;
            p["B->K^*::alpha^A1_0@BSZ2015"]  = 2.89506064e-01;
            p["B->K^*::alpha^A1_1@BSZ2015"]  = 4.60618884e-01;
            p["B->K^*::alpha^A1_2@BSZ2015"]  = 1.21515277e+00;
            p["B->K^*::alpha^A12_1@BSZ2015"] = 5.49963772e-01;
            p["B->K^*::alpha^A12_2@BSZ2015"] = 5.78051325e-01;
            p["B->K^*::alpha^V_0@BSZ2015"]  = 3.63194834e-01;
            p["B->K^*::alpha^V_1@BSZ2015"]  = -1.08676934e+00;
            p["B->K^*::alpha^V_2@BSZ2015"]  = 2.74860257e+00;
            p["B->K^*::alpha^T1_0@BSZ2015"]  = 3.19511637e-01;
            p["B->K^*::alpha^T1_1@BSZ2015"]  = -9.55695695e-01;
            p["B->K^*::alpha^T1_2@BSZ2015"]  = 2.09828722e+00;
            p["B->K^*::alpha^T2_1@BSZ2015"]  = 5.97421164e-01;
            p["B->K^*::alpha^T2_2@BSZ2015"]  = 1.69820398e+00;
            p["B->K^*::alpha^T23_0@BSZ2015"]  = 6.21548438e-01;
            p["B->K^*::alpha^T23_1@BSZ2015"]  = 9.70443767e-01;
            p["B->K^*::alpha^T23_2@BSZ2015"]  = 1.80719475e+00;

            p["B->K^*ccbar::Re_Hhat_at_m7_perp@GRvDV2022"] = 1.78049776e-04;
            p["B->K^*ccbar::Im_Hhat_at_m7_perp@GRvDV2022"] = 5.89058214e-06;
            p["B->K^*ccbar::Re_Hhat_at_m5_perp@GRvDV2022"] = 1.86107972e-04;
            p["B->K^*ccbar::Im_Hhat_at_m5_perp@GRvDV2022"] = 6.91447589e-06;
            p["B->K^*ccbar::Re_Hhat_at_m3_perp@GRvDV2022"] = 1.86697882e-04;
            p["B->K^*ccbar::Im_Hhat_at_m3_perp@GRvDV2022"] = 8.18047904e-06;
            p["B->K^*ccbar::Re_Hhat_at_m1_perp@GRvDV2022"] = 1.77949850e-04;
            p["B->K^*ccbar::Im_Hhat_at_m1_perp@GRvDV2022"] = 9.76018584e-06;
            p["B->K^*ccbar::Abs_Hhat_at_Jpsi_perp@GRvDV2022"] = 6.63880230e-04;
            p["B->K^*ccbar::Arg_Hhat_at_Jpsi_perp_minus_long@GRvDV2022"] = -2.01580763e-01;
            p["B->K^*ccbar::Re_Hhat_at_m7_para@GRvDV2022"] = 1.80068875e-04;
            p["B->K^*ccbar::Im_Hhat_at_m7_para@GRvDV2022"] = 5.87427274e-06;
            p["B->K^*ccbar::Re_Hhat_at_m5_para@GRvDV2022"] = 1.88287075e-04;
            p["B->K^*ccbar::Im_Hhat_at_m5_para@GRvDV2022"] = 6.88597062e-06;
            p["B->K^*ccbar::Re_Hhat_at_m3_para@GRvDV2022"] = 1.88491234e-04;
            p["B->K^*ccbar::Im_Hhat_at_m3_para@GRvDV2022"] = 8.14322578e-06;
            p["B->K^*ccbar::Re_Hhat_at_m1_para@GRvDV2022"] = 1.78614668e-04;
            p["B->K^*ccbar::Im_Hhat_at_m1_para@GRvDV2022"] = 9.72861253e-06;
            p["B->K^*ccbar::Abs_Hhat_at_Jpsi_para@GRvDV2022"] = 7.05500789e-04;
            p["B->K^*ccbar::Arg_Hhat_at_Jpsi_para_minus_long@GRvDV2022"] = 2.01643825e-01;
            p["B->K^*ccbar::Re_Hhat_at_m7_long@GRvDV2022"] = -9.11586797e-07;
            p["B->K^*ccbar::Im_Hhat_at_m7_long@GRvDV2022"] = -7.90619548e-07;
            p["B->K^*ccbar::Re_Hhat_at_m5_long@GRvDV2022"] = 5.10304440e-06;
            p["B->K^*ccbar::Im_Hhat_at_m5_long@GRvDV2022"] = -7.11624621e-07;
            p["B->K^*ccbar::Re_Hhat_at_m3_long@GRvDV2022"] = 7.06738258e-06;
            p["B->K^*ccbar::Im_Hhat_at_m3_long@GRvDV2022"] = -5.47414002e-07;
            p["B->K^*ccbar::Re_Hhat_at_m1_long@GRvDV2022"] = 3.86577153e-06;
            p["B->K^*ccbar::Im_Hhat_at_m1_long@GRvDV2022"] = -2.38953430e-07;
            p["B->K^*ccbar::Abs_Hhat_at_Jpsi_long@GRvDV2022"] = 6.57206533e-04;

            p["B->K^*ccbar::Abs_Hhat_at_psi2S_perp@GRvDV2022"] = 0.;
            p["B->K^*ccbar::Arg_Hhat_at_psi2S_perp_minus_long@GRvDV2022"] = 0.;
            p["B->K^*ccbar::Abs_Hhat_at_psi2S_para@GRvDV2022"] = 0.;
            p["B->K^*ccbar::Arg_Hhat_at_psi2S_para_minus_long@GRvDV2022"] = 0.;
            p["B->K^*ccbar::Abs_Hhat_at_psi2S_long@GRvDV2022"] = 0.;
            p["B->K^*ccbar::Arg_Hhat_at_Jpsi_long@GRvDV2022"] = 0.;
            p["B->K^*ccbar::Arg_Hhat_at_psi2S_long@GRvDV2022"] = 0.;

            p["b->sccbar::t_0"] = 4.0;
            p["b->sccbar::t_s"] = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"] = 0.000181;

            p["b->s::Re{c7}"] = 0.3;
            p["b->s::Im{c7}"] = 0.;
            p["b->s::c8"] = 0;
            p["b->smumu::Re{c9}"] = 4;
            p["b->smumu::Im{c9}"] = 0.;
            p["b->smumu::Re{c10}"] = -4;
            p["b->smumu::Im{c10}"] = 0.;
            p["b->s::c1"] = 0.;
            p["b->s::c2"] = 0.;
            p["b->s::c3"] = 0.;
            p["b->s::c4"] = 0.;
            p["b->s::c5"] = 0.;
            p["b->s::c6"] = 0.;
            p["b->smumu::Re{cS}"] = 0.;
            p["b->smumu::Im{cS}"] = 0.;
            p["b->smumu::Re{cS'}"] = 0.;
            p["b->smumu::Im{cS'}"] = 0.;
            p["b->smumu::Re{cT}"] = 0.;
            p["b->smumu::Im{cT}"] = 0.;
            p["b->smumu::Re{cT5}"] = 0.;
            p["b->smumu::Im{cT5}"] = 0.;

            Options oo
            {
                {"model",                   "WET"},
                {"tag",                     "GvDV2020"},
                {"nonlocal-formfactor",     "GRvDV2022order5"},
                {"form-factors",            "BSZ2015"},
                {"l",                       "mu"},
                {"q",                       "d"}
            };

            static const double eps = 5e-4;
            static const double q2 = 1.0;

            auto nff = NonlocalFormFactor<PToV>::make("B->K^*::GRvDV2022order5", p, oo);
            TEST_CHECK_RELATIVE_ERROR(real(nff->Hhat_perp(q2)),  0.000159636,   eps);
            TEST_CHECK_RELATIVE_ERROR(imag(nff->Hhat_perp(q2)),  0.0000113099,  eps);
            TEST_CHECK_RELATIVE_ERROR(real(nff->H_perp(q2)),     0.000115032,   eps);
            TEST_CHECK_RELATIVE_ERROR(imag(nff->H_perp(q2)),     8.14984e-06,   eps);
            TEST_CHECK_RELATIVE_ERROR(real(nff->Hhat_para(q2)),  0.0001583608,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(nff->Hhat_para(q2)),  1.205929e-05,  eps);
            TEST_CHECK_RELATIVE_ERROR(real(nff->H_para(q2)),     0.000114113,   eps);
            TEST_CHECK_RELATIVE_ERROR(imag(nff->H_para(q2)),     8.68978e-06,   eps);
            TEST_CHECK_RELATIVE_ERROR(real(nff->Hhat_long(q2)), -3.68996e-06,   eps);
            TEST_CHECK_RELATIVE_ERROR(imag(nff->Hhat_long(q2)),  2.74982e-07,   eps);
            TEST_CHECK_RELATIVE_ERROR(real(nff->H_long(q2)),     6.89722e-08,   eps);
            TEST_CHECK_RELATIVE_ERROR(imag(nff->H_long(q2)),    -5.13967e-09,   eps);


            Kinematics k
            {
                {"q2", q2},
            };

            // Javier's results for c7 = 1.0
            // J's = { 9.88734e-21, 1.60842e-20,
            // 3.11517e-21, -1.4707e-20, -2.6687e-22,
            // 6.23113e-21, -1.12587e-20, -6.47412e-21, 0,
            // 1.00589e-22, -5.21581e-23, -1.72649e-24 }

            BToKstarDilepton c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.differential_j_1s(q2),  3.0/4.0*9.88734e-21,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_1c(q2),  3.0/4.0*1.60842e-20,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_2s(q2),  3.0/4.0*3.11517e-21,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_2c(q2), -3.0/4.0*1.4707e-20,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_3(q2),  -3.0/4.0*2.6687e-22,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_4(q2),   3.0/4.0*6.23113e-21,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_5(q2),  -3.0/4.0*1.12587e-20,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_6s(q2), -3.0/4.0*6.47412e-21,  eps);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_6c(q2),  0.0,  1e-22);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_7(q2),   3.0/4.0*1.00589e-22,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_8(q2),  -3.0/4.0*5.21581e-23,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_9(q2),  -3.0/4.0*1.72649e-24,  eps);


            // Javier's results
            // { BR, FL, AFB,
            //   P1, P2, P3, P4p, P5p, P6p, P8p,
            //   S3, S4, S5, S7, S8, S9} =
            // { 0.669443,    0.542507,  0.167357,
            //  -0.0428339,  -0.259782,  0.000138555, 0.920583,   -0.831674,   -0.00743048,  0.0077058,
            //  -0.00919817,  0.214768, -0.388051,    0.00346699, -0.00179773, -0.0000595069 }

            TEST_CHECK_RELATIVE_ERROR(
                 Observable::make("B->K^*ll::dBR/ds", p, k, oo)->evaluate(),    0.669443e-07, eps );
            TEST_CHECK_RELATIVE_ERROR(
                 Observable::make("B->K^*ll::F_L(q2)", p, k, oo)->evaluate(),   0.542507,     eps );
            TEST_CHECK_RELATIVE_ERROR(
                 Observable::make("B->K^*ll::A_FB(q2)", p, k, oo)->evaluate(), -0.167357,     eps );

            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::P'_4(q2)", p, k, oo)->evaluate(),  0.920583,    eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::P'_5(q2)", p, k, oo)->evaluate(), -0.831674,    eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::P'_6(q2)", p, k, oo)->evaluate(), -0.00743048,  eps );

            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::S_3(q2)", p, k, oo)->evaluate(), -0.00919817,   eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::S_4(q2)", p, k, oo)->evaluate(),  0.214768,     eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::S_5(q2)", p, k, oo)->evaluate(), -0.388051,     eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::S_7(q2)", p, k, oo)->evaluate(),  0.00346699,   eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::S_8(q2)", p, k, oo)->evaluate(), -0.00179773,   eps );
            TEST_CHECK_RELATIVE_ERROR(
                Observable::make("B->K^*ll::S_9(q2)", p, k, oo)->evaluate(), -0.0000595069, eps );

        }
} b_to_kstar_dilepton_Javier_GRvDV2022_test;
