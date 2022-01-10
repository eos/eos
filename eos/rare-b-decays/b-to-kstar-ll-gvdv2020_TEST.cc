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
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

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

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_left),  -1.37402e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_left),  -2.96243e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_right),  7.46919e-12, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_right), -1.74463e-12, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_left),   5.67827e-12, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_left),   7.31160e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_right),  1.07717e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_right),  9.27529e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_left),   1.62003e-11, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_left),  -2.63585e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_right),  6.33040e-11, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_right),  7.44445e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_time),       -1.45544e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_time),       -2.80092e-11, eps);
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

            auto nff = NonlocalFormFactor<nff::PToV>::make("B->K^*::GvDV2020", p, oo);
            TEST_CHECK_NEARLY_EQUAL(real(nff->H_perp(q2)),  0.0063082395,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_perp(q2)),  0.,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(nff->H_para(q2)),  0.0063082395,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_para(q2)),  0.,            eps);
            TEST_CHECK_NEARLY_EQUAL(real(nff->H_long(q2)), -0.0001722913,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_long(q2)),  0.,            eps);

            BToKstarDilepton c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.differential_j_1s(q2),  2.17500e-20,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_1c(q2),  1.35530e-19,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_2s(q2),  7.25000e-21,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_2c(q2), -1.35530e-19,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_3(q2),  -5.82559e-21,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_4(q2),   3.70203e-20,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_5(q2),  -4.24896e-20,  eps);
            TEST_CHECK_RELATIVE_ERROR(c.differential_j_6s(q2), -2.23092e-20,  eps);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_6c(q2),  0.0,  1e-20);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_7(q2),   0.0,  1e-20);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_8(q2),   0.0,  1e-20);
            TEST_CHECK_NEARLY_EQUAL(c.differential_j_9(q2),   0.0,  1e-20);

        }
} b_to_kstar_dilepton_Javier_test;
