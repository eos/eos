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
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/utils/complex.hh>

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

            p["mass::B_d"]                               = 5.27942;
            p["mass::K_d"]                               = 0.49761;
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 9.0;
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

            p["b->s::Re{c7}"] = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"] = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"] = -0.1827530948;
            p["b->smumu::Re{c9}"] = 4.294489364 + 1;
            p["b->smumu::Im{c9}"] = 0.5;
            p["b->smumu::Re{c9'}"] = 2;
            p["b->smumu::Im{c9'}"] = 1.5;
            p["b->smumu::Re{c10}"] = -4.196294696 + 3;
            p["b->smumu::Im{c10}"] = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;

            Options oo;
            oo.set("model",                   "WilsonScan");
            oo.set("tag",                     "GvDV2020");
            oo.set("nonlocal-formfactors",    "GvDV2020");
            oo.set("form-factors",            "BSZ2015");
            oo.set("l",                       "mu");
            oo.set("q",                       "d");

            static const double eps = 1e-5;
            static const double q2 = 6.0;

            BToKstarDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_left),  -8.09775e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_left),  -2.42046e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_right), -4.46981e-11, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_right), -9.48108e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_left),   1.90312e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_left),   5.54593e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_right),  7.2903e-10,  eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_right),  6.58267e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_left),   8.60191e-13, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_left),  -2.31604e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_right),  2.49618e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_right),  3.00743e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_time),       -7.68624e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_time),       -1.47917e-10, eps);
       }
    }
} b_to_kstar_dilepton_GvDV2020_test;
