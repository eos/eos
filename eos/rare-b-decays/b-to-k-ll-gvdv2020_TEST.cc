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
#include <eos/rare-b-decays/b-to-k-ll.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class BToKDileptonGvDV2020Test :
    public TestCase
{
    public:
    BToKDileptonGvDV2020Test() :
            TestCase("b_to_k_dilepton_GvDV2020_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]                               = 5.279;
            p["mass::K_d"]                               = 0.492;
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 4.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;

            p["B->Kccbar::Re{alpha_0^plus}@GvDV2020"]  =  0.02;
            p["B->Kccbar::Im{alpha_0^plus}@GvDV2020"]  =  0.03;
            p["B->Kccbar::Re{alpha_1^plus}@GvDV2020"]  = -0.04;
            p["B->Kccbar::Im{alpha_1^plus}@GvDV2020"]  = -0.05;
            p["B->Kccbar::Re{alpha_2^plus}@GvDV2020"]  =  0.06;
            p["B->Kccbar::Im{alpha_2^plus}@GvDV2020"]  =  0.07;

            p["b->s::c3"] = -0.005233499106;
            p["b->s::c4"] = -0.08829686414;
            p["b->s::c5"] = 0.0003601965805;
            p["b->s::c6"] = 0.001020749573;
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
            p["b->smumu::Re{cS}"] = 0.5;
            p["b->smumu::Im{cS}"] = 1;
            p["b->smumu::Re{cS'}"] = 0.6;
            p["b->smumu::Im{cS'}"] = 1.1;
            p["b->smumu::Re{cP}"] = 0.7;
            p["b->smumu::Im{cP}"] = 1.2;
            p["b->smumu::Re{cP'}"] = 0.8;
            p["b->smumu::Im{cP'}"] = 1.3;
            p["b->smumu::Re{cT}"] = 0.9;
            p["b->smumu::Im{cT}"] = 1.4;
            p["b->smumu::Re{cT5}"] = 1.0;
            p["b->smumu::Im{cT5}"] = 1.5;
            p["K::a_1@1GeV"] = 0.1;
            p["K::a_2@1GeV"] = 0.1;
            p["B::1/lambda_B_p"] = 1.0 / 0.485;

            Options oo
            {
                {"model"_ok,                   "WET"},
                {"tag"_ok,                     "GvDV2020"},
                {"nonlocal-formfactors"_ok,    "GvDV2020"},
                {"form-factors"_ok,            "BSZ2015"},
                {"l"_ok,                       "mu"},
                {"q"_ok,                       "d"}
            };

            static const double eps = 1e-5;
            static const double q2 = 6.0;

            BToKDilepton c(p, oo);
            auto amps = c.amplitudes(q2);


            TEST_CHECK_RELATIVE_ERROR_C(amps.F_A,  complex<double>(2.803705304, 6.000000000),  eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_V,  complex<double>(116.5855166, 136.0359514),  eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_S,  complex<double>(3.128079910, 5.971788919),  eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_P,  complex<double>(3.752453111, 6.011203332),  eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_T,  complex<double>(6.062177880, 9.430054481),  eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_T5, complex<double>(6.735753201, 10.10362980),  eps);

        }
} b_to_k_dilepton_GvDV2020_test;

class BToKDileptonJavierTest :
    public TestCase
{
    public:
    BToKDileptonJavierTest() :
            TestCase("b_to_k_dilepton_Javier_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["mass::B_d"]      = 5;
            p["mass::K_d"]      = 0.5;
            p["mass::mu"]       = 1e-15;

            p["B->Kccbar::Re{alpha_0^plus}@GvDV2020"]  =  0.01;
            p["B->Kccbar::Im{alpha_0^plus}@GvDV2020"]  =  0.;
            p["B->Kccbar::Re{alpha_1^plus}@GvDV2020"]  =  0.;
            p["B->Kccbar::Im{alpha_1^plus}@GvDV2020"]  =  0.;
            p["B->Kccbar::Re{alpha_2^plus}@GvDV2020"]  =  0.;
            p["B->Kccbar::Im{alpha_2^plus}@GvDV2020"]  =  0.;

            p["B->K::alpha^f+_0@BSZ2015"] = 1.0;
            p["B->K::alpha^f+_1@BSZ2015"] = 0.0;
            p["B->K::alpha^f+_2@BSZ2015"] = 0.0;
            p["B->K::alpha^fT_0@BSZ2015"] = 1.0;
            p["B->K::alpha^fT_1@BSZ2015"] = 0.0;
            p["B->K::alpha^fT_2@BSZ2015"] = 0.0;

            p["b->s::Re{c7}"] = 1.;
            p["b->s::Im{c7}"] = 0.;
            p["b->s::c8"] = 0;
            p["b->smumu::Re{c9}"] = 4;
            p["b->smumu::Im{c9}"] = 0.;
            p["b->smumu::Re{c10}"] = -4;
            p["b->smumu::Im{c10}"] = 0.;
            p["b->s::c3"] = 0.;
            p["b->s::c4"] = 0.;
            p["b->s::c5"] = 0.;
            p["b->s::c6"] = 0.;

            Options oo
            {
                {"model"_ok,                   "WET"},
                {"tag"_ok,                     "GvDV2020"},
                {"nonlocal-formfactors"_ok,    "GvDV2020"},
                {"form-factors"_ok,            "BSZ2015"},
                {"l"_ok,                       "mu"},
                {"q"_ok,                       "d"}
            };

            static const double eps = 1e-5;
            static const double q2 = 1.0;

            auto nff = NonlocalFormFactor<PToP>::make("B->K::GvDV2020", p, oo);
            TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(q2)), -0.0001717492,  eps);
            TEST_CHECK_NEARLY_EQUAL(imag(nff->H_plus(q2)),  0.,        eps);

            BToKDilepton c(p, oo);
            TEST_CHECK_RELATIVE_ERROR(c.two_differential_decay_width(q2, 0),  1.498599e-19,  eps);
            TEST_CHECK_EQUAL(         c.two_differential_decay_width(q2, 1),  0.                );
        }
} b_to_k_dilepton_Javier_test;
