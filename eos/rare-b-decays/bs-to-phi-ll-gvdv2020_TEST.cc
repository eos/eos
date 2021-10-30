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
#include <eos/rare-b-decays/bs-to-phi-ll.hh>
#include <eos/utils/complex.hh>

#include <iostream>

using namespace test;
using namespace eos;

class BsToPhiDileptonGvDV2020Test :
    public TestCase
{
    public:
    BsToPhiDileptonGvDV2020Test() :
        TestCase("bs_to_phi_dilepton_GvDV2020_test")
    {
    }

    virtual void run() const
    {
        {
            Parameters p = Parameters::Defaults();

            p["mass::B_s"]                               = 5.36688;
            p["mass::phi"]                               = 1.019461;
            p["mass::J/psi"]                             = 3.0969;
            p["mass::psi(2S)"]                           = 3.6860;
            p["mass::D^0"]                               = 1.86723;
            p["b->sccbar::t_0"]                          = 9.0;
            p["b->sccbar::t_s"]                          = -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;

            p["B_s->phiccbar::Re{alpha_0^perp}@GvDV2020"]  = 0.0002;
            p["B_s->phiccbar::Im{alpha_0^perp}@GvDV2020"]  = 0.0003;
            p["B_s->phiccbar::Re{alpha_1^perp}@GvDV2020"]  = 0.0004;
            p["B_s->phiccbar::Im{alpha_1^perp}@GvDV2020"]  = 0.0005;
            p["B_s->phiccbar::Re{alpha_2^perp}@GvDV2020"]  = 0.0006;
            p["B_s->phiccbar::Im{alpha_2^perp}@GvDV2020"]  = 0.0007;
            p["B_s->phiccbar::Re{alpha_0^para}@GvDV2020"]  = 0.0008;
            p["B_s->phiccbar::Im{alpha_0^para}@GvDV2020"]  = 0.0009;
            p["B_s->phiccbar::Re{alpha_1^para}@GvDV2020"]  = 0.0010;
            p["B_s->phiccbar::Im{alpha_1^para}@GvDV2020"]  = 0.0011;
            p["B_s->phiccbar::Re{alpha_2^para}@GvDV2020"]  = 0.0012;
            p["B_s->phiccbar::Im{alpha_2^para}@GvDV2020"]  = 0.0013;
            p["B_s->phiccbar::Re{alpha_0^long}@GvDV2020"]  = 0.0014;
            p["B_s->phiccbar::Im{alpha_0^long}@GvDV2020"]  = 0.0015;
            p["B_s->phiccbar::Re{alpha_1^long}@GvDV2020"]  = 0.0016;
            p["B_s->phiccbar::Im{alpha_1^long}@GvDV2020"]  = 0.0017;
            p["B_s->phiccbar::Re{alpha_2^long}@GvDV2020"]  = 0.0018;
            p["B_s->phiccbar::Im{alpha_2^long}@GvDV2020"]  = 0.0019;

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

            Options oo;
            oo.declare("model",                   "WET");
            oo.declare("tag",                     "GvDV2020");
            oo.declare("nonlocal-formfactors",    "GvDV2020");
            oo.declare("form-factors",            "BZ2004");
            oo.declare("l",                       "mu");

            static const double eps = 1e-5;
            static const double q2 = 6.0;

            BsToPhiDilepton d(p, oo);
            auto amps = d.amplitudes(q2);


            std::cout << real(amps.a_long_left)   << std::endl;
            std::cout << imag(amps.a_long_left)   << std::endl;
            std::cout << real(amps.a_long_right)  << std::endl;
            std::cout << imag(amps.a_long_right)  << std::endl;
            std::cout << real(amps.a_para_left)   << std::endl;
            std::cout << imag(amps.a_para_left)   << std::endl;
            std::cout << real(amps.a_para_right)  << std::endl;
            std::cout << imag(amps.a_para_right)  << std::endl;
            std::cout << real(amps.a_perp_left)   << std::endl;
            std::cout << imag(amps.a_perp_left)   << std::endl;
            std::cout << real(amps.a_perp_right)  << std::endl;
            std::cout << imag(amps.a_perp_right)  << std::endl;



            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_left),  -1.04871e-09, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_left),  -2.41479e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_right),  5.26771e-11, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_right), -2.95225e-11, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_left),   1.79215e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_left),   5.94937e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_right),  8.0494e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_right),  7.15354e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_left),   3.17172e-11, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_left),  -2.52525e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_right),  3.30606e-10, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_right),  3.87104e-10, eps);

            TEST_CHECK_RELATIVE_ERROR(real(amps.a_time),       -1.10463e-09, eps);
            TEST_CHECK_RELATIVE_ERROR(imag(amps.a_time),       -2.12580e-10, eps);
       }
    }
} bs_to_phi_dilepton_GvDV2020_test;
