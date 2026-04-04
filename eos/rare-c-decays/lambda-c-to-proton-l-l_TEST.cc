/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2026 Dominik Suelmann
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

#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/rare-c-decays/lambda-c-to-proton-l-l.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class LambdaCToProtonLeptonLeptonTest : public TestCase
{
    public:
        LambdaCToProtonLeptonLeptonTest() :
            TestCase("lambdac_to_proton_l_l_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p = Parameters::Defaults();

                p["Lambda_c->neutron::a^(t,V)_1@BMRvD2022"]     = +0.689012;
                p["Lambda_c->neutron::a^(t,V)_2@BMRvD2022"]     = +0.399494;
                p["Lambda_c->neutron::a^(t,A)_1@BMRvD2022"]     = +0.246169;
                p["Lambda_c->neutron::a^(t,A)_2@BMRvD2022"]     = +0.284617;
                p["Lambda_c->neutron::a^(0,V)_0@BMRvD2022"]     = +0.0423961;
                p["Lambda_c->neutron::a^(0,V)_1@BMRvD2022"]     = +0.736818;
                p["Lambda_c->neutron::a^(0,V)_2@BMRvD2022"]     = +0.413784;
                p["Lambda_c->neutron::a^(0,A)_0@BMRvD2022"]     = -0.021022;
                p["Lambda_c->neutron::a^(0,A)_1@BMRvD2022"]     = +0.189028;
                p["Lambda_c->neutron::a^(0,A)_2@BMRvD2022"]     = +0.13737;
                p["Lambda_c->neutron::a^(0,T)_0@BMRvD2022"]     = -0.294244;
                p["Lambda_c->neutron::a^(0,T)_1@BMRvD2022"]     = +0.0525338;
                p["Lambda_c->neutron::a^(0,T)_2@BMRvD2022"]     = +0.228615;
                p["Lambda_c->neutron::a^(0,T5)_1@BMRvD2022"]    = +0.181135;
                p["Lambda_c->neutron::a^(0,T5)_2@BMRvD2022"]    = +0.204599;
                p["Lambda_c->neutron::a^(perp,V)_0@BMRvD2022"]  = -0.353174;
                p["Lambda_c->neutron::a^(perp,V)_1@BMRvD2022"]  = +0.321992;
                p["Lambda_c->neutron::a^(perp,V)_2@BMRvD2022"]  = +0.431261;
                p["Lambda_c->neutron::a^(perp,A)_1@BMRvD2022"]  = +0.0665121;
                p["Lambda_c->neutron::a^(perp,A)_2@BMRvD2022"]  = +0.201844;
                p["Lambda_c->neutron::a^(perp,T)_1@BMRvD2022"]  = +0.569999;
                p["Lambda_c->neutron::a^(perp,T)_2@BMRvD2022"]  = +0.393516;
                p["Lambda_c->neutron::a^(perp,T5)_0@BMRvD2022"] = -0.0229303;
                p["Lambda_c->neutron::a^(perp,T5)_1@BMRvD2022"] = +0.327394;
                p["Lambda_c->neutron::a^(perp,T5)_2@BMRvD2022"] = +0.223455;

                p["Lambda_c->proton::res_a_rho@GHM2021"]             = +0.54;
                p["Lambda_c->proton::res_a_omega@GHM2021"]           = +0.074;
                p["Lambda_c->proton::res_a_phi@GHM2021"]             = +0.106;
                p["Lambda_c->proton::res_delta_rho@GHM2021"]         = +0.00;
                p["Lambda_c->proton::res_delta_omega_m_rho@GHM2021"] = +M_PI;
                p["Lambda_c->proton::res_delta_phi_m_rho@GHM2021"]   = +M_PI;

                p["QED::alpha_e(m_c)"] = 7.606410e-03;
                p["mass::rho^0"]       = 7.752600e-01;
                p["mass::omega"]       = 7.826600e-01;
                p["mass::phi"]         = 1.019461e+00;
                p["life_time::rho^0"]  = 4.4654813900949794e-24;
                p["life_time::omega"]  = 7.583087061059907e-23;
                p["life_time::phi"]    = 1.5490985100023533e-22;

                p["mass::Lambda_c"]      = 2.286460;
                p["mass::proton"]        = 0.9382721;
                p["life_time::Lambda_c"] = 2.015e-13;

                p["uc::Re{c7}"]     = 0.0;
                p["uc::Im{c7}"]     = 0.0;
                p["ucmumu::Re{c9}"] = 0.0;
                p["ucmumu::Im{c9}"] = 0.0;


                Options oo{
                    {        "model"_ok,       "WET" },
                    { "form-factors"_ok, "BMRvD2022" },
                    {            "l"_ok,        "mu" }
                };

                LambdaCToProtonLeptonLepton d(p, oo);

                const double eps = 1e-6;

                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(0.75), 2.807053787938076e-07, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(0.9), 1.8000240396675447e-07, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(1.0), 5.383860135156469e-07, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(1.1), 2.893673145090639e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(1.5), 2.489875903462182e-09, eps);

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(0.959, 1.122), 3.0625302416679804e-07, eps);

                const complex<double> ckm_factor_ds(-5.383319178425827e-05, 1.427709577227314e-04);

                p["uc::Re{c7}"]  = real(complex<double>(2.2, -1.1) / ckm_factor_ds);
                p["uc::Im{c7}"]  = imag(complex<double>(2.2, -1.1) / ckm_factor_ds);
                p["uc::Re{c7'}"] = real(complex<double>(1.2, 0.5) / ckm_factor_ds);
                p["uc::Im{c7'}"] = imag(complex<double>(1.2, 0.5) / ckm_factor_ds);

                p["ucmumu::Re{c9}"]  = real(complex<double>(0.3, -1.0) / ckm_factor_ds);
                p["ucmumu::Im{c9}"]  = imag(complex<double>(0.3, -1.0) / ckm_factor_ds);
                p["ucmumu::Re{c9'}"] = real(complex<double>(-0.6, 0.2) / ckm_factor_ds);
                p["ucmumu::Im{c9'}"] = imag(complex<double>(-0.6, 0.2) / ckm_factor_ds);

                p["ucmumu::Re{c10}"]  = real(complex<double>(2.0, 1.4) / ckm_factor_ds);
                p["ucmumu::Im{c10}"]  = imag(complex<double>(2.0, 1.4) / ckm_factor_ds);
                p["ucmumu::Re{c10'}"] = real(complex<double>(-1.2, 0.4) / ckm_factor_ds);
                p["ucmumu::Im{c10'}"] = imag(complex<double>(-1.2, 0.4) / ckm_factor_ds);

                Options oobar{
                    {        "model"_ok,       "WET" },
                    { "form-factors"_ok, "BMRvD2022" },
                    {            "l"_ok,        "mu" },
                    { "cp-conjugate"_ok,      "true" }
                };

                LambdaCToProtonLeptonLepton dNP(p, oo);
                LambdaCToProtonLeptonLepton dNPbar(p, oobar);

                // check NP contributions
                TEST_CHECK_RELATIVE_ERROR(dNP.differential_a_fb_leptonic(1.5), 0.14922822534563807, eps);
                TEST_CHECK_RELATIVE_ERROR(dNPbar.differential_a_fb_leptonic(1.5), 0.1565144188529763, eps);
                TEST_CHECK_RELATIVE_ERROR(dNP.differential_f_l(1.5), 0.3114006688420309, eps);
                TEST_CHECK_RELATIVE_ERROR(dNP.differential_branching_ratio(1.5), 1.2337106994356664e-06, eps);
                TEST_CHECK_RELATIVE_ERROR(dNP.integrated_decay_width(0.4, 0.9), 7.547379484042708e-18, eps);

                // Sigma
                TEST_CHECK_RELATIVE_ERROR(0.5
                                                  * (dNP.integrated_a_fb_leptonic_num(0.4, 0.9) / dNP.integrated_decay_width(0.4, 0.9)
                                                     + dNPbar.integrated_a_fb_leptonic_num(0.4, 0.9) / dNPbar.integrated_decay_width(0.4, 0.9)),
                                          0.11172676105989615,
                                          eps);
                // Delta
                TEST_CHECK_RELATIVE_ERROR(0.5
                                                  * (dNP.integrated_a_fb_leptonic_num(0.4, 0.9) / dNP.integrated_decay_width(0.4, 0.9)
                                                     - dNPbar.integrated_a_fb_leptonic_num(0.4, 0.9) / dNPbar.integrated_decay_width(0.4, 0.9)),
                                          -0.028616987399586513,
                                          eps);
                // FL
                TEST_CHECK_RELATIVE_ERROR(dNP.integrated_f_l_num(0.4, 0.9) / dNP.integrated_decay_width(0.4, 0.9), 0.2657365230011504, eps);
            }
        }
} lambdac_to_proton_l_l_test;
