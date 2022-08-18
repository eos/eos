/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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
#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class LambdaBToLambda1520DileptonNaiveTest :
    public TestCase
{
    public:
    LambdaBToLambda1520DileptonNaiveTest() :
        TestCase("lambda_b_to_lambda1520_ll_naive_test")
    {
    }

    virtual void run() const
    {
        {
            Parameters p = Parameters::Defaults();
            p["Lambda_b->Lambda(1520)::a^(t12,V)_1@ABR2022"]     =  0.1;
            p["Lambda_b->Lambda(1520)::a^(012,V)_1@ABR2022"]     =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp12,V)_1@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,V)_0@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,V)_1@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(t12,A)_1@ABR2022"]     =  0.1;
            p["Lambda_b->Lambda(1520)::a^(012,A)_1@ABR2022"]     =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp12,A)_1@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,A)_0@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,A)_1@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(012,T)_1@ABR2022"]     =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp12,T)_1@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,T)_0@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,T)_1@ABR2022"]  =  0.1;
            p["Lambda_b->Lambda(1520)::a^(012,T5)_1@ABR2022"]    =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp12,T5)_1@ABR2022"] =  0.1;
            p["Lambda_b->Lambda(1520)::a^(perp32,T5)_1@ABR2022"] =  0.1;

            p["b->s::c3"] = 0.;
            p["b->s::c4"] = 0.;
            p["b->s::c5"] = 0.;
            p["b->s::c6"] = 0.;
            p["sb::mu"] = 4.2;
            p["b->s::Re{c7}"] = -0.29;
            p["b->s::Im{c7}"] = 0.;
            p["b->s::Re{c7'}"] = 0.;
            p["b->s::Im{c7'}"] = 0.;
            p["b->s::c8"] = 0.;
            p["sbmumu::mu"] = 4.2;
            p["b->smumu::Re{c9}"] = 4.2;
            p["b->smumu::Im{c9}"] = 0.;
            p["b->smumu::Re{c9'}"] = 0.;
            p["b->smumu::Im{c9'}"] = 0.;
            p["b->smumu::Re{c10}"] = -4.3;
            p["b->smumu::Im{c10}"] = 0.;
            p["b->smumu::Re{c10'}"] = 0.;
            p["b->smumu::Im{c10'}"] = 0.;

            p["CKM::abs(V_ub)"] = 0.;
            p["CKM::arg(V_ub)"] = 0.;
            p["CKM::abs(V_cb)"] = 0.;
            p["CKM::arg(V_cb)"] = 0.;
            p["CKM::abs(V_tb)"] = 1.;
            p["CKM::arg(V_tb)"] = 0.;
            p["CKM::abs(V_us)"] = 0.;
            p["CKM::arg(V_us)"] = 0.;
            p["CKM::abs(V_cs)"] = 0.;
            p["CKM::arg(V_cs)"] = 0.;
            p["CKM::abs(V_ts)"] = 1.;
            p["CKM::arg(V_ts)"] = 0.;

            p["QED::alpha_e(m_b)"] = 1.;
            p["WET::G_Fermi"] = 1.;
            p["mass::Lambda_b"] = 5.62;
            p["mass::Lambda(1520)"] = 1.52;

            Options oo
            {
                { "model", "WET" },
                { "tag", "Naive" },
                { "form-factors", "ABR2022" },
                { "l", "mu" },
            };

            static const double eps = 1e-5;
            static const double q2 = 1.0;

            LambdaBToLambda1520Dilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR_C(amps.b_perp1_right,  complex<double>(-0.0416087,0.00161803)   , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.b_perp1_left,   complex<double>(-0.0256064,0.00161803)   , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.b_para1_right,  complex<double>(-0.042223,-0.000368419)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.b_para1_left,   complex<double>(-0.0458667,-0.000368419) , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp1_right,  complex<double>(0.0105526,-0.000410356)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp1_left,   complex<double>(0.00649415,-0.000410356) , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para1_right,  complex<double>(0.010993,-0.000347111)   , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para1_left,   complex<double>(0.00756008,-0.000347111) , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp0_right,  complex<double>(0.000820407,-0.0127529)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp0_left,   complex<double>(-0.125306,-0.0127529)    , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para0_right,  complex<double>(-0.00451516,0.00240086)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para0_left,   complex<double>(0.0192293,0.00240086)    , eps);
        }
    }
} lambda_b_to_lambda1520_ll_naive_test;
