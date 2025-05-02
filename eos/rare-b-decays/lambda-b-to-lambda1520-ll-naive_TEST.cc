/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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
            p["b->smumu::Im{c9}"] = 0.1;
            p["b->smumu::Re{c9'}"] = 0.;
            p["b->smumu::Im{c9'}"] = 0.;
            p["b->smumu::Re{c10}"] = -4.3;
            p["b->smumu::Im{c10}"] = 0.2;
            p["b->smumu::Re{c10'}"] = 0.;
            p["b->smumu::Im{c10'}"] = 0.;
            p["b->smumu::Re{cS}"] = 0.3;
            p["b->smumu::Im{cS}"] = 0.4;
            p["b->smumu::Re{cS'}"] = 0.;
            p["b->smumu::Im{cS'}"] = 0.;
            p["b->smumu::Re{cP}"] = -0.5;
            p["b->smumu::Im{cP}"] = -0.6;
            p["b->smumu::Re{cP'}"] = 0.;
            p["b->smumu::Im{cP'}"] = 0.;

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
                { "model"_ok, "WET" },
                { "tag"_ok, "Naive" },
                { "form-factors"_ok, "ABR2022" },
                { "l"_ok, "mu" },
            };

            static const double eps = 1e-5;
            static const double q2 = 1.0;

            LambdaBToLambda1520Dilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR_C(amps.b_perp1_right,  complex<double>(-0.0411362,0.00215154), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.b_perp1_left,   complex<double>(-0.0253152,0.0014157), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.b_para1_right,  complex<double>(-0.0417436,-0.000489896), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.b_para1_left,   complex<double>(-0.0453458,-0.000322348), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp1_right,  complex<double>(0.0104327,-0.000545661), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp1_left,   complex<double>(0.00642029,-0.000359041), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para1_right,  complex<double>(0.0108680,-0.000461563), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para1_left,   complex<double>(0.00747410,-0.000303705), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp0_right,  complex<double>(0.00081105,-0.0169579), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp0_left,   complex<double>(-0.123883,-0.0111582), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para0_right,  complex<double>(-0.00446388,0.00319248), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para0_left,   complex<double>(0.0190109,0.00210063), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perpt_right,  complex<double>(-0.0629889,0.00292971), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perpt_left,   complex<double>(0.0629889,-0.00292971), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_parat_right,  complex<double>(0.0121739,-0.000566226), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_parat_left,   complex<double>(-0.0121739,0.000566226), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perpS_right,  complex<double>(-0.000716257,-0.000716257), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perpS_left,   complex<double>(0.00286503,0.00358129), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_paraS_right,  complex<double>(0.000689419,0.000689419), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_paraS_left,   complex<double>(-0.00275768,-0.0034471), eps);
        }
    }
} lambda_b_to_lambda1520_ll_naive_test;
