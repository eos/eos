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
#include <eos/rare-b-decays/lambda-b-to-lambda1520-gamma.hh>

using namespace test;
using namespace eos;

class LambdaBToLambda1520GammaNaiveTest :
    public TestCase
{
    public:
        LambdaBToLambda1520GammaNaiveTest() :
            TestCase("lambdab_to_lambda1520_gamma_naive_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
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

            LambdaBToLambda1520Gamma d(p, oo);

            TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(),  847315588516.,   eps);
        }
} lambdab_to_lambda1520_gamma_naive_test;
