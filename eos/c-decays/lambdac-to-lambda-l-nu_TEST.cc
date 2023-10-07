/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
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
#include <eos/c-decays/lambdac-to-lambda-l-nu.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class LambdaCToLambdaLeptonNeutrinoTest :
    public TestCase
{
    public:
        LambdaCToLambdaLeptonNeutrinoTest() :
            TestCase("lambdac_to_lambda_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();

                p["Lambda_c->Lambda::a^(t,V)_1@BMRvD2022"]     = +0.01;
                p["Lambda_c->Lambda::a^(0,V)_0@BMRvD2022"]     = -0.12;
                p["Lambda_c->Lambda::a^(0,V)_1@BMRvD2022"]     = +0.03;
                p["Lambda_c->Lambda::a^(perp,V)_0@BMRvD2022"]  = -0.04;
                p["Lambda_c->Lambda::a^(perp,V)_1@BMRvD2022"]  = +0.05;
                p["Lambda_c->Lambda::a^(t,A)_1@BMRvD2022"]     = -0.06;
                p["Lambda_c->Lambda::a^(0,A)_0@BMRvD2022"]     = +0.07;
                p["Lambda_c->Lambda::a^(0,A)_1@BMRvD2022"]     = -0.08;
                p["Lambda_c->Lambda::a^(perp,A)_1@BMRvD2022"]  = +0.09;
                p["Lambda_c->Lambda::a^(0,T)_0@BMRvD2022"]     = -0.10;
                p["Lambda_c->Lambda::a^(0,T)_1@BMRvD2022"]     = +0.11;
                p["Lambda_c->Lambda::a^(perp,T)_1@BMRvD2022"]  = -0.02;
                p["Lambda_c->Lambda::a^(0,T5)_1@BMRvD2022"]    = +0.13;
                p["Lambda_c->Lambda::a^(perp,T5)_0@BMRvD2022"] = -0.14;
                p["Lambda_c->Lambda::a^(perp,T5)_1@BMRvD2022"] = +0.15;

                Options oo
                {
                    { "model",        "WET"        },
                    { "form-factors", "BMRvD2022"  },
                    { "l",            "mu"         }
                };

                LambdaCToLambdaLeptonNeutrino d(p, oo);

                const double eps = 1e-4;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_decay_width(0.0112, 1.37), 1.2409565e-11, eps);
            }
        }
} lambdac_to_lambda_l_nu_test;
