/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2018 Danny van Dyk
 * Copyright (c) 2018 Elena Graverini
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
#include <eos/b-decays/lambdab-to-lambdac2595-l-nu.hh>
#include <eos/maths/complex.hh>

#include <iostream>

using namespace test;
using namespace eos;

class LambdaBToLambdaC2595LeptonNeutrinoTest :
    public TestCase
{
    public:
        LambdaBToLambdaC2595LeptonNeutrinoTest() :
            TestCase("lambda_b_to_lambda_c_2595_l_nu_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["Lambda_b->Lambda_c^*::zeta(q^2_max)@HQET"] =  1.00;
            p["Lambda_b->Lambda_c^*::delta_3b@HQET"]      = -0.14;
            p["Lambda_b->Lambda_c^*::rho@HQET"]           =  0.25;
            p["Lambda_b->Lambda_c^*::rho_3b@HQET"]        =  0.25;

            Options o;

            LambdaBToLambdaC2595LeptonNeutrino d(p, o);

            static const double eps = 7.0e-3;
            static constexpr double s_max = 9.1643031076;
            static constexpr double s_min = 0.011163612964000001;
            TEST_CHECK_RELATIVE_ERROR(d.a_l(s_max - 0.1),                                     0.3170446650730, eps);
            TEST_CHECK_RELATIVE_ERROR(d.b_l(s_max - 0.1),                                     0.2657218590415, eps);
            TEST_CHECK_RELATIVE_ERROR(d.c_l(s_max - 0.1),                                    -0.0122146493426, eps);
            TEST_CHECK_RELATIVE_ERROR(d.a_l(s_max - 3.0),                                     7.6237858620049, eps);
            TEST_CHECK_RELATIVE_ERROR(d.b_l(s_max - 3.0),                                    -1.7643518917646, eps);
            TEST_CHECK_RELATIVE_ERROR(d.c_l(s_max - 3.0),                                    -3.02057556788,   eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(s_min, s_max),             0.0436467078537, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(s_min, s_max), -0.0824034043085, eps);

            Kinematics k
            {
                { "q2_mu_min",  0.0111 }, { "q2_mu_max",  9.164 },
                { "q2_tau_min", 3.1570 }, { "q2_tau_max", 9.164 }
            };
            auto obs_R = Observable::make("Lambda_b->Lambda_c(2595)lnu::R_Lambda_c(2595)", p, k, o);
            TEST_CHECK_RELATIVE_ERROR(obs_R->evaluate(), 0.08896965, eps);
        }
} lambda_b_to_lambda_c_2595_l_nu_test;
