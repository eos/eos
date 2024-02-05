/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
 * Copyright (c) 2017 Elena Graverini
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
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/maths/complex.hh>

#include <iostream>

using namespace test;
using namespace eos;

class LambdaBToLambdaC2625LeptonNeutrinoTest :
    public TestCase
{
    public:
        LambdaBToLambdaC2625LeptonNeutrinoTest() :
            TestCase("lambda_b_to_lambda_c_2625_l_nu_test")
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

            LambdaBToLambdaC2625LeptonNeutrino d(p, o);

            static const double eps = 5.0e-3;
            static constexpr double s_max = 8.948473960000001;
            static constexpr double s_min = 0.011163612964000001;
            TEST_CHECK_RELATIVE_ERROR(d.a_l(s_max - 0.1),                                     1.2718441467069,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.b_l(s_max - 0.1),                                     1.4537066913760,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.c_l(s_max - 0.1),                                     0.5008229019503,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.a_l(s_max - 3.0),                                     9.3694820191993,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.b_l(s_max - 3.0),                                     2.0707523203345,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.c_l(s_max - 3.0),                                    -1.9517657097361,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(s_min, s_max),             0.0443817800606,  eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(s_min, s_max),  0.0392696772213,  eps);

            Kinematics k
            {
                { "q2_mu_min",  0.0111 }, { "q2_mu_max",  8.948 },
                { "q2_tau_min", 3.1570 }, { "q2_tau_max", 8.948 }
            };
            auto obs_R = Observable::make("Lambda_b->Lambda_c(2625)lnu::R_Lambda_c(2625)", p, k, o);
            TEST_CHECK_RELATIVE_ERROR(obs_R->evaluate(), 0.0994558945773, eps);
        }
} lambda_b_to_lambda_c_2625_l_nu_test;
