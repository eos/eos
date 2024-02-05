/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2019 Danny van Dyk
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
#include <eos/b-decays/b-to-psd-l-nu.hh>
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

class BToPiLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToPiLeptonNeutrinoTest() :
            TestCase("b_to_pi_l_nu_test")
        {
        }

        virtual void run() const
        {
            // Using IKMvD2014 inputs for V_ub and form factors,
            // from the combined fit to B->pilnu data and LCSR.
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_ub)"]        =  3.32e-3;
                p["B->pi::f_+(0)@BCL2008"] =  0.290;
                p["B->pi::b_+^1@BCL2008"]  = -1.930;
                p["B->pi::b_+^2@BCL2008"]  = -0.441;
                p["mass::B_d"]             =  5.2796;
                p["mass::pi^+"]            =  1.3957e-1;

                Options oo
                {
                    { "model",        "CKM" },
                    { "form-factors", "BCL2008" },
                    { "U",            "u"       },
                    { "q",            "d"       },
                    { "l",            "e"       },
                    { "I",            "1"       }
                };

                BToPseudoscalarLeptonNeutrino d(p, oo);

                const double eps = 1e-8;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 0.01,  2.00), 1.44047e-05, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 2.00,  4.00), 1.43046e-05, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 4.00,  6.00), 1.40803e-05, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 6.00,  8.00), 1.37941e-05, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 8.00, 10.00), 1.34323e-05, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(10.00, 12.00), 1.29770e-05, eps);

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 0.01, 12.00), 8.29930e-5,  eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 0.01, 25.00), 1.43035e-4,  eps);
            }

            // Consistency check for R_pi
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_ub)"]        =  3.32e-3;
                p["B->pi::f_+(0)@BCL2008"] =  0.290;
                p["B->pi::b_+^1@BCL2008"]  = -1.930;
                p["B->pi::b_+^2@BCL2008"]  = -0.441;
                p["mass::B_d"]             =  5.2796;
                p["mass::pi^+"]            =  1.3957e-1;

                Options oo
                {
                    { "model",        "CKM" },
                    { "form-factors", "BCL2008" },
                    { "U",            "u"       },
                    { "q",            "d"       },
                    { "l",            "tau"     },
                    { "I",            "1"     }
                };
                BToPseudoscalarLeptonNeutrino dtau(p, oo);

                oo.declare("l", "mu");
                BToPseudoscalarLeptonNeutrino dmu(p, oo);

                oo =
                {
                    { "model",        "CKM" },
                    { "form-factors", "BCL2008" },
                    { "U",            "u"       },
                    { "q",            "d"       },
                    { "I",            "1"       }
                };
                Kinematics k
                {
                    { "q2_mu_min",   0.011 }, { "q2_mu_max",  10.00 },
                    { "q2_tau_min",  3.154 }, { "q2_tau_max", 10.00 },
                };

                auto obs_Rpi  = Observable::make("B->pilnu::R_pi",   p, k, oo);
                auto obs_Rpip = Observable::make("B->pilnu::R_pi_p", p, k, oo);
                auto obs_Rpi0 = Observable::make("B->pilnu::R_pi_0", p, k, oo);

                const double eps = 1e-5;
                TEST_CHECK_RELATIVE_ERROR(
                    dtau.integrated_branching_ratio(3.154, 10.00) / dmu.integrated_branching_ratio(0.011, 10.00),
                    obs_Rpi->evaluate(),
                    eps
                );
                TEST_CHECK_RELATIVE_ERROR(obs_Rpi->evaluate(),  0.352166, eps);
                TEST_CHECK_RELATIVE_ERROR(obs_Rpip->evaluate(), 0.204647, eps);
                TEST_CHECK_RELATIVE_ERROR(obs_Rpi0->evaluate(), 0.147519, eps);
            }

            // Consistency check for isospin
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_ub)"]        =  3.32e-3;
                p["B->pi::f_+(0)@BCL2008"] =  0.290;
                p["B->pi::b_+^1@BCL2008"]  = -1.930;
                p["B->pi::b_+^2@BCL2008"]  = -0.441;
                p["mass::B_u"]             =  5.2793;
                p["mass::pi^0"]            =  1.3498e-1;
                p["life_time::B_u"]        =  1.519e-12;

                Options oo
                {
                    { "model",        "CKM" },
                    { "form-factors", "BCL2008" },
                    { "U",            "u"       },
                    { "q",            "u"       },
                    { "l",            "e"       },
                    { "I",            "1"       }
                };

                BToPseudoscalarLeptonNeutrino d(p, oo);

                const double eps = 1e-9;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 0.01,  2.00), 1.44047e-05 / 2., eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 2.00,  4.00), 1.43046e-05 / 2., eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 4.00,  6.00), 1.40803e-05 / 2., eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 6.00,  8.00), 1.37941e-05 / 2., eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 8.00, 10.00), 1.34323e-05 / 2., eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(10.00, 12.00), 1.29770e-05 / 2., eps);

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 0.01, 12.00), 8.29930e-5 / 2.,  eps);
            }
        }
} b_to_pi_l_nu_test;
