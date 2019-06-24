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
#include <eos/utils/complex.hh>
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
                    { "model",        "CKMScan" },
                    { "form-factors", "BCL2008" },
                    { "U",            "u"       },
                    { "q",            "d"       },
                    { "l",            "e"       }
                };

                BToPseudoscalarLeptonNeutrino d(p, oo);

                const double eps = 1e-8;

                TEST_CHECK_NEARLY_EQUAL(1.44047e-05, d.integrated_branching_ratio( 0.01,  2.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.43046e-05, d.integrated_branching_ratio( 2.00,  4.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.40803e-05, d.integrated_branching_ratio( 4.00,  6.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.37941e-05, d.integrated_branching_ratio( 6.00,  8.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.34323e-05, d.integrated_branching_ratio( 8.00, 10.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.29770e-05, d.integrated_branching_ratio(10.00, 12.00), eps);

                TEST_CHECK_NEARLY_EQUAL(8.29930e-5,  d.integrated_branching_ratio( 0.01, 12.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.43035e-4,  d.integrated_branching_ratio( 0.01, 25.00), eps);
            }
        }
} b_to_pi_l_nu_test;
