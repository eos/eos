/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
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
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <limits>

using namespace test;
using namespace eos;

class BToPiPiLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToPiPiLeptonNeutrinoTest() :
            TestCase("b_to_pi_pi_l_nu_test")
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
                p["mass::d(2GeV)"]         =  0.0048;
                p["mass::u(2GeV)"]         =  0.0032;
                p["pi::a2@1GeV"]           = +0.16;
                p["pi::a4@1GeV"]           = +0.04;
                p["B->pipi::mu@BFvD2016"]  = +1.5;
                Options oo;
                oo.declare("model"_ok, "CKM");
                oo.declare("form-factors"_ok, "BFvD2016");

                BToPiPiLeptonNeutrino d(p, oo);

                const double eps = 2e-3;

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(0.02, 0.95, 18.60, 26.40, -1.0, +1.0), 5.7200653e-13, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(0.02, 1.95, 15.00, 26.40, -1.0, +1.0), 8.8931904e-12, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(0.02, 0.95, 18.60, 26.40), -0.22715,       eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(0.02, 1.95, 15.00, 26.40), -0.10447,       eps);
            }
        }
} b_to_pi_pi_l_nu_test;
