/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Danny van Dyk
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
#include <eos/b-decays/b-to-pi-l-nu.hh>
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

                Options oo;
                oo.set("model", "CKMScan");
                oo.set("form-factors", "BCL2008");

                BToPiLeptonNeutrino d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 0.0,  2.0), 1.42880e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 2.0,  4.0), 1.41174e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 4.0,  6.0), 1.38961e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 6.0,  8.0), 1.36135e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 8.0, 10.0), 1.32564e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(10.0, 12.0), 1.28071e-5, eps);
            }
        }
} b_to_pi_l_nu_test;
