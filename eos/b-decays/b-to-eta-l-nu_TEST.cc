/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Méril Reboud
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

class BToEtaLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToEtaLeptonNeutrinoTest() :
            TestCase("b_to_eta_l_nu_test")
        {
        }

        virtual void run() const
        {
            // Using IKMvD2014 inputs for V_ub and form factors,
            // from the combined fit to B->pilnu data and LCSR.
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_ub)"]        =  3.32e-3;
                p["mass::B_u"]             =  5.27934;
                p["mass::eta"]             =  0.54786;

                p["B->eta::alpha^f+_0@BSZ2015"] = 1.0;

                Options oo
                {
                    { "model"_ok,        "CKM"     },
                    { "form-factors"_ok, "BSZ2015" },
                    { "P"_ok,            "eta"     },
                    { "q"_ok,            "u"       },
                    { "l"_ok,            "e"       },
                };

                BToPseudoscalarLeptonNeutrino d(p, oo);

                const double eps = 1e-8;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio( 0.01,  2.00), 1.72975e-04, eps);
            }
        }
} b_to_eta_l_nu_test;
