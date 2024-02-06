/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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
#include <eos/rare-b-decays/b-to-psd-nu-nu.hh>
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

class BToKDineutrinoTest :
    public TestCase
{
    public:
        BToKDineutrinoTest() :
            TestCase("b_to_k_nu_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_tb)"]           =  1.00;
                p["CKM::abs(V_ts)"]           =  4.00e-2;
                p["B->K::alpha^f+_0@BSZ2015"] = +3.2909e-01;
                p["B->K::alpha^f+_1@BSZ2015"] = -8.6695e-01;
                p["B->K::alpha^f+_2@BSZ2015"] = +6.0957e-03;
                p["mass::B_u"]                =  5.2796;
                p["mass::K_u"]                =  4.9368e-01;

                Options oo
                {
                    { "model",        "CKM"     },
                    { "form-factors", "BSZ2015" },
                    { "D",            "s"       },
                    { "q",            "u"       },
                    { "I",            "1/2"     }
                };

                BToPseudoscalarDineutrino d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 0.00,  8.00), 2.05845e-06, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 8.00, 16.00), 1.68211e-06, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(16.00, 22.90), 0.59978e-06, eps);

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 0.00, 22.90), 4.34034e-06, eps);
            }
        }
} b_to_k_nu_nu_test;
