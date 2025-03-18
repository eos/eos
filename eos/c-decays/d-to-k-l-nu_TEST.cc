/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
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
#include <eos/c-decays/d-to-psd-l-nu.hh>
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

class DToKLeptonNeutrinoTest :
    public TestCase
{
    public:
        DToKLeptonNeutrinoTest() :
            TestCase("d_to_k_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["CKM::abs(V_cs)"] =  0.9734;
                p["mass::D_u"]      =  1.86483;
                p["mass::K_u"]      =  0.493677;

                p["D->K::alpha^f+_0@BSZ2015"] = 1.;
                p["D->K::alpha^f+_1@BSZ2015"] = 0.;
                p["D->K::alpha^f+_2@BSZ2015"] = 0.;
                p["D->K::alpha^f0_1@BSZ2015"] = 0.;
                p["D->K::alpha^f0_2@BSZ2015"] = 0.;
                p["D->K::alpha^fT_0@BSZ2015"] = 1.;
                p["D->K::alpha^fT_1@BSZ2015"] = 0.;
                p["D->K::alpha^fT_2@BSZ2015"] = 0.;

                Options oo
                {
                    { "model"_ok,        "CKM"     },
                    { "form-factors"_ok, "BSZ2015" },
                    { "P"_ok,            "K"       },
                    { "q"_ok,            "u"       },
                    { "l"_ok,            "e"       },
                };

                DToPseudoscalarLeptonNeutrino d(p, oo);

                const double eps = 1e-5;

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 0.01,  1.  ), 0.04916404, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio( 0.  ,  1.88), 0.06059538, eps);
            }
        }
} d_to_k_l_nu_test;
