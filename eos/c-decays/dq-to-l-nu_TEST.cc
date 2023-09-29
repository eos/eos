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
#include <eos/c-decays/dq-to-l-nu.hh>
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

class DqToLeptonNeutrinoTest :
    public TestCase
{
    public:
        DqToLeptonNeutrinoTest() :
            TestCase("dq_to_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_cs)"]      = 1.000;
                p["mass::D_s"]           = 2.000;
                p["decay-constant::D_s"] = 1.000;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::D_s"]      = 1.000;
                p["scnuee::Re{cVL}"]     = 1.000;

                Options oo
                {
                    { "model", "CKM" },
                    { "l",     "e"   },
                    { "q",     "s"   }
                };

                DqToLeptonNeutrino d(p, oo);

                const double eps = 1e-6;

                TEST_CHECK_NEARLY_EQUAL(
                        d.branching_ratio(),
                        std::pow(1.01033 * 1., 2) * std::pow((1.-std::pow(1./2., 2)), 2) * 2. / (8. * M_PI),
                        eps);
            }
        }
} dq_to_l_nu_test;
