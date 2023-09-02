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
#include <eos/b-decays/b-to-l-nu.hh>
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

class BToLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToLeptonNeutrinoTest() :
            TestCase("b_to_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_ub)"]      = 1.000;
                p["CKM::abs(V_cb)"]      = 2.000;
                p["mass::B_u"]           = 2.000;
                p["mass::B_c"]           = 3.000;
                p["decay-constant::B_u"] = 1.000;
                p["decay-constant::B_c"] = 1.500;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::B_u"]      = 1.000;
                p["life_time::B_c"]      = 0.100;
                p["ubenue::Re{cVL}"]     = 1.000;

                Options oo
                {
                    { "model", "CKM" },
                    { "l",     "e"       },
                    { "q",     "u"       }
                };

                BToLeptonNeutrino d(p, oo);

                const double eps = 1e-6;

                TEST_CHECK_NEARLY_EQUAL(
                        std::pow(1.0066 * 1., 2) * std::pow((1.-std::pow(1./2., 2)), 2) * 2. / (8. * M_PI),
                        d.branching_ratio(),
                        eps);
            }

            {
                Parameters p = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_ub)"]      = 1.000;
                p["CKM::abs(V_cb)"]      = 2.000;
                p["mass::B_u"]           = 2.000;
                p["mass::B_c"]           = 3.000;
                p["decay-constant::B_u"] = 1.000;
                p["decay-constant::B_c"] = 1.500;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::B_u"]      = 1.000;
                p["life_time::B_c"]      = 0.100;
                p["ubenue::Re{cVL}"]     = 1.000;

                Options oo
                {
                    { "model", "CKM" },
                    { "l",     "e"       },
                    { "q",     "c"       }
                };

                BToLeptonNeutrino d(p, oo);

                const double eps = 1e-12;

                TEST_CHECK_NEARLY_EQUAL(
                        std::pow(1.0066 * 1., 2) * 3. / (8. * M_PI) * 0.1
                        * std::pow(1.5 * 2. * (1.-std::pow(1./3., 2)), 2),
                        d.branching_ratio(),
                        eps);
            }

            {
                Parameters p = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_ub)"]      = 1.000;
                p["mass::B_u"]           = 2.000;
                p["decay-constant::B_u"] = 1.000;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::B_u"]      = 1.000;
                p["ubenue::Re{cVL}"]     = 1.000;
                p["ubenue::Re{cVR}"]     = 0.500;
                p["ubenue::Re{cSL}"]     = 0.000;
                p["ubenue::Re{cSR}"]     = 0.000;


                Options oo
                {
                    { "model", "WET" },
                    { "l",     "e"          },
                    { "q",     "u"          }
                };

                BToLeptonNeutrino d(p, oo);

                const double eps = 1e-12;

                // eta factor corrections not yet implemented
                TEST_CHECK_NEARLY_EQUAL(
                        std::pow(0.5, 2) * std::pow((1.-std::pow(1./2., 2)), 2) * 2. / (8. * M_PI),
                        d.branching_ratio(),
                        eps);
            }
        }
} b_to_l_nu_test;
