/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/s-decays/k-to-l-nu.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class KToLeptonNeutrinoTest : public TestCase
{
    public:
        KToLeptonNeutrinoTest() :
            TestCase("k_to_l_nu_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p             = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_us)"]      = 1.000;
                p["mass::K_u"]           = 2.000;
                p["decay-constant::K_u"] = 1.000;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::K_u"]      = 1.000;

                Options oo{
                    { "model", "CKM" },
                    {     "l",   "e" },
                };

                KToLeptonNeutrino d(p, oo);

                const double eps = 1e-12;

                TEST_CHECK_NEARLY_EQUAL(std::pow(1.009653 * 1., 2) * std::pow((1. - std::pow(1. / 2., 2)), 2) * 2. / (8. * M_PI), d.branching_ratio(), eps);
            }

            {
                Parameters p             = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_us)"]      = 2.000;
                p["mass::K_u"]           = 3.000;
                p["decay-constant::K_u"] = 1.500;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::K_u"]      = 0.100;

                Options oo{
                    { "model", "CKM" },
                    {     "l",   "e" },
                };

                KToLeptonNeutrino d(p, oo);

                const double eps = 1e-12;

                TEST_CHECK_NEARLY_EQUAL(std::pow(1.009653 * 1., 2) * 3. / (8. * M_PI) * 0.1 * std::pow(1.5 * 2. * (1. - std::pow(1. / 3., 2)), 2), d.branching_ratio(), eps);
            }

            {
                Parameters p             = Parameters::Defaults();
                p["WET::G_Fermi"]        = 1.000;
                p["CKM::abs(V_us)"]      = 1.000;
                p["mass::K_u"]           = 2.000;
                p["decay-constant::K_u"] = 1.000;
                p["mass::e"]             = 1.000;
                p["QM::hbar"]            = 1.000;
                p["life_time::K_u"]      = 1.000;
                p["usenue::Re{cVL}"]     = 1.000;
                p["usenue::Re{cVR}"]     = 0.500;
                p["usenue::Re{cSL}"]     = 0.000;
                p["usenue::Re{cSR}"]     = 0.000;


                Options oo{
                    { "model", "WET" },
                    {     "l",   "e" },
                };

                KToLeptonNeutrino d(p, oo);

                const double eps = 1e-12;

                // eta factor corrections not yet implemented
                TEST_CHECK_NEARLY_EQUAL(std::pow(0.5, 2) * std::pow((1. - std::pow(1. / 2., 2)), 2) * 2. / (8. * M_PI), d.branching_ratio(), eps);
            }
        }
} k_to_l_nu_test;
