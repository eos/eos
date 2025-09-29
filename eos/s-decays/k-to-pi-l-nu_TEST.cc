/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Matthew Kirk
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
#include <eos/s-decays/k-to-pi-l-nu.hh>
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

class KToPiLeptonNeutrinoTest : public TestCase
{
    public:
        KToPiLeptonNeutrinoTest() :
            TestCase("k_to_pi_l_nu_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p                      = Parameters::Defaults();
                // Parameter values taken to roughly reproduce tau->Kpinu differential rate
                p["0->Kpi::t_0@KSvD2025"]         = -4.0;
                p["0->Kpi::M_(+,0)@KSvD2025"]     = 0.892;
                p["0->Kpi::Gamma_(+,0)@KSvD2025"] = 0.046;
                p["0->Kpi::M_(0,0)@KSvD2025"]     = 0.680;
                p["0->Kpi::Gamma_(0,0)@KSvD2025"] = 0.600;
                p["0->Kpi::M_(+,1)@KSvD2025"]     = 1.368;
                p["0->Kpi::Gamma_(+,1)@KSvD2025"] = 0.212;
                p["0->Kpi::b_+^1@KSvD2025"]       = 0.0088;
                p["0->Kpi::b_+^2@KSvD2025"]       = -0.0052;
                p["0->Kpi::b_+^3@KSvD2025"]       = -0.0273;
                p["0->Kpi::b_+^4@KSvD2025"]       = -0.0094;
                p["0->Kpi::b_0^1@KSvD2025"]       = 0.2653;
                p["0->Kpi::b_0^2@KSvD2025"]       = -0.0588;
                p["0->Kpi::b_0^3@KSvD2025"]       = -0.4228;
                p["0->Kpi::b_0^4@KSvD2025"]       = -0.1646;

                const double eps = 1e-8;

                // K- -> pi0 mu- nubar
                KToPiLeptonNeutrino d(p,
                                      Options{
                                          { "n-resonances-1m"_ok, "2" },
                                          { "n-resonances-0m"_ok, "1" }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 3.3314331e-2, eps); // Compare with PDG 3.24e-2

                // K- -> pi0 e- nubar
                d = KToPiLeptonNeutrino(p,
                                        Options{
                                            {               "l"_ok, "e" },
                                            { "n-resonances-1m"_ok, "2" },
                                            { "n-resonances-0m"_ok, "1" }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 4.8686715e-2, eps); // Compare with PDG 5.07e-2

                // K_S -> pi+ mu- nubar
                d = KToPiLeptonNeutrino(p,
                                        Options{
                                            {               "K"_ok, "K_S" },
                                            { "n-resonances-1m"_ok,   "2" },
                                            { "n-resonances-0m"_ok,   "1" }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 2.43594e-4, eps); // Compare with PDG 2.28e-4

                // K_L -> pi+ mu- nubar
                d = KToPiLeptonNeutrino(p,
                                        Options{
                                            {               "K"_ok, "K_L" },
                                            { "n-resonances-1m"_ok,   "2" },
                                            { "n-resonances-0m"_ok,   "1" }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 0.13918086, eps); // Compare with PDG 0.135
            }
        }
} k_to_pi_lepton_neutrino_test;
