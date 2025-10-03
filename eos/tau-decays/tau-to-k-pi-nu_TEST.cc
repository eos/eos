/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/tau-decays/tau-to-k-pi-nu.hh>
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

class TauToKPiNeutrinoTest : public TestCase
{
    public:
        TauToKPiNeutrinoTest() :
            TestCase("tau_to_k_pi_nu_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p                      = Parameters::Defaults();
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

                // tau -> K- pi0 nu
                TauToKPiNeutrino d(p,
                                   Options{
                                       { "n-resonances-1m"_ok, "2" },
                                       { "n-resonances-0m"_ok, "1" }
                });

                TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.0), 0.0017916075697680681, eps);
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 0.004430198678120239, eps); // Compare with PDG 0.433%

                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_q2(1.0), 0.0017916075697680681 / 0.004430198678120239, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_q2(0.4060, 3.1571), 0.3631445295018548, eps); // Compare with 1 / (3.1571 - 0.4060)

                // tau -> K_S pi- nu
                d = TauToKPiNeutrino(p,
                                     Options{
                                         {               "K"_ok, "K_S" },
                                         { "n-resonances-1m"_ok,   "2" },
                                         { "n-resonances-0m"_ok,   "1" }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 4.25288e-3, eps); // Compare with [Belle:2014B] 4.16e-3

                // tau -> K_L pi- nu
                d = TauToKPiNeutrino(p,
                                     Options{
                                         {               "K"_ok, "K_L" },
                                         { "n-resonances-1m"_ok,   "2" },
                                         { "n-resonances-0m"_ok,   "1" }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 4.25288e-3, eps);
            }
        }
} tau_to_k_pi_nu_test;
