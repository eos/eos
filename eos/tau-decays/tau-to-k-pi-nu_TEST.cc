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
                p["0->Kpi::M_(+,0)@KSvD2025"]     = 0.890;
                p["0->Kpi::Gamma_(+,0)@KSvD2025"] = 0.026;
                p["0->Kpi::M_(0,0)@KSvD2025"]     = 0.680;
                p["0->Kpi::Gamma_(0,0)@KSvD2025"] = 0.300;
                p["0->Kpi::M_(+,1)@KSvD2025"]     = 1.368;
                p["0->Kpi::Gamma_(+,1)@KSvD2025"] = 0.106;
                p["0->Kpi::M_(0,1)@KSvD2025"]     = 1.431;
                p["0->Kpi::Gamma_(0,1)@KSvD2025"] = 0.110;
                p["0->Kpi::b_+^1@KSvD2025"]       = 0.05;
                p["0->Kpi::b_0^1@KSvD2025"]       = 0.1;

                TauToKPiNeutrino d(p, Options{});

                const double eps = 1e-8;

                TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.0), 0.01475307, eps);
                TEST_CHECK_NEARLY_EQUAL(d.branching_ratio(0.4060, 3.1574802249), 0.05127342, eps);
            }
        }
} tau_to_k_pi_nu_test;
