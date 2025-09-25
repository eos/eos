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
#include <eos/tau-decays/tau-to-k-nu.hh>
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

class TauToKNeutrinoTest : public TestCase
{
    public:
        TauToKNeutrinoTest() :
            TestCase("tau_to_k_nu_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p = Parameters::Defaults();

                TauToKNeutrino d(p, Options{});

                const double eps = 1e-8;

                TEST_CHECK_NEARLY_EQUAL(d.branching_ratio(), 0.7192811e-2, eps); // Compare with PDG 0.685%
            }
        }
} tau_to_k_nu_test;
