/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <src/rare-b-decays/inclusive-b-to-s-gamma.hh>

#include <cmath>
#include <limits>

#include <iostream>

using namespace test;
using namespace eos;

class BToXsGammaLargeRecoilTest :
    public TestCase
{
    public:
        BToXsGammaLargeRecoilTest() :
            TestCase("b_to_x_s_gamma_minimal_test")
        {
        }

        virtual void run() const
        {
            /* Minimal */

            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["Abs{c7}"] = +0.331;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;

                Options oo;

                BToXsGamma<Minimal> d(p, oo);

                const double eps = 1e-9;

                TEST_CHECK_NEARLY_EQUAL(3.84894e-4, d.integrated_branching_ratio(), eps);
            }

            // Zero test
            {
                Parameters p = Parameters::Defaults();
                p["Abs{c7}"] = +0.3;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;

                Options oo;

                BToXsGamma<Minimal> d(p, oo);

                const double eps = 1e-9;

                TEST_CHECK_NEARLY_EQUAL(3.15e-4, d.integrated_branching_ratio(), eps);
            }
        }
} b_to_x_s_gamma_large_recoil_test;
