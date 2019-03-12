/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Ahmet Kokulu
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
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
#include <eos/utils/complex.hh>
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

class LambdaBToLambdaCLeptonNeutrinoTest :
    public TestCase
{
    public:
        LambdaBToLambdaCLeptonNeutrinoTest() :
            TestCase("lambdab_to_lambdac_l_nu_test")
        {
        }

        virtual void run() const
        {

            // tests for SM observables, Re{cVL}=1.0 in the SM and all other couplings are zero
            {
                Parameters p1 = Parameters::Defaults();
                // the parameters are fixed as EOS default values

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "DKMR2017");
                oo.set("l", "mu");

                LambdaBToLambdaCLeptonNeutrino d(p1, oo);

                const double eps = 1e-3;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011, 11.1), -0.2028, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_hadronic(0.011, 11.1),  0.3274, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_combined(0.011, 11.1), -0.1177, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_fzero(0.011, 11.1),  0.5869, eps);
            }
            // tests for NP observables
            {
                Parameters p2 = Parameters::Defaults();
                // the rest of input is fixed to default EOS values
                p2["b->cmunumu::Re{cVL}"]   =  1.0;
                p2["b->cmunumu::Im{cVL}"]   = -1.0;
                p2["b->cmunumu::Re{cVR}"]   =  2.0;
                p2["b->cmunumu::Im{cVR}"]   = -2.0;
                p2["b->cmunumu::Re{cSL}"]   =  3.0;
                p2["b->cmunumu::Im{cSL}"]   = -3.0;
                p2["b->cmunumu::Re{cSR}"]   =  4.0;
                p2["b->cmunumu::Im{cSR}"]   = -4.0;
                p2["b->cmunumu::Re{cT}"]    =  0.0;
                p2["b->cmunumu::Im{cT}"]    =  0.0;
                // fix the scale
                p2["mu"]                    =  4.18;
                p2["mass::b(MSbar)"]        =  4.18;
                p2["mass::c"]               =  1.275;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "DKMR2017");
                oo.set("l", "mu");

                LambdaBToLambdaCLeptonNeutrino d(p2, oo);

                const double eps = 1e-2;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011, 11.1),   0.0465, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_hadronic(0.011, 11.1),  -0.0179, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_combined(0.011, 11.1),  -0.0150, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_fzero(0.011, 11.1),  0.4016, eps);
            }
        }
} lambdab_to_lambdac_l_nu_test;
