/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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
#include <eos/b-decays/b-to-psd-psd-l-nu.hh>
#include <eos/maths/complex.hh>
#include <iostream>

using namespace test;
using namespace eos;
using namespace std;

class BToPiPiPWLeptonNeutrinoTest :
    public TestCase
{
    public:
        BToPiPiPWLeptonNeutrinoTest() :
            TestCase("b_to_pi_pi_pw_l_nu_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();
                p["pipi->pipi::Gamman0_K@HKvT2025"] =  5.78022130e-01;
                p["B->pipi::a^F1_2_1_0@HKvT2025"]   = -1.66817478e-02;
                p["B->pipi::a^F1_2_1_1@HKvT2025"]   = -4.85913680e-02;
                p["B->pipi::a^F1_2_1_5@HKvT2025"]   =  1.12230134e-04;
                p["B->pipi::a^f_2_1_0@HKvT2025"]    = -3.99693932e-02;
                p["B->pipi::a^f_2_1_1@HKvT2025"]    = -1.46214509e-02;
                p["B->pipi::a^f_2_1_5@HKvT2025"]    = -1.28391535e-02;
                p["B->pipi::a^g_2_1_0@HKvT2025"]    =  2.39237676e-02;
                p["B->pipi::a^g_2_1_1@HKvT2025"]    =  2.71707816e-02;
                p["B->pipi::a^g_2_1_5@HKvT2025"]    =  2.31612080e-04;
                p["B->pipi::a^F1_1_2_0@HKvT2025"]   =  6.59969573e-03;
                p["B->pipi::a^F1_1_2_1@HKvT2025"]   =  1.63960580e-02;
                p["B->pipi::a^F1_1_2_5@HKvT2025"]   = -1.15792204e-03;
                p["B->pipi::a^f_1_2_0@HKvT2025"]    =  8.43339186e-03;
                p["B->pipi::a^f_1_2_1@HKvT2025"]    =  8.56825464e-03;
                p["B->pipi::a^f_1_2_5@HKvT2025"]    = -3.41325494e-03;
                p["B->pipi::a^g_1_2_0@HKvT2025"]    =  3.33162140e-02;
                p["B->pipi::a^g_1_2_1@HKvT2025"]    =  4.38104301e-02;
                p["B->pipi::a^g_1_2_5@HKvT2025"]    = -3.15353047e-03;
                p["B->pipi::a^F1_1_0_0@HKvT2025"]   = -6.13947145e-03;
                p["B->pipi::a^F1_1_0_1@HKvT2025"]   = -1.69777024e-02;
                p["B->pipi::a^F1_1_0_5@HKvT2025"]   = -1.63231083e-02;

                Options oo
                {
                    { "model",        "CKM" },
                    { "form-factors", "HKvT2025" },
                    { "scattering-amplitudes", "HKvT2025" },
                    { "integration-points", "4096" },
                    { "U",            "u"       },
                    { "q",            "u"       },
                    { "l",            "e"       },
                    { "I1",           "1"       },
                    { "I2",           "1"       },
                    { "C",           "+-"       }
                };

                BToPPLeptonNeutrino test(p, oo);

                TEST_CHECK_NEARLY_EQUAL(test.fully_integrated_branching_ratio_S(),    1.59995e-05, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(test.fully_integrated_branching_ratio_P(),    2.18146e-04, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(test.fully_integrated_branching_ratio_D(),    2.53129e-05, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(test.saturation_1_p(),   0.077110, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(test.saturation_1_m(),   0.298979, 1e-5);
                TEST_CHECK_NEARLY_EQUAL(test.q2_integrated_branching_ratio_S(0.27914, 1.02), 3.97362e-06, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(test.q2_integrated_branching_ratio_P(0.27914, 1.02), 1.84595e-04, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(test.q2_integrated_branching_ratio_D(0.27914, 1.02), 1.49765e-06, 1e-8);

            }
        }
} b_to_pi_pi_pw_l_nu_test;
