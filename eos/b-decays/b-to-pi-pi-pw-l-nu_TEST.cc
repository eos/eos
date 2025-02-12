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
                p["B->pipi::chi_1p_a@HKvT2025"]  = 5.742e-4;
                p["B->pipi::chi_1m_v@HKvT2025"]  = 5.742e-4;
                p["pipi->pipi::Gamman0_K@HKvT2025"] =  5.78022130e-01;
                p["B->pipi::a^F1_1_1_0_0@HKvT2025"]   = -1.66817478e-02;
                p["B->pipi::a^F1_1_1_1_0@HKvT2025"]   = -4.85913680e-02;
                p["B->pipi::a^F1_1_1_0_1@HKvT2025"]   =  1.12230134e-04;
                p["B->pipi::a^f_1_1_0_0@HKvT2025"]    = -3.99693932e-02;
                p["B->pipi::a^f_1_1_1_0@HKvT2025"]    = -1.46214509e-02;
                p["B->pipi::a^f_1_1_0_1@HKvT2025"]    = -1.28391535e-02;
                p["B->pipi::a^g_1_1_0_0@HKvT2025"]    =  2.39237676e-02;
                p["B->pipi::a^g_1_1_1_0@HKvT2025"]    =  2.71707816e-02;
                p["B->pipi::a^g_1_1_0_1@HKvT2025"]    =  2.31612080e-04;
                p["B->pipi::a^F1_0_2_0_0@HKvT2025"]   =  6.59969573e-03;
                p["B->pipi::a^F1_0_2_1_0@HKvT2025"]   =  1.63960580e-02;
                p["B->pipi::a^F1_0_2_0_1@HKvT2025"]   = -1.15792204e-03;
                p["B->pipi::a^f_0_2_0_0@HKvT2025"]    =  8.43339186e-03;
                p["B->pipi::a^f_0_2_1_0@HKvT2025"]    =  8.56825464e-03;
                p["B->pipi::a^f_0_2_0_1@HKvT2025"]    = -3.41325494e-03;
                p["B->pipi::a^g_0_2_0_0@HKvT2025"]    =  3.33162140e-02;
                p["B->pipi::a^g_0_2_1_0@HKvT2025"]    =  4.38104301e-02;
                p["B->pipi::a^g_0_2_0_1@HKvT2025"]    = -3.15353047e-03;
                p["B->pipi::a^F1_0_0_0_0@HKvT2025"]   = -6.13947145e-03;
                p["B->pipi::a^F1_0_0_1_0@HKvT2025"]   = -1.69777024e-02;
                p["B->pipi::a^F1_0_0_0_1@HKvT2025"]   = -1.63231083e-02;

                Options oo
                {
                    { "model"_ok,        "CKM" },
                    { "form-factors"_ok, "HKvT2025" },
                    { "scattering-amplitudes"_ok, "HKvT2025" },
                    { "integration-points"_ok, "4096" },
                    { "U"_ok,            "u"       },
                    { "q"_ok,            "u"       },
                    { "l"_ok,            "e"       },
                    { "I1"_ok,           "1"       },
                    { "I2"_ok,           "1"       },
                    { "C"_ok,           "+-"       },
                    { "I"_ok,           "0|1"      }
                };

                oo.declare("L"_ok, "S");
                BToPPLeptonNeutrino testS(p, oo);
                oo.declare("L"_ok, "P");
                BToPPLeptonNeutrino testP(p, oo);
                oo.declare("L"_ok, "D");
                BToPPLeptonNeutrino testD(p, oo);

                TEST_CHECK_NEARLY_EQUAL(testS.fully_integrated_branching_ratio(),    1.59995e-05, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(testS.q2_integrated_branching_ratio(0.27914, 1.02), 3.97362e-06, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(testP.fully_integrated_branching_ratio(),    2.18146e-04, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(testP.q2_integrated_branching_ratio(0.27914, 1.02), 1.84595e-04, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(testD.fully_integrated_branching_ratio(),    2.53129e-05, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(testD.q2_integrated_branching_ratio(0.27914, 1.02), 1.49765e-06, 1e-8);

            }
        }
} b_to_pi_pi_pw_l_nu_test;
