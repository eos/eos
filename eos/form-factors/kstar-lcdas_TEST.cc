/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Stefan Meiser
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
#include <eos/form-factors/k-star-lcdas.hh>

#include <eos/models/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class KStarLCDAsTest :
    public TestCase
{
    public:
        KStarLCDAsTest() :
            TestCase("kstar_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"] = 0.1176;
            p["K^*::a1para@1GeV"] = 0.03;
            p["K^*::a2para@1GeV"] = 0.11;
            p["K^*::a3para@1GeV"] = 0.21;
            p["K^*::a4para@1GeV"] = 0.14;
            p["K^*::a1perp@1GeV"] = 0.04;
            p["K^*::a2perp@1GeV"] = 0.10;
            p["K^*::a3perp@1GeV"] = 0.15;
            p["K^*::a4perp@1GeV"] = 0.19;
            p["K^*::fperp@1GeV"]  = 0.159;

            /* Diagnostics */
            {
                AntiKStarLCDAs kstar(p, Options{ });
                Diagnostics diagnostics = kstar.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94850, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92874, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91708, 1e-5), // c_rge(mu = 4.0 GeV)
                    std::make_pair(+0.90893, 1e-5), // c_rge(mu = 5.0 GeV)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* Twist 2 */
            {
                AntiKStarLCDAs kstar(p, Options{ });

                // coefficients at mu = 1.0 GeV, and 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.03000,  kstar.a1para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.02486,  kstar.a1para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.11000,  kstar.a2para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.08200,  kstar.a2para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.21000,  kstar.a3para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14521,  kstar.a3para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14000,  kstar.a4para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.09128,  kstar.a4para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.04000,  kstar.a1perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.03238,  kstar.a1perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.10000,  kstar.a2perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.07368,  kstar.a2perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15000,  kstar.a3perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.10299,  kstar.a3perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.19000,  kstar.a4perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.12330,  kstar.a4perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15900,  kstar.fperp(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14637,  kstar.fperp(2.0),    eps);

                // phipara LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL( 0.45242, kstar.phipara(0.1, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.43819, kstar.phipara(0.3, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.64625, kstar.phipara(0.5, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.53401, kstar.phipara(0.7, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.20151, kstar.phipara(0.9, 1.0),  eps);

                // phipara LCDA at various u values for mu = 2.0
                TEST_CHECK_NEARLY_EQUAL( 0.48111, kstar.phipara(0.1, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.38391, kstar.phipara(0.3, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.57223, kstar.phipara(0.5, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.77114, kstar.phipara(0.7, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.00975, kstar.phipara(0.9, 2.0),  eps);

                // phiperp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL( 0.55003, kstar.phiperp(0.1, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.20175, kstar.phiperp(0.3, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.80938, kstar.phiperp(0.5, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.61207, kstar.phiperp(0.7, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.13323, kstar.phiperp(0.9, 1.0),  eps);

                // phiperp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL( 0.54481, kstar.phiperp(0.1, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.22256, kstar.phiperp(0.3, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.68102, kstar.phiperp(0.5, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.83252, kstar.phiperp(0.7, 2.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.95797, kstar.phiperp(0.9, 2.0),  eps);
            }
        }
} kstar_lcdas_test;
