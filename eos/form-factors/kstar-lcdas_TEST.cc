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

class AntiKStarLCDAsTest :
    public TestCase
{
    public:
        AntiKStarLCDAsTest() :
            TestCase("kstar_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"]           = 0.1176;
            p["mass::s(2GeV)"]              = 0.095;
            p["mass::u(2GeV)"]              = 0.0032;
            p["mass::d(2GeV)"]              = 0.0032;   // we use m_ud/2 = m_d = m_u = 3.2 MeV
            p["K^*::a1para@1GeV"]           = 0.03;
            p["K^*::a2para@1GeV"]           = 0.11;
            p["K^*::a3para@1GeV"]           = 0.21;
            p["K^*::a4para@1GeV"]           = 0.14;
            p["K^*::a1perp@1GeV"]           = 0.04;
            p["K^*::a2perp@1GeV"]           = 0.10;
            p["K^*::a3perp@1GeV"]           = 0.15;
            p["K^*::a4perp@1GeV"]           = 0.19;
            p["K^*::fperp@1GeV"]            = 0.159;
            p["K^*::zeta3para@1GeV"]        = 0.023;
            p["K^*::lambda3paratilde@1GeV"] = 0.035;
            p["K^*::omega3paratilde@1GeV"]  = -0.07;
            p["K^*::kappa3para@1GeV"]       = 0.000;
            p["K^*::omega3para@1GeV"]       = 0.1;
            p["K^*::lambda3para@1GeV"]      = -0.008;
            p["K^*::kappa3perp@1GeV"]       = 0.003;
            p["K^*::omega3perp@1GeV"]       = 0.3;
            p["K^*::lambda3perp@1GeV"]      = -0.025;
            p["K^*::fpara"]                 = 0.204;
            p["mass::K_u^*"]                = 0.89166;

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
                TEST_CHECK_NEARLY_EQUAL(kstar.a1para(1.0), 0.03000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a1para(2.0), 0.02486, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2para(1.0), 0.11000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2para(2.0), 0.08200, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3para(1.0), 0.21000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3para(2.0), 0.14521, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4para(1.0), 0.14000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4para(2.0), 0.09128, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a1perp(1.0), 0.04000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a1perp(2.0), 0.03238, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2perp(1.0), 0.10000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2perp(2.0), 0.07368, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3perp(1.0), 0.15000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3perp(2.0), 0.10299, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4perp(1.0), 0.19000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4perp(2.0), 0.12330, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.fperp(1.0),  0.15900, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.fperp(2.0),  0.14818, eps);

                // phipara LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.1, 1.0), 0.45242, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.3, 1.0), 1.43819, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.5, 1.0), 1.64625, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.7, 1.0), 0.53401, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.9, 1.0), 1.20151, eps);

                // phipara LCDA at various u values for mu = 2.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.1, 2.0), 0.48111, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.3, 2.0), 1.38391, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.5, 2.0), 1.57223, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.7, 2.0), 0.77114, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.9, 2.0), 1.00975, eps);

                // phiperp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.1, 1.0), 0.55003, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.3, 1.0), 1.20175, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.5, 1.0), 1.80938, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.7, 1.0), 0.61207, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.9, 1.0), 1.13323, eps);

                // phiperp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.1, 2.0), 0.54481, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.3, 2.0), 1.22256, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.5, 2.0), 1.68102, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.7, 2.0), 0.83252, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.9, 2.0), 0.95797, eps);
            }

            /* Twist 3 */
            {
                AntiKStarLCDAs kstar(p, Options{ });

                // parameters at mu = 1.0, 2.0, 3.0, 4.0, 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(1.0),        0.0230000,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(2.0),        0.0155724,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(3.0),        0.0133572,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(4.0),        0.0121881,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(5.0),        0.0114277,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(1.0), 0.035,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(2.0), 0.0185017,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(3.0), 0.014129 ,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(4.0), 0.0119487,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(5.0), 0.0105829,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(1.0), -0.07,        eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(2.0), -0.0362296,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(3.0), -0.0279061,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(4.0), -0.0238674,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(5.0), -0.0213772,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(1.0),      +0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(2.0),      -0.000882678, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(3.0),      -0.0010668,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(4.0),      -0.00114477,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(5.0),      -0.00118739,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(1.0),       0.1,         eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(2.0),       0.0655062,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(3.0),       0.0552781,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(4.0),       0.0499135,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(5.0),       0.0464411,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(1.0),     -0.008,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(2.0),     -0.00467196,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(3.0),     -0.00377474,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(4.0),     -0.00332192,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(5.0),     -0.00303565,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(1.0),      +0.003,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(2.0),      -0.00109164,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(3.0),      -0.00231601,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(4.0),      -0.00295769,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(5.0),      -0.00337196,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(1.0),       0.3,         eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(2.0),       0.220453,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(3.0),       0.195552,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(4.0),       0.182125,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(5.0),       0.173271,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(1.0),     -0.025,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(2.0),     -0.0156331,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(3.0),     -0.0130251,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(4.0),     -0.0116894,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(5.0),     -0.0108369,   eps);
            }
        }
} anti_kstar_lcdas_test;

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
            p["QCD::alpha_s(MZ)"]           = 0.1176;
            p["mass::s(2GeV)"]              = 0.0032;
            p["mass::u(2GeV)"]              = 0.095;
            p["mass::d(2GeV)"]              = 0.095;
            p["K^*::a1para@1GeV"]           = -0.03;
            p["K^*::a2para@1GeV"]           = 0.11;
            p["K^*::a3para@1GeV"]           = -0.21;
            p["K^*::a4para@1GeV"]           = 0.14;
            p["K^*::a1perp@1GeV"]           = -0.04;
            p["K^*::a2perp@1GeV"]           = 0.10;
            p["K^*::a3perp@1GeV"]           = -0.15;
            p["K^*::a4perp@1GeV"]           = 0.19;
            p["K^*::fperp@1GeV"]            = 0.159;
            p["K^*::zeta3para@1GeV"]        = 0.023;
            p["K^*::lambda3paratilde@1GeV"] = -0.035;
            p["K^*::omega3paratilde@1GeV"]  = -0.07;
            p["K^*::kappa3para@1GeV"]       = -0.000;
            p["K^*::omega3para@1GeV"]       = 0.1;
            p["K^*::lambda3para@1GeV"]      = 0.008;
            p["K^*::kappa3perp@1GeV"]       = -0.003;
            p["K^*::omega3perp@1GeV"]       = 0.3;
            p["K^*::lambda3perp@1GeV"]      = 0.025;
            p["K^*::fpara"]                 = 0.204;
            p["mass::K_u^*"]                = 0.89166;

            /* Diagnostics */
            {
                KStarLCDAs kstar(p, Options{ });
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
                KStarLCDAs kstar(p, Options{ });

                // coefficients at mu = 1.0 GeV, and 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.a1para(1.0), 0.03000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a1para(2.0), 0.02486, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2para(1.0), 0.11000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2para(2.0), 0.08200, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3para(1.0), 0.21000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3para(2.0), 0.14521, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4para(1.0), 0.14000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4para(2.0), 0.09128, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a1perp(1.0), 0.04000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a1perp(2.0), 0.03238, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2perp(1.0), 0.10000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a2perp(2.0), 0.07368, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3perp(1.0), 0.15000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a3perp(2.0), 0.10299, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4perp(1.0), 0.19000, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.a4perp(2.0), 0.12330, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.fperp(1.0),  0.15900, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.fperp(2.0),  0.14818, eps);

                // phipara LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.1, 1.0), 0.45242, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.3, 1.0), 1.43819, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.5, 1.0), 1.64625, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.7, 1.0), 0.53401, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.9, 1.0), 1.20151, eps);

                // phipara LCDA at various u values for mu = 2.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.1, 2.0), 0.48111, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.3, 2.0), 1.38391, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.5, 2.0), 1.57223, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.7, 2.0), 0.77114, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phipara(0.9, 2.0), 1.00975, eps);

                // phiperp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.1, 1.0), 0.55003, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.3, 1.0), 1.20175, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.5, 1.0), 1.80938, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.7, 1.0), 0.61207, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.9, 1.0), 1.13323, eps);

                // phiperp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.1, 2.0), 0.54481, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.3, 2.0), 1.22256, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.5, 2.0), 1.68102, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.7, 2.0), 0.83252, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phiperp(0.9, 2.0), 0.95797, eps);
            }

            /* Twist 3 */
            {
                KStarLCDAs kstar(p, Options{ });

                // parameters at mu = 1.0, 2.0, 3.0, 4.0, 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(1.0),        0.0230000,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(2.0),        0.0155724,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(3.0),        0.0133572,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(4.0),        0.0121881,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta3para(5.0),        0.0114277,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(1.0), 0.035,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(2.0), 0.0185017,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(3.0), 0.014129 ,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(4.0), 0.0119487,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3paratilde(5.0), 0.0105829,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(1.0), -0.07,        eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(2.0), -0.0362296,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(3.0), -0.0279061,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(4.0), -0.0238674,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3paratilde(5.0), -0.0213772,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(1.0),      +0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(2.0),      -0.000882678, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(3.0),      -0.0010668,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(4.0),      -0.00114477,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3para(5.0),      -0.00118739,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(1.0),       0.1,         eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(2.0),       0.0655062,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(3.0),       0.0552781,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(4.0),       0.0499135,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3para(5.0),       0.0464411,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(1.0),     -0.008,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(2.0),     -0.00467196,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(3.0),     -0.00377474,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(4.0),     -0.00332192,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3para(5.0),     -0.00303565,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(1.0),      +0.003,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(2.0),      -0.00109164,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(3.0),      -0.00231601,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(4.0),      -0.00295769,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa3perp(5.0),      -0.00337196,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(1.0),       0.3,         eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(2.0),       0.220453,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(3.0),       0.195552,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(4.0),       0.182125,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega3perp(5.0),       0.173271,    eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(1.0),     -0.025,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(2.0),     -0.0156331,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(3.0),     -0.0130251,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(4.0),     -0.0116894,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.lambda3perp(5.0),     -0.0108369,   eps);
            }
        }
} kstar_lcdas_test;
