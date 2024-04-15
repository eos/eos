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
            p["K^*::zeta4para@1GeV"]        = 0.02;
            p["K^*::omega4paratilde@1GeV"]  = -0.02;
            p["K^*::zeta4perp@1GeV"]        = -0.01;
            p["K^*::zeta4perptilde@1GeV"]   = -0.05;
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

                // two particle LCDAs at scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.1, 1.0), 0.480777, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.3, 1.0), 0.931977, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.5, 1.0), 1.034104, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.7, 1.0), 1.072796, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.9, 1.0), 0.675659, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.1, 1.0), 1.35473,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.3, 1.0), 0.525137, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.5, 1.0), 0.795833, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.7, 1.0), 0.357182, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.9, 1.0), 1.53821,  eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.1, 1.0), 0.499909, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.3, 1.0), 1.04631,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.5, 1.0), 1.14539,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.7, 1.0), 1.14576,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.9, 1.0), 0.730355, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.1, 1.0), 1.10177,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.3, 1.0), 0.844701, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.5, 1.0), 0.883863, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.7, 1.0), 0.739837, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.9, 1.0), 1.32029,  eps);

                // two particle LCDAs at scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.1, 2.0), 0.488303, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.3, 2.0), 1.00108,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.5, 2.0), 1.13694,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.7, 2.0), 1.10367,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.9, 2.0), 0.629763, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.1, 2.0), 1.46711,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.3, 2.0), 0.509224, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.5, 2.0), 0.612806, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.7, 2.0), 0.422571, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.9, 2.0), 1.62548,  eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.1, 2.0), 0.507961, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.3, 2.0), 1.10493,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.5, 2.0), 1.25996,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.7, 2.0), 1.18558,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.9, 2.0), 0.662752, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.1, 2.0), 1.13684,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.3, 2.0), 0.84257,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.5, 2.0), 0.835106, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.7, 2.0), 0.79929,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.9, 2.0), 1.29942,  eps);

                // three particle LCDA scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.1, 0.8, 1.0), -0.0239616, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.3, 0.6, 1.0), -0.0964224, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.5, 0.4, 1.0), -0.112896,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.7, 0.2, 1.0), -0.0540288, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.9, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.1, 0.6, 1.0),  0.0590976, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.3, 0.4, 1.0),  0.0041472, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.5, 0.2, 1.0), -0.029376,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.7, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.1, 0.4, 1.0),  0.117504,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.3, 0.2, 1.0),  0.057024,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.5, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.1, 0.2, 1.0),  0.0669312, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.3, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.0, 0.1, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.1, 0.0, 1.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.1, 0.8, 1.0), -0.156672, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.3, 0.6, 1.0), -0.101088, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.5, 0.4, 1.0),  0.04608,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.7, 0.2, 1.0),  0.058464, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.9, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.1, 0.6, 1.0), -0.046656, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.3, 0.4, 1.0),  0.15552,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.5, 0.2, 1.0),  0.15552,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.7, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.1, 0.4, 1.0),  0.12672,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.3, 0.2, 1.0),  0.18576,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.5, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.1, 0.2, 1.0),  0.1008,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.3, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.0, 0.1, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.1, 0.0, 1.0),  0.0,      eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.1, 0.8, 1.0), -0.067968, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.3, 0.6, 1.0), -0.279936, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.5, 0.4, 1.0), -0.32976,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.7, 0.2, 1.0), -0.158256, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.9, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.1, 0.6, 1.0),  0.186624, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.3, 0.4, 1.0),  0.028512, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.5, 0.2, 1.0), -0.07992,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.7, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.1, 0.4, 1.0),  0.36144,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.3, 0.2, 1.0),  0.17928,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.5, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.1, 0.2, 1.0),  0.204624, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.3, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.0, 0.1, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.1, 0.0, 1.0),  0.0,      eps);

                // three particle LCDA scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.1, 0.8, 2.0), -0.0160271, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.3, 0.6, 2.0), -0.0652682, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.5, 0.4, 2.0), -0.0766598, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.7, 0.2, 2.0), -0.0367404, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.1, 0.6, 2.0),  0.0366071, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.3, 0.4, 2.0), -0.0021538, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.5, 0.2, 2.0), -0.0221321, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.1, 0.4, 2.0),  0.0742666, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.3, 0.2, 2.0),  0.0344653, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.1, 0.2, 2.0),  0.0424959, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.1, 0.8, 2.0), -0.072636,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.3, 0.6, 2.0), -0.0383577, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.5, 0.4, 2.0),  0.0339688, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.7, 0.2, 2.0),  0.0337227, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.1, 0.6, 2.0), -0.0095838, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.3, 0.4, 2.0),  0.0995089, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.5, 0.2, 2.0),  0.0882485, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.1, 0.4, 2.0),  0.0765967, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.3, 0.2, 2.0),  0.104234,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.1, 0.2, 2.0),  0.0561024, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.1, 0.8, 2.0), -0.049339, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.3, 0.6, 2.0), -0.212137, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.5, 0.4, 2.0), -0.252603, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.7, 0.2, 2.0), -0.121824, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.9, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.1, 0.6, 2.0),  0.130711, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.3, 0.4, 2.0),  0.002445, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.5, 0.2, 2.0), -0.070580, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.7, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.1, 0.4, 2.0),  0.25532,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.3, 0.2, 2.0),  0.119891, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.5, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.1, 0.2, 2.0),  0.144836, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.3, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.0, 0.1, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.1, 0.0, 2.0),  0.0,      eps);
            }

            /* Twist 4 */
            {
                AntiKStarLCDAs kstar(p, Options{ });

                // parameters at mu = 1.0, 2.0, 3.0, 4.0, 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(1.0),        0.02,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(2.0),        0.0165725,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(3.0),        0.0153772,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(4.0),        0.0147015,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(5.0),        0.0142425,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(1.0), -0.02,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(2.0), -0.0117872,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(3.0), -0.00954933, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(4.0), -0.00841563, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(5.0), -0.00769734, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(1.0),       -0.01,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(2.0),       -0.00843717, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(3.0),       -0.00784189, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(4.0),       -0.00749527, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(5.0),       -0.00725593, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(1.0),  -0.05,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(2.0),  -0.0365548,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(3.0),  -0.0322774,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(4.0),  -0.0299564,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(5.0),  -0.0284201,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(1.0),      -0.0210942,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(2.0),      -0.017223,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(3.0),      -0.0158359,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(4.0),      -0.0150461,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(5.0),      -0.0145079,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(1.0),       0.0135855,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(2.0),       0.0128504,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(3.0),       0.0124729,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(4.0),       0.0122315,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(5.0),       0.0120558,  eps);
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
            // switch up and strange mass and flip sign of odd parameters
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
            p["K^*::zeta4para@1GeV"]        = 0.02;
            p["K^*::omega4paratilde@1GeV"]  = -0.02;
            p["K^*::zeta4perp@1GeV"]        = -0.01;
            p["K^*::zeta4perptilde@1GeV"]   = -0.05;
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

                // phiperp LCDA at various u values for mu = 2.0
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

                // two particle LCDAs at scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.1, 1.0), 0.480777, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.3, 1.0), 0.931977, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.5, 1.0), 1.034104, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.7, 1.0), 1.072796, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.9, 1.0), 0.675659, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.1, 1.0), 1.35473,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.3, 1.0), 0.525137, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.5, 1.0), 0.795833, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.7, 1.0), 0.357182, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.9, 1.0), 1.53821,  eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.1, 1.0), 0.499909, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.3, 1.0), 1.04631,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.5, 1.0), 1.14539,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.7, 1.0), 1.14576,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.9, 1.0), 0.730355, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.1, 1.0), 1.10177,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.3, 1.0), 0.844701, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.5, 1.0), 0.883863, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.7, 1.0), 0.739837, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.9, 1.0), 1.32029,  eps);

                // two particle LCDAs at scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.1, 2.0), 0.488303, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.3, 2.0), 1.00108,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.5, 2.0), 1.13694,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.7, 2.0), 1.10367,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3para(0.9, 2.0), 0.629763, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.1, 2.0), 1.46711,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.3, 2.0), 0.509224, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.5, 2.0), 0.612806, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.7, 2.0), 0.422571, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3para(0.9, 2.0), 1.62548,  eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.1, 2.0), 0.507961, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.3, 2.0), 1.10493,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.5, 2.0), 1.25996,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.7, 2.0), 1.18558,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.psi3perp(0.9, 2.0), 0.662752, eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.1, 2.0), 1.13684,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.3, 2.0), 0.84257,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.5, 2.0), 0.835106, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.7, 2.0), 0.79929,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.phi3perp(0.9, 2.0), 1.29942,  eps);

                // three particle LCDA scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.1, 0.8, 1.0), -0.0239616, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.3, 0.6, 1.0), -0.0964224, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.5, 0.4, 1.0), -0.112896,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.7, 0.2, 1.0), -0.0540288, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.9, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.1, 0.6, 1.0),  0.0590976, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.3, 0.4, 1.0),  0.0041472, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.5, 0.2, 1.0), -0.029376,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.7, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.1, 0.4, 1.0),  0.117504,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.3, 0.2, 1.0),  0.057024,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.5, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.1, 0.2, 1.0),  0.0669312, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.3, 0.0, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.0, 0.1, 1.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.1, 0.0, 1.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.1, 0.8, 1.0), -0.156672, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.3, 0.6, 1.0), -0.101088, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.5, 0.4, 1.0),  0.04608,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.7, 0.2, 1.0),  0.058464, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.9, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.1, 0.6, 1.0), -0.046656, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.3, 0.4, 1.0),  0.15552,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.5, 0.2, 1.0),  0.15552,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.7, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.1, 0.4, 1.0),  0.12672,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.3, 0.2, 1.0),  0.18576,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.5, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.1, 0.2, 1.0),  0.1008,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.3, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.0, 0.1, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.1, 0.0, 1.0),  0.0,      eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.1, 0.8, 1.0), -0.067968, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.3, 0.6, 1.0), -0.279936, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.5, 0.4, 1.0), -0.32976,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.7, 0.2, 1.0), -0.158256, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.9, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.1, 0.6, 1.0),  0.186624, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.3, 0.4, 1.0),  0.028512, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.5, 0.2, 1.0), -0.07992,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.7, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.1, 0.4, 1.0),  0.36144,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.3, 0.2, 1.0),  0.17928,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.5, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.1, 0.2, 1.0),  0.204624, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.3, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.0, 0.1, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.1, 0.0, 1.0),  0.0,      eps);

                // three particle LCDA scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.1, 0.8, 2.0), -0.0160271, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.3, 0.6, 2.0), -0.0652682, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.5, 0.4, 2.0), -0.0766598, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.7, 0.2, 2.0), -0.0367404, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.1, 0.6, 2.0),  0.0366071, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.3, 0.4, 2.0), -0.0021538, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.5, 0.2, 2.0), -0.0221321, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.1, 0.4, 2.0),  0.0742666, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.3, 0.2, 2.0),  0.0344653, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.1, 0.2, 2.0),  0.0424959, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3para(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.1, 0.8, 2.0), -0.072636,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.3, 0.6, 2.0), -0.0383577, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.5, 0.4, 2.0),  0.0339688, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.7, 0.2, 2.0),  0.0337227, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.1, 0.6, 2.0), -0.0095838, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.3, 0.4, 2.0),  0.0995089, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.5, 0.2, 2.0),  0.0882485, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.1, 0.4, 2.0),  0.0765967, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.3, 0.2, 2.0),  0.104234,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.1, 0.2, 2.0),  0.0561024, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3paratilde(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.1, 0.8, 2.0), -0.049339, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.3, 0.6, 2.0), -0.212137, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.5, 0.4, 2.0), -0.252603, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.7, 0.2, 2.0), -0.121824, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.1, 0.9, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.1, 0.6, 2.0),  0.130711, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.3, 0.4, 2.0),  0.002445, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.5, 0.2, 2.0), -0.070580, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.3, 0.7, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.1, 0.4, 2.0),  0.25532,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.3, 0.2, 2.0),  0.119891, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.5, 0.5, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.1, 0.2, 2.0),  0.144836, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.7, 0.3, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.0, 0.1, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.Phi3perp(0.9, 0.1, 0.0, 2.0),  0.0,      eps);
            }

            /* Twist 4 */
            {
                KStarLCDAs kstar(p, Options{ });

                // parameters at mu = 1.0, 2.0, 3.0, 4.0, 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(1.0),        0.02,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(2.0),        0.0165725,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(3.0),        0.0153772,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(4.0),        0.0147015,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4para(5.0),        0.0142425,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(1.0), -0.02,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(2.0), -0.0117872,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(3.0), -0.00954933, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(4.0), -0.00841563, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.omega4paratilde(5.0), -0.00769734, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(1.0),       -0.01,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(2.0),       -0.00843717, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(3.0),       -0.00784189, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(4.0),       -0.00749527, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perp(5.0),       -0.00725593, eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(1.0),  -0.05,       eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(2.0),  -0.0365548,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(3.0),  -0.0322774,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(4.0),  -0.0299564,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.zeta4perptilde(5.0),  -0.0284201,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(1.0),      -0.0210942,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(2.0),      -0.017223,   eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(3.0),      -0.0158359,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(4.0),      -0.0150461,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4para(5.0),      -0.0145079,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(1.0),       0.0135855,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(2.0),       0.0128504,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(3.0),       0.0124729,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(4.0),       0.0122315,  eps);
                TEST_CHECK_NEARLY_EQUAL(kstar.kappa4perp(5.0),       0.0120558,  eps);
            }
        }
} kstar_lcdas_test;
