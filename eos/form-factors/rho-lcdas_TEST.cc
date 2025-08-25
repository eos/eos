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
#include <eos/form-factors/rho-lcdas.hh>

#include <eos/models/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class RhoLCDAsTest :
    public TestCase
{
    public:
        RhoLCDAsTest() :
            TestCase("rho_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"]          = 0.1176;
            p["mass::d(2GeV)"]             = 0.0048;
            p["mass::u(2GeV)"]             = 0.0032;
            p["rho::a2para@1GeV"]          = 0.22;
            p["rho::a4para@1GeV"]          = 0.16;
            p["rho::a2perp@1GeV"]          = 0.14;
            p["rho::a4perp@1GeV"]          = 0.25;
            p["rho::fperp@1GeV"]           = 0.16;
            p["rho::zeta3para@1GeV"]       = 0.03;
            p["rho::omega3paratilde@1GeV"] = -0.09;
            p["rho::omega3para@1GeV"]      = 0.15;
            p["rho::omega3perp@1GeV"]      = 0.55;
            p["rho::zeta4para@1GeV"]       = 0.07;
            p["rho::omega4paratilde@1GeV"] = -0.03;
            p["rho::zeta4perp@1GeV"]       = -0.03;
            p["rho::zeta4perptilde@1GeV"]  = -0.08;

            /* Diagnostics */
            {
                RhoLCDAs rho(p, Options{ });
                Diagnostics diagnostics = rho.diagnostics();
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
                RhoLCDAs rho(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV, 3.0 GeV, 4.0 GeV and 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.a1para(1.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1para(2.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1para(3.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1para(4.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1para(5.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2para(1.0),  0.22,        eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2para(2.0),  0.164004,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2para(3.0),  0.145901,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2para(4.0),  0.136008,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2para(5.0),  0.129431,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3para(1.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3para(2.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3para(3.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3para(4.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3para(5.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4para(1.0),  0.16,        eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4para(2.0),  0.1043234,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4para(3.0),  0.0879874,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4para(4.0),  0.0794371,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4para(5.0),  0.0739064,   eps);

                TEST_CHECK_NEARLY_EQUAL(rho.a1perp(1.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1perp(2.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1perp(3.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1perp(4.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a1perp(5.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2perp(1.0),  0.14,        eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2perp(2.0),  0.103147,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2perp(3.0),  0.0913332,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2perp(4.0),  0.0849017,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a2perp(5.0),  0.0806361,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3perp(1.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3perp(2.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3perp(3.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3perp(4.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a3perp(5.0),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4perp(1.0),  0.25,        eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4perp(2.0),  0.16224,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4perp(3.0),  0.13658,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4perp(4.0),  0.12317,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.a4perp(5.0),  0.11450,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.fperp(1.0),   0.16,        eps);
                TEST_CHECK_NEARLY_EQUAL(rho.fperp(2.0),   0.149109,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.fperp(3.0),   0.144981,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.fperp(4.0),   0.142559,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.fperp(5.0),   0.140873,    eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.3, 1.0), 0.911333, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.5, 1.0), 1.455,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.7, 1.0), 0.911333, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(1.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.3, 1.0), 0.792225, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.5, 1.0), 1.88813,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.7, 1.0), 0.792225, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(1.0, 1.0), 0.0,      eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.3, 2.0), 1.02489,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.5, 2.0), 1.4244,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(0.7, 2.0), 1.02489,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2para(1.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.3, 2.0), 0.951784, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.5, 2.0), 1.72422,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(0.7, 2.0), 0.951784, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi2perp(1.0, 2.0), 0.0,      eps);
            }

            /* Twist 3 */
            {
                RhoLCDAs rho(p, Options{ });

                // parameters at mu = 1.0 GeV, 2.0 GeV, 3.0 GeV, 4.0 GeV and 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.zeta3para(1.0),        0.03     , eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta3para(2.0),        0.0190839, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta3para(3.0),        0.0159382, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta3para(4.0),        0.0143048, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta3para(5.0),        0.0132535, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3paratilde(1.0), -0.09,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3paratilde(2.0), -0.0471515, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3paratilde(3.0), -0.0365076, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3paratilde(4.0), -0.0313243, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3paratilde(5.0), -0.0281206, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3para(1.0),       0.15     , eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3para(2.0),       0.0932984, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3para(3.0),       0.0770507, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3para(4.0),       0.0686532, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3para(5.0),       0.063268,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3perp(1.0),       0.55,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3perp(2.0),       0.384353,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3perp(3.0),       0.33324,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3perp(4.0),       0.305884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega3perp(5.0),       0.287936,  eps);

                // two particle LCDAs at scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.1, 1.0), 0.85383, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.3, 1.0), 1.19343, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.5, 1.0), 1.10375, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.7, 1.0), 1.19343, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.9, 1.0), 0.85383, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.1, 1.0), 1.30613, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.3, 1.0), 0.09533, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.5, 1.0), 1.03125, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.7, 1.0), 0.09533, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.9, 1.0), 1.30613, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.1, 1.0), 0.80953, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.3, 1.0), 1.20283, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.5, 1.0), 1.15969, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.7, 1.0), 1.20283, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.9, 1.0), 0.80953, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.1, 1.0), 1.25141, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.3, 1.0), 0.64436, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.5, 1.0), 0.83297, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.7, 1.0), 0.64436, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.9, 1.0), 1.25141, eps);

                // two particle LCDAs at scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.1, 2.0), 0.76089, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.3, 2.0), 1.21314, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.5, 2.0), 1.2211,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.7, 2.0), 1.21314, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3para(0.9, 2.0), 0.76089, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.1, 2.0), 1.49203, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.3, 2.0), 0.20838, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.5, 2.0), 0.72066, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.7, 2.0), 0.20838, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3para(0.9, 2.0), 1.49203, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.1, 2.0), 0.71328, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.3, 2.0), 1.22324, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.5, 2.0), 1.28122, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.7, 2.0), 1.22324, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi3perp(0.9, 2.0), 0.71328, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.1, 2.0), 1.25577, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.3, 2.0), 0.72239, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.5, 2.0), 0.78818, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.7, 2.0), 0.72239, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi3perp(0.9, 2.0), 1.25577, eps);

                // three particle LCDA scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.3, 0.6, 1.0), -0.11664, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.5, 0.4, 1.0), -0.1728,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.7, 0.2, 1.0), -0.09072, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.9, 0.0, 1.0),  0.0    , eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.1, 0.6, 1.0),  0.11664, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.5, 0.2, 1.0), -0.0648,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.5, 0.1, 0.4, 1.0),  0.1728,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.5, 0.3, 0.2, 1.0),  0.0648,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.7, 0.1, 0.2, 1.0),  0.09072, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.1, 0.8, 1.0), -0.200448, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.3, 0.6, 1.0), -0.093312, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.5, 0.4, 1.0),  0.11232,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.7, 0.2, 1.0),  0.102816, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.9, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.1, 0.6, 1.0), -0.093312, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.3, 0.4, 1.0),  0.202176, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.5, 0.2, 1.0),  0.22032,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.7, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.5, 0.1, 0.4, 1.0),  0.11232,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.5, 0.3, 0.2, 1.0),  0.22032,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.5, 0.5, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.7, 0.1, 0.2, 1.0),  0.102816, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.7, 0.3, 0.0, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.9, 0.0, 0.1, 1.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.9, 0.1, 0.0, 1.0),  0.0,      eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.3, 0.6, 1.0), -0.42768, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.5, 0.4, 1.0), -0.6336,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.7, 0.2, 1.0), -0.33264, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.1, 0.6, 1.0),  0.42768, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.5, 0.2, 1.0), -0.2376,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.5, 0.1, 0.4, 1.0),  0.6336,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.5, 0.3, 0.2, 1.0),  0.2376,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.7, 0.1, 0.2, 1.0),  0.33264, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                // three particle LCDA scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.1, 0.8, 2.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.3, 0.6, 2.0), -0.0725489, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.5, 0.4, 2.0), -0.10748,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.7, 0.2, 2.0), -0.0564269, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.1, 0.9, 0.0, 2.0),  0.0    , eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.1, 0.6, 2.0),  0.0725489, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.3, 0.4, 2.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.5, 0.2, 2.0), -0.0403049,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.3, 0.7, 0.0, 2.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.5, 0.1, 0.4, 2.0),  0.10748,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.5, 0.3, 0.2, 2.0),  0.0403049,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.5, 0.5, 0.0, 2.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.7, 0.1, 0.2, 2.0),  0.0564269, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.7, 0.3, 0.0, 2.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.9, 0.0, 0.1, 2.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3para(0.9, 0.1, 0.0, 2.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.1, 0.8, 2.0), -0.0972587, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.3, 0.6, 2.0), -0.0357967, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.5, 0.4, 2.0),  0.0685413, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.7, 0.2, 2.0),  0.0572595, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.1, 0.6, 2.0), -0.0357967, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.3, 0.4, 2.0),  0.123374,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.5, 0.2, 2.0),  0.122699,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.5, 0.1, 0.4, 2.0),  0.0685413, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.5, 0.3, 0.2, 2.0),  0.122699,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.7, 0.1, 0.2, 2.0),  0.0572595, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3paratilde(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.1, 0.8, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.3, 0.6, 2.0), -0.298873, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.5, 0.4, 2.0), -0.442775, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.7, 0.2, 2.0), -0.232457, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.1, 0.9, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.1, 0.6, 2.0),  0.298873, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.3, 0.4, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.5, 0.2, 2.0), -0.166041, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.3, 0.7, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.5, 0.1, 0.4, 2.0),  0.442775, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.5, 0.3, 0.2, 2.0),  0.166041, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.5, 0.5, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.7, 0.1, 0.2, 2.0),  0.232457, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.7, 0.3, 0.0, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.9, 0.0, 0.1, 2.0),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi3perp(0.9, 0.1, 0.0, 2.0),  0.0,      eps);
            }

                 /* Twist 4 */
            {
                RhoLCDAs rho(p, Options{ });

                // parameters at mu = 1.0 GeV, 2.0 GeV, 3.0 GeV, 4.0 GeV and 5.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4para(1.0),        0.07,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4para(2.0),        0.0580036, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4para(3.0),        0.0538201, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4para(4.0),        0.0514552, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4para(5.0),        0.0498486, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega4paratilde(1.0), -0.03,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega4paratilde(2.0), -0.0176807, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega4paratilde(3.0), -0.014324,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega4paratilde(4.0), -0.0126234, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.omega4paratilde(5.0), -0.011546,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perp(1.0),       -0.03,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perp(2.0),       -0.0236691, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perp(3.0),       -0.0215038, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perp(4.0),       -0.0202925, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perp(5.0),       -0.0194754, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perptilde(1.0),  -0.08,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perptilde(2.0),  -0.0588162, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perptilde(3.0),  -0.0520483, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perptilde(4.0),  -0.0483689, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.zeta4perptilde(5.0),  -0.0459306, eps);

                // Two-particle LCDAs at mu = 1 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.1, 1.0), 0.56698, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.3, 1.0), 1.33566, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.5, 1.0), 1.42544, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.7, 1.0), 1.33566, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.9, 1.0), 0.56698, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.1, 1.0), 0.39426, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.3, 1.0), 1.33561, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.5, 1.0), 1.52617, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.7, 1.0), 1.33561, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.9, 1.0), 0.39426, eps);

                // Two-particle LCDAs at mu = 2 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.1, 2.0), 0.51937, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.3, 2.0), 1.33651, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.5, 2.0), 1.49017, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.7, 2.0), 1.33651, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.psi4para(0.9, 2.0), 0.51937, eps);

                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.1, 2.0), 0.32926, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.3, 2.0), 1.27599, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.5, 2.0), 1.57309, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.7, 2.0), 1.27599, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phi4para(0.9, 2.0), 0.32926, eps);

                // Three-particle LCDAs at mu = 1 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.3, 0.6, 1.0), -0.01764, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.5, 0.4, 1.0), -0.03920, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.7, 0.2, 1.0), -0.04116, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.1, 0.6, 1.0),  0.01764, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.5, 0.2, 1.0), -0.02940, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.5, 0.1, 0.4, 1.0),  0.03920, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.5, 0.3, 0.2, 1.0),  0.02940, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.7, 0.1, 0.2, 1.0),  0.04116, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.1, 0.8, 1.0), -0.03248, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.3, 0.6, 1.0), -0.02016, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.5, 0.4, 1.0),  0.03640, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.7, 0.2, 1.0),  0.06664, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.1, 0.6, 1.0), -0.02016, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.3, 0.4, 1.0),  0.06552, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.5, 0.2, 1.0),  0.1428,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.5, 0.1, 0.4, 1.0),  0.0364,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.5, 0.3, 0.2, 1.0),  0.1428,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.7, 0.1, 0.2, 1.0),  0.06664, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.3, 0.6, 1.0),  0.01512, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.5, 0.4, 1.0),  0.07616, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.7, 0.2, 1.0),  0.05208, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.1, 0.6, 1.0), -0.01512, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.5, 0.2, 1.0),  0.01736, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.5, 0.1, 0.4, 1.0), -0.07616, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.5, 0.3, 0.2, 1.0), -0.01736, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.7, 0.1, 0.2, 1.0), -0.05208, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.9, 0.0, 0.1, 1.0), -0.02394, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.1, 0.8, 1.0), -0.09856, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.3, 0.6, 1.0),  0.03024, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.5, 0.4, 1.0),  0.11424, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.7, 0.2, 1.0),  0.06944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.1, 0.6, 1.0),  0.03024, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.3, 0.4, 1.0),  0.11424, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.5, 0.2, 1.0),  0.06944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.5, 0.1, 0.4, 1.0),  0.11424, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.5, 0.3, 0.2, 1.0),  0.06944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.7, 0.1, 0.2, 1.0),  0.06944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.9, 0.0, 0.1, 1.0),  0.02394, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.1, 0.8, 1.0),  0.12672, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.3, 0.6, 1.0), -0.03888, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.5, 0.4, 1.0), -0.14688, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.7, 0.2, 1.0), -0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.1, 0.6, 1.0), -0.03888, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.3, 0.4, 1.0), -0.14688, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.5, 0.2, 1.0), -0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.5, 0.1, 0.4, 1.0), -0.14688, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.5, 0.3, 0.2, 1.0), -0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.7, 0.1, 0.2, 1.0), -0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.9, 0.0, 0.1, 1.0), -0.03078, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.1, 0.8, 1.0), -0.12672, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.3, 0.6, 1.0),  0.03888, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.5, 0.4, 1.0),  0.14688, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.7, 0.2, 1.0),  0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.1, 0.6, 1.0),  0.03888, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.3, 0.4, 1.0),  0.14688, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.5, 0.2, 1.0),  0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.5, 0.1, 0.4, 1.0),  0.14688, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.5, 0.3, 0.2, 1.0),  0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.7, 0.1, 0.2, 1.0),  0.08928, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.9, 0.0, 0.1, 1.0),  0.03078, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.3, 0.6, 1.0),  0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.5, 0.4, 1.0),  0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.7, 0.2, 1.0),  0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.1, 0.6, 1.0), -0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.5, 0.2, 1.0),  0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.5, 0.1, 0.4, 1.0), -0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.5, 0.3, 0.2, 1.0), -0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.7, 0.1, 0.2, 1.0), -0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.3, 0.6, 1.0), -0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.5, 0.4, 1.0), -0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.7, 0.2, 1.0), -0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.1, 0.6, 1.0),  0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.5, 0.2, 1.0), -0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.5, 0.1, 0.4, 1.0),  0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.5, 0.3, 0.2, 1.0),  0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.7, 0.1, 0.2, 1.0),  0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.9, 0.0, 0.1, 1.0),  0.03078, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.3, 0.6, 1.0), -0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.5, 0.4, 1.0), -0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.7, 0.2, 1.0), -0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.1, 0.6, 1.0),  0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.5, 0.2, 1.0), -0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.5, 0.1, 0.4, 1.0),  0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.5, 0.3, 0.2, 1.0),  0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.7, 0.1, 0.2, 1.0),  0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.3, 0.6, 1.0), -0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.5, 0.4, 1.0), -0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.7, 0.2, 1.0), -0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.1, 0.6, 1.0),  0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.5, 0.2, 1.0), -0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.5, 0.1, 0.4, 1.0),  0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.5, 0.3, 0.2, 1.0),  0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.7, 0.1, 0.2, 1.0),  0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.9, 0.0, 0.1, 1.0),  0.03078, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.3, 0.6, 1.0),  0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.5, 0.4, 1.0),  0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.7, 0.2, 1.0),  0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.1, 0.6, 1.0), -0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.5, 0.2, 1.0),  0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.5, 0.1, 0.4, 1.0), -0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.5, 0.3, 0.2, 1.0), -0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.7, 0.1, 0.2, 1.0), -0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.3, 0.6, 1.0),  0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.5, 0.4, 1.0),  0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.7, 0.2, 1.0),  0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.1, 0.6, 1.0), -0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.5, 0.2, 1.0),  0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.5, 0.1, 0.4, 1.0), -0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.5, 0.3, 0.2, 1.0), -0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.7, 0.1, 0.2, 1.0), -0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.9, 0.0, 0.1, 1.0), -0.03078, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.3, 0.6, 1.0), -0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.5, 0.4, 1.0), -0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.7, 0.2, 1.0), -0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.1, 0.6, 1.0),  0.02268, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.5, 0.2, 1.0), -0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.5, 0.1, 0.4, 1.0),  0.0504,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.5, 0.3, 0.2, 1.0),  0.0378,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.7, 0.1, 0.2, 1.0),  0.05292, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.9, 0.0, 0.1, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.1, 0.8, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.3, 0.6, 1.0),  0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.5, 0.4, 1.0),  0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.7, 0.2, 1.0),  0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.9, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.1, 0.6, 1.0), -0.01944, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.3, 0.4, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.5, 0.2, 1.0),  0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.7, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.5, 0.1, 0.4, 1.0), -0.09792, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.5, 0.3, 0.2, 1.0), -0.02232, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.5, 0.5, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.7, 0.1, 0.2, 1.0), -0.06696, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.7, 0.3, 0.0, 1.0),  0.0,     eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.9, 0.0, 0.1, 1.0), -0.03078, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.9, 0.1, 0.0, 1.0),  0.0,     eps);

                // Three-particle LCDAs at mu = 2 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.3, 0.6, 2.0), -0.0146169, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.5, 0.4, 2.0), -0.032482,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.7, 0.2, 2.0), -0.0341061, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.1, 0.6, 2.0),  0.0146169, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.5, 0.2, 2.0), -0.0243615, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.5, 0.1, 0.4, 2.0),  0.032482,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.5, 0.3, 0.2, 2.0),  0.0243615, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.7, 0.1, 0.2, 2.0),  0.0341061, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4para(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.1, 0.8, 2.0), -0.0269137, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.3, 0.6, 2.0), -0.016705,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.5, 0.4, 2.0),  0.0301619, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.7, 0.2, 2.0),  0.0552194, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.1, 0.6, 2.0), -0.016705,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.3, 0.4, 2.0),  0.0542914, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.5, 0.2, 2.0),  0.118327,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.5, 0.1, 0.4, 2.0),  0.0301619, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.5, 0.3, 0.2, 2.0),  0.118327,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.7, 0.1, 0.2, 2.0),  0.0552194, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4paratilde(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.3, 0.6, 2.0),  0.0125288, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.5, 0.4, 2.0),  0.0631079, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.7, 0.2, 2.0),  0.0431547, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.1, 0.6, 2.0), -0.0125288, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.5, 0.2, 2.0),  0.0143849, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.5, 0.1, 0.4, 2.0), -0.0631079, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.5, 0.3, 0.2, 2.0), -0.0143849, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.7, 0.1, 0.2, 2.0), -0.0431547, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.9, 0.0, 0.1, 2.0), -0.0198372, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4para(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.1, 0.8, 2.0), -0.0816691, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.3, 0.6, 2.0),  0.0250576, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.5, 0.4, 2.0),  0.0946619, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.7, 0.2, 2.0),  0.0575396, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.1, 0.6, 2.0),  0.0250576, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.3, 0.4, 2.0),  0.0946619, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.5, 0.2, 2.0),  0.0575396, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.5, 0.1, 0.4, 2.0),  0.0946619, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.5, 0.3, 0.2, 2.0),  0.0575396, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.7, 0.1, 0.2, 2.0),  0.0575396, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.9, 0.0, 0.1, 2.0),  0.0198372, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4paratilde(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.1, 0.8, 2.0),  0.0999784, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.3, 0.6, 2.0), -0.0306752, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.5, 0.4, 2.0), -0.115884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.7, 0.2, 2.0), -0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.1, 0.6, 2.0), -0.0306752, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.3, 0.4, 2.0), -0.115884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.5, 0.2, 2.0), -0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.5, 0.1, 0.4, 2.0), -0.115884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.5, 0.3, 0.2, 2.0), -0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.7, 0.1, 0.2, 2.0), -0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.9, 0.0, 0.1, 2.0), -0.0242845, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perp(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.1, 0.8, 2.0), -0.0999784, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.3, 0.6, 2.0),  0.0306752, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.5, 0.4, 2.0),  0.115884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.7, 0.2, 2.0),  0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.1, 0.6, 2.0),  0.0306752, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.3, 0.4, 2.0),  0.115884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.5, 0.2, 2.0),  0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.5, 0.1, 0.4, 2.0),  0.115884,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.5, 0.3, 0.2, 2.0),  0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.7, 0.1, 0.2, 2.0),  0.0704393, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.9, 0.0, 0.1, 2.0),  0.0242845, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Psi4perptilde(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.3, 0.6, 2.0),  0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.5, 0.4, 2.0),  0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.7, 0.2, 2.0),  0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.1, 0.6, 2.0), -0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.5, 0.2, 2.0),  0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.5, 0.1, 0.4, 2.0), -0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.5, 0.3, 0.2, 2.0), -0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.7, 0.1, 0.2, 2.0), -0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp1(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.3, 0.6, 2.0), -0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.5, 0.4, 2.0), -0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.7, 0.2, 2.0), -0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.1, 0.6, 2.0),  0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.5, 0.2, 2.0), -0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.5, 0.1, 0.4, 2.0),  0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.5, 0.3, 0.2, 2.0),  0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.7, 0.1, 0.2, 2.0),  0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.9, 0.0, 0.1, 2.0),  0.0242845, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp2(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.3, 0.6, 2.0), -0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.5, 0.4, 2.0), -0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.7, 0.2, 2.0), -0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.1, 0.6, 2.0),  0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.5, 0.2, 2.0), -0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.5, 0.1, 0.4, 2.0),  0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.5, 0.3, 0.2, 2.0),  0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.7, 0.1, 0.2, 2.0),  0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp3(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.3, 0.6, 2.0), -0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.5, 0.4, 2.0), -0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.7, 0.2, 2.0), -0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.1, 0.6, 2.0),  0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.5, 0.2, 2.0), -0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.5, 0.1, 0.4, 2.0),  0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.5, 0.3, 0.2, 2.0),  0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.7, 0.1, 0.2, 2.0),  0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.9, 0.0, 0.1, 2.0),  0.0242845, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perp4(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.3, 0.6, 2.0),  0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.5, 0.4, 2.0),  0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.7, 0.2, 2.0),  0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.1, 0.6, 2.0), -0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.5, 0.2, 2.0),  0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.5, 0.1, 0.4, 2.0), -0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.5, 0.3, 0.2, 2.0), -0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.7, 0.1, 0.2, 2.0), -0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde1(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.3, 0.6, 2.0),  0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.5, 0.4, 2.0),  0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.7, 0.2, 2.0),  0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.1, 0.6, 2.0), -0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.5, 0.2, 2.0),  0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.5, 0.1, 0.4, 2.0), -0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.5, 0.3, 0.2, 2.0), -0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.7, 0.1, 0.2, 2.0), -0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.9, 0.0, 0.1, 2.0), -0.0242845, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde2(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.3, 0.6, 2.0), -0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.5, 0.4, 2.0), -0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.7, 0.2, 2.0), -0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.1, 0.6, 2.0),  0.0178939, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.5, 0.2, 2.0), -0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.5, 0.1, 0.4, 2.0),  0.0397641, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.5, 0.3, 0.2, 2.0),  0.0298231, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.7, 0.1, 0.2, 2.0),  0.0417523, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.9, 0.0, 0.1, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde3(0.9, 0.1, 0.0, 2.0),  0.0,       eps);

                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.1, 0.8, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.3, 0.6, 2.0),  0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.5, 0.4, 2.0),  0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.7, 0.2, 2.0),  0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.1, 0.9, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.1, 0.6, 2.0), -0.0153376, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.3, 0.4, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.5, 0.2, 2.0),  0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.3, 0.7, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.5, 0.1, 0.4, 2.0), -0.077256,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.5, 0.3, 0.2, 2.0), -0.0176098, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.5, 0.5, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.7, 0.1, 0.2, 2.0), -0.0528295, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.7, 0.3, 0.0, 2.0),  0.0,       eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.9, 0.0, 0.1, 2.0), -0.0242845, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.Phi4perptilde4(0.9, 0.1, 0.0, 2.0),  0.0,       eps);
            }
        }
} rho_lcdas_test;
