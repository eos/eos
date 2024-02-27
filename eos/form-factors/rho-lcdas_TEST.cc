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
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.3, 1.0), 0.911333, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.5, 1.0), 1.455,    eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.7, 1.0), 0.911333, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(1.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.0, 1.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.3, 1.0), 0.792225, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.5, 1.0), 1.88813,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.7, 1.0), 0.792225, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(1.0, 1.0), 0.0,      eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.3, 2.0), 1.02489,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.5, 2.0), 1.4244,   eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(0.7, 2.0), 1.02489,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phipara(1.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.0, 2.0), 0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.3, 2.0), 0.951784, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.5, 2.0), 1.72422,  eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(0.7, 2.0), 0.951784, eps);
                TEST_CHECK_NEARLY_EQUAL(rho.phiperp(1.0, 2.0), 0.0,      eps);
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
            }
        }
} rho_lcdas_test;
