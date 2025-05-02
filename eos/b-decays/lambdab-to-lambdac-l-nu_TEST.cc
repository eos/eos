/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019      Ahmet Kokulu
 * Copyright (c) 2019-2025 Danny van Dyk
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
#include <eos/maths/complex.hh>
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
            // tests for SM observables, Re{cVL}=1.0 in the SM and all other couplings are zero, l = mu
            {
                Parameters p = Parameters::Defaults();
                p["Lambda_c::alpha"]       = -0.78;
                p["mass::Lambda_b"]        = 5.6194;
                p["mass::Lambda_c"]        = 2.2865;
                p["CKM::abs(V_cb)"]        =  0.041996951916414726;

                // the parameters are fixed as EOS default values

                Options oo
                {
                    { "model"_ok,        "WET" },
                    { "form-factors"_ok, "DKMR2017"   },
                    { "l"_ok,            "mu"         }
                };

                LambdaBToLambdaCLeptonNeutrino d(p, oo);

                const double eps = 1e-4;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011, 11.1), -0.20167, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_hadronic(0.011, 11.1),  0.32745, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_combined(0.011, 11.1), -0.11727, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_fzero(0.011, 11.1),          0.58742, eps);
            }

            // tests for SM observables, Re{cVL}=1.0 in the SM and all other couplings are zero, l = mu
            {
                Parameters p = Parameters::Defaults();
                p["Lambda_c::alpha"]       = -0.78;
                p["mass::Lambda_b"]        = 5.6194;
                p["mass::Lambda_c"]        = 2.2865;
                p["CKM::abs(V_cb)"]        =  0.041996951916414726;

                // the parameters are fixed as EOS default values

                Options oo
                {
                    { "model"_ok,        "WET" },
                    { "form-factors"_ok, "DKMR2017"   },
                    { "l"_ok,            "tau"        }
                };

                LambdaBToLambdaCLeptonNeutrino d(p, oo);

                const double eps = 1e-4;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(3.154, 11.1), +0.02447,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_hadronic(3.154, 11.1),  0.29600,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_combined(3.154, 11.1), -0.022086, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_fzero(3.154, 11.1),          0.38041,  eps);
            }

            // Consistency check for R_lambda
            {
                Parameters p = Parameters::Defaults();
                p["Lambda_c::alpha"]       = -0.78;
                p["mass::Lambda_b"]        = 5.6194;
                p["mass::Lambda_c"]        = 2.2865;
                p["CKM::abs(V_cb)"]        =  0.041996951916414726;

                Options oo
                {
                    { "model"_ok,        "WET" },
                    { "form-factors"_ok, "DKMR2017"   },
                    { "l"_ok,            "mu"        }
                };
                LambdaBToLambdaCLeptonNeutrino dmu(p, oo);

                oo.declare("l"_ok, "tau");
                LambdaBToLambdaCLeptonNeutrino dtau(p, oo);

                oo =
                {
                    { "model"_ok,        "WET" },
                    { "form-factors"_ok, "DKMR2017"   },
                };
                Kinematics k
                {
                    { "q2_mu_min",   0.011 }, { "q2_mu_max",  11.1 },
                    { "q2_tau_min",  3.154 }, { "q2_tau_max", 11.1 },
                };

                auto obs_Rlambda = Observable::make("Lambda_b->Lambda_clnu::R(Lambda_c)", p, k, oo);
                TEST_CHECK_RELATIVE_ERROR(
                    dtau.integrated_branching_ratio(3.154, 11.1) / dmu.integrated_branching_ratio(0.011, 11.1),
                    obs_Rlambda->evaluate(),
                    1e-5
                );

            }

            // tests for NP observables (no tensors)
            {
                Parameters p = Parameters::Defaults();
                // the rest of input is fixed to default EOS values
                p["CKM::abs(V_cb)"]      =  0.041996951916414726;
                p["cbmunumu::Re{cVL}"]   =  1.0;
                p["cbmunumu::Im{cVL}"]   = -1.0;
                p["cbmunumu::Re{cVR}"]   =  2.0;
                p["cbmunumu::Im{cVR}"]   = -2.0;
                p["cbmunumu::Re{cSL}"]   =  3.0;
                p["cbmunumu::Im{cSL}"]   = -3.0;
                p["cbmunumu::Re{cSR}"]   =  4.0;
                p["cbmunumu::Im{cSR}"]   = -4.0;
                p["cbmunumu::Re{cT}"]    =  0.0;
                p["cbmunumu::Im{cT}"]    =  0.0;
                // fix the scale
                p["cbmunumu::mu"]          =  4.18;
                p["mass::b(MSbar)"]        =  4.18;
                p["mass::c"]               =  1.275;
                p["Lambda_c::alpha"]       = -0.78;
                p["mass::Lambda_b"]        = 5.6194;
                p["mass::Lambda_c"]        = 2.2865;

                Options oo;
                oo.declare("model"_ok, "WET");
                oo.declare("form-factors"_ok, "DKMR2017");
                oo.declare("l"_ok, "mu");

                LambdaBToLambdaCLeptonNeutrino d(p, oo);

                const double eps = 1e-4;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011, 11.1),   0.04665,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_hadronic(0.011, 11.1),  -0.01808,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_combined(0.011, 11.1),  -0.015045, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_fzero(0.011, 11.1),           0.401858, eps);
            }

            // tests for NP observables (no tensors)
            {
                Parameters p = Parameters::Defaults();
                // the rest of input is fixed to default EOS values
                p["CKM::abs(V_cb)"]      =  0.041996951916414726;
                p["cbmunumu::Re{cVL}"]   =  1.0;
                p["cbmunumu::Im{cVL}"]   = -1.0;
                p["cbmunumu::Re{cVR}"]   =  2.0;
                p["cbmunumu::Im{cVR}"]   = -2.0;
                p["cbmunumu::Re{cSL}"]   =  3.0;
                p["cbmunumu::Im{cSL}"]   = -3.0;
                p["cbmunumu::Re{cSR}"]   =  4.0;
                p["cbmunumu::Im{cSR}"]   = -4.0;
                p["cbmunumu::Re{cT}"]    =  1.0;
                p["cbmunumu::Im{cT}"]    = -2.0;
                // fix the scale
                p["cbmunumu::mu"]          =  4.18;
                p["mass::b(MSbar)"]        =  4.18;
                p["mass::c"]               =  1.275;
                p["Lambda_c::alpha"]       = -0.78;
                p["mass::Lambda_b"]        = 5.6194;
                p["mass::Lambda_c"]        = 2.2865;

                Options oo;
                oo.declare("model"_ok, "WET");
                oo.declare("form-factors"_ok, "DKMR2017");
                oo.declare("l"_ok, "mu");

                LambdaBToLambdaCLeptonNeutrino d(p, oo);

                const double eps = 1e-2;

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_leptonic(0.011, 11.1),   0.1336, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_hadronic(0.011, 11.1),  -0.0147, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_a_fb_combined(0.011, 11.1),  -0.1180, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_fzero(0.011, 11.1),           0.3742, eps);
            }
        }
} lambdab_to_lambdac_l_nu_test;
