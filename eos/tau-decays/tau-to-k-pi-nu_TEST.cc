/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025-2026 Danny van Dyk
 * Copyright (c) 2025      Matthew Kirk
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

#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/tau-decays/tau-to-k-pi-nu.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class TauToKPiNeutrinoTest : public TestCase
{
    public:
        TauToKPiNeutrinoTest() :
            TestCase("tau_to_k_pi_nu_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p                      = Parameters::Defaults();
                p["0->Kpi::t_0@KSvD2025"]         = -4.0;
                p["0->Kpi::M_(+,0)@KSvD2025"]     = 0.892;
                p["0->Kpi::Gamma_(+,0)@KSvD2025"] = 0.046;
                p["0->Kpi::M_(0,0)@KSvD2025"]     = 0.680;
                p["0->Kpi::Gamma_(0,0)@KSvD2025"] = 0.600;
                p["0->Kpi::M_(+,1)@KSvD2025"]     = 1.368;
                p["0->Kpi::Gamma_(+,1)@KSvD2025"] = 0.212;
                p["0->Kpi::b_+^1@KSvD2025"]       = 0.0088;
                p["0->Kpi::b_+^2@KSvD2025"]       = -0.0052;
                p["0->Kpi::b_+^3@KSvD2025"]       = -0.0273;
                p["0->Kpi::b_+^4@KSvD2025"]       = -0.0094;
                p["0->Kpi::b_0^1@KSvD2025"]       = 0.2653;
                p["0->Kpi::b_0^2@KSvD2025"]       = -0.0588;
                p["0->Kpi::b_0^3@KSvD2025"]       = -0.4228;
                p["0->Kpi::b_0^4@KSvD2025"]       = -0.1646;

                const double eps = 1e-8;

                // tau -> K- pi0 nu
                TauToKPiNeutrino d(p,
                                   Options{
                                       { "n-resonances-1m"_ok, "2"_ov },
                                       { "n-resonances-0m"_ok, "1"_ov }
                });

                TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.0), 0.0017916075697680681, eps);
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 0.004430198678120239, eps); // Compare with PDG 0.433%

                TEST_CHECK_NEARLY_EQUAL(d.differential_pdf_q2(1.0), 0.0017916075697680681 / 0.004430198678120239, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_pdf_q2(0.4060, 3.1571), 0.3631445295018548, eps); // Compare with 1 / (3.1571 - 0.4060)

                // tau -> K_S pi- nu
                d = TauToKPiNeutrino(p,
                                     Options{
                                         {               "K"_ok, "K_S"_ov },
                                         { "n-resonances-1m"_ok,   "2"_ov },
                                         { "n-resonances-0m"_ok,   "1"_ov }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 4.25288e-3, eps); // Compare with [Belle:2014B] 4.16e-3

                // tau -> K_L pi- nu
                d = TauToKPiNeutrino(p,
                                     Options{
                                         {               "K"_ok, "K_L"_ov },
                                         { "n-resonances-1m"_ok,   "2"_ov },
                                         { "n-resonances-0m"_ok,   "1"_ov }
                });
                TEST_CHECK_NEARLY_EQUAL(d.total_branching_ratio(), 4.25288e-3, eps);

                // Regression test for the sign of the pseudoscalar amplitude, cf. gh#1203.
                // A non-zero scalar NP coefficient must produce a *linear* shift in the rate; with the
                // former sign bug the |S|^2 + |P|^2 combination cancelled the linear term, so that
                // the rate was identical under Re{cSL} -> -Re{cSL}. Requires the WET model, as the SM
                // model ignores the ustaunutau Wilson coefficients.
                const Options wet_opts{
                    {           "model"_ok, "WET"_ov },
                    { "n-resonances-1m"_ok,   "2"_ov },
                    { "n-resonances-0m"_ok,   "1"_ov }
                };
                p["ustaunutau::Re{cSL}"] = +0.1;
                TauToKPiNeutrino d_cSL_plus(p, wet_opts);
                TEST_CHECK_NEARLY_EQUAL(d_cSL_plus.total_branching_ratio(), 0.0046365212390183843, eps);
                p["ustaunutau::Re{cSL}"] = -0.1;
                TauToKPiNeutrino d_cSL_minus(p, wet_opts);
                TEST_CHECK_NEARLY_EQUAL(d_cSL_minus.total_branching_ratio(), 0.0042855354180517651, eps);
                p["ustaunutau::Re{cSL}"] = 0.0;

                // Two-fold differential rate d^2Gamma/dq2/dcos(theta_K) (here for tau -> K_L pi- nu,
                // the last instance constructed above). The distribution is a quadratic in
                // z = cos(theta_K); Simpson's rule is exact for a quadratic, so
                //     (f(-1) + 4 f(0) + f(+1)) / 3 == dBR/dq2 .
                const double f_m1 = d.double_differential_branching_ratio(1.0, -1.0);
                const double f_0  = d.double_differential_branching_ratio(1.0, 0.0);
                const double f_p1 = d.double_differential_branching_ratio(1.0, +1.0);
                TEST_CHECK_NEARLY_EQUAL((f_m1 + 4.0 * f_0 + f_p1) / 3.0, d.differential_branching_ratio(1.0), eps);
                TEST_CHECK_NEARLY_EQUAL(f_m1, 0.0013654922519557517, eps);
                TEST_CHECK_NEARLY_EQUAL(f_0, 0.00051133788806815986, eps);
                TEST_CHECK_NEARLY_EQUAL(f_p1, 0.0018269935336002955, eps);
                // the S-P interference (linear in cos(theta_K)) makes the distribution forward-backward
                // asymmetric: f(+1) != f(-1).
                TEST_CHECK(std::abs(f_p1 - f_m1) > 1.0e-2 * f_0);
            }

            // Single- and double-differential rates (channel tau -> K_S pi- nu) evaluated at the
            // mode of the FF_2_2_3_3 posterior of PR #1204 (prs/1204/{mode,analysis}.yaml), for two
            // Wilson-coefficient sets and cos(theta_K) in {-0.5, +0.1, +0.8}. The reference values
            // were produced from the f_+ and f_0 form factors of the mode and fed into an independent
            // Mathematica transcription of the rate.
            {
                const double eps = 1e-8;

                Parameters pm                      = Parameters::Defaults();
                // FF_2_2_3_3 fixed parameters
                pm["CKM::abs(V_us)"]               = 0.22515;
                pm["0->Kpi::t_0@KSvD2025"]         = -4.0;
                pm["0->Kpi::M_(+,1)@KSvD2025"]     = 1.368;
                pm["0->Kpi::Gamma_(+,1)@KSvD2025"] = 0.212;
                pm["0->Kpi::M_(0,0)@KSvD2025"]     = 0.680;
                pm["0->Kpi::Gamma_(0,0)@KSvD2025"] = 0.600;
                pm["0->Kpi::M_(0,1)@KSvD2025"]     = 1.431;
                pm["0->Kpi::Gamma_(0,1)@KSvD2025"] = 0.220;
                // FF_2_2_3_3 posterior mode
                pm["0->Kpi::b_+^1@KSvD2025"]       = 0.0661020107264485;
                pm["0->Kpi::b_+^2@KSvD2025"]       = 0.07227068199414584;
                pm["0->Kpi::b_+^3@KSvD2025"]       = 0.01799556162000065;
                pm["0->Kpi::b_0^1@KSvD2025"]       = 0.31182977560257646;
                pm["0->Kpi::b_0^2@KSvD2025"]       = 0.2718566992904784;
                pm["0->Kpi::b_0^3@KSvD2025"]       = -0.012034090375745321;
                pm["0->Kpi::M_(+,0)@KSvD2025"]     = 0.8917120142957791;
                pm["0->Kpi::Gamma_(+,0)@KSvD2025"] = 0.04513720400819739;
                pm["decay-constant::K_u"]          = 0.15569999791603167;
                pm["decay-constant::pi"]           = 0.12969354237876146;

                const Options mode_opts{
                    {           "model"_ok, "WET"_ov },
                    {               "K"_ok, "K_S"_ov },
                    { "n-resonances-1m"_ok,   "2"_ov },
                    { "n-resonances-0p"_ok,   "2"_ov }
                };

                // WC set 1: Standard Model (Re{cVL} = 1.009653 encodes the electroweak correction)
                pm["ustaunutau::Re{cVL}"] = 1.009653;
                pm["ustaunutau::Re{cVR}"] = 0.0;
                pm["ustaunutau::Im{cVR}"] = 0.0;
                pm["ustaunutau::Re{cSL}"] = 0.0;
                pm["ustaunutau::Im{cSL}"] = 0.0;
                pm["ustaunutau::Re{cSR}"] = 0.0;
                pm["ustaunutau::Im{cSR}"] = 0.0;
                {
                    TauToKPiNeutrino d(pm, mode_opts);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(0.5), 0.00038675505491877913, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, -0.5), 6.4551233505561254e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, +0.1), 0.00018221129092269793, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, +0.8), 0.00041832741893284617, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.0), 0.0017914256796805357, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, -0.5), 0.00075821001062914735, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, +0.1), 0.00054835709578110253, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, +0.8), 0.0013068205833905338, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.8), 0.00012653597457715373, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, -0.5), 3.8954451930347461e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, +0.1), 5.6589353783935947e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, +0.8), 0.00010805346454769402, eps);
                }

                // WC set 2: complex vector and scalar new physics
                pm["ustaunutau::Re{cVL}"] = 1.0;
                pm["ustaunutau::Re{cVR}"] = 0.10;
                pm["ustaunutau::Im{cVR}"] = 0.05;
                pm["ustaunutau::Re{cSL}"] = 0.20;
                pm["ustaunutau::Im{cSL}"] = 0.10;
                pm["ustaunutau::Re{cSR}"] = -0.08;
                pm["ustaunutau::Im{cSR}"] = 0.03;
                {
                    TauToKPiNeutrino d(pm, mode_opts);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(0.5), 0.00076669281338883515, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, -0.5), 0.00015656563956361816, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, +0.1), 0.00038477527535802172, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, +0.8), 0.0007685899361395093, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.0), 0.0022271696505980986, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, -0.5), 0.0010436006454148028, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, +0.1), 0.00068171777940277776, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, +0.8), 0.0014528626685785728, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.8), 0.00024388025527597514, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, -0.5), 7.4529812209619498e-05, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, +0.1), 0.00011769461891696896, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, +0.8), 0.00020479500226404309, eps);
                }

                // WC set 3: the couplings of WC set 2 plus a complex tensor coefficient. Since no
                // form factor implementation provides f_T yet, this requires debug-tensor = 1, which
                // substitutes B_T = 1. The vector and scalar couplings are unchanged with respect to
                // WC set 2, so that the difference between the two isolates the tensor sector: the
                // tensor-tensor term, the vector-tensor interference in a and c, and the
                // scalar-tensor interference in b.
                const Options mode_opts_tensor{
                    {           "model"_ok, "WET"_ov },
                    {               "K"_ok, "K_S"_ov },
                    { "n-resonances-1m"_ok,   "2"_ov },
                    { "n-resonances-0p"_ok,   "2"_ov },
                    {    "debug-tensor"_ok,   "1"_ov }
                };

                pm["ustaunutau::Re{cVL}"] = 1.0;
                pm["ustaunutau::Re{cVR}"] = 0.10;
                pm["ustaunutau::Im{cVR}"] = 0.05;
                pm["ustaunutau::Re{cSL}"] = 0.20;
                pm["ustaunutau::Im{cSL}"] = 0.10;
                pm["ustaunutau::Re{cSR}"] = -0.08;
                pm["ustaunutau::Im{cSR}"] = 0.03;
                pm["ustaunutau::Re{cT}"]  = 0.15;
                pm["ustaunutau::Im{cT}"]  = -0.06;
                {
                    TauToKPiNeutrino d(pm, mode_opts_tensor);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(0.5), 0.00073574631204537380, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, -0.5), 0.00015177467973205986, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, +0.1), 0.00036836581045930140, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(0.5, +0.8), 0.00073541313105551930, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.0), 0.0027900709232341150, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, -0.5), 0.0013327140281895065, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, +0.1), 0.00097079255254210580, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.0, +0.8), 0.0017173847960770617, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.differential_branching_ratio(1.8), 0.00030612793006224910, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, -0.5), 0.00010263108846277380, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, +0.1), 0.00015614154875673390, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.double_differential_branching_ratio(1.8, +0.8), 0.00023732977893832068, eps);
                }
            }
        }
} tau_to_k_pi_nu_test;
