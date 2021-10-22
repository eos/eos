/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/maths/complex.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

// Uncomment the following #define to generate new test data for the Bobeth compatibility tests
//#define EOS_GENERATE_TEST_DATA
#ifdef EOS_GENERATE_TEST_DATA
#include <gsl/gsl_rng.h>
#endif

using namespace test;
using namespace eos;

class BToKstarDileptonLowRecoilTest :
    public TestCase
{
    public:
        BToKstarDileptonLowRecoilTest() :
            TestCase("b_to_kstar_dilepton_low_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil (SM) */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.530e-12;
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["sb::mu"] = 4.2;
                p["b->s::Re{c7}"] = -0.331;
                p["b->s::c8"] = -0.18100000;
                p["sbmumu::mu"] = 4.2;
                p["b->smumu::Re{c9}"] = +4.27;
                p["b->smumu::Re{c10}"] = -4.173;
                // PDG 2008 CKM parameters
                p["CKM::A"]         =  0.814;
                p["CKM::lambda"]    =  0.2257;
                p["CKM::rhobar"]    =  0.135;
                p["CKM::etabar"]    =  0.349;
                p["CKM::abs(V_ub)"] =  0.00359255775926898;
                p["CKM::arg(V_ub)"] = -1.2023040533144056;
                p["CKM::abs(V_cb)"] =  0.04146529127297828;
                p["CKM::arg(V_cb)"] =  0.0;
                p["CKM::abs(V_tb)"] =  0.9991334809397352;
                p["CKM::arg(V_tb)"] =  0.0;
                p["CKM::abs(V_us)"] =  0.22569854350471902;
                p["CKM::arg(V_us)"] =  0.0;
                p["CKM::abs(V_cs)"] =  0.973346862850555;
                p["CKM::arg(V_cs)"] = -3.222382085887583e-05;
                p["CKM::abs(V_ts)"] =  0.040694467854567457;
                p["CKM::arg(V_ts)"] = -3.1230200317017145;
                // Kaon mass
                p["mass::K_d^*"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // muon mass near zero to avoid artificial divergence
                p["mass::mu"] = 1e-5;

                Options oo
                {
                    {"model", "WET"},
                    {"tag", "GP2004"},
                    {"l", "mu"},
                    {"form-factors", "BZ2004"}
                };

                BToKstarDilepton d(p, oo);

                /* q^2 = [14.00, 19.21] */
                {
                    const double eps = 1e-4;

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(14.00, 19.21), -0.4093, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(14.00, 19.21),  +0.3497, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(14.00, 19.21),     -0.4835, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(14.00, 19.21),     +1.6893, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(14.00, 19.21),     +0.5758, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_5(14.00, 19.21),      0.1244, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_re(14.00, 19.21),    -0.8391, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_im(14.00, 19.21),     0.0,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(14.00, 19.21),                        +0.9967, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(14.00, 19.21),                        -0.9727, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(14.00, 19.21),                        -0.9587, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_4(14.00, 19.21),                         0.0,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_5(14.00, 19.21),                         0.0,    eps);

                    double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(14.00, 19.21) / d.integrated_branching_ratio(14.00, 19.21);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(14.00, 19.21), a_fb,    eps);
                }

                /* q^2 = [16.00, 19.21] */
                {
                    const double eps = 1e-4;

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(16.00, 19.21), -0.381708, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(16.00, 19.21),  +0.337697, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(16.00, 19.21),     -0.599389, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(16.00, 19.21),     +1.99535,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(16.00, 19.21),     +0.486256, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_5(16.00, 19.21),      0.112158, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_re(16.00, 19.21),    -0.768382, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_im(16.00, 19.21),     0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(16.00, 19.21),                        +0.998622, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(16.00, 19.21),                        -0.970214, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(16.00, 19.21),                        -0.959887, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_4(16.00, 19.21),                         0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_5(16.00, 19.21),                         0.0,      eps);

                    double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(16.00, 19.21) / d.integrated_branching_ratio(16.00, 19.21);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(16.00, 19.21), a_fb,    eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -9.860564941316e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -2.941484608501e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), +8.071641897174e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), -2.941484608501e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),  +7.179697602811e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),  +2.141760651448e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right), -5.877142772730e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right), +2.141760651448e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -1.139839686524e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -3.400232049605e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), +9.330477335285e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), -3.400232049605e-12, eps);
                }
            }

            /* Low Recoil (Benchmark Point) */
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["sb::mu"] = 4.2;
                p["b->s::Re{c7}"] = 0.0;
                p["b->s::Im{c7}"] = -0.3;
                p["b->s::c8"] = -0.181;
                p["sbmumu::mu"] = 4.2;
                p["b->smumu::Re{c9}"] = 0.0;
                p["b->smumu::Im{c9}"] = 4.2;
                p["b->smumu::Re{c10}"] = 0.0;
                p["b->smumu::Im{c10}"] = -4.2;
                // PDG 2008 CKM parameters
                p["CKM::A"]         = 0.814;
                p["CKM::lambda"]    = 0.2257;
                p["CKM::rhobar"]    = 0.135;
                p["CKM::etabar"]    = 0.349;
                p["CKM::abs(V_ub)"] =  0.00359255775926898;
                p["CKM::arg(V_ub)"] = -1.2023040533144056;
                p["CKM::abs(V_cb)"] =  0.04146529127297828;
                p["CKM::arg(V_cb)"] =  0.0;
                p["CKM::abs(V_tb)"] =  0.9991334809397352;
                p["CKM::arg(V_tb)"] =  0.0;
                p["CKM::abs(V_us)"] =  0.22569854350471902;
                p["CKM::arg(V_us)"] =  0.0;
                p["CKM::abs(V_cs)"] =  0.973346862850555;
                p["CKM::arg(V_cs)"] = -3.222382085887583e-05;
                p["CKM::abs(V_ts)"] =  0.040694467854567457;
                p["CKM::arg(V_ts)"] = -3.1230200317017145;
                // Kaon mass
                p["mass::K_d^*"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // mu mass
                p["mass::mu"] = 1e-5;

                Options oo
                {
                    {"model", "WET"},
                    {"tag", "GP2004"},
                    {"l", "mu"},
                    {"form-factors", "BZ2004"}
                };

                BToKstarDilepton d(p, oo);

                /* observables */
                {
                    static const double eps = 1e-4;
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(14.18, 19.21),                         2.77767e-7, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(14.18, 19.21),             -4.08565e-1, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation(14.18, 19.21),              +0.34841,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2(14.18, 19.21),                 -4.92697e-1, eps);

                    Kinematics k_mu  = Kinematics({{"q2_min", 14.18}, {"q2_max", 19.21}});
                    auto obs_BR   = Observable::make("B->K^*ll::BR", p, k_mu, oo);
                    auto obs_FL   = Observable::make("B->K^*ll::F_L",  p, k_mu, oo);
                    auto obs_AT2  = Observable::make("B->K^*ll::A_T^2",  p, k_mu, oo);
                    auto obs_AFB  = Observable::make("B->K^*ll::A_FB",  p, k_mu, oo);

                    TEST_CHECK_RELATIVE_ERROR(obs_BR->evaluate(),   2.63504e-7, eps);
                    TEST_CHECK_RELATIVE_ERROR(obs_FL->evaluate(),   0.34851,    eps);
                    TEST_CHECK_RELATIVE_ERROR(obs_AT2->evaluate(), -4.91581e-1, eps);
                    TEST_CHECK_RELATIVE_ERROR(obs_AFB->evaluate(), -4.02902e-1, eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -2.41522826885e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -9.96805582174e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), -2.41522826886e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), +7.68695280669e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),  +1.75858165484e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),  +7.25796411402e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right), +1.75858165484e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right), -5.59704205262e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -2.79190193386e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -1.15226517859e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), -2.79190193386e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), +8.88579298412e-12, eps);
                }
            }

            /* Low Recoil (Zero Point for C_7 = C_9 = C_10 = 0) */
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["sb::mu"] = 4.2;
                p["b->s::Re{c7}"] = 0.0;
                p["b->s::c8"] = -0.181;
                p["sbmumu::mu"] = 4.2;
                p["b->smumu::Re{c9}"] = 0.0;
                p["b->smumu::Re{c10}"] = 0.0;
                // PDG 2008 CKM parameters
                p["CKM::A"]         = 0.814;
                p["CKM::lambda"]    = 0.2257;
                p["CKM::rhobar"]    = 0.135;
                p["CKM::etabar"]    = 0.349;
                p["CKM::abs(V_ub)"] =  0.00359255775926898;
                p["CKM::arg(V_ub)"] = -1.2023040533144056;
                p["CKM::abs(V_cb)"] =  0.04146529127297828;
                p["CKM::arg(V_cb)"] =  0.0;
                p["CKM::abs(V_tb)"] =  0.9991334809397352;
                p["CKM::arg(V_tb)"] =  0.0;
                p["CKM::abs(V_us)"] =  0.22569854350471902;
                p["CKM::arg(V_us)"] =  0.0;
                p["CKM::abs(V_cs)"] =  0.973346862850555;
                p["CKM::arg(V_cs)"] = -3.222382085887583e-05;
                p["CKM::abs(V_ts)"] =  0.040694467854567457;
                p["CKM::arg(V_ts)"] = -3.1230200317017145;
                // Kaon mass
                p["mass::K_d^*"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;

                Options oo
                {
                    {"model", "WET"},
                    {"tag", "GP2004"},
                    {"l", "mu"},
                    {"form-factors", "BZ2004"}
                };

                BToKstarDilepton d(p, oo);

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -2.413541335202e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -2.939430107299e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), -2.413541335202e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), -2.939430107299e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),  +1.757353360762e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),  +2.140264723229e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right), +1.757353360762e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right), +2.140264723229e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -2.789951909754e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -3.397857132935e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), -2.789951909754e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), -3.397857132935e-12, eps);
                }
            }
        }
} b_to_kstar_dilepton_low_recoil_test;

class BToKstarDileptonLowRecoilPolynomialTest :
    public TestCase
{
    public:
        BToKstarDileptonLowRecoilPolynomialTest() :
            TestCase("b_to_kstar_dilepton_low_recoil_polynomial_test")
        {
        }

        void run_one(const ObservablePtr & o, const WilsonPolynomial & p, const std::array<double, 6> & values) const
        {
            Parameters parameters = o->parameters();
            Parameter abs_c7(parameters["b->s::Re{c7}"]);
            Parameter arg_c7(parameters["b->s::Im{c7}"]);
            Parameter abs_c9(parameters["b->smumu::Re{c9}"]);
            Parameter arg_c9(parameters["b->smumu::Im{c9}"]);
            Parameter abs_c10(parameters["b->smumu::Re{c10}"]);
            Parameter arg_c10(parameters["b->smumu::Im{c10}"]);

            abs_c7 = values[0];
            arg_c7 = values[1];
            abs_c9 = values[2];
            arg_c9 = values[3];
            abs_c10 = values[4];
            arg_c10 = values[5];

            static const double eps = 3e-14;
            WilsonPolynomialEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(), p.accept_returning<double>(evaluator), eps);
        }

        virtual void run() const
        {
            // Test make_polynomial
            {
                static const std::vector<std::string> names
                {
                    "B->K^*ll::BR",
                    "B->K^*ll::Abar_FB",
                };
                static const std::vector<std::array<double, 6>> inputs
                {
                    std::array<double, 6>{{0.0,       0.0,       0.0,       0.0,       0.0,       0.0      }},
                    std::array<double, 6>{{1.0,       0.0,       1.0,       0.0,       1.0,       0.0      }},
                    std::array<double, 6>{{0.7808414, 0.8487257, 0.7735165, 0.5383695, 0.6649164, 0.7235497}},
                    std::array<double, 6>{{0.5860642, 0.9830907, 0.7644369, 0.8330194, 0.4935018, 0.4492084}},
                    std::array<double, 6>{{0.2177456, 0.5062894, 0.6463376, 0.3624364, 0.6770480, 0.0718421}},
                    std::array<double, 6>{{0.0088306, 0.9441413, 0.8721501, 0.2984633, 0.2961408, 0.9145809}},
                    std::array<double, 6>{{0.7967655, 0.2427081, 0.8403112, 0.3351082, 0.6477823, 0.5569495}},
                    std::array<double, 6>{{0.7607454, 0.5025871, 0.5877762, 0.5516025, 0.2930899, 0.4882813}},
                };

                Parameters parameters = Parameters::Defaults();
                parameters["CKM::abs(V_ub)"] =  0.003631275231633653;
                parameters["CKM::arg(V_ub)"] = -1.210765774253535;
                parameters["CKM::abs(V_cb)"] =  0.041996951916414726;
                parameters["CKM::arg(V_cb)"] =  0.0;
                parameters["CKM::abs(V_tb)"] =  0.9991111344469873;
                parameters["CKM::arg(V_tb)"] =  0.0;
                parameters["CKM::abs(V_us)"] =  0.22534851424944366;
                parameters["CKM::arg(V_us)"] =  0.0;
                parameters["CKM::abs(V_cs)"] =  0.9734061815416853;
                parameters["CKM::arg(V_cs)"] = -3.304199362533668e-05;
                parameters["CKM::abs(V_ts)"] =  0.04121212396309175;
                parameters["CKM::arg(V_ts)"] = -3.1230250224697222;
                Kinematics kinematics
                {
                    { "q2_min", 14.18 },
                    { "q2_max", 19.21 }
                };
                Options options
                {
                    { "model",        "WET" },
                    { "l",            "mu"         },
                    { "tag",          "GP2004"     },
                    { "form-factors", "BZ2004"     }
                };

                for (const auto & name : names)
                {
                    ObservablePtr observable = Observable::make(name, parameters, kinematics, options);
                    WilsonPolynomial polynomial = make_polynomial(observable, std::list<std::string>{ "b->s::Re{c7}", "b->s::Im{c7}", "b->smumu::Re{c9}", "b->smumu::Im{c9}", "b->smumu::Re{c10}", "b->smumu::Im{c10}" });

                    for (const auto & input : inputs)
                    {
                        run_one(observable, polynomial, input);
                    }
                }
            }

            // Test ratios
            {
                static const double eps = 1e-7;
                Kinematics kinematics
                {
                    { "q2_min", 14.18 },
                    { "q2_max", 19.21 }
                };

                Parameters parameters = Parameters::Defaults();
                parameters["life_time::B_d"] = 1.530e-12;
                Parameter lambda = parameters["CKM::lambda"];
                Parameter A = parameters["CKM::A"];

                Options options;
                options.declare("model", "WET");
                options.declare("tag", "GP2004");

                ObservablePtr numerator = Observable::make("B->K^*ll::Abar_FB", parameters, kinematics, options);
                ObservablePtr denominator = Observable::make("B->K^*ll::BR", parameters, kinematics, options);
                ObservablePtr observable = Observable::make("B->K^*ll::A_FB", parameters, kinematics, options);

                TEST_CHECK_NEARLY_EQUAL(numerator->evaluate() / denominator->evaluate(), observable->evaluate(), eps);

                // vary CKM::lambda
                {
                    lambda = lambda.max();
                    TEST_CHECK_NEARLY_EQUAL(numerator->evaluate() / denominator->evaluate(), observable->evaluate(), eps);
                    lambda = lambda.min();
                    TEST_CHECK_NEARLY_EQUAL(numerator->evaluate() / denominator->evaluate(), observable->evaluate(), eps);
                    lambda = lambda.central();
                }

                // vary CKM::A
                {
                    A = A.max();
                    TEST_CHECK_NEARLY_EQUAL(numerator->evaluate() / denominator->evaluate(), observable->evaluate(), eps);
                    A = A.min();
                    TEST_CHECK_NEARLY_EQUAL(numerator->evaluate() / denominator->evaluate(), observable->evaluate(), eps);
                    A = A.central();
                }

                std::list<std::string> coefficients{"b->s::Re{c7}", "b->s::Im{c7}", "b->smumu::Re{c9}", "b->smumu::Im{c9}", "b->smumu::Re{c10}", "b->smumu::Im{c10}"};

                // central ratio
                {
                    ObservablePtr ratio = make_polynomial_ratio(
                            make_polynomial(numerator, coefficients),
                            make_polynomial(denominator, coefficients),
                            parameters);
                    TEST_CHECK_NEARLY_EQUAL(ratio->evaluate(), observable->evaluate(), eps);
                }

                // lambda ratios
                {
                    lambda = lambda.max();
                    ObservablePtr ratio = make_polynomial_ratio(
                            make_polynomial(numerator, coefficients),
                            make_polynomial(denominator, coefficients),
                            parameters);
                    TEST_CHECK_NEARLY_EQUAL(ratio->evaluate(), observable->evaluate(), eps);
                }
                {
                    lambda = lambda.min();
                    ObservablePtr ratio = make_polynomial_ratio(
                            make_polynomial(numerator, coefficients),
                            make_polynomial(denominator, coefficients),
                            parameters);
                    TEST_CHECK_NEARLY_EQUAL(ratio->evaluate(), observable->evaluate(), eps);
                }

                // A ratios
                {
                    A = A.max();
                    ObservablePtr ratio = make_polynomial_ratio(
                            make_polynomial(numerator, coefficients),
                            make_polynomial(denominator, coefficients),
                            parameters);
                    TEST_CHECK_NEARLY_EQUAL(ratio->evaluate(), observable->evaluate(), eps);
                }
                {
                    A = A.min();
                    ObservablePtr ratio = make_polynomial_ratio(
                            make_polynomial(numerator, coefficients),
                            make_polynomial(denominator, coefficients),
                            parameters);
                    TEST_CHECK_NEARLY_EQUAL(ratio->evaluate(), observable->evaluate(), eps);
                }
            }
        }
} b_to_kstar_dilepton_low_recoil_polynomial_test;

class BToKstarDileptonLowRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKstarDileptonLowRecoilBobethCompatibilityTest() :
            TestCase("b_to_kstar_dilepton_low_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            static const std::vector<std::string> variation_names
            {
                "b->s::Re{c7}",      "b->s::Im{c7}",      "b->s::Re{c7'}",      "b->s::Im{c7'}",
                "b->smumu::Re{c9}",  "b->smumu::Im{c9}",  "b->smumu::Re{c9'}",  "b->smumu::Im{c9'}",
                "b->smumu::Re{c10}", "b->smumu::Im{c10}", "b->smumu::Re{c10'}", "b->smumu::Im{c10'}",
            };

            Parameters p = Parameters::Defaults();
            // comparison done for zero lepton mass
            // but this leads to a NaN in the timelike transversity amplitude
            // so make the mass very small
            p["mass::mu"] = 1e-5;
            p["mass::B_d"] = 5.27958;
            p["mass::K_d^*"] = 0.89594;
            p["CKM::abs(V_ub)"] =  0.003631275231633653;
            p["CKM::arg(V_ub)"] = -1.210765774253535;
            p["CKM::abs(V_cb)"] =  0.041996951916414726;
            p["CKM::arg(V_cb)"] =  0.0;
            p["CKM::abs(V_tb)"] =  0.9991111344469873;
            p["CKM::arg(V_tb)"] =  0.0;
            p["CKM::abs(V_us)"] =  0.22534851424944366;
            p["CKM::arg(V_us)"] =  0.0;
            p["CKM::abs(V_cs)"] =  0.9734061815416853;
            p["CKM::arg(V_cs)"] = -3.304199362533668e-05;
            p["CKM::abs(V_ts)"] =  0.04121212396309175;
            p["CKM::arg(V_ts)"] = -3.1230250224697222;

            Options o;
            o.declare("model", "WET");
            o.declare("l", "mu");
            o.declare("q", "d");
            o.declare("tag", "GP2004");
            o.declare("form-factors", "BZ2004");

            std::vector<Parameter> variations;

            for (const auto & variation_name : variation_names)
            {
                variations.push_back(p[variation_name]);
            }

            Kinematics k
            {
                { "q2_min", 14.18 },
                { "q2_max", 19.21 }
            };

            std::vector<ObservablePtr> observables;
            std::vector<std::string> observable_names = {
                    "B->K^*ll::BR_CP_specific",
                    "B->K^*ll::A_FB_CP_specific",
                    "B->K^*ll::F_L_CP_specific"
            };
            for (auto & s : observable_names)
            {
                observables.push_back(Observable::make(s, p, k, o));
                TEST_CHECK_MSG(observables.back(), "Could not create '" + s + "'");
            }

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil_TEST-btokstarll.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->K^*ll at low recoil --" << std::endl;
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus2);
                std::fstream file(filename.c_str(), std::fstream::out);
                file.precision(17);

                for (int i = 0 ; i < 1000 ; ++i)
                {
                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        *v = v->min() + (v->max() - v->min()) * gsl_rng_uniform(rng);

                        file << *v << '\t';
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        file << (*o)->evaluate() << '\t';
                    }
                    file << std::endl;
                }
            }
#else
            // Verify the test case data
            {
                std::cout << "-- Verifying test case data for B->K^*ll at low recoil --" << std::endl;
                std::fstream file(filename.c_str(), std::fstream::in);
                TEST_CHECK_MSG(file, "'" + filename + "' does not exist");
                std::string line;
                while (file)
                {
                    std::getline(file, line);
                    if (line.empty())
                        break;

                    std::stringstream ss(line);

                    for (auto & variation : variations)
                    {
                        double value;
                        ss >> value;
                        variation = value;
                    }

                    for (const auto & observable : observables)
                    {
                        double reference;
                        ss >> reference;

                        TEST_CHECK_RELATIVE_ERROR(reference, observable->evaluate(), 1e-3);
                    }
                }
            }
#endif
        }
} b_to_kstar_dilepton_low_recoil_bobeth_compatibility_test;

class BToKstarDileptonTensorLowRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKstarDileptonTensorLowRecoilBobethCompatibilityTest() :
            TestCase("b_to_kstar_tensor_dilepton_low_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            // Christoph uses \Delta C instead of C for C9, C10
            // important to agree to alpha_s, can change values by 1%
            Parameters p = Parameters::Defaults();
            p["CKM::abs(V_ub)"] =  0.003631275231633653;
            p["CKM::arg(V_ub)"] = -1.210765774253535;
            p["CKM::abs(V_cb)"] =  0.041996951916414726;
            p["CKM::arg(V_cb)"] =  0.0;
            p["CKM::abs(V_tb)"] =  0.9991111344469873;
            p["CKM::arg(V_tb)"] =  0.0;
            p["CKM::abs(V_us)"] =  0.22534851424944366;
            p["CKM::arg(V_us)"] =  0.0;
            p["CKM::abs(V_cs)"] =  0.9734061815416853;
            p["CKM::arg(V_cs)"] = -3.304199362533668e-05;
            p["CKM::abs(V_ts)"] =  0.04121212396309175;
            p["CKM::arg(V_ts)"] = -3.1230250224697222;
            p["b->s::c1"] = -0.3231323312;
            p["b->s::c2"] = 1.009301831;
            p["b->s::c3"] = -0.005233499106;
            p["b->s::c4"] = -0.08829686414;
            p["b->s::c5"] = 0.0003601965805;
            p["b->s::c6"] = 0.001020749573;
            p["sb::mu"] = 4.2;
            p["b->s::Re{c7}"] = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"] = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"] = -0.1827530948;
            p["sbmumu::mu"] = 4.2;
            p["b->smumu::Re{c9}"] = 4.294489364 + 1;
            p["b->smumu::Im{c9}"] = 0.5;
            p["b->smumu::Re{c9'}"] = 2;
            p["b->smumu::Im{c9'}"] = 1.5;
            p["b->smumu::Re{c10}"] = -4.196294696 + 3;
            p["b->smumu::Im{c10}"] = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;
            p["b->smumu::Re{cS}"] = 0.5;
            p["b->smumu::Im{cS}"] = 1;
            p["b->smumu::Re{cS'}"] = 0.6;
            p["b->smumu::Im{cS'}"] = 1.1;
            p["b->smumu::Re{cP}"] = 0.7;
            p["b->smumu::Im{cP}"] = 1.2;
            p["b->smumu::Re{cP'}"] = 0.8;
            p["b->smumu::Im{cP'}"] = 1.3;
            p["b->smumu::Re{cT}"] = 0.9;
            p["b->smumu::Im{cT}"] = 1.4;
            p["b->smumu::Re{cT5}"] = -1.0;
            p["b->smumu::Im{cT5}"] = -1.5;

            p["mass::s(2GeV)"] = 0.095;

            // increase sensitivity to m_l^2/q^2 terms
            p["mass::mu"] = 1.5;

            Options oo;
            oo.declare("model", "WET");
            oo.declare("scan-mode", "cartesian");
            oo.declare("tag", "GP2004");
            oo.declare("form-factors", "KMPW2010");
            oo.declare("l", "mu");
            oo.declare("q", "d");

            static const double q2 = 14.0;
            static const double q2_max = 19.0;

            {
                double eps = 7.5e-3;

                BToKstarDilepton d(p, oo);
                auto amps = d.amplitudes(q2);

                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.021407965e-11, -1.564297789e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 8.645626526e-11,  8.331646455e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-7.455049449e-11,  4.565517978e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 3.565928422e-11,  2.577481906e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-7.541145186e-11,  4.618243535e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 3.607110071e-11,  2.607248335e-11), eps);

                // nearly identically implemented, only difference from alpha_s
                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-1.677697256e-10, -3.507403558e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.767698228e-12,  2.767698228e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 2.38060e-11,      3.70316e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>( 2.64511e-11,      3.96767e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 1.46932e-11,      2.28561e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>( 1.63258e-11,      2.44887e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>( 3.12340e-11,      4.6851e-11 ),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.81106e-11,      4.37276e-11),     eps);

                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1s(q2),  6.080153751e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1c(q2),  4.48478951e-20,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2s(q2), -6.418495462e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2c(q2),  6.966335387e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_3(q2),   1.032048382e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_4(q2),  -9.218261443e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_5(q2),  -5.01299633e-21,  10*eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6s(q2), -2.407918046e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6c(q2), -1.457760738e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_7(q2),  -2.413388446e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_8(q2),  -3.305335877e-22, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_9(q2),   4.621083562e-22, eps);

                TEST_CHECK_RELATIVE_ERROR(d.differential_forward_backward_asymmetry(q2),      -0.1862620012, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(q2, q2_max), -0.1855108351, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_1(q2, q2_max), -1.004538219,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_2(q2, q2_max), -0.6513218755, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_3(q2, q2_max), -1.553614345,  eps);
            }

            {
                oo.declare("cp-conjugate", "true");
                BToKstarDilepton d(p, oo);

                auto amps = d.amplitudes(q2);

                double eps = 7e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.021407965e-11,  1.843164004e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 8.645626526e-11, -8.05278024e-11 ), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-7.455049449e-11, -8.452349138e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 3.565928422e-11, -2.966165022e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-7.541145186e-11, -8.549962337e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 3.607110071e-11, -3.000420215e-11), eps);

                // nearly identically implemented, only difference from alpha_s
                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-1.677697256e-10,  3.507403558e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.767698228e-12, -2.767698228e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 2.3806e-11,      -3.70316e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>( 2.64511e-11,     -3.96767e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 1.46932e-11,     -2.28561e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>( 1.63258e-11,     -2.44887e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>( 3.1234e-11,      -4.6851e-11 ),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.81106e-11,     -4.37276e-11),     eps);

                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1s(q2),  6.154137843e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1c(q2),  4.668428684e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2s(q2), -6.424911528e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2c(q2),  6.893410893e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_3(q2),   1.023638282e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_4(q2),  -9.167284751e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_5(q2),  -5.227165173e-21, 10*eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6s(q2), -2.437095943e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6c(q2), -1.456196508e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_7(q2),   2.417311002e-21, eps);

                eps = 3.5e-3;
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_8(q2),   2.569833906e-22, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_9(q2),  -3.592801961e-22, eps);

                // j6c needed, nonzero only with tensor contributions
                TEST_CHECK_RELATIVE_ERROR(d.differential_forward_backward_asymmetry(q2),      -0.184288085,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(q2, q2_max), -0.1816297646, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_1(q2, q2_max), -1.004829116,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_2(q2, q2_max), -0.6687996564, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_3(q2, q2_max), -1.532088081,  eps);
            }

            {
                Options o{
                    {"model", "WET"},
                    {"scan-mode", "cartesian"},
                    {"tag", "GP2004"},
                    {"form-factors", "KMPW2010"},
                    {"l", "mu"},
                    {"q", "d"}
                };
                Kinematics k{{"q2", q2}};

                auto observables = Observables();
                observables.insert("B->K^*ll::J_6c_cp_averaged(q2)", R"()", Unit::None(), Options(),
                R"(
                0.5 * (<<B->K^*ll::J_6c(q2);cp-conjugate=false>> + <<B->K^*ll::J_6c(q2);cp-conjugate=true>>)
                )");
                observables.insert("B->K^*ll::J1c_plus_J2c_cp_averaged(q2)", R"()", Unit::None(), Options(),
                R"(
                0.5 * (<<B->K^*ll::J_1c(q2);cp-conjugate=false>> + <<B->K^*ll::J_1c(q2);cp-conjugate=true>>
                     + <<B->K^*ll::J_2c(q2);cp-conjugate=false>> + <<B->K^*ll::J_2c(q2);cp-conjugate=true>>)
                )");
                observables.insert("B->K^*ll::J1s_minus_3J2s_cp_averaged(q2)", R"()", Unit::None(), Options(),
                R"(
                0.5 * (<<B->K^*ll::J_1s(q2);cp-conjugate=false>> + <<B->K^*ll::J_1s(q2);cp-conjugate=true>>)
                -
                1.5 * (<<B->K^*ll::J_2s(q2);cp-conjugate=false>> + <<B->K^*ll::J_2s(q2);cp-conjugate=true>>)
                )");
                auto J6c_avg            = Observable::make("B->K^*ll::J_6c_cp_averaged(q2)", p, k, o);
                auto J1c_plus_J2c_avg   = Observable::make("B->K^*ll::J1c_plus_J2c_cp_averaged(q2)", p, k, o);
                auto J1s_minus_3J2s_avg = Observable::make("B->K^*ll::J1s_minus_3J2s_cp_averaged(q2)", p, k, o);

                double eps = 7e-4;
                TEST_CHECK_RELATIVE_ERROR(J6c_avg->evaluate(), 0.5 * (-1.457760738e-20 -1.456196508e-20), eps);
                TEST_CHECK_RELATIVE_ERROR(J1c_plus_J2c_avg->evaluate(),
                        0.5 * (4.48478951e-20 + 4.668428684e-20 + 6.966335387e-21 + 6.893410893e-21), eps);
                TEST_CHECK_RELATIVE_ERROR(J1s_minus_3J2s_avg->evaluate(),
                        0.5 * (6.080153751e-20 + 6.154137843e-20 - 3.0 * (-6.418495462e-21 - 6.424911528e-21)), eps);
            }
        }
} b_to_kstar_dilepton_tensor_low_recoil_bobeth_compatibility_test;
