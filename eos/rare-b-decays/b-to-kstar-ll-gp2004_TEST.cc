/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2021      Méril Reboud
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
#include <eos/rare-b-decays/b-to-kstar-ll-impl.hh>
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
                    {"model"_ok,        "WET"     },
                    {"tag"_ok,          "GP2004"  },
                    {"l"_ok,            "mu"      },
                    {"form-factors"_ok, "BSZ2015" }
                };

                BToKstarDilepton d(p, oo);

                /* q^2 = [14.00, 19.21] */
                {
                    const double eps = 1e-4;
                    auto ir = d.prepare(14.00, 19.21);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(ir), -0.410151, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(ir),   0.315794, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(ir),     -0.548440, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(ir),      1.847569, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(ir),      0.524309, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_5(ir),      0.122853, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_re(ir),    -0.799275, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_im(ir),     0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(ir),                         0.997726, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(ir),                        -0.968699, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(ir),                        -0.955853, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_4(ir),                         0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_5(ir),                         0.0,      eps);

                    double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(ir) / d.integrated_branching_ratio(ir);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(ir), a_fb,    eps);
                }

                /* q^2 = [16.00, 19.21] */
                {
                    const double eps = 1e-4;
                    auto ir = d.prepare(16.00, 19.21);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(ir), -0.374292, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(ir),   0.308243, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(ir),     -0.657588, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(ir),      2.198274, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(ir),      0.439617, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_5(ir),      0.108524, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_re(ir),    -0.721433, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_im(ir),     0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(ir),                         0.999119, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(ir),                        -0.966399, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(ir),                        -0.957599, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_4(ir),                         0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_5(ir),                         0.0,      eps);

                    double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(ir) / d.integrated_branching_ratio(ir);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(ir), a_fb,    eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-7;

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_left),  -8.7364199e-11,   eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_left),  -2.606141641e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_right),  7.15128034e-12,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_right), -2.606141641e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_left),   6.474223364e-11, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_left),   1.931310892e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_right), -5.299537604e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_right),  1.931310892e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_left),  -1.115956981e-10, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_left),  -3.328985969e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_right),  9.134772862e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_right), -3.328985969e-12, eps);
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
                    {"model"_ok,        "WET"     },
                    {"tag"_ok,          "GP2004"  },
                    {"l"_ok,            "mu"      },
                    {"form-factors"_ok, "BSZ2015" }
                };

                BToKstarDilepton d(p, oo);

                /* observables */
                {
                    static const double eps = 1e-4;
                    auto ir = d.prepare(14.18, 19.21);

                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(ir),               2.458835412e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(ir),   -0.4087203569,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation(ir),     0.3149021373,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2(ir),       -0.5572039303,    eps);

                    Kinematics k_mu  = Kinematics({{"q2_min", 14.18}, {"q2_max", 19.21}});
                    auto obs_BR   = Observable::make("B->K^*ll::BR", p, k_mu, oo);
                    auto obs_FL   = Observable::make("B->K^*ll::F_L",  p, k_mu, oo);
                    auto obs_AT2  = Observable::make("B->K^*ll::A_T^2",  p, k_mu, oo);
                    auto obs_AFB  = Observable::make("B->K^*ll::A_FB",  p, k_mu, oo);

                    TEST_CHECK_RELATIVE_ERROR(obs_BR->evaluate(),   2.331870232e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(obs_FL->evaluate(),   0.3149648891,    eps);
                    TEST_CHECK_RELATIVE_ERROR(obs_AT2->evaluate(), -0.556111894,     eps);
                    TEST_CHECK_RELATIVE_ERROR(obs_AFB->evaluate(), -0.4031418964,    eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-7;

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_left),  -2.140027659e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_left),  -8.831641416e-11, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_long_right), -2.140027659e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_long_right),  6.810594153e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_left),   1.585891826e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_left),   6.544788278e-11, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_perp_right),  1.585891826e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_perp_right), -5.047068226e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_left),  -2.733589740e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_left),  -1.128120202e-10, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(amps.a_para_right), -2.733589740e-12, eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(amps.a_para_right),  8.699593306e-12, eps);
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
                    {"model"_ok,        "WET"     },
                    {"tag"_ok,          "GP2004"  },
                    {"l"_ok,            "mu"      },
                    {"form-factors"_ok, "BSZ2015" }
                };

                BToKstarDilepton d(p, oo);

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-7;

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -2.138384054e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -2.604318549e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), -2.138384054e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), -2.604318549e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),   1.584673814e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),   1.929959869e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right),  1.584673814e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right),  1.929959869e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -2.731490261e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -3.326657221e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), -2.731490261e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), -3.326657221e-12, eps);
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
                    { "model"_ok,        "WET"     },
                    { "l"_ok,            "mu"      },
                    { "tag"_ok,          "GP2004"  },
                    { "form-factors"_ok, "BSZ2015" }
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
                options.declare("model"_ok, "WET");
                options.declare("tag"_ok, "GP2004");

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
            o.declare("model"_ok, "WET");
            o.declare("l"_ok, "mu");
            o.declare("q"_ok, "d");
            o.declare("tag"_ok, "GP2004");
            o.declare("form-factors"_ok, "BSZ2015");

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
            oo.declare("model"_ok, "WET");
            oo.declare("scan-mode"_ok, "cartesian");
            oo.declare("tag"_ok, "GP2004");
            oo.declare("form-factors"_ok, "KMPW2010");
            oo.declare("l"_ok, "mu");
            oo.declare("q"_ok, "d");

            static const double q2 = 14.0;
            static const double q2_max = 19.0;

            {
                double eps = 7.5e-3;

                BToKstarDilepton d(p, oo);
                auto amps = d.amplitudes(q2);

                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.022237305e-11, -1.565547054e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 8.648694084e-11,  8.335187033e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-7.452197142e-11,  4.596972947e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 3.56948263e-11,   2.58076248e-11 ), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-7.540711414e-11,  4.651573987e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 3.611879543e-11,  2.611415763e-11), eps);

                // nearly identically implemented, only difference from alpha_s
                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-1.678438392e-10, -3.508952983e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.768920882e-12,  2.768920882e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 2.381028679e-11,  3.703822389e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>( 2.645587421e-11,  3.968381131e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 1.469979946e-11,  2.286635471e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>( 1.633311051e-11,  2.449966576e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>( 3.124004517e-11,  4.686006776e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.811604066e-11,  4.373606324e-11), eps);

                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1s(q2),  6.085019472e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1c(q2),  4.489187687e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2s(q2), -6.421382622e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2c(q2),  6.967884897e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_3(q2),   1.032360198e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_4(q2),  -9.221388938e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_5(q2),  -5.020415785e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6s(q2), -2.409650924e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6c(q2), -1.458618418e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_7(q2),  -2.414995565e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_8(q2),  -3.305986949e-22, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_9(q2),   4.620491199e-22, eps);

                TEST_CHECK_RELATIVE_ERROR(d.differential_forward_backward_asymmetry(q2),      -0.1862325546, eps);

                auto ir = d.prepare(q2, q2_max);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(ir), -0.1855329818, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_1(ir), -1.004548102,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_2(ir), -0.6518372271, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_3(ir), -1.553829809,  eps);
            }

            {
                oo.declare("cp-conjugate"_ok, "true");
                BToKstarDilepton d(p, oo);

                auto amps = d.amplitudes(q2);

                double eps = 7e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.022237305e-11,  1.840549748e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 8.648694084e-11, -8.060184339e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-7.452197142e-11, -8.434606297e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 3.56948263e-11,  -2.964525815e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-7.540711414e-11, -8.534789239e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 3.611879543e-11, -2.999737289e-11), eps);

                // nearly identically implemented, only difference from alpha_s
                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-1.678438392e-10,  3.508952983e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.768920882e-12, -2.768920882e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 2.381028679e-11, -3.703822389e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>( 2.645587421e-11, -3.968381131e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 1.469979946e-11, -2.286635471e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>( 1.633311051e-11, -2.449966576e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>( 3.124004517e-11, -4.686006776e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.811604066e-11, -4.373606324e-11), eps);

                eps = 1e-3;
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1s(q2),  6.158136673e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_1c(q2),  4.670611082e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2s(q2), -6.427703401e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_2c(q2),  6.895804905e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_3(q2),   1.024056268e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_4(q2),  -9.171018962e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_5(q2),  -5.232104087e-21, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6s(q2), -2.43848199e-20,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_6c(q2), -1.457072802e-20, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_7(q2),   2.418865485e-21, eps);

                eps = 3.5e-3;
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_8(q2),   2.57868544e-22,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_j_9(q2),  -3.604004966e-22, eps);

                // j6c needed, nonzero only with tensor contributions
                TEST_CHECK_RELATIVE_ERROR(d.differential_forward_backward_asymmetry(q2),      -0.1842839266,  eps);

                auto ir = d.prepare(q2, q2_max);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(ir), -0.1816844542, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_1(ir), -1.004836959,  eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_2(ir), -0.6691776451, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_h_3(ir), -1.53250009,  eps);
            }

            {
                Options o{
                    {"model"_ok, "WET"},
                    {"scan-mode"_ok, "cartesian"},
                    {"tag"_ok, "GP2004"},
                    {"form-factors"_ok, "KMPW2010"},
                    {"l"_ok, "mu"},
                    {"q"_ok, "d"}
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
