/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

// Uncomment the following #define to generate new test data for the Bobeth compatibility tests
//#define EOS_GENERATE_TEST_DATA

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
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.18100000;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // s quark mass
                p["mass::mu"] = 0.0;


                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LowRecoil> d(p, oo);

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
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  16.00)), -9.860564941316e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  16.00)), -2.941484608501e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 16.00)), +8.071641897174e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 16.00)), -2.941484608501e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  16.00)), +7.179697602811e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  16.00)), +2.141760651448e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 16.00)), -5.877142772730e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 16.00)), +2.141760651448e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   16.00)), -1.139839686524e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   16.00)), -3.400232049605e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  16.00)), +9.330477335285e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  16.00)), -3.400232049605e-12, eps);
                }
            }

            /* Low Recoil (Benchmark Point) */
            {
                Parameters p = Parameters::Defaults();
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.3;
                p["Arg{c7}"] = -M_PI / 2.0;
                p["c8"] = -0.181;
                p["Abs{c9}"] = 4.2;
                p["Arg{c9}"] = +M_PI / 2.0;
                p["Abs{c10}"] = 4.2;
                p["Arg{c10}"] = -M_PI / 2.0;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // mu mass
                p["mass::mu"] = 0.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LowRecoil> d(p, oo);

                /* observables */
                {
                    static const double eps = 1e-4;
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(14.18, 19.21),                         2.77767e-7, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(14.18, 19.21),             2.63504e-7, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(14.18, 19.21),             -4.08565e-1, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry_cp_averaged(14.18, 19.21), -4.02902e-1, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation(14.18, 19.21),              +0.34841,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation_cp_averaged(14.18, 19.21),  +0.34851,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2(14.18, 19.21),                 -4.92697e-1, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2_cp_averaged(14.18, 19.21),     -4.91581e-1, eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  16.00)), -2.41522826885e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  16.00)), -9.96805582174e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 16.00)), -2.41522826886e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 16.00)), +7.68695280669e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  16.00)), +1.75858165484e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  16.00)), +7.25796411402e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 16.00)), +1.75858165484e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 16.00)), -5.59704205262e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   16.00)), -2.79190193386e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   16.00)), -1.15226517859e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  16.00)), -2.79190193386e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  16.00)), +8.88579298412e-12, eps);
                }
            }

            /* Low Recoil (Zero Point for C_7 = C_9 = C_10 = 0) */
            {
                Parameters p = Parameters::Defaults();
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.0;
                p["c8"] = -0.181;
                p["Abs{c9}"] = 0.0;
                p["Abs{c10}"] = 0.0;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LowRecoil> d(p, oo);

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  16.00)), -2.413541335202e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  16.00)), -2.939430107299e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 16.00)), -2.413541335202e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 16.00)), -2.939430107299e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  16.00)), +1.757353360762e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  16.00)), +2.140264723229e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 16.00)), +1.757353360762e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 16.00)), +2.140264723229e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   16.00)), -2.789951909754e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   16.00)), -3.397857132935e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  16.00)), -2.789951909754e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  16.00)), -3.397857132935e-12, eps);
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
            Parameter abs_c7(parameters["Abs{c7}"]);
            Parameter arg_c7(parameters["Arg{c7}"]);
            Parameter abs_c9(parameters["Abs{c9}"]);
            Parameter arg_c9(parameters["Arg{c9}"]);
            Parameter abs_c10(parameters["Abs{c10}"]);
            Parameter arg_c10(parameters["Arg{c10}"]);

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
                    "B->K^*ll::BR@LowRecoil,model=WilsonScan",
                    "B->K^*ll::Abar_FB@LowRecoil,model=WilsonScan",
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
                Kinematics kinematics;
                kinematics.declare("s_min");
                kinematics.set("s_min", 14.18);
                kinematics.declare("s_max");
                kinematics.set("s_max", 19.21);
                Options options;
                options.set("model", "WilsonScan");
                options.set("form-factors", "BZ2004");

                for (auto n = names.cbegin(), n_end = names.cend() ; n != n_end ; ++n)
                {
                    ObservablePtr observable = Observable::make(*n, parameters, kinematics, options);
                    WilsonPolynomial polynomial = make_polynomial(observable, std::list<std::string>{ "c7", "c9", "c10" });

                    for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i)
                    {
                        run_one(observable, polynomial, *i);
                    }
                }
            }

            // Test ratios
            {
                static const double eps = 1e-7;
                Kinematics kinematics;
                kinematics.declare("s_min");
                kinematics.set("s_min", 14.18);
                kinematics.declare("s_max");
                kinematics.set("s_max", 19.21);

                Parameters parameters = Parameters::Defaults();
                parameters["life_time::B_d"] = 1.530e-12;
                Parameter lambda = parameters["CKM::lambda"];
                Parameter A = parameters["CKM::A"];

                Options options;
                options.set("model", "WilsonScan");

                ObservablePtr numerator = Observable::make("B->K^*ll::Abar_FB@LowRecoil", parameters, kinematics, options);
                ObservablePtr denominator = Observable::make("B->K^*ll::BR@LowRecoil", parameters, kinematics, options);
                ObservablePtr observable = Observable::make("B->K^*ll::A_FB@LowRecoil", parameters, kinematics, options);

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

                std::list<std::string> coefficients{"c7", "c9", "c10"};

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
                "Abs{c7}",  "Arg{c7}",  "Abs{c7'}",  "Arg{c7'}",
                "Abs{c9}",  "Arg{c9}",  "Abs{c9'}",  "Arg{c9'}",
                "Abs{c10}", "Arg{c10}", "Abs{c10'}", "Arg{c10'}",
            };

            Parameters p = Parameters::Defaults();
            p["mass::mu"] = 0.0;
            Options o;
            o.set("model", "WilsonScan");
            o.set("form-factors", "BZ2004");

            std::vector<Parameter> variations;

            for (auto n = variation_names.cbegin(), n_end = variation_names.cend() ; n != n_end ; ++n)
            {
                variations.push_back(p[*n]);
            }

            Kinematics k;
            k.declare("s_min"); k.set("s_min", 14.18);
            k.declare("s_max"); k.set("s_max", 19.21);

            std::vector<ObservablePtr> observables;
            observables.push_back(Observable::make("B->K^*ll::BR@LowRecoil,q=d,l=mu",   p, k, o));
            observables.push_back(Observable::make("B->K^*ll::A_FB@LowRecoil,q=d,l=mu", p, k, o));
            observables.push_back(Observable::make("B->K^*ll::F_L@LowRecoil,q=d,l=mu",  p, k, o));

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil_TEST-btokstarll.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->K^*ll@LowRecoil --" << std::endl;
                RandomNumberGenerator rng;
                std::fstream file(filename.c_str(), std::fstream::out);
                file.precision(17);

                for (int i = 0 ; i < 1000 ; ++i)
                {
                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        *v = v->min() + (v->max() - v->min()) * rng();

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
                std::cout << "-- Verifying test case data for B->K^*ll@LowRecoil --" << std::endl;
                std::fstream file(filename.c_str(), std::fstream::in);

                std::string line;
                while (file)
                {
                    std::getline(file, line);
                    if (line.empty())
                        break;

                    std::stringstream ss(line);

                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        double value;
                        ss >> value;
                        *v = value;
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        double reference;
                        ss >> reference;

                        TEST_CHECK_RELATIVE_ERROR(reference, (*o)->evaluate(), 1e-3);
                    }
                }
            }
#endif
        }
} b_to_kstar_dilepton_low_recoil_bobeth_compatibility_test;

class BToKDileptonLowRecoilTest :
    public TestCase
{
    public:
        BToKDileptonLowRecoilTest() :
            TestCase("b_to_k_dilepton_low_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil (SM) */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.530e-12;
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.331;
                p["c8"] = -0.181;
                p["Abs{c9}"] = 4.27;
                p["Abs{c10}"] = 4.17;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K_d"] = 0.49761;
                // B mass
                p["mass::B_d"] = 5.27953;
                p["mass::b(MSbar)"] = 4.2;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004v2");

                BToKDilepton<LowRecoil> d(p, oo);

                /* q^2 = [14.18, 22.8] */
                {
                    const double eps = 1e-5;
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(14.18  ), 2.498607492e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(15.2575), 2.454634936e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(16.335 ), 2.374070832e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(17.4125), 2.244476570e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(18.49  ), 2.046333084e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(19.5675), 1.749827352e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(20.645 ), 1.312607210e-08, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(21.7225), 6.929281810e-09, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(22.8   ), 1.579971652e-10, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.differential_flat_term(15.0), 0.005562348378, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.differential_flat_term(22.0), 0.008213626852, eps);

                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(14.18, 22.8),       1.5267386e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(14.18, 22.8),             5.4236817e-03, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(14.18, 22.8), 1.0015589,     eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(14.18, 22.8),          2.2706273e-05, eps);
                }
            }

            /* Benchmark Point */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.530e-12;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // B mass
                p["mass::B_d"] = 5.27953;
                // Kaon mass
                p["mass::K_d"] = 0.49761;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.331;
                p["c8"] = -0.181;
                p["Abs{c9}"] = 4.27;
                p["Abs{c10}"] = 4.17;
                p["Arg{c7}"] = -M_PI / 2.0;
                p["Arg{c9}"] = +M_PI / 2.0;
                p["Arg{c10}"] = -M_PI / 2.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004v2");

                BToKDilepton<LowRecoil> d(p, oo);

                /* q^2 = [14.18, 22.8] */
                {
                    const double eps = 1e-5;

                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(14.18, 22.8),              1.5520940e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(14.18, 22.8),  1.4629637e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(14.18, 22.8),                    5.3935506e-03, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(14.18, 22.8),        1.0015315,     eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(14.18, 22.8),                 0.0609245,     eps);
                }
            }
        }
} b_to_k_dilepton_low_recoil_test;

class BToKDileptonLowRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKDileptonLowRecoilBobethCompatibilityTest() :
            TestCase("b_to_k_dilepton_low_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            static const std::vector<std::string> variation_names
            {
                "Abs{c7}",  "Arg{c7}",  "Abs{c7'}",  "Arg{c7'}",
                "Abs{c9}",  "Arg{c9}",  "Abs{c9'}",  "Arg{c9'}",
                "Abs{c10}", "Arg{c10}", "Abs{c10'}", "Arg{c10'}",
            };

            Parameters p = Parameters::Defaults();
            // old test data generated for K^+ mass set to K0 mass
            p["mass::K_u"] = 0.497614;
            Options o;
            o.set("model", "WilsonScan");
            o.set("form-factors", "KMPW2010");

            std::vector<Parameter> variations;

            for (auto n = variation_names.cbegin(), n_end = variation_names.cend() ; n != n_end ; ++n)
            {
                variations.push_back(p[*n]);
            }

            Kinematics k;
            k.declare("s_min"); k.set("s_min", 14.18);
            k.declare("s_max"); k.set("s_max", 22.86);

            std::vector<ObservablePtr> observables;
            std::vector<std::string> observable_names = {
                    "B->Kll::BR@LowRecoil,q=u,l=mu",
                    "B->Kll::F_H@LowRecoil,q=u,l=mu",
            };
            for (auto & s : observable_names)
            {
                observables.push_back(Observable::make(s, p, k, o));
                TEST_CHECK_MSG(observables.back(), "Could not create '" + s + "'");
            }

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil_TEST-btokll.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->Kll@LowRecoil --" << std::endl;
                RandomNumberGenerator rng;
                std::fstream file(filename.c_str(), std::fstream::out);
                file.precision(17);

                for (int i = 0 ; i < 1000 ; ++i)
                {
                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        *v = v->min() + (v->max() - v->min()) * rng();

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
                std::cout << "-- Verifying test case data for B->Kll@LowRecoil --" << std::endl;
                std::fstream file(filename.c_str(), std::fstream::in);
                TEST_CHECK_MSG(file, "'" + filename + "' does not exist");
                std::string line;
                while (file)
                {
                    std::getline(file, line);
                    if (line.empty())
                        break;

                    std::stringstream ss(line);

                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        double value;
                        ss >> value;
                        *v = value;
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        double reference;
                        ss >> reference;

                        TEST_CHECK_RELATIVE_ERROR(reference, (*o)->evaluate(), 1e-3);
                    }
                }
            }
#endif
        }
} b_to_k_dilepton_low_recoil_bobeth_compatibility_test;
