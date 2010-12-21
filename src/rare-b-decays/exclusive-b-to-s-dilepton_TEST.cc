/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/factory.hh>
#include <src/utils/complex.hh>
#include <src/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <iostream>

using namespace test;
using namespace eos;

class BToKstarDileptonLargeRecoilTest :
    public TestCase
{
    public:
        BToKstarDileptonLargeRecoilTest() :
            TestCase("b_to_kstar_dilepton_large_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Large Recoil */

            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;

                ObservableOptions oo;
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LargeRecoil> d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(2.00, 4.30) * 1e7,      +1.1122, eps); // 1.1125 for m_b_PS = 4.6
                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), +0.0857, eps); // 0.0911 for m_b_PS = 4.6
                TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(2.00, 4.30),  +0.7827, eps); // 0.7800 for m_b_PS = 4.6

                double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(2.00, 4.30) / d.integrated_branching_ratio(2.00, 4.30);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), a_fb,    eps);

                TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(1.00, 6.00) * 1e7,      +2.6105, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), +0.0471, eps);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(1.00, 6.00),  +0.7310, eps);

                a_fb = d.integrated_unnormalized_forward_backward_asymmetry(1.00, 6.00) / d.integrated_branching_ratio(1.00, 6.00);
                TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), a_fb,    eps);
            }
        }
} b_to_kstar_dilepton_large_recoil_test;

class BToKstarDileptonLargeRecoilPolynomialTest :
    public TestCase
{
    public:
        BToKstarDileptonLargeRecoilPolynomialTest() :
            TestCase("b_to_kstar_dilepton_large_recoil_polynomial_test")
        {
        }

        void run_one(const ObservablePtr & o, const Kinematics & k, const WilsonPolynomial & p, const std::array<double, 6> & values) const
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

            static const double eps = 2e-10;
            WilsonPolynomialEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(k), p.accept_returning<double>(evaluator), eps);
        }

        virtual void run() const
        {
            static const std::vector<std::string> names
            {
                "B->K^*ll::BR@LargeRecoil",
                "B->K^*ll::Abar_FB@LargeRecoil",
                "B->K^*ll::Fbar_L@LargeRecoil",
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
            kinematics.set("s_min", 1.0);
            kinematics.declare("s_max");
            kinematics.set("s_max", 6.0);
            ObservableOptions options;
            options.set("form-factors", "BZ2004");

            for (auto n = names.cbegin(), n_end = names.cend() ; n != n_end ; ++n)
            {
                ObservablePtr observable = RareBFactory::make(*n, parameters, options);
                TEST_CHECK(0 != observable);

                WilsonPolynomial polynomial = make_polynomial(observable, kinematics, std::list<std::string>{ "c7", "c9", "c10" });

                for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i)
                {
                    run_one(observable, kinematics, polynomial, *i);
                }
            }
        }
} b_to_kstar_dilepton_large_recoil_polynomial_test;

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
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;

                ObservableOptions oo;
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LowRecoil> d(p, oo);

                /* q^2 = [14.00, 19.21] */
                {
                    const double eps = 1e-4;

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(14.00, 19.21), -0.4093, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(14.00, 19.21),  +0.3497, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(14.00, 19.21),     -0.4835, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(14.00, 19.21),     +1.6892, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(14.00, 19.21),     -0.5758, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(14.00, 19.21),                        +0.9967, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(14.00, 19.21),                        -0.9727, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(14.00, 19.21),                        -0.9586, eps);

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
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(16.00, 19.21),     -0.486256, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(16.00, 19.21),                        +0.998622, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(16.00, 19.21),                        -0.970214, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(16.00, 19.21),                        -0.959887, eps);

                    double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(16.00, 19.21) / d.integrated_branching_ratio(16.00, 19.21);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(16.00, 19.21), a_fb,    eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  16.00)), +3.106224560e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  16.00)), -1.005941560e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 16.00)), +3.106224560e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 16.00)), +8.197948310e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  16.00)), -2.261711490e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  16.00)), +7.324485210e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 16.00)), -2.261711490e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 16.00)), -5.969109290e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   16.00)), +3.590664490e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   16.00)), -1.162825983e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  16.00)), +3.590664490e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  16.00)), +9.476482210e-12, eps);
                }
            }

            /* Low Recoil (Benchmark Point) */
            {
                Parameters p = Parameters::Defaults();
                p["Abs{c7}"] = 0.3;
                p["Arg{c7}"] = -M_PI / 2.0;
                p["c8"] = -0.181;
                p["Abs{c9}"] = 4.2;
                p["Arg{c9}"] = +M_PI / 2.0;
                p["Abs{c10}"] = 4.2;
                p["Arg{c10}"] = -M_PI / 2.0;

                ObservableOptions oo;
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LowRecoil> d(p, oo);

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  16.00)), +1.017630920e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  16.00)), -2.496817743e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 16.00)), -7.732915218e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 16.00)), -2.496817743e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  16.00)), -7.409598052e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  16.00)), +1.817988774e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 16.00)), +5.630508311e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 16.00)), +1.817988774e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   16.00)), +1.176338390e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   16.00)), -2.886215922e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  16.00)), -8.938923592e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  16.00)), -2.886215922e-12, eps);
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

        void run_one(const ObservablePtr & o, const Kinematics & k, const WilsonPolynomial & p, const std::array<double, 6> & values) const
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
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(k), p.accept_returning<double>(evaluator), eps);
        }

        virtual void run() const
        {
            static const std::vector<std::string> names
            {
                "B->K^*ll::BR@LowRecoil",
                "B->K^*ll::Abar_FB@LowRecoil",
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
            ObservableOptions options;
            options.set("form-factors", "BZ2004");

            for (auto n = names.cbegin(), n_end = names.cend() ; n != n_end ; ++n)
            {
                ObservablePtr observable = RareBFactory::make(*n, parameters, options);
                WilsonPolynomial polynomial = make_polynomial(observable, kinematics, std::list<std::string>{ "c7", "c9", "c10" });

                for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i)
                {
                    run_one(observable, kinematics, polynomial, *i);
                }
            }
        }
} b_to_kstar_dilepton_low_recoil_polynomial_test;
