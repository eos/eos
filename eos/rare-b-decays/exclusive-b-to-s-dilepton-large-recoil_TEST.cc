/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
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

#include <iomanip>

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
                p["life_time::B_d"] = 1.530e-12;
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"] = +4.27;
                p["Arg{c9}"] = 0.0;
                p["Abs{c10}"] = +4.173;
                p["Arg{c10}"] = M_PI;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // BHvD2010 parameters
                p["decay-constant::B_d"] = 0.200;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.2795;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["mass::W"] = 80.398;
                p["mass::tau"] = 1.77684;
                p["mass::mu"] = 0.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LargeRecoil> d(p, oo);

                const double eps = 1e-4;

                /* integrated observables */
                {
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(2.00, 4.30) * 1e7,      +1.0646, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), +0.0740, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(2.00, 4.30),  +0.7892, eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(1.00, 6.00) * 1e7,      +2.4788, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), +0.0394, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(1.00, 6.00),  +0.7387, eps);
                }

                /* inverse observables */
                {
                    TEST_CHECK_NEARLY_EQUAL(d.a_fb_zero_crossing(), +3.9370, eps);
                }

                /* transversity amplitudes at q^2 = 6.00 GeV^2 */
                {
                    static const double eps = 1e-17; // 1e-7..1e-4 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  6.00)),-1.26337489744e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  6.00)),-8.74365094218e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 6.00)),+7.65932494892e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 6.00)),-8.74365094218e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  6.00)),+5.38369467677e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  6.00)),-2.69948269691e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 6.00)),-2.44576143945e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 6.00)),-2.69948269691e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   6.00)),-5.79174565150e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   6.00)),+2.78267668675e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  6.00)),+2.51835429593e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  6.00)),+2.78267668675e-12, eps);
                }
            }

            // Benchmark Point (CPV)
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
                // b quark mass
                p["mass::b(MSbar)"] = 4.20;
                p["mass::W"] = 80.398;
                p["mass::mu"] = 0.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LargeRecoil> d(p, oo);

                /* observables */
                {
                    static const double eps = 1e-4;
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(1.0, 6.0),                        2.42908e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(1.0, 6.0),            2.38139e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(1.0, 6.0),             +0.0160992,  eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry_cp_averaged(1.0, 6.0), -0.0114422,  eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation(1.0, 6.0),              +0.755664,   eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation_cp_averaged(1.0, 6.0),  +0.757228,   eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2(1.0, 6.0),                 -0.0373159,  eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2_cp_averaged(1.0, 6.0),     -0.0408299,  eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_3(1.0, 6.0),                 +0.591474,   eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_4(1.0, 6.0),                 +0.75243,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_5(1.0, 6.0),                 +0.497716,   eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_re(1.0, 6.0),                +0.0878531,  eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_im(1.0, 6.0),                -0.00206701, eps);
                }

                /* transversity amplitudes at q^2 = 6.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  6.00)),-1.75667412993e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  6.00)),-1.25792842651e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 6.00)),-1.75667412993e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 6.00)),+9.07262247124e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  6.00)),-7.21683758673e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  6.00)),+5.39530540830e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 6.00)),-7.21683758673e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 6.00)),-2.48488890209e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   6.00)),+7.43924297713e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   6.00)),-5.80247546021e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  6.00)),+7.43924297713e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  6.00)),+2.56146713733e-11, eps);
                }
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

            static const double eps = 2e-10;
            WilsonPolynomialEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(), p.accept_returning<double>(evaluator), eps);
        }

        virtual void run() const
        {
            static const std::vector<std::string> names
            {
                "B->K^*ll::BR@LargeRecoil,model=WilsonScan",
                "B->K^*ll::J_6s@LargeRecoil,model=WilsonScan",
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
            // Kaon mass
            parameters["mass::K_d"] = 0.49761;
            // b quark mass
            parameters["mass::b(MSbar)"] = 4.2;
            Kinematics kinematics;
            kinematics.declare("s_min");
            kinematics.set("s_min", 1.0);
            kinematics.declare("s_max");
            kinematics.set("s_max", 6.0);
            Options options;
            options.set("form-factors", "BZ2004");

            for (auto n = names.cbegin(), n_end = names.cend() ; n != n_end ; ++n)
            {
                ObservablePtr observable = Observable::make(*n, parameters, kinematics, options);
                TEST_CHECK(ObservablePtr() != observable);

                WilsonPolynomial polynomial = make_polynomial(observable, std::list<std::string>{ "c7", "c9", "c10" });

                for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i)
                {
                    run_one(observable, polynomial, *i);
                }
            }
        }
} b_to_kstar_dilepton_large_recoil_polynomial_test;

class BToKDileptonLargeRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKDileptonLargeRecoilBobethCompatibilityTest() :
            TestCase("b_to_k_dilepton_large_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            // Christoph uses \Delta C instead of C for C9, C10
            Parameters p = Parameters::Defaults();
            p["c1"] = -0.321961806;
            p["c2"] = 1.009246233;
            p["c3"] = -0.005191551838;
            p["c4"] = -0.08787376815;
            p["c5"] = 0.0003573490117;
            p["c6"] = 0.001011323348;
            p["Re{c7}"] = -0.3367243768 + 0.1;
            p["Im{c7}"] = 0.2;
            p["Re{c7'}"] = 0.3;
            p["Im{c7'}"] = 0.4;
            p["c8"] = -0.1825051462;
            p["Re{c9}"] = 4.295481065 + 1;
            p["Im{c9}"] = 0.5;
            p["Re{c9'}"] = 2;
            p["Im{c9'}"] = 1.5;
            p["Re{c10}"] = -4.196943811 + 3;
            p["Im{c10}"] = 2.5;
            p["Re{c10'}"] = 4;
            p["Im{c10'}"] = 3.5;
            p["Re{cS}"] = 0.5;
            p["Im{cS}"] = 1;
            p["Re{cS'}"] = 0.6;
            p["Im{cS'}"] = 1.1;
            p["Re{cP}"] = 0.7;
            p["Im{cP}"] = 1.2;
            p["Re{cP'}"] = 0.8;
            p["Im{cP'}"] = 1.3;
            p["Re{cT}"] = 0.9;
            p["Im{cT}"] = 1.4;
            p["Re{cT5}"] = 1.0;
            p["Im{cT5}"] = 1.5;

            Options oo;
            oo.set("model", "WilsonScan");
            oo.set("scan-mode", "cartesian");
            oo.set("form-factors", "KMPW2010");
            oo.set("l", "mu");
            oo.set("q", "u");

            BToKDilepton<LargeRecoil> d(p, oo);
            static const double eps = 1e-3;
            static const double s = 6.0;

            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_A(s)), 2.803056189, 1e-14);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_A(s)), 6, 1e-14);
            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_S(s)), 3.277235546, eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_S(s)), 6.256540588, eps);
            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_T(s)), 7.695315895, eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_T(s)), 11.97049139, eps);
            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_T5(s)), 8.550350995, eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_T5(s)), 12.82552649, eps);
            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_P(s)), 4.010598621, eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_P(s)), 6.467135768 , eps);

            /* difference comes from cal_T, F_V affects everything below */
            TEST_CHECK_RELATIVE_ERROR(std::real(d.F_V(s)), 7.787757339, 7 * eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(d.F_V(s)), 3.226118406, 7 * eps);

            TEST_CHECK_RELATIVE_ERROR(d.a_l(s),  3.935609789e-20, 8 * eps);
            TEST_CHECK_RELATIVE_ERROR(d.b_l(s),  9.695777426e-21, 2 * eps);
            TEST_CHECK_RELATIVE_ERROR(d.c_l(s), -2.771498043e-20, 10 * eps);

            const double tau_over_hbar = p["life_time::B_u"] / p["hbar"];
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(1, 6),
                                      2.910731397e-19 * tau_over_hbar, 7 * eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(1, 6),
                                      2.887189571e-19 * tau_over_hbar, 5 * eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(1, 6), 0.1093626586 , 6 * eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(1, 6),  0.2778727304, 5 * eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(1, 6), 1.07282598, eps);
            // todo debug cal_T to find out why difference is 100% or more
//            TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(1, 6), -0.003672979605, eps);
        }
} b_to_k_dilepton_large_recoil_bobeth_compatibility_test;
