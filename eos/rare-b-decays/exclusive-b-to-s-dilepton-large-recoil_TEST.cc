/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <iostream>

// Uncomment the following #define to generate new test data for the Bobeth compatibility tests
//#define EOS_GENERATE_TEST_DATA

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
                p["mass::K^*0"] = 0.896;
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
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(2.00, 4.30) * 1e7,      +1.0694, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), +0.0857, eps); // 0.0911 for m_b_PS = 4.6
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(2.00, 4.30),  +0.7827, eps); // 0.7800 for m_b_PS = 4.6

                    //double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(2.00, 4.30) / d.integrated_branching_ratio(2.00, 4.30);
                    //TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), a_fb,    eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(1.00, 6.00) * 1e7,      +2.50997, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), +0.0471, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(1.00, 6.00),  +0.7310, eps);

                    //a_fb = d.integrated_unnormalized_forward_backward_asymmetry(1.00, 6.00) / d.integrated_branching_ratio(1.00, 6.00);
                    //TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), a_fb,    eps);
                }

                /* inverse observables */
                {
                    TEST_CHECK_NEARLY_EQUAL(d.a_fb_zero_crossing(), +4.03841, eps);
                }

                /* transversity amplitudes at q^2 = 6.00 GeV^2 */
                {
                    static const double eps = 1e-18; // 1e-6 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  6.00)), -1.28943624e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  6.00)), +1.80137030e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 6.00)), +5.05219138e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 6.00)), +1.80137030e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  6.00)), +5.47777162e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  6.00)), -2.71163350e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 6.00)), -2.35162609e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 6.00)), -2.71163350e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   6.00)), -5.88872010e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   6.00)), +2.79520190e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  6.00)), +2.42131786e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  6.00)), +2.79520190e-12, eps);
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
                p["mass::K^*0"] = 0.896;
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
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(1.0, 6.0),                         2.41382e-7, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(1.0, 6.0),             2.40579e-7, eps);
          //          TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(1.0, 6.0),             +2.39239e-2, eps);
          //          TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry_cp_averaged(1.0, 6.0), -9.16268e-3, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation(1.0, 6.0),              +0.74575,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_longitudinal_polarisation_cp_averaged(1.0, 6.0),  +0.75409,    eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2(1.0, 6.0),                 -3.57972e-2, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_transverse_asymmetry_2_cp_averaged(1.0, 6.0),     -4.01245e-2, eps);
                }

                /* transversity amplitudes at q^2 = 6.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  6.00)), -4.314625142e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  6.00)), -1.230969071e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 6.00)), -4.314625142e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 6.00)), +1.176755204e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  6.00)), +1.805581664e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  6.00)), +5.371595755e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 6.00)), +1.805581664e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 6.00)), -2.508539778e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   6.00)), -1.861225302e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   6.00)), -5.778033336e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  6.00)), -1.861225302e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  6.00)), +2.585846875e-11, eps);
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
    //            "B->K^*ll::Fbar_L@LargeRecoil,model=WilsonScan",
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
            parameters["mass::K0"] = 0.49761;
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

class BToKDileptonLargeRecoilTest :
    public TestCase
{
    public:
        BToKDileptonLargeRecoilTest() :
            TestCase("b_to_k_dilepton_large_recoil_test")
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
                // B decay constant
                p["decay-constant::B_d"] = 0.200;
                // Kaon mass
                p["mass::K0"] = 0.49761;
                // B mass
                p["mass::B_d"] = 5.27953;
                // b quark mass
                p["mass::b(MSbar)"] = 4.20;
                p["mass::W"] = 80.398;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004v2");
                oo.set("l", "mu");

                BToKDilepton<LargeRecoil> d(p, oo);

                const double eps = 1e-7;
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(1.0625), 3.145791429e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(2.125 ), 3.124835257e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(3.1875), 3.119065711e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(4.25  ), 3.127702599e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(5.3125), 3.163614767e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(6.375 ), 3.345862737e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(7.4375), 3.297253470e-08, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(8.5   ), 3.149745610e-08, eps);

                TEST_CHECK_RELATIVE_ERROR(d.differential_flat_term(1.0), 0.06529669784, eps);
                TEST_CHECK_RELATIVE_ERROR(d.differential_flat_term(8.0), 0.009094268877, eps);

                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(0.1, 8.68), 2.7329633e-7, eps);
                TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(0.1, 8.68), 3.325593091e-2, eps);

                // Standard bin 1.0 <= s <= 6.0 GeV^2
                TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(1.0, 6.0), 1.0003573, eps);
            }

            /* Benchmark Point */
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
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // B decay constant
                p["decay-constant::B_d"] = 0.200;
                // B mass
                p["mass::B_d"] = 5.27953;
                // Kaon mass
                p["mass::K0"] = 0.49761;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                p["Abs{c9}"] = 4.27;
                p["Abs{c10}"] = 4.17;
                p["Arg{c7}"] = -M_PI / 2.0;
                p["Arg{c9}"] = +M_PI / 2.0;
                p["Arg{c10}"] = -M_PI / 2.0;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004v2");

                BToKDilepton<LargeRecoil> d(p, oo);

                /* q^2 = [14.18, 22.8] */
                {
                    const double eps = 1e-5;

                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(1.0, 6.0),              1.5242329e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(1.0, 6.0),  1.5217576e-07, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(1.0, 6.0),                    2.4491886e-02, eps);
                    TEST_CHECK_RELATIVE_ERROR(d.integrated_ratio_muons_electrons(1.0, 6.0),        1.0003682,     eps);
                }
            }
        }
} b_to_k_dilepton_large_recoil_test;

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
            static const std::vector<std::string> variation_names
            {
                "Abs{c7}",  "Arg{c7}",  "Abs{c7'}",  "Arg{c7'}",
                "Abs{c9}",  "Arg{c9}",  "Abs{c9'}",  "Arg{c9'}",
                "Abs{c10}", "Arg{c10}", "Abs{c10'}", "Arg{c10'}",
            };

            Parameters p = Parameters::Defaults();
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
            observables.push_back(Observable::make("B->Kll::BR@LowRecoil,q=u,l=mu",  p, k, o));
            observables.push_back(Observable::make("B->Kll::F_H@LowRecoil,q=u,l=mu", p, k, o));

            std::string filename(EOS_BUILDDIR "/eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil_TEST-btokll.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->Kll@LargeRecoil --" << std::endl;
                RandomNumberGenerator rng;
                std::fstream file(filename.c_str(), std::fstream::out);
                file.precision(17);

                for (int i = 0 ; i < 1000 ; ++i)
                {
                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        *v = v->sample(rng);
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
                std::cout << "-- Verifying test case data for B->Kll@LargeRecoil --" << std::endl;
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
} b_to_k_dilepton_large_recoil_bobeth_compatibility_test;
