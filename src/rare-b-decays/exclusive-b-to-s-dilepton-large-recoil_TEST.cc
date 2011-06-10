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
#include <src/observable.hh>
#include <src/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
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
                p["life_time::B_d"] = 1.530e-12;
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

                    double a_fb = d.integrated_unnormalized_forward_backward_asymmetry(2.00, 4.30) / d.integrated_branching_ratio(2.00, 4.30);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(2.00, 4.30), a_fb,    eps);

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_branching_ratio(1.00, 6.00) * 1e7,      +2.50997, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), +0.0471, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(1.00, 6.00),  +0.7310, eps);

                    a_fb = d.integrated_unnormalized_forward_backward_asymmetry(1.00, 6.00) / d.integrated_branching_ratio(1.00, 6.00);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(1.00, 6.00), a_fb,    eps);
                }

                /* inverse observables */
                {
                    TEST_CHECK_NEARLY_EQUAL(d.a_fb_zero_crossing(), +4.03824, eps);
                }

                /* transversity amplitudes at q^2 = 6.00 GeV^2 */
                {
                    static const double eps = 1e-18; // 1e-6 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  6.00)), -1.28943854e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  6.00)), +1.8013703e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 6.00)), +5.0519606e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 6.00)), +1.8013703e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  6.00)), +5.4778926e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  6.00)), -2.7116335e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 6.00)), -2.3515051e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 6.00)), -2.7116335e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   6.00)), -5.8888449e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   6.00)), +2.7952019e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  6.00)), +2.4211931e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(right_handed,  6.00)), +2.7952019e-12, eps);
                }
            }

            // Benchmark Point (CPV)
            {
                Parameters p = Parameters::Defaults();
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
                // B mass
                p["mass::B_d"] = 5.27953;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton<LargeRecoil> d(p, oo);

                /* transversity amplitudes at q^2 = 6.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(left_handed,  6.00)), -4.314869773e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(left_handed,  6.00)), -1.230969071e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_long(right_handed, 6.00)), -4.314869773e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_long(right_handed, 6.00)), +1.176755204e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(left_handed,  6.00)), +1.818410935e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(left_handed,  6.00)), +5.371595755e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_perp(right_handed, 6.00)), +1.818410935e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_perp(right_handed, 6.00)), -2.508539778e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(left_handed,   6.00)), -1.874449940e-13, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(d.a_par(left_handed,   6.00)), -5.778033336e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(d.a_par(right_handed,  6.00)), -1.874449940e-13, eps);
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
                "B->K^*ll::Abar_FB@LargeRecoil,model=WilsonScan",
                "B->K^*ll::Fbar_L@LargeRecoil,model=WilsonScan",
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
                // B mass
                p["mass::B_d"] = 5.27953;
                // B decay constant
                p["decay-constant::B_d"] = 0.200;

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
        }
} b_to_k_dilepton_large_recoil_test;
