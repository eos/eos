/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BToXsDileptonLargeRecoilTest :
    public TestCase
{
    public:
        BToXsDileptonLargeRecoilTest() :
            TestCase("b_to_x_s_dilepton_HLMW2005_test")
        {
        }

        virtual void run() const
        {
            /* HLMW2005 */

            // Standard Model
            {
                Parameters p = Parameters::Defaults();
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

                // quark masses
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.27;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("l", "mu");

                BToXsDilepton<HLMW2005> d_mu(p, oo);

                oo.set("l", "e");
                BToXsDilepton<HLMW2005> d_e(p, oo);

                const double eps = 1e-11;

                {
                    Diagnostics diagnostics = d_mu.diagnostics();
                    static const std::vector<std::pair<double, double>> reference
                    {
                        /* phi_ll */
                        std::make_pair(+1.01995e-05, 1e-9), // phi_ll(s = 1.0GeV^2)
                        std::make_pair(+6.91719e-06, 1e-9), // phi_ll(s = 6.0GeV^2)
                    };

                    TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
                }

                TEST_CHECK_NEARLY_EQUAL(1.47462e-06, d_mu.integrated_branching_ratio(1.00, 6.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.53382e-06, d_e.integrated_branching_ratio(1.00, 6.00), eps);
            }

            // Benchmark Point (C_7,9,10 = 0, C_7',9',10' = C_7,9,10^SM)
            {
                Parameters p = Parameters::Defaults();
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"]  = 0.0;
                p["Abs{c7'}"] = 0.331;
                p["Arg{c7'}"] = M_PI;
                p["c8"] = -0.181;
                p["Abs{c9}"]  = 0.0;
                p["Abs{c9'}"] = +4.27;
                p["Arg{c9'}"] = 0.0;
                p["Abs{c10}"]  = 0.0;
                p["Abs{c10'}"] = +4.173;
                p["Arg{c10'}"] = M_PI;

                // quark masses
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.27;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("l", "mu");

                BToXsDilepton<HLMW2005> d_mu(p, oo);

                oo.set("l", "e");
                BToXsDilepton<HLMW2005> d_e(p, oo);

                const double eps = 1e-11;

                {
                    Diagnostics diagnostics = d_mu.diagnostics();
                    static const std::vector<std::pair<double, double>> reference
                    {
                        /* phi_ll */
                        std::make_pair(+0.93156e-05, 1e-9), // phi_ll(s = 1.0GeV^2)
                        std::make_pair(+5.84538e-06, 1e-9), // phi_ll(s = 6.0GeV^2)
                    };

                    TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
                }

                TEST_CHECK_NEARLY_EQUAL(1.35795e-06, d_mu.integrated_branching_ratio(1.00, 6.00), eps);
                TEST_CHECK_NEARLY_EQUAL(1.40541e-06, d_e.integrated_branching_ratio(1.00, 6.00), eps);
            }
        }
} b_to_x_s_dilepton_large_recoil_test;

class BToXsDileptonLargeRecoilPolynomialTest :
    public TestCase
{
    public:
        BToXsDileptonLargeRecoilPolynomialTest() :
            TestCase("b_to_x_s_dilepton_large_recoil_polynomial_test")
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

            static const double eps = 5e-8;
            WilsonPolynomialEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(), p.accept_returning<double>(evaluator), eps);
        }

        virtual void run() const
        {
            static const std::vector<std::string> names
            {
                "B->X_sll::BR@HLMW2005,model=WilsonScan",
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

            for (auto n = names.cbegin(), n_end = names.cend() ; n != n_end ; ++n)
            {
                ObservablePtr observable = Observable::make(*n, parameters, kinematics, Options());
                TEST_CHECK(ObservablePtr() != observable);

                WilsonPolynomial polynomial = make_polynomial(observable, std::list<std::string>{ "c7", "c9", "c10" });

                for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i)
                {
                    run_one(observable, polynomial, *i);
                }
            }
        }
} b_to_x_s_dilepton_large_recoil_polynomial_test;
