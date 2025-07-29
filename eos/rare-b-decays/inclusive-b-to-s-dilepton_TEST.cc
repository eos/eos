/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
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
#include <eos/utils/expression-evaluator.hh>
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
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["sb::mu"] = 4.2;
                p["b->s::Re{c7}"] = -0.331;
                p["b->s::c8"] = -0.181;
                p["sbmumu::mu"] = 4.2;
                p["b->smumu::Re{c9}"]  = +4.27;
                p["b->smumu::Re{c10}"] = -4.173;
                p["CKM::abs(V_ub)"]    =  0.003631275231633653;
                p["CKM::arg(V_ub)"]    = -1.210765774253535;
                p["CKM::abs(V_cb)"]    =  0.041996951916414726;
                p["CKM::arg(V_cb)"]    =  0.0;
                p["CKM::abs(V_tb)"]    =  0.9991111344469873;
                p["CKM::arg(V_tb)"]    =  0.0;
                p["CKM::abs(V_us)"]    =  0.22534851424944366;
                p["CKM::arg(V_us)"]    =  0.0;
                p["CKM::abs(V_cs)"]    =  0.9734061815416853;
                p["CKM::arg(V_cs)"]    = -3.304199362533668e-05;
                p["CKM::abs(V_ts)"]    =  0.04121212396309175;
                p["CKM::arg(V_ts)"]    = -3.1230250224697222;

                // quark masses
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.27;

                Options oo;
                oo.declare("model"_ok, "WET");
                oo.declare("l"_ok, "mu");

                BToXsDilepton<HLMW2005> d_mu(p, oo);

                oo.declare("l"_ok, "e");
                BToXsDilepton<HLMW2005> d_e(p, oo);

                {
                    Diagnostics diagnostics = d_mu.diagnostics();
                    static const std::vector<std::pair<double, double>> reference
                    {
                        /* phi_ll */
                        std::make_pair(+1.01380e-05, 1e-9), // phi_ll(s = 1.0GeV^2)
                        std::make_pair(+5.90542e-06, 1e-9), // phi_ll(s = 6.0GeV^2)
                    };

                    TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
                }

                static const double eps = 1e-11;

                TEST_CHECK_NEARLY_EQUAL(d_mu.integrated_branching_ratio(1.00, 6.00), 1.40781e-06, eps);
                TEST_CHECK_NEARLY_EQUAL(d_e.integrated_branching_ratio(1.00, 6.00),  1.46487e-06, eps);
            }

            // Benchmark Point (C_7,9,10 = 0, C_7',9',10' = C_7,9,10^SM)
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["sb::mu"] = 4.2;
                p["b->s::Re{c7}"]  = 0.0;
                p["b->s::Re{c7'}"] = -0.331;
                p["b->s::c8"] = -0.181;
                p["sbmumu::mu"] = 4.2;
                p["b->smumu::Re{c9}"]   = 0.0;
                p["b->smumu::Re{c9'}"]  = +4.27;
                p["b->smumu::Re{c10}"]  = 0.0;
                p["b->smumu::Re{c10'}"] = -4.173;
                p["CKM::abs(V_ub)"]     =  0.003631275231633653;
                p["CKM::arg(V_ub)"]     = -1.210765774253535;
                p["CKM::abs(V_cb)"]     =  0.041996951916414726;
                p["CKM::arg(V_cb)"]     =  0.0;
                p["CKM::abs(V_tb)"]     =  0.9991111344469873;
                p["CKM::arg(V_tb)"]     =  0.0;
                p["CKM::abs(V_us)"]     =  0.22534851424944366;
                p["CKM::arg(V_us)"]     =  0.0;
                p["CKM::abs(V_cs)"]     =  0.9734061815416853;
                p["CKM::arg(V_cs)"]     = -3.304199362533668e-05;
                p["CKM::abs(V_ts)"]     =  0.04121212396309175;
                p["CKM::arg(V_ts)"]     = -3.1230250224697222;

                // quark masses
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.27;

                Options oo;
                oo.declare("model"_ok, "WET");
                oo.declare("l"_ok, "mu");

                BToXsDilepton<HLMW2005> d_mu(p, oo);

                oo.declare("l"_ok, "e");
                BToXsDilepton<HLMW2005> d_e(p, oo);

                const double eps = 1e-11;

                {
                    Diagnostics diagnostics = d_mu.diagnostics();
                    static const std::vector<std::pair<double, double>> reference
                    {
                        /* phi_ll */
                        std::make_pair(+0.93156e-05, 1e-9), // phi_ll(s = 1.0GeV^2)
                        std::make_pair(+5.71634e-06, 1e-9), // phi_ll(s = 6.0GeV^2)
                    };

                    TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
                }

                TEST_CHECK_NEARLY_EQUAL(d_mu.integrated_branching_ratio(1.00, 6.00), 1.35152e-06, eps);
                TEST_CHECK_NEARLY_EQUAL(d_e.integrated_branching_ratio(1.00, 6.00),  1.39879e-06, eps);
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

        void run_one(const ObservablePtr & o, const exp::Expression & p, const std::array<double, 6> & values) const
        {
            Parameters parameters = o->parameters();
            Parameter re_c7(parameters["b->s::Re{c7}"]);
            Parameter im_c7(parameters["b->s::Im{c7}"]);
            Parameter re_c9(parameters["b->smumu::Re{c9}"]);
            Parameter im_c9(parameters["b->smumu::Im{c9}"]);
            Parameter re_c10(parameters["b->smumu::Re{c10}"]);
            Parameter im_c10(parameters["b->smumu::Im{c10}"]);

            re_c7 = values[0];
            im_c7 = values[1];
            re_c9 = values[2];
            im_c9 = values[3];
            re_c10 = values[4];
            im_c10 = values[5];

            static const double eps = 5e-8;
            exp::ExpressionEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(), std::visit(evaluator, p), eps);
        }

        virtual void run() const
        {
            static const std::vector<std::string> names
            {
                "B->X_sll::BR@HLMW2005;model=WET",
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
            parameters["CKM::abs(V_ub)"]     =  0.003631275231633653;
            parameters["CKM::arg(V_ub)"]     = -1.210765774253535;
            parameters["CKM::abs(V_cb)"]     =  0.041996951916414726;
            parameters["CKM::arg(V_cb)"]     =  0.0;
            parameters["CKM::abs(V_tb)"]     =  0.9991111344469873;
            parameters["CKM::arg(V_tb)"]     =  0.0;
            parameters["CKM::abs(V_us)"]     =  0.22534851424944366;
            parameters["CKM::arg(V_us)"]     =  0.0;
            parameters["CKM::abs(V_cs)"]     =  0.9734061815416853;
            parameters["CKM::arg(V_cs)"]     = -3.304199362533668e-05;
            parameters["CKM::abs(V_ts)"]     =  0.04121212396309175;
            parameters["CKM::arg(V_ts)"]     = -3.1230250224697222;
            Kinematics kinematics
            {
                { "q2_min", 1.0 },
                { "q2_max", 6.0 }
            };

            for (const auto & name : names)
            {
                ObservablePtr observable = Observable::make(name, parameters, kinematics, Options());
                TEST_CHECK(ObservablePtr() != observable);

                auto polynomial = make_polynomial(observable, std::list<std::string>{ "b->s::Re{c7}", "b->s::Im{c7}", "b->smumu::Re{c9}", "b->smumu::Im{c9}", "b->smumu::Re{c10}", "b->smumu::Im{c10}" });

                for (const auto & input : inputs)
                {
                    run_one(observable, polynomial, input);
                }
            }
        }
} b_to_x_s_dilepton_large_recoil_polynomial_test;
