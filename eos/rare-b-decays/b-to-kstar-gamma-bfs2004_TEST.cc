/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2025 Danny van Dyk
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
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
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

class BToKstarGammaTest :
    public TestCase
{
    public:
        BToKstarGammaTest() :
            TestCase("b_to_kstar_gamma_test")
        {
        }

        virtual void run() const
        {
            /* QCDF */

            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = -0.331;
                p["b->s::Re{c7'}"] = -0.00659; // m_s(m_b) / m_b(m_b) * Abs{c7} = 85 / 4200 * Abs{c7}
                p["b->s::c8"] = -0.181;
                // PDG 2010 CKM parameters
                p["CKM::A"]         =  0.812;
                p["CKM::lambda"]    =  0.22543;
                p["CKM::rhobar"]    =  0.144;
                p["CKM::etabar"]    =  0.342;
                p["CKM::abs(V_ub)"] =  0.003540950873054711;
                p["CKM::arg(V_ub)"] = -1.1728563751359748;
                p["CKM::abs(V_cb)"] =  0.04126451344307112;
                p["CKM::arg(V_cb)"] =  0.0;
                p["CKM::abs(V_tb)"] =  0.9991419776905534;
                p["CKM::arg(V_tb)"] =  0.0;
                p["CKM::abs(V_td)"] =  0.008576901910577167;
                p["CKM::arg(V_td)"] = -0.37951557931964897;
                p["CKM::abs(V_us)"] =  0.22542858674178629;
                p["CKM::arg(V_us)"] =  0.0;
                p["CKM::abs(V_cs)"] =  0.9734167680132911;
                p["CKM::arg(V_cs)"] = -3.119448393424795e-05;
                p["CKM::abs(V_ts)"] =  0.04051834255894421;
                p["CKM::arg(V_ts)"] = -3.123445879630718;
                p["decay-constant::B_d"] = 0.200;
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.27;
                p["mass::B_d"] = 5.27958;
                p["mass::K_d^*"] = 0.89594;
                p["K^*::a_1_para@1GeV"] = 0.1;
                p["K^*::a_1_perp@1GeV"] = 0.1;
                p["K^*::a_2_para@1GeV"] = 0.1;
                p["K^*::a_2_perp@1GeV"] = 0.1;
                p["B::1/lambda_B_p"] = 1.0 / 0.485;

                Options oo
                {
                    {"model"_ok, "WET"},
                    {"tag"_ok, "BFS2004"},
                    {"form-factors"_ok, "KMPW2010"}
                };

                BToKstarGamma d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(),             +5.45306e-5, eps);

                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::S_K^*gamma", p, Kinematics(), oo)->evaluate(),  -3.94778e-2, eps);
                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::C_K^*gamma", p, Kinematics(), oo)->evaluate(),  3.66320e-3, eps);
                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::BR",         p, Kinematics(), oo)->evaluate(),  5.47311e-5, eps);
            }

            // Benchmark Point (CPV)
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = 0.0;
                p["b->s::Im{c7}"] = -0.331;
                p["b->s::Re{c7'}"] = 0.0;
                p["b->s::Im{c7'}"] = -0.00659; // m_s(m_b) / m_b(m_b) * Abs{c7} = 85 / 4200 * Abs{c7}
                p["b->s::c8"] = -0.181;
                // PDG 2010 CKM parameters
                p["CKM::A"]         =  0.812;
                p["CKM::lambda"]    =  0.22543;
                p["CKM::rhobar"]    =  0.144;
                p["CKM::etabar"]    =  0.342;
                p["CKM::abs(V_ub)"] =  0.003540950873054711;
                p["CKM::arg(V_ub)"] = -1.1728563751359748;
                p["CKM::abs(V_cb)"] =  0.04126451344307112;
                p["CKM::arg(V_cb)"] =  0.0;
                p["CKM::abs(V_tb)"] =  0.9991419776905534;
                p["CKM::arg(V_tb)"] =  0.0;
                p["CKM::abs(V_td)"] =  0.008576901910577167;
                p["CKM::arg(V_td)"] = -0.37951557931964897;
                p["CKM::abs(V_us)"] =  0.22542858674178629;
                p["CKM::arg(V_us)"] =  0.0;
                p["CKM::abs(V_cs)"] =  0.9734167680132911;
                p["CKM::arg(V_cs)"] = -3.119448393424795e-05;
                p["CKM::abs(V_ts)"] =  0.04051834255894421;
                p["CKM::arg(V_ts)"] = -3.123445879630718;
                p["decay-constant::B_d"] = 0.200;
                p["mass::b(MSbar)"] = 4.2;
                p["mass::B_d"] = 5.27958;
                p["mass::K_d^*"] = 0.89594;
                p["K^*::a_1_para@1GeV"] = 0.1;
                p["K^*::a_1_perp@1GeV"] = 0.1;
                p["K^*::a_2_para@1GeV"] = 0.1;
                p["K^*::a_2_perp@1GeV"] = 0.1;
                p["B::1/lambda_B_p"] = 1.0 / 0.485;

                Options oo
                {
                    {"model"_ok, "WET"},
                    {"tag"_ok, "BFS2004"},
                    {"form-factors"_ok, "KMPW2010"}
                };

                BToKstarGamma d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(),             +5.65584e-5, eps);

                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::S_K^*gamma", p, Kinematics(), oo)->evaluate(), +4.72504e-2, eps);
                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::C_K^*gamma", p, Kinematics(), oo)->evaluate(), -4.13944e-1, eps);
                TEST_CHECK_RELATIVE_ERROR(Observable::make("B->K^*gamma::BR",         p, Kinematics(), oo)->evaluate(),  4.00005e-5, eps);
            }
        }
} b_to_kstar_gamma_test;

class BToKstarGammaBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKstarGammaBobethCompatibilityTest() :
            TestCase("b_to_kstar_gamma_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            static const std::vector<std::string> variation_names
            {
                "b->s::Re{c7}",  "b->s::Im{c7}",  "b->s::Re{c7'}",  "b->s::Im{c7'}",
            };

            Parameters p = Parameters::Defaults();
            p["CKM::abs(V_ub)"] =  0.003631275231633653;
            p["CKM::arg(V_ub)"] = -1.210765774253535;
            p["CKM::abs(V_cb)"] =  0.041996951915501936;
            p["CKM::arg(V_cb)"] =  0.0;
            p["CKM::abs(V_tb)"] =  0.9991111398988599;
            p["CKM::arg(V_tb)"] =  0.0;
            p["CKM::abs(V_us)"] =  0.22534851424454577;
            p["CKM::arg(V_us)"] =  0.0;
            p["CKM::abs(V_cs)"] =  0.9734061890640963;
            p["CKM::arg(V_cs)"] = -0.0000330419933699906;
            p["CKM::abs(V_ts)"] =  0.04121211253368258;
            p["CKM::arg(V_ts)"] = -3.1230250245535283;
            p["CKM::abs(V_td)"] =  0.008859566045351227;
            p["CKM::arg(V_td)"] = -0.38266;
            p["decay-constant::B_d"] = 0.1906;
            p["mass::B_d"] = 5.27958;
            p["mass::K_d^*"] = 0.89594;
            p["K^*::a_1_para@1GeV"] = 0.1;
            p["K^*::a_1_perp@1GeV"] = 0.1;
            p["K^*::a_2_para@1GeV"] = 0.1;
            p["K^*::a_2_perp@1GeV"] = 0.1;
            p["B::1/lambda_B_p"] = 1.0 / 0.485;

            Options o
            {
                {"model"_ok, "WET"},
                {"tag"_ok, "BFS2004"},
                {"form-factors"_ok, "KMPW2010"}
            };

            std::vector<Parameter> variations;

            for (const auto & variation_name : variation_names)
            {
                variations.push_back(p[variation_name]);
            }

            Kinematics k;

            std::vector<ObservablePtr> observables;
            observables.push_back(Observable::make("B->K^*gamma::BR_CP_specific;q=d",  p, k, o));
            observables.push_back(Observable::make("B->K^*gamma::S_K^*gamma;q=d", p, k, o));
            observables.push_back(Observable::make("B->K^*gamma::C_K^*gamma;q=d", p, k, o));
            observables.push_back(Observable::make("B->K^*gamma::A_I", p, k, o));

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-gamma_TEST-btokstargamma.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->K^*gamma --" << std::endl;
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
                std::cout << "-- Verifying test case data for B->K^*gamma --" << std::endl;
                std::fstream file(filename.c_str(), std::fstream::in);

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
} b_to_kstar_gamma_bobeth_compatibility_test;
