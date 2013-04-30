/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <eos/utils/complex.hh>
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
                p["c1"] = -0.32300000;
                p["c2"] = +1.00931000;
                p["c3"] = -0.00522869;
                p["c4"] = -0.08794730;
                p["c5"] = +0.00037476;
                p["c6"] = +0.00105859;
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = M_PI;
                p["Abs{c7'}"] = 0.00659; // m_s(m_b) / m_b(m_b) * Abs{c7} = 85 / 4200 * Abs{c7}
                p["Arg{c7'}"] = M_PI;
                p["c8"] = -0.181;
                // PDG 2010 CKM parameters
                p["CKM::A"] = 0.812;
                p["CKM::lambda"] = 0.22543;
                p["CKM::rhobar"] = 0.144;
                p["CKM::etabar"] = 0.342;
                p["decay-constant::B_d"] = 0.200;
                p["mass::b(MSbar)"] = 4.2;
                p["mass::c"] = 1.27;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "KMPW2010");

                BToKstarGamma d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(),             +5.45306e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_cp_averaged(), +5.47311e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.s_kstar_gamma(),               -3.94778e-2, eps);
                TEST_CHECK_RELATIVE_ERROR(d.c_kstar_gamma(),               +3.66320e-3, eps);
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
                p["Abs{c7}"] = 0.331;
                p["Arg{c7}"] = -M_PI / 2.0;
                p["Abs{c7'}"] = 0.00659; // m_s(m_b) / m_b(m_b) * Abs{c7} = 85 / 4200 * Abs{c7}
                p["Arg{c7'}"] = -M_PI / 2.0;
                p["c8"] = -0.181;
                // PDG 2010 CKM parameters
                p["CKM::A"] = 0.812;
                p["CKM::lambda"] = 0.22543;
                p["CKM::rhobar"] = 0.144;
                p["CKM::etabar"] = 0.342;
                p["decay-constant::B_d"] = 0.200;
                p["mass::b(MSbar)"] = 4.2;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("form-factors", "KMPW2010");

                BToKstarGamma d(p, oo);

                const double eps = 1e-4;

                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio(),             +5.65584e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.branching_ratio_cp_averaged(), +4.00005e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(d.s_kstar_gamma(),               +4.72504e-2, eps);
                TEST_CHECK_RELATIVE_ERROR(d.c_kstar_gamma(),               -4.13944e-1, eps);
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
                "Abs{c7}",  "Arg{c7}",  "Abs{c7'}",  "Arg{c7'}",
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

            std::vector<ObservablePtr> observables;
            observables.push_back(Observable::make("B->K^*gamma::BR,q=d",  p, k, o));
            observables.push_back(Observable::make("B->K^*gamma::S_K^*gamma,q=d", p, k, o));
            observables.push_back(Observable::make("B->K^*gamma::C_K^*gamma,q=d", p, k, o));
            observables.push_back(Observable::make("B->K^*gamma::A_I", p, k, o));

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-gamma_TEST-btokstargamma.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->K^*gamma --" << std::endl;
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
                std::cout << "-- Verifying test case data for B->K^*gamma --" << std::endl;
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
} b_to_kstar_gamma_bobeth_compatibility_test;
