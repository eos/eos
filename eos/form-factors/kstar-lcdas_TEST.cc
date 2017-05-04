/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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
#include <eos/form-factors/kstar-lcdas.hh>

#include <eos/utils/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class KstarLCDAsTest :
    public TestCase
{
    public:
        KstarLCDAsTest() :
            TestCase("kstar_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"] = 0.1176;
            p["K^*::a_1_para@1GeV"] = 0.03; // taken from [BBL2007], p. 23, table 1
            p["K^*::a_2_para@1GeV"] = 0.11;
            p["K^*::a_1_perp@1GeV"] = 0.04;
            p["K^*::a_2_perp@1GeV"] = 0.10;
            p["K^*::f_perp@1GeV"]   = 0.159; // taken from [BSZ2015], table 1

            /* Diagnostics */
            {
                KstarLCDAs kstar(p, Options{ });
                Diagnostics diagnostics = kstar.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94850, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92874, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91708, 1e-5), // c_rge(mu = 4.0 GeV)
                    std::make_pair(+0.90893, 1e-5), // c_rge(mu = 5.0 GeV)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* Twist 2 */
            {
                KstarLCDAs kstar(p, Options{ });

                // coefficients at mu = 1.0 GeV, and 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.03000,  kstar.a_1_para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.02486,  kstar.a_1_para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.11000,  kstar.a_2_para(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.08200,  kstar.a_2_para(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.04000,  kstar.a_1_perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.03238,  kstar.a_1_perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.10000,  kstar.a_2_perp(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.07368,  kstar.a_2_perp(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.15900,  kstar.f_perp(1.0),     eps);
                TEST_CHECK_NEARLY_EQUAL( 0.14637,  kstar.f_perp(2.0),     eps);

                // phi_2_para LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL( 0.69714, kstar.phi_2_para(0.1, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.17306, kstar.phi_2_para(0.3, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.25250, kstar.phi_2_para(0.5, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.26378, kstar.phi_2_para(0.7, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.77490, kstar.phi_2_para(0.9, 1.0),  eps);

                // phi_2_perp LCDA at various u values for mu = 1.0
                TEST_CHECK_NEARLY_EQUAL( 0.66636, kstar.phi_2_perp(0.1, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.16172, kstar.phi_2_perp(0.3, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.27500, kstar.phi_2_perp(0.5, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 1.28268, kstar.phi_2_perp(0.7, 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.77004, kstar.phi_2_perp(0.9, 1.0),  eps);
            }
        }
} kstar_lcdas_test;
