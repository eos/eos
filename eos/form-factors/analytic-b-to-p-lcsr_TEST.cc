/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
 * Copyright (c) 2018 Nico Gubernari
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
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/mesonic.hh>

#include <vector>
#include <utility>

using namespace test;
using namespace eos;

class LCSRFormFactorsTest :
    public TestCase
{
    public:
        LCSRFormFactorsTest() :
            TestCase("lcsr_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> pi diagnostic values */
            {
                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.3174;
                p["B::lambda_H^2"]            = 1.2696;
                p["mass::d(2GeV)"]            = 0.0048;
                p["mass::u(2GeV)"]            = 0.0032;
                p["mass::B_d"]                = 5.2795;
                p["mass::pi^+"]               = 0.13957;
                p["decay-constant::B_d"]      = 0.180;
                p["decay-constant::pi"]       = 0.1302;
                p["B->pi::mu@B-LCSR"]         = 1.0;
                p["B->pi::s_0^+,0@B-LCSR"]    = 0.7;
                p["B->pi::s_0^+,1@B-LCSR"]    = 0.0;
                p["B->pi::s_0^+/-,0@B-LCSR"]  = 0.7;
                p["B->pi::s_0^+/-,1@B-LCSR"]  = 0.0;
                p["B->pi::s_0^T,0@B-LCSR"]    = 0.7;
                p["B->pi::s_0^T,1@B-LCSR"]    = 0.0;
                p["B->pi::M^2@B-LCSR"]        = 1.0;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "zero" }
                };

                AnalyticFormFactorBToPLCSR<lcsr::BToPi> ff{ p, o };
                auto diagnostics = ff.diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    //ATTENTION: the test cases are evaluated with mathematica using m_v = m_d
                    std::make_pair(0.005340672,  1.0e-7), // m_v(mu) in the MSbar scheme

                    std::make_pair(0.7,          1.0e-7), // s_0 value for fT

                    /* f_+ */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-6.306061e-1,  1.0e-6), // I1_fp_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.306061e-1,  1.0e-6), // I1_fp_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.306061e-1,  1.0e-6), // I1_fp_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.969100e-1,  1.0e-6), // I1_fp_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.969100e-1,  1.0e-6), // I1_fp_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.969100e-1,  1.0e-6), // I1_fp_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(2.121437,  1.0e-6), // I2_fp_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.121437,  1.0e-6), // I2_fp_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.121437,  1.0e-6), // I2_fp_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.458780,  1.0e-6), // I2_fp_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.458780,  1.0e-6), // I2_fp_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.458780,  1.0e-6), // I2_fp_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(24.29358,  1.0e-5), // I2d1_fp_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(24.29358,  1.0e-5), // I2d1_fp_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(24.29358,  1.0e-5), // I2d1_fp_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.038187,  1.0e-6), // I2d1_fp_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.038187,  1.0e-6), // I2d1_fp_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.038187,  1.0e-6), // I2d1_fp_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_fp_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_fp_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_fp_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_fp_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_fp_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_fp_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-4.813004,  1.0e-6), // I2d1_fp_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.813004,  1.0e-6), // I2d1_fp_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.813004,  1.0e-6), // I2d1_fp_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.871988,  1.0e-6), // I2d1_fp_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.871988,  1.0e-6), // I2d1_fp_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.871988,  1.0e-6), // I2d1_fp_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(7.087783e-6,  1.0e-6), // I3_fp_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.087783e-6,  1.0e-6), // I3_fp_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.087783e-6,  1.0e-6), // I3_fp_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.999202e-5,  1.0e-6), // I3_fp_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.999202e-5,  1.0e-6), // I3_fp_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.999202e-5,  1.0e-6), // I3_fp_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(2.934497e-4,  1.0e-6), // I3d1_fp_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.934497e-4,  1.0e-6), // I3d1_fp_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.934497e-4,  1.0e-6), // I3d1_fp_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.239029e-4,  1.0e-6), // I3d1_fp_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(3.239029e-4,  1.0e-6), // I3d1_fp_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(3.239029e-4,  1.0e-6), // I3d1_fp_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(3.229436e-3,  1.0e-6),     // I3d2_fp_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.229436e-3,  1.0e-6),     // I3d2_fp_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.229436e-3,  1.0e-6),     // I3d2_fp_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.055932e-3,  1.0e-6), // I3d2_fp_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.055932e-3,  1.0e-6), // I3d2_fp_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.055932e-3,  1.0e-6), // I3d2_fp_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-9.86142e-2,  1.0e-6), // I3_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-9.86142e-2,  1.0e-6), // I3_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-9.86142e-2,  1.0e-6), // I3_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.040432e-1,  1.0e-6), // I3_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.040432e-1,  1.0e-6), // I3_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.040432e-1,  1.0e-6), // I3_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-6.741819,  1.0e-6), // I3d1_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.741819,  1.0e-6), // I3d1_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.741819,  1.0e-6), // I3d1_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.862556e1,  1.0e-5), // I3d1_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.862556e1,  1.0e-5), // I3d1_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.862556e1,  1.0e-5), // I3d1_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar */
                    std::make_pair(-2.754366e2,  1.0e-4), // I3d2_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.754366e2,  1.0e-4), // I3d2_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.754366e2,  1.0e-4), // I3d2_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.925539e2,  1.0e-4), // I3d2_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.925539e2,  1.0e-4), // I3d2_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.925539e2,  1.0e-4), // I3d2_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(8.792167e-6,  1.0e-6), // I4_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.792167e-6,  1.0e-6), // I4_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(8.792167e-6,  1.0e-6), // I4_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.619631e-5,  1.0e-6), // I4_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.619631e-5,  1.0e-6), // I4_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.619631e-5,  1.0e-6), // I4_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(6.102403e-4,  1.0e-6), // I4d1_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.102403e-4,  1.0e-6), // I4d1_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(6.102403e-4,  1.0e-6), // I4d1_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.793886e-3,  1.0e-6), // I4d1_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.793886e-3,  1.0e-6), // I4d1_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.793886e-3,  1.0e-6), // I4d1_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(2.582168e-2,  1.0e-6), // I4d2_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.582168e-2,  1.0e-6), // I4d2_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.582168e-2,  1.0e-6), // I4d2_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.110889e-2,  1.0e-6), // I4d2_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(3.110889e-2,  1.0e-6), // I4d2_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(3.110889e-2,  1.0e-6), // I4d2_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar */
                    std::make_pair(3.26455e-1,   1.0e-6), // I4d3_fp_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.26455e-1,   1.0e-6), // I4d3_fp_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.26455e-1,   1.0e-6), // I4d3_fp_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.188097e-2,  1.0e-6), // I4d3_fp_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.188097e-2,  1.0e-6), // I4d3_fp_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.188097e-2,  1.0e-6),  // I4d3_fp_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair(-6.0901e-2,  1.0e-6),  // I2_fp_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.2750e-1,  1.0e-5),  // I2_fp_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-8.0096e-3,  1.0e-7),  // I2_fp_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6780e-2,  1.0e-6),  // I2_fp_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair( 1.4154e-2,  1.0e-6),  // I2_fp_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.9663e-2,  1.0e-6),  // I2_fp_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.4227e-2,  1.0e-6),  // I2_fp_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5555e-1,  1.0e-5),  // I2_fp_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(-9.1040e-1,  1.0e-5),  // I3_fp_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.9047   ,  1.0e-4),  // I3_fp_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.7756   ,  1.0e-4),  // I3_fp_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0003e+1,  1.0e-3),  // I3_fp_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair( 4.5024,     1.0e-4),  // I3d1A_fp_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.4422,     1.0e-4),  // I3d1A_fp_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.3800,     1.0e-4),  // I3d1A_fp_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.9979,     1.0e-4),  // I3d1A_fp_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair( 3.7021e2,  1.0e-2),   // I3d1B_fp_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.5009e5,  1.0e+1),   // I3d1B_fp_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 1.5620e-3, 1.0e-7),   // I3d1C_fp_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 3.2735e-3, 1.0e-7),   // I3d1C_fp_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-8.8243e-4,  1.0e-8),  // I4_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.1837e-2,  1.0e-6),  // I4_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.5300e-3,  1.0e-7),  // I4_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.1816e-2,  1.0e-6),  // I4_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.7024e-3,  1.0e-7),  // I4d1A_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.0107e-2,  1.0e-6),  // I4d1A_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0310e-3,  1.0e-7),  // I4d1A_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2811e-2,  1.0e-6),  // I4d1A_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(-1.9188   ,  1.0e-4),  // I4d1B_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.3916e+3,  1.0e-1),  // I4d1B_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fp_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fp_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair( 1.7587e-2,  1.0e-6),  // I4d2A_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.3309e-1,  1.0e-5),  // I4d2A_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.5760e-3,  1.0e-7),  // I4d2A_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.9196e-2,  1.0e-6),  // I4d2A_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair( 4.6189e1,  1.0e-3),   // I4d2B_fp_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.6190e4,  1.0e-0),   // I4d2B_fp_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 2.2997e-2,  1.0e-5),  // I4d2C_fp_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 6.1266e-2,  1.0e-4),  // I4d2C_fp_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_fp_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair( 5.2280e-2,  1.0e-6),  // I2_fp_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5632e-1,  1.0e-5),  // I2_fp_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1192e-3,  1.0e-7),  // I2_fp_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.0459e-3,  1.0e-7),  // I2_fp_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-3.6625e-2,  1.0e-6),  // I2_fp_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0951e-1,  1.0e-5),  // I2_fp_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2931e-1,  1.0e-5),  // I2_fp_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.0790e-1,  1.0e-5),  // I2_fp_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.3893   ,  1.0e-4),  // I3_fp_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.5699   ,  1.0e-4),  // I3_fp_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.4140   ,  1.0e-4),  // I3_fp_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0961e+1,  1.0e-3),  // I3_fp_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.3641   ,  1.0e-4),  // I3d1A_fp_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8582e+1,  1.0e-3),  // I3d1A_fp_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.3273   ,  1.0e-4),  // I3d1A_fp_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6543e+1,  1.0e-3),  // I3d1A_fp_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair( 7.1751e2,  1.0e-2),   // I3d1B_fp_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1178e5,  1.0e+1),   // I3d1B_fp_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair( 7.2331e-1,  1.0e-5),  // I3d1C_fp_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.5158   ,  1.0e-4),  // I3d1C_fp_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 2.1594e-1,  1.0e-5),  // I3_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.0469e-1,  1.0e-5),  // I3_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6829   ,  1.0e-4),  // I3_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.6415   ,  1.0e-4),  // I3_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-2.2010   ,  1.0e-4),  // I3d1A_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.4128   ,  1.0e-4),  // I3d1A_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.3618   ,  1.0e-4),  // I3d1A_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1521e+1,  1.0e-3),  // I3d1A_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-2.6259e+1,  1.0e-3),  // I3d1B_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.3864e+4,  1.0e-0),  // I3d1B_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fp_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fp_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 1.4158e+1,  1.0e-3),  // I4_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.6203e+1,  1.0e-3),  // I4_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7591e+2,  1.0e-2),  // I4_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.0102e+2,  1.0e-2),  // I4_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.6581e+2,  1.0e-2),  // I4d1A_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.9064e+2,  1.0e-2),  // I4d1A_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.5317e+2,  1.0e-2),  // I4d1A_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5164e+3,  1.0e-1),  // I4d1A_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(-1.7217e+3,  1.0e-1),  // I4d1B_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-9.0902e+5,  1.0e+1),  // I4d1B_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fp_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fp_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 8.6802e+2,  1.0e-2),  // I4d2A_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.2511e+3,  1.0e-1),  // I4d2A_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.2257e+2,  1.0e-2),  // I4d2A_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0477e+3,  1.0e-1),  // I4d2A_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 6.8324e4,  1.0e-0),   // I4d2B_fp_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.6709e7,  1.0e+3),   // I4d2B_fp_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-8.9097e1,   1.0e-3),  // I4d2C_fp_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.3736e2,   1.0e-2),  // I4d2C_fp_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_fp_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 9.1997e-2,  1.0e-6),  // I2_fp_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.7250e-2,  1.0e-6),  // I2_fp_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.6328e-1,  1.0e-5),  // I2_fp_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0763e-1,  1.0e-5),  // I2_fp_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 4.0212   ,  1.0e-4),  // I3_fp_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.0653   ,  1.0e-4),  // I3_fp_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.0250e+1,  1.0e-3),  // I3_fp_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 9.0754   ,  1.0e-4),  // I3_fp_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.1736e+1,  1.0e-3),  // I3d1A_fp_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-8.9034   ,  1.0e-4),  // I3d1A_fp_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.5998e+1,  1.0e-3),  // I3d1A_fp_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0576e+1,  1.0e-3),  // I3d1A_fp_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 9.1786e1,  1.0e-3),   // I3d1B_fp_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0981e4,  1.0e-0),   // I3d1B_fp_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 9.2825e-1,  1.0e-5),  // I3d1C_fp_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 3.8906e-1,  1.0e-5),  // I3d1C_fp_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.0736e-2,  1.0e-6),  // I2_fp_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.3456e-2,  1.0e-6),  // I2_fp_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.2091e-2,  1.0e-6),  // I2_fp_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.7252e-2,  1.0e-6),  // I2_fp_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 1.9663e+1,  1.0e-3),  // I3_fp_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.7532   ,  1.0e-4),  // I3_fp_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.9783e+1,  1.0e-3),  // I3_fp_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.4162e+1,  1.0e-3),  // I3_fp_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-9.8456e+1,  1.0e-3),  // I3d1A_fp_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.0455e+1,  1.0e-3),  // I3d1A_fp_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-8.0809e+1,  1.0e-3),  // I3d1A_fp_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.2614e+1,  1.0e-3),  // I3d1A_fp_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair( 5.4237e2,  1.0e-2),   // I3d1B_fp_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.4890e4,  1.0e-0),   // I3d1B_fp_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 3.6938,  1.0e-4),     // I3d1C_fp_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.5482,  1.0e-4),     // I3d1C_fp_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* f_± */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_fpm_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_fpm_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_fpm_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_fpm_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_fpm_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_fpm_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(-1.744575e-1,  1.0e-6), // I2_fpm_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.744575e-1,  1.0e-6), // I2_fpm_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.744575e-1,  1.0e-6), // I2_fpm_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.246749e-1,  1.0e-6), // I2_fpm_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.246749e-1,  1.0e-6), // I2_fpm_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.246749e-1,  1.0e-6), // I2_fpm_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(-6.596762,  1.0e-6), // I2d1_fpm_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.596762,  1.0e-6), // I2d1_fpm_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.596762,  1.0e-6), // I2d1_fpm_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.278834,  1.0e-6), // I2d1_fpm_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.278834,  1.0e-6), // I2d1_fpm_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-5.278834,  1.0e-6), // I2d1_fpm_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_fpm_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_fpm_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_fpm_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_fpm_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_fpm_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_fpm_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_fpm_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_fpm_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_fpm_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.067506,  1.0e-6), // I2d1_fpm_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.067506,  1.0e-6), // I2d1_fpm_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.067506,  1.0e-6), // I2d1_fpm_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(6.783834e-6,  1.0e-6), // I3_fpm_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.783834e-6,  1.0e-6), // I3_fpm_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(6.783834e-6,  1.0e-6), // I3_fpm_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.825358e-5,  1.0e-6), // I3_fpm_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.825358e-5,  1.0e-6), // I3_fpm_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.825358e-5,  1.0e-6), // I3_fpm_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(2.735417e-4,  1.0e-6), // I3d1_fpm_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.735417e-4,  1.0e-6), // I3d1_fpm_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.735417e-4,  1.0e-6), // I3d1_fpm_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.721174e-4,  1.0e-6), // I3d1_fpm_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.721174e-4,  1.0e-6), // I3d1_fpm_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.721174e-4,  1.0e-6), // I3d1_fpm_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(2.442048e-3,  1.0e-6),     // I3d2_fpm_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.442048e-3,  1.0e-6),     // I3d2_fpm_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.442048e-3,  1.0e-6),     // I3d2_fpm_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.780826e-3,  1.0e-6), // I3d2_fpm_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.780826e-3,  1.0e-6), // I3d2_fpm_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.780826e-3,  1.0e-6), // I3d2_fpm_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.050510e-1,  1.0e-6), // I3_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.050510e-1,  1.0e-6), // I3_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.050510e-1,  1.0e-6), // I3_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.66655,  1.0e-6), // I3d1_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.66655,  1.0e-6), // I3d1_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.66655,  1.0e-6), // I3d1_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar */
                    std::make_pair(5.266027e1,  1.0e-5), // I3d2_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.266027e1,  1.0e-5), // I3d2_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.266027e1,  1.0e-5), // I3d2_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.420043e2,  1.0e-4), // I3d2_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.420043e2,  1.0e-4), // I3d2_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.420043e2,  1.0e-4), // I3d2_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(-7.230287e-7,  1.0e-6), // I4_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-7.230287e-7,  1.0e-6), // I4_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-7.230287e-7,  1.0e-6), // I4_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-9.706098e-6,  1.0e-6), // I4_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-9.706098e-6,  1.0e-6), // I4_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-9.706098e-6,  1.0e-6), // I4_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(-6.924356e-5,  1.0e-6), // I4d1_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.924356e-5,  1.0e-6), // I4d1_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.924356e-5,  1.0e-6), // I4d1_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.424789e-4,  1.0e-6), // I4d1_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.424789e-4,  1.0e-6), // I4d1_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.424789e-4,  1.0e-6), // I4d1_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(-4.808265e-3,  1.0e-6), // I4d2_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.808265e-3,  1.0e-6), // I4d2_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.808265e-3,  1.0e-6), // I4d2_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.412735e-2,  1.0e-6), // I4d2_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.412735e-2,  1.0e-6), // I4d2_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.412735e-2,  1.0e-6), // I4d2_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar */
                    std::make_pair(-2.031646e-1,  1.0e-6), // I4d3_fpm_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.031646e-1,  1.0e-6), // I4d3_fpm_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.031646e-1,  1.0e-6), // I4d3_fpm_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.467625e-1,  1.0e-6), // I4d3_fpm_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.467625e-1,  1.0e-6), // I4d3_fpm_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.467625e-1,  1.0e-6), // I4d3_fpm_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 6.6758e-2,  1.0e-6),  // I2_fpm_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4016e-1,  1.0e-5),  // I2_fpm_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.7754e-3,  1.0e-7),  // I2_fpm_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.8400e-2,  1.0e-6),  // I2_fpm_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-9.0011e-4,  1.0e-8),  // I2_fpm_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.8863e-3,  1.0e-7),  // I2_fpm_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.7203e-3,  1.0e-7),  // I2_fpm_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.8921e-3,  1.0e-7),  // I2_fpm_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 6.9665e-1,  1.0e-5),  // I3_fpm_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4630   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.6519   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.6580   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-1.3654,     1.0e-4),  // I3d1A_fpm_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.8585,     1.0e-4),  // I3d1A_fpm_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.0672,     1.0e-4),  // I3d1A_fpm_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 19.006,     1.0e-3),  // I3d1A_fpm_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-2.8078e2,  1.0e-2),   // I3d1B_fpm_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1384e5,  1.0e+1),   // I3d1B_fpm_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 1.5123e-3, 1.0e-7),   // I3d1C_fpm_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 3.1694e-3, 1.0e-7),   // I3d1C_fpm_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair( 8.1649e-5,  1.0e-9),  // I4_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.3611e-4,  1.0e-8),  // I4_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.2379e-4,  1.0e-8),  // I4_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.1722e-3,  1.0e-7),  // I4_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 1.4358e-3,  1.0e-7),  // I4d1A_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.1431e-2,  1.0e-6),  // I4d1A_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.9018e-3,  1.0e-7),  // I4d1A_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.3312e-1,  1.0e-5),  // I4d1A_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 1.1892e-1,  1.0e-5),  // I4d1B_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.6249e+1,  1.0e-3),  // I4d1B_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fpm_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fpm_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.3131e-2,  1.0e-6),  // I4d2A_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.7699e-1,  1.0e-5),  // I4d2A_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.2914e-2,  1.0e-6),  // I4d2A_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.1887e-1,  1.0e-5),  // I4d2A_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair( 5.3349  ,  1.0e-4),   // I4d2B_fpm_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 4.9216e3,  1.0e-1),   // I4d2B_fpm_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(-1.4994e-3,  1.0e-7),  // I4d2C_fpm_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-3.9945e-3,  1.0e-7),  // I4d2C_fpm_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_fpm_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair( 1.5850e-1,  1.0e-5),  // I2_fpm_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.7393e-1,  1.0e-5),  // I2_fpm_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.4250e-3,  1.0e-7),  // I2_fpm_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5298e-2,  1.0e-6),  // I2_fpm_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 2.3291e-3,  1.0e-7),  // I2_fpm_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 6.9640e-3,  1.0e-7),  // I2_fpm_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.2234e-3,  1.0e-7),  // I2_fpm_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.9580e-2,  1.0e-6),  // I2_fpm_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.5628   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.7683   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.4870   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3181e+1,  1.0e-3),  // I3_fpm_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 9.1284   ,  1.0e-4),  // I3d1A_fpm_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6025e+1,  1.0e-3),  // I3d1A_fpm_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.2352e+1,  1.0e-3),  // I3d1A_fpm_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2388e+1,  1.0e-3),  // I3d1A_fpm_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-4.3315e1,  1.0e-3),   // I3d1B_fpm_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.2785e4,  1.0e+0),   // I3d1B_fpm_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-3.1613   ,  1.0e-4),  // I3d1C_fpm_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-6.6250   ,  1.0e-4),  // I3d1C_fpm_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair(-4.5941e-1,  1.0e-5),  // I3_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.4999   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7073   ,  1.0e-4),  // I3_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6258e+1,  1.0e-3),  // I3_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair( 3.7625   ,  1.0e-4),  // I3d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0642e+1,  1.0e-3),  // I3d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1471   ,  1.0e-4),  // I3d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-8.0410   ,  1.0e-4),  // I3d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair( 5.5717e+1,  1.0e-3),  // I3d1B_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.9418e+4,  1.0e-0),  // I3d1B_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fpm_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fpm_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(-9.0291e-1,  1.0e-5),  // I4_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.9725   ,  1.0e-4),  // I4_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.1196e+1,  1.0e-3),  // I4_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.1992e+1,  1.0e-3),  // I4_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.9592e+1,  1.0e-3),  // I4d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.7049e+1,  1.0e-3),  // I4d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.3937e+2,  1.0e-2),  // I4d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.7044e+2,  1.0e-2),  // I4d1A_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.0393e+2,  1.0e-2),  // I4d1B_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.4875e+4,  1.0e+0),  // I4d1B_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fpm_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fpm_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 5.8876e+2,  1.0e-2),  // I4d2A_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7440e+3,  1.0e-1),  // I4d2A_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5559e+3,  1.0e-1),  // I4d2A_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.1896e+3,  1.0e-1),  // I4d2A_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 3.1145e3,  1.0e-1),   // I4d2B_fpm_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.8394e6,  1.0e+2),   // I4d2B_fpm_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair( 5.9543  ,   1.0e-4),  // I4d2C_fpm_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.5863e1,   1.0e-3),  // I4d2C_fpm_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_fpm_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-5.8503e-3,  1.0e-7),  // I2_fpm_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0048e-3,  1.0e-7),  // I2_fpm_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.9461e-2,  1.0e-6),  // I2_fpm_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3204e-2,  1.0e-6),  // I2_fpm_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 3.8927   ,  1.0e-4),  // I3_fpm_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9984   ,  1.0e-4),  // I3_fpm_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9605e+1,  1.0e-3),  // I3_fpm_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.7849   ,  1.0e-4),  // I3_fpm_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.5328e+1,  1.0e-3),  // I3d1A_fpm_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0823e+1,  1.0e-3),  // I3d1A_fpm_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.6733e+1,  1.0e-3),  // I3d1A_fpm_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.9908e+1,  1.0e-3),  // I3d1A_fpm_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 8.9111e1,  1.0e-3),   // I3d1B_fpm_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0662e4,  1.0e-0),   // I3d1B_fpm_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 8.9626e-1,  1.0e-5),  // I3d1C_fpm_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 3.7565e-1,  1.0e-5),  // I3d1C_fpm_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair(-1.3187e-3,  1.0e-7),  // I2_fpm_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.7635e-3,  1.0e-7),  // I2_fpm_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0407e-3,  1.0e-7),  // I2_fpm_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.2767e-3,  1.0e-7),  // I2_fpm_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair(-3.2005e-1,  1.0e-5),  // I3_fpm_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3295   ,  1.0e-4),  // I3_fpm_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.9056   ,  1.0e-4),  // I3_fpm_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0903e-1,  1.0e-5),  // I3_fpm_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-3.4147e+1,  1.0e-3),  // I3d1A_fpm_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5136e+1,  1.0e-3),  // I3d1A_fpm_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0506e+2,  1.0e-2),  // I3d1A_fpm_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-8.7223e+1,  1.0e-3),  // I3d1A_fpm_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-4.1135e2,  1.0e-2),   // I3d1B_fpm_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-4.9215e4,  1.0e-0),   // I3d1B_fpm_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 3.5764,  1.0e-4),     // I3d1C_fpm_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.4990,  1.0e-4),     // I3d1C_fpm_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* f_T */

                    /* 2 particle */

                    /* I_1 phi_bar */
                    std::make_pair(-7.928177e-2,  1.0e-6), // I1_fT_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-7.928177e-2,  1.0e-6), // I1_fT_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-7.928177e-2,  1.0e-6), // I1_fT_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-9.588403e-2,  1.0e-6), // I1_fT_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-9.588403e-2,  1.0e-6), // I1_fT_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-9.588403e-2,  1.0e-6), // I1_fT_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(2.501327,  1.0e-6), // I2_fT_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.121435,  1.0e-6), // I2_fT_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.741543,  1.0e-6), // I2_fT_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.896509,  1.0e-6), // I2_fT_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.458777,  1.0e-6), // I2_fT_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.021046,  1.0e-6), // I2_fT_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(2.860946e1,  1.0e-4), // I2d1_fT_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.429355e1,  1.0e-4), // I2d1_fT_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.997764e1,  1.0e-4), // I2d1_fT_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.66970,  1.0e-5), // I2d1_fT_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.03819,  1.0e-5), // I2d1_fT_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.40668,  1.0e-5), // I2d1_fT_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(3.685383e-3,  1.0e-6), // I2_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.685383e-3,  1.0e-6), // I2_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.685383e-3,  1.0e-6), // I2_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.355562e-2,  1.0e-6), // I2_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.355562e-2,  1.0e-6), // I2_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.355562e-2,  1.0e-6), // I2_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_bar */
                    std::make_pair(2.557924e-1,  1.0e-6), // I2d1_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.557924e-1,  1.0e-6), // I2d1_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.557924e-1,  1.0e-6), // I2d1_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.519371e-1,  1.0e-6), // I2d1_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.519371e-1,  1.0e-6), // I2d1_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.519371e-1,  1.0e-6), // I2d1_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-1.162735e-1,  1.0e-6), // I3_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-9.861442e-2,  1.0e-6), // I3_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-8.095529e-2,  1.0e-6), // I3_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.115810e-1,  1.0e-6), // I3_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.040445e-1,  1.0e-6), // I3_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.965080e-1,  1.0e-6), // I3_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-7.947512,  1.0e-6), // I3d1_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.741835,  1.0e-6), // I3d1_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-5.536158,  1.0e-6), // I3d1_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.191921e1,  1.0e-5), // I3d1_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.862561e1,  1.0e-5), // I3d1_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.533200e1,  1.0e-5), // I3d1_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar */
                    std::make_pair(-3.244968e2,  1.0e-4), // I3d2_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.754373e2,  1.0e-4), // I3d2_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.263777e2,  1.0e-4), // I3d2_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.429136e2,  1.0e-4), // I3d2_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.925548e2,  1.0e-4), // I3d2_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.421960e2,  1.0e-4), // I3d2_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(1.036660e-5,  1.0e-6), // I4_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.792157e-6,  1.0e-6), // I4_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.217719e-6,  1.0e-6), // I4_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(6.620075e-5,  1.0e-6), // I4_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.619624e-5,  1.0e-6), // I4_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.619173e-5,  1.0e-6), // I4_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(7.193745e-4,  1.0e-6), // I4d1_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.102396e-4,  1.0e-6), // I4d1_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.011048e-4,  1.0e-6), // I4d1_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.111174e-3,  1.0e-6), // I4d1_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.793884e-3,  1.0e-6), // I4d1_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.476593e-3,  1.0e-6), // I4d1_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(3.042181e-2,  1.0e-6), // I4d2_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.582165e-2,  1.0e-6), // I4d2_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.122149e-2,  1.0e-6), // I4d2_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.648225e-2,  1.0e-5), // I4d2_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(3.111705e-2,  1.0e-5), // I4d2_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.574222e-2,  1.0e-5), // I4d2_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar */
                    std::make_pair(3.827981e-1,  1.0e-6), // I4d3_fT_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.264544e-1,  1.0e-6), // I4d3_fT_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.701108e-1,  1.0e-6), // I4d3_fT_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.075999e-2,  1.0e-6), // I4d3_fT_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.188125e-2,  1.0e-6), // I4d3_fT_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.003306e-3,  1.0e-6), // I4d3_fT_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_1 phi_3 */
                    std::make_pair( 4.5100e-3,  1.0e-7),  // I1_fT_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.4514e-3,  1.0e-7),  // I1_fT_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.9303e-4,  1.0e-8),  // I1_fT_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2428e-3,  1.0e-7),  // I1_fT_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_3 */
                    std::make_pair(-1.0000e-1,  1.0e-5),  // I2_fT_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0957e-1,  1.0e-5),  // I2_fT_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.3149e-2,  1.0e-6),  // I2_fT_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.7557e-2,  1.0e-6),  // I2_fT_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair( 5.6610e-2,  1.0e-6),  // I2_fT_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1858e-1,  1.0e-5),  // I2_fT_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9690e-1,  1.0e-5),  // I2_fT_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.2210e-1,  1.0e-5),  // I2_fT_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(-1.2552   ,  1.0e-4),  // I3_fT_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.6292   ,  1.0e-4),  // I3_fT_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-6.5831   ,  1.0e-4),  // I3_fT_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3794e+1,  1.0e-3),  // I3_fT_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair( 6.6021,     1.0e-4),  // I3d1A_fT_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3838e+1,  1.0e-3),  // I3d1A_fT_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.3604,     1.0e-4),  // I3d1A_fT_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1238e+1,  1.0e-3),  // I3d1A_fT_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3*/
                    std::make_pair( 5.0903e2,  1.0e-2),   // I3d1B_fT_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.0637e5,  1.0e+1),   // I3d1B_fT_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 6.4106e-4, 1.0e-8),   // I3d1C_fT_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.3434e-3, 1.0e-7),   // I3d1C_fT_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 3.2970e-5,  1.0e-9),  // I3_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.3918e-4,  1.0e-8),  // I3_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7290e-4,  1.0e-8),  // I3_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.3031e-3,  1.0e-7),  // I3_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-1.0583e-4,  1.0e-8),  // I3d1A_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.4097e-3,  1.0e-7),  // I3d1A_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1347e-4,  1.0e-8),  // I3d1A_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.8435e-3,  1.0e-7),  // I3d1A_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 7.0995e-2,  1.0e-6),  // I3d1B_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.1487e+1,  1.0e-3),  // I3d1B_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fT_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fT_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-7.3106e-4,  1.0e-8),  // I4_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.7380e-3,  1.0e-7),  // I4_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.8338e-3,  1.0e-7),  // I4_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.1067e-2,  1.0e-6),  // I4_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.0901e-3,  1.0e-7),  // I4d1A_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.1162e-2,  1.0e-6),  // I4d1A_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-8.3439e-4,  1.0e-8),  // I4d1A_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1114e-2,  1.0e-6),  // I4d1A_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(-1.5741   ,  1.0e-4),  // I4d1B_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1416e+3,  1.0e-1),  // I4d1B_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fT_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fT_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair( 1.4062e-2,  1.0e-6),  // I4d2A_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8731e-1,  1.0e-5),  // I4d2A_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9140e-3,  1.0e-7),  // I4d2A_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.8816e-2,  1.0e-6),  // I4d2A_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair( 3.7837e1,  1.0e-3),   // I4d2B_fT_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.3247e4,  1.0e-0),   // I4d2B_fT_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 1.8883e-2,  1.0e-6),  // I4d2C_fT_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 5.0306e-2,  1.0e-6),  // I4d2C_fT_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_fT_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-7.8150e-3,  1.0e-7),  // I2_fT_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-8.1889e-2,  1.0e-6),  // I2_fT_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-8.8049e-3,  1.0e-7),  // I2_fT_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.2261e-2,  1.0e-6),  // I2_fT_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair( 1.7328e-1,  1.0e-5),  // I3_fT_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8157   ,  1.0e-4),  // I3_fT_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9523e-1,  1.0e-5),  // I3_fT_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0457   ,  1.0e-4),  // I3_fT_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 1.8136e-1,  1.0e-4),  // I3d1A_fT_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9003   ,  1.0e-3),  // I3d1A_fT_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.0433e-1,  1.0e-4),  // I3d1A_fT_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.1410   ,  1.0e-3),  // I3d1A_fT_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair( 5.8771e2,  1.0e-2),   // I3d1B_fT_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.7347e5,  1.0e+1),   // I3d1B_fT_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair( 3.0751  ,  1.0e-4),   // I3d1C_fT_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 6.4444  ,  1.0e-4),   // I3d1C_fT_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair(-1.5989e-2,  1.0e-6),  // I2_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.2197e-2,  1.0e-6),  // I2_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.9864e-1,  1.0e-5),  // I2_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.6584e-1,  1.0e-5),  // I2_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair(-1.6964e-1,  1.0e-5),  // I3_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.5379e-1,  1.0e-5),  // I3_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1074   ,  1.0e-4),  // I3_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.0033   ,  1.0e-4),  // I3_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair( 1.9895   ,  1.0e-4),  // I3d1A_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.8887   ,  1.0e-4),  // I3d1A_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.6632   ,  1.0e-4),  // I3d1A_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.8270e+1,  1.0e-3),  // I3d1A_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair( 2.0583e+1,  1.0e-3),  // I3d1B_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0868e+4,  1.0e-0),  // I3d1B_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fT_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_fT_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 1.1622e+1,  1.0e-3),  // I4_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.7942e+1,  1.0e-3),  // I4_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.4439e+2,  1.0e-2),  // I4_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.1130e+2,  1.0e-2),  // I4_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.3593e+2,  1.0e-2),  // I4d1A_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.0224e+2,  1.0e-2),  // I4d1A_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.5190e+2,  1.0e-2),  // I4d1A_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2386e+3,  1.0e-1),  // I4d1A_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(-1.4102e+3,  1.0e-1),  // I4d1B_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-7.4457e+5,  1.0e+1),  // I4d1B_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fT_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_fT_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 7.1453e+2,  1.0e-2),  // I4d2A_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8554e+3,  1.0e-1),  // I4d2A_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.0933e+2,  1.0e-2),  // I4d2A_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0392e+3,  1.0e-1),  // I4d2A_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 5.5964e4,  1.0e-0),   // I4d2B_fT_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.3673e7,  1.0e+3),   // I4d2B_fT_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-7.3283e+1,  1.0e-3),  // I4d2C_fT_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.9523e+2,  1.0e-2),  // I4d2C_fT_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_fT_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-2.1643e-5,  1.0e-9),  // I2_fT_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.5357e-5,  1.0e-9),  // I2_fT_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.3494e-5,  1.0e-9),  // I2_fT_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.0194e-5,  1.0e-9),  // I2_fT_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 4.7997e-4,  1.0e-7),  // I3_fT_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0059e-3,  1.0e-6),  // I3_fT_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.4278e-4,  1.0e-7),  // I3_fT_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5566e-3,  1.0e-6),  // I3_fT_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair( 9.9742e-4,  1.0e-8),  // I3d1A_fT_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.0902e-3,  1.0e-7),  // I3d1A_fT_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5435e-3,  1.0e-7),  // I3d1A_fT_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.2348e-3,  1.0e-7),  // I3d1A_fT_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(-1.9442e-1, 1.0e-3),   // I3d1B_fT_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.3260e+1, 1.0e-1),   // I3d1B_fT_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 1.9659e-3,  1.0e-7),  // I3d1C_fT_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 8.2398e-4,  1.0e-8),  // I3d1C_fT_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair(-6.9442e-1,  1.0e-5),  // I2_fT_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.9091e-1,  1.0e-5),  // I2_fT_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.6419   ,  1.0e-4),  // I2_fT_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5262   ,  1.0e-4),  // I2_fT_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 1.5397e+1,  1.0e-3),  // I3_fT_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 6.4504   ,  1.0e-4),  // I3_fT_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.0753e+1,  1.0e-3),  // I3_fT_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.3841e+1,  1.0e-3),  // I3_fT_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-8.0985e+1,  1.0e-3),  // I3d1A_fT_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.3950e+1,  1.0e-3),  // I3d1A_fT_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-6.5754e+1,  1.0e-3),  // I3d1A_fT_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.7570e+1,  1.0e-3),  // I3d1A_fT_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair( 7.4575e2,  1.0e-2),   // I3d1B_fT_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.9223e4,  1.0e-0),   // I3d1B_fT_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair(-7.8647e-3,  1.0e-1),  // I3d1C_fT_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-3.2964e-3,  1.0e-1),  // I3d1C_fT_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_p(-5.0), 0.286545, 1.0e-4);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_p( 0.0), 0.282318, 1.0e-4);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_p( 5.0), 0.277101, 1.0e-4);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_pm(-5.0), 0.292937, 1.0e-4);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_pm( 0.0), 0.284811, 1.0e-4);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_pm( 5.0), 0.272936, 1.0e-4);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_t(-5.0), 0.279865, 1.0e-4);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_t( 0.0), 0.276469, 1.0e-4);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_f_t( 5.0), 0.272724, 1.0e-4);
            }


            /* B -> pi form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.3174;
                p["B::lambda_H^2"]            = 1.2696;
                p["mass::d(2GeV)"]            = 0.0048;
                p["mass::u(2GeV)"]            = 0.0032;
                p["mass::B_d"]                = 5.2795;
                p["mass::pi^+"]               = 0.13957;
                p["decay-constant::B_d"]      = 0.180;
                p["decay-constant::pi"]       = 0.1302;
                p["B->pi::mu@B-LCSR"]         = 1.0;
                p["B->pi::s_0^+,0@B-LCSR"]    = 0.7;
                p["B->pi::s_0^+,1@B-LCSR"]    = 0.0;
                p["B->pi::s_0^+/-,0@B-LCSR"]  = 0.7;
                p["B->pi::s_0^+/-,1@B-LCSR"]  = 0.0;
                p["B->pi::s_0^T,0@B-LCSR"]    = 0.7;
                p["B->pi::s_0^T,1@B-LCSR"]    = 0.0;
                p["B->pi::M^2@B-LCSR"]        = 1.0;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "zero" }
                };

                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->f_p(-5.0), 0.270388,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p( 0.0), 0.356854,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p(+5.0), 0.494302,  eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_0(-5.0), 0.304492,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0( 0.0), 0.356854,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(+5.0), 0.431392,  eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_t(-5.0), 0.227664,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t( 0.0), 0.301374,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t(+5.0), 0.419634,  eps);

            }


            /* B -> K form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.03;
                p["B::lambda_H^2"]            = 0.06;
                p["mass::B_d"]                = 5.27958;
                p["mass::K_d"]                = 0.497614;
                p["decay-constant::B_d"]      = 0.1905;
                p["decay-constant::K_d"]      = 0.1561;
                p["B->K::mu@B-LCSR"]          = 1.0;
                p["B->K::s_0^+,0@B-LCSR"]     = 1.05;
                p["B->K::s_0^+,1@B-LCSR"]     = 0.0;
                p["B->K::s_0^+/-,0@B-LCSR"]   = 1.05;
                p["B->K::s_0^+/-,1@B-LCSR"]   = 0.0;
                p["B->K::s_0^T,0@B-LCSR"]     = 1.05;
                p["B->K::s_0^T,1@B-LCSR"]     = 0.0;
                p["B->K::M^2@B-LCSR"]         = 1.0;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->K::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->f_p(-5.0), 0.208620,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p( 0.0), 0.267282,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p(+5.0), 0.354006,  eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_0(-5.0), 0.240493,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0( 0.0), 0.267282,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(+5.0), 0.299113,  eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_t(-5.0), 0.196080,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t( 0.0), 0.252352,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t(+5.0), 0.336631,  eps);

            }


            /* B -> D form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                /*charm mass = 1.066273     */
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.03;
                p["B::lambda_H^2"]            = 0.06;
                p["mass::B_d"]                = 5.27958;
                p["mass::D^+"]                = 1.86959;
                p["decay-constant::B_d"]      = 0.1905;
                p["decay-constant::D_d"]      = 0.2127;
                p["B->D::mu@B-LCSR"]          = 2.1213;
                p["B->D::s_0^+,0@B-LCSR"]     = 6.0;
                p["B->D::s_0^+,1@B-LCSR"]     = 0.0;
                p["B->D::s_0^+/-,0@B-LCSR"]   = 6.0;
                p["B->D::s_0^+/-,1@B-LCSR"]   = 0.0;
                p["B->D::s_0^T,0@B-LCSR"]     = 6.0;
                p["B->D::s_0^T,1@B-LCSR"]     = 0.0;
                p["B->D::M^2@B-LCSR"]         = 4.5;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->D::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->f_p(-5.0), 0.628668,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p( 0.0), 0.745726,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p(+5.0), 0.917246,       eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_0(-5.0), 0.692330,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0( 0.0), 0.745726,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(+5.0), 0.810308,       eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_t(-5.0), 0.501645,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t( 0.0), 0.616377, 3.0 * eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t(+5.0), 0.823555, 9.0 * eps);

            }


            /* B_s -> K form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B_s::1/lambda_B_p"]          = 1.69348;
                p["B_s::lambda_E^2"]            = 0.03;
                p["B_s::lambda_H^2"]            = 0.06;
                p["mass::B_s"]                  = 5.36677;
                p["mass::K_u"]                  = 0.493677;
                p["decay-constant::B_s"]        = 0.2307;
                p["decay-constant::K_u"]        = 0.1561;
                p["B_s->K::mu@B-LCSR"]          = 1.0;
                p["B_s->K::s_0^+,0@B-LCSR"]     = 1.05;
                p["B_s->K::s_0^+,1@B-LCSR"]     = 0.0;
                p["B_s->K::s_0^+/-,0@B-LCSR"]   = 1.05;
                p["B_s->K::s_0^+/-,1@B-LCSR"]   = 0.0;
                p["B_s->K::s_0^T,0@B-LCSR"]     = 1.05;
                p["B_s->K::s_0^T,1@B-LCSR"]     = 0.0;
                p["B_s->K::M^2@B-LCSR"]         = 1.0;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B_s->K::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->f_p(-5.0), 0.189587,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p( 0.0), 0.239226,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p(+5.0), 0.308755,  eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_0(-5.0), 0.219599,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0( 0.0), 0.239226,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(+5.0), 0.259143,  eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_t(-5.0), 0.184833,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t( 0.0), 0.234306,  eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t(+5.0), 0.304495,  eps);

            }


            /* B_s -> D_s form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                /*charm mass = 1.066273     */
                p["B_s::1/lambda_B_p"]            = 1.69348;
                p["B_s::lambda_E^2"]              = 0.03;
                p["B_s::lambda_H^2"]              = 0.06;
                p["mass::B_s"]                    = 5.36677;
                p["mass::D_s"]                    = 1.96828;
                p["decay-constant::B_s"]          = 0.2307;
                p["decay-constant::D_s"]          = 0.2499;
                p["B_s->D_s::mu@B-LCSR"]          = 2.1213;
                p["B_s->D_s::s_0^+,0@B-LCSR"]     = 6.0;
                p["B_s->D_s::s_0^+,1@B-LCSR"]     = 0.0;
                p["B_s->D_s::s_0^+/-,0@B-LCSR"]   = 6.0;
                p["B_s->D_s::s_0^+/-,1@B-LCSR"]   = 0.0;
                p["B_s->D_s::s_0^T,0@B-LCSR"]     = 6.0;
                p["B_s->D_s::s_0^T,1@B-LCSR"]     = 0.0;
                p["B_s->D_s::M^2@B-LCSR"]         = 4.5;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B_s->D_s::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->f_p(-5.0), 0.539744,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p( 0.0), 0.642184,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_p(+5.0), 0.787744,       eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_0(-5.0), 0.600434,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0( 0.0), 0.642184,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_0(+5.0), 0.688616,       eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_t(-5.0), 0.518002,       eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t( 0.0), 0.643604, 2.0 * eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_t(+5.0), 0.862329, 5.0 * eps);

            }
        }
} kmo2006_form_factors_test;
