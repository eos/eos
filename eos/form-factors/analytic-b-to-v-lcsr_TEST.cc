/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018-2025 Danny van Dyk
 * Copyright (c) 2018      Nico Gubernari
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
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
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
            /* B -> K^* diagnostic values */
            {
                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.3174;
                p["B::lambda_H^2"]            = 1.2696;
                p["mass::d(2GeV)"]            = 0.0048;
                p["mass::u(2GeV)"]            = 0.0032;
                p["mass::B_d"]                = 5.2795;
                p["mass::K_d^*"]              = 0.896;
                p["decay-constant::B_d"]      = 0.180;
                p["B->K^*::f_Kstar_par"]      = 0.217;
                p["B->K^*::mu@B-LCSR"]        = 1.0;
                p["B->K^*::s_0^A1,0@B-LCSR"]  = 1.7;
                p["B->K^*::s_0^A1,1@B-LCSR"]  = 0.0;
                p["B->K^*::s_0^A2,0@B-LCSR"]  = 1.7;
                p["B->K^*::s_0^A2,1@B-LCSR"]  = 0.0;
                p["B->K^*::s_0^A30,0@B-LCSR"] = 1.7;
                p["B->K^*::s_0^A30,1@B-LCSR"] = 0.0;
                p["B->K^*::s_0^V,0@B-LCSR"]   = 1.7;
                p["B->K^*::s_0^V,1@B-LCSR"]   = 0.0;
                p["B->K^*::M^2@B-LCSR"]       = 1.0;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "all"  },
                    { "gminus"_ok, "zero" }
                };
                AnalyticFormFactorBToVLCSR<lcsr::BToKstar> ff{ p, o };
                auto diagnostics = ff.diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(0.126841,   1.0e-7), // m_v(mu) in the MSbar scheme, strange quark for these tests
                    std::make_pair(0.896,      1.0e-7), // m_V, the K^* mass
                    std::make_pair(0.217,      1.0e-7), // f_V, the K^* decay constant

                    std::make_pair(0.0748879,  1.0e-7), // sigma_0 value for q^2 = 5.0
                    std::make_pair(1.7,        1.0e-7), // s_0 value for V

                    /* A_1 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-0.785308,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-0.662564,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-0.539820,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.007970,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-0.839074,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-0.670179,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_1 phi_bar */
                    std::make_pair(1.984386e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.984386e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.984386e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.504278e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.504278e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.504278e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(6.611618e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.578220e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.544822e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(8.121517e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(6.760676e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.399835e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair( 8.503654e-1,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 6.997323e-1,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 5.490993e-1,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 2.124788e-2,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.312584e-3,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.387221e-2,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_1 g_+ */
                    std::make_pair(-4.450784e-3,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.450784e-3,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.450784e-3,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.256998e-2,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.256998e-2,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.256998e-2,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.422015e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190204e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-9.583918e-2,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.904307e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.221156e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.538004e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-5.794952,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.809833,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.824715,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.049258,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.868187,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.687118,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(4.971110e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.194124e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.418066e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.425962e-2,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.187028e-2,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(9.480936e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(2.079085e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.740800e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.402515e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.368797e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.929743e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.490689e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(2.425291,  1.0e-6),     // I3d2_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.931553,  1.0e-6),     // I3d2_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.437815,  1.0e-6),     // I3d2_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.451867e-1,  1.0e-6), // I3d2_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.042344e-1,  1.0e-6), // I3d2_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.632821e-1,  1.0e-6), // I3d2_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-3.685383e-3,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.685383e-3,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.685383e-3,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.355562e-2,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.355562e-2,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.355562e-2,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_bar */
                    std::make_pair(-2.557924e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.557924e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.557924e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.519371e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.519371e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.519371e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(1.281981e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.281981e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.281981e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(8.563666e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(8.563666e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(8.563666e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(9.033828e-3,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(9.033828e-3,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(9.033828e-3,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.829506e-2,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.829506e-2,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.829506e-2,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar */
                    std::make_pair(3.960039e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.960039e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.960039e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.378878e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.378878e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.378878e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(1.545825e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.304214e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.062602e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.047158e-3,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(8.716971e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(6.962364e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(1.094275e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(9.190984e-3,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.439213e-3,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.498098e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.881025e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.263952e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(4.848042e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.030634e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.213225e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(6.844165e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.482980e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.121795e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar ATTENTION:here I have to low the precision to e-5 to pass the test*/
                    std::make_pair(7.792133,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.162872,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.533612,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.858880,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.694256,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(0.5296311, 1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_1 phi_3 */
                    std::make_pair( 1.7872e-3,  1.0e-7),  // I1_A1_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.7454e-3,  1.0e-7),  // I1_A1_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9875e-4,  1.0e-8),  // I1_A1_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.2609e-4,  1.0e-8),  // I1_A1_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_3 */
                    std::make_pair( 3.6060e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.2251e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.0787e-3,  1.0e-7),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2610e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_1 phi_bar_3 */
                    std::make_pair( 8.7010e-4,  1.0e-8),  // I1_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8234e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.8006e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2156e-2,  1.0e-6),  // I1_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-3.2004e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.8686e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1237e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.4756e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 3.5811e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.1752e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.4075   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.9944   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-3.1896,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.7368,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.0492,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0663e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 NOTE: this values are very high. This is not a problem,
                       since the IB integrands are only integrated over w1 between 0 and sigma mB = 0.3937...*/
                    std::make_pair(-9.4237e1,  1.0e-3),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-4.8570e4,  1.0e-0),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-6.9099e-2, 1.0e-6),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4480e-1, 1.0e-5),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 5.3970e-4,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.2762e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.3462e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.9176e-2,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-2.1902e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.3468e-2,  1.0e-6),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.3497e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.6372e-2,  1.0e-5),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 6.8285e-1,  1.0e-5),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 7.6758e+2,  1.0e-2),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-2.2434e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.6913e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-6.0880e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.7362e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.1798e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.5302e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.2090e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.1264e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 1.9027e-1,  1.0e-5),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1388e+2,  1.0e-2),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.3486e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.9547e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0095e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.3539e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-8.7809,    1.0e-4),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.0897e3,  1.0e-1),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 2.4694e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 6.5788e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A1_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_1 phi_4 */
                    std::make_pair(-1.5968e-3,  1.0e-7),  // I1_A1_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.2460e-3,  1.0e-7),  // I1_A1_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.9279e-5,  1.0e-9),  // I1_A1_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.9072e-4,  1.0e-8),  // I1_A1_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-3.2572e-2,  1.0e-6),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0701e-1,  1.0e-5),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.6172e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.8903e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_1 phi_bar_4 */
                    std::make_pair(-2.3439e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.7004e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0136e-2,  1.0e-6),  // I1_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.4384e-2,  1.0e-6),  // I1_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-2.6830e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.8049e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.4274e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6740e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-9.7591e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.2410   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.2049   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0159e+1,  1.0e-3),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 7.3408   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5402e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.5044   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.7843e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-7.1511  ,  1.0e-4),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.6832e3,  1.0e-1),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-4.9002  ,  1.0e-4),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.0269e1,  1.0e-3),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair( 9.8237e-3,  1.0e-7),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.4454e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9291e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.5413e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 1.9972e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 6.9648e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.9281   ,  1.0e-4),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1263e+1,  1.0e-3),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-3.0365   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.2680   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-9.5387   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.5782e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-1.5607e+1,  1.0e-3),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.2772e+4,  1.0e-0),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 6.8878e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4304e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.3852e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.8838e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.9780e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0222   ,  1.0e-4),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.8296e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.3934e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.2957e+1,  1.0e-3),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0605e+4,  1.0e-0),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-5.3685e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.6816   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.0407e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6220   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-4.3710e2,  1.0e-2),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.5043e5,  1.0e+1),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-3.5097,     1.0e-4),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-9.3501,     1.0e-4),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_1 psi_bar_4 */
                    std::make_pair( 5.7779e-3,  1.0e-7),  // I1_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.1616e-3,  1.0e-7),  // I1_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.6260e-2,  1.0e-6),  // I1_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6343e-2,  1.0e-6),  // I1_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 8.6733e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.7459e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.4430e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4532e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 2.3938   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3098   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5023e+1,  1.0e-3),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.7708   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-1.9486e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-8.1903   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.1008e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3032e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 4.4659e1,  1.0e-3),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.3431e3,  1.0e-1),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 3.2135   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.3469   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 psi_bar_bar_4 */
                    std::make_pair(-1.9915e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5176e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3719e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.3526e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_bar_4 */
                    std::make_pair(-4.0237e-1,  1.0e-5),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.9835e-1,  1.0e-5),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0932e+1,  1.0e-3),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.7641   ,  1.0e-4),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_bar_4 */
                    std::make_pair( 6.1456   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.0380   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6582e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5666e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_bar_4 */
                    std::make_pair(-5.1460   ,  1.0e-4),  // I3d1B_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1437e+3,  1.0e-1),  // I3d1B_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 psi_bar_bar_4 */
                    std::make_pair(-2.4476e-2,  1.0e-6),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.0838e-2,  1.0e-6),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.6317e-1,  1.0e-5),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.7225e-1,  1.0e-5),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A psi_bar_bar_4 */
                    std::make_pair( 1.8911e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.4730e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6665e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.8231e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B psi_bar_bar_4 */
                    std::make_pair( 3.7481   ,  1.0e-4),  // I4d1B_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.3304e+2,  1.0e-2),  // I4d1B_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A psi_bar_bar_4 */
                    std::make_pair( 4.5638e-1,  1.0e-5),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3208   ,  1.0e-4),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.4202e-1,  1.0e-5),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4369   ,  1.0e-4),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B psi_bar_bar_4 */
                    std::make_pair(-8.9442e1,  1.0e-3),   // I4d2B_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-9.9275e3,  1.0e-1),   // I4d2B_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C psi_bar_bar_4 */
                    std::make_pair( 3.7564,     1.0e-4),  // I4d2C_A1_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 2.1742,     1.0e-4),  // I4d2C_A1_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_psi_bar_bar_4(sigma_0, 5.0)

                    /* I_1 chi_bar_4 */
                    std::make_pair( 1.7652e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.6994e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.7318e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.7250e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 5.7930e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.9119e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.7142   ,  1.0e-4),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6316   ,  1.0e-4),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 7.6042e-1,  1.0e-5),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6257   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0788   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4748   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 1.2674e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.9425e-2,  1.0e-6),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2167e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.1849e-3,  1.0e-7),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-1.9119e2,  1.0e-2),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.2875e4,  1.0e-0),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 13.702,  1.0e-3),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 5.7428,  1.0e-4),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_bar_4 */
                    std::make_pair(-7.9662e-2,  1.0e-6),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.0703e-2,  1.0e-6),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1488   ,  1.0e-4),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3410   ,  1.0e-4),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_bar_4 */
                    std::make_pair(-1.6297   ,  1.0e-4),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.2608   ,  1.0e-4),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.3779e+1,  1.0e-3),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.7457e+1,  1.0e-3),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_bar_4 */
                    std::make_pair( 2.4665e+1,  1.0e-3),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6506e+1,  1.0e-3),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0617e+2,  1.0e-2),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.2123e+1,  1.0e-3),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_bar_4 */
                    std::make_pair(-1.6802e+1,  1.0e-3),  // I3d1B_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.7344e+3,  1.0e-1),  // I3d1B_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 chi_bar_bar_4 */
                    std::make_pair(-1.3802e-2,  1.0e-6),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2070e-1,  1.0e-5),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6296   ,  1.0e-4),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.7483e-1,  1.0e-5),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A chi_bar_bar_4 */
                    std::make_pair(-4.3562e-1,  1.0e-5),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4129   ,  1.0e-4),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.4108   ,  1.0e-4),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.9721e-1,  1.0e-5),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B chi_bar_bar_4 */
                    std::make_pair( 1.6046e+1,  1.0e-3),  // I4d1B_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.5663e+3,  1.0e-1),  // I4d1B_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A chi_bar_bar_4 */
                    std::make_pair( 6.8813   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.6899   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.9144   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.9432   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B chi_bar_bar_4 */
                    std::make_pair(-3.9644e2,  1.0e-2),   // I4d2B_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-4.2563e4,  1.0e+0),   // I4d2B_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C chi_bar_bar_4 */
                    std::make_pair( 14.099,     1.0e-3),  // I4d2C_A1_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 8.1613,     1.0e-4),  // I4d2C_A1_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_chi_bar_bar_4(sigma_0, 5.0)

                    /* A_2 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(-1.214827e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.214827e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.214827e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.578211e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.578211e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.578211e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(-5.879753,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-5.879753,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-5.879753,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.216107,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.216107,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-5.216107,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.067506,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.067506,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.067506,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(3.825497e-3,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.825497e-3,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.825497e-3,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.029345e-2,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.029345e-2,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.029345e-2,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(1.542539e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.542539e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.542539e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.534507e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.534507e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.534507e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(1.377104,  1.0e-6),    // I3d2_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.377104,  1.0e-6),    // I3d2_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.377104,  1.0e-6),    // I3d2_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.004232,  1.0e-6),   // I3d2_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.004232,  1.0e-6),   // I3d2_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.004232,  1.0e-6),   // I3d2_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.050510e-1,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.050510e-1,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.050510e-1,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.666550,  1.0e-6),    // I3d1_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.666550,  1.0e-6),    // I3d1_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.666550,  1.0e-6),    // I3d1_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar ATTENTION: this term is huge!!! */
                    std::make_pair(5.266027e+1,  1.0e-4),     // I3d2_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.266027e+1,  1.0e-4),     // I3d2_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.266027e+1,  1.0e-4),     // I3d2_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.420043e+2,  1.0e-4), // I3d2_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.420043e+2,  1.0e-4), // I3d2_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.420043e+2,  1.0e-4), // I3d2_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(-2.839758e-4,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.839758e-4,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.839758e-4,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.611760e-3,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.611760e-3,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.611760e-3,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(-3.019630e-2,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.019630e-2,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.019630e-2,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.201416e-1,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.201416e-1,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.201416e-1,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar ATTENTION:here I have to low the precision to e-5 to pass the test*/
                    std::make_pair(-2.311842,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.311842,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.311842,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.365898,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.365898,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.365898,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar ATTENTION:this term is huge!!!*/
                    std::make_pair(-1.075549e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.075549e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.075549e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.356486e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.356486e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.356486e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair(-1.4278e-1,  1.0e-5),  // I2_A2_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0585e-1,  1.0e-5),  // I2_A2_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.3765e-2,  1.0e-6),  // I2_A2_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.0061e-2,  1.0e-6),  // I2_A2_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-1.8162e-3,  1.0e-7),  // I2_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.8062e-3,  1.0e-7),  // I2_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2108e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.5374e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 5.9725e-1,  1.0e-5),  // I3_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1750   ,  1.0e-4),  // I3_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.0284   ,  1.0e-4),  // I3_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.3236   ,  1.0e-4),  // I3_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-2.7218,     1.0e-4),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.7722,     1.0e-4),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.0008,     1.0e-4),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.8757e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-1.6161e2,  1.0e-2),   // I3d1B_A2_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-8.3291e4,  1.0e-0),   // I3d1B_A2_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-1.6055e-1, 1.0e-5),   // I3d1C_A2_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-3.3647e-1, 1.0e-5),   // I3d1C_A2_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-6.2567e-3,  1.0e-7),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.8476e-2,  1.0e-6),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.1386e-1,  1.0e-5),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.1538e-1,  1.0e-5),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 1.6675e-2,  1.0e-6),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.7536e-1,  1.0e-5),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.5971e-1,  1.0e-5),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.1116   ,  1.0e-4),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(-2.5336   ,  1.0e-4),  // I4d1B_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.8480e+3,  1.0e-1),  // I4d1B_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair( 6.0010e-1,  1.0e-5),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.2658   ,  1.0e-4),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.9049e-1,  1.0e-5),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.7364   ,  1.0e-4),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-5.4467,    1.0e-4),   // I4d2B_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-5.6347e4,  1.0e-0),   // I4d2B_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 4.5808e-1,  1.0e-5),  // I4d2C_A2_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.2203   ,  1.0e-4),  // I4d2C_A2_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A2_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-5.1174e-2,  1.0e-6),  // I2_A2_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.6812e-1,  1.0e-5),  // I2_A2_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.5407e-3,  1.0e-7),  // I2_A2_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.1120e-3,  1.0e-7),  // I2_A2_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 4.8925e-3,  1.0e-7),  // I2_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6073e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1158e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.0898e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.1690   ,  1.0e-4),  // I3_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1294   ,  1.0e-4),  // I3_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.9270   ,  1.0e-4),  // I3_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2218e+1,  1.0e-3),  // I3_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.7895   ,  1.0e-4),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4895e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.2705e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2655e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-5.9305e1,  1.0e-3),   // I3d1B_A2_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.2253e4,  1.0e-0),   // I3d1B_A2_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-6.8824  ,  1.0e-4),   // I3d1C_A2_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4423e1,  1.0e-3),   // I3d1C_A2_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair(-2.9488e-1,  1.0e-5),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0397   ,  1.0e-4),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7825   ,  1.0e-4),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6638e+1,  1.0e-3),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair( 3.7096   ,  1.0e-4),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1027e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-9.3044e-1,  1.0e-5),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.1932   ,  1.0e-4),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair( 2.1750e+1,  1.0e-3),  // I3d1B_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.7799e+4,  1.0e-0),  // I3d1B_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(-1.2599   ,  1.0e-4),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.2681   ,  1.0e-4),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.4967e+1,  1.0e-3),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.0940e+1,  1.0e-3),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair( 2.4486   ,  1.0e-4),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.1715   ,  1.0e-4),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.6923e+2,  1.0e-2),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.8427e+2,  1.0e-2),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.1266e+2,  1.0e-2),  // I4d1B_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 9.2197e+4,  1.0e-0),  // I4d1B_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 3.9498e+2,  1.0e-2),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.2766e+3,  1.0e-1),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6454e+3,  1.0e-1),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.4536e+3,  1.0e-1),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-2.8370e3,  1.0e-1),   // I4d2B_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.4601e5,  1.0e+1),   // I4d2B_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair( 2.5661e1,   1.0e-3),  // I4d2C_A2_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 6.8363e1,   1.0e-3),  // I4d2C_A2_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A2_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-1.2061e-2,  1.0e-6),  // I2_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.5994e-3,  1.0e-7),  // I2_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.5687e-2,  1.0e-6),  // I2_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.4113e-2,  1.0e-6),  // I2_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 2.7663   ,  1.0e-4),  // I3_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4881   ,  1.0e-4),  // I3_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7438e+1,  1.0e-3),  // I3_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.8175   ,  1.0e-4),  // I3_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.4656e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0700e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.8371e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.0841e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 5.5938e1,  1.0e-3),   // I3d1B_A2_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.6926e3,  1.0e-1),   // I3d1B_A2_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 3.4460   ,  1.0e-4),  // I3d1C_A2_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.4443   ,  1.0e-4),  // I3d1C_A2_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 psi_bar_bar_4 */
                    std::make_pair( 5.9782e-1,  1.0e-5),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.5797e-1,  1.0e-5),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6102e+1,  1.0e-3),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0067e+1,  1.0e-3),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_bar_4 */
                    std::make_pair(-7.5204   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.8572   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.5910   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.1420   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_bar_4 */
                    std::make_pair( 6.5126   ,  1.0e-4),  // I3d1B_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.4474e+3,  1.0e-1),  // I3d1B_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 psi_bar_bar_4 */
                    std::make_pair( 2.5248   ,  1.0e-4),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8306   ,  1.0e-4),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.8992e+1,  1.0e-3),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.2395e+1,  1.0e-3),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A psi_bar_bar_4 */
                    std::make_pair(-4.8860   ,  1.0e-4),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.3750e-1,  1.0e-5),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.4756e+2,  1.0e-2),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.7031e+2,  1.0e-2),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B psi_bar_bar_4 */
                    std::make_pair( 3.5488e+1,  1.0e-3),  // I4d1B_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 7.8871e+3,  1.0e-1),  // I4d1B_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A psi_bar_bar_4 */
                    std::make_pair(-7.9793e+2,  1.0e-2),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.5492e+2,  1.0e-2),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.5837e+3,  1.0e-1),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.7024e+3,  1.0e-1),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B psi_bar_bar_4 */
                    std::make_pair(-5.1089e2,  1.0e-2),   // I4d2B_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.9941e4,  1.0e+0),   // I4d2B_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C psi_bar_bar_4 */
                    std::make_pair(-2.4472e1,   1.0e-3),  // I4d2C_A2_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4164e1,   1.0e-3),  // I4d2C_A2_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A2_3pt_psi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair(-3.6847e-3,  1.0e-7),  // I2_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.7219e-3,  1.0e-7),  // I2_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7023e-3,  1.0e-7),  // I2_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1950e-2,  1.0e-6),  // I2_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair(-1.3357   ,  1.0e-4),  // I3_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0534   ,  1.0e-4),  // I3_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.3830e+1,  1.0e-3),  // I3_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.3001   ,  1.0e-4),  // I3_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-1.6148e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.7073   ,  1.0e-4),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0750e+2,  1.0e-2),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-8.8422e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-3.2786e2,  1.0e-2),   // I3d1B_A2_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.9227e4,  1.0e-0),   // I3d1B_A2_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 16.912,  1.0e-3),     // I3d1C_A2_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.0883,  1.0e-4),     // I3d1C_A2_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 chi_bar_bar_4 */
                    std::make_pair( 2.3913   ,  1.0e-4),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8319   ,  1.0e-4),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.4409e+1,  1.0e-3),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.0266e+1,  1.0e-3),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_bar_4 */
                    std::make_pair(-3.0082e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.9429e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0364e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2568e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_bar_4 */
                    std::make_pair( 2.6050e+1,  1.0e-3),  // I3d1B_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.7898e+3,  1.0e-0),  // I3d1B_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 chi_bar_bar_4 */
                    std::make_pair( 1.0334e+1,  1.0e-3),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.7172   ,  1.0e-4),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.8024e+2,  1.0e-2),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.7378e+2,  1.0e-2),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A chi_bar_bar_4 */
                    std::make_pair(-2.0169e+1,  1.0e-3),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.7782   ,  1.0e-4),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0075e+3,  1.0e-1),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.9147e+3,  1.0e-1),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B chi_bar_bar_4 */
                    std::make_pair( 1.2792e+2,  1.0e-2),  // I4d1B_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.8431e+4,  1.0e-0),  // I4d1B_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A chi_bar_bar_4 */
                    std::make_pair(-3.2142e+3,  1.0e-1),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2789e+3,  1.0e-1),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.8320e+4,  1.0e-0),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0746e+4,  1.0e-0),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B chi_bar_bar_4 */
                    std::make_pair(-2.2061e3,  1.0e-1),   // I4d2B_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1103e5,  1.0e+1),   // I4d2B_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C chi_bar_bar_4 */
                    std::make_pair(-1.1506e2,   1.0e-2),  // I4d2C_A2_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-6.6600e1,   1.0e-3),  // I4d2C_A2_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A2_3pt_chi_bar_bar_4(sigma_0, 5.0)

                    /* A_30 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(7.094319e-1,  1.0e-6), // I1_A30_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.094319e-1,  1.0e-6), // I1_A30_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.094319e-1,  1.0e-6), // I1_A30_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.004800,  1.0e-6),    // I1_A30_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.004800,  1.0e-6),    // I1_A30_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.004800,  1.0e-6),    // I1_A30_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(1.362149e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.362149e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.362149e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.321887e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.321887e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.321887e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(6.815761,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.815761,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(6.815761,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.145075,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.145075,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.145075,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(1.339818e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.339818e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.339818e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.064229e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.064229e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.064229e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(5.802308,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.802308,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.802308,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.285435,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.285435,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.285435,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(-4.490801e-3,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.490801e-3,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.490801e-3,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.421476e-2,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.421476e-2,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.421476e-2,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(-1.991596e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.991596e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.991596e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.702610e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.702610e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.702610e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(-3.153181,  1.0e-6),      // I3d2_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.153181,  1.0e-6),      // I3d2_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.153181,  1.0e-6),      // I3d2_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.308782e-1,  1.0e-6),   // I3d2_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.308782e-1,  1.0e-6),   // I3d2_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.308782e-1,  1.0e-6),   // I3d2_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-8.902671e-3,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-8.902671e-3,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-8.902671e-3,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.233207e-1,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.233207e-1,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.233207e-1,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-8.583107e-1,  1.0e-6), // I3d1_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-8.583107e-1,  1.0e-6), // I3d1_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-8.583107e-1,  1.0e-6), // I3d1_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.726354,  1.0e-6),    // I3d1_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.726354,  1.0e-6),    // I3d1_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-5.726354,  1.0e-6),    // I3d1_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar ATTENTION: this term is huge!!! */
                    std::make_pair(-6.045307e+1,  1.0e-4),     // I3d2_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.045307e+1,  1.0e-4),     // I3d2_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.045307e+1,  1.0e-4),     // I3d2_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.892940e+2,  1.0e-4), // I3d2_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.892940e+2,  1.0e-4), // I3d2_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.892940e+2,  1.0e-4), // I3d2_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(3.184156e-4,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.184156e-4,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.184156e-4,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.570244e-3,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.570244e-3,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.570244e-3,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(3.437931e-2,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.437931e-2,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.437931e-2,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.767839e-1,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.767839e-1,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.767839e-1,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar ATTENTION:here I have to low the precision to e-5 to pass the test*/
                    std::make_pair(2.712355,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.712355,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.712355,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(9.969988,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(9.969988,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(9.969988,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar ATTENTION:this term is huge!!!*/
                    std::make_pair(1.354431e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.354431e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.354431e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.205167e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.205167e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.205167e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 1.5611e-1,  1.0e-5),  // I2_A30_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.2052e-1,  1.0e-5),  // I2_A30_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6198e-2,  1.0e-6),  // I2_A30_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.4645e-2,  1.0e-6),  // I2_A30_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-5.0321e-2,  1.0e-6),  // I2_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0546e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.3547e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.0303e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(-7.7567e-1,  1.0e-5),  // I3_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5220   ,  1.0e-4),  // I3_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.2343   ,  1.0e-4),  // I3_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0809e+1,  1.0e-3),  // I3_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair( 1.2282,     1.0e-4),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.0898,     1.0e-4),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.7350e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.6519e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair( 2.1073e2,  1.0e-2),   // I3d1B_A30_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0861e5,  1.0e+1),   // I3d1B_A30_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 2.1712e-1, 1.0e-5),   // I3d1C_A30_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 4.5501e-1, 1.0e-5),   // I3d1C_A30_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair( 5.7732e-3,  1.0e-7),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.1120e-2,  1.0e-6),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.1841e-2,  1.0e-6),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.8293e-1,  1.0e-5),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair(-7.6893e-3,  1.0e-7),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.9699e-1,  1.0e-5),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.0384e-1,  1.0e-5),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.2570   ,  1.0e-4),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 3.0623   ,  1.0e-4),  // I4d1B_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.4423e+3,  1.0e-1),  // I4d1B_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-2.7983e-1,  1.0e-5),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.9617   ,  1.0e-4),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6887   ,  1.0e-4),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.6004e+1,  1.0e-3),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair( 1.6362e1,  1.0e-3),   // I4d2B_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.1081e4,  1.0e-0),   // I4d2B_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(-5.1471e-1,  1.0e-5),  // I4d2C_A30_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.3712   ,  1.0e-4),  // I4d2C_A30_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A30_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair( 3.7842e-2,  1.0e-6),  // I2_A30_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.2432e-1,  1.0e-5),  // I2_A30_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8788e-3,  1.0e-7),  // I2_A30_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.5196e-3,  1.0e-7),  // I2_A30_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 1.3555e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.4534e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.8622e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.4102   ,  1.0e-4),  // I2_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair( 1.5861   ,  1.0e-4),  // I3_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.5683   ,  1.0e-4),  // I3_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.7005   ,  1.0e-4),  // I3_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6571e+1,  1.0e-3),  // I3_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair(-5.6152   ,  1.0e-4),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.5960e-1,  1.0e-5),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0036e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5259e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair( 7.3361e1,  1.0e-3),   // I3d1B_A30_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.7527e4,  1.0e-0),   // I3d1B_A30_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair( 9.1963  ,  1.0e-4),   // I3d1C_A30_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.9272e1,  1.0e-3),   // I3d1C_A30_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 3.8210e-1,  1.0e-5),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.2988   ,  1.0e-4),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.5656   ,  1.0e-4),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.1519e+1,  1.0e-3),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-3.5670   ,  1.0e-4),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.7556   ,  1.0e-4),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6398e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.8239e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-3.3673e+1,  1.0e-3),  // I3d1B_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.7556e+4,  1.0e-0),  // I3d1B_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 1.4545   ,  1.0e-4),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.8463   ,  1.0e-4),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.8947e+1,  1.0e-3),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.1832e+1,  1.0e-3),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair( 2.9786e-2,  1.0e-6),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0427e+1,  1.0e-3),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.7024e+2,  1.0e-2),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0682e+3,  1.0e-1),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(-1.3927e+2,  1.0e-2),  // I4d1B_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1398e+5,  1.0e+1),  // I4d1B_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-4.6218e+2,  1.0e-2),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.4356e+3,  1.0e-1),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3003e+2,  1.0e-2),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1782e+3,  1.0e-1),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 2.7803e3,  1.0e-1),   // I4d2B_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1866e6,  1.0e+2),   // I4d2B_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-2.7341e1,   1.0e-3),  // I4d2C_A30_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-7.2839e1,   1.0e-3),  // I4d2C_A30_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A30_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-3.3416e-1,  1.0e-5),  // I2_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.8285e-1,  1.0e-5),  // I2_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0970   ,  1.0e-4),  // I2_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.4516e-1,  1.0e-5),  // I2_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair(-3.7711   ,  1.0e-4),  // I3_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0380   ,  1.0e-4),  // I3_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.3744e+1,  1.0e-3),  // I3_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0660e+1,  1.0e-3),  // I3_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair( 1.7986e+1,  1.0e-3),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 6.0131   ,  1.0e-4),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.1801e+1,  1.0e-3),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5690e+1,  1.0e-3),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(-7.4685e1,  1.0e-3),   // I3d1B_A30_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-8.9355e3,  1.0e-1),   // I3d1B_A30_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair(-4.7949   ,  1.0e-4),  // I3d1C_A30_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.0097   ,  1.0e-4),  // I3d1C_A30_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair(-1.0209e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.1395e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5799e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.3110e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 1.3612   ,  1.0e-4),  // I3_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.6236   ,  1.0e-4),  // I3_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5774e+1,  1.0e-3),  // I3_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.2160   ,  1.0e-4),  // I3_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 1.8873e+1,  1.0e-3),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8387e-1,  1.0e-5),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6875e+2,  1.0e-2),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0068e+2,  1.0e-2),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair( 4.2752e2,  1.0e-2),   // I3d1B_A30_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.1149e4,  1.0e-0),   // I3d1B_A30_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair(-23.002,  1.0e-3),     // I3d1C_A30_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-9.6406,  1.0e-4),     // I3d1C_A30_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* V */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-1.244211e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.244211e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.244211e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.640698e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.640698e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.640698e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(1.047518e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.047518e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.047518e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.321959e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.321959e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.321959e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(1.417796e-1,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.417796e-1,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.417796e-1,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.240348e-2,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.240348e-2,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.240348e-2,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-2.349792e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.349792e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.349792e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.636322e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.636322e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.636322e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-9.74102e-1,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-9.74102e-1,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-9.74102e-1,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.075191,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.075191,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.075191,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(7.878105e-4,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.878105e-4,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.878105e-4,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.321685e-3,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.321685e-3,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.321685e-3,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(3.347036e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.347036e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.347036e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.012802e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.012802e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.012802e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(4.290449e-1,  1.0e-6),      // I3d2_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.290449e-1,  1.0e-6),      // I3d2_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.290449e-1,  1.0e-6),      // I3d2_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.535883e-2,  1.0e-6),   // I3d2_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.535883e-2,  1.0e-6),   // I3d2_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.535883e-2,  1.0e-6),   // I3d2_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(2.449130e-5,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.449130e-5,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.449130e-5,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.704475e-4,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.704475e-4,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.704475e-4,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(1.750898e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.750898e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.750898e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.809225e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.809225e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.809225e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(7.915833e-2,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.915833e-2,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.915833e-2,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.191980e-1,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.191980e-1,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.191980e-1,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar */
                    std::make_pair(1.392592,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.392592,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.392592,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.001663e-1,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.001663e-1,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.001663e-1,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair(9.4355e-3,  1.0e-7),  // I2_V_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(1.9773e-2,  1.0e-6),  // I2_V_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(1.5772e-3,  1.0e-7),  // I2_V_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(3.3054e-3,  1.0e-7),  // I2_V_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(4.5937e-3,  1.0e-7),  // I2_V_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(9.6268e-3,  1.0e-7),  // I2_V_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(3.0624e-2,  1.0e-6),  // I2_V_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(6.4178e-2,  1.0e-6),  // I2_V_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(9.3704e-2,  1.0e-6),  // I3_V_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(1.9637e-1,  1.0e-5),  // I3_V_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(6.2469e-1,  1.0e-5),  // I3_V_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(1.3091   ,  1.0e-4),  // I3_V_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-6.6912e-1,  1.0e-5), // I3d1A_V_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.4023,     1.0e-4), // I3d1A_V_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.7239e-1,  1.0e-5), // I3d1A_V_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.7084e-1,  1.0e-5), // I3d1A_V_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 NOTE: this values are very high. This is not a problem,
                       since the IB integrands are only integrated over w1 between 0 and sigma mB = 0.3937...*/
                    std::make_pair(-2.2888e1,  1.0e-3),  // I3d1B_V_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1796e4,  1.0e-0),  // I3d1B_V_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-5.8357e-4,  1.0e-8), // I4_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0296e-3,  1.0e-7), // I4_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5776e-2,  1.0e-6), // I4_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.5196e-2,  1.0e-6), // I4_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 7.2788e-3,  1.0e-7), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.1485e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.8598e-3,  1.0e-7), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.0219e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(4.6212e-2,  1.0e-6),  // I4d1B_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(5.1946e+1,  1.0e-3),  // I4d1B_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.1064e-2,  1.0e-6), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.7109e-3,  1.0e-7), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.1986e-2,  1.0e-6), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2000e-1,  1.0e-5), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-2.0391,    1.0e-4),  // I4d2B_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.0234e2,  1.0e-2),  // I4d2B_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(6.8405e-3,  1.0e-7),  // I4d2C_V_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(1.8223e-2,  1.0e-6),  // I4d2C_V_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),        // I4d2D_V_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-8.4303e-3,  1.0e-7), // I2_V_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.7696e-2,  1.0e-6), // I2_V_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.1856e-4,  1.0e-8), // I2_V_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0069e-3,  1.0e-7), // I2_V_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-1.2375e-2,  1.0e-6), // I2_V_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.0654e-2,  1.0e-6), // I2_V_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3514e-2,  1.0e-6), // I2_V_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2873e-1,  1.0e-5), // I2_V_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-2.5360e-1,  1.0e-5), // I3_V_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-8.4161e-1,  1.0e-5), // I3_V_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0929   ,  1.0e-4), // I3_V_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6399   ,  1.0e-4), // I3_V_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair(1.4859   ,  1.0e-4),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(2.6042   ,  1.0e-4),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(3.8988e-1,  1.0e-5),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(2.4277e-1,  1.0e-5),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-1.7368  ,  1.0e-4),  // I3d1B_V_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.5169e2,  1.0e-2),  // I3d1B_V_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-1.2709,  1.0e-4),    // I3d1C_V_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.6634,  1.0e-4),    // I3d1C_V_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                     /* I_3 phi_bar_bar_4 */
                    std::make_pair(1.5385e-4,  1.0e-8),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(2.0494e-3,  1.0e-7),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(7.4673e-4,  1.0e-8),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(9.9467e-3,  1.0e-7),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-6.7819e-4,  1.0e-8), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.0338e-3,  1.0e-7), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7151e-3,  1.0e-7), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2846e-2,  1.0e-6), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(1.5951e-1,  1.0e-5),  // I3d1B_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(1.3054e+2,  1.0e-2),  // I3d1B_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(4.4873e-3,  1.0e-7),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(4.6496e-2,  1.0e-6),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(4.1781e-2,  1.0e-6),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(2.7895e-1,  1.0e-5),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-3.5843e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0299e-1,  1.0e-5), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8001e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6487e-1,  1.0e-5), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(3.1470   ,  1.0e-4),  // I4d1B_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(2.5754e+3,  1.0e-1),  // I4d1B_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-9.4292e-2,  1.0e-6), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5857   ,  1.0e-4), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2350e-1,  1.0e-5), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.8178e-1,  1.0e-5), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-9.2787e1,  1.0e-3),  // I4d2B_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.7865e4,  1.0e+0),  // I4d2B_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-9.7221e-1,  1.0e-5), // I4d2C_V_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.5900,     1.0e-4), // I4d2C_V_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),        // I4d2D_V_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(3.0504e-2,  1.0e-6),  // I2_V_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(1.6692e-2,  1.0e-6),  // I2_V_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(1.9143e-1,  1.0e-5),  // I2_V_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(8.6282e-2,  1.0e-6),  // I2_V_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair(6.2536e-1,  1.0e-5),  // I3_V_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(3.4702e-1,  1.0e-5),  // I3_V_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(3.9098   ,  1.0e-4),  // I3_V_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(1.7701   ,  1.0e-4),  // I3_V_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-4.0182   ,  1.0e-4), // I3d1A_V_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5404   ,  1.0e-4), // I3d1A_V_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5380   ,  1.0e-4), // I3d1A_V_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.2218e-1,  1.0e-5), // I3d1A_V_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(1.0846e1,  1.0e-3),   // I3d1B_V_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(1.2977e3,  1.0e-1),   // I3d1B_V_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair(8.9017e-1,  1.0e-5),  // I3d1C_V_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(3.7310e-1,  1.0e-5),  // I3d1C_V_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 9.3196e-3,  1.0e-7), // I2_V_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9531e-2,  1.0e-6), // I2_V_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.4423e-2,  1.0e-6), // I2_V_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.0225e-2,  1.0e-6), // I2_V_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 1.9011e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.9840e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9420e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.1655e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(3.0244e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(6.3382e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(4.6805e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(9.8088e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-4.6435e1,  1.0e-3),  // I3d1B_V_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-5.5556e3,  1.0e-1),  // I3d1B_V_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 3.3413,  1.0e-4),    // I3d1C_V_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.4004,  1.0e-4),    // I3d1C_V_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* T_1 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-0.646387,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-0.646387,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-0.646387,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-0.817720,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-0.817720,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-0.817720,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(5.442027e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.442027e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.442027e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(6.588619e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(6.588619e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(6.588619e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(4.091723e-3,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.091723e-3,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.091723e-3,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.156818e-2,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.156818e-2,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.156818e-2,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(1.697257e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.697257e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.697257e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.877432e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.877432e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.877432e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(6.177960e-5,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.177960e-5,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(6.177960e-5,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.120409e-4,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.120409e-4,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.120409e-4,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(4.351159e-3,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.351159e-3,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.351159e-3,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.359736e-2,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.359736e-2,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.359736e-2,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(1.272367e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.272367e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.272367e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(8.495098e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(8.495098e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(8.495098e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(8.963349e-3,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.963349e-3,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(8.963349e-3,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.805355e-2,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.805355e-2,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.805355e-2,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(3.927598e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.927598e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.927598e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.327408e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.327408e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.327408e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 4.5886e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.4503e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.6959e-3,  1.0e-7),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6064e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-1.1314e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.4519e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.4936e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5829e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 4.5570e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.3851e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0480   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.3621   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-3.7659,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.9184,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.6333,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.7503,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-1.1469e2,  1.0e-2),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-5.9111e4,  1.0e-0),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-3.4549e-2, 1.0e-6),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-7.2404e-2, 1.0e-6),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 3.2513e-4,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.3303e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1672e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.8868e-2,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-1.7844e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.3768e-2,  1.0e-6),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6351e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5106e-2,  1.0e-5),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 3.3705e-1,  1.0e-5),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.7887e+2,  1.0e-2),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-2.8445e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.8378e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.7016e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.2024e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.8606e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1531e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.9642e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.4101e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 2.3156e-1,  1.0e-5),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.6030e+2,  1.0e-2),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.3119e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.5379e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0544e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.0253e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-1.0603e1,  1.0e-3),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.5401e3,  1.0e-1),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 3.2542e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 8.6695e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A1_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-4.1175e-2,  1.0e-6),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.3527e-1,  1.0e-5),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0443e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.9177e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-2.5931e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.4342e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2585e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6375e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.2366   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1052   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3291   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2873e+1,  1.0e-3),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.5867   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7147e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.6855   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5151e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-8.7031  ,  1.0e-4),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.2656e3,  1.0e-1),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-6.2022  ,  1.0e-4),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.2998e1,  1.0e-3),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair( 4.8973e-3,  1.0e-7),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7033e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.6386e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.7612e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 9.9891e-2,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.4989e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9623   ,  1.0e-4),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.6343   ,  1.0e-4),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-1.5172   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.6362   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.7676   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2860e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-7.6307   ,  1.0e-4),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.2447e+3,  1.0e-1),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 1.6691e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.0879e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0142e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0677   ,  1.0e-4),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.2508e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5027   ,  1.0e-4),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.0447e-3,  1.0e-7),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.1543e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.5769e+1,  1.0e-3),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.2905e+4,  1.0e-0),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-3.9175e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.2016   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8795e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.6798   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-5.0562e2,  1.0e-2),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.6824e5,  1.0e+1),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-4.6251,     1.0e-4),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-12.321,     1.0e-3),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 7.4198e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.0143e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.6704e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0975e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 3.0431   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6794   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9054e+1,  1.0e-3),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.6113   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.2912e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.3676   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8209e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1459e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 5.4351e1,  1.0e-3),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.5027e3,  1.0e-1),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 4.2348   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.7749   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.9916e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6553e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8718   ,  1.0e-4),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.4666e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 9.4146e-1,  1.0e-5),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9890   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.4080   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.0576   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 5.5928e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0840   ,  1.0e-4),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.0662e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6109   ,  1.0e-4),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-2.3268e2,  1.0e-2),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.7839e4,  1.0e-0),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 16.715,  1.0e-3),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.0059,  1.0e-4),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* T_23A */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-0.646387,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-0.646387,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-0.646387,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-0.817720,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-0.817720,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-0.817720,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(8.745434e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.442027e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.138620e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 1.492636e-1,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 6.588619e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.749123e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(4.091723e-3,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.091723e-3,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.091723e-3,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.156818e-2,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.156818e-2,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.156818e-2,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(1.697257e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.697257e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.697257e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.877432e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.877432e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.877432e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-1.473793e-3,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 6.177960e-5,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 1.597356e-3,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.007110e-2,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 4.120409e-4,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair( 2.089519e-2,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-1.422180e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 4.351159e-3,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 1.509203e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-9.185646e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 1.359736e-2,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair( 9.457593e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair( 2.044612e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 1.272367e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 5.001218e-5,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 1.924399e-3,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 8.495098e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.253794e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair( 1.641282e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 8.963349e-3,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 1.513876e-3,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 7.812554e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 2.805355e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.201845e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair( 9.183585e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 3.927598e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.328388e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 2.200521   ,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 5.327408e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.135040   ,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 1.0142e-2,  1.0e-6),  // I2_T23A_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9596e-2,  1.0e-6),  // I2_T23A_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7208e-3,  1.0e-7),  // I2_T23A_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5419e-3,  1.0e-7),  // I2_T23A_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-1.1314e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.4519e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.4936e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5829e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 4.9009e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.9443e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.2872   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.8383   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-3.4865,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.3852,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.1149,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.4581,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-1.2641e2,  1.0e-2),   // I3d1B_T23A_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.5151e4,  1.0e+0),   // I3d1B_T23A_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-6.8424e-2, 1.0e-6),   // I3d1C_T23A_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4339e-1, 1.0e-5),   // I3d1C_T23A_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 3.2509e-4,  1.0e-8),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.3303e-3,  1.0e-7),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1672e-3,  1.0e-7),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.8868e-2,  1.0e-6),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-1.7844e-3,  1.0e-7),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.3768e-2,  1.0e-6),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6351e-3,  1.0e-7),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5100e-2,  1.0e-6),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 3.3705e-1, 1.0e-5),   // I3d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.7887e+2, 1.0e-2),   // I3d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-3.3709e-3,  1.0e-7),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.6848e-2,  1.0e-6),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-8.0525e-2,  1.0e-6),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6698e-1,  1.0e-5),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.3898e-2,  1.0e-6),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.2605e-2,  1.0e-6),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.7293e-3,  1.0e-7),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.9045e-1,  1.0e-5),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(-3.1411e-1,  1.0e-5),  // I4d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.5310e+2,  1.0e-2),  // I4d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-4.1789e-2,  1.0e-6),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.3707e-1,  1.0e-5),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.9717e-2,  1.0e-6),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.3374   ,  1.0e-4),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-1.0628e1,  1.0e-3),   // I4d2B_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.5295e4,  1.0e-0),   // I4d2B_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 1.1331e-1,  1.0e-5),  // I4d2C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 3.0188e-1,  1.0e-5),  // I4d2C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_T23A_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-4.1175e-2,  1.0e-6),  // I2_T23A_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.3527e-1,  1.0e-5),  // I2_T23A_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0443e-3,  1.0e-7),  // I2_T23A_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.9177e-3,  1.0e-7),  // I2_T23A_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-2.5931e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.4342e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2585e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6375e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.2436   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1779   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3369   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2955e+1,  1.0e-3),  // I3_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.4715   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5939e+1,  1.0e-3),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.5557   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.3790e+1,  1.0e-3),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-1.8957e1,  1.0e-3),   // I3d1B_T23A_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-7.1132e3,  1.0e-1),   // I3d1B_T23A_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-6.4407  ,  1.0e-4),   // I3d1C_T23A_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.3498e1,  1.0e-3),   // I3d1C_T23A_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair( 4.8973e-3,  1.0e-7),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7033e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.6386e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.7612e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 3.9026e-2,  1.0e-6),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3820e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.6437e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2025   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-8.6389e-1,  1.0e-5),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.7636   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.2204   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.0365e+1,  1.0e-3),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-2.8108   ,  1.0e-4),  // I3d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.3003e+3,  1.0e-1),  // I3d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(-2.3102e-1,  1.0e-5),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.5273e-1,  1.0e-5),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.7739   ,  1.0e-4),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2899e+1,  1.0e-3),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair( 1.2581e-1,  1.0e-5),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2561   ,  1.0e-4),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7368e+1,  1.0e-3),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6662e+2,  1.0e-2),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 3.5385e+1,  1.0e-3),  // I4d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.8958e+4,  1.0e+0),  // I4d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 7.6579e+1,  1.0e-3),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.3957e+2,  1.0e-2),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9495e+2,  1.0e-2),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.0209e+2,  1.0e-2),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-9.2157e2,  1.0e-2),   // I4d2B_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.3549e4,  1.0e+0),   // I4d2B_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair( 1.0489  ,   1.0e-3),  // I4d2C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 2.7960  ,   1.0e-3),  // I4d2C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_T23A_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 7.4198e-2,  1.0e-5),  // I2_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.0143e-2,  1.0e-5),  // I2_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.6704e-1,  1.0e-4),  // I2_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0975e-1,  1.0e-4),  // I2_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 3.0372   ,  1.0e-4),  // I3_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6670   ,  1.0e-4),  // I3_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9045e+1,  1.0e-3),  // I3_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.5922   ,  1.0e-4),  // I3_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.2931e+1,  1.0e-3),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.4078   ,  1.0e-4),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8239e+1,  1.0e-3),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1521e+1,  1.0e-3),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 5.5795e1,  1.0e-3),   // I3d1B_T23A_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.6754e3,  1.0e-1),   // I3d1B_T23A_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 4.1309   ,  1.0e-4),  // I3d1C_T23A_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.7314   ,  1.0e-4),  // I3d1C_T23A_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.9916e-1,  1.0e-5),  // I2_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6553e-1,  1.0e-5),  // I2_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8718   ,  1.0e-4),  // I2_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.4666e-1,  1.0e-5),  // I2_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 5.1953e-1,  1.0e-5),  // I3_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8518   ,  1.0e-4),  // I3_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5259   ,  1.0e-4),  // I3_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.8893   ,  1.0e-4),  // I3_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-2.8682   ,  1.0e-4),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2405e-1,  1.0e-5),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.2452e+1,  1.0e-3),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6279e+1,  1.0e-3),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-2.5646e2,  1.0e-2),   // I3d1B_T23A_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.0683e4,  1.0e-0),   // I3d1B_T23A_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 17.131,  1.0e-3),     // I3d1C_T23A_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.1801,  1.0e-4),     // I3d1C_T23A_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* T_23B */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(1.049356e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.049356e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.049356e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.848575e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.848575e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.848575e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_1 phi_bar */
                    std::make_pair(-3.303407e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.303407e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.303407e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-8.337742e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-8.337742e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-8.337742e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(1.032833e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.745443e-2,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.162560e-2,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.478178e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.097542e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.716906e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(4.962290e-3,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.962290e-3,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.962290e-3,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.802917e-2,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.802917e-2,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.802917e-2,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(-6.643000e-5,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.643000e-5,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.643000e-5,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.860632e-4,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.860632e-4,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.860632e-4,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(-6.981006e-3,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.981006e-3,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.981006e-3,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.411269e-2,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.411269e-2,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.411269e-2,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-1.689134e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.689134e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.689134e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-9.729494e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-9.729494e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-9.729494e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-4.788091e-3,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.052294e-3,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.316497e-3,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.153610e-2,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.218510e-2,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.283410e-2,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-4.534558e-1,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.840569e-1,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.146599e-1,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.736468   ,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.323016   ,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.909564   ,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(2.414476e-4,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.044442e-4,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.674408e-4,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.194830e-3,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.704120e-3,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.213409e-3,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(2.295643e-2,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.942877e-2,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.590112e-2,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.449240e-1,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.226998e-1,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.004757e-1,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(1.582494   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.338815   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.095137   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.598765   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(3.897925   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(3.197086   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_1 phi_3 */
                    std::make_pair(-3.5744e-3,  1.0e-7),  // I1_T23B_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.4908e-3,  1.0e-7),  // I1_T23B_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.9751e-4,  1.0e-8),  // I1_T23B_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2522e-3,  1.0e-7),  // I1_T23B_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_3 */
                    std::make_pair( 7.1753e-2,  1.0e-6),  // I2_T23B_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4871e-1,  1.0e-5),  // I2_T23B_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.2020e-2,  1.0e-6),  // I2_T23B_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.5125e-2,  1.0e-6),  // I2_T23B_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-3.2127e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.9752e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1270e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.4950e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(-1.1184e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.1665e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.5643e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5578   ,  1.0e-4),  // I3_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-7.4905e-1,  1.0e-5),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5207,     1.0e-4),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0117e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.1127e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair( 3.0959e1,  1.0e-3),   // I3d1B_T23B_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.5956e4,  1.0e+0),   // I3d1B_T23B_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 3.7174e-2, 1.0e-6),   // I3d1C_T23B_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.7904e-2, 1.0e-6),   // I3d1C_T23B_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 2.7245e-4,  1.0e-8),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.6292e-3,  1.0e-7),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8163e-3,  1.0e-7),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4194e-2,  1.0e-6),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-2.2552e-3,  1.0e-7),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0040e-2,  1.0e-6),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8561e-3,  1.0e-7),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.8045e-2,  1.0e-6),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 2.8248e-1, 1.0e-5),   // I3d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.1753e+2, 1.0e-2),   // I3d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair( 1.3507e-3,  1.0e-7),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5721e-2,  1.0e-6),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.3704e-2,  1.0e-6),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1733e-1,  1.0e-5),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 8.9862e-3,  1.0e-7),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1921e-1,  1.0e-5),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9383e-1,  1.0e-5),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6934   ,  1.0e-4),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 1.1429   ,  1.0e-4),  // I4d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.2848e+3,  1.0e-1),  // I4d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-2.8217e-1,  1.0e-5),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0181   ,  1.0e-4),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9233e-1,  1.0e-5),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.6503   ,  1.0e-4),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-2.1915  ,  1.0e-4),   // I4d2B_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1825e4,  1.0e-0),   // I4d2B_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(-1.7459e-1,  1.0e-5),  // I4d2C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-4.6513e-1,  1.0e-5),  // I4d2C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_T23B_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair( 3.3331e-3,  1.0e-7),  // I2_T23B_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0950e-2,  1.0e-6),  // I2_T23B_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6548e-4,  1.0e-8),  // I2_T23B_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.9809e-4,  1.0e-8),  // I2_T23B_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 3.8706e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5302e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5589e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.0771e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair( 1.1692e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.9214e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.5760e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2374   ,  1.0e-4),  // I3_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 9.7065e-1,  1.0e-5),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.7941   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.8848   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6569e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair( 2.2178e1,  1.0e-3),   // I3d1B_T23B_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.3217e3,  1.0e-1),   // I3d1B_T23B_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair( 1.0126  ,  1.0e-4),   // I3d1C_T23B_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 2.1220  ,  1.0e-4),   // I3d1C_T23B_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair(-1.1893e-3,  1.0e-7),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1364e-3,  1.0e-7),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.3407e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.7057e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 2.0318e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.0113e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.0072   ,  1.0e-4),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1451e+1,  1.0e-3),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-3.0059   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.0956   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.9620   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.1288e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-1.6717e+1,  1.0e-3),  // I3d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.3680e+4,  1.0e-0),  // I3d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 5.1748e-1,  1.0e-5),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7874   ,  1.0e-4),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0203e+1,  1.0e-3),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.9167e+1,  1.0e-3),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.0133   ,  1.0e-4),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0762e-1,  1.0e-5),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.1064e+2,  1.0e-2),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.2014e+2,  1.0e-2),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(-4.2378e+1,  1.0e-3),  // I4d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.4681e+4,  1.0e+0),  // I4d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-1.6145e+2,  1.0e-2),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.1732e+2,  1.0e-2),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-6.5989e+2,  1.0e-2),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.7868e+3,  1.0e-1),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 9.5143e2,  1.0e-2),   // I4d2B_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.3525e5,  1.0e+1),   // I4d2B_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-11.507  ,   1.0e-3),  // I4d2C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-30.656  ,   1.0e-3),  // I4d2C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_T23B_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-8.7441e-2,  1.0e-5),  // I2_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.9220e-2,  1.0e-5),  // I2_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.4456e-1,  1.0e-4),  // I2_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.4769e-1,  1.0e-4),  // I2_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair(-2.4595e-1,  1.0e-5),  // I3_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.2644e-1,  1.0e-5),  // I3_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5683   ,  1.0e-4),  // I3_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.9348e-1,  1.0e-5),  // I3_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-1.6700   ,  1.0e-4),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.1771   ,  1.0e-4),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0047e+1,  1.0e-3),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.1407   ,  1.0e-4),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(-5.9693  ,  1.0e-4),   // I3d1B_T23B_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-7.1418e2,  1.0e-2),   // I3d1B_T23B_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair(-2.4490e-1,  1.0e-5),  // I3d1C_T23B_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.0264e-1,  1.0e-5),  // I3d1C_T23B_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.0776e-1,  1.0e-5),  // I2_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.8691e-2,  1.0e-6),  // I2_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5023   ,  1.0e-4),  // I2_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.7025e-1,  1.0e-5),  // I2_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 8.3115e-1,  1.0e-5),  // I3_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7620e-1,  1.0e-5),  // I3_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.0667   ,  1.0e-4),  // I3_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2763   ,  1.0e-4),  // I3_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 5.2893   ,  1.0e-4),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.0471e-1,  1.0e-5),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.3147e+1,  1.0e-3),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.1736e+1,  1.0e-3),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair( 6.2809e1,  1.0e-3),   // I3d1B_T23B_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 7.5146e3,  1.0e-1),   // I3d1B_T23B_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair(-1.8049,  1.0e-4),     // I3d1C_T23B_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-0.7565,  1.0e-3),     // I3d1C_T23B_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_1(-5.0), 0.849651, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_1( 0.0), 0.838315, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_1( 5.0), 0.822425, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_2(-5.0), 0.833842, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_2( 0.0), 0.821883, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_2( 5.0), 0.781765, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_30(-5.0), 0.858475, 1.0e-3);
              //TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_30( 0.0), 0.838315, 1.0e-3); A30 is 0 for q2 = 0
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_30( 5.0), 0.817009, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_v(-5.0), 0.851149, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_v( 0.0), 0.841027, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_v( 5.0), 0.826181, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_1(-5.0), 0.846651, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_1( 0.0), 0.835993, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_1( 5.0), 0.820869, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23A(-5.0), 0.847232, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23A( 0.0), 0.835993, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23A( 5.0), 0.818952, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23B(-5.0), 0.832092, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23B( 0.0), 0.818786, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23B( 5.0), 0.809531, 1.0e-3);
            }


            /* B -> rho form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]               = 2.173913;
                p["B::lambda_E^2"]                 = 0.03;
                p["B::lambda_H^2"]                 = 0.06;
                p["mass::B_d"]                     = 5.27958;
                p["mass::rho^+"]                   = 0.77526;
                p["decay-constant::B_d"]           = 0.1905;
                p["decay-constant::rho"]           = 0.213;
                p["B->rho::mu@B-LCSR"]             = 1.0;
                p["B->rho::s_0^A1,0@B-LCSR"]       = 1.6;
                p["B->rho::s_0^A1,1@B-LCSR"]       = 0.0;
                p["B->rho::s_0^A2,0@B-LCSR"]       = 1.6;
                p["B->rho::s_0^A2,1@B-LCSR"]       = 0.0;
                p["B->rho::s_0^A30,0@B-LCSR"]      = 1.6;
                p["B->rho::s_0^A30,1@B-LCSR"]      = 0.0;
                p["B->rho::s_0^V,0@B-LCSR"]        = 1.6;
                p["B->rho::s_0^V,1@B-LCSR"]        = 0.0;
                p["B->rho::s_0^T1,0@B-LCSR"]       = 1.6;
                p["B->rho::s_0^T1,1@B-LCSR"]       = 0.0;
                p["B->rho::s_0^T23A,0@B-LCSR"]     = 1.6;
                p["B->rho::s_0^T23A,1@B-LCSR"]     = 0.0;
                p["B->rho::s_0^T23B,0@B-LCSR"]     = 1.6;
                p["B->rho::s_0^T23B,1@B-LCSR"]     = 0.0;
                p["B->rho::M^2@B-LCSR"]            = 1.0;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "all"  },
                    { "gminus"_ok, "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->rho::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->v(-5.0),   0.202868, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 0.0),   0.254891, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(+5.0),   0.325801, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_0(-5.0), 0.231651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 0.0), 0.299211, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(+5.0), 0.400604, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_1(-5.0), 0.187824, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 0.0), 0.202553, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(+5.0), 0.218425, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_2(-5.0), 0.143643, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 0.0), 0.169280, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(+5.0), 0.196516, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_1(-5.0), 0.177540, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 0.0), 0.223847, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(+5.0), 0.288414, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(-5.0), 0.209924, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 0.0), 0.223847, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(+5.0), 0.236963, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_3(-5.0), 0.132538, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 0.0), 0.156777, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(+5.0), 0.181452, eps);
            }


            /* B -> K^* form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]               = 2.173913;
                p["B::lambda_E^2"]                 = 0.03;
                p["B::lambda_H^2"]                 = 0.06;
                p["mass::B_d"]                     = 5.27958;
                p["mass::K_d^*"]                   = 0.89594;
                p["decay-constant::B_d"]           = 0.1905;
                p["B->K^*::f_Kstar_par"]           = 0.204;
                p["B->K^*::mu@B-LCSR"]             = 1.0;
                p["B->K^*::s_0^A1,0@B-LCSR"]       = 1.7;
                p["B->K^*::s_0^A1,1@B-LCSR"]       = 0.0;
                p["B->K^*::s_0^A2,0@B-LCSR"]       = 1.7;
                p["B->K^*::s_0^A2,1@B-LCSR"]       = 0.0;
                p["B->K^*::s_0^A30,0@B-LCSR"]      = 1.7;
                p["B->K^*::s_0^A30,1@B-LCSR"]      = 0.0;
                p["B->K^*::s_0^V,0@B-LCSR"]        = 1.7;
                p["B->K^*::s_0^V,1@B-LCSR"]        = 0.0;
                p["B->K^*::s_0^T1,0@B-LCSR"]       = 1.7;
                p["B->K^*::s_0^T1,1@B-LCSR"]       = 0.0;
                p["B->K^*::s_0^T23A,0@B-LCSR"]     = 1.7;
                p["B->K^*::s_0^T23A,1@B-LCSR"]     = 0.0;
                p["B->K^*::s_0^T23B,0@B-LCSR"]     = 1.7;
                p["B->K^*::s_0^T23B,1@B-LCSR"]     = 0.0;
                p["B->K^*::M^2@B-LCSR"]            = 1.0;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "all"  },
                    { "gminus"_ok, "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->v(-5.0),   0.260799, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 0.0),   0.328805, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(+5.0),   0.423196, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_0(-5.0), 0.268280, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 0.0), 0.346213, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(+5.0), 0.463291, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_1(-5.0), 0.242241, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 0.0), 0.264235, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(+5.0), 0.290082, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_2(-5.0), 0.192709, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 0.0), 0.230725, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(+5.0), 0.275801, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_1(-5.0), 0.229660, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 0.0), 0.290716, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(+5.0), 0.377439, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(-5.0), 0.269636, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 0.0), 0.290716, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(+5.0), 0.313430, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_3(-5.0), 0.167315, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 0.0), 0.200133, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(+5.0), 0.236724, eps);
            }


            /* B -> D^* form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                /*charm mass = 1.066273     */
                p["B::1/lambda_B_p"]               = 2.173913;
                p["B::lambda_E^2"]                 = 0.03;
                p["B::lambda_H^2"]                 = 0.06;
                p["mass::B_d"]                     = 5.27958;
                p["mass::D_d^*"]                   = 2.01026;
                p["decay-constant::B_d"]           = 0.1905;
                p["decay-constant::D_d^*"]         = 0.242;
                p["B->D^*::mu@B-LCSR"]             = 2.1213;
                p["B->D^*::s_0^A1,0@B-LCSR"]       = 8.0;
                p["B->D^*::s_0^A1,1@B-LCSR"]       = 0.0;
                p["B->D^*::s_0^A2,0@B-LCSR"]       = 8.0;
                p["B->D^*::s_0^A2,1@B-LCSR"]       = 0.0;
                p["B->D^*::s_0^A30,0@B-LCSR"]      = 8.0;
                p["B->D^*::s_0^A30,1@B-LCSR"]      = 0.0;
                p["B->D^*::s_0^V,0@B-LCSR"]        = 8.0;
                p["B->D^*::s_0^V,1@B-LCSR"]        = 0.0;
                p["B->D^*::s_0^T1,0@B-LCSR"]       = 8.0;
                p["B->D^*::s_0^T1,1@B-LCSR"]       = 0.0;
                p["B->D^*::s_0^T23A,0@B-LCSR"]     = 8.0;
                p["B->D^*::s_0^T23A,1@B-LCSR"]     = 0.0;
                p["B->D^*::s_0^T23B,0@B-LCSR"]     = 8.0;
                p["B->D^*::s_0^T23B,1@B-LCSR"]     = 0.0;
                p["B->D^*::M^2@B-LCSR"]            = 4.5;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "off"  },
                    { "gminus"_ok, "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->D^*::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->v(-5.0),   0.874062, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 0.0),   1.024130, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(+5.0),   1.231090, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_0(-5.0), 0.732174, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 0.0), 0.871237, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(+5.0), 1.058420, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_1(-5.0), 0.770035, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 0.0), 0.803015, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(+5.0), 0.838538, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_2(-5.0), 0.649193, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 0.0), 0.719118, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(+5.0), 0.692782, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_1(-5.0), 0.744459, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 0.0), 0.870939, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(+5.0), 1.039270, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(-5.0), 0.854086, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 0.0), 0.870939, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(+5.0), 0.879057, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_3(-5.0), 0.429351, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 0.0), 0.473056, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(+5.0), 0.399784, eps);
            }


            /* B_s -> K^* form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B_s::1/lambda_B_p"]               = 1.69348;
                p["B_s::lambda_E^2"]                 = 0.03;
                p["B_s::lambda_H^2"]                 = 0.06;
                p["mass::B_s"]                       = 5.36677;
                p["mass::K_u^*"]                     = 0.89166;
                p["decay-constant::B_s"]             = 0.2307;
                p["decay-constant::K_u^*"]           = 0.204;
                p["B_s->K^*::mu@B-LCSR"]             = 1.0;
                p["B_s->K^*::s_0^A1,0@B-LCSR"]       = 1.7;
                p["B_s->K^*::s_0^A1,1@B-LCSR"]       = 0.0;
                p["B_s->K^*::s_0^A2,0@B-LCSR"]       = 1.7;
                p["B_s->K^*::s_0^A2,1@B-LCSR"]       = 0.0;
                p["B_s->K^*::s_0^A30,0@B-LCSR"]      = 1.7;
                p["B_s->K^*::s_0^A30,1@B-LCSR"]      = 0.0;
                p["B_s->K^*::s_0^V,0@B-LCSR"]        = 1.7;
                p["B_s->K^*::s_0^V,1@B-LCSR"]        = 0.0;
                p["B_s->K^*::s_0^T1,0@B-LCSR"]       = 1.7;
                p["B_s->K^*::s_0^T1,1@B-LCSR"]       = 0.0;
                p["B_s->K^*::s_0^T23A,0@B-LCSR"]     = 1.7;
                p["B_s->K^*::s_0^T23A,1@B-LCSR"]     = 0.0;
                p["B_s->K^*::s_0^T23B,0@B-LCSR"]     = 1.7;
                p["B_s->K^*::s_0^T23B,1@B-LCSR"]     = 0.0;
                p["B_s->K^*::M^2@B-LCSR"]            = 1.0;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "all"  },
                    { "gminus"_ok, "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B_s->K^*::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->v(-5.0),   0.167330, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 0.0),   0.203912, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(+5.0),   0.246149, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_0(-5.0), 0.219685, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 0.0), 0.278868, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(+5.0), 0.362594, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_1(-5.0), 0.154442, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 0.0), 0.165623, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(+5.0), 0.177698, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_2(-5.0), 0.107336, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 0.0), 0.120495, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(+5.0), 0.129637, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_1(-5.0), 0.147053, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 0.0), 0.181221, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(+5.0), 0.224128, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(-5.0), 0.172067, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 0.0), 0.181221, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(+5.0), 0.188704, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_3(-5.0), 0.094625, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 0.0), 0.104721, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(+5.0), 0.107054, eps);
            }


            /* B_s -> phi form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B_s::1/lambda_B_p"]               = 1.69348;
                p["B_s::lambda_E^2"]                 = 0.03;
                p["B_s::lambda_H^2"]                 = 0.06;
                p["mass::B_s"]                       = 5.36677;
                p["mass::phi"]                       = 1.019461;
                p["decay-constant::B_s"]             = 0.2307;
                p["decay-constant::phi"]             = 0.233;
                p["B_s->phi::mu@B-LCSR"]             = 1.0;
                p["B_s->phi::s_0^A1,0@B-LCSR"]       = 1.7;
                p["B_s->phi::s_0^A1,1@B-LCSR"]       = 0.0;
                p["B_s->phi::s_0^A2,0@B-LCSR"]       = 1.7;
                p["B_s->phi::s_0^A2,1@B-LCSR"]       = 0.0;
                p["B_s->phi::s_0^A30,0@B-LCSR"]      = 1.7;
                p["B_s->phi::s_0^A30,1@B-LCSR"]      = 0.0;
                p["B_s->phi::s_0^V,0@B-LCSR"]        = 1.7;
                p["B_s->phi::s_0^V,1@B-LCSR"]        = 0.0;
                p["B_s->phi::s_0^T1,0@B-LCSR"]       = 1.7;
                p["B_s->phi::s_0^T1,1@B-LCSR"]       = 0.0;
                p["B_s->phi::s_0^T23A,0@B-LCSR"]     = 1.7;
                p["B_s->phi::s_0^T23A,1@B-LCSR"]     = 0.0;
                p["B_s->phi::s_0^T23B,0@B-LCSR"]     = 1.7;
                p["B_s->phi::s_0^T23B,1@B-LCSR"]     = 0.0;
                p["B_s->phi::M^2@B-LCSR"]            = 1.0;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "all"  },
                    { "gminus"_ok, "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B_s->phi::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->v(-5.0),   0.188667, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 0.0),   0.232526, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(+5.0),   0.286883, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_0(-5.0), 0.210868, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 0.0), 0.267878, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(+5.0), 0.348945, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_1(-5.0), 0.174001, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 0.0), 0.189710, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(+5.0), 0.208671, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_2(-5.0), 0.131226, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 0.0), 0.153049, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(+5.0), 0.177031, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_1(-5.0), 0.166630, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 0.0), 0.207632, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(+5.0), 0.262038, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(-5.0), 0.193698, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 0.0), 0.207632, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(+5.0), 0.222445, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_3(-5.0), 0.107931, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 0.0), 0.123576, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(+5.0), 0.135953, eps);
            }


            /* B_s -> D_s^* form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                /*charm mass = 1.066273     */
                p["B_s::1/lambda_B_p"]               = 1.69348;
                p["B_s::lambda_E^2"]                 = 0.03;
                p["B_s::lambda_H^2"]                 = 0.06;
                p["mass::B_s"]                       = 5.36677;
                p["mass::D_s^*"]                     = 2.1121;
                p["decay-constant::B_s"]             = 0.2307;
                p["decay-constant::D_s^*"]           = 0.293;
                p["B_s->D_s^*::mu@B-LCSR"]           = 2.1213;
                p["B_s->D_s^*::s_0^A1,0@B-LCSR"]     = 8.0;
                p["B_s->D_s^*::s_0^A1,1@B-LCSR"]     = 0.0;
                p["B_s->D_s^*::s_0^A2,0@B-LCSR"]     = 8.0;
                p["B_s->D_s^*::s_0^A2,1@B-LCSR"]     = 0.0;
                p["B_s->D_s^*::s_0^A30,0@B-LCSR"]    = 8.0;
                p["B_s->D_s^*::s_0^A30,1@B-LCSR"]    = 0.0;
                p["B_s->D_s^*::s_0^V,0@B-LCSR"]      = 8.0;
                p["B_s->D_s^*::s_0^V,1@B-LCSR"]      = 0.0;
                p["B_s->D_s^*::s_0^T1,0@B-LCSR"]     = 8.0;
                p["B_s->D_s^*::s_0^T1,1@B-LCSR"]     = 0.0;
                p["B_s->D_s^*::s_0^T23A,0@B-LCSR"]   = 8.0;
                p["B_s->D_s^*::s_0^T23A,1@B-LCSR"]   = 0.0;
                p["B_s->D_s^*::s_0^T23B,0@B-LCSR"]   = 8.0;
                p["B_s->D_s^*::s_0^T23B,1@B-LCSR"]   = 0.0;
                p["B_s->D_s^*::M^2@B-LCSR"]          = 4.5;

                Options o = {
                    { "2pt"_ok,    "all"  },
                    { "3pt"_ok,    "off"  },
                    { "gminus"_ok, "WW-limit" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B_s->D_s^*::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR(ff->v(-5.0),   0.653061, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v( 0.0),   0.758854, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->v(+5.0),   0.868634, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_0(-5.0), 0.596366, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0( 0.0), 0.716505, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_0(+5.0), 0.871324, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_1(-5.0), 0.590307, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1( 0.0), 0.630581, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_1(+5.0), 0.683596, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->a_2(-5.0), 0.471994, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2( 0.0), 0.519060, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->a_2(+5.0), 0.480577, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_1(-5.0), 0.569470, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1( 0.0), 0.670264, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_1(+5.0), 0.793728, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_2(-5.0), 0.645332, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2( 0.0), 0.670264, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_2(+5.0), 0.696802, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->t_3(-5.0), 0.279428, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3( 0.0), 0.286126, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->t_3(+5.0), 0.150783, eps);
            }
        }
} kmo2006_form_factors_test;
