/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
 * Copyright (c) 2018      Danny van Dyk
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
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/maths/integrate.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class AnalyticFormFactorBToPiPiBFvD2016Test :
    public TestCase
{
    public:
        AnalyticFormFactorBToPiPiBFvD2016Test() :
            TestCase("analytic_form_factor_b_to_pi_pi_BFvD2016_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.1e-5;

            Parameters p = Parameters::Defaults();
            p["mass::B_d"] = 5.2795;
            p["mass::pi^+"] = 0.13957;
            p["mass::d(2GeV)"] = 0.0048;
            p["mass::u(2GeV)"] = 0.0032;
            p["pi::a2@1GeV"] = 0.17;
            p["pi::a4@1GeV"] = 0.06;
            p["decay-constant::pi"] = 0.1304;
            p["B->pipi::mu@BFvD2016"] = 1.5;
            p["B->pi::f_+(0)@BCL2008"] = +2.6000e-01;
            p["B->pi::b_+^1@BCL2008"]  = +0.0000e+00;
            p["B->pi::b_+^2@BCL2008"]  = +0.0000e+00;
            p["B->pi::b_0^1@BCL2008"]  = +0.0000e+00;
            p["B->pi::b_0^2@BCL2008"]  = +0.0000e+00;

            /* Factory */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToPP>> ff = FormFactorFactory<PToPP>::create("B->pipi::BFvD2016", p, Options{ });

                TEST_CHECK(0 != ff.get());
            }

            /* Internal Diagonostics: Integrals */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    // LO, tw2, I_1: integral over f_1
                    std::make_pair(+0.682961, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(+0.735448, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(+0.799572, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    std::make_pair(+0.703571, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(+0.742255, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(+0.786845, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    // LO, tw2, I_2: integral over f_2
                    std::make_pair(+3.391770, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(+3.239920, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(+3.061500, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    std::make_pair(+3.354750, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(+3.244530, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(+3.120950, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    // LO, tw3, I_{sigma,1}: integral over f_{sigma,1}
                    std::make_pair(-0.431354, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(-0.476524, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(-0.534146, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    std::make_pair(-0.447129, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(-0.480559, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(-0.520284, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    // LO, tw3, I_{sigma,2}: integral over f_{sigma,2}
                    std::make_pair(-0.427658, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(-0.449110, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(-0.474262, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    std::make_pair(-0.437860, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(-0.453609, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(-0.471241, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    // LO, tw3, I_{finite}: integral over f_{finite} @ tw3
                    std::make_pair(+1.274055, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(+1.372867, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(+1.495634, eps), // q2 = 0.0112245 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                    std::make_pair(+1.297684, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = -1)
                    std::make_pair(+1.369796, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z =  0)
                    std::make_pair(+1.453899, eps), // q2 = 0.0224490 M_B^2, k2 = 0.6666667 M_B^2, z = +1)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* F_perp, fixed scale */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.30, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.30, 18.60, -1.0)), +0.001599203, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.30, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.30, 18.60,  0.0)), +0.002372988, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.30, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.30, 18.60, +1.0)), +0.003660134, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.60, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.60, 18.60, -1.0)), +0.001327286, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.60, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.60, 18.60,  0.0)), +0.001771697, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.60, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.60, 18.60, +1.0)), +0.002414450, eps);
            }

            /* F_para, fixed scale */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.30, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.30, 18.60, -1.0)), -0.014711514, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.30, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.30, 18.60,  0.0)), -0.021592105, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.30, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.30, 18.60, +1.0)), -0.032743147, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.60, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.60, 18.60, -1.0)), -0.016570605, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.60, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.60, 18.60,  0.0)), -0.021951495, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.60, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.60, 18.60, +1.0)), -0.029590914, eps);
            }

            /* F_long, fixed scale */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 18.60, -1.0)), +0.117508425, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 18.60,  0.0)), +0.102828699, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 18.60, +1.0)),  0.059754354, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60, -1.0)), +0.069270859, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60,  0.0)), +0.052968246, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60, +1.0)), +0.021346673, eps);
            }

            /* F_time, fixed scale */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 18.60, -1.0)), +0.126396815, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 18.60,  0.0)), +0.117632744, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 18.60, +1.0)), +0.086289383, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60, -1.0)), +0.084516600, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60,  0.0)), +0.078294218, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60, +1.0)), +0.063515338, eps);
            }

            /* F_long, variable scale */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ { "scale"_ok, "variable" } });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 18.60, -1.0)), +0.157064905, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 18.60,  0.0)), +0.133953013, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 18.60, +1.0)), +0.070974841, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60, -1.0)), +0.093299574, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60,  0.0)), +0.068992421, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60, +1.0)), +0.023277423, eps);
            }

            /* F_time, variable scale */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ { "scale"_ok, "variable" } });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 18.60, -1.0)), +0.167697744, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 18.60,  0.0)), +0.153005539, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 18.60, +1.0)), +0.106492190, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60, -1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60, -1.0)), +0.111584771, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60,  0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60,  0.0)), +0.101649520, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60, +1.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60, +1.0)), +0.079671685, eps);
            }

            p["B->pipi::mu@BFvD2016"] = 5.2795 / 2.0; /* for tests with Thorsten */

            /* F_perp, comparison with Thorsten */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.30, 20.00, -0.25)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.30, 20.00, -0.25)), +0.000950503, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_perp(0.60, 18.60, +0.50)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_perp(0.60, 18.60, +0.50)), +0.001190890, eps);
            }

            /* F_para, comparison with Thorsten */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.30, 20.00, -0.25)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.30, 20.00, -0.25)), -0.013625800, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_para(0.60, 18.60, +0.50)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_para(0.60, 18.60, +0.50)), -0.018046200, eps);
            }

            /* F_long, comparison with Thorsten */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.30, 20.00, -0.25)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.30, 20.00, -0.25)), +0.057997292, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_long(0.60, 18.60, +0.50)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_long(0.60, 18.60, +0.50)), +0.032167500, eps);
            }

            /* F_time, comparison with Thorsten */
            {
                AnalyticFormFactorBToPiPiBFvD2016 ff(p, Options{ });

                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.30, 20.00, -0.25)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.30, 20.00, -0.25)), +0.070875500, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(ff.f_time(0.60, 18.60, +0.50)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(ff.f_time(0.60, 18.60, +0.50)), +0.057476100, eps);
            }
        }
} analytic_form_factor_b_to_pi_pi_BFvD2016_test;

class AnalyticFormFactorBToPiPiFvDV2018Test :
    public TestCase
{
    public:
        AnalyticFormFactorBToPiPiFvDV2018Test() :
            TestCase("analytic_form_factor_b_to_pi_pi_FvDV2018_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-3;

            /* Factory */
            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"] = 5.27958;
                p["mass::B_d^*"] = 5.32465;

                std::shared_ptr<FormFactors<PToPP>> ff = FormFactorFactory<PToPP>::create("B->pipi::FvDV2018-Dispersive", p, Options{ });
                TEST_CHECK(0 != ff.get());
            }

            /* Internal Diagnostics */
            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"] = 5.27958;
                p["mass::B_d^*"] = 5.32465;

                AnalyticFormFactorBToPiPiFvDV2018 ff(p, Options{ });
                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  2910.308,  eps), // f_time_im_res_qhat2(q2 = 0.05, k2 = 13.0)
                    std::make_pair(  2927.843,  eps), // f_long_im_res_qhat2(q2 = 0.05, k2 = 13.0)
                    std::make_pair(   -46.067,  eps), // f_perp_im_res_qhat2(q2 = 0.05, k2 = 13.0)
                    std::make_pair(   129.103,  eps), // f_para_im_res_qhat2(q2 = 0.05, k2 = 13.0)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }
        }
} analytic_form_factor_b_to_pi_pi_FvDV2018_test;
