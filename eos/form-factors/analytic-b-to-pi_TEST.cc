/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2023 Danny van Dyk
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
#include <eos/form-factors/analytic-b-to-psd-dkmmo2008.hh>
#include <eos/form-factors/mesonic.hh>

#include <cmath>
#include <limits>
#include <vector>

#include <iostream>
#include <array>

using namespace test;
using namespace eos;

class AnalyticFormFactorBToPiDKMMO2008Test :
    public TestCase
{
    public:
        AnalyticFormFactorBToPiDKMMO2008Test() :
            TestCase("analytic_form_factor_b_to_pi_DKMMO2008_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            /* Factory */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B->pi::DKMMO2008", p, Options{ });

                TEST_CHECK(0 != ff.get());
            }

            /* Decay Constant */
            {
                Parameters p = Parameters::Defaults();
                Options o{
                    { "decay-constant", "sum-rule" }
                };
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff(p, o);
                p["mass::B_d"] = 5.2795;
                p["mass::b(MSbar)"] = 4.2;
                p["B->pi::mu@DKMMO2008"] = 4.2;
                p["B->pi::Mp^2@DKMMO2008"] = 5.0;       // decay constant
                p["B->pi::sp_0^B@DKMMO2008"] = 35.75;   // decay constant
                p["B->pi::s_0^+(0)@DKMMO2008"] = 37.5;  // f_+
                p["B->pi::s_0^+'(0)@DKMMO2008"] = 0.0;  // f_+
                p["B->pi::s_0^0(0)@DKMMO2008"] = 37.5;  // f_0
                p["B->pi::s_0^0'(0)@DKMMO2008"] = 0.0;  // f_0
                p["B->pi::s_0^T(0)@DKMMO2008"] = 37.5;  // f_T
                p["B->pi::s_0^T'(0)@DKMMO2008"] = 0.0;  // f_T
                p["QCD::m_0^2"] = 0.8;
                p["QCD::cond_GG"] = 0.012;
                p["QCD::r_vac"] = 1.0;

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(-5.05150,  eps), // rho_1(s = 19.60, m_b = 4.16, mu = 4.16)
                    std::make_pair(-4.62757,  eps), // rho_1(s = 22.05, m_b = 4.16, mu = 4.16)
                    std::make_pair(+0.67764,  eps), // rho_1(s = 25.20, m_b = 4.16, mu = 4.16)
                    std::make_pair( 0.22315, 1e-3), // f_B
                    std::make_pair( 5.33019,  eps), // M_B for SVZ
                    std::make_pair( 1.00000,  eps), // rescale factor for f_+ at s =  0.0 GeV^2
                    std::make_pair( 1.09380,  eps), // rescale factor for f_+ at s = 10.0 GeV^2
                    std::make_pair( 1.00000,  eps), // rescale factor for f_0 at s =  0.0 GeV^2
                    std::make_pair( 1.14083,  eps), // rescale factor for f_0 at s = 10.0 GeV^2
                    std::make_pair( 1.00000,  eps), // rescale factor for f_T at s =  0.0 GeV^2
                    std::make_pair( 1.07377,  eps), // rescale factor for f_T at s = 10.0 GeV^2
                    std::make_pair( 5.30187,  eps), // M_B for f_+ at s =  0.0 GeV^2
                    std::make_pair( 5.32078,  eps), // M_B for f_+ at s = 10.0 GeV^2
                    std::make_pair( 5.30187,  eps), // M_B for f_0 at s =  0.0 GeV^2
                    std::make_pair( 5.35957,  eps), // M_B for f_0 at s = 10.0 GeV^2
                    std::make_pair( 5.30246,  eps), // M_B for f_T at s =  0.0 GeV^2
                    std::make_pair( 5.34903,  eps), // M_B for f_T at s = 10.0 GeV^2
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /*
             * B -> pi f_+ Form Factor at test scale mu = 3.0 GeV.
             * These test values are in reasonably agreement with values
             * derived from the Mathematica notebook graciously
             * provided by I. Sentitemsu Imsong.
             */
            {
                static const double eps = 1e-4;

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"] = 0.13957;
                p["mass::b(MSbar)"] = 4.18;
                p["mass::d(2GeV)"] = 0.0048;
                p["mass::u(2GeV)"] = 0.0032;
                p["pi::a2@1GeV"] = 0.17;
                p["pi::a4@1GeV"] = 0.06;
                p["pi::f3@1GeV"] = 0.0045;
                p["pi::omega3@1GeV"] = -1.5;
                p["pi::omega4@1GeV"] = 0.2;
                p["pi::delta4@1GeV"] = 0.18;
                p["B->pi::M^2@DKMMO2008"] = 12.0;
                p["B->pi::Mp^2@DKMMO2008"] = 4.5;
                p["B->pi::mu@DKMMO2008"] = 3.0;
                p["B->pi::s_0^+(0)@DKMMO2008"] = 37.5;
                p["B->pi::s_0^+'(0)@DKMMO2008"] = 0.0;
                p["B->pi::s_0^0(0)@DKMMO2008"] = 37.5;
                p["B->pi::s_0^0'(0)@DKMMO2008"] = 0.0;
                p["B->pi::s_0^T(0)@DKMMO2008"] = 37.5;
                p["B->pi::s_0^T'(0)@DKMMO2008"] = 0.0;
                p["B->pi::sp_0^B@DKMMO2008"] = 36.5;
                p["QCD::m_0^2"] = 0.8;
                p["QCD::cond_GG"] = 0.012;
                p["QCD::r_vac"] = 1.0;
                p["QCD::alpha_s(MZ)"] = 0.1184;

                Options o{
                    { "decay-constant", "sum-rule" }
                };
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff(p, o);

                // LO, tw2
                TEST_CHECK_NEARLY_EQUAL( 0.1167, ff.F_lo_tw2(-5.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1484, ff.F_lo_tw2(-1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1584, ff.F_lo_tw2( 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1696, ff.F_lo_tw2( 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2290, ff.F_lo_tw2( 5.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.3604, ff.F_lo_tw2(10.0),  eps);

                // LO, tw3
                TEST_CHECK_NEARLY_EQUAL( 0.1261, ff.F_lo_tw3(-5.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1628, ff.F_lo_tw3(-1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1746, ff.F_lo_tw3( 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1876, ff.F_lo_tw3( 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2580, ff.F_lo_tw3( 5.0),  eps);
                TEST_CHECK_NEARLY_EQUAL( 0.4214, ff.F_lo_tw3(10.0),  eps);

                // LO, tw4
                TEST_CHECK_NEARLY_EQUAL(-0.0013, ff.F_lo_tw4( 0.0),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0016, ff.F_lo_tw4( 1.0),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0034, ff.F_lo_tw4( 5.0),  eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0087, ff.F_lo_tw4(10.0),  eps);

                // NLO, tw2
                const auto nlo_eps = 400 * eps;
                TEST_CHECK_NEARLY_EQUAL(+0.7706, ff.F_nlo_tw2( 0.0), nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(+0.8190, ff.F_nlo_tw2( 1.0), nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(+1.0609, ff.F_nlo_tw2( 5.0), nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(+1.4741, ff.F_nlo_tw2(10.0), nlo_eps);

                // NLO, tw3
                TEST_CHECK_NEARLY_EQUAL(-0.9221, ff.F_nlo_tw3( 0.0), nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(-0.9963, ff.F_nlo_tw3( 1.0), nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(-1.4371, ff.F_nlo_tw3( 5.0), nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(-2.7571, ff.F_nlo_tw3(10.0), nlo_eps);

                // fp form factor @ mu = 3.0
                TEST_CHECK_NEARLY_EQUAL( 0.2831, ff.f_p( 0.0), 10 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2988, ff.f_p( 1.0), 10 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.3777, ff.f_p( 5.0), 10 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.5346, ff.f_p(10.0), 10 * eps);

                o = Options{
                    { "decay-constant", "sum-rule" },
                    { "rescale-borel",  "0"        }
                };
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff_no_rescale(p, o);

                // Ftil LO, tw3
                TEST_CHECK_NEARLY_EQUAL( 0.0283, ff_no_rescale.Ftil_lo_tw3(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0452, ff_no_rescale.Ftil_lo_tw3( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0480, ff_no_rescale.Ftil_lo_tw3(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0512, ff_no_rescale.Ftil_lo_tw3(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0677, ff_no_rescale.Ftil_lo_tw3(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1058, ff_no_rescale.Ftil_lo_tw3( 10.0), 1. * eps);

                // Ftil LO, tw4
                TEST_CHECK_NEARLY_EQUAL( 0.0010, ff_no_rescale.Ftil_lo_tw4(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0012, ff_no_rescale.Ftil_lo_tw4( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0012, ff_no_rescale.Ftil_lo_tw4(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0013, ff_no_rescale.Ftil_lo_tw4(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0012, ff_no_rescale.Ftil_lo_tw4(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0006, ff_no_rescale.Ftil_lo_tw4( 10.0), 1. * eps);

                // Ftil NLO, tw2
                TEST_CHECK_NEARLY_EQUAL( 0.1980, ff_no_rescale.Ftil_nlo_tw2(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2397, ff_no_rescale.Ftil_nlo_tw2( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2454, ff_no_rescale.Ftil_nlo_tw2( 1e-5), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2513, ff_no_rescale.Ftil_nlo_tw2(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2775, ff_no_rescale.Ftil_nlo_tw2(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.3147, ff_no_rescale.Ftil_nlo_tw2( 10.0), 1. * eps);

                // Ftil NLO, tw3
                TEST_CHECK_NEARLY_EQUAL(-0.1072, ff_no_rescale.Ftil_nlo_tw3(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1772, ff_no_rescale.Ftil_nlo_tw3( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1907, ff_no_rescale.Ftil_nlo_tw3( 1e-5), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.2064, ff_no_rescale.Ftil_nlo_tw3(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.3023, ff_no_rescale.Ftil_nlo_tw3(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.6126, ff_no_rescale.Ftil_nlo_tw3( 10.0), 1. * eps);

                // f0 form factor @ mu = 3.0
                TEST_CHECK_NEARLY_EQUAL( 0.2234, ff_no_rescale.f_0(-10.0), 10. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2757, ff_no_rescale.f_0( -1.0), 10. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2835, ff_no_rescale.f_0(  0.0), 10. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2918, ff_no_rescale.f_0(  1.0), 10. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.3318, ff_no_rescale.f_0(  5.0), 10. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.4057, ff_no_rescale.f_0( 10.0), 10. * eps);

                // FT LO, tw2
                TEST_CHECK_NEARLY_EQUAL( 0.0225, ff_no_rescale.FT_lo_tw2(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0336, ff_no_rescale.FT_lo_tw2( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0354, ff_no_rescale.FT_lo_tw2(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0373, ff_no_rescale.FT_lo_tw2(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0473, ff_no_rescale.FT_lo_tw2(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0680, ff_no_rescale.FT_lo_tw2( 10.0), 1. * eps);

                // FT LO, tw3
                TEST_CHECK_NEARLY_EQUAL( 0.0137, ff_no_rescale.FT_lo_tw3(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0219, ff_no_rescale.FT_lo_tw3( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0233, ff_no_rescale.FT_lo_tw3(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0248, ff_no_rescale.FT_lo_tw3(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0330, ff_no_rescale.FT_lo_tw3(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0520, ff_no_rescale.FT_lo_tw3( 10.0), 1. * eps);

                // FT LO, tw4
                TEST_CHECK_NEARLY_EQUAL(-0.0008, ff_no_rescale.FT_lo_tw4(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0014, ff_no_rescale.FT_lo_tw4( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0016, ff_no_rescale.FT_lo_tw4(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0017, ff_no_rescale.FT_lo_tw4(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0025, ff_no_rescale.FT_lo_tw4(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0044, ff_no_rescale.FT_lo_tw4( 10.0), 1. * eps);

                // FT NLO, tw2
                TEST_CHECK_NEARLY_EQUAL( 0.1014, ff_no_rescale.FT_nlo_tw2(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1443, ff_no_rescale.FT_nlo_tw2(- 1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1506, ff_no_rescale.FT_nlo_tw2(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1573, ff_no_rescale.FT_nlo_tw2(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.1870, ff_no_rescale.FT_nlo_tw2(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2211, ff_no_rescale.FT_nlo_tw2( 10.0), 1. * eps);

                // FT NLO, tw3
                TEST_CHECK_NEARLY_EQUAL(-0.0314, ff_no_rescale.FT_nlo_tw3(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0603, ff_no_rescale.FT_nlo_tw3(- 1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0665, ff_no_rescale.FT_nlo_tw3(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.0740, ff_no_rescale.FT_nlo_tw3(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.1220, ff_no_rescale.FT_nlo_tw3(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(-0.2879, ff_no_rescale.FT_nlo_tw3( 10.0), 1. * eps);

                // fT form factor @ mu = 3.0
                TEST_CHECK_NEARLY_EQUAL( 0.1750, ff_no_rescale.f_t(-10.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2169, ff_no_rescale.f_t( -5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2636, ff_no_rescale.f_t( -1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2779, ff_no_rescale.f_t(  0.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.2935, ff_no_rescale.f_t(  1.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.3725, ff_no_rescale.f_t(  5.0), 1. * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.5322, ff_no_rescale.f_t( 10.0), 1. * eps);
            }

            {
                // Comparison with Blazenka's notebook underlysing the [DKKMO:2008A] results
                static const double eps = 1e-4;

                Parameters p = Parameters::Defaults();
                p["decay-constant::pi"] = 0.1307;
                p["mass::B_d"] = 5.279;
                p["mass::pi^+"] = 0.13957;
                p["mass::b(MSbar)"] = 4.164;
                p["mass::d(2GeV)"] = 0.006;
                p["mass::u(2GeV)"] = 0.003;
                p["pi::a2@1GeV"] = 0.161995;
                p["pi::a4@1GeV"] = 0.038004;
                p["pi::f3@1GeV"] = 0.0045;
                p["pi::omega3@1GeV"] = -1.5;
                p["pi::omega4@1GeV"] = 0.2;
                p["pi::delta4@1GeV"] = 0.18;
                p["B->pi::M^2@DKMMO2008"] = 18.0;
                p["B->pi::Mp^2@DKMMO2008"] = 5.;
                p["B->pi::mu@DKMMO2008"] = 3.0;
                p["B->pi::s_0^+(0)@DKMMO2008"] = 35.75;
                p["B->pi::s_0^+'(0)@DKMMO2008"] = 0.0;
                p["B->pi::s_0^0(0)@DKMMO2008"] = 35.75;
                p["B->pi::s_0^0'(0)@DKMMO2008"] = 0.0;
                p["B->pi::s_0^T(0)@DKMMO2008"] = 35.75;
                p["B->pi::s_0^T'(0)@DKMMO2008"] = 0.0;
                p["B->pi::sp_0^B@DKMMO2008"] = 35.6;
                p["QCD::m_0^2"] = 0.8;
                p["QCD::cond_GG"] = 0.012;
                p["QCD::r_vac"] = 1.0;
                p["QCD::alpha_s(MZ)"] = 0.1176;

                Options o{
                    { "decay-constant", "sum-rule" },
                    { "rescale-borel",  "0"        }
                };
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff_no_rescale(p, o);

                TEST_CHECK_NEARLY_EQUAL( 0.2641, ff_no_rescale.f_p(  0.0),   2 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.4964, ff_no_rescale.f_p( 10.0),  15 * eps);
                // f_0(0) = f_+(0)
                TEST_CHECK_NEARLY_EQUAL( 0.3725, ff_no_rescale.f_0( 10.0),   7 * eps);

                // The values for f_T used here differe from the published manuscript due to a typeo
                // in the formulas for the leading-order expression. The shift is ~2%, and the values
                // below are taken from an updated Mathematica notebook free of this typo.
                TEST_CHECK_NEARLY_EQUAL( 0.2606, ff_no_rescale.f_t(  0.0),  10 * eps);
                TEST_CHECK_NEARLY_EQUAL( 0.4990, ff_no_rescale.f_t( 10.0),  19 * eps);
            }
        }
} analytic_form_factor_b_to_pi_DKMMO2008_test;
