/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2025 Danny van Dyk
 * Copyright (c) 2023 Carolina Bolognani
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

class AnalyticFormFactorBsToKDKMMO2008Test :
    public TestCase
{
    public:
        AnalyticFormFactorBsToKDKMMO2008Test() :
            TestCase("analytic_form_factor_bs_to_k_DKMMO2008_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            /* Factory */
            {
                Parameters p = Parameters::Defaults();
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("B_s->K::DKMMO2008", p, Options{ });

                TEST_CHECK(0 != ff.get());
            }

            /* Decay Constant */
            {
                Parameters p = Parameters::Defaults();
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff(p, Options{ });
                p["mass::B_s"] = 5.3667;
                p["mass::b(MSbar)"] = 4.2;
                p["B_s->K::mu@DKMMO2008"] = 4.2;
                p["B_s->K::Mp^2@DKMMO2008"] = 5.0;       // decay constant
                p["B_s->K::sp_0^B@DKMMO2008"] = 35.75;   // decay constant
                p["B_s->K::s_0^+(0)@DKMMO2008"] = 37.5;  // f_+
                p["B_s->K::s_0^+'(0)@DKMMO2008"] = 0.0;  // f_+
                p["B_s->K::s_0^0(0)@DKMMO2008"] = 37.5;  // f_0
                p["B_s->K::s_0^0'(0)@DKMMO2008"] = 0.0;  // f_0
                p["B_s->K::s_0^T(0)@DKMMO2008"] = 37.5;  // f_T
                p["B_s->K::s_0^T'(0)@DKMMO2008"] = 0.0;  // f_T
                p["QCD::m_0^2"] = 0.8;
                p["QCD::cond_ss@2GeV"] = 0.0;
                p["QCD::cond_GG"] = 0.012;
                p["QCD::r_vac"] = 1.0;

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(-5.05150,  eps), // rho_1(s = 19.60, m_b = 4.16, mu = 4.16)
                    std::make_pair(-4.62757,  eps), // rho_1(s = 22.05, m_b = 4.16, mu = 4.16)
                    std::make_pair(+0.67764,  eps), // rho_1(s = 25.20, m_b = 4.16, mu = 4.16)
                    std::make_pair( 0.20216, 1e-3), // f_B
                    std::make_pair( 5.30431,  eps), // M_B for SVZ
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

                //TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /*
             * B_s -> K f_+ Form Factor at test scale mu = 3.0 GeV.
             * These test values are in reasonably agreement with values
             * derived from the Mathematica notebook graciously
             * provided by Domagoj Leljak.
             */
            {
                static const double eps = 1e-4;

                Parameters p = Parameters::Defaults();
                p["mass::K_u"] = 0.49368;
                p["mass::b(MSbar)"] = 4.18;
                p["mass::d(2GeV)"] = 0.0048;
                p["mass::u(2GeV)"] = 0.0032;
                p["K::a1@1GeV"] =  0.06;
                p["K::a2@1GeV"] =  0.25;
                p["K::a3@1GeV"] =  0.00;
                p["K::a4@1GeV"] = -0.15;
                p["K::f3@1GeV"] = 0.0045;
                p["K::omega3@1GeV"] = -1.5;
                p["K::omega4@1GeV"] = 0.2;
                p["K::delta4@1GeV"] = 0.18;
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

                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff(p, Options{ { "decay-constant"_ok, "sum-rule" } });

                // LO, tw2
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw2(-5.0), 0.1167663129, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw2(-1.0), 0.1484394092, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw2( 0.0), 0.1584577215, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw2( 1.0), 0.169560937, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw2( 5.0), 0.2285685098, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw2(10.0), 0.3595046485, eps);

                // LO, tw3
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw3(-5.0), 0.1261998773, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw3(-1.0), 0.1628625335, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw3( 0.0), 0.1745584606, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw3( 1.0), 0.1875771, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw3( 5.0), 0.2578032862, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw3(10.0), 0.4208530654, eps);

                // LO, tw4
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw4( 0.0), -0.001347845547, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw4( 1.0), -0.001631829059, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw4( 5.0), -0.003401789525, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_lo_tw4(10.0), -0.008687803229, eps);

                // NLO, tw2
                const auto nlo_eps = 400 * eps;
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw2( 0.0), 0.7744550115, nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw2( 1.0), 0.821673966, nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw2( 5.0), 1.055216673, nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw2(10.0), 1.451914987, nlo_eps);

                // NLO, tw3
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw3( 0.0), -0.9050878798, nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw3( 1.0), -0.9780325217, nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw3( 5.0), -1.412757924, nlo_eps);
                TEST_CHECK_NEARLY_EQUAL(ff.F_nlo_tw3(10.0), -2.727613208, nlo_eps);

                // fp form factor @ mu = 3.0
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 0.0), 0.2835562036, 10 * eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 1.0), 0.2992200556, 10 * eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 5.0), 0.3780079263, 10 * eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p(10.0), 0.5345373344, 10 * eps);

                Options o
                {
                    { "rescale-borel"_ok,  "0"        },
                    { "decay-constant"_ok, "sum-rule" }
                };
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff_no_rescale(p, o);

                // Ftil LO, tw3
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw3(-10.0), 0.02832490463, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw3( -1.0), 0.04520072704, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw3( 0.0), 0.04803228036, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw3( 1.0), 0.05115992429, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw3( 5.0), 0.06773084243, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw3( 10.0), 0.1057816771, 1. * eps);

                // Ftil LO, tw4
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw4(-10.0), 0.001036568783, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw4( -1.0), 0.001235548506, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw4( 0.0), 0.001247981061, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw4( 1.0), 0.001255592717, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw4( 5.0), 0.001200905175, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_lo_tw4( 10.0), 0.000621844038, 1. * eps);

                // Ftil NLO, tw2
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw2(-10.0), 0.1980712141, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw2( -1.0), 0.2399580895, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw2( 1e-5), 0.2455887772, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw2( 1.0), 0.2514958573, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw2( 5.0), 0.2774746283, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw2( 10.0), 0.3145996418, 1. * eps);

                // Ftil NLO, tw3
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw3(-10.0), -0.1072388289, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw3( -1.0), -0.1771775591, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw3( 1e-5), -0.1907704872, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw3( 1.0), -0.206426068, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw3( 5.0), -0.3023306927, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.Ftil_nlo_tw3( 10.0), -0.6125574901, 1. * eps);

                // f0 form factor @ mu = 3.0
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0(-10.0), 0.2233836095, 10. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0( -1.0), 0.275827994, 10. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0( 0.0), 0.2835562036, 10. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0( 1.0), 0.2918197449, 10. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0( 5.0), 0.331639046, 10. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0( 10.0), 0.4053014335, 10. * eps);

                // FT LO, tw2
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw2(-10.0), 0.02250363594, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw2( -1.0), 0.03357643699, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw2( 0.0), 0.03535972602, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw2( 1.0), 0.03730682314, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw2( 5.0), 0.04724556234, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw2( 10.0), 0.06789349568, 1. * eps);

                // FT LO, tw3
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw3(-10.0), 0.01371275199, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw3( -1.0), 0.02192018299, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw3( 0.0), 0.02330854566, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw3( 1.0), 0.0248455308, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw3( 5.0), 0.03303953084, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw3( 10.0), 0.05201170974, 1. * eps);

                // FT LO, tw4
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw4(-10.0), -0.0007575148009, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw4( -1.0), -0.001441505804, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw4( 0.0), -0.001565969574, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw4( 1.0), -0.00170610462, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw4( 5.0), -0.002487417233, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_lo_tw4( 10.0), -0.004441057184, 1. * eps);

                // FT NLO, tw2
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw2(-10.0), 0.102793122, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw2(- 1.0), 0.1455997986, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw2( 0.0), 0.151722538, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw2( 1.0), 0.1581471538, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw2( 5.0), 0.1864415337, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw2( 10.0), 0.2181614612, 1. * eps);

                // FT NLO, tw3
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw3(-10.0), -0.03143549821, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw3(- 1.0), -0.06027476052, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw3( 0.0), -0.06654060366, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw3( 1.0), -0.07395675912, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw3( 5.0), -0.1220326809, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.FT_nlo_tw3( 10.0), -0.2878912548, 1. * eps);

                // fT form factor @ mu = 3.0
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t(-10.0), 0.1750151738, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( -5.0), 0.2170830843, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( -1.0), 0.2638103349, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( 0.0), 0.278050279, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( 1.0), 0.2935706186, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( 5.0), 0.3722379609, 1. * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( 10.0), 0.5313809249, 1. * eps);
            }

            {
                // Comparison with Blazenka's notebook underlysing the [DKMMO:2008A] results
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
                    { "rescale-borel"_ok,  "0"        },
                    { "decay-constant"_ok, "sum-rule" }
                };
                AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down> ff_no_rescale(p, o);

                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_p( 0.0), 0.264200304, 2 * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_p( 10.0), 0.4975396448, 15 * eps);
                // f_0(0) = f_+(0)
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_0( 10.0), 0.373123391, 7 * eps);

                // The values for f_T used here differe from the published manuscript due to a typeo
                // in the formulas for the leading-order expression. The shift is ~2%, and the values
                // below are taken from an updated Mathematica notebook free of this typo.
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( 0.0), 0.2612949722, 10 * eps);
                TEST_CHECK_NEARLY_EQUAL(ff_no_rescale.f_t( 10.0), 0.4973622224, 19 * eps);

                //TEST_CHECK(false);
            }
        }
} analytic_form_factor_bs_to_k_DKMMO2008_test;
