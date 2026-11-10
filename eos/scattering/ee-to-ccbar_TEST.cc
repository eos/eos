/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2025-2026 Danny van Dyk
 * Copyright (c) 2026 Simon Mutke
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

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/utils/kmatrix-impl.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <memory>
#include <string>


using namespace test;
using namespace eos;

class eetoccbarTest : public TestCase
{
    public:
        eetoccbarTest() :
            TestCase("ee->ccbar scattering tests")
        {
        }

        virtual void
        run() const
        {
            constexpr double eps = 1e-5;

            {
                Parameters p                            = Parameters::Defaults();
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]      = 10.0;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]    = 12.0;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"]   = 0.3;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 0.4;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"]      = 0.5;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]    = 0.6;

                // Set all cst to zero
                p["ee->ccbar::c(e^+e^-,e^+e^-)"]       = 0.0;
                p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"]    = 0.0;
                p["ee->ccbar::c(e^+e^-,D^+D^-)"]       = 0.0;
                p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"] = 0.0;
                p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"]    = 0.0;
                p["ee->ccbar::c(D^+D^-,D^+D^-)"]       = 0.0;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                // Masses
                p["mass::psi(2S)"]   = 3.68610;
                p["mass::psi(3770)"] = 3.7737;
                p["mass::e"]         = 5.10999e-4;
                p["mass::D^0"]       = 1.86483;
                p["mass::D^+"]       = 1.86965;

                // Build K Matrix
                auto psi2S   = std::make_shared<CharmoniumResonance<3, 2>>("psi2S", p["mass::psi(2S)"]);
                auto psi3770 = std::make_shared<CharmoniumResonance<3, 2>>("psi3770", p["mass::psi(3770)"]);

                std::array<std::array<Parameter, 3>, 3> bkgcst{
                    p["ee->ccbar::c(e^+e^-,e^+e^-)"],    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"],    p["ee->ccbar::c(e^+e^-,D^+D^-)"],
                    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],
                    p["ee->ccbar::c(e^+e^-,D^+D^-)"],    p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],    p["ee->ccbar::c(D^+D^-,D^+D^-)"],
                };

                std::array<Parameter, 2> ee_g0s{
                    { p["ee->ccbar::g0(psi(2S),e^+e^-)"], p["ee->ccbar::g0(psi(3770),e^+e^-)"] }
                };
                std::array<Parameter, 2> D0Dbar0_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"], p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] }
                };
                std::array<Parameter, 2> DpDm_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^+D^-)"], p["ee->ccbar::g0(psi(3770),D^+D^-)"] }
                };

                auto ee      = std::make_shared<EEChannel<3, 2>>("ee", p["mass::e"], p["mass::e"], p["ee->ccbar::q_0"], ee_g0s);
                auto D0Dbar0 = std::make_shared<PWavePPChannel<3, 2>>("D0Dbar0", p["mass::D^0"], p["mass::D^0"], p["ee->ccbar::q_0"], D0Dbar0_g0s);
                auto DpDm    = std::make_shared<PWavePPChannel<3, 2>>("DpDm", p["mass::D^+"], p["mass::D^+"], p["ee->ccbar::q_0"], DpDm_g0s);

                KMatrix<3, 2> KMatrix32({ ee, D0Dbar0, DpDm }, { psi2S, psi3770 }, bkgcst, "ee->ccbar");

                TEST_CHECK_RELATIVE_ERROR(std::abs(ee->rho(20)), 0.0198944, eps);
                TEST_CHECK_RELATIVE_ERROR(std::abs(DpDm->rho(20)), 0.0109126, eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->rho(complex<double>(20, 5))), 0.0119492, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(DpDm->rho(complex<double>(20, 5))), 0.00272429, eps);

                TEST_CHECK_RELATIVE_ERROR(real(ee->chew_mandelstam(20)), -0.1149620, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ee->chew_mandelstam(20)), 0.0198944, eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->chew_mandelstam(-20)), -0.0125291, eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->chew_mandelstam(2)), -0.0081324, eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->chew_mandelstam(20)), 0.00024738, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(DpDm->chew_mandelstam(20)), 0.0093575, eps);

                auto tMatrix32s0 = KMatrix32.tmatrix_row(0, 0.0001);
                auto tMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
                auto tMatrix32s2 = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].real(), 10.1176, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].imag(), 2.11433, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].real(), 0.0534574, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].imag(), 0.00902802, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].real(), 0.313608, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].imag(), 0.0691847, eps);


                // Test the implentation of background constants
                p["ee->ccbar::c(e^+e^-,D^+D^-)"] = 0.5;

                tMatrix32s0 = KMatrix32.tmatrix_row(0, 0.0001);
                tMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
                tMatrix32s2 = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_NEARLY_EQUAL(tMatrix32s0[0].real(), 10.1145, eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s0[0].imag(), 2.11295, eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s1[1].real(), 0.0534576, eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s1[1].imag(), 0.00902805, eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s2[1].real(), 0.313633, eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s2[1].imag(), 0.0690969, eps);
            }

            {
                Parameters p                            = Parameters::Defaults();
                p["mass::psi(3770)"]                    = 3.796443282051135;
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]      = -0.02077753547690923;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]    = -0.001999916489715092;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"]   = -3.3693829070086214;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 8.38428874933062;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"]      = 3.5692138012280807;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]    = -8.391199026268827;
                p["ee->ccbar::g0(psi(2S),eff(2S))"]     = 0.09237874214140328;
                p["ee->ccbar::g0(psi(3770),eff(3770))"] = -1.5409316160476978;

                p["ee->ccbar::effective_mass"] = 0.1349768;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                // Masses
                p["mass::psi(2S)"] = 3.68610;
                p["mass::e"]       = 5.10999e-4;
                p["mass::D^0"]     = 1.86483;
                p["mass::D^+"]     = 1.86965;

                // Build K Matrix
                auto psi2S   = std::make_shared<CharmoniumResonance<5, 2>>("psi2S", p["mass::psi(2S)"]);
                auto psi3770 = std::make_shared<CharmoniumResonance<5, 2>>("psi3770", p["mass::psi(3770)"]);

                std::array<std::array<Parameter, 5>, 5> bkgcst{
                    p["ee->ccbar::c(e^+e^-,e^+e^-)"],       p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"],    p["ee->ccbar::c(e^+e^-,D^+D^-)"],       p["ee->ccbar::c(e^+e^-,eff(2S))"],
                    p["ee->ccbar::c(e^+e^-,eff(3770))"],    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"],    p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],
                    p["ee->ccbar::c(D^0Dbar^0,eff(2S))"],   p["ee->ccbar::c(D^0Dbar^0,eff(3770))"], p["ee->ccbar::c(e^+e^-,D^+D^-)"],       p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],
                    p["ee->ccbar::c(D^+D^-,D^+D^-)"],       p["ee->ccbar::c(D^+D^-,eff(2S))"],      p["ee->ccbar::c(D^+D^-,eff(3770))"],    p["ee->ccbar::c(e^+e^-,eff(2S))"],
                    p["ee->ccbar::c(D^0Dbar^0,eff(2S))"],   p["ee->ccbar::c(D^+D^-,eff(2S))"],      p["ee->ccbar::c(eff(2S),eff(2S))"],     p["ee->ccbar::c(eff(2S),eff(3770))"],
                    p["ee->ccbar::c(e^+e^-,eff(3770))"],    p["ee->ccbar::c(D^0Dbar^0,eff(3770))"], p["ee->ccbar::c(D^+D^-,eff(3770))"],    p["ee->ccbar::c(eff(2S),eff(3770))"],
                    p["ee->ccbar::c(eff(3770),eff(3770))"],
                };

                std::array<Parameter, 2> ee_g0s{
                    { p["ee->ccbar::g0(psi(2S),e^+e^-)"], p["ee->ccbar::g0(psi(3770),e^+e^-)"] }
                };
                std::array<Parameter, 2> D0Dbar0_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"], p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] }
                };
                std::array<Parameter, 2> DpDm_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^+D^-)"], p["ee->ccbar::g0(psi(3770),D^+D^-)"] }
                };
                std::array<Parameter, 2> eff2S_g0s{
                    { p["ee->ccbar::g0(psi(2S),eff(2S))"], p["ee->ccbar::g0(psi(3770),eff(2S))"] }
                };
                std::array<Parameter, 2> eff3770_g0s{
                    { p["ee->ccbar::g0(psi(2S),eff(3770))"], p["ee->ccbar::g0(psi(3770),eff(3770))"] }
                };

                auto ee      = std::make_shared<EEChannel<5, 2>>("ee", p["mass::e"], p["mass::e"], p["ee->ccbar::q_0"], ee_g0s);
                auto D0Dbar0 = std::make_shared<PWavePPChannel<5, 2>>("D0Dbar0", p["mass::D^0"], p["mass::D^0"], p["ee->ccbar::q_0"], D0Dbar0_g0s);
                auto DpDm    = std::make_shared<PWavePPChannel<5, 2>>("DpDm", p["mass::D^+"], p["mass::D^+"], p["ee->ccbar::q_0"], DpDm_g0s);
                auto eff2S   = std::make_shared<EffChannel<5, 2>>("eff(2S)", p["ee->ccbar::effective_mass"], p["ee->ccbar::effective_mass"], p["ee->ccbar::q_0"], eff2S_g0s);
                auto eff3770 = std::make_shared<EffChannel<5, 2>>("eff(3770)", p["ee->ccbar::effective_mass"], p["ee->ccbar::effective_mass"], p["ee->ccbar::q_0"], eff3770_g0s);

                KMatrix<5, 2> KMatrix32({ ee, D0Dbar0, DpDm, eff2S, eff3770 }, { psi2S, psi3770 }, bkgcst, "ee->ccbar");

                auto tMatrix32s0 = KMatrix32.tmatrix_row(0, 0.0001);
                auto tMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
                auto tMatrix32s2 = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_RELATIVE_ERROR(real(eff2S->chew_mandelstam(20)), -0.0172797, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(eff2S->chew_mandelstam(20)), 0.0189092, eps);
                TEST_CHECK_RELATIVE_ERROR(real(eff3770->chew_mandelstam(20)), -0.0172797, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(eff3770->chew_mandelstam(20)), 0.0189092, eps);

                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].real(), 3.17782e-5, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].imag(), 1.99853e-11, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].real(), -0.0391188, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].imag(), -0.000556894, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].real(), -0.0067422, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].imag(), 0.00187893, eps);
            }

            {
                // Test the full K matrix
                Parameters p                            = Parameters::Defaults();
                p["mass::psi(3770)"]                    = 3.796443282051135;
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]      = -0.02077753547690923;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]    = -0.001999916489715092;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"]   = -3.3693829070086214;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 8.38428874933062;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"]      = 3.5692138012280807;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]    = -8.391199026268827;
                p["ee->ccbar::g0(psi(2S),eff(2S))"]     = 0.09237874214140328;
                p["ee->ccbar::g0(psi(3770),eff(3770))"] = -1.5409316160476978;

                p["ee->ccbar::effective_mass"] = 0.1349768;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                // Masses
                p["mass::psi(2S)"] = 3.68610;
                p["mass::e"]       = 5.10999e-4;
                p["mass::D^0"]     = 1.86483;
                p["mass::D^+"]     = 1.86965;

                Options   oo;
                EEToCCBar c(p, oo);

                auto ir = c.prepare(3.78);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbar0(ir), 0.872205, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDpDm(ir), 0.676807, eps);

                auto irc = c.prepare_complex(3.78, -0.001);

                // Test the Chew Mandelstam function on the first and second Riemann sheets
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_ee(irc), -0.112832, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_ee(irc), -0.019891, eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_II_ee(irc), -0.112832, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_II_ee(irc), 0.0198977, eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_DpDm(irc), 0.000985784, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_DpDm(irc), -0.0006997, eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_II_DpDm(irc), 0.00102809, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_II_DpDm(irc), 0.000664733, eps);

                // Test the amplitude on the first and second Riemann sheets
                TEST_CHECK_RELATIVE_ERROR(c.re_T_eetoDpDm(irc), 0.0275484, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_eetoDpDm(irc), -0.096347, eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_II_eetoDpDm(irc), 0.0238397, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_II_eetoDpDm(irc), 0.111547, eps);

                // Set ee -> DD cst to .5
                p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"] = 0.5;
                p["ee->ccbar::c(e^+e^-,D^+D^-)"]    = 0.5;

                irc = c.prepare_complex(3.78, -0.001);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_eetoDpDm(irc), 0.288606, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_eetoDpDm(irc), -0.0809925, eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_II_eetoDpDm(irc), 0.288325, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_II_eetoDpDm(irc), 0.0890536, eps);

                irc = c.prepare_complex(3.8037, -0.03923);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_eetoDpDm(irc), 0.330701, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_eetoDpDm(irc), -0.0672643, eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_II_eetoDpDm(irc), 0.296099, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_II_eetoDpDm(irc), -0.0682443, eps);
            }

            {
                // The D^* Dbar + h.c. channels, coupled to the resonances.
                Parameters p                             = Parameters::Defaults();
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]       = 10.0;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]     = 12.0;
                p["ee->ccbar::g0(psi(2S),D^*0Dbar^0)"]   = 0.3;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^0)"] = 0.4;
                p["ee->ccbar::g0(psi(2S),D^*+D^-)"]      = 0.5;
                p["ee->ccbar::g0(psi(3770),D^*+D^-)"]    = 0.6;

                // Set all cst to zero
                p["ee->ccbar::c(e^+e^-,e^+e^-)"]         = 0.0;
                p["ee->ccbar::c(e^+e^-,D^*0Dbar^0)"]     = 0.0;
                p["ee->ccbar::c(e^+e^-,D^*+D^-)"]        = 0.0;
                p["ee->ccbar::c(D^*0Dbar^0,D^*0Dbar^0)"] = 0.0;
                p["ee->ccbar::c(D^*0Dbar^0,D^*+D^-)"]    = 0.0;
                p["ee->ccbar::c(D^*+D^-,D^*+D^-)"]       = 0.0;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                // Masses
                p["mass::psi(2S)"]   = 3.68610;
                p["mass::psi(3770)"] = 3.7737;
                p["mass::e"]         = 5.10999e-4;
                p["mass::D^0"]       = 1.86483;
                p["mass::D^+"]       = 1.86965;
                p["mass::D_u^*"]     = 2.00685;
                p["mass::D_d^*"]     = 2.01026;

                auto psi2S   = std::make_shared<CharmoniumResonance<3, 2>>("psi2S", p["mass::psi(2S)"]);
                auto psi3770 = std::make_shared<CharmoniumResonance<3, 2>>("psi3770", p["mass::psi(3770)"]);

                std::array<std::array<Parameter, 3>, 3> bkgcst{
                    p["ee->ccbar::c(e^+e^-,e^+e^-)"],     p["ee->ccbar::c(e^+e^-,D^*0Dbar^0)"],     p["ee->ccbar::c(e^+e^-,D^*+D^-)"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^0)"], p["ee->ccbar::c(D^*0Dbar^0,D^*0Dbar^0)"], p["ee->ccbar::c(D^*0Dbar^0,D^*+D^-)"],
                    p["ee->ccbar::c(e^+e^-,D^*+D^-)"],    p["ee->ccbar::c(D^*0Dbar^0,D^*+D^-)"],    p["ee->ccbar::c(D^*+D^-,D^*+D^-)"],
                };

                std::array<Parameter, 2> ee_g0s{
                    { p["ee->ccbar::g0(psi(2S),e^+e^-)"], p["ee->ccbar::g0(psi(3770),e^+e^-)"] }
                };
                std::array<Parameter, 2> Dst0D0bar_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^*0Dbar^0)"], p["ee->ccbar::g0(psi(3770),D^*0Dbar^0)"] }
                };
                std::array<Parameter, 2> DstpDm_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^*+D^-)"], p["ee->ccbar::g0(psi(3770),D^*+D^-)"] }
                };

                auto ee        = std::make_shared<EEChannel<3, 2>>("ee", p["mass::e"], p["mass::e"], p["ee->ccbar::q_0"], ee_g0s);
                auto Dst0D0bar = std::make_shared<PWaveVPChannel<3, 2>>("Dst0D0bar", p["mass::D_u^*"], p["mass::D^0"], p["ee->ccbar::q_0"], Dst0D0bar_g0s);
                auto DstpDm    = std::make_shared<PWaveVPChannel<3, 2>>("DstpDm", p["mass::D_d^*"], p["mass::D^+"], p["ee->ccbar::q_0"], DstpDm_g0s);

                KMatrix<3, 2> KMatrix32({ ee, Dst0D0bar, DstpDm }, { psi2S, psi3770 }, bkgcst, "ee->ccbar");

                // s = 20 lies above both thresholds, s = 14 below them
                const auto rho_open     = Dst0D0bar->rho(20.0);
                const auto cm_open      = Dst0D0bar->chew_mandelstam(20.0);
                const auto cm_closed    = Dst0D0bar->chew_mandelstam(14.0);
                const auto tmatrix_open = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_RELATIVE_ERROR(std::abs(rho_open), 0.0099522, eps);
                TEST_CHECK_RELATIVE_ERROR(real(cm_open), 0.000793319, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(cm_open), 0.00829489, eps);
                TEST_CHECK_RELATIVE_ERROR(real(cm_closed), -0.00190795, eps);
                TEST_CHECK_RELATIVE_ERROR(tmatrix_open[1].real(), 0.3089342, eps);
                TEST_CHECK_RELATIVE_ERROR(tmatrix_open[1].imag(), 0.0681416, eps);
            }

            {
                // The D^*0 Dbar^*0 channels, coupled to the resonances.
                Parameters p                                   = Parameters::Defaults();
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]             = 10.0;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]           = 12.0;
                p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0(1P1))"]   = 0.3;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(1P1))"] = 0.4;
                p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0(5P1))"]   = 0.5;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(5P1))"] = 0.6;
                p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0(5F1))"]   = 0.7;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(5F1))"] = 0.8;

                // Set all cst to zero
                p["ee->ccbar::c(e^+e^-,e^+e^-)"]                     = 0.0;
                p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(1P1))"]           = 0.0;
                p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(5P1))"]           = 0.0;
                p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(5F1))"]           = 0.0;
                p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(1P1))"] = 0.0;
                p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(5P1))"] = 0.0;
                p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(5F1))"] = 0.0;
                p["ee->ccbar::c(D^*0Dbar^*0(5P1),D^*0Dbar^*0(5P1))"] = 0.0;
                p["ee->ccbar::c(D^*0Dbar^*0(5P1),D^*0Dbar^*0(5F1))"] = 0.0;
                p["ee->ccbar::c(D^*0Dbar^*0(5F1),D^*0Dbar^*0(5F1))"] = 0.0;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                // Masses
                p["mass::psi(2S)"]   = 3.68610;
                p["mass::psi(3770)"] = 3.7737;
                p["mass::e"]         = 5.10999e-4;
                p["mass::D_u^*"]     = 2.00685;

                auto psi2S   = std::make_shared<CharmoniumResonance<4, 2>>("psi2S", p["mass::psi(2S)"]);
                auto psi3770 = std::make_shared<CharmoniumResonance<4, 2>>("psi3770", p["mass::psi(3770)"]);

                std::array<std::array<Parameter, 4>, 4> bkgcst{
                    p["ee->ccbar::c(e^+e^-,e^+e^-)"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(1P1))"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(5P1))"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(5F1))"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(1P1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(1P1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(5P1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(5F1))"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(5P1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(5P1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(5P1),D^*0Dbar^*0(5P1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(5P1),D^*0Dbar^*0(5F1))"],
                    p["ee->ccbar::c(e^+e^-,D^*0Dbar^*0(5F1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(1P1),D^*0Dbar^*0(5F1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(5P1),D^*0Dbar^*0(5F1))"],
                    p["ee->ccbar::c(D^*0Dbar^*0(5F1),D^*0Dbar^*0(5F1))"],
                };

                std::array<Parameter, 2> ee_g0s{
                    { p["ee->ccbar::g0(psi(2S),e^+e^-)"], p["ee->ccbar::g0(psi(3770),e^+e^-)"] }
                };
                std::array<Parameter, 2> Dst0Dst0bar1P1_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0(1P1))"], p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(1P1))"] }
                };
                std::array<Parameter, 2> Dst0Dst0bar5P1_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0(5P1))"], p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(5P1))"] }
                };
                std::array<Parameter, 2> Dst0Dst0bar5F1_g0s{
                    { p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0(5F1))"], p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(5F1))"] }
                };

                auto ee             = std::make_shared<EEChannel<4, 2>>("ee", p["mass::e"], p["mass::e"], p["ee->ccbar::q_0"], ee_g0s);
                auto Dst0Dst0bar1P1 = std::make_shared<PWaveVVChannel<4, 2>>("Dst0Dst0bar1P1", p["mass::D_u^*"], p["mass::D_u^*"], p["ee->ccbar::q_0"], Dst0Dst0bar1P1_g0s);
                auto Dst0Dst0bar5P1 = std::make_shared<PWaveVVChannel<4, 2>>("Dst0Dst0bar5P1", p["mass::D_u^*"], p["mass::D_u^*"], p["ee->ccbar::q_0"], Dst0Dst0bar5P1_g0s);
                auto Dst0Dst0bar5F1 = std::make_shared<FWaveVVChannel<4, 2>>("Dst0Dst0bar5F1", p["mass::D_u^*"], p["mass::D_u^*"], p["ee->ccbar::q_0"], Dst0Dst0bar5F1_g0s);

                KMatrix<4, 2> KMatrix42({ ee, Dst0Dst0bar1P1, Dst0Dst0bar5P1, Dst0Dst0bar5F1 }, { psi2S, psi3770 }, bkgcst, "ee->ccbar");

                // s = 20 lies above both thresholds, s = 14 below them
                const auto rho_open     = Dst0Dst0bar1P1->rho(20.0);
                const auto cm_open_P    = Dst0Dst0bar1P1->chew_mandelstam(20.0);
                const auto cm_closed_P  = Dst0Dst0bar1P1->chew_mandelstam(14.0);
                const auto cm_open_F    = Dst0Dst0bar5F1->chew_mandelstam(20.0);
                const auto cm_closed_F  = Dst0Dst0bar5F1->chew_mandelstam(14.0);
                const auto tmatrix_open = KMatrix42.tmatrix_row(0, 20.0);

                TEST_CHECK_RELATIVE_ERROR(std::abs(rho_open), 0.008774084, eps);
                TEST_CHECK_RELATIVE_ERROR(real(cm_open_P), 0.001314708, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(cm_open_P), 0.006979871, eps);
                TEST_CHECK_RELATIVE_ERROR(real(cm_closed_P), -0.00305087, eps);
                TEST_CHECK_RELATIVE_ERROR(real(cm_open_F), 0.002703345, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(cm_open_F), 0.000939654, eps);
                TEST_CHECK_RELATIVE_ERROR(real(cm_closed_F), -0.000929367, eps);
                TEST_CHECK_RELATIVE_ERROR(tmatrix_open[1].real(), 0.30186, eps);
                TEST_CHECK_RELATIVE_ERROR(tmatrix_open[1].imag(), 0.0665914, eps);
            }

            {
                // Narrow-resonance limit: at the peak, sigma(e^+e^- -> f) = 12 pi / s * B(e^+e^-) * B(f)
                Parameters p                            = Parameters::Defaults();
                p["mass::psi(3770)"]                    = 3.7737;
                p["ee->ccbar::q_0"]                     = 0.5;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]    = 0.05;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 1.0;

                EEToCCBar c(p, Options{});

                // Locate the peak through a ternary search
                double lo = 3.76, hi = 3.79;
                for (unsigned i = 0; i < 200; ++i)
                {
                    const double E1     = lo + (hi - lo) / 3.0;
                    const double E2     = hi - (hi - lo) / 3.0;
                    const double sigma1 = c.sigma_eetoD0Dbar0(c.prepare(E1));
                    const double sigma2 = c.sigma_eetoD0Dbar0(c.prepare(E2));
                    if (sigma1 < sigma2)
                    {
                        lo = E1;
                    }
                    else
                    {
                        hi = E2;
                    }
                }
                const double E = (lo + hi) / 2.0, s = E * E;

                const double m_e      = p["mass::e"].evaluate();
                const double m_D      = p["mass::D^0"].evaluate();
                const double q2       = s / 4.0 - m_D * m_D;
                const double rho_ee   = std::sqrt(1.0 - 4.0 * m_e * m_e / s) / (16.0 * M_PI);
                const double rho_D    = std::sqrt(4.0 * q2 / s) / (16.0 * M_PI);
                const double Gamma_ee = power_of<2>(0.05) * rho_ee;
                // P-wave barrier factor n^2 = q^2 / (q^2 + q_0^2)
                const double Gamma_D  = power_of<2>(1.0) * rho_D * q2 / (q2 + 0.25);
                const double hbarc2   = 0.3893793721e6; // (hbar c)^2 in GeV^2 nb
                const double expected = hbarc2 * 12.0 * M_PI / s * Gamma_ee * Gamma_D / power_of<2>(Gamma_ee + Gamma_D) / (16.0 * M_PI * rho_ee);

                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbar0(c.prepare(E)), expected, eps);
            }

            {
                // An unknown value of assume-symmetry is rejected
                Parameters p = Parameters::Defaults();
                Options    o{
                       { "assume-symmetry"_ok, "su3"_ov }
                };

                TEST_CHECK_THROWS(InvalidOptionValueError, EEToCCBar(p, o));
            }

            {
                // assume-symmetry replaces the couplings of related channels; their own couplings are then ignored
                Parameters p                             = Parameters::Defaults();
                p["ee->ccbar::q_0"]                      = 0.5;
                p["mass::psi(3770)"]                     = 3.7737;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]     = 0.05;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]  = 1.0;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^0)"] = 1.0;

                // cross sections into D^+ D^-, D^*+ D^- + h.c., and D_s^+ D_s^-
                const auto sigmas = [&](const qnp::OptionValue & symmetry)
                {
                    EEToCCBar  c(p,
                                Options{
                                     { "assume-symmetry"_ok, symmetry }
                    });
                    const auto ir = c.prepare(4.3);

                    return std::array<double, 3>{ c.sigma_eetoDpDm(ir), c.sigma_eetoDstpDm(ir), c.sigma_eetoDsDs(ir) };
                };

                const auto none    = sigmas("none"_ov);
                const auto isospin = sigmas("isospin"_ov);
                const auto su3f    = sigmas("su3f"_ov);

                TEST_CHECK_NEARLY_EQUAL(none[0], 0.0, 1.0e-15);
                TEST_CHECK_NEARLY_EQUAL(none[1], 0.0, 1.0e-15);
                TEST_CHECK_NEARLY_EQUAL(none[2], 0.0, 1.0e-15);
                TEST_CHECK(isospin[0] > 0.0);
                TEST_CHECK(isospin[1] > 0.0);
                TEST_CHECK_NEARLY_EQUAL(isospin[2], 0.0, 1.0e-15);
                TEST_CHECK(su3f[0] > 0.0);
                TEST_CHECK(su3f[1] > 0.0);
                TEST_CHECK(su3f[2] > 0.0);

                // coupling the isospin partners changes nothing under either symmetry
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]  = 2.0;
                p["ee->ccbar::g0(psi(3770),D^*+D^-)"] = 2.0;

                const auto none2    = sigmas("none"_ov);
                const auto isospin2 = sigmas("isospin"_ov);
                const auto su3f2    = sigmas("su3f"_ov);

                TEST_CHECK(none2[0] > 0.0);
                TEST_CHECK(none2[1] > 0.0);
                TEST_CHECK_RELATIVE_ERROR(isospin2[0], isospin[0], 1.0e-14);
                TEST_CHECK_RELATIVE_ERROR(isospin2[1], isospin[1], 1.0e-14);
                TEST_CHECK_RELATIVE_ERROR(su3f2[0], su3f[0], 1.0e-14);
                TEST_CHECK_RELATIVE_ERROR(su3f2[1], su3f[1], 1.0e-14);

                // coupling D_s^+ D_s^- opens it under isospin, and changes nothing under SU(3)_F
                p["ee->ccbar::g0(psi(3770),D_s^+D_s^-)"] = 2.0;

                const auto isospin3 = sigmas("isospin"_ov);
                const auto su3f3    = sigmas("su3f"_ov);

                TEST_CHECK(isospin3[2] > 0.0);
                TEST_CHECK_RELATIVE_ERROR(su3f3[0], su3f[0], 1.0e-14);
                TEST_CHECK_RELATIVE_ERROR(su3f3[1], su3f[1], 1.0e-14);
                TEST_CHECK_RELATIVE_ERROR(su3f3[2], su3f[2], 1.0e-14);
            }

            {
                // The D_s^+ D_s^- channel is closed below threshold
                Parameters p                             = Parameters::Defaults();
                p["ee->ccbar::q_0"]                      = 0.5;
                p["mass::D_s"]                           = 1.96835;
                p["mass::psi(3770)"]                     = 3.7737;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]     = 0.05;
                p["ee->ccbar::g0(psi(3770),D_s^+D_s^-)"] = 1.0;

                EEToCCBar c(p, Options{});

                const double threshold = 2.0 * 1.96835;
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDsDs(c.prepare(threshold - 0.01)), 0.0, 1.0e-15);
                TEST_CHECK(c.sigma_eetoDsDs(c.prepare(threshold + 0.01)) > 0.0);
            }

            {
                // VV helicity cross sections: each LS channel alone has fixed helicity fractions TT : TL : LL
                Parameters p                         = Parameters::Defaults();
                p["ee->ccbar::q_0"]                  = 0.5;
                p["mass::D_u^*"]                     = 2.00685;
                p["mass::D_d^*"]                     = 2.01026;
                p["mass::psi(3770)"]                 = 3.7737;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"] = 0.05;

                const std::array<std::string, 3>           waves{ "1P1", "5P1", "5F1" };
                const std::array<std::array<double, 3>, 3> fractions{
                    { { 2.0 / 3.0, 0.0, 1.0 / 3.0 }, { 2.0 / 15.0, 3.0 / 5.0, 4.0 / 15.0 }, { 1.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0 } }
                };

                for (unsigned i = 0; i < waves.size(); ++i)
                {
                    for (unsigned j = 0; j < waves.size(); ++j)
                    {
                        p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(" + waves[j] + "))"] = (i == j) ? 1.0 : 0.0;
                    }

                    EEToCCBar  c(p, Options{});
                    const auto ir    = c.prepare(4.3);
                    const auto total = c.sigma_eetoDst0Dst0bar(ir);

                    TEST_CHECK(total > 0.0);
                    TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDst0Dst0bar_TT(ir) / total, fractions[i][0], 1.0e-12);
                    TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDst0Dst0bar_TL(ir) / total, fractions[i][1], 1.0e-12);
                    TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDst0Dst0bar_LL(ir) / total, fractions[i][2], 1.0e-12);
                }

                // with all LS channels coupled, the helicity cross sections still sum to the total
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(1P1))"] = 1.0;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(5P1))"] = 0.5;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0(5F1))"] = 0.7;

                for (const auto & symmetry : { "none"_ov, "isospin"_ov })
                {
                    EEToCCBar  c(p,
                                Options{
                                     { "assume-symmetry"_ok, symmetry }
                    });
                    const auto ir = c.prepare(4.3);

                    const auto total0 = c.sigma_eetoDst0Dst0bar(ir);
                    TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0Dst0bar_1P1(ir) + c.sigma_eetoDst0Dst0bar_5P1(ir) + c.sigma_eetoDst0Dst0bar_5F1(ir), total0, 1.0e-12);
                    TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0Dst0bar_TT(ir) + c.sigma_eetoDst0Dst0bar_TL(ir) + c.sigma_eetoDst0Dst0bar_LL(ir), total0, 1.0e-12);

                    // the charged channels have no couplings of their own, but inherit the neutral ones under isospin
                    const auto totalp = c.sigma_eetoDstpDstm(ir);
                    if (symmetry == "none"_ov)
                    {
                        TEST_CHECK_NEARLY_EQUAL(totalp, 0.0, 1.0e-15);
                    }
                    else
                    {
                        TEST_CHECK(totalp > 0.0);
                        TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstpDstm_1P1(ir) + c.sigma_eetoDstpDstm_5P1(ir) + c.sigma_eetoDstpDstm_5F1(ir), totalp, 1.0e-12);
                        TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstpDstm_TT(ir) + c.sigma_eetoDstpDstm_TL(ir) + c.sigma_eetoDstpDstm_LL(ir), totalp, 1.0e-12);
                    }
                }
            }

            {
                // The D_s^*+ D_s^- + h.c. channel: closed below threshold, and its SU(3)_F partner D^*0 Dbar^0 + h.c.
                Parameters p                             = Parameters::Defaults();
                p["ee->ccbar::q_0"]                      = 0.5;
                p["mass::D_s"]                           = 1.96835;
                p["mass::D_s^*"]                         = 2.1122;
                p["mass::psi(3770)"]                     = 3.7737;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]     = 0.05;
                p["ee->ccbar::g0(psi(3770),D^*0Dbar^0)"] = 1.0;

                EEToCCBar none(p, Options{});
                EEToCCBar su3f(p,
                               Options{
                                   { "assume-symmetry"_ok, "su3f"_ov }
                });

                const double threshold = 2.1122 + 1.96835;
                TEST_CHECK_NEARLY_EQUAL(none.sigma_eetoDstsDs(none.prepare(4.3)), 0.0, 1.0e-15);
                TEST_CHECK_NEARLY_EQUAL(su3f.sigma_eetoDstsDs(su3f.prepare(threshold - 0.01)), 0.0, 1.0e-15);
                TEST_CHECK(su3f.sigma_eetoDstsDs(su3f.prepare(threshold + 0.01)) > 0.0);

                // its own coupling is ignored under SU(3)_F
                const double sigma                        = su3f.sigma_eetoDstsDs(su3f.prepare(4.3));
                p["ee->ccbar::g0(psi(3770),D_s^*+D_s^-)"] = 2.0;
                TEST_CHECK_RELATIVE_ERROR(su3f.sigma_eetoDstsDs(su3f.prepare(4.3)), sigma, 1.0e-14);
                TEST_CHECK(none.sigma_eetoDstsDs(none.prepare(4.3)) > 0.0);
            }

            {
                // Narrow-resonance limit for psi(4040) and psi(4160) decaying into their effective channels, and their widths
                struct Case
                {
                        std::string resonance, channel;
                        double      mass;
                        double (EEToCCBar::*eff_width)() const;
                        double (EEToCCBar::*total_width)() const;
                };

                for (const auto & [resonance, channel, M, eff_width, total_width] : {
                         Case{ "psi(4040)", "eff(4040)", 4.040, &EEToCCBar::psi4040_eff_width, &EEToCCBar::psi4040_total_width },
                         Case{ "psi(4160)", "eff(4160)", 4.191, &EEToCCBar::psi4160_eff_width, &EEToCCBar::psi4160_total_width }
                })
                {
                    Parameters p                                          = Parameters::Defaults();
                    p["mass::" + resonance]                               = M;
                    p["ee->ccbar::q_0"]                                   = 0.5;
                    p["ee->ccbar::effective_mass"]                        = 0.1;
                    p["ee->ccbar::g0(" + resonance + ",e^+e^-)"]          = 0.05;
                    p["ee->ccbar::g0(" + resonance + "," + channel + ")"] = 1.0;

                    EEToCCBar c(p, Options{});

                    // Locate the peak through a ternary search
                    double lo = M - 0.03, hi = M + 0.03;
                    for (unsigned i = 0; i < 200; ++i)
                    {
                        const double E1     = lo + (hi - lo) / 3.0;
                        const double E2     = hi - (hi - lo) / 3.0;
                        const double sigma1 = c.sigma_eetoeff(c.prepare(E1));
                        const double sigma2 = c.sigma_eetoeff(c.prepare(E2));
                        if (sigma1 < sigma2)
                        {
                            lo = E1;
                        }
                        else
                        {
                            hi = E2;
                        }
                    }
                    const double E = (lo + hi) / 2.0, s = E * E;

                    const double m_e       = p["mass::e"].evaluate();
                    const double q2        = s / 4.0 - 0.01;
                    const double rho_ee    = std::sqrt(1.0 - 4.0 * m_e * m_e / s) / (16.0 * M_PI);
                    const double rho_eff   = std::sqrt(4.0 * q2 / s) / (16.0 * M_PI);
                    const double Gamma_ee  = power_of<2>(0.05) * rho_ee;
                    // P-wave barrier factor n^2 = q^2 / (q^2 + q_0^2)
                    const double Gamma_eff = power_of<2>(1.0) * rho_eff * q2 / (q2 + 0.25);
                    const double hbarc2    = 0.3893793721e6; // (hbar c)^2 in GeV^2 nb
                    const double expected  = hbarc2 * 12.0 * M_PI / s * Gamma_ee * Gamma_eff / power_of<2>(Gamma_ee + Gamma_eff) / (16.0 * M_PI * rho_ee);

                    TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoeff(c.prepare(E)), expected, eps);

                    // partial widths g0^2 rho(M^2) / M, without barrier factor
                    const double width_ee  = power_of<2>(0.05) * std::sqrt(1.0 - 4.0 * m_e * m_e / (M * M)) / (16.0 * M_PI) / M;
                    const double width_eff = power_of<2>(1.0) * std::sqrt(1.0 - 4.0 * 0.01 / (M * M)) / (16.0 * M_PI) / M;
                    TEST_CHECK_RELATIVE_ERROR((c.*eff_width)(), width_eff, eps);
                    TEST_CHECK_RELATIVE_ERROR((c.*total_width)(), width_ee + width_eff, eps);
                }
            }

            {
                // Coupled-channel cross sections above all open-charm thresholds. The reference values stem from an
                // independent evaluation that computes the Chew-Mandelstam functions as dispersion integrals.
                Parameters p                                   = Parameters::Defaults();
                p["ee->ccbar::q_0"]                            = 0.5;
                p["ee->ccbar::exclusive_norm"]                 = 1.0;
                p["ee->ccbar::Rconstant"]                      = 0.0;
                p["mass::e"]                                   = 5.10999e-4;
                p["mass::D^0"]                                 = 1.86483;
                p["mass::D^+"]                                 = 1.86965;
                p["mass::D_u^*"]                               = 2.00685;
                p["mass::D_d^*"]                               = 2.01026;
                p["mass::D_s"]                                 = 1.96835;
                p["mass::D_s^*"]                               = 2.1122;
                p["mass::psi(4040)"]                           = 4.040;
                p["mass::psi(4160)"]                           = 4.191;
                p["ee->ccbar::g0(psi(4040),e^+e^-)"]           = 0.05;
                p["ee->ccbar::g0(psi(4040),D^*0Dbar^0)"]       = 1.0;
                p["ee->ccbar::g0(psi(4040),D^*+D^-)"]          = 0.8;
                p["ee->ccbar::g0(psi(4040),D_s^+D_s^-)"]       = 0.6;
                p["ee->ccbar::g0(psi(4040),D^*+D^*-(1P1))"]    = 0.3;
                p["ee->ccbar::g0(psi(4040),D^*+D^*-(5F1))"]    = 0.6;
                p["ee->ccbar::g0(psi(4160),e^+e^-)"]           = 0.04;
                p["ee->ccbar::g0(psi(4160),D^*0Dbar^0)"]       = 0.5;
                p["ee->ccbar::g0(psi(4160),D_s^+D_s^-)"]       = 0.3;
                p["ee->ccbar::g0(psi(4160),D_s^*+D_s^-)"]      = 0.9;
                p["ee->ccbar::g0(psi(4160),D^*0Dbar^*0(1P1))"] = 0.7;
                p["ee->ccbar::g0(psi(4160),D^*0Dbar^*0(5P1))"] = 0.4;
                p["ee->ccbar::g0(psi(4160),D^*0Dbar^*0(5F1))"] = 0.5;

                EEToCCBar c(p, Options{});

                const auto ir = c.prepare(4.2);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0D0bar(ir), 7.786533313, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstpDm(ir), 0.08112239587, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDsDs(ir), 2.368817589, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstsDs(ir), 8.011028937, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0Dst0bar(ir), 9.807427087, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0Dst0bar_TT(ir), 3.119090068, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0Dst0bar_TL(ir), 1.871259624, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDst0Dst0bar_LL(ir), 4.817077395, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstpDstm(ir), 0.007626808885, eps);

                // R sums all channels, unaffected by the normalisation of the exclusive cross sections
                p["ee->ccbar::exclusive_norm"] = 2.0;
                p["ee->ccbar::Rconstant"]      = 0.3;

                const auto   ir2 = c.prepare(4.2);
                const double sum = c.sigma_eetoeff(ir2) + c.sigma_eetoD0Dbar0(ir2) + c.sigma_eetoDpDm(ir2) + c.sigma_eetoDst0D0bar(ir2) + c.sigma_eetoDstpDm(ir2)
                                   + c.sigma_eetoDsDs(ir2) + c.sigma_eetoDst0Dst0bar(ir2) + c.sigma_eetoDstpDstm(ir2) + c.sigma_eetoDstsDs(ir2);
                const double GeVtonb    = 10.0 * power_of<2>(1.0e18 * p["QM::hbar"].evaluate() * 299792458.0);
                const double sigma_mumu = GeVtonb * 4.0 * M_PI * power_of<2>(p["QED::alpha_e(0)"].evaluate()) / (3.0 * 4.2 * 4.2);
                TEST_CHECK_RELATIVE_ERROR(c.R(ir2), sum / 2.0 / sigma_mumu + 0.3, 1.0e-12);
            }
        }
} eetoccbar_test;
