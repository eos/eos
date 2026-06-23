/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/maths/complex.hh>
#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/utils/kmatrix-impl.hh>

#include <array>
#include <cmath>
#include <memory>
#include <vector>


using namespace test;
using namespace eos;

class eetoccbarTest :
    public TestCase
{
    public:
    eetoccbarTest() :
        TestCase("ee->ccbar scattering tests")
        {
        }

        virtual void run() const
        {

            constexpr double eps = 1e-5;

            {
                // Standalone validation of the unequal-mass S-wave PV channel (D Dbar^*).
                // Reference values were computed independently with mpmath and the closed
                // form was checked against the equal-mass limit (EEChannel) and a numerical
                // once-subtracted dispersion integral of rho; see PHASE2-HARD-SCOPE.md and
                // analytic/chew-mandelstam.py. CM is threshold-subtracted (CM(s_th) = 0),
                // matching EEChannel/PWavePPChannel; Im CM is unchanged by the subtraction.
                Parameters p = Parameters::Defaults();
                p["ee->ccbar::q_0"] = 0.5;
                std::array<Parameter, 2> g0s {{ p["ee->ccbar::q_0"], p["ee->ccbar::q_0"] }};

                // D^0 Dbar^*0: m1 = 1.86483, m2 = 2.00685 GeV, threshold s = 14.989 GeV^2
                auto ch = std::make_shared<SWavePVChannel<3, 2>>("D^0Dbar^*0", p["mass::D^0"], p["mass::D_u^*"], p["ee->ccbar::q_0"], g0s);

                // rho: zero below threshold, sqrt(lambda)/(16 pi s) above
                TEST_CHECK_NEARLY_EQUAL(  std::abs(ch->rho(14.0)),    0.0,            eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(16.5)),        0.006014846605, eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(25.0)),        0.01258357249,  eps);

                // Chew-Mandelstam: real below threshold; Im = rho above threshold
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(14.0)), -0.004412244375, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(ch->chew_mandelstam(14.0)),  0.0,            eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(16.5)), -0.001195053117, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(16.5)),  0.006014846605, eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(25.0)), -0.005972998193, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(25.0)),  0.01258357249,  eps);

                // complex s = 16 - 0.1 i
                const complex<double> sc(16.0, -0.1);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(sc)), -0.001047685124, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(sc)), -0.004924905559, eps);
            }

            {
                // High-accuracy validation of the S-wave (l = 0) Chew-Mandelstam function
                // against the reference values obtained from a Python script using 40-digit
                // mpmath. See EOS-ANALYSIS-2026-03 for the values.
                // Three mass configurations are covered:
                // - P1 equal masses 1/sqrt(2) (the EEChannel limit),
                // - P2 equal masses 1, and
                // - P3 the unequal pair (1/sqrt(2), sqrt(2)).
                // The points span negative s, the pseudothreshold, the threshold cusp and
                // the asymptotic region. CM is real below threshold (Im CM = 0), Im CM = rho
                // above threshold (unitarity), and CM(s_th) = 0 (threshold subtraction).
                Parameters p = Parameters::Defaults();
                p["ee->ccbar::q_0"] = 0.5;
                std::array<Parameter, 2> g0s {{ p["ee->ccbar::q_0"], p["ee->ccbar::q_0"] }};

                const double eps_rel  = 1.0e-10; // Re CM, and Im CM above threshold
                const double eps_cusp = 1.0e-8;  // within ~0.15 of the sqrt cusp at s_th
                const double eps_abs  = 1.0e-10; // Im CM below threshold (analytically 0)

                // each row: { s, Re CM(s), Im CM(s) }. m1, m2 reuse two distinct mass
                // parameters set to the exact double literals used to compute the oracle.
                auto check_pair = [&](double m1v, double m2v,
                        const std::vector<std::array<double, 3>> & rows)
                {
                    p["mass::D^0"] = m1v;
                    p["mass::D^+"] = m2v;
                    auto ch = std::make_shared<SWavePVChannel<3, 2>>("CM-ref",
                            p["mass::D^0"], p["mass::D^+"], p["ee->ccbar::q_0"], g0s);
                    const double sth = power_of<2>(m1v + m2v);

                    for (const auto & row : rows)
                    {
                        const double s = row[0], re_ref = row[1], im_ref = row[2];
                        const complex<double> cm = ch->chew_mandelstam(s);
                        const double eps = (std::abs(s - sth) < 0.15) ? eps_cusp : eps_rel;

                        TEST_CHECK_RELATIVE_ERROR(real(cm), re_ref, eps);
                        if (s < sth)
                        {
                            // below threshold Im CM and rho are analytically zero
                            TEST_CHECK_NEARLY_EQUAL(imag(cm),            0.0, eps_abs);
                            TEST_CHECK_NEARLY_EQUAL(std::abs(ch->rho(s)), 0.0, eps_abs);
                        }
                        else
                        {
                            // above threshold unitarity fixes Im CM = rho
                            TEST_CHECK_RELATIVE_ERROR(imag(cm), im_ref,          eps);
                            TEST_CHECK_RELATIVE_ERROR(imag(cm), real(ch->rho(s)), eps);
                        }
                    }

                    // threshold subtraction: CM(s_th) = 0 by construction. The exact
                    // threshold sits on the sqrt cusp, where the +i*1e-15 prescription leaves
                    // a residue ~3e-10, so this uses the looser cusp tolerance (cf. §5).
                    TEST_CHECK_NEARLY_EQUAL(real(ch->chew_mandelstam(sth)), 0.0, eps_cusp);
                    TEST_CHECK_NEARLY_EQUAL(imag(ch->chew_mandelstam(sth)), 0.0, eps_cusp);

                    return ch;
                };

                // P1: m1 = m2 = 1/sqrt(2)
                auto ch_P1 = check_pair(1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), {
                    {       -5.0,   -0.01856633106339601,                    0.0 },
                    {    -3.3273,   -0.01718330725143701,                    0.0 },
                    {    -1.6545,   -0.01535843245361532,                    0.0 },
                    {     0.0182,   -0.01262658976809803,                    0.0 },
                    {     1.2727,  -0.008841779331811922,                    0.0 },
                    {        1.9,  -0.003908832816629901,                    0.0 },
                    {        2.1, -0.0006129585132922205,   0.004341306987767856 },
                    {     6.7719,   -0.01296294970584483,    0.01670015659585684 },
                    {    13.3125,   -0.01869243450631893,    0.01833916805856752 },
                    {    20.7875,   -0.02213977623780263,    0.01891313446938620 },
                    {    28.2625,   -0.02440736344346382,    0.01917753966534671 },
                    {       32.0,   -0.02530388085671864,    0.01926263887692187 },
                });

                // P2: m1 = m2 = 1
                check_pair(1.0, 1.0, {
                    {       -5.0,   -0.01635357875576924,                    0.0 },
                    {    -2.8424,   -0.01504971686221515,                    0.0 },
                    {    -0.6848,   -0.01334276145479074,                    0.0 },
                    {     1.4727,   -0.01081749037523438,                    0.0 },
                    {     3.0909,  -0.007376014124534837,                    0.0 },
                    {        3.9,  -0.002863634700957743,                    0.0 },
                    {        4.1, -0.0003114548843090558,   0.003106978273228536 },
                    {     8.7719,  -0.008828618422788801,    0.01467333982471145 },
                    {    15.3125,   -0.01405927798112908,    0.01709961726772085 },
                    {    22.7875,   -0.01743411331848962,    0.01806409882882630 },
                    {    30.2625,   -0.01970511522029685,    0.01853300215040734 },
                    {       34.0,   -0.02060982647817307,    0.01868750463945365 },
                });

                // P3: m1 = 1/sqrt(2), m2 = sqrt(2) (the genuine unequal-mass branch,
                // pseudothreshold at s_- = 0.5)
                check_pair(1.0 / std::sqrt(2.0), std::sqrt(2.0), {
                    {       -5.0,   -0.01549899689444485,                    0.0 },
                    {    -2.7212,   -0.01423049248800826,                    0.0 },
                    {    -0.4424,   -0.01257581665976170,                    0.0 },
                    {     0.4121,   -0.01179276938694260,                    0.0 },
                    {     1.8364,   -0.01014143015504747,                    0.0 },
                    {     3.5455,  -0.006852433318057166,                    0.0 },
                    {        4.4,  -0.002571474740987491,                    0.0 },
                    {        4.6, -0.0002457938041429069,   0.002769263243529780 },
                    {     9.2719,  -0.007599502018633736,    0.01388205972976265 },
                    {    15.8125,   -0.01256796688621726,    0.01655891790747672 },
                    {    23.2875,   -0.01587938608409122,    0.01767624914919503 },
                    {    30.7625,   -0.01813573212282312,    0.01823177497201795 },
                    {       34.5,   -0.01903884180507123,    0.01841667123888186 },
                });

                // Equal masses reduce exactly to the equal-mass EEChannel: cross-check P1
                // directly against an EEChannel built with the same masses (independent of
                // the oracle table above). ch_P1 holds live references to both mass
                // parameters, which the P2/P3 calls mutated, so restore the P1 masses first.
                p["mass::D^0"] = 1.0 / std::sqrt(2.0);
                p["mass::D^+"] = 1.0 / std::sqrt(2.0);
                auto ee = std::make_shared<EEChannel<3, 2>>("ee-ref",
                        p["mass::D^0"], p["mass::D^0"], p["ee->ccbar::q_0"], g0s);
                for (const double s : { -3.3273, 1.2727, 6.7719, 20.7875 })
                {
                    TEST_CHECK_RELATIVE_ERROR(real(ch_P1->chew_mandelstam(s)),
                            real(ee->chew_mandelstam(s)), 1.0e-10);
                    TEST_CHECK_NEARLY_EQUAL(imag(ch_P1->chew_mandelstam(s)),
                            imag(ee->chew_mandelstam(s)), 1.0e-10);
                }
            }

            {
                // Standalone validation of the unequal-mass D-wave (l = 2) PV channel.
                // Reference values were computed independently with mpmath; the closed form
                // (l = 0 dispersive function plus the four Blatt-Weisskopf poles, threshold-
                // subtracted) was checked against a numerical once-subtracted dispersion
                // integral of rho*n^2 and against EOS's own l = 1 form; see
                // analytic/chew-mandelstam.py (checks D0-D5).
                Parameters p = Parameters::Defaults();
                p["ee->ccbar::q_0"] = 0.5;
                std::array<Parameter, 2> g0s {{ p["ee->ccbar::q_0"], p["ee->ccbar::q_0"] }};

                // D^0 Dbar^*0 masses (a PV pair): threshold s = 14.989 GeV^2, q0 = 0.5 GeV
                auto ch = std::make_shared<DWavePVChannel<3, 2>>("D^0Dbar^*0 (D-wave)", p["mass::D^0"], p["mass::D_u^*"], p["ee->ccbar::q_0"], g0s);

                // rho is the same bare phase space as the S-wave channel
                TEST_CHECK_NEARLY_EQUAL(  std::abs(ch->rho(14.0)),    0.0,            eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(18.0)),        0.008130934478, eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(28.0)),        0.01355610326,  eps);

                // Chew-Mandelstam: real below threshold; Im = rho * n^2 above threshold
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(14.0)), -0.0008931396197, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(ch->chew_mandelstam(14.0)),  0.0,             eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(18.0)),  0.003115315952,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(18.0)),  0.002716379008,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(28.0)),  0.001411672509,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(28.0)),  0.01055767226,   eps);

                // complex s = 18 - 0.1 i
                const complex<double> sc(18.0, -0.1);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(sc)),  0.0029893857,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(sc)), -0.002760542364, eps);

                // threshold subtraction: CM(s_th) = 0 by construction
                const double sth = power_of<2>(double(p["mass::D^0"]) + double(p["mass::D_u^*"]));
                TEST_CHECK_NEARLY_EQUAL(real(ch->chew_mandelstam(sth)), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ch->chew_mandelstam(sth)), 0.0, eps);
            }

            {
                // Standalone validation of the equal-mass F-wave (l = 3) channel (D*Dbar*).
                // Reference values were computed independently with mpmath; the closed form
                // (l = 0 dispersive function plus the six Blatt-Weisskopf poles, threshold-
                // subtracted) was checked against a numerical once-subtracted dispersion
                // integral of rho*n^2; see analytic/chew-mandelstam.py (checks F0-F5).
                Parameters p = Parameters::Defaults();
                p["ee->ccbar::q_0"] = 0.5;
                std::array<Parameter, 2> g0s {{ p["ee->ccbar::q_0"], p["ee->ccbar::q_0"] }};

                // D^*0 Dbar^*0 masses (equal-mass VV pair): threshold s = 16.110 GeV^2, q0 = 0.5 GeV
                auto ch = std::make_shared<FWavePPChannel<3, 2>>("D^*0Dbar^*0 (F-wave)", p["mass::D_u^*"], p["mass::D_u^*"], p["ee->ccbar::q_0"], g0s);

                // rho is the bare equal-mass phase space
                TEST_CHECK_NEARLY_EQUAL(  std::abs(ch->rho(14.0)),    0.0,            eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(17.0)),        0.004552526492, eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(25.0)),        0.01186359211,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->rho(35.0)),        0.0146155291,   eps);

                // Chew-Mandelstam: real below threshold; Im = rho * n^2 above threshold
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(14.0)), -0.0009293672174, eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(ch->chew_mandelstam(14.0)),  0.0,             eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(17.0)),  0.0005132195813, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(17.0)),  1.1872278e-5,    eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(25.0)),  0.004436822172,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(25.0)),  0.004626120868,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(35.0)),  0.003104964563,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(35.0)),  0.009894674852,  eps);

                // complex s = 25 - 0.1 i
                const complex<double> sc(25.0, -0.1);
                TEST_CHECK_RELATIVE_ERROR(real(ch->chew_mandelstam(sc)),  0.00436416736,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ch->chew_mandelstam(sc)), -0.004632750481,  eps);

                // threshold subtraction: CM(s_th) = 0 by construction
                const double sth = power_of<2>(2.0 * double(p["mass::D_u^*"]));
                TEST_CHECK_NEARLY_EQUAL(real(ch->chew_mandelstam(sth)), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ch->chew_mandelstam(sth)), 0.0, eps);
            }

            {
                Parameters p = Parameters::Defaults();
                p["ee->ccbar::g0(psi(2S),e^+e^-)"] = 10.0;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"] = 12.0;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"] = 0.3;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 0.4;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"] = 0.5;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"] = 0.6;

                // Set all cst to zero
                p["ee->ccbar::c(e^+e^-,e^+e^-)"] = 0.0;
                p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"] = 0.0;
                p["ee->ccbar::c(e^+e^-,D^+D^-)"] = 0.0;
                p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"] = 0.0;
                p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"] = 0.0;
                p["ee->ccbar::c(D^+D^-,D^+D^-)"] = 0.0;

                // Channels effective momentum
                p["ee->ccbar::q_0"]    = 0.5;

                // Build K Matrix
                auto psi2S = std::make_shared<CharmoniumResonance<3, 2>>("psi2S", p["mass::psi(2S)"]);
                auto psi3770 = std::make_shared<CharmoniumResonance<3, 2>>("psi3770", p["mass::psi(3770)"]);

                std::array<std::array<Parameter, 3>, 3> bkgcst {
                    p["ee->ccbar::c(e^+e^-,e^+e^-)"],    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"],    p["ee->ccbar::c(e^+e^-,D^+D^-)"],
                    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],
                    p["ee->ccbar::c(e^+e^-,D^+D^-)"],    p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],    p["ee->ccbar::c(D^+D^-,D^+D^-)"],
                };

                std::array<Parameter, 2> ee_g0s       {{p["ee->ccbar::g0(psi(2S),e^+e^-)"],       p["ee->ccbar::g0(psi(3770),e^+e^-)"]}};
                std::array<Parameter, 2> D0Dbar0_g0s  {{p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"],    p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]}};
                std::array<Parameter, 2> DpDm_g0s     {{p["ee->ccbar::g0(psi(2S),D^+D^-)"],       p["ee->ccbar::g0(psi(3770),D^+D^-)"]}};

                auto ee        = std::make_shared<EEChannel<3, 2>>("ee", p["mass::e"], p["mass::e"], p["ee->ccbar::q_0"], ee_g0s);
                auto D0Dbar0   = std::make_shared<PWavePPChannel<3, 2>>("D0Dbar0", p["mass::D^0"], p["mass::D^0"], p["ee->ccbar::q_0"], D0Dbar0_g0s);
                auto DpDm      = std::make_shared<PWavePPChannel<3, 2>>("DpDm", p["mass::D^+"], p["mass::D^+"], p["ee->ccbar::q_0"], DpDm_g0s);

                KMatrix<3, 2> KMatrix32({ee, D0Dbar0, DpDm}, {psi2S, psi3770}, bkgcst, "ee->ccbar");

                TEST_CHECK_RELATIVE_ERROR(std::abs(ee->rho(20)),    0.0198944,  eps);
                TEST_CHECK_RELATIVE_ERROR(std::abs(DpDm->rho(20)),  0.0109126,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->rho(complex<double>(20, 5))),  0.0119492,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(DpDm->rho(complex<double>(20, 5))),  0.00272429,  eps);

                TEST_CHECK_RELATIVE_ERROR(real(ee->chew_mandelstam(20)),    -0.1149620,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(ee->chew_mandelstam(20)),     0.0198944,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->chew_mandelstam(-20)), -0.0125291,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->chew_mandelstam(2)),   -0.0081324,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(DpDm->chew_mandelstam(20)),   0.00024738, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(DpDm->chew_mandelstam(20)),   0.0093575,  eps);

                auto tMatrix32s0 = KMatrix32.tmatrix_row(0,  0.0001);
                auto tMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
                auto tMatrix32s2 = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].real(),  10.1176,     eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].imag(),  2.11433,     eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].real(),  0.0534574,   eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].imag(),  0.00902802,  eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].real(),  0.313608,    eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].imag(),  0.0691847,   eps);


                // Test the implentation of background constants
                p["ee->ccbar::c(e^+e^-,D^+D^-)"] = 0.5;

                tMatrix32s0 = KMatrix32.tmatrix_row(0,  0.0001);
                tMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
                tMatrix32s2 = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_NEARLY_EQUAL(tMatrix32s0[0].real(),  10.1145,    eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s0[0].imag(),  2.11295,    eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s1[1].real(),  0.0534576,  eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s1[1].imag(),  0.00902805, eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s2[1].real(),  0.313633,   eps);
                TEST_CHECK_NEARLY_EQUAL(tMatrix32s2[1].imag(),  0.0690969,  eps);
            }

            {
                Parameters p = Parameters::Defaults();
                p["mass::psi(3770)"]                      =  3.796443282051135;
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]        = -0.02077753547690923;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]      = -0.001999916489715092;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"]     = -3.3693829070086214;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]   =  8.38428874933062;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"]        =  3.5692138012280807;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]      = -8.391199026268827;
                p["ee->ccbar::g0(psi(2S),eff(2S))"]       =  0.09237874214140328;
                p["ee->ccbar::g0(psi(3770),eff(3770))"]   = -1.5409316160476978;

                p["ee->ccbar::effective_mass"] = 0.1349768;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                // Build K Matrix
                auto psi2S = std::make_shared<CharmoniumResonance<5, 2>>("psi2S", p["mass::psi(2S)"]);
                auto psi3770 = std::make_shared<CharmoniumResonance<5, 2>>("psi3770", p["mass::psi(3770)"]);

                std::array<std::array<Parameter, 5>, 5> bkgcst {
                    p["ee->ccbar::c(e^+e^-,e^+e^-)"],    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"],     p["ee->ccbar::c(e^+e^-,D^+D^-)"],      p["ee->ccbar::c(e^+e^-,eff(2S))"],    p["ee->ccbar::c(e^+e^-,eff(3770))"],
                    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"],  p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],   p["ee->ccbar::c(D^0Dbar^0,eff(2S))"], p["ee->ccbar::c(D^0Dbar^0,eff(3770))"],
                    p["ee->ccbar::c(e^+e^-,D^+D^-)"],    p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],     p["ee->ccbar::c(D^+D^-,D^+D^-)"],      p["ee->ccbar::c(D^+D^-,eff(2S))"],    p["ee->ccbar::c(D^+D^-,eff(3770))"],
                    p["ee->ccbar::c(e^+e^-,eff(2S))"],   p["ee->ccbar::c(D^0Dbar^0,eff(2S))"],    p["ee->ccbar::c(D^+D^-,eff(2S))"],     p["ee->ccbar::c(eff(2S),eff(2S))"],   p["ee->ccbar::c(eff(2S),eff(3770))"],
                    p["ee->ccbar::c(e^+e^-,eff(3770))"], p["ee->ccbar::c(D^0Dbar^0,eff(3770))"],  p["ee->ccbar::c(D^+D^-,eff(3770))"],   p["ee->ccbar::c(eff(2S),eff(3770))"], p["ee->ccbar::c(eff(3770),eff(3770))"],
                };

                std::array<Parameter, 2> ee_g0s       {{p["ee->ccbar::g0(psi(2S),e^+e^-)"],       p["ee->ccbar::g0(psi(3770),e^+e^-)"]}};
                std::array<Parameter, 2> D0Dbar0_g0s  {{p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"],    p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]}};
                std::array<Parameter, 2> DpDm_g0s     {{p["ee->ccbar::g0(psi(2S),D^+D^-)"],       p["ee->ccbar::g0(psi(3770),D^+D^-)"]}};
                std::array<Parameter, 2> eff2S_g0s    {{p["ee->ccbar::g0(psi(2S),eff(2S))"],      p["ee->ccbar::g0(psi(3770),eff(2S))"]}};
                std::array<Parameter, 2> eff3770_g0s  {{p["ee->ccbar::g0(psi(2S),eff(3770))"],    p["ee->ccbar::g0(psi(3770),eff(3770))"]}};

                auto ee        = std::make_shared<EEChannel<5, 2>>("ee", p["mass::e"], p["mass::e"], p["ee->ccbar::q_0"], ee_g0s);
                auto D0Dbar0   = std::make_shared<PWavePPChannel<5, 2>>("D0Dbar0", p["mass::D^0"], p["mass::D^0"], p["ee->ccbar::q_0"], D0Dbar0_g0s);
                auto DpDm      = std::make_shared<PWavePPChannel<5, 2>>("DpDm", p["mass::D^+"], p["mass::D^+"], p["ee->ccbar::q_0"], DpDm_g0s);
                auto eff2S     = std::make_shared<EffChannel<5, 2>>("eff(2S)", p["ee->ccbar::effective_mass"], p["ee->ccbar::effective_mass"], p["ee->ccbar::q_0"], eff2S_g0s);
                auto eff3770   = std::make_shared<EffChannel<5, 2>>("eff(3770)", p["ee->ccbar::effective_mass"], p["ee->ccbar::effective_mass"], p["ee->ccbar::q_0"], eff3770_g0s);

                KMatrix<5, 2> KMatrix32({ee, D0Dbar0, DpDm, eff2S, eff3770}, {psi2S, psi3770}, bkgcst, "ee->ccbar");

                auto tMatrix32s0 = KMatrix32.tmatrix_row(0,  0.0001);
                auto tMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
                auto tMatrix32s2 = KMatrix32.tmatrix_row(0, 20.0);

                TEST_CHECK_RELATIVE_ERROR(real(eff2S->chew_mandelstam(20)),   -0.0172797,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(eff2S->chew_mandelstam(20)),    0.0189092,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(eff3770->chew_mandelstam(20)), -0.0172797,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(eff3770->chew_mandelstam(20)),  0.0189092,  eps);

                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].real(),  3.17782e-5,  eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s0[0].imag(),  1.99853e-11, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].real(), -0.0391188,   eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s1[1].imag(), -0.000556894, eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].real(), -0.0067422,   eps);
                TEST_CHECK_RELATIVE_ERROR(tMatrix32s2[1].imag(),  0.00187893,  eps);
            }

            {
                // Test the full K matrix
                Parameters p = Parameters::Defaults();
                p["mass::psi(3770)"]                      =  3.796443282051135;
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]        = -0.02077753547690923;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]      = -0.001999916489715092;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"]     = -3.3693829070086214;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]   =  8.38428874933062;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"]        =  3.5692138012280807;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]      = -8.391199026268827;
                p["ee->ccbar::g0(psi(2S),eff(2S))"]       =  0.09237874214140328;
                p["ee->ccbar::g0(psi(3770),eff(3770))"]   = -1.5409316160476978;

                p["ee->ccbar::effective_mass"] = 0.1349768;

                // Channels effective momentum
                p["ee->ccbar::q_0"] = 0.5;

                Options oo;
                EEToCCBar c(p, oo);

                auto ir = c.prepare(3.78);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbar0(ir), 3.48882, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDpDm(ir), 2.70723,    eps);

                auto irc = c.prepare_complex(3.78, -0.001);

                // Test the Chew Mandelstam function on the first and second Riemann sheets
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_ee(irc),      -0.112832,    eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_ee(irc),      -0.019891,    eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_II_ee(irc),   -0.112832,    eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_II_ee(irc),    0.0198977,   eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_DpDm(irc),     0.000985784, eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_DpDm(irc),    -0.0006997,   eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_chew_mandelstam_II_DpDm(irc),  0.00102809,  eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_chew_mandelstam_II_DpDm(irc),  0.000664733, eps);

                // Test the amplitude on the first and second Riemann sheets
                TEST_CHECK_RELATIVE_ERROR(c.re_T_eetoDpDm(irc),   0.0275484,  eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_eetoDpDm(irc),  -0.096347,   eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_II_eetoDpDm(irc),  0.0238397,  eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_II_eetoDpDm(irc),  0.111547,  eps);

                // Set ee -> DD cst to .5
                p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"] = 0.5;
                p["ee->ccbar::c(e^+e^-,D^+D^-)"]    = 0.5;

                irc = c.prepare_complex(3.78, -0.001);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_eetoDpDm(irc),   0.288606,   eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_eetoDpDm(irc),  -0.0809925,  eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_II_eetoDpDm(irc),  0.288325,    eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_II_eetoDpDm(irc),  0.0890536,  eps);

                irc = c.prepare_complex(3.8037, -0.03923);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_eetoDpDm(irc),   0.330701,   eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_eetoDpDm(irc),  -0.0672643,  eps);
                TEST_CHECK_RELATIVE_ERROR(c.re_T_II_eetoDpDm(irc),  0.296099,   eps);
                TEST_CHECK_RELATIVE_ERROR(c.im_T_II_eetoDpDm(irc), -0.0682443,  eps);
            }

            {
                // Phase-1/2 extension: the full EEToCCBar now carries 4 resonances and
                // 18 channels (psi(4040) + eff(4040); the equal-mass open-charm channels
                // D_s^+D_s^-, D^*0Dbar^*0, D^*+D^*-; the S- and D-wave D Dbar^* channels;
                // and the second P-wave (5P1) and F-wave (5F1) of D^*Dbar^*). With every
                // new coupling left at its default 0, the new resonance and channels must
                // decouple completely, so the N_C=5 results are reproduced exactly.
                Parameters p = Parameters::Defaults();
                p["mass::psi(3770)"]                      =  3.796443282051135;
                p["ee->ccbar::g0(psi(2S),e^+e^-)"]        = -0.02077753547690923;
                p["ee->ccbar::g0(psi(3770),e^+e^-)"]      = -0.001999916489715092;
                p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"]     = -3.3693829070086214;
                p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]   =  8.38428874933062;
                p["ee->ccbar::g0(psi(2S),D^+D^-)"]        =  3.5692138012280807;
                p["ee->ccbar::g0(psi(3770),D^+D^-)"]      = -8.391199026268827;
                p["ee->ccbar::g0(psi(2S),eff(2S))"]       =  0.09237874214140328;
                p["ee->ccbar::g0(psi(3770),eff(3770))"]   = -1.5409316160476978;

                p["ee->ccbar::effective_mass"] = 0.1349768;
                p["ee->ccbar::q_0"] = 0.5;

                Options oo;
                EEToCCBar c(p, oo);

                // Backward compatibility at the extended (18x4) default size: identical to
                // the N_C=5 value pinned in the block above.
                auto ir = c.prepare(3.78);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbar0(ir), 3.48882, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDpDm(ir),    2.70723, eps);

                // psi(4040) is inert while its couplings are at their default 0.
                TEST_CHECK_NEARLY_EQUAL(c.psi4040_total_width(), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.psi4040_eff_width(),   0.0, eps);

                // The new open-charm channels are inert at default couplings. The equal-mass
                // ones are evaluated above all thresholds (D_sD̄_s 3.937, D*0D̄*0 4.014,
                // D*+D*- 4.021 GeV); the unequal-mass D Dbar^* channels (thresholds 3.872 /
                // 3.880 GeV) likewise.
                auto ir_hi = c.prepare(4.05);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDsDs(ir_hi),            0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0(ir_hi), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstarpDstarm(ir_hi),    0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoD0Dbarstar0(ir_hi),     0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDpDstarm(ir_hi),        0.0, eps);
                // ... including the new D-wave partial waves separately.
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoD0Dbarstar0_D(ir_hi),   0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDpDstarm_D(ir_hi),      0.0, eps);
                // ... and the three D^*Dbar^* partial waves (1P1, 5P1, 5F1) separately.
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0_1P1(ir_hi), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0_5P1(ir_hi), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0_5F1(ir_hi), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstarpDstarm_5F1(ir_hi),    0.0, eps);

                // Positive control: switching the couplings on makes psi(4040) acquire a
                // width and feed the new open-charm channels, confirming they are genuinely
                // wired in (not dead code). The D-wave D Dbar^* coupling is switched on too.
                p["ee->ccbar::g0(psi(4040),e^+e^-)"]            = 5.0;
                p["ee->ccbar::g0(psi(4040),eff(4040))"]         = 2.0;
                p["ee->ccbar::g0(psi(4040),D_s^+D_s^-)"]        = 2.0;
                p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0)"]       = 2.0;
                p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0(5P1))"]  = 2.0;
                p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0(5F1))"]  = 2.0;
                p["ee->ccbar::g0(psi(4040),D^0Dbar^*0)"]        = 2.0;
                p["ee->ccbar::g0(psi(4040),D^0Dbar^*0(D))"]     = 2.0;
                TEST_CHECK(c.psi4040_eff_width()   > 0.0);
                TEST_CHECK(c.psi4040_total_width() > 0.0);

                ir_hi = c.prepare(4.05);
                TEST_CHECK(c.sigma_eetoDsDs(ir_hi)            > 0.0);
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0(ir_hi) > 0.0);
                TEST_CHECK(c.sigma_eetoD0Dbarstar0(ir_hi)     > 0.0);
                // The D-wave partial wave now contributes, and the total is the sum S + D.
                TEST_CHECK(c.sigma_eetoD0Dbarstar0_S(ir_hi)   > 0.0);
                TEST_CHECK(c.sigma_eetoD0Dbarstar0_D(ir_hi)   > 0.0);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbarstar0(ir_hi),
                        c.sigma_eetoD0Dbarstar0_S(ir_hi) + c.sigma_eetoD0Dbarstar0_D(ir_hi), eps);
                // The three D^*Dbar^* partial waves each contribute, and the total is their
                // incoherent sum 1P1 + 5P1 + 5F1.
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_1P1(ir_hi) > 0.0);
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_5P1(ir_hi) > 0.0);
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_5F1(ir_hi) > 0.0);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstar0Dstarbar0(ir_hi),
                          c.sigma_eetoDstar0Dstarbar0_1P1(ir_hi)
                        + c.sigma_eetoDstar0Dstarbar0_5P1(ir_hi)
                        + c.sigma_eetoDstar0Dstarbar0_5F1(ir_hi), eps);
                // The transverse/longitudinal components are non-negative coherent
                // combinations of the same waves, and the orthonormal recoupling makes
                // their sum equal to the total: sigma_TT + sigma_TL + sigma_LL = total.
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_TT(ir_hi) > 0.0);
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_TL(ir_hi) > 0.0);
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_LL(ir_hi) > 0.0);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstar0Dstarbar0(ir_hi),
                          c.sigma_eetoDstar0Dstarbar0_TT(ir_hi)
                        + c.sigma_eetoDstar0Dstarbar0_TL(ir_hi)
                        + c.sigma_eetoDstar0Dstarbar0_LL(ir_hi), eps);

                // Selection-rule check: the mixed component sigma_TL is fed only by the
                // two S=2 waves (5P1, 5F1) -- the spin-singlet 1P1 cannot produce it.
                // With both S=2 couplings off, sigma_TL must vanish and the surviving
                // 1P1 splits as sigma_TT : sigma_LL = 2 : 1 (= 2/3 : 1/3 of the total).
                p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0(5P1))"] = 0.0;
                p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0(5F1))"] = 0.0;
                ir_hi = c.prepare(4.05);
                TEST_CHECK(c.sigma_eetoDstar0Dstarbar0_1P1(ir_hi) > 0.0);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0_5P1(ir_hi), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0_5F1(ir_hi), 0.0, eps);
                TEST_CHECK_NEARLY_EQUAL(c.sigma_eetoDstar0Dstarbar0_TL(ir_hi), 0.0, eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstar0Dstarbar0_TT(ir_hi),
                        2.0 / 3.0 * c.sigma_eetoDstar0Dstarbar0(ir_hi), eps);
                TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoDstar0Dstarbar0_LL(ir_hi),
                        1.0 / 3.0 * c.sigma_eetoDstar0Dstarbar0(ir_hi), eps);
            }
        }
} eetoccbar_test;
