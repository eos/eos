/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_B_LCDAS_PARAM_HH
#define EOS_GUARD_EOS_FORM_FACTORS_B_LCDAS_PARAM_HH 1

#include <eos/form-factors/b-lcdas.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <array>
#include <cmath>

namespace eos
{
    namespace b_lcdas
    {
        namespace aux
        {
            // Implement mathematical expressions as pure functions

            inline double mPow(const double x, const double a) {return std::pow(x, a);}
            inline double mExp(const double x) {return std::exp(x);}

            // Leading-order pieces

            double L0_phi_plus(const double omega0, const std::array<const double, 9> a){
                return 1/omega0 * (a[0] + a[2]/3. + a[4]/5. + a[6]/7. + a[8]/9.);
            }

            double L0inc_phi_plus(const double Omega, const double omega0, const std::array<const double, 9> a) {
                const double xOm = Omega / omega0;
                const double xOm2 = std::pow(xOm, 2);
                const double xOm3 = std::pow(xOm, 3);
                const double xOm4 = std::pow(xOm, 4);
                const double xOm5 = std::pow(xOm, 5);
                const double xOm6 = std::pow(xOm, 6);
                const double xOm7 = std::pow(xOm, 7);
                const double xOm8 = std::pow(xOm, 8);

                const double mexp = std::exp(-xOm);

                return 1/omega0 * (
                    +a[0]*(1. - mexp)
                    +a[1]*(mexp*xOm)
                    +a[2]*(0.3333333333333333 + mexp*(-0.3333333333333333 + 0.6666666666666666*xOm - 0.6666666666666666*xOm2))
                    +a[3]*(mexp*xOm*(1 - xOm + 0.3333333333333333*xOm2))
                    +a[4]*(0.2 + mexp*(-0.2 + 0.8*xOm - 1.6*xOm2 + 0.8*xOm3 - 0.13333333333333333*xOm4))
                    +a[5]*(mexp*xOm*(1 - 2.*xOm + 1.5555555555555556*xOm2 - 0.4444444444444444*xOm3 + 0.044444444444444446*xOm4))
                    +a[6]*(0.14285714285714285 + mexp*(-0.14285714285714285 + 0.8571428571428571*xOm - 2.5714285714285716*xOm2 + 2.4761904761904763*xOm3 - 1.0476190476190477*xOm4 + 0.19047619047619047*xOm5 - 0.012698412698412698*xOm6))
                    +a[7]*(mexp*xOm*(1 - 3.*xOm + 3.6666666666666665*xOm2 - 2.*xOm3 + 0.5333333333333333*xOm4 - 0.06666666666666667*xOm5 + 0.0031746031746031746*xOm6))
                    +a[8]*(0.1111111111111111 + mexp*(-0.1111111111111111 + 0.8888888888888888*xOm - 3.5555555555555554*xOm2 + 5.037037037037037*xOm3 - 3.4074074074074074*xOm4 + 1.1851851851851851*xOm5 - 0.21728395061728395*xOm6 + 0.019753086419753086*xOm7 - 0.0007054673721340388*xOm8))
                );
            }

            double B_phi_plus(const double Omega, const double sigma, const double omega0, const std::array<double, 9> a){
                const double xOm = Omega / omega0;
                const double xOm2 = std::pow(xOm, 2);
                const double xOm3 = std::pow(xOm, 3);
                const double xOm4 = std::pow(xOm, 4);
                const double xOm5 = std::pow(xOm, 5);
                const double xOm6 = std::pow(xOm, 6);
                const double xOm7 = std::pow(xOm, 7);
                const double xOm8 = std::pow(xOm, 8);
                const double xsg = sigma * omega0;
                const double xsg2 = std::pow(xsg, 2);
                const double xsg3 = std::pow(xsg, 3);
                const double xsg4 = std::pow(xsg, 4);
                const double xsg5 = std::pow(xsg, 5);
                const double xsg6 = std::pow(xsg, 6);
                const double xsg7 = std::pow(xsg, 7);
                const double xsg8 = std::pow(xsg, 8);

                return (
                    +a[0]*((1 + (-1. + xOm*(-1. - xsg))*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-2.))
                    +a[1]*((-1. + xsg + (1 + xOm - xsg - xOm*xsg2 + xOm2*(1 + 2.*xsg + xsg2))*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-3.))
                    +a[2]*((1 - 2.*xsg + xsg2 + (-1. + 2.*xsg + 2.*xOm2*xsg - xsg2 + 4.*xOm2*xsg2 + xOm*(-1. + xsg + xsg2 - xsg3) + xOm3*(-0.6666666666666666 - 2.*xsg - 2.*xsg2 - 0.6666666666666666*xsg3) + 2.*xOm2*xsg3)*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-4.))
                    +a[3]*((-1. + 3.*xsg - 3.*xsg2 + xsg3 + (1 - 0.6666666666666666*xOm3 + 0.3333333333333333*xOm4 - 3.*xsg - 4.*xOm3*xsg + 1.3333333333333333*xOm4*xsg + 3.*xsg2 - 8.*xOm3*xsg2 + 2.*xOm4*xsg2 - xsg3 - 6.666666666666667*xOm3*xsg3 + 1.3333333333333333*xOm4*xsg3 + xOm*(1 - 2.*xsg + 2.*xsg3 - xsg4) - 2.*xOm3*xsg4 + 0.3333333333333333*xOm4*xsg4 + xOm2*(1 + 2.*xsg + 4.*xsg2 + 6.*xsg3 + 3.*xsg4))*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-5.))
                    +a[4]*((1 - 4.*xsg + 6.*xsg2 - 4.*xsg3 + xsg4 + (-1. + 0.6666666666666666*xOm4 - 0.13333333333333333*xOm5 + 4.*xsg + 4.*xOm2*xsg + 4.*xOm4*xsg - 0.6666666666666666*xOm5*xsg - 6.*xsg2 + 8.*xOm2*xsg2 + 9.333333333333334*xOm4*xsg2 - 1.3333333333333333*xOm5*xsg2 + 4.*xsg3 + 8.*xOm2*xsg3 + 10.666666666666666*xOm4*xsg3 - 1.3333333333333333*xOm5*xsg3 - xsg4 + 8.*xOm2*xsg4 + 6.*xOm4*xsg4 - 0.6666666666666666*xOm5*xsg4 + xOm3*(-1.3333333333333333 - 6.666666666666667*xsg - 16.*xsg2 - 21.333333333333332*xsg3 - 14.666666666666666*xsg4 - 4.*xsg5) + xOm*(-1. + 3.*xsg - 2.*xsg2 - 2.*xsg3 + 3.*xsg4 - xsg5) + 4.*xOm2*xsg5 + 1.3333333333333333*xOm4*xsg5 - 0.13333333333333333*xOm5*xsg5)*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-6.))
                    +a[5]*((-1. + 5.*xsg - 10.*xsg2 + 10.*xsg3 - 5.*xsg4 + xsg5 + (1 - 1.3333333333333333*xOm3 + 1.3333333333333333*xOm4 - 0.4*xOm5 + 0.044444444444444446*xOm6 - 5.*xsg - 10.666666666666666*xOm3*xsg + 8.666666666666666*xOm4*xsg - 2.6666666666666665*xOm5*xsg + 0.26666666666666666*xOm6*xsg + 10.*xsg2 - 30.666666666666668*xOm3*xsg2 + 24.666666666666668*xOm4*xsg2 - 7.333333333333333*xOm5*xsg2 + 0.6666666666666666*xOm6*xsg2 - 10.*xsg3 - 48.*xOm3*xsg3 + 38.666666666666664*xOm4*xsg3 - 10.666666666666666*xOm5*xsg3 + 0.8888888888888888*xOm6*xsg3 + 5.*xsg4 - 46.666666666666664*xOm3*xsg4 + 34.666666666666664*xOm4*xsg4 - 8.666666666666666*xOm5*xsg4 + 0.6666666666666666*xOm6*xsg4 - xsg5 - 26.666666666666668*xOm3*xsg5 + 16.666666666666668*xOm4*xsg5 - 3.7333333333333334*xOm5*xsg5 + 0.26666666666666666*xOm6*xsg5 + xOm*(1 - 4.*xsg + 5.*xsg2 - 5.*xsg4 + 4.*xsg5 - xsg6) - 6.666666666666667*xOm3*xsg6 + 3.3333333333333335*xOm4*xsg6 - 0.6666666666666666*xOm5*xsg6 + 0.044444444444444446*xOm6*xsg6 + xOm2*(1 + 2.*xsg + 11.*xsg2 + 20.*xsg3 + 15.*xsg4 + 10.*xsg5 + 5.*xsg6))*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-7.))
                    +a[6]*((1 - 6.*xsg + 15.*xsg2 - 20.*xsg3 + 15.*xsg4 - 6.*xsg5 + xsg6 + (-1. + 2.*xOm4 - 0.9333333333333333*xOm5 + 0.17777777777777778*xOm6 - 0.012698412698412698*xOm7 + 6.*xsg + 6.*xOm2*xsg + 16.*xOm4*xsg - 7.066666666666666*xOm5*xsg + 1.3333333333333333*xOm6*xsg - 0.08888888888888889*xOm7*xsg - 15.*xsg2 + 12.*xOm2*xsg2 + 54.*xOm4*xsg2 - 23.333333333333332*xOm5*xsg2 + 4.266666666666667*xOm6*xsg2 - 0.26666666666666666*xOm7*xsg2 + 20.*xsg3 + 26.*xOm2*xsg3 + 102.66666666666667*xOm4*xsg3 - 43.333333333333336*xOm5*xsg3 + 7.555555555555555*xOm6*xsg3 - 0.4444444444444444*xOm7*xsg3 - 15.*xsg4 + 40.*xOm2*xsg4 + 120.66666666666667*xOm4*xsg4 - 48.666666666666664*xOm5*xsg4 + 8.*xOm6*xsg4 - 0.4444444444444444*xOm7*xsg4 + 6.*xsg5 + 26.*xOm2*xsg5 + 88.*xOm4*xsg5 - 32.93333333333333*xOm5*xsg5 + 5.066666666666666*xOm6*xsg5 - 0.26666666666666666*xOm7*xsg5 - xsg6 + 12.*xOm2*xsg6 + 36.666666666666664*xOm4*xsg6 - 12.4*xOm5*xsg6 + 1.7777777777777777*xOm6*xsg6 - 0.08888888888888889*xOm7*xsg6 + xOm3*(-2. - 14.*xsg - 50.*xsg2 - 99.33333333333333*xsg3 - 118.*xsg4 - 90.*xsg5 - 43.333333333333336*xsg6 - 10.*xsg7) + xOm*(-1. + 5.*xsg - 9.*xsg2 + 5.*xsg3 + 5.*xsg4 - 9.*xsg5 + 5.*xsg6 - xsg7) + 6.*xOm2*xsg7 + 6.666666666666667*xOm4*xsg7 - 2.*xOm5*xsg7 + 0.26666666666666666*xOm6*xsg7 - 0.012698412698412698*xOm7*xsg7)*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-8.))
                    +a[7]*((-1. + 7.*xsg - 21.*xsg2 + 35.*xsg3 - 35.*xsg4 + 21.*xsg5 - 7.*xsg6 + xsg7 + (1 - 2.*xOm3 + 3.*xOm4 - 1.7333333333333334*xOm5 + 0.4888888888888889*xOm6 - 0.06349206349206349*xOm7 + 0.0031746031746031746*xOm8 - 7.*xsg - 20.*xOm3*xsg + 26.*xOm4*xsg - 15.2*xOm5*xsg + 4.177777777777778*xOm6*xsg - 0.5333333333333333*xOm7*xsg + 0.025396825396825397*xOm8*xsg + 21.*xsg2 - 76.*xOm3*xsg2 + 102.*xOm4*xsg2 - 58.4*xOm5*xsg2 + 15.733333333333333*xOm6*xsg2 - 1.9555555555555555*xOm7*xsg2 + 0.08888888888888889*xOm8*xsg2 - 35.*xsg3 - 174.66666666666666*xOm3*xsg3 + 231.33333333333334*xOm4*xsg3 - 129.33333333333334*xOm5*xsg3 + 34.044444444444444*xOm6*xsg3 - 4.088888888888889*xOm7*xsg3 + 0.17777777777777778*xOm8*xsg3 + 35.*xsg4 - 261.3333333333333*xOm3*xsg4 + 332.*xOm4*xsg4 - 181.33333333333334*xOm5*xsg4 + 46.22222222222222*xOm6*xsg4 - 5.333333333333333*xOm7*xsg4 + 0.2222222222222222*xOm8*xsg4 - 21.*xsg5 - 252.*xOm3*xsg5 + 312.6666666666667*xOm4*xsg5 - 165.06666666666666*xOm5*xsg5 + 40.266666666666666*xOm6*xsg5 - 4.444444444444445*xOm7*xsg5 + 0.17777777777777778*xOm8*xsg5 + 7.*xsg6 - 158.66666666666666*xOm3*xsg6 + 191.33333333333334*xOm4*xsg6 - 95.2*xOm5*xsg6 + 21.955555555555556*xOm6*xsg6 - 2.311111111111111*xOm7*xsg6 + 0.08888888888888889*xOm8*xsg6 - xsg7 - 65.33333333333333*xOm3*xsg7 + 70.*xOm4*xsg7 - 31.733333333333334*xOm5*xsg7 + 6.844444444444444*xOm6*xsg7 - 0.6857142857142857*xOm7*xsg7 + 0.025396825396825397*xOm8*xsg7 + xOm*(1 - 6.*xsg + 14.*xsg2 - 14.*xsg3 + 14.*xsg5 - 14.*xsg6 + 6.*xsg7 - xsg8) - 14.*xOm3*xsg8 + 11.666666666666666*xOm4*xsg8 - 4.666666666666667*xOm5*xsg8 + 0.9333333333333333*xOm6*xsg8 - 0.08888888888888889*xOm7*xsg8 + 0.0031746031746031746*xOm8*xsg8 + xOm2*(1 + 2.*xsg + 22.*xsg2 + 42.*xsg3 + 56.*xsg4 + 70.*xsg5 + 42.*xsg6 + 14.*xsg7 + 7.*xsg8))*mExp(xOm*(-1. - xsg)))*mPow(1. + xsg,-9.))
                    +a[8]*((1 - 8.*xsg + 28.*xsg2 - 56.*xsg3 + 70.*xsg4 - 56.*xsg5 + 28.*xsg6 - 8.*xsg7 + xsg8 + mExp(xOm*(-1. - xsg))*(-1. - xOm - 2.6666666666666665*xOm3 + 4.*xOm4 - 2.933333333333333*xOm5 + 1.0666666666666667*xOm6 - 0.20317460317460317*xOm7 + 0.01904761904761905*xOm8 + 8.*xsg + 7.*xOm*xsg + 8.*xOm2*xsg - 24.*xOm3*xsg + 40.*xOm4*xsg - 28.533333333333335*xOm5*xsg + 10.311111111111112*xOm6*xsg - 1.9301587301587302*xOm7*xsg + 0.17777777777777778*xOm8*xsg - 28.*xsg2 - 20.*xOm*xsg2 + 16.*xOm2*xsg2 - 112.*xOm3*xsg2 + 176.*xOm4*xsg2 - 124.8*xOm5*xsg2 + 44.44444444444444*xOm6*xsg2 - 8.177777777777777*xOm7*xsg2 + 0.7365079365079366*xOm8*xsg2 + 56.*xsg3 + 28.*xOm*xsg3 + 64.*xOm2*xsg3 - 293.3333333333333*xOm3*xsg3 + 458.6666666666667*xOm4*xsg3 - 321.06666666666666*xOm5*xsg3 + 112.35555555555555*xOm6*xsg3 - 20.266666666666666*xOm7*xsg3 + 1.7777777777777777*xOm8*xsg3 - 70.*xsg4 - 14.*xOm*xsg4 + 112.*xOm2*xsg4 - 501.3333333333333*xOm3*xsg4 + 781.3333333333334*xOm4*xsg4 - 536.*xOm5*xsg4 + 183.82222222222222*xOm6*xsg4 - 32.355555555555554*xOm7*xsg4 + 2.7555555555555555*xOm8*xsg4 + 56.*xsg5 - 14.*xOm*xsg5 + 112.*xOm2*xsg5 - 597.3333333333334*xOm3*xsg5 + 901.3333333333334*xOm4*xsg5 - 604.2666666666667*xOm5*xsg5 + 201.95555555555555*xOm6*xsg5 - 34.48888888888889*xOm7*xsg5 + 2.8444444444444446*xOm8*xsg5 - 28.*xsg6 + 28.*xOm*xsg6 + 112.*xOm2*xsg6 - 485.3333333333333*xOm3*xsg6 + 709.3333333333334*xOm4*xsg6 - 461.8666666666667*xOm5*xsg6 + 148.9777777777778*xOm6*xsg6 - 24.533333333333335*xOm7*xsg6 + 1.9555555555555555*xOm8*xsg6 + 8.*xsg7 - 20.*xOm*xsg7 + 64.*xOm2*xsg7 - 261.3333333333333*xOm3*xsg7 + 373.3333333333333*xOm4*xsg7 - 231.46666666666667*xOm5*xsg7 + 71.11111111111111*xOm6*xsg7 - 11.225396825396825*xOm7*xsg7 + 0.8634920634920635*xOm8*xsg7 - xsg8 + 7.*xOm*xsg8 + 16.*xOm2*xsg8 - 93.33333333333333*xOm3*xsg8 + 121.33333333333333*xOm4*xsg8 - 69.06666666666666*xOm5*xsg8 + 19.91111111111111*xOm6*xsg8 - 2.996825396825397*xOm7*xsg8 + 0.2222222222222222*xOm8*xsg8 + mPow(xOm,9.)*(-0.0007054673721340388 - 0.006349206349206349*xsg - 0.025396825396825397*xsg2 - 0.05925925925925926*xsg3 - 0.08888888888888889*xsg4 - 0.08888888888888889*xsg5 - 0.05925925925925926*xsg6 - 0.025396825396825397*xsg7 - 0.006349206349206349*xsg8 - 0.0007054673721340388*mPow(xsg,9.)) + (-1.*xOm + 8.*xOm2 - 18.666666666666668*xOm3 + 18.666666666666668*xOm4 - 9.333333333333334*xOm5 + 2.488888888888889*xOm6 - 0.35555555555555557*xOm7 + 0.025396825396825397*xOm8)*mPow(xsg,9.)))*mPow(1. + xsg,-10.))
                );
            }

            // Next-to-leading-order pieces

            double L0_Dphi_plus_eff_1(const double Egamma, const double mu, const double omega0, const std::array<const double, 9> a) {
                // Sum of all contributions 1a, 1b, 1c and 1d
                const double mlog = -0.11593151565841242 + std::log(std::pow(mu,2)/(Egamma*omega0));
                const double mlog2 = std::pow(mlog, 2);

                return 1/omega0 * (
                    +a[0]*(-1. + mlog2)
                    +a[1]*(2.*mlog)
                    +a[2]*(1. + 0.3333333333333333*mlog2)
                    +a[3]*(1.3333333333333333*mlog)
                    +a[4]*(1.1333333333333333 + 0.2*mlog2)
                    +a[5]*(1.0222222222222221*mlog)
                    +a[6]*(1.1015873015873017 + 0.14285714285714285*mlog2)
                    +a[7]*(0.8380952380952381*mlog)
                    +a[8]*(1.0430335097001764 + 0.1111111111111111*mlog2)
                );
            }
        }
    }
}

namespace eos
{
    namespace b_lcdas
    {
        /*!
        * Parametrization of the B-meson LCDAs
        */
        class Param:
            public BMesonLCDAs,
            public PrivateImplementationPattern<Param>
        {
            public:
                Param(const Parameters &, const Options &);
                ~Param();

                /*!
                *  B to gamma l nu
                */
                double L0() const;
                double L0inc(const double & Omega) const;
                double Binc(const double & Omega, const double & sigma) const;

                /*!
                * Leading twist two-particle LCDAs
                * omega: plus-component of the spectator momentum
                */
                double phi_plus(const double & omega) const;
                double phi_minus(const double & omega) const;
                double phi_bar(const double & omega) const;
                double phi_bar_d1(const double & omega) const;

                /*!
                * Next-to-leading twist two-particle LCDAs
                * omega: plus-component of the spectator momentum
                */
                double g_plus(const double & omega) const;
                double g_plus_d1(const double & omega) const;
                double g_plus_d2(const double & omega) const;

                double g_minusWW(const double & omega) const;
                double g_minusWW_d1(const double & omega) const;
                double g_minusWW_d2(const double & omega) const;

                double g_bar(const double & omega) const;
                double g_bar_d1(const double & omega) const;
                double g_bar_d2(const double & omega) const;
                double g_bar_d3(const double & omega) const;

                /*!
                * Leading power three-particle LCDAs
                * omega_1: plus-component of the spectator momentum
                * omega_2: plus-component of the gluon momentum
                * */
                double phi_3(const double & omega_1, const double & omega_2) const;
                double phi_4(const double & omega_1, const double & omega_2) const;

                double phi_bar_3(const double & omega_1, const double & omega_2) const;
                double phi_bar_4(const double & omega_1, const double & omega_2) const;

                double phi_bar2_3(const double & omega_1, const double & omega_2) const;
                double phi_bar2_4(const double & omega_1, const double & omega_2) const;

                double phi_bar_bar_3(const double & omega_1, const double & omega_2) const;
                double phi_bar_bar_4(const double & omega_1, const double & omega_2) const;

                double psi_bar_4(const double & omega_1, const double & omega_2) const;
                double chi_bar_4(const double & omega_1, const double & omega_2) const;

                double psi_bar_bar_4(const double & omega_1, const double & omega_2) const;
                double chi_bar_bar_4(const double & omega_1, const double & omega_2) const;
                /*!
                * Pseudo observables for the two-particle LCDAs
                */
                double inverse_lambda_plus() const;

                /*!
                * Leading power three-particle LCDAs
                * omega: plus-component of the spectator momentum
                * xi:    plus-component of the gluon momentum
                * */
                double psi_A(const double & omega, const double & xi) const;
                double psi_V(const double & omega, const double & xi) const;
                double X_A(const double & omega, const double & xi) const;
                double Y_A(const double & omega, const double & xi) const;

                /*!
                * Auxiliary functions for the three-particle LCDAs
                * See [KMO2006], below eq. (72), p. 28 for their definition.
                */
                double Xbar_A(const double & omega, const double & xi) const;
                double Ybar_A(const double & omega, const double & xi) const;

                /* Internal diagnostics */
                Diagnostics diagnostics() const;
        };
    }
}

#endif
