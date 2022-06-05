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
                    +a[0]*(1 - mexp)
                    +a[1]*(mexp*xOm)
                    +a[2]*(0.3333333333333333 + (mexp*(-1 + 2*xOm - 2*xOm2))/3.)
                    +a[3]*((mexp*xOm*(3 - 3*xOm + xOm2))/3.)
                    +a[4]*(0.2 + (mexp*(-3 + 12*xOm - 24*xOm2 + 12*xOm3 - 2*xOm4))/15.)
                    +a[5]*((mexp*xOm*(45 - 90*xOm + 70*xOm2 - 20*xOm3 + 2*xOm4))/45.)
                    +a[6]*(0.14285714285714285 + (mexp*(-45 + 270*xOm - 810*xOm2 + 780*xOm3 - 330*xOm4 + 60*xOm5 - 4*xOm6))/315.)
                    +a[7]*((mexp*xOm*(315 - 945*xOm + 1155*xOm2 - 630*xOm3 + 168*xOm4 - 21*xOm5 + xOm6))/315.)
                    +a[8]*(0.1111111111111111 + (mexp*(-315 + 2520*xOm - 10080*xOm2 + 14280*xOm3 - 9660*xOm4 + 3360*xOm5 - 616*xOm6 + 56*xOm7 - 2*xOm8))/2835.)
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
                    +a[0]*((1 - (1 + xOm + xOm*xsg)*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-2))
                    +a[1]*((-1 + xsg + (1 + xOm - xsg - xOm*xsg2 + xOm2*(1 + 2*xsg + xsg2))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-3))
                    +a[2]*(((3 - 6*xsg + 3*xsg2 + (-3 + 6*xsg + 6*xOm2*xsg - 3*xsg2 + 12*xOm2*xsg2 + 3*xOm*(-1 + xsg + xsg2 - xsg3) + 6*xOm2*xsg3 - 2*xOm3*(1 + 3*xsg + 3*xsg2 + xsg3))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-4))/3.)
                    +a[3]*(((-3 + 9*xsg - 9*xsg2 + 3*xsg3 + (3 - 2*xOm3 + xOm4 - 9*xsg - 12*xOm3*xsg + 4*xOm4*xsg + 9*xsg2 - 24*xOm3*xsg2 + 6*xOm4*xsg2 - 3*xsg3 - 20*xOm3*xsg3 + 4*xOm4*xsg3 - 6*xOm3*xsg4 + xOm4*xsg4 - 3*xOm*(-1 + 2*xsg - 2*xsg3 + xsg4) + 3*xOm2*(1 + 2*xsg + 4*xsg2 + 6*xsg3 + 3*xsg4))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-5))/3.)
                    +a[4]*(((15 - 60*xsg + 90*xsg2 - 60*xsg3 + 15*xsg4 + (-15 + 10*xOm4 - 2*xOm5 + 60*xsg + 60*xOm2*xsg + 60*xOm4*xsg - 10*xOm5*xsg - 90*xsg2 + 120*xOm2*xsg2 + 140*xOm4*xsg2 - 20*xOm5*xsg2 + 60*xsg3 + 120*xOm2*xsg3 + 160*xOm4*xsg3 - 20*xOm5*xsg3 - 15*xsg4 + 120*xOm2*xsg4 + 90*xOm4*xsg4 - 10*xOm5*xsg4 + 15*xOm*(-1 + 3*xsg - 2*xsg2 - 2*xsg3 + 3*xsg4 - xsg5) + 60*xOm2*xsg5 + 20*xOm4*xsg5 - 2*xOm5*xsg5 - 20*xOm3*(1 + 5*xsg + 12*xsg2 + 16*xsg3 + 11*xsg4 + 3*xsg5))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-6))/15.)
                    +a[5]*(((-45 + 225*xsg - 450*xsg2 + 450*xsg3 - 225*xsg4 + 45*xsg5 + (45 - 60*xOm3 + 60*xOm4 - 18*xOm5 + 2*xOm6 - 225*xsg - 480*xOm3*xsg + 390*xOm4*xsg - 120*xOm5*xsg + 12*xOm6*xsg + 450*xsg2 - 1380*xOm3*xsg2 + 1110*xOm4*xsg2 - 330*xOm5*xsg2 + 30*xOm6*xsg2 - 450*xsg3 - 2160*xOm3*xsg3 + 1740*xOm4*xsg3 - 480*xOm5*xsg3 + 40*xOm6*xsg3 + 225*xsg4 - 2100*xOm3*xsg4 + 1560*xOm4*xsg4 - 390*xOm5*xsg4 + 30*xOm6*xsg4 - 45*xsg5 - 1200*xOm3*xsg5 + 750*xOm4*xsg5 - 168*xOm5*xsg5 + 12*xOm6*xsg5 - 300*xOm3*xsg6 + 150*xOm4*xsg6 - 30*xOm5*xsg6 + 2*xOm6*xsg6 - 45*xOm*(-1 + 4*xsg - 5*xsg2 + 5*xsg4 - 4*xsg5 + xsg6) + 45*xOm2*(1 + 2*xsg + 11*xsg2 + 20*xsg3 + 15*xsg4 + 10*xsg5 + 5*xsg6))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-7))/45.)
                    +a[6]*(((-315*(-1 + 6*xsg - 15*xsg2 + 20*xsg3 - 15*xsg4 + 6*xsg5 - xsg6) + (-315 + 630*xOm4 - 294*xOm5 + 56*xOm6 - 4*xOm7 + 1890*xsg + 1890*xOm2*xsg + 5040*xOm4*xsg - 2226*xOm5*xsg + 420*xOm6*xsg - 28*xOm7*xsg - 4725*xsg2 + 3780*xOm2*xsg2 + 17010*xOm4*xsg2 - 7350*xOm5*xsg2 + 1344*xOm6*xsg2 - 84*xOm7*xsg2 + 6300*xsg3 + 8190*xOm2*xsg3 + 32340*xOm4*xsg3 - 13650*xOm5*xsg3 + 2380*xOm6*xsg3 - 140*xOm7*xsg3 - 4725*xsg4 + 12600*xOm2*xsg4 + 38010*xOm4*xsg4 - 15330*xOm5*xsg4 + 2520*xOm6*xsg4 - 140*xOm7*xsg4 + 1890*xsg5 + 8190*xOm2*xsg5 + 27720*xOm4*xsg5 - 10374*xOm5*xsg5 + 1596*xOm6*xsg5 - 84*xOm7*xsg5 - 315*xsg6 + 3780*xOm2*xsg6 + 11550*xOm4*xsg6 - 3906*xOm5*xsg6 + 560*xOm6*xsg6 - 28*xOm7*xsg6 + 315*xOm*(-1 + 5*xsg - 9*xsg2 + 5*xsg3 + 5*xsg4 - 9*xsg5 + 5*xsg6 - xsg7) + 1890*xOm2*xsg7 + 2100*xOm4*xsg7 - 630*xOm5*xsg7 + 84*xOm6*xsg7 - 4*xOm7*xsg7 - 210*xOm3*(3 + 21*xsg + 75*xsg2 + 149*xsg3 + 177*xsg4 + 135*xsg5 + 65*xsg6 + 15*xsg7))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-8))/315.)
                    +a[7]*(((-315 + 2205*xsg - 6615*xsg2 + 11025*xsg3 - 11025*xsg4 + 6615*xsg5 - 2205*xsg6 + 315*xsg7 + (315 - 630*xOm3 + 945*xOm4 - 546*xOm5 + 154*xOm6 - 20*xOm7 + xOm8 - 2205*xsg - 6300*xOm3*xsg + 8190*xOm4*xsg - 4788*xOm5*xsg + 1316*xOm6*xsg - 168*xOm7*xsg + 8*xOm8*xsg + 6615*xsg2 - 23940*xOm3*xsg2 + 32130*xOm4*xsg2 - 18396*xOm5*xsg2 + 4956*xOm6*xsg2 - 616*xOm7*xsg2 + 28*xOm8*xsg2 - 11025*xsg3 - 55020*xOm3*xsg3 + 72870*xOm4*xsg3 - 40740*xOm5*xsg3 + 10724*xOm6*xsg3 - 1288*xOm7*xsg3 + 56*xOm8*xsg3 + 11025*xsg4 - 82320*xOm3*xsg4 + 104580*xOm4*xsg4 - 57120*xOm5*xsg4 + 14560*xOm6*xsg4 - 1680*xOm7*xsg4 + 70*xOm8*xsg4 - 6615*xsg5 - 79380*xOm3*xsg5 + 98490*xOm4*xsg5 - 51996*xOm5*xsg5 + 12684*xOm6*xsg5 - 1400*xOm7*xsg5 + 56*xOm8*xsg5 + 2205*xsg6 - 49980*xOm3*xsg6 + 60270*xOm4*xsg6 - 29988*xOm5*xsg6 + 6916*xOm6*xsg6 - 728*xOm7*xsg6 + 28*xOm8*xsg6 - 315*xsg7 - 20580*xOm3*xsg7 + 22050*xOm4*xsg7 - 9996*xOm5*xsg7 + 2156*xOm6*xsg7 - 216*xOm7*xsg7 + 8*xOm8*xsg7 - 4410*xOm3*xsg8 + 3675*xOm4*xsg8 - 1470*xOm5*xsg8 + 294*xOm6*xsg8 - 28*xOm7*xsg8 + xOm8*xsg8 - 315*xOm*(-1 + 6*xsg - 14*xsg2 + 14*xsg3 - 14*xsg5 + 14*xsg6 - 6*xsg7 + xsg8) + 315*xOm2*(1 + 2*xsg + 22*xsg2 + 42*xsg3 + 56*xsg4 + 70*xsg5 + 42*xsg6 + 14*xsg7 + 7*xsg8))*mExp(-(xOm*(1 + xsg))))*mPow(1 + xsg,-9))/315.)
                    +a[8]*(((-2835*(-1 + 8*xsg - 28*xsg2 + 56*xsg3 - 70*xsg4 + 56*xsg5 - 28*xsg6 + 8*xsg7 - xsg8) + mExp(-(xOm*(1 + xsg)))*(-2*mPow(xOm,9)*(1 + 9*xsg + 36*xsg2 + 84*xsg3 + 126*xsg4 + 126*xsg5 + 84*xsg6 + 36*xsg7 + 9*xsg8 + mPow(xsg,9)) + 9*(-315 - 840*xOm3 + 1260*xOm4 - 924*xOm5 + 336*xOm6 - 64*xOm7 + 6*xOm8 + 2520*xsg + 2520*xOm2*xsg - 7560*xOm3*xsg + 12600*xOm4*xsg - 8988*xOm5*xsg + 3248*xOm6*xsg - 608*xOm7*xsg + 56*xOm8*xsg - 8820*xsg2 + 5040*xOm2*xsg2 - 35280*xOm3*xsg2 + 55440*xOm4*xsg2 - 39312*xOm5*xsg2 + 14000*xOm6*xsg2 - 2576*xOm7*xsg2 + 232*xOm8*xsg2 + 17640*xsg3 + 20160*xOm2*xsg3 - 92400*xOm3*xsg3 + 144480*xOm4*xsg3 - 101136*xOm5*xsg3 + 35392*xOm6*xsg3 - 6384*xOm7*xsg3 + 560*xOm8*xsg3 - 22050*xsg4 + 35280*xOm2*xsg4 - 157920*xOm3*xsg4 + 246120*xOm4*xsg4 - 168840*xOm5*xsg4 + 57904*xOm6*xsg4 - 10192*xOm7*xsg4 + 868*xOm8*xsg4 + 17640*xsg5 + 35280*xOm2*xsg5 - 188160*xOm3*xsg5 + 283920*xOm4*xsg5 - 190344*xOm5*xsg5 + 63616*xOm6*xsg5 - 10864*xOm7*xsg5 + 896*xOm8*xsg5 - 8820*xsg6 + 35280*xOm2*xsg6 - 152880*xOm3*xsg6 + 223440*xOm4*xsg6 - 145488*xOm5*xsg6 + 46928*xOm6*xsg6 - 7728*xOm7*xsg6 + 616*xOm8*xsg6 + 2520*xsg7 + 20160*xOm2*xsg7 - 82320*xOm3*xsg7 + 117600*xOm4*xsg7 - 72912*xOm5*xsg7 + 22400*xOm6*xsg7 - 3536*xOm7*xsg7 + 272*xOm8*xsg7 - 315*xsg8 + 5040*xOm2*xsg8 - 29400*xOm3*xsg8 + 38220*xOm4*xsg8 - 21756*xOm5*xsg8 + 6272*xOm6*xsg8 - 944*xOm7*xsg8 + 70*xOm8*xsg8 + 315*xOm*(-1 + 7*xsg - 20*xsg2 + 28*xsg3 - 14*xsg4 - 14*xsg5 + 28*xsg6 - 20*xsg7 + 7*xsg8) + (-315*xOm + 4*(630*xOm2 - 1470*xOm3 + 1470*xOm4 - 735*xOm5 + 196*xOm6 - 28*xOm7 + 2*xOm8))*mPow(xsg,9))))*mPow(1 + xsg,-10))/2835.)
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
