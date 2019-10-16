/* vim: set sw=4 sts=4 tw=140 et foldmethod=marker : */

/*
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2019 Nico Gubernari
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

#include <eos/form-factors/unitarity-bounds.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>

#include <iostream>

namespace eos
{
    template <> struct Implementation<BGLCoefficients>
    {
        //B^(*)->D^(*)
        // parameters for the leading Isgur-Wise function xi
        UsedParameter xipone, xippone, xipppone;

        // parameters for the subleading Isgur-Wise function chi2
        UsedParameter chi2one, chi2pone, chi2ppone;

        // parameters for the subleading Isgur-Wise function chi3
        UsedParameter chi3pone, chi3ppone;

        // parameters for the subleading Isgur-Wise function eta
        UsedParameter etaone, etapone, etappone;

        // parameters for subsubleading 1/mc corrections in h+ (B->D), equal to delta{h+}
        UsedParameter l1one, l1pone, l1ppone;

        // parameters for subsubleading 1/mc corrections in hA1 (B->D^*), equal to delta{A1}
        UsedParameter l2one, l2pone, l2ppone;

        // parameters for subsubleading 1/m_c corrections
        UsedParameter l3one, l3pone, l3ppone;
        UsedParameter l4one, l4pone, l4ppone;
        UsedParameter l5one, l5pone, l5ppone;
        UsedParameter l6one, l6pone, l6ppone;

        //B^_s(*)->D_s^(*)
        // parameters for the leading Isgur-Wise function xi
        UsedParameter xispone, xisppone, xispppone;

        // parameters for the subleading Isgur-Wise function chi2
        UsedParameter chi2sone, chi2spone, chi2sppone;

        // parameters for the subleading Isgur-Wise function chi3
        UsedParameter chi3spone, chi3sppone;

        // parameters for the subleading Isgur-Wise function eta
        UsedParameter etasone, etaspone, etasppone;

        // parameters for subsubleading 1/mc corrections in h+ (B->D), equal to delta{h+}
        UsedParameter l1sone, l1spone, l1sppone;

        // parameters for subsubleading 1/mc corrections in hA1 (B->D^*), equal to delta{A1}
        UsedParameter l2sone, l2spone, l2sppone;

        // parameters for subsubleading 1/m_c corrections
        UsedParameter l3sone, l3spone, l3sppone;
        UsedParameter l4sone, l4spone, l4sppone;
        UsedParameter l5sone, l5spone, l5sppone;
        UsedParameter l6sone, l6spone, l6sppone;

        Implementation(const Parameters & p, const Options & /*o*/, ParameterUser & u) :
            xipone(p["B(*)->D(*)::xi'(1)@HQET"], u),
            xippone(p["B(*)->D(*)::xi''(1)@HQET"], u),
            xipppone(p["B(*)->D(*)::xi'''(1)@HQET"], u),
            chi2one(p["B(*)->D(*)::chi_2(1)@HQET"], u),
            chi2pone(p["B(*)->D(*)::chi_2'(1)@HQET"], u),
            chi2ppone(p["B(*)->D(*)::chi_2''(1)@HQET"], u),
            chi3pone(p["B(*)->D(*)::chi_3'(1)@HQET"], u),
            chi3ppone(p["B(*)->D(*)::chi_3''(1)@HQET"], u),
            etaone(p["B(*)->D(*)::eta(1)@HQET"], u),
            etapone(p["B(*)->D(*)::eta'(1)@HQET"], u),
            etappone(p["B(*)->D(*)::eta''(1)@HQET"], u),
            l1one(p["B(*)->D(*)::l_1(1)@HQET"], u),
            l1pone(p["B(*)->D(*)::l_1'(1)@HQET"], u),
            l1ppone(p["B(*)->D(*)::l_1''(1)@HQET"], u),
            l2one(p["B(*)->D(*)::l_2(1)@HQET"], u),
            l2pone(p["B(*)->D(*)::l_2'(1)@HQET"], u),
            l2ppone(p["B(*)->D(*)::l_2''(1)@HQET"], u),
            l3one(p["B(*)->D(*)::l_3(1)@HQET"], u),
            l3pone(p["B(*)->D(*)::l_3'(1)@HQET"], u),
            l3ppone(p["B(*)->D(*)::l_3''(1)@HQET"], u),
            l4one(p["B(*)->D(*)::l_4(1)@HQET"], u),
            l4pone(p["B(*)->D(*)::l_4'(1)@HQET"], u),
            l4ppone(p["B(*)->D(*)::l_4''(1)@HQET"], u),
            l5one(p["B(*)->D(*)::l_5(1)@HQET"], u),
            l5pone(p["B(*)->D(*)::l_5'(1)@HQET"], u),
            l5ppone(p["B(*)->D(*)::l_5''(1)@HQET"], u),
            l6one(p["B(*)->D(*)::l_6(1)@HQET"], u),
            l6pone(p["B(*)->D(*)::l_6'(1)@HQET"], u),
            l6ppone(p["B(*)->D(*)::l_6''(1)@HQET"], u),
            xispone(p["B_s(*)->D_s(*)::xi'(1)@HQET"], u),
            xisppone(p["B_s(*)->D_s(*)::xi''(1)@HQET"], u),
            xispppone(p["B_s(*)->D_s(*)::xi'''(1)@HQET"], u),
            chi2sone(p["B_s(*)->D_s(*)::chi_2(1)@HQET"], u),
            chi2spone(p["B_s(*)->D_s(*)::chi_2'(1)@HQET"], u),
            chi2sppone(p["B_s(*)->D_s(*)::chi_2''(1)@HQET"], u),
            chi3spone(p["B_s(*)->D_s(*)::chi_3'(1)@HQET"], u),
            chi3sppone(p["B_s(*)->D_s(*)::chi_3''(1)@HQET"], u),
            etasone(p["B_s(*)->D_s(*)::eta(1)@HQET"], u),
            etaspone(p["B_s(*)->D_s(*)::eta'(1)@HQET"], u),
            etasppone(p["B_s(*)->D_s(*)::eta''(1)@HQET"], u),
            l1sone(p["B_s(*)->D_s(*)::l_1(1)@HQET"], u),
            l1spone(p["B_s(*)->D_s(*)::l_1'(1)@HQET"], u),
            l1sppone(p["B_s(*)->D_s(*)::l_1''(1)@HQET"], u),
            l2sone(p["B_s(*)->D_s(*)::l_2(1)@HQET"], u),
            l2spone(p["B_s(*)->D_s(*)::l_2'(1)@HQET"], u),
            l2sppone(p["B_s(*)->D_s(*)::l_2''(1)@HQET"], u),
            l3sone(p["B_s(*)->D_s(*)::l_3(1)@HQET"], u),
            l3spone(p["B_s(*)->D_s(*)::l_3'(1)@HQET"], u),
            l3sppone(p["B_s(*)->D_s(*)::l_3''(1)@HQET"], u),
            l4sone(p["B_s(*)->D_s(*)::l_4(1)@HQET"], u),
            l4spone(p["B_s(*)->D_s(*)::l_4'(1)@HQET"], u),
            l4sppone(p["B_s(*)->D_s(*)::l_4''(1)@HQET"], u),
            l5sone(p["B_s(*)->D_s(*)::l_5(1)@HQET"], u),
            l5spone(p["B_s(*)->D_s(*)::l_5'(1)@HQET"], u),
            l5sppone(p["B_s(*)->D_s(*)::l_5''(1)@HQET"], u),
            l6sone(p["B_s(*)->D_s(*)::l_6(1)@HQET"], u),
            l6spone(p["B_s(*)->D_s(*)::l_6'(1)@HQET"], u),
            l6sppone(p["B_s(*)->D_s(*)::l_6''(1)@HQET"], u)
        {
        }

        ~Implementation() = default;

        /*
         * HQET parameters following [BLPR2017]
         */
        inline double _mu() const { return 2.31; } // mu^2 = m_b * m_c
        inline double _alpha_s() const { return 0.26; }
        inline double _m_b_1S() const { return 4.71; }
        inline double _m_b_pole() const { return _m_b_1S() * (1 + 2.0 / 9.0 * power_of<2>(_alpha_s())); }
        inline double _m_c_pole() const { return _m_b_pole() - 3.40; }
        inline double _lambda_1() const { return -0.30; }
        // q = u, d
        inline double _LambdaBar() const { return 5.313 - _m_b_pole() + _lambda_1() / (2.0 * _m_b_1S()); }
        inline double _eps_b() const { return _LambdaBar() / (2.0 * _m_b_pole()); }
        inline double _eps_c() const { return _LambdaBar() / (2.0 * _m_c_pole()); }
        // q = s
        inline double _Lambda_sBar() const { return 5.403 - _m_b_pole() + _lambda_1() / (2.0 * _m_b_1S()); }
        inline double _eps_b_s() const { return _Lambda_sBar() / (2.0 * _m_b_pole()); }
        inline double _eps_c_s() const { return _Lambda_sBar() / (2.0 * _m_c_pole()); }

        // B -> D form factors
        // {{{
        double V1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.008703312831206976 + 0.003212569766727432*as + -0.004155423451479522*epsb + 0.004155423451479522*epsc + 0.008310846902959045*epsb*etaone + -0.008310846902959045*epsc*etaone + 0.008703312831206976*l1one*pow(epsc,2) + -0.004155423451479522*l4one*pow(epsc,2);
        }

        double V1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.060227992698537264 + 0.029238681829861973*as + -0.02875604015951836*epsb + -0.27850601059862323*chi2one*epsb + 0.8355180317958697*chi3pone*epsb + 0.02875604015951836*epsc + -0.27850601059862323*chi2one*epsc + 0.8355180317958697*chi3pone*epsc + 0.05751208031903672*epsb*etaone + -0.05751208031903672*epsc*etaone + 0.06648677522367237*epsb*etapone + -0.06648677522367237*epsc*etapone + 0.06962650264965581*xipone + 0.025700558133819457*as*xipone + -0.033243387611836185*epsb*xipone + 0.033243387611836185*epsc*xipone + 0.06648677522367237*epsb*etaone*xipone + -0.06648677522367237*epsc*etaone*xipone + 0.060227992698537264*l1one*pow(epsc,2) + 0.06962650264965581*l1pone*pow(epsc,2) + -0.02875604015951836*l4one*pow(epsc,2) + -0.033243387611836185*l4pone*pow(epsc,2) + 0.06962650264965581*l1one*xipone*pow(epsc,2) + -0.033243387611836185*l4one*xipone*pow(epsc,2);
        }

        double V1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();


            return 0.1306576264897831 + 0.015976186987051347*as + -0.062382885202461706*epsb + -2.4843077875504394*chi2one*epsb + -2.2280480847889863*chi2pone*epsb + 7.452923362651318*chi3pone*epsb + 3.3420721271834792*chi3ppone*epsb + 0.062382885202461706*epsc + -2.4843077875504394*chi2one*epsc + -2.2280480847889863*chi2pone*epsc + 7.452923362651318*chi3pone*epsc + 3.3420721271834792*chi3ppone*epsc + 0.12476577040492341*epsb*etaone + -0.12476577040492341*epsc*etaone + 0.5930701929996385*epsb*etapone + -0.5930701929996385*epsc*etapone + 0.2659471008946895*epsb*etappone + -0.2659471008946895*epsc*etappone + 0.6210769468876098*xipone + 0.2853105709065347*as*xipone + -0.29653509649981924*epsb*xipone + -2.2280480847889863*chi2one*epsb*xipone + 6.6841442543669585*chi3pone*epsb*xipone + 0.29653509649981924*epsc*xipone + -2.2280480847889863*chi2one*epsc*xipone + 6.6841442543669585*chi3pone*epsc*xipone + 0.5930701929996385*epsb*etaone*xipone + -0.5930701929996385*epsc*etaone*xipone + 0.531894201789379*epsb*etapone*xipone + -0.531894201789379*epsc*etapone*xipone + 0.2785060105986233*xippone + 0.10280223253527784*as*xippone + -0.13297355044734474*epsb*xippone + 0.13297355044734474*epsc*xippone + 0.2659471008946895*epsb*etaone*xippone + -0.2659471008946895*epsc*etaone*xippone + 0.1306576264897831*l1one*pow(epsc,2) + 0.48182394158829817*l1pone*pow(epsc,2) + -0.062382885202461706*l4one*pow(epsc,2) + -0.23004832127614688*l4pone*pow(epsc,2) + 0.6210769468876098*l1one*xipone*pow(epsc,2) + 0.5570120211972466*l1pone*xipone*pow(epsc,2) + -0.29653509649981924*l4one*xipone*pow(epsc,2) + -0.2659471008946895*l4pone*xipone*pow(epsc,2) + 0.2785060105986233*l1one*xippone*pow(epsc,2) + -0.13297355044734474*l4one*xippone*pow(epsc,2);
        }

        double S1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.04733475740152116 + 0.011600723452028892*as + 0.04733475740152116*l1one*pow(epsc,2);
        }

        double S1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();


            return 0.2843488842022793 + 0.21589749983616183*as + -0.39656050100795237*epsb + -1.5147122368486774*chi2one*epsb + 4.544136710546032*chi3pone*epsb + 0.39656050100795237*epsc + -1.5147122368486774*chi2one*epsc + 4.544136710546032*chi3pone*epsc + 0.7931210020159047*epsb*etaone + -0.7931210020159047*epsc*etaone + 0.37867805921216935*xipone + 0.09280578761623114*as*xipone + 0.2843488842022793*l1one*pow(epsc,2) + 0.37867805921216935*l1pone*pow(epsc,2) + -0.39656050100795237*l4one*pow(epsc,2) + 0.37867805921216935*l1one*xipone*pow(epsc,2);
        }

        double S1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();


            return 0.6559703043190044 + 0.3219759165015938*as + -1.5890932982243737*epsb + -12.128588768170292*chi2one*epsb + -12.117697894789417*chi2pone*epsb + 36.385766304510874*chi3pone*epsb + 18.176546842184127*chi3ppone*epsb + 1.5890932982243737*epsc + -12.128588768170292*chi2one*epsc + -12.117697894789417*chi2pone*epsc + 36.385766304510874*chi3pone*epsc + 18.176546842184127*chi3ppone*epsc + 3.1781865964487475*epsb*etaone + -3.1781865964487475*epsc*etaone + 6.344968016127237*epsb*etapone + -6.344968016127237*epsc*etapone + 3.032147192042573*xipone + 1.9127915739217567*as*xipone + -3.1724840080636185*epsb*xipone + -12.117697894789417*chi2one*epsb*xipone + 36.353093684368254*chi3pone*epsb*xipone + 3.1724840080636185*epsc*xipone + -12.117697894789417*chi2one*epsc*xipone + 36.353093684368254*chi3pone*epsc*xipone + 6.344968016127237*epsb*etaone*xipone + -6.344968016127237*epsc*etaone*xipone + 1.5147122368486772*xippone + 0.37122315046492455*as*xippone + 0.6559703043190044*l1one*pow(epsc,2) + 2.274791073618234*l1pone*pow(epsc,2) + -1.5890932982243737*l4one*pow(epsc,2) + -3.1724840080636185*l4pone*pow(epsc,2) + 3.032147192042573*l1one*xipone*pow(epsc,2) + 3.0294244736973543*l1pone*xipone*pow(epsc,2) + -3.1724840080636185*l4one*xipone*pow(epsc,2) + 1.5147122368486772*l1one*xippone*pow(epsc,2);
        }
        // }}}

        // B -> D^* form factors
        // {{{
        double A1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.008363045281441904 + -0.003525762321913306*as + 0.008363045281441904*l2one*pow(epsc,2);
        }

        double A1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.0812300174749896 + 0.033452181125767616*epsb + -0.26761744900614093*chi2one*epsb + 0.8028523470184228*chi3pone*epsb + 0.033452181125767616*epsc + -0.26761744900614093*chi3pone*epsc + -0.06690436225153523*epsb*etaone + as*(0.00850018218230231 - 0.028206098575306444*xipone) + 0.06690436225153523*xipone + (0.06690436225153523*l2pone - 0.033452181125767616*l5one + l2one*(0.0812300174749896 + 0.06690436225153523*xipone))*pow(epsc,2);
        }

        double A1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.31721556277653973 + 0.25801570764842313*epsb + -3.134595457211949*chi2one*epsb + -2.1409395920491274*chi2pone*epsb + 9.403786371635848*chi3pone*epsb + 3.211409388073691*chi3ppone*epsb + 0.25801570764842313*epsc + -3.134595457211949*chi3pone*epsc + -1.0704697960245637*chi3ppone*epsc + -0.5160314152968463*epsb*etaone + -0.5352348980122819*epsb*etapone + 0.7836488643029873*xipone + 0.26761744900614093*epsb*xipone + -2.1409395920491274*chi2one*epsb*xipone + 6.422818776147382*chi3pone*epsb*xipone + 0.26761744900614093*epsc*xipone + -2.1409395920491274*chi3pone*epsc*xipone + -0.5352348980122819*epsb*etaone*xipone + as*(0.17093446696648967 + 0.011589260307805627*xipone - 0.11282439430122577*xippone) + 0.26761744900614093*xippone + (0.6498401397999168*l2pone - 0.25801570764842313*l5one - 0.26761744900614093*l5pone + 0.5352348980122819*l2pone*xipone - 0.26761744900614093*l5one*xipone + l2one*(0.31721556277653973 + 0.7836488643029873*xipone + 0.26761744900614093*xippone))*pow(epsc,2);
        }

        double A5_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.005602202787090022 + -0.0023618233360843633*as + 0.005602202787090022*l2one*pow(epsc,2);
        }

        double A5_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.043523865623533844 + 0.012992295758625651*as + -0.0499621450168773*epsb + -0.17927048918688068*chi2one*epsb + 0.537811467560642*chi3pone*epsb + 0.04996214501687729*epsc + 0.17927048918688068*chi2one*epsc + -0.17927048918688068*chi3pone*epsc + 0.0999242900337546*epsb*etaone + 0.0999242900337546*epsc*etaone + 0.04481762229672017*xipone + -0.018894586688674903*as*xipone + 0.043523865623533844*l2one*pow(epsc,2) + 0.04481762229672017*l2pone*pow(epsc,2) + 0.04481762229672017*l3one*pow(epsc,2) + 0.0499621450168773*l5one*pow(epsc,2) + -0.0999242900337546*l6one*pow(epsc,2) + 0.04481762229672017*l2one*xipone*pow(epsc,2);
        }

        double A5_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.11731273819106425 + 0.11394596757274658*as + -0.28823475536695886*epsb + -1.7513046783268447*chi2one*epsb + -1.4341639134950457*chi2pone*epsb + 5.253914034980533*chi3pone*epsb + 2.1512458702425685*chi3ppone*epsb + 0.28823475536695875*epsc + 1.7513046783268447*chi2one*epsc + 1.4341639134950457*chi2pone*epsc + -1.7513046783268447*chi3pone*epsc + -0.7170819567475228*chi3ppone*epsc + 0.5764695107339177*epsb*etaone + 0.5764695107339177*epsc*etaone + 0.7993943202700368*epsb*etapone + 0.7993943202700368*epsc*etapone + 0.43782616958171117*xipone + 0.06614919269165541*as*xipone + -0.3996971601350184*epsb*xipone + -1.4341639134950457*chi2one*epsb*xipone + 4.302491740485137*chi3pone*epsb*xipone + 0.3996971601350183*epsc*xipone + 1.4341639134950457*chi2one*epsc*xipone + -1.4341639134950457*chi3pone*epsc*xipone + 0.7993943202700368*epsb*etaone*xipone + 0.7993943202700368*epsc*etaone*xipone + 0.1792704891868807*xippone + -0.07557834675469961*as*xippone + 0.11731273819106425*l2one*pow(epsc,2) + 0.34819092498827076*l2pone*pow(epsc,2) + 0.43782616958171117*l3one*pow(epsc,2) + 0.3585409783737614*l3pone*pow(epsc,2) + 0.28823475536695886*l5one*pow(epsc,2) + 0.3996971601350184*l5pone*pow(epsc,2) + -0.976166670868936*l6one*pow(epsc,2) + -0.7993943202700368*l6pone*pow(epsc,2) + 0.43782616958171117*l2one*xipone*pow(epsc,2) + 0.3585409783737614*l2pone*xipone*pow(epsc,2) + 0.3585409783737614*l3one*xipone*pow(epsc,2) + 0.3996971601350184*l5one*xipone*pow(epsc,2) + -0.7993943202700368*l6one*xipone*pow(epsc,2) + 0.1792704891868807*l2one*xippone*pow(epsc,2);
        }

        double V4_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.007502586789045188 + 0.006840445991871229*as + 0.007502586789045188*epsb + 0.007502586789045188*epsc + -0.015005173578090376*epsb*etaone + 0.007502586789045188*l2one*pow(epsc,2) + -0.007502586789045188*l5one*pow(epsc,2);
        }

        double V4_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.05963010511442658 + 0.06603935741648106*as + 0.05963010511442658*epsb + -0.24008277724944604*chi2one*epsb + 0.720248331748338*chi3pone*epsb + 0.05963010511442658*epsc + -0.24008277724944604*chi3pone*epsc + -0.11926021022885316*epsb*etaone + -0.12004138862472302*epsb*etapone + 0.060020694312361504*xipone + 0.05472356793496983*as*xipone + 0.060020694312361504*epsb*xipone + 0.060020694312361504*epsc*xipone + -0.12004138862472301*epsb*etaone*xipone + 0.05963010511442658*l2one*pow(epsc,2) + 0.06002069431236151*l2pone*pow(epsc,2) + -0.05963010511442658*l5one*pow(epsc,2) + -0.06002069431236151*l5pone*pow(epsc,2) + 0.060020694312361504*l2one*xipone*pow(epsc,2) + -0.060020694312361504*l5one*xipone*pow(epsc,2);
        }

        double V4_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.1543810939148951 + 0.1663860668436065*as + 0.1543810939148951*epsb + -2.3883289181605427*chi2one*epsb + -1.9206622179955681*chi2pone*epsb + 7.164986754481628*chi3pone*epsb + 2.880993326993352*chi3ppone*epsb + 0.1543810939148951*epsc + -2.3883289181605427*chi3pone*epsc + -0.9603311089977841*chi3ppone*epsc + -0.3087621878297902*epsb*etaone + -1.1941644590802714*epsb*etapone + -0.48016555449889203*epsb*etappone + 0.5970822295401357*xipone + 0.6377619952017881*as*xipone + 0.5970822295401357*epsb*xipone + -1.9206622179955684*chi2one*epsb*xipone + 5.761986653986704*chi3pone*epsb*xipone + 0.5970822295401357*epsc*xipone + -1.9206622179955684*chi3pone*epsc*xipone + -1.1941644590802714*epsb*etaone*xipone + -0.9603311089977842*epsb*etapone*xipone + 0.24008277724944602*xippone + 0.21889427173987933*as*xippone + 0.24008277724944602*epsb*xippone + 0.24008277724944602*epsc*xippone + -0.48016555449889203*epsb*etaone*xippone + 0.1543810939148951*l2one*pow(epsc,2) + 0.4770408409154127*l2pone*pow(epsc,2) + -0.1543810939148951*l5one*pow(epsc,2) + -0.4770408409154127*l5pone*pow(epsc,2) + 0.5970822295401357*l2one*xipone*pow(epsc,2) + 0.4801655544988921*l2pone*xipone*pow(epsc,2) + -0.5970822295401357*l5one*xipone*pow(epsc,2) + -0.4801655544988921*l5pone*xipone*pow(epsc,2) + 0.24008277724944602*l2one*xippone*pow(epsc,2) + -0.24008277724944602*l5one*xippone*pow(epsc,2);
        }

        double P1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.03148856751270911 + -0.0022676800849227457*as + -0.014123119863775367*epsb + 0.014123119863775367*epsc + 0.028246239727550734*epsb*etaone + 0.028246239727550734*epsc*etaone + (0.03148856751270911*l2one + 0.014123119863775367*l5one - 0.028246239727550734*l6one)*pow(epsc,2);
        }

        double P1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.26142777033909864 + -0.11725448401958674*epsb + -1.0076341604066916*chi2one*epsb + 3.022902481220075*chi3pone*epsb + 0.11725448401958674*epsc + 1.0076341604066916*chi2one*epsc + -1.0076341604066916*chi3pone*epsc + 0.23450896803917348*epsb*etaone + 0.23450896803917348*epsc*etaone + 0.22596991782040585*epsb*etapone + 0.22596991782040585*epsc*etapone + as*(-0.04568809944563187 - 0.018141440679381962*xipone) + 0.2519085401016729*xipone + -0.11298495891020292*epsb*xipone + 0.11298495891020292*epsc*xipone + 0.22596991782040585*epsb*etaone*xipone + 0.22596991782040585*epsc*etaone*xipone + (0.2519085401016729*l2pone + 0.2519085401016729*l3one + 0.11725448401958674*l5one + 0.11298495891020292*l5pone - 0.34749392694937636*l6one - 0.22596991782040585*l6pone + l2one*(0.26142777033909864 + 0.2519085401016729*xipone) + 0.11298495891020292*l5one*xipone - 0.22596991782040585*l6one*xipone)*pow(epsc,2);
        }

        double P1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.7662176967652994 + -0.34366073873619646*epsb + -10.38095697166454*chi2one*epsb + -8.061073283253533*chi2pone*epsb + 31.142870914993626*chi3pone*epsb + 12.091609924880299*chi3ppone*epsb + 0.34366073873619657*epsc + 10.38095697166454*chi2one*epsc + 8.061073283253533*chi2pone*epsc + -10.38095697166454*chi3pone*epsc + -4.0305366416267665*chi3ppone*epsc + 0.6873214774723929*epsb*etaone + 0.687321477472393*epsc*etaone + 2.3280115799542*epsb*etapone + 2.3280115799542*epsc*etapone + 0.9038796712816234*epsb*etappone + 0.9038796712816234*epsc*etappone + 2.595239242916135*xipone + -1.1640057899771*epsb*xipone + -8.061073283253533*chi2one*epsb*xipone + 24.183219849760597*chi3pone*epsb*xipone + 1.1640057899771*epsc*xipone + 8.061073283253533*chi2one*epsc*xipone + -8.061073283253533*chi3pone*epsc*xipone + 2.3280115799542*epsb*etaone*xipone + 2.3280115799542*epsc*etaone*xipone + 1.8077593425632468*epsb*etapone*xipone + 1.8077593425632468*epsc*etapone*xipone + as*(-0.5977425693276348 - 0.401787676923819*xipone - 0.07256576271752786*xippone) + 1.0076341604066916*xippone + -0.4519398356408117*epsb*xippone + 0.4519398356408117*epsc*xippone + 0.9038796712816234*epsb*etaone*xippone + 0.9038796712816234*epsc*etaone*xippone + (2.091422162712789*l2pone + 2.595239242916135*l3one + 2.0152683208133833*l3pone + 0.34366073873619646*l5one + 0.9380358721566939*l5pone - 1.8513272674494927*l6one - 2.779951415595011*l6pone + 2.0152683208133833*l2pone*xipone + 2.0152683208133833*l3one*xipone + 1.1640057899771*l5one*xipone + 0.9038796712816234*l5pone*xipone - 3.2318912512358233*l6one*xipone - 1.8077593425632468*l6pone*xipone + 0.4519398356408117*l5one*xippone - 0.9038796712816234*l6one*xippone + l2one*(0.7662176967652994 + 2.595239242916135*xipone + 1.0076341604066916*xippone))*pow(epsc,2);
        }
        // }}}

        // B^* -> D form factors
        // {{{
        double P2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.0438486230520249 + -0.002055679483844757*as + -0.021080868288500247*epsb + 0.021080868288500247*epsc + -0.042161736577000494*epsb*etaone + -0.042161736577000494*epsc*etaone + (0.0438486230520249*l1one - 0.021080868288500247*l4one)*pow(epsc,2);
        }

        double P2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.3296001313229333 + -0.1584600945860323*epsb + 1.4031559376647968*chi2one*epsb + -1.4031559376647968*chi3pone*epsb + 0.1584600945860323*epsc + -1.4031559376647968*chi2one*epsc + 4.209467812994389*chi3pone*epsc + -0.3169201891720646*epsb*etaone + -0.3169201891720646*epsc*etaone + -0.33729389261600395*epsb*etapone + -0.33729389261600395*epsc*etapone + as*(-0.05126982691808179 - 0.016445435870758057*xipone) + 0.3507889844161992*xipone + -0.16864694630800198*epsb*xipone + 0.16864694630800198*epsc*xipone + -0.33729389261600395*epsb*etaone*xipone + -0.33729389261600395*epsc*etaone*xipone + (0.3507889844161992*l1pone - 0.1584600945860323*l4one - 0.16864694630800198*l4pone + l1one*(0.3296001313229333 + 0.3507889844161992*xipone) - 0.16864694630800198*l4one*xipone)*pow(epsc,2);
        }

        double P2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.8529477601473068 + -0.41006713864886873*epsb + 13.35351607766346*chi2one*epsb + 11.225247501318373*chi2pone*epsb + -13.35351607766346*chi3pone*epsb + -5.6126237506591865*chi3ppone*epsb + 0.41006713864886873*epsc + -13.35351607766346*chi2one*epsc + -11.225247501318373*chi2pone*epsc + 40.06054823299038*chi3pone*epsc + 16.83787125197756*chi3ppone*epsc + -0.8201342772977375*epsb*etaone + -0.8201342772977375*epsc*etaone + -3.209949298608525*epsb*etapone + -3.209949298608525*epsc*etapone + -1.3491755704640158*epsb*etappone + -1.3491755704640158*epsc*etappone + 3.338379019415865*xipone + -1.6049746493042625*epsb*xipone + 11.225247501318373*chi2one*epsb*xipone + -11.225247501318373*chi3pone*epsb*xipone + 1.6049746493042625*epsc*xipone + -11.225247501318373*chi2one*epsc*xipone + 33.67574250395512*chi3pone*epsc*xipone + -3.209949298608525*epsb*etaone*xipone + -3.209949298608525*epsc*etaone*xipone + -2.6983511409280316*epsb*etapone*xipone + -2.6983511409280316*epsc*etapone*xipone + as*(-0.7534023367765568 - 0.44304948708617065*xipone - 0.06578174348303224*xippone) + 1.4031559376647966*xippone + -0.6745877852320079*epsb*xippone + 0.6745877852320079*epsc*xippone + -1.3491755704640158*epsb*etaone*xippone + -1.3491755704640158*epsc*etaone*xippone + (2.6368010505834665*l1pone - 0.41006713864886873*l4one - 1.2676807566882586*l4pone + 2.8063118753295933*l1pone*xipone - 1.6049746493042625*l4one*xipone - 1.3491755704640158*l4pone*xipone - 0.6745877852320079*l4one*xippone + l1one*(0.8529477601473068 + 3.338379019415865*xipone + 1.4031559376647966*xippone))*pow(epsc,2);
        }

        double V5_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.009099008915484857 + 0.008295975883519244*as + 0.009099008915484857*epsb + 0.009099008915484857*epsc + -0.018198017830969714*epsc*etaone + 0.009099008915484857*l1one*pow(epsc,2) + -0.009099008915484857*l4one*pow(epsc,2);
        }

        double V5_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.06695557218376086 + 0.07520189758366089*as + 0.06695557218376086*epsb + -0.2911682852955154*chi3pone*epsb + 0.06695557218376086*epsc + -0.2911682852955154*chi2one*epsc + 0.8735048558865462*chi3pone*epsc + -0.13391114436752172*epsc*etaone + -0.1455841426477577*epsc*etapone + 0.07279207132387885*xipone + 0.06636780706815396*as*xipone + 0.07279207132387885*epsb*xipone + 0.07279207132387885*epsc*xipone + -0.1455841426477577*epsc*etaone*xipone + 0.06695557218376086*l1one*pow(epsc,2) + 0.07279207132387885*l1pone*pow(epsc,2) + -0.06695557218376086*l4one*pow(epsc,2) + -0.07279207132387885*l4pone*pow(epsc,2) + 0.07279207132387885*l1one*xipone*pow(epsc,2) + -0.07279207132387885*l4one*xipone*pow(epsc,2);
        }

        double V5_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.1574600063604854 + 0.1663038473068933*as + 0.1574600063604854*epsb + -2.7249148804713785*chi3pone*epsb + -1.1646731411820617*chi3ppone*epsb + 0.1574600063604854*epsc + -2.7249148804713785*chi2one*epsc + -2.3293462823641233*chi2pone*epsc + 8.174744641414136*chi3pone*epsc + 3.494019423546185*chi3ppone*epsc + -0.3149200127209708*epsc*etaone + -1.3624574402356893*epsc*etapone + -0.5823365705910308*epsc*etappone + 0.6812287201178446*xipone + 0.7343507948055951*as*xipone + 0.6812287201178446*epsb*xipone + -2.3293462823641238*chi3pone*epsb*xipone + 0.6812287201178446*epsc*xipone + -2.3293462823641238*chi2one*epsc*xipone + 6.988038847092371*chi3pone*epsc*xipone + -1.3624574402356893*epsc*etaone*xipone + -1.1646731411820619*epsc*etapone*xipone + 0.2911682852955154*xippone + 0.2654712282726158*as*xippone + 0.2911682852955154*epsb*xippone + 0.2911682852955154*epsc*xippone + -0.5823365705910308*epsc*etaone*xippone + 0.1574600063604854*l1one*pow(epsc,2) + 0.535644577470087*l1pone*pow(epsc,2) + -0.1574600063604854*l4one*pow(epsc,2) + -0.535644577470087*l4pone*pow(epsc,2) + 0.6812287201178446*l1one*xipone*pow(epsc,2) + 0.5823365705910309*l1pone*xipone*pow(epsc,2) + -0.6812287201178446*l4one*xipone*pow(epsc,2) + -0.5823365705910309*l4pone*xipone*pow(epsc,2) + 0.2911682852955154*l1one*xippone*pow(epsc,2) + -0.2911682852955154*l4one*xippone*pow(epsc,2);
        }

        double A2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.013497592386244992 + -0.005690427478321802*as + 0.013497592386244992*l1one*pow(epsc,2);
        }

        double A2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.11544627665098033 + 0.05399036954497997*epsb + -0.43192295635983974*chi3pone*epsb + 0.05399036954497997*epsc + -0.43192295635983974*chi2one*epsc + 1.2957688690795193*chi3pone*epsc + -0.10798073908995993*epsc*etaone + as*(0.020319079912870634 - 0.045523419826574416*xipone) + 0.10798073908995993*xipone + (0.10798073908995993*l1pone - 0.05399036954497997*l4one + l1one*(0.11544627665098033 + 0.10798073908995993*xipone))*pow(epsc,2);
        }

        double A2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.392953572653608 + 0.3538043675139614*epsb + -4.55812676555105*chi3pone*epsb + -1.727691825439359*chi3ppone*epsb + 0.3538043675139614*epsc + -4.55812676555105*chi2one*epsc + -3.455383650878718*chi2pone*epsc + 13.67438029665315*chi3pone*epsc + 5.183075476318077*chi3ppone*epsc + -0.7076087350279228*epsc*etaone + -0.8638459127196795*epsc*etapone + 1.1395316913877624*xipone + 0.43192295635983974*epsb*xipone + -3.455383650878718*chi3pone*epsb*xipone + 0.43192295635983974*epsc*xipone + -3.455383650878718*chi2one*epsc*xipone + 10.366150952636154*chi3pone*epsc*xipone + -0.8638459127196795*epsc*etaone*xipone + as*(0.2460383882047104 + 0.07150579964981622*xipone - 0.18209367930629766*xippone) + 0.43192295635983974*xippone + (0.9235702132078426*l1pone - 0.3538043675139614*l4one - 0.43192295635983974*l4pone + 0.8638459127196795*l1pone*xipone - 0.43192295635983974*l4one*xipone + l1one*(0.392953572653608 + 1.1395316913877624*xipone + 0.43192295635983974*xippone))*pow(epsc,2);
        }

        double A6_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.00977922297982081 + -0.004122806317496859*as + 0.00977922297982081*l1one*pow(epsc,2);
        }

        double A6_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.06472590811759746 + -0.0813639089842952*epsb + 0.31293513535426587*chi2one*epsb + -0.31293513535426587*chi3pone*epsb + 0.0813639089842952*epsc + -0.31293513535426587*chi2one*epsc + 0.9388054060627977*chi3pone*epsc + -0.16272781796859037*epsb*etaone + -0.1627278179685904*epsc*etaone + as*(0.02286243511586715 - 0.03298245053997487*xipone) + 0.07823378383856647*xipone + (0.07823378383856647*l1pone - 0.0813639089842952*l4one + l1one*(0.06472590811759746 + 0.07823378383856647*xipone))*pow(epsc,2);
        }

        double A6_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.14122128288513747 + -0.375796859081217*epsb + 2.697099330471651*chi2one*epsb + 2.503481082834127*chi2pone*epsb + -2.697099330471651*chi3pone*epsb + -1.2517405414170635*chi3ppone*epsb + 0.37579685908121685*epsc + -2.697099330471651*chi2one*epsc + -2.503481082834127*chi2pone*epsc + 8.091297991414953*chi3pone*epsc + 3.7552216242511904*chi3ppone*epsc + -0.7515937181624335*epsb*etaone + -0.7515937181624337*epsc*etaone + -1.301822543748723*epsb*etapone + -1.3018225437487234*epsc*etapone + 0.6742748326179128*xipone + -0.6509112718743617*epsb*xipone + 2.503481082834127*chi2one*epsb*xipone + -2.503481082834127*chi3pone*epsb*xipone + 0.6509112718743617*epsc*xipone + -2.503481082834127*chi2one*epsc*xipone + 7.510443248502381*chi3pone*epsc*xipone + -1.301822543748723*epsb*etaone*xipone + -1.3018225437487234*epsc*etaone*xipone + as*(0.1351402435746817 + 0.11693457984698744*xipone - 0.13192980215989944*xippone) + 0.31293513535426587*xippone + (0.5178072649407798*l1pone - 0.37579685908121685*l4one - 0.6509112718743617*l4pone + 0.6258702707085317*l1pone*xipone - 0.6509112718743617*l4one*xipone + l1one*(0.14122128288513747 + 0.6742748326179128*xipone + 0.31293513535426587*xippone))*pow(epsc,2);
        }
        // }}}

        // B^* -> D^* form factors
        // {{{
        double S2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.03876949172451259 + 0.009501562415472311*as + 0.03876949172451259*l2one*pow(epsc,2);
        }

        double S2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.2783004793038407 + 0.19272149041023498*as + -0.3431363158965911*epsb + -1.240623735184403*chi3pone*epsb + 0.3431363158965911*epsc + -1.240623735184403*chi3pone*epsc + 0.31015593379610074*xipone + 0.0760124993237785*as*xipone + 0.2783004793038407*l2one*pow(epsc,2) + 0.31015593379610074*l2pone*pow(epsc,2) + -0.3431363158965911*l5one*pow(epsc,2) + 0.31015593379610074*l2one*xipone*pow(epsc,2);
        }

        double S2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.7635413763985075 + 0.4799658571364192*as + -1.776875501773842*epsb + -11.386862808091708*chi3pone*epsb + -4.962494940737613*chi3ppone*epsb + 1.776875501773842*epsc + -11.386862808091708*chi3pone*epsc + -4.962494940737613*chi3ppone*epsc + 2.846715702022927*xipone + 1.6937969219294364*as*xipone + -2.7450905271727284*epsb*xipone + -9.924989881475225*chi3pone*epsb*xipone + 2.7450905271727284*epsc*xipone + -9.924989881475225*chi3pone*epsc*xipone + 1.2406237351844032*xippone + 0.304049997295114*as*xippone + 0.7635413763985075*l2one*pow(epsc,2) + 2.2264038344307258*l2pone*pow(epsc,2) + -1.776875501773842*l5one*pow(epsc,2) + -2.7450905271727284*l5pone*pow(epsc,2) + 2.846715702022927*l2one*xipone*pow(epsc,2) + 2.4812474703688063*l2pone*xipone*pow(epsc,2) + -2.7450905271727284*l5one*xipone*pow(epsc,2) + 1.2406237351844032*l2one*xippone*pow(epsc,2);
        }

        double S3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.027414170501558595 + 0.006718619215847705*as + 0.027414170501558595*l2one*pow(epsc,2);
        }

        double S3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.1967881561232122 + 0.13627467274945537*as + -0.24263401584184888*epsb + 0.8772534560498749*chi2one*epsb + -0.8772534560498749*chi3pone*epsb + 0.24263401584184888*epsc + 0.8772534560498749*chi2one*epsc + -0.8772534560498749*chi3pone*epsc + -0.48526803168369775*epsb*etaone + 0.48526803168369775*epsc*etaone + 0.21931336401246873*xipone + 0.05374895372678164*as*xipone + 0.1967881561232122*l2one*pow(epsc,2) + 0.21931336401246873*l2pone*pow(epsc,2) + 0.21931336401246873*l3one*pow(epsc,2) + 0.24263401584184888*l5one*pow(epsc,2) + 0.21931336401246873*l2one*xipone*pow(epsc,2) + -0.48526803168369775*l6one*pow(epsc,2);
        }

        double S3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.539905284967895 + 0.33938711231917595*as + -1.256440716628533*epsb + 8.051727908042539*chi2one*epsb + 7.018027648399001*chi2pone*epsb + -8.051727908042539*chi3pone*epsb + -3.5090138241995006*chi3ppone*epsb + 1.256440716628533*epsc + 8.051727908042539*chi2one*epsc + 7.018027648399001*chi2pone*epsc + -8.051727908042539*chi3pone*epsc + -3.5090138241995006*chi3ppone*epsc + -2.512881433257066*epsb*etaone + 2.512881433257066*epsc*etaone + -3.882144253469582*epsb*etapone + 3.882144253469582*epsc*etapone + 2.0129319770106346*xipone + 1.1976952894492057*as*xipone + -1.941072126734791*epsb*xipone + 7.018027648399001*chi2one*epsb*xipone + -7.018027648399001*chi3pone*epsb*xipone + 1.941072126734791*epsc*xipone + 7.018027648399001*chi2one*epsc*xipone + -7.018027648399001*chi3pone*epsc*xipone + -3.882144253469582*epsb*etaone*xipone + 3.882144253469582*epsc*etaone*xipone + 0.8772534560498751*xippone + 0.21499581490712652*as*xippone + 0.539905284967895*l2one*pow(epsc,2) + 1.5743052489856972*l2pone*pow(epsc,2) + 2.0129319770106346*l3one*pow(epsc,2) + 1.7545069120997503*l3pone*pow(epsc,2) + 1.2564407166285327*l5one*pow(epsc,2) + 1.941072126734791*l5pone*pow(epsc,2) + -3.882144253469582*l6pone*pow(epsc,2) + 2.0129319770106346*l2one*xipone*pow(epsc,2) + 1.7545069120997503*l2pone*xipone*pow(epsc,2) + 1.7545069120997503*l3one*xipone*pow(epsc,2) + 1.941072126734791*l5one*xipone*pow(epsc,2) + 0.8772534560498751*l2one*xippone*pow(epsc,2) + -4.453953559991857*l6one*pow(epsc,2) + -3.882144253469582*xipone*l6one*pow(epsc,2);
        }

        double P3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.03684032631637572 + -0.0025546933818376485*as + -0.016649718028464284*epsb + 0.016649718028464284*epsc + (0.03684032631637572*l2one - 0.016649718028464284*l5one)*pow(epsc,2);
        }

        double P3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.3236923482372957 + -0.14629040687207026*epsb + -1.178890442124023*chi3pone*epsb + 0.14629040687207026*epsc + -1.178890442124023*chi3pone*epsc + as*(-0.053731167406765354 - 0.02043754705470119*xipone) + 0.29472261053100574*xipone + -0.13319774422771427*epsb*xipone + 0.13319774422771427*epsc*xipone + (0.29472261053100574*l2pone - 0.14629040687207026*l5one - 0.13319774422771427*l5pone + l2one*(0.3236923482372957 + 0.29472261053100574*xipone) - 0.13319774422771427*l5one*xipone)*pow(epsc,2);
        }

        double P3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 1.0139125139266254 + -0.4582303999545689*epsb + -12.715936027841508*chi3pone*epsb + -4.715561768496092*chi3ppone*epsb + 0.4582303999545689*epsc + -12.715936027841508*chi3pone*epsc + -4.715561768496092*chi3ppone*epsc + 3.178984006960377*xipone + -1.4367187434319904*epsb*xipone + -9.431123536992184*chi3pone*epsb*xipone + 1.4367187434319904*epsc*xipone + -9.431123536992184*chi3pone*epsc*xipone + as*(-0.7189809143133374 - 0.47072443336352526*xipone - 0.08175018821880474*xippone) + 1.178890442124023*xippone + -0.5327909769108571*epsb*xippone + 0.5327909769108571*epsc*xippone + (2.5895387858983656*l2pone - 0.4582303999545689*l5one - 1.170323254976562*l5pone + 2.357780884248046*l2pone*xipone - 1.4367187434319904*l5one*xipone - 1.0655819538217142*l5pone*xipone - 0.5327909769108571*l5one*xippone + l2one*(1.0139125139266254 + 3.178984006960377*xipone + 1.178890442124023*xippone))*pow(epsc,2);
        }

        double V2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.007335151340155247 + 0.002658939321365706*as + -0.0033150683970844777*epsb + 0.0033150683970844777*epsc + 0.007335151340155247*l2one*pow(epsc,2) + -0.0033150683970844777*l5one*pow(epsc,2);
        }

        double V2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.06109473379049465 + 0.028094157378472037*as + -0.027611321406338075*epsb + -0.23472484288496792*chi3pone*epsb + 0.027611321406338075*epsc + -0.23472484288496792*chi3pone*epsc + 0.05868121072124198*xipone + 0.02127151457092565*as*xipone + -0.026520547176675822*epsb*xipone + 0.026520547176675822*epsc*xipone + 0.06109473379049465*l2one*pow(epsc,2) + 0.05868121072124198*l2pone*pow(epsc,2) + -0.027611321406338075*l5one*pow(epsc,2) + -0.026520547176675822*l5pone*pow(epsc,2) + 0.05868121072124198*l2one*xipone*pow(epsc,2) + -0.026520547176675822*l5one*xipone*pow(epsc,2);
        }

        double V2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.16740988852041994 + 0.04231490821343695*as + -0.07565968376894236*epsb + -2.4244811670657653*chi3pone*epsb + -0.9388993715398717*chi3ppone*epsb + 0.07565968376894236*epsc + -2.4244811670657653*chi3pone*epsc + -0.9388993715398717*chi3ppone*epsc + 0.6061202917664413*xipone + 0.2672962881696276*as*xipone + -0.2739316656040563*epsb*xipone + -1.8777987430797434*chi3pone*epsb*xipone + 0.2739316656040563*epsc*xipone + -1.8777987430797434*chi3pone*epsc*xipone + 0.23472484288496792*xippone + 0.08508605828370261*as*xippone + -0.1060821887067033*epsb*xippone + 0.1060821887067033*epsc*xippone + 0.16740988852041994*l2one*pow(epsc,2) + 0.4887578703239573*l2pone*pow(epsc,2) + -0.07565968376894236*l5one*pow(epsc,2) + -0.22089057125070463*l5pone*pow(epsc,2) + 0.6061202917664413*l2one*xipone*pow(epsc,2) + 0.46944968576993584*l2pone*xipone*pow(epsc,2) + -0.2739316656040563*l5one*xipone*pow(epsc,2) + -0.2121643774134066*l5pone*xipone*pow(epsc,2) + 0.23472484288496792*l2one*xippone*pow(epsc,2) + -0.1060821887067033*l5one*xippone*pow(epsc,2);
        }

        double V3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.005186735253653369 + 0.0018801540249012481*as + -0.0023441073436756533*epsb + 0.0023441073436756533*epsc + -0.0046882146873513065*epsb*etaone + 0.0046882146873513065*epsc*etaone + 0.005186735253653369*l2one*pow(epsc,2) + 0.0023441073436756533*l5one*pow(epsc,2) + -0.0046882146873513065*l6one*pow(epsc,2);
        }

        double V3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.043200500558045696 + 0.019865569194039662*as + -0.01952415260394294*epsb + 0.1659755281169078*chi2one*epsb + -0.1659755281169078*chi3pone*epsb + 0.01952415260394294*epsc + 0.1659755281169078*chi2one*epsc + -0.1659755281169078*chi3pone*epsc + -0.03904830520788588*epsb*etaone + 0.03904830520788588*epsc*etaone + -0.03750571749881045*epsb*etapone + 0.03750571749881045*epsc*etapone + 0.04149388202922695*xipone + 0.015041232199209985*as*xipone + -0.018752858749405226*epsb*xipone + 0.018752858749405226*epsc*xipone + -0.03750571749881045*epsb*etaone*xipone + 0.03750571749881045*epsc*etaone*xipone + 0.043200500558045696*l2one*pow(epsc,2) + 0.04149388202922695*l2pone*pow(epsc,2) + 0.04149388202922695*l3one*pow(epsc,2) + 0.01952415260394294*l5one*pow(epsc,2) + 0.018752858749405226*l5pone*pow(epsc,2) + -0.057801163957291096*l6one*pow(epsc,2) + -0.03750571749881045*l6pone*pow(epsc,2) + 0.04149388202922695*l2one*xipone*pow(epsc,2) + 0.018752858749405226*l5one*xipone*pow(epsc,2) + -0.03750571749881045*l6one*xipone*pow(epsc,2);
        }

        double V3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.11837666741047298 + 0.029921158543007614*as + -0.05349947545544891*epsb + 1.714367074091277*chi2one*epsb + 1.327804224935262*chi2pone*epsb + -1.7143670740912778*chi3pone*epsb + -0.6639021124676311*chi3ppone*epsb + 0.05349947545544891*epsc + 1.714367074091277*chi2one*epsc + 1.327804224935262*chi2pone*epsc + -1.7143670740912778*chi3pone*epsc + -0.6639021124676311*chi3ppone*epsc + -0.10699895091089782*epsb*etaone + 0.10699895091089782*epsc*etaone + -0.38739787666070785*epsb*etapone + 0.38739787666070785*epsc*etapone + -0.15002286999524184*epsb*etappone + 0.15002286999524184*epsc*etappone + 0.42859176852281944*xipone + 0.18900701795073727*as*xipone + -0.19369893833035393*epsb*xipone + 1.327804224935262*chi2one*epsb*xipone + -1.3278042249352622*chi3pone*epsb*xipone + 0.19369893833035393*epsc*xipone + 1.327804224935262*chi2one*epsc*xipone + -1.3278042249352622*chi3pone*epsc*xipone + -0.3873978766607078*epsb*etaone*xipone + 0.38739787666070785*epsc*etaone*xipone + -0.3000457399904837*epsb*etapone*xipone + 0.3000457399904837*epsc*etapone*xipone + 0.16597552811690777*xippone + 0.060164928796839934*as*xippone + -0.07501143499762092*epsb*xippone + 0.07501143499762092*epsc*xippone + -0.15002286999524184*epsb*etaone*xippone + 0.15002286999524184*epsc*etaone*xippone + 0.11837666741047298*l2one*pow(epsc,2) + 0.34560400446436557*l2pone*pow(epsc,2) + 0.42859176852281927*l3one*pow(epsc,2) + 0.3319510562338155*l3pone*pow(epsc,2) + 0.05349947545544891*l5one*pow(epsc,2) + 0.1561932208315435*l5pone*pow(epsc,2) + -0.30069788924125174*l6one*pow(epsc,2) + -0.46240931165832877*l6pone*pow(epsc,2) + 0.42859176852281944*l2one*xipone*pow(epsc,2) + 0.33195105623381554*l2pone*xipone*pow(epsc,2) + 0.3319510562338155*l3one*xipone*pow(epsc,2) + 0.19369893833035393*l5one*xipone*pow(epsc,2) + 0.15002286999524184*l5pone*xipone*pow(epsc,2) + -0.5374207466559496*l6one*xipone*pow(epsc,2) + -0.3000457399904837*l6pone*xipone*pow(epsc,2) + 0.16597552811690777*l2one*xippone*pow(epsc,2) + 0.07501143499762092*l5one*xippone*pow(epsc,2) + -0.15002286999524184*l6one*xippone*pow(epsc,2);
        }

        double V6_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.006543299374599414 + 0.005965820488201019*as + 0.006543299374599414*epsb + 0.006543299374599414*epsc + 0.013086598749198828*epsc*etaone + 0.006543299374599414*l2one*pow(epsc,2) + 0.006543299374599414*l5one*pow(epsc,2) + -0.013086598749198828*l6one*pow(epsc,2);
        }

        double V6_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.05449937088744518 + 0.05986905570747457*as + 0.05449937088744518*epsb + -0.20938557998718124*chi3pone*epsb + 0.05449937088744518*epsc + 0.20938557998718124*chi2one*epsc + -0.20938557998718124*chi3pone*epsc + 0.10899874177489036*epsc*etaone + 0.10469278999359062*epsc*etapone + 0.05234639499679531*xipone + 0.04772656390560815*as*xipone + 0.05234639499679531*epsb*xipone + 0.05234639499679531*epsc*xipone + 0.10469278999359062*epsc*etaone*xipone + 0.05449937088744518*l2one*pow(epsc,2) + 0.05234639499679531*l2pone*pow(epsc,2) + 0.05234639499679531*l3one*pow(epsc,2) + 0.05449937088744518*l5one*pow(epsc,2) + 0.05234639499679531*l5pone*pow(epsc,2) + -0.16134513677168566*l6one*pow(epsc,2) + -0.10469278999359062*l6pone*pow(epsc,2) + 0.05234639499679531*l2one*xipone*pow(epsc,2) + 0.05234639499679531*l5one*xipone*pow(epsc,2) + -0.10469278999359062*l6one*xipone*pow(epsc,2);
        }

        double V6_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.1493374803135605 + 0.16238989364956619*as + 0.1493374803135605*epsb + -2.1627510283726084*chi3pone*epsb + -0.837542319948725*chi3ppone*epsb + 0.1493374803135605*epsc + 2.1627510283726084*chi2one*epsc + 1.67508463989745*chi2pone*epsc + -2.1627510283726084*chi3pone*epsc + -0.837542319948725*chi3ppone*epsc + 0.298674960627121*epsc*etaone + 1.0813755141863042*epsc*etapone + 0.4187711599743625*epsc*etappone + 0.5406877570931521*xipone + 0.5744055734710128*as*xipone + 0.5406877570931521*epsb*xipone + -1.67508463989745*chi3pone*epsb*xipone + 0.5406877570931521*epsc*xipone + 1.67508463989745*chi2one*epsc*xipone + -1.67508463989745*chi3pone*epsc*xipone + 1.0813755141863042*epsc*etaone*xipone + 0.837542319948725*epsc*etapone*xipone + 0.20938557998718124*xippone + 0.1909062556224326*as*xippone + 0.20938557998718124*epsb*xippone + 0.20938557998718124*epsc*xippone + 0.4187711599743625*epsc*etaone*xippone + 0.1493374803135605*l2one*pow(epsc,2) + 0.4359949670995614*l2pone*pow(epsc,2) + 0.5406877570931521*l3one*pow(epsc,2) + 0.4187711599743625*l3pone*pow(epsc,2) + 0.1493374803135605*l5one*pow(epsc,2) + 0.4359949670995614*l5pone*pow(epsc,2) + -0.839362717720273*l6one*pow(epsc,2) + -1.2907610941734853*l6pone*pow(epsc,2) + 0.5406877570931521*l2one*xipone*pow(epsc,2) + 0.4187711599743625*l2pone*xipone*pow(epsc,2) + 0.4187711599743625*l3one*xipone*pow(epsc,2) + 0.5406877570931521*l5one*xipone*pow(epsc,2) + 0.4187711599743625*l5pone*xipone*pow(epsc,2) + -1.5001466741606666*l6one*xipone*pow(epsc,2) + -0.837542319948725*l6pone*xipone*pow(epsc,2) + 0.20938557998718124*l2one*xippone*pow(epsc,2) + 0.20938557998718124*l5one*xippone*pow(epsc,2) + -0.4187711599743625*l6one*xippone*pow(epsc,2);
        }

        double V7_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.006543299374599414 + 0.005965820488201019*as + 0.006543299374599414*epsb + 0.006543299374599414*epsc + 0.013086598749198828*epsb*etaone + 0.006543299374599414*l2one*pow(epsc,2) + -0.006543299374599414*l5one*pow(epsc,2);
        }

        double V7_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.05449937088744518 + 0.05986905570747457*as + 0.05449937088744518*epsb + 0.20938557998718124*chi2one*epsb + -0.20938557998718124*chi3pone*epsb + 0.05449937088744518*epsc + -0.20938557998718124*chi3pone*epsc + 0.10899874177489036*epsb*etaone + 0.10469278999359062*epsb*etapone + 0.05234639499679531*xipone + 0.04772656390560815*as*xipone + 0.05234639499679531*epsb*xipone + 0.05234639499679531*epsc*xipone + 0.10469278999359062*epsb*etaone*xipone + 0.05449937088744518*l2one*pow(epsc,2) + 0.05234639499679531*l2pone*pow(epsc,2) + -0.05449937088744518*l5one*pow(epsc,2) + -0.05234639499679531*l5pone*pow(epsc,2) + 0.05234639499679531*l2one*xipone*pow(epsc,2) + -0.05234639499679531*l5one*xipone*pow(epsc,2);
        }

        double V7_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.1493374803135605 + 0.16238989364956619*as + 0.1493374803135605*epsb + 2.1627510283726084*chi2one*epsb + 1.67508463989745*chi2pone*epsb + -2.1627510283726084*chi3pone*epsb + -0.837542319948725*chi3ppone*epsb + 0.1493374803135605*epsc + -2.1627510283726084*chi3pone*epsc + -0.837542319948725*chi3ppone*epsc + 0.298674960627121*epsb*etaone + 1.0813755141863042*epsb*etapone + 0.4187711599743625*epsb*etappone + 0.5406877570931521*xipone + 0.5744055734710128*as*xipone + 0.5406877570931521*epsb*xipone + 1.67508463989745*chi2one*epsb*xipone + -1.67508463989745*chi3pone*epsb*xipone + 0.5406877570931521*epsc*xipone + -1.67508463989745*chi3pone*epsc*xipone + 1.0813755141863042*epsb*etaone*xipone + 0.837542319948725*epsb*etapone*xipone + 0.20938557998718124*xippone + 0.1909062556224326*as*xippone + 0.20938557998718124*epsb*xippone + 0.20938557998718124*epsc*xippone + 0.4187711599743625*epsb*etaone*xippone + 0.1493374803135605*l2one*pow(epsc,2) + 0.4359949670995614*l2pone*pow(epsc,2) + -0.1493374803135605*l5one*pow(epsc,2) + -0.4359949670995614*l5pone*pow(epsc,2) + 0.5406877570931521*l2one*xipone*pow(epsc,2) + 0.4187711599743625*l2pone*xipone*pow(epsc,2) + -0.5406877570931521*l5one*xipone*pow(epsc,2) + -0.4187711599743625*l5pone*xipone*pow(epsc,2) + 0.20938557998718124*l2one*xippone*pow(epsc,2) + -0.20938557998718124*l5one*xippone*pow(epsc,2);
        }

        double A3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.006711208710683109 + -0.0028293672950845494*as + 0.006711208710683109*l2one*pow(epsc,2);
        }

        double A3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.06881919420215554 + 0.026844834842732437*epsb + 0.2147586787418595*chi2one*epsb + -0.2147586787418595*chi3pone*epsb + 0.026844834842732437*epsc + -0.2147586787418595*chi3pone*epsc + 0.053689669685464875*epsb*etaone + as*(0.005289452358091395 - 0.022634938360676395*xipone) + 0.053689669685464875*xipone + (0.053689669685464875*l2pone - 0.026844834842732437*l5one + l2one*(0.06881919420215554 + 0.053689669685464875*xipone))*pow(epsc,2);
        }

        double A3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.28459648348892946 + 0.2215871071231573*epsb + 2.6317315719526966*chi2one*epsb + 1.718069429934876*chi2pone*epsb + -2.6317315719526966*chi3pone*epsb + -0.859034714967438*chi3ppone*epsb + 0.2215871071231573*epsc + -2.6317315719526966*chi3pone*epsc + -0.859034714967438*chi3ppone*epsc + 0.4431742142463146*epsb*etaone + 0.429517357483719*epsb*etapone + 0.6579328929881741*xipone + 0.2147586787418595*epsb*xipone + 1.718069429934876*chi2one*epsb*xipone + -1.718069429934876*chi3pone*epsb*xipone + 0.2147586787418595*epsc*xipone + -1.718069429934876*chi3pone*epsc*xipone + 0.429517357483719*epsb*etaone*xipone + as*(0.14308066076970227 - 0.002954257856621647*xipone - 0.09053975344270558*xippone) + 0.2147586787418595*xippone + (0.5505535536172445*l2pone - 0.2215871071231573*l5one - 0.2147586787418595*l5pone + 0.429517357483719*l2pone*xipone - 0.2147586787418595*l5one*xipone + l2one*(0.28459648348892946 + 0.6579328929881741*xipone + 0.2147586787418595*xippone))*pow(epsc,2);
        }

        double A4_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.006711208710683109 + -0.0028293672950845494*as + 0.006711208710683109*l2one*pow(epsc,2);
        }

        double A4_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.06881919420215554 + 0.026844834842732437*epsb + -0.2147586787418595*chi3pone*epsb + 0.026844834842732437*epsc + 0.2147586787418595*chi2one*epsc + -0.2147586787418595*chi3pone*epsc + 0.053689669685464875*epsc*etaone + as*(0.005289452358091395 - 0.022634938360676395*xipone) + 0.053689669685464875*xipone + (0.053689669685464875*l2pone + 0.053689669685464875*l3one + 0.026844834842732437*l5one - 0.053689669685464875*l6one + l2one*(0.06881919420215554 + 0.053689669685464875*xipone))*pow(epsc,2);
        }

        double A4_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.28459648348892946 + 0.2215871071231573*epsb + -2.6317315719526966*chi3pone*epsb + -0.859034714967438*chi3ppone*epsb + 0.2215871071231573*epsc + 2.6317315719526966*chi2one*epsc + 1.718069429934876*chi2pone*epsc + -2.6317315719526966*chi3pone*epsc + -0.859034714967438*chi3ppone*epsc + 0.4431742142463146*epsc*etaone + 0.429517357483719*epsc*etapone + 0.6579328929881741*xipone + 0.2147586787418595*epsb*xipone + -1.718069429934876*chi3pone*epsb*xipone + 0.2147586787418595*epsc*xipone + 1.718069429934876*chi2one*epsc*xipone + -1.718069429934876*chi3pone*epsc*xipone + 0.429517357483719*epsc*etaone*xipone + as*(0.14308066076970227 - 0.002954257856621647*xipone - 0.09053975344270558*xippone) + 0.2147586787418595*xippone + (0.5505535536172445*l2pone + 0.6579328929881741*l3one + 0.429517357483719*l3pone + 0.2215871071231573*l5one + 0.2147586787418595*l5pone - 0.6579328929881741*l6one - 0.429517357483719*l6pone + 0.429517357483719*l2pone*xipone + 0.429517357483719*l3one*xipone + 0.2147586787418595*l5one*xipone - 0.429517357483719*l6one*xipone + l2one*(0.28459648348892946 + 0.6579328929881741*xipone + 0.2147586787418595*xippone))*pow(epsc,2);
        }

        double A7_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();

            return 0.006412276517110343 + -0.0027033409698122754*as + 0.006412276517110343*l2one*pow(epsc,2);
        }

        double A7_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.053295141150212194 + -0.05675299940030692*epsb + -0.205192848547531*chi3pone*epsb + 0.05675299940030692*epsc + -0.205192848547531*chi3pone*epsc + as*(0.01306681716595138 - 0.021626727758498204*xipone) + 0.05129821213688275*xipone + (0.05129821213688275*l2pone - 0.05675299940030692*l5one + l2one*(0.053295141150212194 + 0.05129821213688275*xipone))*pow(epsc,2);
        }

        double A7_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();

            return 0.15627762643821283 + -0.3581921735514904*epsb + -2.1158302139018526*chi3pone*epsb + -0.820771394190124*chi3ppone*epsb + 0.3581921735514904*epsc + -2.1158302139018526*chi3pone*epsc + -0.820771394190124*chi3ppone*epsc + 0.5289575534754631*xipone + -0.4540239952024554*epsb*xipone + -1.641542788380248*chi3pone*epsb*xipone + 0.4540239952024554*epsc*xipone + -1.641542788380248*chi3pone*epsc*xipone + as*(0.13798340585773658 + 0.061281081810614674*xipone - 0.08650691103399281*xippone) + 0.205192848547531*xippone + (0.4263611292016976*l2pone - 0.3581921735514904*l5one - 0.4540239952024554*l5pone + 0.410385697095062*l2pone*xipone - 0.4540239952024554*l5one*xipone + l2one*(0.15627762643821283 + 0.5289575534754631*xipone + 0.205192848547531*xippone))*pow(epsc,2);
        }
        // }}}

        // B_s -> D_s form factors
        // {{{
        double V1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004856223859143635 + 0.0017747026012450353*as + -0.002249995326965876*epsb + 0.002249995326965876*epsc + 0.004499990653931752*epsb*etasone + -0.004499990653931752*epsc*etasone + (0.004856223859143635*l1sone - 0.002249995326965876*l4sone)*pow(epsc,2);
        }

        double V1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04087217915418119 + -0.018936979588918035*epsb + -0.15539916349259633*chi2sone*epsb + 0.466197490477789*chi3spone*epsb + 0.018936979588918035*epsc + -0.15539916349259633*chi2sone*epsc + 0.466197490477789*chi3spone*epsc + 0.03787395917783607*epsb*etasone + -0.03787395917783607*epsc*etasone + 0.03599992523145402*epsb*etaspone + -0.03599992523145402*epsc*etaspone + as*(0.01886199951407103 + 0.014197620809960282*xispone) + 0.03884979087314908*xispone + -0.01799996261572701*epsb*xispone + 0.01799996261572701*epsc*xispone + 0.03599992523145402*epsb*etasone*xispone + -0.03599992523145402*epsc*etasone*xispone + (0.03884979087314908*l1spone - 0.018936979588918035*l4sone - 0.01799996261572701*l4spone + l1sone*(0.04087217915418119 + 0.03884979087314908*xispone) - 0.01799996261572701*l4sone*xispone)*pow(epsc,2);
        }

        double V1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();


            return 0.11343397437390446 + -0.052556455316592324*epsb + -1.618708059918991*chi2sone*epsb + -1.2431933079407707*chi2spone*epsb + 4.856124179756971*chi3spone*epsb + 1.8647899619111559*chi3sppone*epsb + 0.052556455316592324*epsc + -1.618708059918991*chi2sone*epsc + -1.2431933079407707*chi2spone*epsc + 4.856124179756971*chi3spone*epsc + 1.8647899619111559*chi3sppone*epsc + 0.10511291063318465*epsb*etasone + -0.10511291063318465*epsc*etasone + 0.3749915238855965*epsb*etaspone + -0.3749915238855965*epsc*etaspone + 0.14399970092581607*epsb*etasppone + -0.14399970092581607*epsc*etasppone + 0.40467701497974773*xispone + -0.18749576194279824*epsb*xispone + -1.2431933079407707*chi2sone*epsb*xispone + 3.7295799238223117*chi3spone*epsb*xispone + 0.18749576194279824*epsc*xispone + -1.2431933079407707*chi2sone*epsc*xispone + 3.7295799238223117*chi3spone*epsc*xispone + 0.3749915238855965*epsb*etasone*xispone + -0.3749915238855965*epsc*etasone*xispone + 0.28799940185163214*epsb*etaspone*xispone + -0.28799940185163214*epsc*etaspone*xispone + as*(0.029490534535912005 + 0.17929123773248878*xispone + 0.05679048323984112*xisppone) + 0.15539916349259633*xisppone + -0.07199985046290804*epsb*xisppone + 0.07199985046290804*epsc*xisppone + 0.14399970092581607*epsb*etasone*xisppone + -0.14399970092581607*epsc*etasone*xisppone + (0.3269774332334495*l1spone - 0.052556455316592324*l4sone - 0.15149583671134426*l4spone + 0.31079832698519266*l1spone*xispone - 0.18749576194279824*l4sone*xispone - 0.14399970092581607*l4spone*xispone + l1sone*(0.11343397437390446 + 0.40467701497974773*xispone + 0.15539916349259633*xisppone) - 0.07199985046290804*l4sone*xisppone)*pow(epsc,2);
        }

        double S1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.027578719782544572 + 0.006758946679380605*as + 0.027578719782544572*l1sone*pow(epsc,2);
        }

        double S1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();


            return 0.1992219965165159 + 0.1358421677238942*as + -0.2380954936350598*epsb + -0.8825190330414263*chi2sone*epsb + 2.6475570991242794*chi3spone*epsb + 0.2380954936350598*epsc + -0.8825190330414263*chi2sone*epsc + 2.6475570991242794*chi3spone*epsc + 0.4761909872701196*epsb*etasone + -0.4761909872701196*epsc*etasone + 0.22062975826035658*xispone + 0.054071573435044845*as*xispone + 0.1992219965165159*l1sone*pow(epsc,2) + 0.22062975826035658*l1spone*pow(epsc,2) + -0.2380954936350598*l4sone*pow(epsc,2) + 0.22062975826035658*l1sone*xispone*pow(epsc,2);
        }

        double S1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();


            return 0.5502882657016179 + 0.3404075815206873*as + -1.2437532297774159*epsb + -8.14014195461136*chi2sone*epsb + -7.0601522643314105*chi2spone*epsb + 24.42042586383408*chi3spone*epsb + 10.590228396497118*chi3sppone*epsb + 1.2437532297774159*epsc + -8.14014195461136*chi2sone*epsc + -7.0601522643314105*chi2spone*epsc + 24.42042586383408*chi3spone*epsc + 10.590228396497118*chi3sppone*epsc + 2.4875064595548317*epsb*etasone + -2.4875064595548317*epsc*etasone + 3.8095278981609573*epsb*etaspone + -3.8095278981609573*epsc*etaspone + 2.03503548865284*xispone + 1.1948804886612434*as*xispone + -1.9047639490804786*epsb*xispone + -7.0601522643314105*chi2sone*epsb*xispone + 21.180456792994235*chi3spone*epsb*xispone + 1.9047639490804786*epsc*xispone + -7.0601522643314105*chi2sone*epsc*xispone + 21.180456792994235*chi3spone*epsc*xispone + 3.8095278981609573*epsb*etasone*xispone + -3.8095278981609573*epsc*etasone*xispone + 0.8825190330414263*xisppone + 0.21628629374017938*as*xisppone + 0.5502882657016179*l1sone*pow(epsc,2) + 1.593775972132127*l1spone*pow(epsc,2) + -1.2437532297774159*l4sone*pow(epsc,2) + -1.9047639490804786*l4spone*pow(epsc,2) + 2.03503548865284*l1sone*xispone*pow(epsc,2) + 1.7650380660828526*l1spone*xispone*pow(epsc,2) + -1.9047639490804786*l4sone*xispone*pow(epsc,2) + 0.8825190330414263*l1sone*xisppone*pow(epsc,2);
        }
        // }}}

        // B_s -> D_s^* form factors
        // {{{
        double A1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0039057998402766254 + -0.0016466396450514063*as + 0.0039057998402766254*l2sone*pow(epsc,2);
        }

        double A1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04541491752041236 + 0.0008172109651720427*as + 0.015623199361106501*epsb + -0.12498559488885201*chi2sone*epsb + 0.3749567846665561*chi3spone*epsb + 0.015623199361106501*epsc + -0.12498559488885201*chi3spone*epsc + -0.031246398722213003*epsb*etasone + 0.031246398722213003*xispone + -0.01317311716041125*as*xispone + 0.04541491752041236*l2sone*pow(epsc,2) + 0.031246398722213003*l2spone*pow(epsc,2) + -0.015623199361106501*l5sone*pow(epsc,2) + 0.031246398722213003*l2sone*xispone*pow(epsc,2);
        }

        double A1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.2138412950270066 + 0.09035882754434083*as + 0.15041327135943644*epsb + -1.7032485504308996*chi2sone*epsb + -0.9998847591108161*chi2spone*epsb + 5.109745651292698*chi3spone*epsb + 1.4998271386662243*chi3sppone*epsb + 0.15041327135943644*epsc + -1.7032485504308996*chi3spone*epsc + -0.49994237955540805*chi3sppone*epsc + -0.3008265427188729*epsb*etasone + -0.24997118977770402*epsb*etaspone + 0.4258121376077249*xispone + -0.019808546599446152*as*xispone + 0.12498559488885201*epsb*xispone + -0.9998847591108161*chi2sone*epsb*xispone + 2.9996542773324486*chi3spone*epsb*xispone + 0.12498559488885201*epsc*xispone + -0.9998847591108161*chi3spone*epsc*xispone + -0.24997118977770402*epsb*etasone*xispone + 0.12498559488885201*xisppone + -0.052692468641645*as*xisppone + 0.2138412950270066*l2sone*pow(epsc,2) + 0.36331934016329887*l2spone*pow(epsc,2) + -0.15041327135943644*l5sone*pow(epsc,2) + -0.12498559488885201*l5spone*pow(epsc,2) + 0.4258121376077249*l2sone*xispone*pow(epsc,2) + 0.24997118977770402*l2spone*xispone*pow(epsc,2) + -0.12498559488885201*l5sone*xispone*pow(epsc,2) + 0.12498559488885201*l2sone*xisppone*pow(epsc,2);
        }

        double A5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0025298459520406175 + -0.0010665535385468866*as + 0.0025298459520406175*l2sone*pow(epsc,2);
        }

        double A5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.024488910660312665 + -0.023253219521903006*epsb + -0.08095507046529976*chi2sone*epsb + 0.24286521139589928*chi3spone*epsb + 0.023253219521903006*epsc + 0.08095507046529976*chi2sone*epsc + -0.08095507046529976*chi3spone*epsc + 0.04650643904380601*epsb*etasone + 0.04650643904380601*epsc*etasone + as*(0.004367737341141352 - 0.008532428308375093*xispone) + 0.02023876761632494*xispone + (0.02023876761632494*l2spone + 0.02023876761632494*l3sone + 0.023253219521903006*l5sone - 0.04650643904380601*l6sone + l2sone*(0.024488910660312665 + 0.02023876761632494*xispone))*pow(epsc,2);
        }

        double A5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.08602032953147207 + -0.178584742882586*epsb + -0.9455552820606049*chi2sone*epsb + -0.6476405637223982*chi2spone*epsb + 2.8366658461818144*chi3spone*epsb + 0.971460845583597*chi3sppone*epsb + 0.17858474288258597*epsc + 0.9455552820606049*chi2sone*epsc + 0.6476405637223982*chi2spone*epsc + -0.9455552820606049*chi3spone*epsc + -0.3238202818611991*chi3sppone*epsc + 0.357169485765172*epsb*etasone + 0.357169485765172*epsc*etasone + 0.37205151235044814*epsb*etaspone + 0.37205151235044814*epsc*etaspone + 0.23638882051515117*xispone + -0.18602575617522407*epsb*xispone + -0.6476405637223982*chi2sone*epsb*xispone + 1.942921691167194*chi3spone*epsb*xispone + 0.18602575617522407*epsc*xispone + 0.6476405637223982*chi2sone*epsc*xispone + -0.6476405637223982*chi3spone*epsc*xispone + 0.37205151235044814*epsb*etasone*xispone + 0.37205151235044814*epsc*etasone*xispone + as*(0.069484337086071 + 0.017877042112380646*xispone - 0.03412971323350037*xisppone) + 0.08095507046529977*xisppone + (0.1959112852825013*l2spone + 0.23638882051515123*l3sone + 0.16191014093059955*l3spone + 0.178584742882586*l5sone + 0.18602575617522407*l5spone - 0.543195241940396*l6sone - 0.37205151235044814*l6spone + 0.16191014093059955*l2spone*xispone + 0.16191014093059955*l3sone*xispone + 0.18602575617522407*l5sone*xispone - 0.37205151235044814*l6sone*xispone + l2sone*(0.08602032953147207 + 0.23638882051515117*xispone + 0.08095507046529977*xisppone))*pow(epsc,2);
        }

        double V4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004686857964851364 + 0.004273219315096704*as + 0.004686857964851364*epsb + 0.004686857964851364*epsc + -0.009373715929702727*epsb*etasone + 0.004686857964851364*l2sone*pow(epsc,2) + -0.004686857964851364*l5sone*pow(epsc,2);
        }

        double V4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04346831296433728 + 0.04692344500022692*as + 0.04346831296433728*epsb + -0.14997945487524364*chi2sone*epsb + 0.4499383646257309*chi3spone*epsb + 0.04346831296433728*epsc + -0.14997945487524364*chi3spone*epsc + -0.08693662592867456*epsb*etasone + -0.07498972743762182*epsb*etaspone + 0.03749486371881091*xispone + 0.03418575452077363*as*xispone + 0.03749486371881091*epsb*xispone + 0.03749486371881091*epsc*xispone + -0.07498972743762182*epsb*etasone*xispone + 0.04346831296433728*l2sone*pow(epsc,2) + 0.03749486371881091*l2spone*pow(epsc,2) + -0.04346831296433728*l5sone*pow(epsc,2) + -0.03749486371881091*l5spone*pow(epsc,2) + 0.03749486371881091*l2sone*xispone*pow(epsc,2) + -0.03749486371881091*l5sone*xispone*pow(epsc,2);
        }

        double V4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.134960773305362 + 0.148733363687041*as + 0.134960773305362*epsb + -1.6909449246092803*chi2sone*epsb + -1.1998356390019491*chi2spone*epsb + 5.07283477382784*chi3spone*epsb + 1.7997534585029236*chi3sppone*epsb + 0.134960773305362*epsc + -1.6909449246092803*chi3spone*epsc + -0.5999178195009746*chi3sppone*epsc + -0.269921546610724*epsb*etasone + -0.8454724623046401*epsb*etaspone + -0.2999589097504873*epsb*etasppone + 0.42273623115232*xispone + 0.44375906904336254*as*xispone + 0.42273623115232*epsb*xispone + -1.1998356390019491*chi2sone*epsb*xispone + 3.599506917005847*chi3spone*epsb*xispone + 0.42273623115232*epsc*xispone + -1.1998356390019491*chi3spone*epsc*xispone + -0.84547246230464*epsb*etasone*xispone + -0.5999178195009746*epsb*etaspone*xispone + 0.14997945487524364*xisppone + 0.13674301808309453*as*xisppone + 0.14997945487524364*epsb*xisppone + 0.14997945487524364*epsc*xisppone + -0.2999589097504873*epsb*etasone*xisppone + 0.134960773305362*l2sone*pow(epsc,2) + 0.34774650371469823*l2spone*pow(epsc,2) + -0.134960773305362*l5sone*pow(epsc,2) + -0.34774650371469823*l5spone*pow(epsc,2) + 0.42273623115232*l2sone*xispone*pow(epsc,2) + 0.2999589097504873*l2spone*xispone*pow(epsc,2) + -0.42273623115232*l5sone*xispone*pow(epsc,2) + -0.2999589097504873*l5spone*xispone*pow(epsc,2) + 0.14997945487524364*l2sone*xisppone*pow(epsc,2) + -0.14997945487524364*l5sone*xisppone*pow(epsc,2);
        }

        double P1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.017623813657930924 + -0.0014523480573094521*as + -0.007669567407650894*epsb + 0.007669567407650894*epsc + 0.015339134815301789*epsb*etasone + 0.015339134815301789*epsc*etasone + 0.017623813657930924*l2sone*pow(epsc,2) + 0.007669567407650894*l5sone*pow(epsc,2) + -0.015339134815301789*l6sone*pow(epsc,2);
        }

        double P1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.17391618332456305 + -0.02962974998746001*as + -0.07568520169236201*epsb + -0.5639620370537896*chi2sone*epsb + 1.6918861111613686*chi3spone*epsb + 0.07568520169236201*epsc + 0.5639620370537896*chi2sone*epsc + -0.5639620370537896*chi3spone*epsc + 0.15137040338472402*epsb*etasone + 0.15137040338472402*epsc*etasone + 0.12271307852241431*epsb*etaspone + 0.12271307852241431*epsc*etaspone + 0.1409905092634474*xispone + -0.011618784458475612*as*xispone + -0.061356539261207155*epsb*xispone + 0.061356539261207155*epsc*xispone + 0.12271307852241431*epsb*etasone*xispone + 0.12271307852241431*epsc*etasone*xispone + 0.17391618332456305*l2sone*pow(epsc,2) + 0.1409905092634474*l2spone*pow(epsc,2) + 0.1409905092634474*l3sone*pow(epsc,2) + 0.07568520169236201*l5sone*pow(epsc,2) + 0.061356539261207155*l5spone*pow(epsc,2) + -0.21272694264593117*l6sone*pow(epsc,2) + -0.12271307852241431*l6spone*pow(epsc,2) + 0.1409905092634474*l2sone*xispone*pow(epsc,2) + 0.061356539261207155*l5sone*xispone*pow(epsc,2) + -0.12271307852241431*l6sone*xispone*pow(epsc,2);
        }

        double P1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.6200566729823669 + -0.38104285222530004*as + -0.2698375358651133*epsb + -6.6932419404935946*chi2sone*epsb + -4.511696296430316*chi2spone*epsb + 20.079725821480782*chi3spone*epsb + 6.767544444645476*chi3sppone*epsb + 0.2698375358651133*epsc + 6.6932419404935946*chi2sone*epsc + 4.511696296430316*chi2spone*epsc + -6.6932419404935946*chi3spone*epsc + -2.255848148215158*chi3sppone*epsc + 0.5396750717302266*epsb*etasone + 0.5396750717302266*epsc*etasone + 1.4563893841226205*epsb*etaspone + 1.456389384122621*epsc*etaspone + 0.49085231408965724*epsb*etasppone + 0.49085231408965724*epsc*etasppone + 1.6733104851233989*xispone + -0.26027556881663116*as*xispone + -0.7281946920613103*epsb*xispone + -4.511696296430316*chi2sone*epsb*xispone + 13.535088889290952*chi3spone*epsb*xispone + 0.7281946920613105*epsc*xispone + 4.511696296430316*chi2sone*epsc*xispone + -4.511696296430316*chi3spone*epsc*xispone + 1.4563893841226205*epsb*etasone*xispone + 1.4563893841226205*epsc*etasone*xispone + 0.9817046281793145*epsb*etaspone*xispone + 0.9817046281793145*epsc*etaspone*xispone + 0.5639620370537896*xisppone + -0.04647513783390245*as*xisppone + -0.24542615704482862*epsb*xisppone + 0.24542615704482862*epsc*xisppone + 0.49085231408965724*epsb*etasone*xisppone + 0.49085231408965724*epsc*etasone*xisppone + 0.6200566729823669*l2sone*pow(epsc,2) + 1.391329466596504*l2spone*pow(epsc,2) + 1.6733104851233986*l3sone*pow(epsc,2) + 1.127924074107579*l3spone*pow(epsc,2) + 0.2698375358651133*l5sone*pow(epsc,2) + 0.6054816135388961*l5spone*pow(epsc,2) + -1.267869763791537*l6sone*pow(epsc,2) + -1.7018155411674496*l6spone*pow(epsc,2) + 1.6733104851233989*l2sone*xispone*pow(epsc,2) + 1.127924074107579*l2spone*xispone*pow(epsc,2) + 1.127924074107579*l3sone*xispone*pow(epsc,2) + 0.7281946920613103*l5sone*xispone*pow(epsc,2) + 0.49085231408965724*l5spone*xispone*pow(epsc,2) + -1.947241698212278*l6sone*xispone*pow(epsc,2) + -0.9817046281793145*l6spone*xispone*pow(epsc,2) + 0.5639620370537896*l2sone*xisppone*pow(epsc,2) + 0.24542615704482862*l5sone*xisppone*pow(epsc,2) + -0.49085231408965724*l6sone*xisppone*pow(epsc,2);
        }
        // }}}

        // B_s^* -> D_s form factors
        // {{{
        double P2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03295002855689904 + -0.0019019140260784866*as + -0.015382939460954132*epsb + 0.015382939460954132*epsc + -0.030765878921908265*epsb*etasone + -0.030765878921908265*epsc*etasone + 0.03295002855689904*l1sone*pow(epsc,2) + -0.015382939460954132*l4sone*pow(epsc,2);
        }

        double P2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.2781321885136799 + -0.04348367290967678*as + -0.12984785766301843*epsb + 1.0544009138207693*chi2sone*epsb + -1.0544009138207693*chi3spone*epsb + 0.12984785766301843*epsc + -1.0544009138207693*chi2sone*epsc + 3.1632027414623076*chi3spone*epsc + -0.25969571532603686*epsb*etasone + -0.25969571532603686*epsc*etasone + -0.24612703137526612*epsb*etaspone + -0.24612703137526612*epsc*etaspone + 0.2636002284551923*xispone + -0.015215312208627893*as*xispone + -0.12306351568763306*epsb*xispone + 0.12306351568763306*epsc*xispone + -0.24612703137526612*epsb*etasone*xispone + -0.24612703137526612*epsc*etasone*xispone + 0.2781321885136799*l1sone*pow(epsc,2) + 0.2636002284551923*l1spone*pow(epsc,2) + -0.12984785766301843*l4sone*pow(epsc,2) + -0.12306351568763306*l4spone*pow(epsc,2) + 0.2636002284551923*l1sone*xispone*pow(epsc,2) + -0.12306351568763306*l4sone*xispone*pow(epsc,2);
        }

        double P2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.8002963903908009 + -0.6117683772791116*as + -0.3736236799595781*epsb + 11.009031860079295*chi2sone*epsb + 8.435207310566156*chi2spone*epsb + -11.009031860079295*chi3spone*epsb + -4.217603655283078*chi3sppone*epsb + 0.3736236799595781*epsc + -11.009031860079295*chi2sone*epsc + -8.435207310566156*chi2spone*epsc + 33.027095580237884*chi3spone*epsc + 12.65281096584923*chi3sppone*epsc + -0.7472473599191563*epsb*etasone + -0.7472473599191563*epsc*etasone + -2.569819785358827*epsb*etaspone + -2.569819785358827*epsc*etaspone + -0.9845081255010645*epsb*etasppone + -0.9845081255010645*epsc*etasppone + 2.752257965019824*xispone + -0.37830000769467*as*xispone + -1.2849098926794136*epsb*xispone + 8.435207310566156*chi2sone*epsb*xispone + -8.435207310566156*chi3spone*epsb*xispone + 1.2849098926794136*epsc*xispone + -8.435207310566156*chi2sone*epsc*xispone + 25.30562193169846*chi3spone*epsc*xispone + -2.569819785358827*epsb*etasone*xispone + -2.569819785358827*epsc*etasone*xispone + -1.969016251002129*epsb*etaspone*xispone + -1.969016251002129*epsc*etaspone*xispone + 1.0544009138207695*xisppone + -0.06086124883451157*as*xisppone + -0.49225406275053224*epsb*xisppone + 0.49225406275053224*epsc*xisppone + -0.9845081255010645*epsb*etasone*xisppone + -0.9845081255010645*epsc*etasone*xisppone + 0.8002963903908009*l1sone*pow(epsc,2) + 2.225057508109439*l1spone*pow(epsc,2) + -0.3736236799595781*l4sone*pow(epsc,2) + -1.0387828613041474*l4spone*pow(epsc,2) + 2.752257965019824*l1sone*xispone*pow(epsc,2) + 2.108801827641539*l1spone*xispone*pow(epsc,2) + -1.2849098926794136*l4sone*xispone*pow(epsc,2) + -0.9845081255010645*l4spone*xispone*pow(epsc,2) + 1.0544009138207695*l1sone*xisppone*pow(epsc,2) + -0.49225406275053224*l4sone*xisppone*pow(epsc,2);
        }

        double V5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.0052838463097047055 + 0.004817520453566517*as + 0.0052838463097047055*epsb + 0.0052838463097047055*epsc + -0.010567692619409411*epsc*etasone + 0.0052838463097047055*l1sone*pow(epsc,2) + -0.0052838463097047055*l4sone*pow(epsc,2);
        }

        double V5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04663432375816549 + 0.05073878485941332*as + 0.04663432375816549*epsb + -0.16908308191055058*chi3spone*epsb + 0.04663432375816549*epsc + -0.16908308191055058*chi2sone*epsc + 0.5072492457316516*chi3spone*epsc + -0.09326864751633097*epsc*etasone + -0.08454154095527529*epsc*etaspone + 0.042270770477637644*xispone + 0.03854016362853214*as*xispone + 0.042270770477637644*epsb*xispone + 0.042270770477637644*epsc*xispone + -0.08454154095527529*epsc*etasone*xispone + 0.04663432375816549*l1sone*pow(epsc,2) + 0.042270770477637644*l1spone*pow(epsc,2) + -0.04663432375816549*l4sone*pow(epsc,2) + -0.042270770477637644*l4spone*pow(epsc,2) + 0.042270770477637644*l1sone*xispone*pow(epsc,2) + -0.042270770477637644*l4sone*xispone*pow(epsc,2);
        }

        double V5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1368797742680169 + 0.15006620551181418*as + 0.1368797742680169*epsb + -1.8304645240823967*chi3spone*epsb + -0.6763323276422023*chi3sppone*epsb + 0.1368797742680169*epsc + -1.8304645240823967*chi2sone*epsc + -1.3526646552844046*chi2spone*epsc + 5.49139357224719*chi3spone*epsc + 2.028996982926607*chi3sppone*epsc + -0.2737595485360338*epsc*etasone + -0.9152322620411983*epsc*etaspone + -0.33816616382110115*epsc*etasppone + 0.4576161310205992*xispone + 0.4829906061323709*as*xispone + 0.4576161310205992*epsb*xispone + -1.3526646552844046*chi3spone*epsb*xispone + 0.4576161310205992*epsc*xispone + -1.3526646552844046*chi2sone*epsc*xispone + 4.057993965853213*chi3spone*epsc*xispone + -0.9152322620411985*epsc*etasone*xispone + -0.6763323276422023*epsc*etaspone*xispone + 0.16908308191055058*xisppone + 0.15416065451412855*as*xisppone + 0.16908308191055058*epsb*xisppone + 0.16908308191055058*epsc*xisppone + -0.33816616382110115*epsc*etasone*xisppone + 0.1368797742680169*l1sone*pow(epsc,2) + 0.3730745900653239*l1spone*pow(epsc,2) + -0.1368797742680169*l4sone*pow(epsc,2) + -0.3730745900653239*l4spone*pow(epsc,2) + 0.4576161310205992*l1sone*xispone*pow(epsc,2) + 0.33816616382110115*l1spone*xispone*pow(epsc,2) + -0.4576161310205992*l4sone*xispone*pow(epsc,2) + -0.33816616382110115*l4spone*xispone*pow(epsc,2) + 0.16908308191055058*l1sone*xisppone*pow(epsc,2) + -0.16908308191055058*l4sone*xisppone*pow(epsc,2);
        }

        double A2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.005036604965156558 + -0.0021233738929904216*as + 0.005036604965156558*l1sone*pow(epsc,2);
        }

        double A2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.054823780290389645 + 0.002630398208386344*as + 0.020146419860626232*epsb + -0.16117135888500986*chi3spone*epsb + 0.020146419860626232*epsc + -0.16117135888500986*chi2sone*epsc + 0.4835140766550296*chi3spone*epsc + -0.040292839721252464*epsc*etasone + 0.040292839721252464*xispone + -0.01698699114392337*as*xispone + 0.054823780290389645*l1sone*pow(epsc,2) + 0.040292839721252464*l1spone*pow(epsc,2) + -0.020146419860626232*l4sone*pow(epsc,2) + 0.040292839721252464*l1sone*xispone*pow(epsc,2);
        }

        double A2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.2413414334718427 + 0.11191245450655582*as + 0.1790022814403061*epsb + -2.076703687062488*chi3spone*epsb + -0.6446854355400394*chi3sppone*epsb + 0.1790022814403061*epsc + -2.076703687062488*chi2sone*epsc + -1.2893708710800789*chi2spone*epsc + 6.230111061187466*chi3spone*epsc + 1.9340563066201184*chi3sppone*epsc + -0.3580045628806122*epsc*etasone + -0.3223427177700197*epsc*etaspone + 0.519175921765622*xispone + -0.012930796620755983*as*xispone + 0.16117135888500986*epsb*xispone + -1.2893708710800789*chi3spone*epsb*xispone + 0.16117135888500986*epsc*xispone + -1.2893708710800789*chi2sone*epsc*xispone + 3.868112613240237*chi3spone*epsc*xispone + -0.3223427177700197*epsc*etasone*xispone + 0.16117135888500986*xisppone + -0.06794796457569348*as*xisppone + 0.2413414334718427*l1sone*pow(epsc,2) + 0.43859024232311716*l1spone*pow(epsc,2) + -0.1790022814403061*l4sone*pow(epsc,2) + -0.16117135888500986*l4spone*pow(epsc,2) + 0.519175921765622*l1sone*xispone*pow(epsc,2) + 0.3223427177700197*l1spone*xispone*pow(epsc,2) + -0.16117135888500986*l4sone*xispone*pow(epsc,2) + 0.16117135888500986*l1sone*xisppone*pow(epsc,2);
        }

        double A6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0035294640354068636 + -0.001487980860297349*as + 0.0035294640354068636*l1sone*pow(epsc,2);
        }

        double A6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03157617379125103 + 0.0054695955072789635*as + -0.030240238818437362*epsb + 0.11294284913301962*chi2sone*epsb + -0.11294284913301962*chi3spone*epsb + 0.030240238818437362*epsc + -0.11294284913301962*chi2sone*epsc + 0.3388285473990589*chi3spone*epsc + -0.060480477636874724*epsb*etasone + -0.060480477636874724*epsc*etasone + 0.028235712283254905*xispone + -0.01190384688237879*as*xispone + 0.03157617379125103*l1sone*pow(epsc,2) + 0.028235712283254905*l1spone*pow(epsc,2) + -0.030240238818437362*l4sone*pow(epsc,2) + 0.028235712283254905*l1sone*xispone*pow(epsc,2);
        }

        double A6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1012833161041574 + 0.07541197778656607*as + -0.2100623092681177*epsb + 1.2363232595860718*chi2sone*epsb + 0.9035427930641569*chi2spone*epsb + -1.2363232595860718*chi3spone*epsb + -0.45177139653207843*chi3sppone*epsb + 0.2100623092681177*epsc + -1.2363232595860718*chi2sone*epsc + -0.9035427930641569*chi2spone*epsc + 3.708969778758216*chi3spone*epsc + 1.3553141895962353*chi3sppone*epsc + -0.4201246185362354*epsb*etasone + -0.4201246185362354*epsc*etasone + -0.4838438210949978*epsb*etaspone + -0.4838438210949978*epsc*etaspone + 0.30908081489651795*xispone + 0.019949070293474133*as*xispone + -0.2419219105474989*epsb*xispone + 0.9035427930641569*chi2sone*epsb*xispone + -0.9035427930641569*chi3spone*epsb*xispone + 0.2419219105474989*epsc*xispone + -0.9035427930641569*chi2sone*epsc*xispone + 2.7106283791924706*chi3spone*epsc*xispone + -0.4838438210949978*epsb*etasone*xispone + -0.4838438210949978*epsc*etasone*xispone + 0.11294284913301961*xisppone + -0.047615387529515156*as*xisppone + 0.1012833161041574*l1sone*pow(epsc,2) + 0.2526093903300082*l1spone*pow(epsc,2) + -0.2100623092681177*l4sone*pow(epsc,2) + -0.2419219105474989*l4spone*pow(epsc,2) + 0.30908081489651795*l1sone*xispone*pow(epsc,2) + 0.22588569826603921*l1spone*xispone*pow(epsc,2) + -0.2419219105474989*l4sone*xispone*pow(epsc,2) + 0.11294284913301961*l1sone*xisppone*pow(epsc,2);
        }
        // }}}

        // B_s^* -> D_s^* form factors
        // {{{
        double S2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.02706476794041704 + 0.006632988218505706*as + 0.02706476794041704*l2sone*pow(epsc,2);
        }

        double S2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.21900165467601973 + 0.1424559890126617*as + -0.2466988050391903*epsb + -0.8660725740933454*chi3spone*epsb + 0.2466988050391903*epsc + -0.8660725740933454*chi3spone*epsc + 0.21651814352333634*xispone + 0.05306390574804565*as*xispone + 0.21900165467601973*l2sone*pow(epsc,2) + 0.21651814352333634*l2spone*pow(epsc,2) + -0.2466988050391903*l5sone*pow(epsc,2) + 0.21651814352333634*l2sone*xispone*pow(epsc,2);
        }

        double S2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.6697744863595609 + 0.4576953641177782*as + -1.5028303505352347*epsb + -8.740198097819322*chi3spone*epsb + -3.4642902963733815*chi3sppone*epsb + 1.5028303505352347*epsc + -8.740198097819322*chi3spone*epsc + -3.4642902963733815*chi3sppone*epsc + 2.1850495244548305*xispone + 1.245775723597385*as*xispone + -1.9735904403135225*epsb*xispone + -6.928580592746763*chi3spone*epsb*xispone + 1.9735904403135225*epsc*xispone + -6.928580592746763*chi3spone*epsc*xispone + 0.8660725740933454*xisppone + 0.2122556229921826*as*xisppone + 0.6697744863595609*l2sone*pow(epsc,2) + 1.7520132374081578*l2spone*pow(epsc,2) + -1.5028303505352347*l5sone*pow(epsc,2) + -1.9735904403135225*l5spone*pow(epsc,2) + 2.1850495244548305*l2sone*xispone*pow(epsc,2) + 1.7321451481866907*l2spone*xispone*pow(epsc,2) + -1.9735904403135225*l5sone*xispone*pow(epsc,2) + 0.8660725740933454*l2sone*xisppone*pow(epsc,2);
        }

        double S3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.019137680941909156 + 0.0046902309488358615*as + 0.019137680941909156*l2sone*pow(epsc,2);
        }

        double S3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.15485755511248808 + 0.10073159585148939*as + -0.17444239795382943*epsb + 0.6124057901410931*chi2sone*epsb + -0.612405790141093*chi3spone*epsb + 0.17444239795382943*epsc + 0.6124057901410931*chi2sone*epsc + -0.612405790141093*chi3spone*epsc + -0.34888479590765886*epsb*etasone + 0.34888479590765886*epsc*etasone + 0.15310144753527324*xispone + 0.03752184759068689*as*xispone + 0.15485755511248808*l2sone*pow(epsc,2) + 0.15310144753527324*l2spone*pow(epsc,2) + 0.15310144753527327*l3sone*pow(epsc,2) + 0.17444239795382943*l5sone*pow(epsc,2) + -0.34888479590765886*l6sone*pow(epsc,2) + 0.15310144753527324*l2sone*xispone*pow(epsc,2);
        }

        double S3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.473602081170582 + 0.3236394956853269*as + -1.0626615318364205*epsb + 6.180253343881806*chi2sone*epsb + 4.8992463211287465*chi2spone*epsb + -6.1802533438818035*chi3spone*epsb + -2.449623160564373*chi3sppone*epsb + 1.0626615318364205*epsc + 6.180253343881806*chi2sone*epsc + 4.8992463211287465*chi2spone*epsc + -6.1802533438818035*chi3spone*epsc + -2.449623160564373*chi3sppone*epsc + -2.125323063672841*epsb*etasone + 2.125323063672841*epsc*etasone + -2.791078367261272*epsb*etaspone + 2.791078367261272*epsc*etaspone + 1.5450633359704509*xispone + 0.880896461993289*as*xispone + -1.395539183630636*epsb*xispone + 4.8992463211287465*chi2sone*epsb*xispone + -4.899246321128746*chi3spone*epsb*xispone + 1.395539183630636*epsc*xispone + 4.8992463211287465*chi2sone*epsc*xispone + -4.899246321128746*chi3spone*epsc*xispone + -2.791078367261272*epsb*etasone*xispone + 2.791078367261272*epsc*etasone*xispone + 0.6124057901410932*xisppone + 0.1500873903627476*as*xisppone + 0.473602081170582*l2sone*pow(epsc,2) + 1.2388604408999044*l2spone*pow(epsc,2) + 1.5450633359704515*l3sone*pow(epsc,2) + 1.2248115802821866*l3spone*pow(epsc,2) + 1.0626615318364205*l5sone*pow(epsc,2) + 1.395539183630636*l5spone*pow(epsc,2) + -3.5208622473034765*l6sone*pow(epsc,2) + -2.791078367261272*l6spone*pow(epsc,2) + 1.5450633359704509*l2sone*xispone*pow(epsc,2) + 1.2248115802821864*l2spone*xispone*pow(epsc,2) + 1.2248115802821866*l3sone*xispone*pow(epsc,2) + 1.395539183630636*l5sone*xispone*pow(epsc,2) + -2.791078367261272*l6sone*xispone*pow(epsc,2) + 0.6124057901410932*l2sone*xisppone*pow(epsc,2);
        }

        double P3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.021814694186003132 + -0.0017356714764909195*as + -0.009572963042792976*epsb + 0.009572963042792976*epsc + 0.021814694186003132*l2sone*pow(epsc,2) + -0.009572963042792976*l5sone*pow(epsc,2);
        }

        double P3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.22410145124848382 + -0.03667648192922286*as + -0.09834265345853424*epsb + -0.6980702139521002*chi3spone*epsb + 0.09834265345853424*epsc + -0.6980702139521002*chi3spone*epsc + 0.17451755348802506*xispone + -0.013885371811927357*as*xispone + -0.07658370434234381*epsb*xispone + 0.07658370434234381*epsc*xispone + 0.22410145124848382*l2sone*pow(epsc,2) + 0.17451755348802506*l2spone*pow(epsc,2) + -0.09834265345853424*l5sone*pow(epsc,2) + -0.07658370434234381*l5spone*pow(epsc,2) + 0.17451755348802506*l2sone*xispone*pow(epsc,2) + -0.07658370434234381*l5sone*xispone*pow(epsc,2);
        }

        double P3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.8335390223413004 + -0.48142504321106716*as + -0.3657827236798429*epsb + -8.567386867855681*chi3spone*epsb + -2.792280855808401*chi3sppone*epsb + 0.3657827236798429*epsc + -8.567386867855681*chi3spone*epsc + -2.792280855808401*chi3sppone*epsc + 2.1418467169639204*xispone + -0.3211825990576375*as*xispone + -0.9399086363529615*epsb*xispone + -5.584561711616802*chi3spone*epsb*xispone + 0.9399086363529615*epsc*xispone + -5.584561711616802*chi3spone*epsc*xispone + 0.6980702139521002*xisppone + -0.05554148724770942*as*xisppone + -0.30633481736937523*epsb*xisppone + 0.30633481736937523*epsc*xisppone + 0.8335390223413004*l2sone*pow(epsc,2) + 1.7928116099878706*l2spone*pow(epsc,2) + -0.3657827236798429*l5sone*pow(epsc,2) + -0.7867412276682738*l5spone*pow(epsc,2) + 2.1418467169639204*l2sone*xispone*pow(epsc,2) + 1.3961404279042005*l2spone*xispone*pow(epsc,2) + -0.9399086363529615*l5sone*xispone*pow(epsc,2) + -0.6126696347387505*l5spone*xispone*pow(epsc,2) + 0.6980702139521002*l2sone*xisppone*pow(epsc,2) + -0.30633481736937523*l5sone*xisppone*pow(epsc,2);
        }

        double V2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004610740186162548 + 0.0016556537960441093*as + -0.0020233355107207897*epsb + 0.0020233355107207897*epsc + (0.004610740186162548*l2sone - 0.0020233355107207897*l5sone)*pow(epsc,2);
        }

        double V2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04453670166388495 + -0.019544083242286438*epsb + -0.14754368595720155*chi3spone*epsb + 0.019544083242286438*epsc + -0.14754368595720155*chi3spone*epsc + as*(0.01974474150231262 + 0.013245230368352874*xispone) + 0.03688592148930039*xispone + -0.01618668408576632*epsb*xispone + 0.01618668408576632*epsc*xispone + (0.03688592148930039*l2spone - 0.019544083242286438*l5sone - 0.01618668408576632*l5spone + l2sone*(0.04453670166388495 + 0.03688592148930039*xispone) - 0.01618668408576632*l5sone*xispone)*pow(epsc,2);
        }

        double V2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.14480511659548714 + -0.0635449673397373*epsb + -1.7202618251587218*chi3spone*epsb + -0.5901747438288062*chi3sppone*epsb + 0.0635449673397373*epsc + -1.7202618251587218*chi3spone*epsc + -0.5901747438288062*chi3sppone*epsc + 0.43006545628968046*xispone + -0.18872603410982414*epsb*xispone + -1.1803494876576124*chi3spone*epsb*xispone + 0.18872603410982414*epsc*xispone + -1.1803494876576124*chi3spone*epsc*xispone + as*(0.045600901401397734 + 0.18444839275520675*xispone + 0.052980921473411505*xisppone) + 0.14754368595720155*xisppone + -0.06474673634306527*epsb*xisppone + 0.06474673634306527*epsc*xisppone + (0.3562936133110796*l2spone - 0.0635449673397373*l5sone - 0.15635266593829147*l5spone + 0.2950873719144031*l2spone*xispone - 0.18872603410982414*l5sone*xispone - 0.12949347268613054*l5spone*xispone + l2sone*(0.14480511659548714 + 0.43006545628968046*xispone + 0.14754368595720155*xisppone) - 0.06474673634306527*l5sone*xisppone)*pow(epsc,2);
        }

        double V3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.0032602856519248634 + 0.001170724026480039*as + -0.001430714260246217*epsb + 0.001430714260246217*epsc + -0.002861428520492434*epsb*etasone + 0.002861428520492434*epsc*etasone + (0.0032602856519248634*l2sone + 0.001430714260246217*l5sone - 0.002861428520492434*l6sone)*pow(epsc,2);
        }

        double V3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03149220375821526 + -0.013819753792695108*epsb + 0.10432914086159562*chi2sone*epsb + -0.10432914086159563*chi3spone*epsb + 0.013819753792695108*epsc + 0.10432914086159562*chi2sone*epsc + -0.10432914086159563*chi3spone*epsc + -0.027639507585390216*epsb*etasone + 0.027639507585390216*epsc*etasone + -0.02289142816393947*epsb*etaspone + 0.02289142816393947*epsc*etaspone + as*(0.01396164060906072 + 0.009365792211840312*xispone) + 0.026082285215398907*xispone + -0.011445714081969736*epsb*xispone + 0.011445714081969736*epsc*xispone + -0.02289142816393947*epsb*etasone*xispone + 0.02289142816393947*epsc*etasone*xispone + (0.026082285215398907*l2spone + 0.026082285215398904*l3sone + 0.013819753792695108*l5sone + 0.011445714081969736*l5spone - 0.03908522166735995*l6sone - 0.02289142816393947*l6spone + l2sone*(0.03149220375821526 + 0.026082285215398907*xispone) + 0.011445714081969736*l5sone*xispone - 0.02289142816393947*l6sone*xispone)*pow(epsc,2);
        }

        double V3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.10239267989517772 + -0.04493307731620595*epsb + 1.2164088019860795*chi2sone*epsb + 0.8346331268927651*chi2spone*epsb + -1.2164088019860797*chi3spone*epsb + -0.4173165634463825*chi3sppone*epsb + 0.04493307731620595*epsc + 1.2164088019860795*chi2sone*epsc + 0.8346331268927651*chi2spone*epsc + -1.2164088019860797*chi3spone*epsc + -0.4173165634463825*chi3sppone*epsc + -0.0898661546324119*epsb*etasone + 0.0898661546324119*epsc*etasone + -0.2668989170110007*epsb*etaspone + 0.2668989170110007*epsc*etaspone + -0.09156571265575791*epsb*etasppone + 0.09156571265575791*epsc*etasppone + 0.3041022004965199*xispone + -0.13344945850550036*epsb*xispone + 0.8346331268927651*chi2sone*epsb*xispone + -0.834633126892765*chi3spone*epsb*xispone + 0.13344945850550036*epsc*xispone + 0.8346331268927651*chi2sone*epsc*xispone + -0.834633126892765*chi3spone*epsc*xispone + -0.2668989170110007*epsb*etasone*xispone + 0.2668989170110007*epsc*etasone*xispone + -0.18313142531151583*epsb*etaspone*xispone + 0.18313142531151583*epsc*etaspone*xispone + as*(0.032244706609147494 + 0.13042470929616637*xispone + 0.03746316884736126*xisppone) + 0.10432914086159563*xisppone + -0.045782856327878964*epsb*xisppone + 0.045782856327878964*epsc*xisppone + -0.09156571265575791*epsb*etasone*xisppone + 0.09156571265575791*epsc*etasone*xisppone + (0.25193763006572206*l2spone + 0.30410220049651987*l3sone + 0.20865828172319129*l3spone + 0.04493307731620595*l5sone + 0.11055803034156086*l5spone - 0.22331561313791226*l6sone - 0.3126817733388797*l6spone + 0.20865828172319126*l2spone*xispone + 0.20865828172319129*l3sone*xispone + 0.13344945850550038*l5sone*xispone + 0.0915657126557579*l5spone*xispone - 0.3584646296667586*l6sone*xispone - 0.18313142531151583*l6spone*xispone + l2sone*(0.10239267989517772 + 0.3041022004965199*xispone + 0.10432914086159563*xisppone) + 0.04578285632787895*l5sone*xisppone - 0.09156571265575791*l6sone*xisppone)*pow(epsc,2);
        }

        double V6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004143071140512856 + 0.0037774243969479295*as + 0.004143071140512856*epsb + 0.004143071140512856*epsc + 0.008286142281025711*epsc*etasone + 0.004143071140512856*l2sone*pow(epsc,2) + 0.004143071140512856*l5sone*pow(epsc,2) + -0.008286142281025711*l6sone*pow(epsc,2);
        }

        double V6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04001932789686092 + 0.042932871310228105*as + 0.04001932789686092*epsb + -0.13257827649641138*chi3spone*epsb + 0.04001932789686092*epsc + 0.13257827649641138*chi2sone*epsc + -0.13257827649641138*chi3spone*epsc + 0.08003865579372184*epsc*etasone + 0.06628913824820569*epsc*etaspone + 0.033144569124102845*xispone + 0.030219395175583436*as*xispone + 0.033144569124102845*epsb*xispone + 0.033144569124102845*epsc*xispone + 0.06628913824820569*epsc*etasone*xispone + 0.04001932789686092*l2sone*pow(epsc,2) + 0.033144569124102845*l2spone*pow(epsc,2) + 0.033144569124102845*l3sone*pow(epsc,2) + 0.04001932789686092*l5sone*pow(epsc,2) + 0.033144569124102845*l5spone*pow(epsc,2) + -0.11318322491782468*l6sone*pow(epsc,2) + -0.06628913824820569*l6spone*pow(epsc,2) + 0.033144569124102845*l2sone*xispone*pow(epsc,2) + 0.033144569124102845*l5sone*xispone*pow(epsc,2) + -0.06628913824820569*l6sone*xispone*pow(epsc,2);
        }

        double V6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.13011748121610847 + 0.14381800621477134*as + 0.13011748121610847*epsb + -1.5457750456923722*chi3spone*epsb + -0.5303131059856455*chi3sppone*epsb + 0.13011748121610847*epsc + 1.5457750456923722*chi2sone*epsc + 1.060626211971291*chi2spone*epsc + -1.5457750456923722*chi3spone*epsc + -0.5303131059856455*chi3sppone*epsc + 0.26023496243221694*epsc*etasone + 0.7728875228461861*epsc*etaspone + 0.26515655299282276*epsc*etasppone + 0.38644376142309306*xispone + 0.40390176083299173*as*xispone + 0.38644376142309306*epsb*xispone + -1.060626211971291*chi3spone*epsb*xispone + 0.38644376142309306*epsc*xispone + 1.060626211971291*chi2sone*epsc*xispone + -1.060626211971291*chi3spone*epsc*xispone + 0.7728875228461861*epsc*etasone*xispone + 0.5303131059856455*epsc*etaspone*xispone + 0.13257827649641138*xisppone + 0.12087758070233373*as*xisppone + 0.13257827649641138*epsb*xisppone + 0.13257827649641138*epsc*xisppone + 0.26515655299282276*epsc*etasone*xisppone + 0.13011748121610847*l2sone*pow(epsc,2) + 0.32015462317488735*l2spone*pow(epsc,2) + 0.38644376142309306*l3sone*pow(epsc,2) + 0.26515655299282276*l3spone*pow(epsc,2) + 0.13011748121610847*l5sone*pow(epsc,2) + 0.32015462317488735*l5spone*pow(epsc,2) + -0.64667872385531*l6sone*pow(epsc,2) + -0.9054657993425975*l6spone*pow(epsc,2) + 0.38644376142309306*l2sone*xispone*pow(epsc,2) + 0.26515655299282276*l2spone*xispone*pow(epsc,2) + 0.26515655299282276*l3sone*xispone*pow(epsc,2) + 0.38644376142309306*l5sone*xispone*pow(epsc,2) + 0.26515655299282276*l5spone*xispone*pow(epsc,2) + -1.0380440758390088*l6sone*xispone*pow(epsc,2) + -0.5303131059856455*l6spone*xispone*pow(epsc,2) + 0.13257827649641138*l2sone*xisppone*pow(epsc,2) + 0.13257827649641138*l5sone*xisppone*pow(epsc,2) + -0.26515655299282276*l6sone*xisppone*pow(epsc,2);
        }

        double V7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004143071140512856 + 0.0037774243969479295*as + 0.004143071140512856*epsb + 0.004143071140512856*epsc + 0.008286142281025711*epsb*etasone + 0.004143071140512856*l2sone*pow(epsc,2) + -0.004143071140512856*l5sone*pow(epsc,2);
        }

        double V7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04001932789686092 + 0.042932871310228105*as + 0.04001932789686092*epsb + 0.13257827649641138*chi2sone*epsb + -0.13257827649641138*chi3spone*epsb + 0.04001932789686092*epsc + -0.13257827649641138*chi3spone*epsc + 0.08003865579372184*epsb*etasone + 0.06628913824820569*epsb*etaspone + 0.033144569124102845*xispone + 0.030219395175583436*as*xispone + 0.033144569124102845*epsb*xispone + 0.033144569124102845*epsc*xispone + 0.06628913824820569*epsb*etasone*xispone + 0.04001932789686092*l2sone*pow(epsc,2) + 0.033144569124102845*l2spone*pow(epsc,2) + -0.04001932789686092*l5sone*pow(epsc,2) + -0.033144569124102845*l5spone*pow(epsc,2) + 0.033144569124102845*l2sone*xispone*pow(epsc,2) + -0.033144569124102845*l5sone*xispone*pow(epsc,2);
        }

        double V7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.13011748121610847 + 0.14381800621477134*as + 0.13011748121610847*epsb + 1.5457750456923722*chi2sone*epsb + 1.060626211971291*chi2spone*epsb + -1.5457750456923722*chi3spone*epsb + -0.5303131059856455*chi3sppone*epsb + 0.13011748121610847*epsc + -1.5457750456923722*chi3spone*epsc + -0.5303131059856455*chi3sppone*epsc + 0.26023496243221694*epsb*etasone + 0.7728875228461861*epsb*etaspone + 0.26515655299282276*epsb*etasppone + 0.38644376142309306*xispone + 0.40390176083299173*as*xispone + 0.38644376142309306*epsb*xispone + 1.060626211971291*chi2sone*epsb*xispone + -1.060626211971291*chi3spone*epsb*xispone + 0.38644376142309306*epsc*xispone + -1.060626211971291*chi3spone*epsc*xispone + 0.7728875228461861*epsb*etasone*xispone + 0.5303131059856455*epsb*etaspone*xispone + 0.13257827649641138*xisppone + 0.12087758070233373*as*xisppone + 0.13257827649641138*epsb*xisppone + 0.13257827649641138*epsc*xisppone + 0.26515655299282276*epsb*etasone*xisppone + 0.13011748121610847*l2sone*pow(epsc,2) + 0.32015462317488735*l2spone*pow(epsc,2) + -0.13011748121610847*l5sone*pow(epsc,2) + -0.32015462317488735*l5spone*pow(epsc,2) + 0.38644376142309306*l2sone*xispone*pow(epsc,2) + 0.26515655299282276*l2spone*xispone*pow(epsc,2) + -0.38644376142309306*l5sone*xispone*pow(epsc,2) + -0.26515655299282276*l5spone*xispone*pow(epsc,2) + 0.13257827649641138*l2sone*xisppone*pow(epsc,2) + -0.13257827649641138*l5sone*xisppone*pow(epsc,2);
        }

        double A3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0032655729460533954 + -0.001376727455751547*as + 0.0032655729460533954*l2sone*pow(epsc,2);
        }

        double A3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03959827548987528 + -2.9352316138022542e-6*as + 0.013062291784213581*epsb + 0.10449833427370865*chi2sone*epsb + -0.10449833427370865*chi3spone*epsb + 0.013062291784213581*epsc + -0.10449833427370865*chi3spone*epsc + 0.026124583568427163*epsb*etasone + 0.026124583568427163*xispone + -0.011013819646012377*as*xispone + 0.03959827548987528*l2sone*pow(epsc,2) + 0.026124583568427163*l2spone*pow(epsc,2) + -0.013062291784213581*l5sone*pow(epsc,2) + 0.026124583568427163*l2sone*xispone*pow(epsc,2);
        }

        double A3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1946223526057385 + 0.0771916390118414*as + 0.13226851839107395*epsb + 1.4761414842234262*chi2sone*epsb + 0.8359866741896693*chi2spone*epsb + -1.4761414842234262*chi3spone*epsb + -0.41799333709483466*chi3sppone*epsb + 0.13226851839107395*epsc + -1.4761414842234262*chi3spone*epsc + -0.41799333709483466*chi3sppone*epsc + 0.2645370367821479*epsb*etasone + 0.20899666854741733*epsb*etaspone + 0.36903537105585654*xispone + -0.02205112114493516*as*xispone + 0.10449833427370867*epsb*xispone + 0.8359866741896693*chi2sone*epsb*xispone + -0.8359866741896693*chi3spone*epsb*xispone + 0.10449833427370867*epsc*xispone + -0.8359866741896693*chi3spone*epsc*xispone + 0.20899666854741733*epsb*etasone*xispone + 0.10449833427370867*xisppone + -0.044055278584049506*as*xisppone + 0.1946223526057385*l2sone*pow(epsc,2) + 0.31678620391900225*l2spone*pow(epsc,2) + -0.13226851839107395*l5sone*pow(epsc,2) + -0.10449833427370867*l5spone*pow(epsc,2) + 0.36903537105585654*l2sone*xispone*pow(epsc,2) + 0.20899666854741733*l2spone*xispone*pow(epsc,2) + -0.10449833427370867*l5sone*xispone*pow(epsc,2) + 0.10449833427370867*l2sone*xisppone*pow(epsc,2);
        }

        double A4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0032655729460533954 + -0.001376727455751547*as + 0.0032655729460533954*l2sone*pow(epsc,2);
        }

        double A4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03959827548987528 + -2.9352316138022542e-6*as + 0.013062291784213581*epsb + -0.10449833427370865*chi3spone*epsb + 0.013062291784213581*epsc + 0.10449833427370865*chi2sone*epsc + -0.10449833427370865*chi3spone*epsc + 0.026124583568427163*epsc*etasone + 0.026124583568427163*xispone + -0.011013819646012377*as*xispone + 0.03959827548987528*l2sone*pow(epsc,2) + 0.026124583568427163*l2spone*pow(epsc,2) + 0.026124583568427163*l3sone*pow(epsc,2) + 0.013062291784213581*l5sone*pow(epsc,2) + -0.026124583568427163*l6sone*pow(epsc,2) + 0.026124583568427163*l2sone*xispone*pow(epsc,2);
        }

        double A4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1946223526057385 + 0.0771916390118414*as + 0.13226851839107395*epsb + -1.4761414842234262*chi3spone*epsb + -0.41799333709483466*chi3sppone*epsb + 0.13226851839107395*epsc + 1.4761414842234262*chi2sone*epsc + 0.8359866741896693*chi2spone*epsc + -1.4761414842234262*chi3spone*epsc + -0.41799333709483466*chi3sppone*epsc + 0.2645370367821479*epsc*etasone + 0.20899666854741733*epsc*etaspone + 0.36903537105585654*xispone + -0.02205112114493516*as*xispone + 0.10449833427370867*epsb*xispone + -0.8359866741896693*chi3spone*epsb*xispone + 0.10449833427370867*epsc*xispone + 0.8359866741896693*chi2sone*epsc*xispone + -0.8359866741896693*chi3spone*epsc*xispone + 0.20899666854741733*epsc*etasone*xispone + 0.10449833427370867*xisppone + -0.044055278584049506*as*xisppone + 0.1946223526057385*l2sone*pow(epsc,2) + 0.31678620391900225*l2spone*pow(epsc,2) + 0.36903537105585654*l3sone*pow(epsc,2) + 0.20899666854741733*l3spone*pow(epsc,2) + 0.13226851839107395*l5sone*pow(epsc,2) + 0.10449833427370867*l5spone*pow(epsc,2) + -0.36903537105585654*l6sone*pow(epsc,2) + -0.20899666854741733*l6spone*pow(epsc,2) + 0.36903537105585654*l2sone*xispone*pow(epsc,2) + 0.20899666854741733*l2spone*xispone*pow(epsc,2) + 0.20899666854741733*l3sone*xispone*pow(epsc,2) + 0.10449833427370867*l5sone*xispone*pow(epsc,2) + -0.20899666854741733*l6sone*xispone*pow(epsc,2) + 0.10449833427370867*l2sone*xisppone*pow(epsc,2);
        }

        double A7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0030191876237791126 + -0.0012728542783726298*as + 0.0030191876237791126*l2sone*pow(epsc,2);
        }

        double A7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.030733535458538706 + 0.004397061025152454*as + -0.02752027952410895*epsb + -0.09661400396093159*chi3spone*epsb + 0.02752027952410895*epsc + -0.09661400396093159*chi3spone*epsc + 0.024153500990232897*xispone + -0.010182834226981036*as*xispone + 0.030733535458538706*l2sone*pow(epsc,2) + 0.024153500990232897*l2spone*pow(epsc,2) + -0.02752027952410895*l5sone*pow(epsc,2) + 0.024153500990232897*l2sone*xispone*pow(epsc,2);
        }

        double A7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.11439730382532164 + 0.08500224488347237*as + -0.22509952894155016*epsb + -1.176701142595102*chi3spone*epsb + -0.3864560158437264*chi3sppone*epsb + 0.22509952894155016*epsc + -1.176701142595102*chi3spone*epsc + -0.3864560158437264*chi3sppone*epsc + 0.2941752856487755*xispone + 0.014810819747257556*as*xispone + -0.22016223619287162*epsb*xispone + -0.7729120316874528*chi3spone*epsb*xispone + 0.22016223619287162*epsc*xispone + -0.7729120316874528*chi3spone*epsc*xispone + 0.0966140039609316*xisppone + -0.04073133690792415*as*xisppone + 0.11439730382532164*l2sone*pow(epsc,2) + 0.2458682836683097*l2spone*pow(epsc,2) + -0.22509952894155016*l5sone*pow(epsc,2) + -0.22016223619287162*l5spone*pow(epsc,2) + 0.2941752856487755*l2sone*xispone*pow(epsc,2) + 0.1932280079218632*l2spone*xispone*pow(epsc,2) + -0.22016223619287162*l5sone*xispone*pow(epsc,2) + 0.0966140039609316*l2sone*xisppone*pow(epsc,2);
        }
        // }}}
    };

    BGLCoefficients::BGLCoefficients(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BGLCoefficients>(new Implementation<BGLCoefficients>(p, o, *this))
    {
    }

    BGLCoefficients::~BGLCoefficients() = default;

    // B -> D form factors
    // {{{
    double BGLCoefficients::V1_a0() const { return _imp->V1_a0(); }
    double BGLCoefficients::V1_a1() const { return _imp->V1_a1(); }
    double BGLCoefficients::V1_a2() const { return _imp->V1_a2(); }
    double BGLCoefficients::S1_a0() const { return _imp->S1_a0(); }
    double BGLCoefficients::S1_a1() const { return _imp->S1_a1(); }
    double BGLCoefficients::S1_a2() const { return _imp->S1_a2(); }
    // }}}

    // B -> D^* form factors
    // {{{
    double BGLCoefficients::A1_a0() const { return _imp->A1_a0(); }
    double BGLCoefficients::A1_a1() const { return _imp->A1_a1(); }
    double BGLCoefficients::A1_a2() const { return _imp->A1_a2(); }
    double BGLCoefficients::A5_a0() const { return _imp->A5_a0(); }
    double BGLCoefficients::A5_a1() const { return _imp->A5_a1(); }
    double BGLCoefficients::A5_a2() const { return _imp->A5_a2(); }
    double BGLCoefficients::V4_a0() const { return _imp->V4_a0(); }
    double BGLCoefficients::V4_a1() const { return _imp->V4_a1(); }
    double BGLCoefficients::V4_a2() const { return _imp->V4_a2(); }
    double BGLCoefficients::P1_a0() const { return _imp->P1_a0(); }
    double BGLCoefficients::P1_a1() const { return _imp->P1_a1(); }
    double BGLCoefficients::P1_a2() const { return _imp->P1_a2(); }
    // }}}

    // B^* -> D form factors
    // {{{
    double BGLCoefficients::P2_a0() const { return _imp->P2_a0(); }
    double BGLCoefficients::P2_a1() const { return _imp->P2_a1(); }
    double BGLCoefficients::P2_a2() const { return _imp->P2_a2(); }
    double BGLCoefficients::V5_a0() const { return _imp->V5_a0(); }
    double BGLCoefficients::V5_a1() const { return _imp->V5_a1(); }
    double BGLCoefficients::V5_a2() const { return _imp->V5_a2(); }
    double BGLCoefficients::A2_a0() const { return _imp->A2_a0(); }
    double BGLCoefficients::A2_a1() const { return _imp->A2_a1(); }
    double BGLCoefficients::A2_a2() const { return _imp->A2_a2(); }
    double BGLCoefficients::A6_a0() const { return _imp->A6_a0(); }
    double BGLCoefficients::A6_a1() const { return _imp->A6_a1(); }
    double BGLCoefficients::A6_a2() const { return _imp->A6_a2(); }
    // }}}

    // B^* -> D^* form factors
    // {{{
    double BGLCoefficients::S2_a0() const { return _imp->S2_a0(); }
    double BGLCoefficients::S2_a1() const { return _imp->S2_a1(); }
    double BGLCoefficients::S2_a2() const { return _imp->S2_a2(); }
    double BGLCoefficients::S3_a0() const { return _imp->S3_a0(); }
    double BGLCoefficients::S3_a1() const { return _imp->S3_a1(); }
    double BGLCoefficients::S3_a2() const { return _imp->S3_a2(); }
    double BGLCoefficients::P3_a0() const { return _imp->P3_a0(); }
    double BGLCoefficients::P3_a1() const { return _imp->P3_a1(); }
    double BGLCoefficients::P3_a2() const { return _imp->P3_a2(); }
    double BGLCoefficients::V2_a0() const { return _imp->V2_a0(); }
    double BGLCoefficients::V2_a1() const { return _imp->V2_a1(); }
    double BGLCoefficients::V2_a2() const { return _imp->V2_a2(); }
    double BGLCoefficients::V3_a0() const { return _imp->V3_a0(); }
    double BGLCoefficients::V3_a1() const { return _imp->V3_a1(); }
    double BGLCoefficients::V3_a2() const { return _imp->V3_a2(); }
    double BGLCoefficients::V6_a0() const { return _imp->V6_a0(); }
    double BGLCoefficients::V6_a1() const { return _imp->V6_a1(); }
    double BGLCoefficients::V6_a2() const { return _imp->V6_a2(); }
    double BGLCoefficients::V7_a0() const { return _imp->V7_a0(); }
    double BGLCoefficients::V7_a1() const { return _imp->V7_a1(); }
    double BGLCoefficients::V7_a2() const { return _imp->V7_a2(); }
    double BGLCoefficients::A3_a0() const { return _imp->A3_a0(); }
    double BGLCoefficients::A3_a1() const { return _imp->A3_a1(); }
    double BGLCoefficients::A3_a2() const { return _imp->A3_a2(); }
    double BGLCoefficients::A4_a0() const { return _imp->A4_a0(); }
    double BGLCoefficients::A4_a1() const { return _imp->A4_a1(); }
    double BGLCoefficients::A4_a2() const { return _imp->A4_a2(); }
    double BGLCoefficients::A7_a0() const { return _imp->A7_a0(); }
    double BGLCoefficients::A7_a1() const { return _imp->A7_a1(); }
    double BGLCoefficients::A7_a2() const { return _imp->A7_a2(); }
    // }}}

    // B_s -> D_s form factors
    // {{{
    double BGLCoefficients::V1s_a0() const { return _imp->V1s_a0(); }
    double BGLCoefficients::V1s_a1() const { return _imp->V1s_a1(); }
    double BGLCoefficients::V1s_a2() const { return _imp->V1s_a2(); }
    double BGLCoefficients::S1s_a0() const { return _imp->S1s_a0(); }
    double BGLCoefficients::S1s_a1() const { return _imp->S1s_a1(); }
    double BGLCoefficients::S1s_a2() const { return _imp->S1s_a2(); }
    // }}}

    // B_s -> D_s^* form factors
    // {{{
    double BGLCoefficients::A1s_a0() const { return _imp->A1s_a0(); }
    double BGLCoefficients::A1s_a1() const { return _imp->A1s_a1(); }
    double BGLCoefficients::A1s_a2() const { return _imp->A1s_a2(); }
    double BGLCoefficients::A5s_a0() const { return _imp->A5s_a0(); }
    double BGLCoefficients::A5s_a1() const { return _imp->A5s_a1(); }
    double BGLCoefficients::A5s_a2() const { return _imp->A5s_a2(); }
    double BGLCoefficients::V4s_a0() const { return _imp->V4s_a0(); }
    double BGLCoefficients::V4s_a1() const { return _imp->V4s_a1(); }
    double BGLCoefficients::V4s_a2() const { return _imp->V4s_a2(); }
    double BGLCoefficients::P1s_a0() const { return _imp->P1s_a0(); }
    double BGLCoefficients::P1s_a1() const { return _imp->P1s_a1(); }
    double BGLCoefficients::P1s_a2() const { return _imp->P1s_a2(); }
    // }}}

    // B_s^* -> D_s form factors
    // {{{
    double BGLCoefficients::P2s_a0() const { return _imp->P2s_a0(); }
    double BGLCoefficients::P2s_a1() const { return _imp->P2s_a1(); }
    double BGLCoefficients::P2s_a2() const { return _imp->P2s_a2(); }
    double BGLCoefficients::V5s_a0() const { return _imp->V5s_a0(); }
    double BGLCoefficients::V5s_a1() const { return _imp->V5s_a1(); }
    double BGLCoefficients::V5s_a2() const { return _imp->V5s_a2(); }
    double BGLCoefficients::A2s_a0() const { return _imp->A2s_a0(); }
    double BGLCoefficients::A2s_a1() const { return _imp->A2s_a1(); }
    double BGLCoefficients::A2s_a2() const { return _imp->A2s_a2(); }
    double BGLCoefficients::A6s_a0() const { return _imp->A6s_a0(); }
    double BGLCoefficients::A6s_a1() const { return _imp->A6s_a1(); }
    double BGLCoefficients::A6s_a2() const { return _imp->A6s_a2(); }
    // }}}

    // B_s^* -> D_s^* form factors
    // {{{
    double BGLCoefficients::S2s_a0() const { return _imp->S2s_a0(); }
    double BGLCoefficients::S2s_a1() const { return _imp->S2s_a1(); }
    double BGLCoefficients::S2s_a2() const { return _imp->S2s_a2(); }
    double BGLCoefficients::S3s_a0() const { return _imp->S3s_a0(); }
    double BGLCoefficients::S3s_a1() const { return _imp->S3s_a1(); }
    double BGLCoefficients::S3s_a2() const { return _imp->S3s_a2(); }
    double BGLCoefficients::P3s_a0() const { return _imp->P3s_a0(); }
    double BGLCoefficients::P3s_a1() const { return _imp->P3s_a1(); }
    double BGLCoefficients::P3s_a2() const { return _imp->P3s_a2(); }
    double BGLCoefficients::V2s_a0() const { return _imp->V2s_a0(); }
    double BGLCoefficients::V2s_a1() const { return _imp->V2s_a1(); }
    double BGLCoefficients::V2s_a2() const { return _imp->V2s_a2(); }
    double BGLCoefficients::V3s_a0() const { return _imp->V3s_a0(); }
    double BGLCoefficients::V3s_a1() const { return _imp->V3s_a1(); }
    double BGLCoefficients::V3s_a2() const { return _imp->V3s_a2(); }
    double BGLCoefficients::V6s_a0() const { return _imp->V6s_a0(); }
    double BGLCoefficients::V6s_a1() const { return _imp->V6s_a1(); }
    double BGLCoefficients::V6s_a2() const { return _imp->V6s_a2(); }
    double BGLCoefficients::V7s_a0() const { return _imp->V7s_a0(); }
    double BGLCoefficients::V7s_a1() const { return _imp->V7s_a1(); }
    double BGLCoefficients::V7s_a2() const { return _imp->V7s_a2(); }
    double BGLCoefficients::A3s_a0() const { return _imp->A3s_a0(); }
    double BGLCoefficients::A3s_a1() const { return _imp->A3s_a1(); }
    double BGLCoefficients::A3s_a2() const { return _imp->A3s_a2(); }
    double BGLCoefficients::A4s_a0() const { return _imp->A4s_a0(); }
    double BGLCoefficients::A4s_a1() const { return _imp->A4s_a1(); }
    double BGLCoefficients::A4s_a2() const { return _imp->A4s_a2(); }
    double BGLCoefficients::A7s_a0() const { return _imp->A7s_a0(); }
    double BGLCoefficients::A7s_a1() const { return _imp->A7s_a1(); }
    double BGLCoefficients::A7s_a2() const { return _imp->A7s_a2(); }
    // }}}

    template <> struct Implementation<HQETUnitarityBounds>
    {
        // option to determine if we use z^3 terms in the leading-power IW function
        SwitchOption opt_zorder_bound;
        unsigned zorder_bound;

        // number of light flavor multiplets
        UsedParameter nf;
        UsedParameter ns;

        std::shared_ptr<BGLCoefficients> bgl;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_zorder_bound(o, "z-order-bound", { "1", "2" }, "2"),
            nf(p["B(*)->D(*)::n_f@HQET"], u),
            ns(p["B_s(*)->D_s(*)::n_s@HQET"], u),
            bgl(new BGLCoefficients(p, o))
        {
            if ("1" == opt_zorder_bound.value())
            {
                zorder_bound = 1;
            }
            else if ("2" == opt_zorder_bound.value())
            {
                zorder_bound = 2;
            }
            else
            {
                throw InternalError("Only z-order-bound=2 is presently supported");
            }
        }

        ~Implementation() = default;

        // bounds up to z^2
        // {{{
        double bound_0p() const
        {
            // 3 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=u,d
            const std::array<std::array<double, 3>, 3> bgl_coeffs_ud
            {{
                // B -> D S_1
                { bgl->S1_a0(), bgl->S1_a1(), bgl->S1_a2() },
                // B^* -> D^* S_2
                { bgl->S2_a0(), bgl->S2_a1(), bgl->S2_a2() },
                // B^* -> D^* S_3
                { bgl->S3_a0(), bgl->S3_a1(), bgl->S3_a2() }
            }};

            double result = 0.0;
            for (const auto & ff : bgl_coeffs_ud)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * nf; // to account for flavour symmetry
                }
            }

            // 3 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=s
            const std::array<std::array<double, 3>, 3> bgl_coeffs_s
            {{
                // B_s -> D_s S_1
                { bgl->S1s_a0(), bgl->S1s_a1(), bgl->S1s_a2() },
                // B_s^* -> D_s^* S_2
                { bgl->S2s_a0(), bgl->S2s_a1(), bgl->S2s_a2() },
                // B_s^* -> D_s^* S_3
                { bgl->S3s_a0(), bgl->S3s_a1(), bgl->S3s_a2() }
            }};

            for (const auto & ff : bgl_coeffs_s)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * ns;
                }
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 0^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_0m() const
        {
            // 3 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=u,d
            const std::array<std::array<double, 3>, 3> bgl_coeffs_ud
            {{
                // B -> D^* P_1
                { bgl->P1_a0(), bgl->P1_a1(), bgl->P1_a2() },
                // B^* -> D P_2
                { bgl->P2_a0(), bgl->P2_a1(), bgl->P2_a2() },
                // B^* -> D^* P_3
                { bgl->P3_a0(), bgl->P3_a1(), bgl->P3_a2() }
            }};

            double result = 0.0;
            for (const auto & ff : bgl_coeffs_ud)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * nf; // to account for flavour symmetry
                }
            }

            // 3 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=s
            const std::array<std::array<double, 3>, 3> bgl_coeffs_s
            {{
                // B_s -> D_s S_1
                { bgl->P1s_a0(), bgl->P1s_a1(), bgl->P1s_a2() },
                // B_s^* -> D_s^* S_2
                { bgl->P2s_a0(), bgl->P2s_a1(), bgl->P2s_a2() },
                // B_s^* -> D_s^* S_3
                { bgl->P3s_a0(), bgl->P3s_a1(), bgl->P3s_a2() }
            }};

            for (const auto & ff : bgl_coeffs_s)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * ns;
                }
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 0^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // using the same relative uncertainty as for 0^+, cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1p() const
        {
            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=u,d
            const std::array<std::array<double, 3>, 7> bgl_coeffs_ud
            {{
                // B -> D V_1
                { bgl->V1_a0(), bgl->V1_a1(), bgl->V1_a2() },
                // B -> D^* V_2
                { bgl->V2_a0(), bgl->V2_a1(), bgl->V2_a2() },
                // B^* -> D V_3
                { bgl->V3_a0(), bgl->V3_a1(), bgl->V3_a2() },
                // B^* -> D^* V_4
                { bgl->V4_a0(), bgl->V4_a1(), bgl->V4_a2() },
                // B^* -> D^* V_5
                { bgl->V5_a0(), bgl->V5_a1(), bgl->V5_a2() },
                // B^* -> D^* V_6
                { bgl->V6_a0(), bgl->V6_a1(), bgl->V6_a2() },
                // B^* -> D^* V_7
                { bgl->V7_a0(), bgl->V7_a1(), bgl->V7_a2() }
            }};

            double result = 0.0;
            for (const auto & ff : bgl_coeffs_ud)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * nf; // to account for flavour symmetry
                }
            }

            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=s
            const std::array<std::array<double, 3>, 7> bgl_coeffs_s
            {{
                // B_s -> D_s V_1
                { bgl->V1s_a0(), bgl->V1s_a1(), bgl->V1s_a2() },
                // B_s -> D_s^* V_2
                { bgl->V2s_a0(), bgl->V2s_a1(), bgl->V2s_a2() },
                // B_s^* -> D_s V_3
                { bgl->V3s_a0(), bgl->V3s_a1(), bgl->V3s_a2() },
                // B_s^* -> D_s^* V_4
                { bgl->V4s_a0(), bgl->V4s_a1(), bgl->V4s_a2() },
                // B_s^* -> D_s^* V_5
                { bgl->V5s_a0(), bgl->V5s_a1(), bgl->V5s_a2() },
                // B_s^* -> D_s^* V_6
                { bgl->V6s_a0(), bgl->V6s_a1(), bgl->V6s_a2() },
                // B_s^* -> D_s^* V_7
                { bgl->V7s_a0(), bgl->V7s_a1(), bgl->V7s_a2() }
            }};

            for (const auto & ff : bgl_coeffs_s)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * ns;
                }
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 1^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0093549; // cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1m() const
        {
            // 3 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=u,d
            const std::array<std::array<double, 3>, 7> bgl_coeffs_ud
            {{
                // B -> D^* A_1
                { bgl->A1_a0(), bgl->A1_a1(), bgl->A1_a2() },
                // B^* -> D A_2
                { bgl->A2_a0(), bgl->A2_a1(), bgl->A2_a2() },
                // B^* -> D^* A_3
                { bgl->A3_a0(), bgl->A3_a1(), bgl->A3_a2() },
                // B^* -> D^* A_4
                { bgl->A4_a0(), bgl->A4_a1(), bgl->A4_a2() },
                // B -> D^* A_5
                { bgl->A5_a0(), bgl->A5_a1(), bgl->A5_a2() },
                // B^* -> D A_6
                { bgl->A6_a0(), bgl->A6_a1(), bgl->A6_a2() },
                // B^* -> D^* A_7
                { bgl->A7_a0(), bgl->A7_a1(), bgl->A7_a2() }
            }};

            double result = 0.0;
            for (const auto & ff : bgl_coeffs_ud)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * nf; // to account for flavour symmetry
                }
            }

            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=s
            const std::array<std::array<double, 3>, 7> bgl_coeffs_s
            {{
                // B_s -> D_s^* A_1
                { bgl->A1s_a0(), bgl->A1s_a1(), bgl->A1s_a2() },
                // B_s^* -> D_s A_2
                { bgl->A2s_a0(), bgl->A2s_a1(), bgl->A2s_a2() },
                // B_s^* -> D_s^* A_3
                { bgl->A3s_a0(), bgl->A3s_a1(), bgl->A3s_a2() },
                // B_s^* -> D_s^* A_4
                { bgl->A4s_a0(), bgl->A4s_a1(), bgl->A4s_a2() },
                // B_s -> D_s^* A_5
                { bgl->A5s_a0(), bgl->A5s_a1(), bgl->A5s_a2() },
                // B_s^* -> D_s A_6
                { bgl->A6s_a0(), bgl->A6s_a1(), bgl->A6s_a2() },
                // B_s^* -> D_s^* A_7
                { bgl->A7s_a0(), bgl->A7s_a1(), bgl->A7s_a2() }
            }};

            for (const auto & ff : bgl_coeffs_s)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += pow(ff[i], 2) * ns;
                }
            }

            if (result < 0.0)
            {
                throw InternalError("Contribution to 1^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= result) && (result < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0093549; // same relative uncertainty as for 1^-, cf. [BG2016], eq. (2.8), p.5
                return -pow((result - 1.0) / sigma, 2) / 2.0;
            }
        }
        // }}}
    };

    HQETUnitarityBounds::HQETUnitarityBounds(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<HQETUnitarityBounds>(new Implementation<HQETUnitarityBounds>(p, o, *this))
    {
    }

    HQETUnitarityBounds::~HQETUnitarityBounds() = default;

    double
    HQETUnitarityBounds::bound_0p() const
    {
        return _imp->bound_0p();
    }

    double
    HQETUnitarityBounds::bound_0m() const
    {
        return _imp->bound_0m();
    }

    double
    HQETUnitarityBounds::bound_1p() const
    {
        return _imp->bound_1p();
    }

    double
    HQETUnitarityBounds::bound_1m() const
    {
        return _imp->bound_1m();
    }
}
