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

#include <gsl/gsl_sf_dilog.h>

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
        // option to determine if we use the SU3_F-symmetry limit for the subsubleading-power IW functions
        SwitchOption opt_sslp_limit;

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

        static const std::vector<OptionSpecification> options;

        std::string _sslp_prefix()
        {
            if ("1" == opt_sslp_limit.value())
                return "B(*)->D(*)";

            return "B_s(*)->D_s(*)";
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
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
            opt_sslp_limit(o, "SU3F-limit-sslp", { "0", "1" }, "0"),
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
            l1sone(p[_sslp_prefix() + "::l_1(1)@HQET"], u),
            l1spone(p[_sslp_prefix() + "::l_1'(1)@HQET"], u),
            l1sppone(p[_sslp_prefix() + "::l_1''(1)@HQET"], u),
            l2sone(p[_sslp_prefix() + "::l_2(1)@HQET"], u),
            l2spone(p[_sslp_prefix() + "::l_2'(1)@HQET"], u),
            l2sppone(p[_sslp_prefix() + "::l_2''(1)@HQET"], u),
            l3sone(p[_sslp_prefix() + "::l_3(1)@HQET"], u),
            l3spone(p[_sslp_prefix() + "::l_3'(1)@HQET"], u),
            l3sppone(p[_sslp_prefix() + "::l_3''(1)@HQET"], u),
            l4sone(p[_sslp_prefix() + "::l_4(1)@HQET"], u),
            l4spone(p[_sslp_prefix() + "::l_4'(1)@HQET"], u),
            l4sppone(p[_sslp_prefix() + "::l_4''(1)@HQET"], u),
            l5sone(p[_sslp_prefix() + "::l_5(1)@HQET"], u),
            l5spone(p[_sslp_prefix() + "::l_5'(1)@HQET"], u),
            l5sppone(p[_sslp_prefix() + "::l_5''(1)@HQET"], u),
            l6sone(p[_sslp_prefix() + "::l_6(1)@HQET"], u),
            l6spone(p[_sslp_prefix() + "::l_6'(1)@HQET"], u),
            l6sppone(p[_sslp_prefix() + "::l_6''(1)@HQET"], u)
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

            return 0.004853927809740896 + 0.0018073507175946574*as + -0.0022489124823136626*epsb + 0.0022489124823136626*epsc + 0.004497824964627325*epsb*etasone + -0.004497824964627325*epsc*etasone + (0.004853927809740896*l1sone - 0.0022489124823136626*l4sone)*pow(epsc,2);
        }

        double V1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04085920220485458 + -0.018930806855320836*epsb + -0.15532568991170867*chi2sone*epsb + 0.465977069735126*chi3spone*epsb + 0.018930806855320836*epsc + -0.15532568991170867*chi2sone*epsc + 0.465977069735126*chi3spone*epsc + 0.03786161371064167*epsb*etasone + -0.03786161371064167*epsc*etasone + 0.0359825997170186*epsb*etaspone + -0.0359825997170186*epsc*etaspone + as*(0.01923252020198305 + 0.014458805740757258*xispone) + 0.03883142247792717*xispone + -0.0179912998585093*epsb*xispone + 0.0179912998585093*epsc*xispone + 0.0359825997170186*epsb*etasone*xispone + -0.0359825997170186*epsc*etasone*xispone + (0.03883142247792717*l1spone - 0.018930806855320836*l4sone - 0.0179912998585093*l4spone + l1sone*(0.04085920220485458 + 0.03883142247792717*xispone) - 0.0179912998585093*l4sone*xispone)*pow(epsc,2);
        }

        double V1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();


            return 0.11341891835977361 + -0.05254903476684066*epsb + -1.618145850378764*chi2sone*epsb + -1.2426055192936694*chi2spone*epsb + 4.854437551136292*chi3spone*epsb + 1.863908278940504*chi3sppone*epsb + 0.05254903476684066*epsc + -1.618145850378764*chi2sone*epsc + -1.2426055192936694*chi2spone*epsc + 4.854437551136292*chi3spone*epsc + 1.863908278940504*chi3sppone*epsc + 0.10509806953368132*epsb*etasone + -0.10509806953368132*epsc*etasone + 0.3748581091191706*epsb*etaspone + -0.3748581091191706*epsc*etaspone + 0.1439303988680744*epsb*etasppone + -0.1439303988680744*epsc*etasppone + 0.404536462594691*xispone + -0.1874290545595853*epsb*xispone + -1.2426055192936694*chi2sone*epsb*xispone + 3.727816557881008*chi3spone*epsb*xispone + 0.1874290545595853*epsc*xispone + -1.2426055192936694*chi2sone*epsc*xispone + 3.727816557881008*chi3spone*epsc*xispone + 0.3748581091191706*epsb*etasone*xispone + -0.3748581091191706*epsc*etasone*xispone + 0.2878607977361488*epsb*etaspone*xispone + -0.2878607977361488*epsc*etaspone*xispone + as*(0.03101685877692137 + 0.18277777309737894*xispone + 0.05783522296302903*xisppone) + 0.15532568991170867*xisppone + -0.0719651994340372*epsb*xisppone + 0.0719651994340372*epsc*xisppone + 0.1439303988680744*epsb*etasone*xisppone + -0.1439303988680744*epsc*etasone*xisppone + (0.32687361763883666*l1spone - 0.05254903476684066*l4sone - 0.1514464548425667*l4spone + 0.31065137982341734*l1spone*xispone - 0.1874290545595853*l4sone*xispone - 0.1439303988680744*l4spone*xispone + l1sone*(0.11341891835977361 + 0.404536462594691*xispone + 0.15532568991170867*xisppone) - 0.0719651994340372*l4sone*xisppone)*pow(epsc,2);
        }

        double S1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.027568784109983837 + 0.00691165359571278*as + 0.027568784109983837*l1sone*pow(epsc,2);
        }

        double S1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();


            return 0.1991754644782667 + 0.13815036774277073*as + -0.23801173042451831*epsb + -0.8822010915194828*chi2sone*epsb + 2.6466032745584482*chi3spone*epsb + 0.23801173042451831*epsc + -0.8822010915194828*chi2sone*epsc + 2.6466032745584482*chi3spone*epsc + 0.47602346084903663*epsb*etasone + -0.47602346084903663*epsc*etasone + 0.2205502728798707*xispone + 0.05529322876570224*as*xispone + 0.1991754644782667*l1sone*pow(epsc,2) + 0.2205502728798707*l1spone*pow(epsc,2) + -0.23801173042451831*l4sone*pow(epsc,2) + 0.2205502728798707*l1sone*xispone*pow(epsc,2);
        }

        double S1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();


            return 0.5502247632998395 + 0.35019441270807145*as + -1.2435335848826854*epsb + -8.138017046343501*chi2sone*epsb + -7.057608732155861*chi2spone*epsb + 24.414051139030498*chi3spone*epsb + 10.586413098233795*chi3sppone*epsb + 1.2435335848826854*epsc + -8.138017046343501*chi2sone*epsc + -7.057608732155861*chi2spone*epsc + 24.414051139030498*chi3spone*epsc + 10.586413098233795*chi3sppone*epsc + 2.487067169765371*epsb*etasone + -2.487067169765371*epsc*etasone + 3.808187686792293*epsb*etaspone + -3.808187686792293*epsc*etaspone + 2.0345042615858753*xispone + 1.2157893994735705*as*xispone + -1.9040938433961465*epsb*xispone + -7.057608732155861*chi2sone*epsb*xispone + 21.17282619646759*chi3spone*epsb*xispone + 1.9040938433961465*epsc*xispone + -7.057608732155861*chi2sone*epsc*xispone + 21.17282619646759*chi3spone*epsc*xispone + 3.808187686792293*epsb*etasone*xispone + -3.808187686792293*epsc*etasone*xispone + 0.8822010915194827*xisppone + 0.22117291506280895*as*xisppone + 0.5502247632998395*l1sone*pow(epsc,2) + 1.5934037158261338*l1spone*pow(epsc,2) + -1.2435335848826854*l4sone*pow(epsc,2) + -1.9040938433961465*l4spone*pow(epsc,2) + 2.0345042615858753*l1sone*xispone*pow(epsc,2) + 1.7644021830389653*l1spone*xispone*pow(epsc,2) + -1.9040938433961465*l4sone*xispone*pow(epsc,2) + 0.8822010915194827*l1sone*xisppone*pow(epsc,2);
        }
        // }}}

        // B_s -> D_s^* form factors
        // {{{
        double A1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0039029629021226243 + -0.0016234798992337527*as + 0.0039029629021226243*l2sone*pow(epsc,2);
        }

        double A1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.045389619854771975 + 0.0011297414860574153*as + 0.015611851608490495*epsb + -0.12489481286792396*chi2sone*epsb + 0.3746844386037719*chi3spone*epsb + 0.015611851608490495*epsc + -0.12489481286792396*chi3spone*epsc + -0.03122370321698099*epsb*etasone + 0.03122370321698099*xispone + -0.012987839193870022*as*xispone + 0.045389619854771975*l2sone*pow(epsc,2) + 0.03122370321698099*l2spone*pow(epsc,2) + -0.015611851608490495*l5sone*pow(epsc,2) + 0.03122370321698099*l2sone*xispone*pow(epsc,2);
        }

        double A1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.21375918823072088 + 0.09214542757833816*as + 0.15033477620210686*epsb + -1.7022574610885512*chi2sone*epsb + -0.9991585029433917*chi2spone*epsb + 5.106772383265654*chi3spone*epsb + 1.4987377544150877*chi3sppone*epsb + 0.15033477620210686*epsc + -1.7022574610885512*chi3spone*epsc + -0.49957925147169585*chi3sppone*epsc + -0.3006695524042137*epsb*etasone + -0.24978962573584793*epsb*etaspone + 0.4255643652721378*xispone + -0.01693774649928072*as*xispone + 0.12489481286792396*epsb*xispone + -0.9991585029433917*chi2sone*epsb*xispone + 2.9974755088301754*chi3spone*epsb*xispone + 0.12489481286792396*epsc*xispone + -0.9991585029433917*chi3spone*epsc*xispone + -0.24978962573584793*epsb*etasone*xispone + 0.12489481286792396*xisppone + -0.05195135677548009*as*xisppone + 0.21375918823072088*l2sone*pow(epsc,2) + 0.3631169588381758*l2spone*pow(epsc,2) + -0.15033477620210686*l5sone*pow(epsc,2) + -0.12489481286792396*l5spone*pow(epsc,2) + 0.4255643652721378*l2sone*xispone*pow(epsc,2) + 0.24978962573584793*l2spone*xispone*pow(epsc,2) + -0.12489481286792396*l5sone*xispone*pow(epsc,2) + 0.12489481286792396*l2sone*xisppone*pow(epsc,2);
        }

        double A5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.002527938212054363 + -0.0010515234135438955*as + 0.002527938212054363*l2sone*pow(epsc,2);
        }

        double A5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.024475416515784738 + -0.023236265467587044*epsb + -0.08089402278573964*chi2sone*epsb + 0.2426820683572189*chi3spone*epsb + 0.023236265467587044*epsc + 0.08089402278573964*chi2sone*epsc + -0.08089402278573964*chi3spone*epsc + 0.04647253093517409*epsb*etasone + 0.04647253093517409*epsc*etasone + as*(0.004718133656994875 - 0.008412187308351166*xispone) + 0.02022350569643491*xispone + (0.02022350569643491*l2spone + 0.02022350569643491*l3sone + 0.023236265467587044*l5sone - 0.04647253093517409*l6sone + l2sone*(0.024475416515784738 + 0.02022350569643491*xispone))*pow(epsc,2);
        }

        double A5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.08599311230540732 + -0.17850024445887427*epsb + -0.9450013740765909*chi2sone*epsb + -0.6471521822859172*chi2spone*epsb + 2.8350041222297726*chi3spone*epsb + 0.9707282734288758*chi3sppone*epsb + 0.1785002444588743*epsc + 0.9450013740765909*chi2sone*epsc + 0.6471521822859172*chi2spone*epsc + -0.9450013740765909*chi3spone*epsc + -0.3235760911429586*chi3sppone*epsc + 0.35700048891774855*epsb*etasone + 0.3570004889177486*epsc*etasone + 0.3717802474813927*epsb*etaspone + 0.3717802474813927*epsc*etaspone + 0.23625034351914773*xispone + -0.18589012374069636*epsb*xispone + -0.6471521822859172*chi2sone*epsb*xispone + 1.9414565468577516*chi3spone*epsb*xispone + 0.18589012374069636*epsc*xispone + 0.6471521822859172*chi2sone*epsc*xispone + -0.6471521822859172*chi3spone*epsc*xispone + 0.3717802474813927*epsb*etasone*xispone + 0.3717802474813927*epsc*etasone*xispone + as*(0.0719079842270774 + 0.02092069463925666*xispone - 0.03364874923340466*xisppone) + 0.08089402278573965*xisppone + (0.1958033321262779*l2spone + 0.23625034351914773*l3sone + 0.1617880455714793*l3spone + 0.17850024445887427*l5sone + 0.18589012374069636*l5spone - 0.5428906126584448*l6sone - 0.3717802474813927*l6spone + 0.1617880455714793*l2spone*xispone + 0.1617880455714793*l3sone*xispone + 0.18589012374069636*l5sone*xispone - 0.3717802474813927*l6sone*xispone + l2sone*(0.08599311230540732 + 0.23625034351914773*xispone + 0.08089402278573965*xisppone))*pow(epsc,2);
        }

        double V4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004684646612240956 + 0.004297565736635373*as + 0.004684646612240956*epsb + 0.004684646612240956*epsc + -0.009369293224481911*epsb*etasone + 0.004684646612240956*l2sone*pow(epsc,2) + -0.004684646612240956*l5sone*pow(epsc,2);
        }

        double V4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04345453602144384 + 0.04722512441781111*as + 0.04345453602144384*epsb + -0.14990869159171058*chi2sone*epsb + 0.4497260747751317*chi3spone*epsb + 0.04345453602144384*epsc + -0.14990869159171058*chi3spone*epsc + -0.08690907204288768*epsb*etasone + -0.07495434579585529*epsb*etaspone + 0.037477172897927645*xispone + 0.034380525893082985*as*xispone + 0.037477172897927645*epsb*xispone + 0.037477172897927645*epsc*xispone + -0.07495434579585529*epsb*etasone*xispone + 0.04345453602144384*l2sone*pow(epsc,2) + 0.037477172897927645*l2spone*pow(epsc,2) + -0.04345453602144384*l5sone*pow(epsc,2) + -0.037477172897927645*l5spone*pow(epsc,2) + 0.037477172897927645*l2sone*xispone*pow(epsc,2) + -0.037477172897927645*l5sone*xispone*pow(epsc,2);
        }

        double V4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.13494178999195733 + 0.15007087537821776*as + 0.13494178999195733*epsb + -1.6903625358696244*chi2sone*epsb + -1.1992695327336846*chi2spone*epsb + 5.071087607608872*chi3spone*epsb + 1.7989042991005268*chi3sppone*epsb + 0.13494178999195733*epsc + -1.6903625358696244*chi3spone*epsc + -0.5996347663668423*chi3sppone*epsc + -0.26988357998391466*epsb*etasone + -0.8451812679348122*epsb*etaspone + -0.29981738318342116*epsb*etasppone + 0.4225906339674061*xispone + 0.44656204712865494*as*xispone + 0.4225906339674061*epsb*xispone + -1.1992695327336846*chi2sone*epsb*xispone + 3.5978085982010537*chi3spone*epsb*xispone + 0.4225906339674061*epsc*xispone + -1.1992695327336846*chi3spone*epsc*xispone + -0.8451812679348122*epsb*etasone*xispone + -0.5996347663668423*epsb*etaspone*xispone + 0.14990869159171058*xisppone + 0.13752210357233197*as*xisppone + 0.14990869159171058*epsb*xisppone + 0.14990869159171058*epsc*xisppone + -0.29981738318342116*epsb*etasone*xisppone + 0.13494178999195733*l2sone*pow(epsc,2) + 0.34763628817155073*l2spone*pow(epsc,2) + -0.13494178999195733*l5sone*pow(epsc,2) + -0.34763628817155073*l5spone*pow(epsc,2) + 0.4225906339674061*l2sone*xispone*pow(epsc,2) + 0.29981738318342116*l2spone*xispone*pow(epsc,2) + -0.4225906339674061*l5sone*xispone*pow(epsc,2) + -0.29981738318342116*l5spone*xispone*pow(epsc,2) + 0.14990869159171058*l2sone*xisppone*pow(epsc,2) + -0.14990869159171058*l5sone*xisppone*pow(epsc,2);
        }

        double P1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.017614604136229357 + -0.0012894515619258611*as + -0.007665367904889766*epsb + 0.007665367904889766*epsc + 0.015330735809779531*epsb*etasone + 0.015330735809779531*epsc*etasone + 0.017614604136229357*l2sone*pow(epsc,2) + 0.007665367904889766*l5sone*pow(epsc,2) + -0.015330735809779531*l6sone*pow(epsc,2);
        }

        double P1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.17385254509323594 + -0.027729362113336223*as + -0.07565561559229919*epsb + -0.5636673323593394*chi2sone*epsb + 1.6910019970780181*chi3spone*epsb + 0.07565561559229919*epsc + 0.5636673323593394*chi2sone*epsc + -0.5636673323593394*chi3spone*epsc + 0.15131123118459838*epsb*etasone + 0.15131123118459838*epsc*etasone + 0.12264588647823625*epsb*etaspone + 0.12264588647823625*epsc*etaspone + 0.14091683308983485*xispone + -0.01031561249540689*as*xispone + -0.061322943239118126*epsb*xispone + 0.061322943239118126*epsc*xispone + 0.12264588647823625*epsb*etasone*xispone + 0.12264588647823625*epsc*etasone*xispone + 0.17385254509323594*l2sone*pow(epsc,2) + 0.14091683308983485*l2spone*pow(epsc,2) + 0.14091683308983485*l3sone*pow(epsc,2) + 0.07565561559229919*l5sone*pow(epsc,2) + 0.061322943239118126*l5spone*pow(epsc,2) + -0.21263417442371652*l6sone*pow(epsc,2) + -0.12264588647823625*l6spone*pow(epsc,2) + 0.14091683308983485*l2sone*xispone*pow(epsc,2) + 0.061322943239118126*l5sone*xispone*pow(epsc,2) + -0.12264588647823625*l6sone*xispone*pow(epsc,2);
        }

        double P1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.6199338182961519 + -0.37274031827421483*as + -0.269777325517593*epsb + -6.69061610770223*chi2sone*epsb + -4.509338658874715*chi2spone*epsb + 20.071848323106686*chi3spone*epsb + 6.7640079883120725*chi3sppone*epsb + 0.269777325517593*epsc + 6.690616107702229*chi2sone*epsc + 4.509338658874715*chi2spone*epsc + -6.69061610770223*chi3spone*epsc + -2.2546693294373576*chi3sppone*epsc + 0.539554651035186*epsb*etasone + 0.539554651035186*epsc*etasone + 1.45578162243326*epsb*etaspone + 1.4557816224332596*epsc*etaspone + 0.49058354591294506*epsb*etasppone + 0.49058354591294506*epsc*etasppone + 1.6726540269255574*xispone + -0.24246612189750352*as*xispone + -0.72789081121663*epsb*xispone + -4.509338658874715*chi2sone*epsb*xispone + 13.528015976624145*chi3spone*epsb*xispone + 0.7278908112166298*epsc*xispone + 4.509338658874715*chi2sone*epsc*xispone + -4.509338658874715*chi3spone*epsc*xispone + 1.45578162243326*epsb*etasone*xispone + 1.4557816224332598*epsc*etasone*xispone + 0.9811670918258901*epsb*etaspone*xispone + 0.9811670918258901*epsc*etaspone*xispone + 0.5636673323593394*xisppone + -0.04126244998162756*as*xisppone + -0.24529177295647253*epsb*xisppone + 0.24529177295647253*epsc*xisppone + 0.49058354591294506*epsb*etasone*xisppone + 0.49058354591294506*epsc*etasone*xisppone + 0.6199338182961519*l2sone*pow(epsc,2) + 1.3908203607458876*l2spone*pow(epsc,2) + 1.6726540269255572*l3sone*pow(epsc,2) + 1.1273346647186788*l3spone*pow(epsc,2) + 0.269777325517593*l5sone*pow(epsc,2) + 0.6052449247383936*l5spone*pow(epsc,2) + -1.2674454622518159*l6sone*pow(epsc,2) + -1.7010733953897321*l6spone*pow(epsc,2) + 1.6726540269255574*l2sone*xispone*pow(epsc,2) + 1.1273346647186788*l2spone*xispone*pow(epsc,2) + 1.1273346647186788*l3sone*xispone*pow(epsc,2) + 0.72789081121663*l5sone*xispone*pow(epsc,2) + 0.49058354591294506*l5spone*xispone*pow(epsc,2) + -1.9463651683462044*l6sone*xispone*pow(epsc,2) + -0.9811670918258901*l6spone*xispone*pow(epsc,2) + 0.5636673323593394*l2sone*xisppone*pow(epsc,2) + 0.24529177295647253*l5sone*xisppone*pow(epsc,2) + -0.49058354591294506*l6sone*xisppone*pow(epsc,2);
        }
        // }}}

        // B_s^* -> D_s form factors
        // {{{
        double P2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03294794420017363 + -0.0015899447001142385*as + -0.01538157363810894*epsb + 0.01538157363810894*epsc + -0.03076314727621788*epsb*etasone + -0.03076314727621788*epsc*etasone + 0.03294794420017363*l1sone*pow(epsc,2) + -0.01538157363810894*l4sone*pow(epsc,2);
        }

        double P2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.27812193901117405 + -0.040298078998178166*as + -0.12983975750606846*epsb + 1.054334214405556*chi2sone*epsb + -1.054334214405556*chi3spone*epsb + 0.12983975750606844*epsc + -1.054334214405556*chi2sone*epsc + 3.1630026432166685*chi3spone*epsc + -0.2596795150121369*epsb*etasone + -0.2596795150121369*epsc*etasone + -0.24610517820974304*epsb*etaspone + -0.24610517820974304*epsc*etaspone + 0.263583553601389*xispone + -0.012719557600913908*as*xispone + -0.12305258910487152*epsb*xispone + 0.12305258910487152*epsc*xispone + -0.24610517820974304*epsb*etasone*xispone + -0.24610517820974304*epsc*etasone*xispone + 0.27812193901117405*l1sone*pow(epsc,2) + 0.263583553601389*l1spone*pow(epsc,2) + -0.12983975750606844*l4sone*pow(epsc,2) + -0.12305258910487152*l4spone*pow(epsc,2) + 0.263583553601389*l1sone*xispone*pow(epsc,2) + -0.12305258910487152*l4sone*xispone*pow(epsc,2);
        }

        double P2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.8002854055840344 + -0.600240705338798*as + -0.3736090125292201*epsb + 11.008570477168684*chi2sone*epsb + 8.43467371524445*chi2spone*epsb + -11.008570477168682*chi3spone*epsb + -4.217336857622223*chi3sppone*epsb + 0.37360901252921974*epsc + -11.008570477168682*chi2sone*epsc + -8.434673715244447*chi2spone*epsc + 33.02571143150605*chi3spone*epsc + 12.652010572866674*chi3sppone*epsc + -0.7472180250584398*epsb*etasone + -0.7472180250584395*epsc*etasone + -2.5696464765165814*epsb*etaspone + -2.569646476516581*epsc*etaspone + -0.9844207128389721*epsb*etasppone + -0.9844207128389721*epsc*etasppone + 2.7521426192921705*xispone + -0.3478237471872532*as*xispone + -1.2848232382582907*epsb*xispone + 8.43467371524445*chi2sone*epsb*xispone + -8.434673715244447*chi3spone*epsb*xispone + 1.2848232382582905*epsc*xispone + -8.434673715244447*chi2sone*epsc*xispone + 25.304021145733348*chi3spone*epsc*xispone + -2.5696464765165814*epsb*etasone*xispone + -2.569646476516581*epsc*etasone*xispone + -1.9688414256779443*epsb*etaspone*xispone + -1.9688414256779443*epsc*etaspone*xispone + 1.0543342144055559*xisppone + -0.050878230403655675*as*xisppone + -0.4922103564194861*epsb*xisppone + 0.4922103564194861*epsc*xisppone + -0.9844207128389721*epsb*etasone*xisppone + -0.9844207128389721*epsc*etasone*xisppone + 0.8002854055840344*l1sone*pow(epsc,2) + 2.224975512089393*l1spone*pow(epsc,2) + -0.37360901252921974*l4sone*pow(epsc,2) + -1.0387180600485475*l4spone*pow(epsc,2) + 2.7521426192921705*l1sone*xispone*pow(epsc,2) + 2.1086684288111117*l1spone*xispone*pow(epsc,2) + -1.2848232382582905*l4sone*xispone*pow(epsc,2) + -0.9844207128389721*l4spone*xispone*pow(epsc,2) + 1.0543342144055559*l1sone*xisppone*pow(epsc,2) + -0.4922103564194861*l4sone*xisppone*pow(epsc,2);
        }

        double V5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.005283133663761428 + 0.004846601268945023*as + 0.005283133663761428*epsb + 0.005283133663761428*epsc + -0.010566267327522857*epsc*etasone + 0.005283133663761428*l1sone*pow(epsc,2) + -0.005283133663761428*l4sone*pow(epsc,2);
        }

        double V5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.046629963052817115 + 0.05107859492650357*as + 0.046629963052817115*epsb + -0.1690602772403657*chi3spone*epsb + 0.046629963052817115*epsc + -0.1690602772403657*chi2sone*epsc + 0.5071808317210972*chi3spone*epsc + -0.09325992610563423*epsc*etasone + -0.08453013862018285*epsc*etaspone + 0.04226506931009143*xispone + 0.038772810151560186*as*xispone + 0.04226506931009143*epsb*xispone + 0.04226506931009143*epsc*xispone + -0.08453013862018285*epsc*etasone*xispone + 0.046629963052817115*l1sone*pow(epsc,2) + 0.04226506931009143*l1spone*pow(epsc,2) + -0.046629963052817115*l4sone*pow(epsc,2) + -0.04226506931009143*l4spone*pow(epsc,2) + 0.04226506931009143*l1sone*xispone*pow(epsc,2) + -0.04226506931009143*l4sone*xispone*pow(epsc,2);
        }

        double V5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.13687352617980045 + 0.15146661888422702*as + 0.13687352617980045*epsb + -1.8302793721708792*chi3spone*epsb + -0.6762411089614629*chi3sppone*epsb + 0.13687352617980045*epsc + -1.8302793721708792*chi2sone*epsc + -1.3524822179229259*chi2spone*epsc + 5.490838116512637*chi3spone*epsc + 2.0287233268843887*chi3sppone*epsc + -0.2737470523596009*epsc*etasone + -0.9151396860854396*epsc*etaspone + -0.33812055448073147*epsc*etasppone + 0.4575698430427198*xispone + 0.4861743797151489*as*xispone + 0.4575698430427198*epsb*xispone + -1.3524822179229259*chi3spone*epsb*xispone + 0.4575698430427198*epsc*xispone + -1.3524822179229259*chi2sone*epsc*xispone + 4.0574466537687774*chi3spone*epsc*xispone + -0.9151396860854396*epsc*etasone*xispone + -0.6762411089614629*epsc*etaspone*xispone + 0.16906027724036574*xisppone + 0.15509124060624077*as*xisppone + 0.16906027724036574*epsb*xisppone + 0.16906027724036574*epsc*xisppone + -0.33812055448073147*epsc*etasone*xisppone + 0.13687352617980045*l1sone*pow(epsc,2) + 0.3730397044225369*l1spone*pow(epsc,2) + -0.13687352617980045*l4sone*pow(epsc,2) + -0.3730397044225369*l4spone*pow(epsc,2) + 0.4575698430427198*l1sone*xispone*pow(epsc,2) + 0.33812055448073147*l1spone*xispone*pow(epsc,2) + -0.4575698430427198*l4sone*xispone*pow(epsc,2) + -0.33812055448073147*l4spone*xispone*pow(epsc,2) + 0.16906027724036574*l1sone*xisppone*pow(epsc,2) + -0.16906027724036574*l4sone*xisppone*pow(epsc,2);
        }

        double A2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.00503545648285162 + -0.0020945529302699613*as + 0.00503545648285162*l1sone*pow(epsc,2);
        }

        double A2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.05481409325271715 + 0.0030156952718291374*as + 0.02014182593140648*epsb + -0.16113460745125183*chi3spone*epsb + 0.02014182593140648*epsc + -0.16113460745125183*chi2sone*epsc + 0.4834038223537555*chi3spone*epsc + -0.04028365186281296*epsc*etasone + 0.04028365186281296*xispone + -0.016756423442159687*as*xispone + 0.05481409325271715*l1sone*pow(epsc,2) + 0.04028365186281296*l1spone*pow(epsc,2) + -0.02014182593140648*l4sone*pow(epsc,2) + 0.04028365186281296*l1sone*xispone*pow(epsc,2);
        }

        double A2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.24131157329136588 + 0.11401702767769033*as + 0.17897272114805565*epsb + -2.076320198989453*chi3spone*epsb + -0.6445384298050074*chi3sppone*epsb + 0.17897272114805565*epsc + -2.076320198989453*chi2sone*epsc + -1.2890768596100148*chi2spone*epsc + 6.228960596968358*chi3spone*epsc + 1.9336152894150223*chi3sppone*epsc + -0.3579454422961113*epsc*etasone + -0.3222692149025037*epsc*etaspone + 0.5190800497473632*xispone + -0.009387284709686276*as*xispone + 0.16113460745125185*epsb*xispone + -1.2890768596100148*chi3spone*epsb*xispone + 0.16113460745125185*epsc*xispone + -1.2890768596100148*chi2sone*epsc*xispone + 3.8672305788300445*chi3spone*epsc*xispone + -0.3222692149025037*epsc*etasone*xispone + 0.16113460745125185*xisppone + -0.06702569376863876*as*xisppone + 0.24131157329136588*l1sone*pow(epsc,2) + 0.43851274602173723*l1spone*pow(epsc,2) + -0.17897272114805565*l4sone*pow(epsc,2) + -0.16113460745125185*l4spone*pow(epsc,2) + 0.5190800497473632*l1sone*xispone*pow(epsc,2) + 0.3222692149025037*l1spone*xispone*pow(epsc,2) + -0.16113460745125185*l4sone*xispone*pow(epsc,2) + 0.16113460745125185*l1sone*xisppone*pow(epsc,2);
        }

        double A6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0035285573463436095 + -0.0014677418332536486*as + 0.0035285573463436095*l1sone*pow(epsc,2);
        }

        double A6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.031570021713464555 + 0.005931223319652669*as + -0.03023324226499239*epsb + 0.1129138350829955*chi2sone*epsb + -0.1129138350829955*chi3spone*epsb + 0.03023324226499239*epsc + -0.1129138350829955*chi2sone*epsc + 0.33874150524898655*chi3spone*epsc + -0.06046648452998478*epsb*etasone + -0.06046648452998478*epsc*etasone + 0.028228458770748876*xispone + -0.01174193466602919*as*xispone + 0.031570021713464555*l1sone*pow(epsc,2) + 0.028228458770748876*l1spone*pow(epsc,2) + -0.03023324226499239*l4sone*pow(epsc,2) + 0.028228458770748876*l1sone*xispone*pow(epsc,2);
        }

        double A6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1012710007502295 + 0.07834959986071738*as + -0.21003049797287487*epsb + 1.2360683649968567*chi2sone*epsb + 0.903310680663964*chi2spone*epsb + -1.2360683649968567*chi3spone*epsb + -0.451655340331982*chi3sppone*epsb + 0.21003049797287487*epsc + -1.2360683649968567*chi2sone*epsc + -0.903310680663964*chi2spone*epsc + 3.708205094990571*chi3spone*epsc + 1.354966020995946*chi3sppone*epsc + -0.42006099594574975*epsb*etasone + -0.42006099594574975*epsc*etasone + -0.4837318762398783*epsb*etaspone + -0.4837318762398782*epsc*etaspone + 0.30901709124921417*xispone + 0.02396591722516295*as*xispone + -0.24186593811993914*epsb*xispone + 0.903310680663964*chi2sone*epsb*xispone + -0.903310680663964*chi3spone*epsb*xispone + 0.2418659381199391*epsc*xispone + -0.903310680663964*chi2sone*epsc*xispone + 2.709932041991892*chi3spone*epsc*xispone + -0.4837318762398783*epsb*etasone*xispone + -0.4837318762398782*epsc*etasone*xispone + 0.1129138350829955*xisppone + -0.046967738664116764*as*xisppone + 0.1012710007502295*l1sone*pow(epsc,2) + 0.25256017370771644*l1spone*pow(epsc,2) + -0.21003049797287487*l4sone*pow(epsc,2) + -0.2418659381199391*l4spone*pow(epsc,2) + 0.30901709124921417*l1sone*xispone*pow(epsc,2) + 0.225827670165991*l1spone*xispone*pow(epsc,2) + -0.2418659381199391*l4sone*xispone*pow(epsc,2) + 0.1129138350829955*l1sone*xisppone*pow(epsc,2);
        }
        // }}}

        // B_s^* -> D_s^* form factors
        // {{{
        double S2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.02706009741680484 + 0.006784122900199494*as + 0.02706009741680484*l2sone*pow(epsc,2);
        }

        double S2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.2189743421429717 + 0.14491175989743943*as + -0.24666697664657322*epsb + -0.8659231173377548*chi3spone*epsb + 0.24666697664657322*epsc + -0.8659231173377548*chi3spone*epsc + 0.2164807793344387*xispone + 0.05427298320159596*as*xispone + 0.2189743421429717*l2sone*pow(epsc,2) + 0.2164807793344387*l2spone*pow(epsc,2) + -0.24666697664657322*l5sone*pow(epsc,2) + 0.2164807793344387*l2sone*xispone*pow(epsc,2);
        }

        double S2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.6697186753961397 + 0.46925290063727393*as + -1.502731992354025*epsb + -8.739025183250604*chi3spone*epsb + -3.4636924693510194*chi3sppone*epsb + 1.502731992354025*epsc + -8.739025183250604*chi3spone*epsc + -3.4636924693510194*chi3sppone*epsc + 2.184756295812651*xispone + 1.2678400455827072*as*xispone + -1.9733358131725856*epsb*xispone + -6.927384938702039*chi3spone*epsb*xispone + 1.9733358131725856*epsc*xispone + -6.927384938702039*chi3spone*epsc*xispone + 0.8659231173377548*xisppone + 0.2170919328063838*as*xisppone + 0.6697186753961397*l2sone*pow(epsc,2) + 1.7517947371437732*l2spone*pow(epsc,2) + -1.502731992354025*l5sone*pow(epsc,2) + -1.9733358131725856*l5spone*pow(epsc,2) + 2.184756295812651*l2sone*xispone*pow(epsc,2) + 1.7318462346755097*l2spone*xispone*pow(epsc,2) + -1.9733358131725856*l5sone*xispone*pow(epsc,2) + 0.8659231173377548*l2sone*xisppone*pow(epsc,2);
        }

        double S3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.01913437838299128 + 0.00479709930713401*as + 0.01913437838299128*l2sone*pow(epsc,2);
        }

        double S3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.15483824223515844 + 0.1024680880971562*as + -0.17441989188157567*epsb + 0.6123001082557211*chi2sone*epsb + -0.6123001082557209*chi3spone*epsb + 0.17441989188157567*epsc + 0.6123001082557211*chi2sone*epsc + -0.6123001082557209*chi3spone*epsc + -0.34883978376315133*epsb*etasone + 0.34883978376315133*epsc*etasone + 0.15307502706393022*xispone + 0.03837679445707208*as*xispone + 0.15483824223515844*l2sone*pow(epsc,2) + 0.15307502706393022*l2spone*pow(epsc,2) + 0.15307502706393028*l3sone*pow(epsc,2) + 0.17441989188157567*l5sone*pow(epsc,2) + -0.34883978376315133*l6sone*pow(epsc,2) + 0.15307502706393022*l2sone*xispone*pow(epsc,2);
        }

        double S3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.4735626168598823 + 0.3318119081320737*as + -1.0625919820995022*epsb + 6.179423968036516*chi2sone*epsb + 4.898400866045769*chi2spone*epsb + -6.179423968036511*chi3spone*epsb + -2.449200433022884*chi3sppone*epsb + 1.0625919820995022*epsc + 6.179423968036516*chi2sone*epsc + 4.898400866045769*chi2spone*epsc + -6.179423968036511*chi3spone*epsc + -2.449200433022884*chi3sppone*epsc + -2.1251839641990045*epsb*etasone + 2.1251839641990045*epsc*etasone + -2.7907182701052107*epsb*etaspone + 2.7907182701052107*epsc*etaspone + 1.5448559920091278*xispone + 0.8964982936913937*as*xispone + -1.3953591350526053*epsb*xispone + 4.898400866045769*chi2sone*epsb*xispone + -4.898400866045768*chi3spone*epsb*xispone + 1.3953591350526053*epsc*xispone + 4.898400866045769*chi2sone*epsc*xispone + -4.898400866045768*chi3spone*epsc*xispone + -2.7907182701052107*epsb*etasone*xispone + 2.7907182701052107*epsc*etasone*xispone + 0.612300108255721*xisppone + 0.15350717782828832*as*xisppone + 0.4735626168598823*l2sone*pow(epsc,2) + 1.2387059378812677*l2spone*pow(epsc,2) + 1.544855992009129*l3sone*pow(epsc,2) + 1.2246002165114422*l3spone*pow(epsc,2) + 1.0625919820995022*l5sone*pow(epsc,2) + 1.3953591350526053*l5spone*pow(epsc,2) + -3.52054309925161*l6sone*pow(epsc,2) + -2.7907182701052107*l6spone*pow(epsc,2) + 1.5448559920091278*l2sone*xispone*pow(epsc,2) + 1.224600216511442*l2spone*xispone*pow(epsc,2) + 1.2246002165114422*l3sone*xispone*pow(epsc,2) + 1.3953591350526053*l5sone*xispone*pow(epsc,2) + -2.7907182701052107*l6sone*xispone*pow(epsc,2) + 0.612300108255721*l2sone*xisppone*pow(epsc,2);
        }

        double P3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.02181039634593078 + -0.0015340565139036527*as + -0.009570660132031263*epsb + 0.009570660132031263*epsc + 0.02181039634593078*l2sone*pow(epsc,2) + -0.009570660132031263*l5sone*pow(epsc,2);
        }

        double P3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.22407038485915534 + -0.03424640011175238*as + -0.09832473766761804*epsb + -0.697932683069785*chi3spone*epsb + 0.09832473766761804*epsc + -0.697932683069785*chi3spone*epsc + 0.17448317076744624*xispone + -0.012272452111229222*as*xispone + -0.07656528105625009*epsb*xispone + 0.07656528105625009*epsc*xispone + 0.22407038485915534*l2sone*pow(epsc,2) + 0.17448317076744624*l2spone*pow(epsc,2) + -0.09832473766761804*l5sone*pow(epsc,2) + -0.07656528105625009*l5spone*pow(epsc,2) + 0.17448317076744624*l2sone*xispone*pow(epsc,2) + -0.07656528105625009*l5sone*xispone*pow(epsc,2);
        }

        double P3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.8334754182339632 + -0.4704868585026707*as + -0.36573888111887276*epsb + -8.56611768163254*chi3spone*epsb + -2.79173073227914*chi3sppone*epsb + 0.36573888111887276*epsc + -8.56611768163254*chi3spone*epsc + -2.79173073227914*chi3sppone*epsc + 2.141529420408135*xispone + -0.2985161051164776*as*xispone + -0.9397284634534445*epsb*xispone + -5.58346146455828*chi3spone*epsb*xispone + 0.9397284634534445*epsc*xispone + -5.58346146455828*chi3spone*epsc*xispone + 0.697932683069785*xisppone + -0.04908980844491689*as*xisppone + -0.30626112422500035*epsb*xisppone + 0.30626112422500035*epsc*xisppone + 0.8334754182339632*l2sone*pow(epsc,2) + 1.7925630788732423*l2spone*pow(epsc,2) + -0.36573888111887276*l5sone*pow(epsc,2) + -0.7865979013409443*l5spone*pow(epsc,2) + 2.141529420408135*l2sone*xispone*pow(epsc,2) + 1.39586536613957*l2spone*xispone*pow(epsc,2) + -0.9397284634534445*l5sone*xispone*pow(epsc,2) + -0.6125222484500007*l5spone*xispone*pow(epsc,2) + 0.697932683069785*l2sone*xisppone*pow(epsc,2) + -0.30626112422500035*l5sone*xisppone*pow(epsc,2);
        }

        double V2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.004609867550019866 + 0.0016868153974270927*as + -0.0020228644576260183*epsb + 0.0020228644576260183*epsc + (0.004609867550019866*l2sone - 0.0020228644576260183*l5sone)*pow(epsc,2);
        }

        double V2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04453083033743383 + -0.01954065555696522*epsb + -0.1475157616006357*chi3spone*epsb + 0.01954065555696522*epsc + -0.1475157616006357*chi3spone*epsc + as*(0.020136736299265223 + 0.013494523179416741*xispone) + 0.036878940400158926*xispone + -0.01618291566100815*epsb*xispone + 0.01618291566100815*epsc*xispone + (0.036878940400158926*l2spone - 0.01954065555696522*l5sone - 0.01618291566100815*l5spone + l2sone*(0.04453083033743383 + 0.036878940400158926*xispone) - 0.01618291566100815*l5sone*xispone)*pow(epsc,2);
        }

        double V2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.14479527883566448 + -0.06353788259869902*epsb + -1.720018093999154*chi3spone*epsb + -0.5900630464025428*chi3sppone*epsb + 0.06353788259869902*epsc + -1.720018093999154*chi3spone*epsc + -0.5900630464025428*chi3sppone*epsc + 0.4300045234997885*xispone + -0.18869107577773805*epsb*xispone + -1.1801260928050856*chi3spone*epsb*xispone + 0.18869107577773805*epsc*xispone + -1.1801260928050856*chi3spone*epsc*xispone + as*(0.04740742180967111 + 0.18808293675295523*xispone + 0.053978092717666966*xisppone) + 0.1475157616006357*xisppone + -0.06473166264403259*epsb*xisppone + 0.06473166264403259*epsc*xisppone + (0.3562466426994706*l2spone - 0.06353788259869902*l5sone - 0.15632524445572177*l5spone + 0.2950315232012714*l2spone*xispone - 0.18869107577773805*l5sone*xispone - 0.12946332528806517*l5spone*xispone + l2sone*(0.14479527883566448 + 0.4300045234997885*xispone + 0.1475157616006357*xisppone) - 0.06473166264403259*l5sone*xisppone)*pow(epsc,2);
        }

        double V3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.0032596686049908637 + 0.0011927586061305786*as + -0.0014303811754086053*epsb + 0.0014303811754086053*epsc + -0.0028607623508172106*epsb*etasone + 0.0028607623508172106*epsc*etasone + (0.0032596686049908637*l2sone + 0.0014303811754086053*l5sone - 0.0028607623508172106*l6sone)*pow(epsc,2);
        }

        double V3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03148805210346711 + -0.0138173300531607*epsb + 0.10430939535970762*chi2sone*epsb + -0.10430939535970764*chi3spone*epsb + 0.013817330053160699*epsc + 0.10430939535970761*chi2sone*epsc + -0.10430939535970764*chi3spone*epsc + -0.0276346601063214*epsb*etasone + 0.0276346601063214*epsc*etasone + -0.022886098806537684*epsb*etaspone + 0.022886098806537684*epsc*etaspone + as*(0.014238822788175745 + 0.009542068849044629*xispone) + 0.02607734883992691*xispone + -0.011443049403268842*epsb*xispone + 0.011443049403268842*epsc*xispone + -0.022886098806537684*epsb*etasone*xispone + 0.022886098806537684*epsc*etasone*xispone + (0.02607734883992691*l2spone + 0.026077348839926903*l3sone + 0.0138173300531607*l5sone + 0.011443049403268842*l5spone - 0.03907770950959024*l6sone - 0.022886098806537684*l6spone + l2sone*(0.03148805210346711 + 0.02607734883992691*xispone) + 0.011443049403268842*l5sone*xispone - 0.022886098806537684*l6sone*xispone)*pow(epsc,2);
        }

        double V3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.10238572354849537 + -0.044928067647774815*epsb + 1.2162364580303622*chi2sone*epsb + 0.8344751628776611*chi2spone*epsb + -1.2162364580303626*chi3spone*epsb + -0.41723758143883055*chi3sppone*epsb + 0.04492806764777479*epsc + 1.2162364580303622*chi2sone*epsc + 0.834475162877661*chi2spone*epsc + -1.2162364580303626*chi3spone*epsc + -0.41723758143883055*chi3sppone*epsc + -0.08985613529554964*epsb*etasone + 0.08985613529554964*epsc*etasone + -0.2668494784636466*epsb*etaspone + 0.2668494784636466*epsc*etaspone + -0.09154439522615074*epsb*etasppone + 0.09154439522615074*epsc*etasppone + 0.30405911450759066*xispone + -0.1334247392318233*epsb*xispone + 0.8344751628776611*chi2sone*epsb*xispone + -0.8344751628776611*chi3spone*epsb*xispone + 0.13342473923182327*epsc*xispone + 0.834475162877661*chi2sone*epsc*xispone + -0.8344751628776611*chi3spone*epsc*xispone + -0.2668494784636466*epsb*etasone*xispone + 0.2668494784636466*epsc*etasone*xispone + -0.18308879045230148*epsb*etaspone*xispone + 0.18308879045230148*epsc*etaspone*xispone + as*(0.03352210944018954 + 0.13299472000349521*xispone + 0.03816827539617852*xisppone) + 0.10430939535970764*xisppone + -0.04577219761307537*epsb*xisppone + 0.04577219761307537*epsc*xisppone + -0.09154439522615074*epsb*etasone*xisppone + 0.09154439522615074*epsc*etasone*xisppone + (0.2519044168277369*l2spone + 0.30405911450759054*l3sone + 0.20861879071941525*l3spone + 0.04492806764777482*l5sone + 0.11053864042528559*l5spone - 0.22328087452737289*l6sone - 0.3126216760767219*l6spone + 0.20861879071941528*l2spone*xispone + 0.20861879071941525*l3sone*xispone + 0.1334247392318233*l5sone*xispone + 0.09154439522615074*l5spone*xispone - 0.35839387368979736*l6sone*xispone - 0.18308879045230148*l6spone*xispone + l2sone*(0.10238572354849537 + 0.30405911450759066*xispone + 0.10430939535970764*xisppone) + 0.04577219761307537*l5sone*xisppone - 0.09154439522615074*l6sone*xisppone)*pow(epsc,2);
        }

        double V6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.00414233004658001 + 0.003800059460515357*as + 0.00414233004658001*epsb + 0.00414233004658001*epsc + 0.00828466009316002*epsc*etasone + 0.00414233004658001*l2sone*pow(epsc,2) + 0.00414233004658001*l5sone*pow(epsc,2) + -0.00828466009316002*l6sone*pow(epsc,2);
        }

        double V6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04001446777036228 + 0.04321714176889576*as + 0.04001446777036228*epsb + -0.13255456149056033*chi3spone*epsb + 0.04001446777036228*epsc + 0.13255456149056033*chi2sone*epsc + -0.13255456149056033*chi3spone*epsc + 0.08002893554072456*epsc*etasone + 0.06627728074528016*epsc*etaspone + 0.03313864037264008*xispone + 0.030400475684122855*as*xispone + 0.03313864037264008*epsb*xispone + 0.03313864037264008*epsc*xispone + 0.06627728074528016*epsc*etasone*xispone + 0.04001446777036228*l2sone*pow(epsc,2) + 0.03313864037264008*l2spone*pow(epsc,2) + 0.03313864037264008*l3sone*pow(epsc,2) + 0.04001446777036228*l5sone*pow(epsc,2) + 0.03313864037264008*l5spone*pow(epsc,2) + -0.11316757591336464*l6sone*pow(epsc,2) + -0.06627728074528016*l6spone*pow(epsc,2) + 0.03313864037264008*l2sone*xispone*pow(epsc,2) + 0.03313864037264008*l5sone*xispone*pow(epsc,2) + -0.06627728074528016*l6sone*xispone*pow(epsc,2);
        }

        double V6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1301099928828365 + 0.14509523330458124*as + 0.1301099928828365*epsb + -1.5455720916327136*chi3spone*epsb + -0.5302182459622413*chi3sppone*epsb + 0.1301099928828365*epsc + 1.5455720916327136*chi2sone*epsc + 1.0604364919244826*chi2spone*epsc + -1.5455720916327136*chi3spone*epsc + -0.5302182459622413*chi3sppone*epsc + 0.260219985765673*epsc*etasone + 0.7727860458163568*epsc*etaspone + 0.26510912298112066*epsc*etasppone + 0.3863930229081784*xispone + 0.4065380855194118*as*xispone + 0.3863930229081784*epsb*xispone + -1.0604364919244826*chi3spone*epsb*xispone + 0.3863930229081784*epsc*xispone + 1.0604364919244826*chi2sone*epsc*xispone + -1.0604364919244826*chi3spone*epsc*xispone + 0.7727860458163568*epsc*etasone*xispone + 0.5302182459622413*epsc*etaspone*xispone + 0.13255456149056033*xisppone + 0.12160190273649144*as*xisppone + 0.13255456149056033*epsb*xisppone + 0.13255456149056033*epsc*xisppone + 0.26510912298112066*epsc*etasone*xisppone + 0.1301099928828365*l2sone*pow(epsc,2) + 0.32011574216289823*l2spone*pow(epsc,2) + 0.3863930229081784*l3sone*pow(epsc,2) + 0.26510912298112066*l3spone*pow(epsc,2) + 0.1301099928828365*l5sone*pow(epsc,2) + 0.32011574216289823*l5spone*pow(epsc,2) + -0.6466130086738514*l6sone*pow(epsc,2) + -0.9053406073069172*l6spone*pow(epsc,2) + 0.3863930229081784*l2sone*xispone*pow(epsc,2) + 0.26510912298112066*l2spone*xispone*pow(epsc,2) + 0.26510912298112066*l3sone*xispone*pow(epsc,2) + 0.3863930229081784*l5sone*xispone*pow(epsc,2) + 0.26510912298112066*l5spone*xispone*pow(epsc,2) + -1.0378951687974776*l6sone*xispone*pow(epsc,2) + -0.5302182459622413*l6spone*xispone*pow(epsc,2) + 0.13255456149056033*l2sone*xisppone*pow(epsc,2) + 0.13255456149056033*l5sone*xisppone*pow(epsc,2) + -0.26510912298112066*l6sone*xisppone*pow(epsc,2);
        }

        double V7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.00414233004658001 + 0.003800059460515357*as + 0.00414233004658001*epsb + 0.00414233004658001*epsc + 0.00828466009316002*epsb*etasone + 0.00414233004658001*l2sone*pow(epsc,2) + -0.00414233004658001*l5sone*pow(epsc,2);
        }

        double V7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.04001446777036228 + 0.04321714176889576*as + 0.04001446777036228*epsb + 0.13255456149056033*chi2sone*epsb + -0.13255456149056033*chi3spone*epsb + 0.04001446777036228*epsc + -0.13255456149056033*chi3spone*epsc + 0.08002893554072456*epsb*etasone + 0.06627728074528016*epsb*etaspone + 0.03313864037264008*xispone + 0.030400475684122855*as*xispone + 0.03313864037264008*epsb*xispone + 0.03313864037264008*epsc*xispone + 0.06627728074528016*epsb*etasone*xispone + 0.04001446777036228*l2sone*pow(epsc,2) + 0.03313864037264008*l2spone*pow(epsc,2) + -0.04001446777036228*l5sone*pow(epsc,2) + -0.03313864037264008*l5spone*pow(epsc,2) + 0.03313864037264008*l2sone*xispone*pow(epsc,2) + -0.03313864037264008*l5sone*xispone*pow(epsc,2);
        }

        double V7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.1301099928828365 + 0.14509523330458124*as + 0.1301099928828365*epsb + 1.5455720916327136*chi2sone*epsb + 1.0604364919244826*chi2spone*epsb + -1.5455720916327136*chi3spone*epsb + -0.5302182459622413*chi3sppone*epsb + 0.1301099928828365*epsc + -1.5455720916327136*chi3spone*epsc + -0.5302182459622413*chi3sppone*epsc + 0.260219985765673*epsb*etasone + 0.7727860458163568*epsb*etaspone + 0.26510912298112066*epsb*etasppone + 0.3863930229081784*xispone + 0.4065380855194118*as*xispone + 0.3863930229081784*epsb*xispone + 1.0604364919244826*chi2sone*epsb*xispone + -1.0604364919244826*chi3spone*epsb*xispone + 0.3863930229081784*epsc*xispone + -1.0604364919244826*chi3spone*epsc*xispone + 0.7727860458163568*epsb*etasone*xispone + 0.5302182459622413*epsb*etaspone*xispone + 0.13255456149056033*xisppone + 0.12160190273649144*as*xisppone + 0.13255456149056033*epsb*xisppone + 0.13255456149056033*epsc*xisppone + 0.26510912298112066*epsb*etasone*xisppone + 0.1301099928828365*l2sone*pow(epsc,2) + 0.32011574216289823*l2spone*pow(epsc,2) + -0.1301099928828365*l5sone*pow(epsc,2) + -0.32011574216289823*l5spone*pow(epsc,2) + 0.3863930229081784*l2sone*xispone*pow(epsc,2) + 0.26510912298112066*l2spone*xispone*pow(epsc,2) + -0.3863930229081784*l5sone*xispone*pow(epsc,2) + -0.26510912298112066*l5spone*xispone*pow(epsc,2) + 0.13255456149056033*l2sone*xisppone*pow(epsc,2) + -0.13255456149056033*l5sone*xisppone*pow(epsc,2);
        }

        double A3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0032646160303035363 + -0.0013579525700927507*as + 0.0032646160303035363*l2sone*pow(epsc,2);
        }

        double A3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03958930938988248 + 0.0002697117831198897*as + 0.013058464121214145*epsb + 0.10446771296971316*chi2sone*epsb + -0.10446771296971316*chi3spone*epsb + 0.013058464121214145*epsc + -0.10446771296971316*chi3spone*epsc + 0.02611692824242829*epsb*etasone + 0.02611692824242829*xispone + -0.010863620560742008*as*xispone + 0.03958930938988248*l2sone*pow(epsc,2) + 0.02611692824242829*l2spone*pow(epsc,2) + -0.013058464121214145*l5sone*pow(epsc,2) + 0.02611692824242829*l2sone*xispone*pow(epsc,2);
        }

        double A3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.19459138473995258 + 0.07882797120349648*as + 0.1322403093171016*epsb + 1.4757933264156655*chi2sone*epsb + 0.8357417037577054*chi2spone*epsb + -1.4757933264156655*chi3spone*epsb + -0.4178708518788527*chi3sppone*epsb + 0.1322403093171016*epsc + -1.4757933264156655*chi3spone*epsc + -0.4178708518788527*chi3sppone*epsc + 0.2644806186342032*epsb*etasone + 0.20893542593942635*epsb*etaspone + 0.3689483316039164*xispone + -0.019569546856524914*as*xispone + 0.10446771296971318*epsb*xispone + 0.8357417037577054*chi2sone*epsb*xispone + -0.8357417037577054*chi3spone*epsb*xispone + 0.10446771296971318*epsc*xispone + -0.8357417037577054*chi3spone*epsc*xispone + 0.20893542593942635*epsb*etasone*xispone + 0.10446771296971318*xisppone + -0.04345448224296803*as*xisppone + 0.19459138473995258*l2sone*pow(epsc,2) + 0.3167144751190598*l2spone*pow(epsc,2) + -0.1322403093171016*l5sone*pow(epsc,2) + -0.10446771296971318*l5spone*pow(epsc,2) + 0.3689483316039164*l2sone*xispone*pow(epsc,2) + 0.20893542593942635*l2spone*xispone*pow(epsc,2) + -0.10446771296971318*l5sone*xispone*pow(epsc,2) + 0.10446771296971318*l2sone*xisppone*pow(epsc,2);
        }

        double A4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0032646160303035363 + -0.0013579525700927507*as + 0.0032646160303035363*l2sone*pow(epsc,2);
        }

        double A4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03958930938988248 + 0.0002697117831198897*as + 0.013058464121214145*epsb + -0.10446771296971316*chi3spone*epsb + 0.013058464121214145*epsc + 0.10446771296971316*chi2sone*epsc + -0.10446771296971316*chi3spone*epsc + 0.02611692824242829*epsc*etasone + 0.02611692824242829*xispone + -0.010863620560742008*as*xispone + 0.03958930938988248*l2sone*pow(epsc,2) + 0.02611692824242829*l2spone*pow(epsc,2) + 0.02611692824242829*l3sone*pow(epsc,2) + 0.013058464121214145*l5sone*pow(epsc,2) + -0.02611692824242829*l6sone*pow(epsc,2) + 0.02611692824242829*l2sone*xispone*pow(epsc,2);
        }

        double A4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.19459138473995258 + 0.07882797120349648*as + 0.1322403093171016*epsb + -1.4757933264156655*chi3spone*epsb + -0.4178708518788527*chi3sppone*epsb + 0.1322403093171016*epsc + 1.4757933264156655*chi2sone*epsc + 0.8357417037577054*chi2spone*epsc + -1.4757933264156655*chi3spone*epsc + -0.4178708518788527*chi3sppone*epsc + 0.2644806186342032*epsc*etasone + 0.20893542593942635*epsc*etaspone + 0.3689483316039164*xispone + -0.019569546856524914*as*xispone + 0.10446771296971318*epsb*xispone + -0.8357417037577054*chi3spone*epsb*xispone + 0.10446771296971318*epsc*xispone + 0.8357417037577054*chi2sone*epsc*xispone + -0.8357417037577054*chi3spone*epsc*xispone + 0.20893542593942635*epsc*etasone*xispone + 0.10446771296971318*xisppone + -0.04345448224296803*as*xisppone + 0.19459138473995258*l2sone*pow(epsc,2) + 0.3167144751190598*l2spone*pow(epsc,2) + 0.3689483316039164*l3sone*pow(epsc,2) + 0.20893542593942635*l3spone*pow(epsc,2) + 0.1322403093171016*l5sone*pow(epsc,2) + 0.10446771296971318*l5spone*pow(epsc,2) + -0.3689483316039164*l6sone*pow(epsc,2) + -0.20893542593942635*l6spone*pow(epsc,2) + 0.3689483316039164*l2sone*xispone*pow(epsc,2) + 0.20893542593942635*l2spone*xispone*pow(epsc,2) + 0.20893542593942635*l3sone*xispone*pow(epsc,2) + 0.10446771296971318*l5sone*xispone*pow(epsc,2) + -0.20893542593942635*l6sone*xispone*pow(epsc,2) + 0.10446771296971318*l2sone*xisppone*pow(epsc,2);
        }

        double A7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();

            return 0.0030181566001492937 + -0.0012554350876400194*as + 0.0030181566001492937*l2sone*pow(epsc,2);
        }

        double A7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.03072546294993736 + 0.0048262778721295235*as + -0.02751207995069488*epsb + -0.0965810112047774*chi3spone*epsb + 0.02751207995069488*epsc + -0.0965810112047774*chi3spone*epsc + 0.02414525280119435*xispone + -0.010043480701120155*as*xispone + 0.03072546294993736*l2sone*pow(epsc,2) + 0.02414525280119435*l2spone*pow(epsc,2) + -0.02751207995069488*l5sone*pow(epsc,2) + 0.02414525280119435*l2sone*xispone*pow(epsc,2);
        }

        double A7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();

            return 0.11437742691975071 + 0.08810044101490708*as + -0.2250545454775743*epsb + -1.1763768368075502*chi3spone*epsb + -0.3863240448191096*chi3sppone*epsb + 0.2250545454775743*epsc + -1.1763768368075502*chi3spone*epsc + -0.3863240448191096*chi3sppone*epsc + 0.29409420920188756*xispone + 0.01852326157479586*as*xispone + -0.22009663960555897*epsb*xispone + -0.7726480896382192*chi3spone*epsb*xispone + 0.22009663960555897*epsc*xispone + -0.7726480896382192*chi3spone*epsc*xispone + 0.0965810112047774*xisppone + -0.040173922804480615*as*xisppone + 0.11437742691975071*l2sone*pow(epsc,2) + 0.24580370359949882*l2spone*pow(epsc,2) + -0.2250545454775743*l5sone*pow(epsc,2) + -0.22009663960555897*l5spone*pow(epsc,2) + 0.29409420920188756*l2sone*xispone*pow(epsc,2) + 0.1931620224095548*l2spone*xispone*pow(epsc,2) + -0.22009663960555897*l5sone*xispone*pow(epsc,2) + 0.0965810112047774*l2sone*xisppone*pow(epsc,2);
        }
        // }}}
    };

    const std::vector<OptionSpecification>
    Implementation<BGLCoefficients>::options
    {
        { "SU3F-limit-sslp", { "0", "1" }, "0" }
    };

    BGLCoefficients::BGLCoefficients(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BGLCoefficients>(new Implementation<BGLCoefficients>(p, o, *this))
    {
    }

    BGLCoefficients::~BGLCoefficients() = default;

    const std::set<ReferenceName>
    BGLCoefficients::references
    {
    };


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

    std::vector<OptionSpecification>::const_iterator
    BGLCoefficients::begin_options()
    {
        return Implementation<BGLCoefficients>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BGLCoefficients::end_options()
    {
        return Implementation<BGLCoefficients>::options.cend();
    }

    template <> struct Implementation<HQETUnitarityBounds>
    {
        // option to determine if we use z^3 terms in the leading-power IW function
        SwitchOption opt_zorder_bound;
        unsigned zorder_bound;

        // number of light flavor multiplets
        UsedParameter nf;
        UsedParameter ns;

        std::shared_ptr<BGLCoefficients> bgl;

        static const std::vector<OptionSpecification> options;

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

            return result;
        }

        double bound_0p_prior() const
        {
            const double value = bound_0p();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 0^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
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

            return result;
        }

        double bound_0m_prior() const
        {
            const double value = bound_0m();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 0^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // using the same relative uncertainty as for 0^+, cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
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

            return result;
        }

        double bound_1p_prior() const
        {
            const double value = bound_1p();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 1^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0093549; // cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
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

            return result;
        }

        double bound_1m_prior() const
        {
            const double value = bound_1m();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 1^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0093549; // same relative uncertainty as for 1^-, cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }
        // }}}
    };

    const std::vector<OptionSpecification>
    Implementation<HQETUnitarityBounds>::options
    {
        { "SU3F-limit-sslp", { "0", "1" }, "0" },
        { "z-order-bound",   { "1", "2" }, "2" }
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

    double
    HQETUnitarityBounds::bound_0p_prior() const
    {
        return _imp->bound_0p_prior();
    }

    double
    HQETUnitarityBounds::bound_0m_prior() const
    {
        return _imp->bound_0m_prior();
    }

    double
    HQETUnitarityBounds::bound_1p_prior() const
    {
        return _imp->bound_1p_prior();
    }

    double
    HQETUnitarityBounds::bound_1m_prior() const
    {
        return _imp->bound_1m_prior();
    }

    const std::set<ReferenceName>
    HQETUnitarityBounds::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    HQETUnitarityBounds::begin_options()
    {
        return Implementation<HQETUnitarityBounds>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    HQETUnitarityBounds::end_options()
    {
        return Implementation<HQETUnitarityBounds>::options.cend();
    }

    template <> struct Implementation<OPEUnitarityBounds>
    {
        const double cond_qq, cond_G2;
        const double mu;   // TODO maybe declare the scale of alpha_s(mu) as a UsedParameter

        std::shared_ptr<Model> model;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            cond_qq(-0.02/12), // (p["B->D^*::<qq>@BGL1997"], u),         // [BGL1997] quark condensate
            cond_G2(0.02),     // (p["B->D^*::<alS/pi G^2>@BGL1997"], u), // [BGL1997] gluon condensate
            mu(4.2),                                       // TODO remove hard-coded numerical value
            model(Model::make("SM", p, o))
        {
            u.uses(*model);
        }

        ~Implementation() = default;


        double chi_T(const double & u) const
        {
            // [BGL:1997A] eq.(4.1) + (4.2) + (4.5, 4.8)
            // limit for u -> 1 in eq.(4.11)
            // limit for u -> 0 in eq.(4.10)
            const double aS = model->alpha_s(mu) / M_PI;
            const double u2 = power_of<2>(u);
            const double u3 = power_of<3>(u);
            const double u4 = power_of<4>(u);
            const double u5 = power_of<5>(u);
            const double u6 = power_of<6>(u);
            const double lnu = std::log(u2);

            const double m_b = model->m_b_pole();

            const double chi_pert =
                1.0 / 32.0 * (
                   (1.0 - u2) * (3. + 4. * u - 21. * u2 + 40. * u3 - 21. * u4 + 4. * u5 + 3. * u6)
                   + 12. * u3 * (2. - 3. * u + 2. * u2) * lnu
                )
                + aS / (576.0 * (1.0 - u2)) * (
                   power_of<2>(1.0 - u2) * (75. + 360. * u - 1031. * u2 + 1776. * u3 - 1031. * u4 + 360. * u5 + 75. * u6)
                   + 4.0 * u * (1.0 - u2) * (18. - 99. * u + 732. * u2 - 1010. * u3 + 732. * u4 - 99. * u5 + 18. * u6) * lnu
                   + 4.0 * u3 * (108. - 324. * u + 648. * u2 - 456. * u3 + 132. * u4 + 59. * u5 - 12. * u6 - 9. *u6 * u) * lnu * lnu
                   + 8.0 * power_of<3>(1.0 - u2) * (9. + 12. * u - 32. * u2 + 12. * u3 + 9. * u4) * gsl_sf_dilog(1.0 - u2)
                );

            const double chi_cond =
                -1.0 * cond_qq / 2.0 * (2. - 3. * u + 2. * u2)
                -1.0 * cond_G2 / (24.0 * m_b * power_of<2>(1.0 - u2) ) * (
                    (1.0 - u2) * (2. - 104. * u + 148. * u2 - 270. * u3 + 145. * u4 - 104. * u5 + 5. * u6 - 2. * u6 * u)
                   - 12.0 * u * (3. - 5. * u + 17. * u2 - 15. * u3 + 17. * u4 - 5. * u5 + 3. * u6) * lnu
                );

            //  TODO catch u->1 and u->0, but should not be needed for b->c
            return chi_pert / (power_of<2>(m_b * M_PI) * power_of<5>(1.0 - u2))
                 + chi_cond / (power_of<5>(m_b * (1.0 - u2)) );
        }

        double chi_L(const double & u) const
        {
            // [BGL1997] eq.(4.1) + (4.3) + (4.6, 4.9)
            // limit for u -> 1 in eq.(4.11)
            // limit for u -> 0 in eq.(4.10)
            const double aS = model->alpha_s(mu) / M_PI;
            const double u2 = power_of<2>(u);
            const double u3 = power_of<3>(u);
            const double u4 = power_of<4>(u);
            const double u5 = power_of<5>(u);
            const double lnu = std::log(u2);

            const double m_b = model->m_b_pole();

            const double chi_pert =
                1.0 / 8.0 * (
                   (1.0 - u2) * (1. + u + u2) * (1 - 4. * u + u2) - 6. * u3 * lnu
                )
                + aS / (48.0 * (1.0 - u2)) * (
                   power_of<2>(1.0 - u2) * (1. - 36. * u - 22. * u2 - 36. * u3 + u4)
                   - 2.0 * u * (1.0 - u2) * (9. + 4. * u + 66. * u2 + 4. * u3 + 9. * u4) * lnu
                   - 4.0 * u3 * (9. + 18. * u2 - 2. * u3 - 3. * u4 + u5) * lnu * lnu
                   + 8.0 * power_of<3>(1.0 - u2) * (1. - 3. * u + u2) * gsl_sf_dilog(1.0 - u2)
                );

            const double chi_cond =
                cond_qq + cond_G2 / (12.0 * m_b * power_of<2>(1.0 - u2) ) * (
                    (1.0 - u2) * (1. - 21. * u + 10. * u2 - 20. * u3 + u4 - u5)
                   - 3.0 * u * (3. - 2. * u + 8. * u2 - 2. * u3 + 3. * u4) * lnu
                );

            //  TODO catch u->1 and u->0, but should not be needed for b->c
            return chi_pert/(M_PI * M_PI * power_of<3>(1.0 - u2))
                 + chi_cond / (power_of<3>(m_b * (1.0 - u2)) );
        }

        double chi_1m() const
        {
            return chi_T(+1.0 * model->m_c_pole() / model->m_b_pole());
        }

        double chi_1p() const
        {
            return chi_T(-1.0 * model->m_c_pole() / model->m_b_pole());
        }

        double chi_0p() const
        {
            return chi_L(+1.0 * model->m_c_pole() / model->m_b_pole());
        }

        double chi_0m() const
        {
            return chi_L(-1.0 * model->m_c_pole() / model->m_b_pole());
        }
    };

    const std::vector<OptionSpecification>
    Implementation<OPEUnitarityBounds>::options
    {
    };

    OPEUnitarityBounds::OPEUnitarityBounds(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<OPEUnitarityBounds>(new Implementation<OPEUnitarityBounds>(p, o, *this))
    {
    }

    OPEUnitarityBounds::~OPEUnitarityBounds() = default;

    double
    OPEUnitarityBounds::bound_0p() const
    {
        return _imp->chi_0p();
    }

    double
    OPEUnitarityBounds::bound_0m() const
    {
        return _imp->chi_0m();
    }

    double
    OPEUnitarityBounds::bound_1p() const
    {
        return _imp->chi_1p();
    }

    double
    OPEUnitarityBounds::bound_1m() const
    {
        return _imp->chi_1m();
    }

    const std::set<ReferenceName>
    OPEUnitarityBounds::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    OPEUnitarityBounds::begin_options()
    {
        return Implementation<OPEUnitarityBounds>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    OPEUnitarityBounds::end_options()
    {
        return Implementation<OPEUnitarityBounds>::options.cend();
    }

    template <> struct Implementation<BGLUnitarityBounds>
    {
        // B->D^* parameters
        std::array<UsedParameter, 4> _a_g, _a_f, _a_F1, _a_F2;
        // B->D parameters
        std::array<UsedParameter, 4> _a_f_p, _a_f_0, _a_f_t;

        // option to determine if we use z^3 terms in the leading-power IW function
        SwitchOption opt_zorder_bound;
        unsigned zorder_bound;

        // number of light flavor multiplets
        UsedParameter nf;
        UsedParameter ns;

        static const std::vector<OptionSpecification> options;

        std::string _par_name_dstar(const std::string & ff_name)
        {
            return std::string("B->D^*") + std::string("::a^") + ff_name + std::string("@BGL1997");
        }

        std::string _par_name_d(const std::string & ff_name)
        {
            return std::string("B->D") + std::string("::a^") + ff_name + std::string("@BGL1997");
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            // B->D^*
            _a_g{{   UsedParameter(p[_par_name_dstar("g_0")],  u),
                     UsedParameter(p[_par_name_dstar("g_1")],  u),
                     UsedParameter(p[_par_name_dstar("g_2")],  u),
                     UsedParameter(p[_par_name_dstar("g_3")],  u) }},
            _a_f{{   UsedParameter(p[_par_name_dstar("f_0")],  u),
                     UsedParameter(p[_par_name_dstar("f_1")],  u),
                     UsedParameter(p[_par_name_dstar("f_2")],  u),
                     UsedParameter(p[_par_name_dstar("f_3")],  u) }},
            _a_F1{{  UsedParameter(p[_par_name_dstar("F1_0")], u),
                     UsedParameter(p[_par_name_dstar("F1_1")], u),
                     UsedParameter(p[_par_name_dstar("F1_2")], u),
                     UsedParameter(p[_par_name_dstar("F1_3")], u) }},
            _a_F2{{  UsedParameter(p[_par_name_dstar("F2_0")], u),
                     UsedParameter(p[_par_name_dstar("F2_1")], u),
                     UsedParameter(p[_par_name_dstar("F2_2")], u),
                     UsedParameter(p[_par_name_dstar("F2_3")], u) }},
            // B->D
            _a_f_p{{ UsedParameter(p[_par_name_d("f+_0")], u),
                     UsedParameter(p[_par_name_d("f+_1")], u),
                     UsedParameter(p[_par_name_d("f+_2")], u),
                     UsedParameter(p[_par_name_d("f+_3")], u) }},
            _a_f_0{{ UsedParameter(p[_par_name_d("f0_0")], u),
                     UsedParameter(p[_par_name_d("f0_1")], u),
                     UsedParameter(p[_par_name_d("f0_2")], u),
                     UsedParameter(p[_par_name_d("f0_3")], u) }},
            _a_f_t{{ UsedParameter(p[_par_name_d("fT_0")], u),
                     UsedParameter(p[_par_name_d("fT_1")], u),
                     UsedParameter(p[_par_name_d("fT_2")], u),
                     UsedParameter(p[_par_name_d("fT_3")], u) }},
            // further parameters
            opt_zorder_bound(o, "z-order-bound", { "1", "2" }, "2"),
            nf(p["B(*)->D(*)::n_f@BGL"], u),
            ns(p["B_s(*)->D_s(*)::n_s@BGL"], u)
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
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += pow(_a_f_0[i], 2) * nf; // to account for flavour symmetry
            }

            return result;
        }

        double bound_0p_prior() const
        {
            const double value = bound_0p();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 0^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_0m() const
        {
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += pow(_a_F2[i], 2) * nf; // to account for flavour symmetry
            }

            return result;
        }

        double bound_0m_prior() const
        {
            const double value = bound_0m();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 0^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0130561; // using the same relative uncertainty as for 0^+, cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1p() const
        {
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += pow(_a_f[i],  2) * nf; // to account for flavour symmetry
                result += pow(_a_F1[i], 2) * nf; // to account for flavour symmetry
            }

            return result;
        }

        double bound_1p_prior() const
        {
            const double value = bound_1p();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 1^+ unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0093549; // cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }

        double bound_1m() const
        {
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += pow(_a_f_p[i], 2) * nf; // to account for flavour symmetry
                result += pow(_a_g[i],   2) * nf; // to account for flavour symmetry
            }

            return result;
        }

        double bound_1m_prior() const
        {
            const double value = bound_1m();

            if (value < 0.0)
            {
                throw InternalError("Contribution to 1^- unitarity bound must be positive; found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.0093549; // same relative uncertainty as for 1^-, cf. [BG2016], eq. (2.8), p.5
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }
        // }}}
    };

    const std::vector<OptionSpecification>
    Implementation<BGLUnitarityBounds>::options
    {
        { "z-order-bound", { "1", "2" }, "2" }
    };

    BGLUnitarityBounds::BGLUnitarityBounds(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BGLUnitarityBounds>(new Implementation<BGLUnitarityBounds>(p, o, *this))
    {
    }

    BGLUnitarityBounds::~BGLUnitarityBounds() = default;

    double
    BGLUnitarityBounds::bound_0p() const
    {
        return _imp->bound_0p();
    }

    double
    BGLUnitarityBounds::bound_0m() const
    {
        return _imp->bound_0m();
    }

    double
    BGLUnitarityBounds::bound_1p() const
    {
        return _imp->bound_1p();
    }

    double
    BGLUnitarityBounds::bound_1m() const
    {
        return _imp->bound_1m();
    }

    double
    BGLUnitarityBounds::bound_0p_prior() const
    {
        return _imp->bound_0p_prior();
    }

    double
    BGLUnitarityBounds::bound_0m_prior() const
    {
        return _imp->bound_0m_prior();
    }

    double
    BGLUnitarityBounds::bound_1p_prior() const
    {
        return _imp->bound_1p_prior();
    }

    double
    BGLUnitarityBounds::bound_1m_prior() const
    {
        return _imp->bound_1m_prior();
    }

    const std::set<ReferenceName>
    BGLUnitarityBounds::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BGLUnitarityBounds::begin_options()
    {
        return Implementation<BGLUnitarityBounds>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BGLUnitarityBounds::end_options()
    {
        return Implementation<BGLUnitarityBounds>::options.cend();
    }
}
