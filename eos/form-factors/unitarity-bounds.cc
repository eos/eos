/* vim: set sw=4 sts=4 tw=140 et foldmethod=marker : */

/*
 * Copyright (c) 2019-2025 Danny van Dyk
 * Copyright (c) 2019-2024 Nico Gubernari
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
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <gsl/gsl_sf_dilog.h>

#include <cmath>

#include <iostream>

namespace eos
{
    using namespace std::literals::string_literals;

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
            opt_sslp_limit(o, "SU3F-limit-sslp"_ok, { "0", "1" }, "0"),
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
         * HQET parameters following [BLPR:2017]
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
            const double epsc2 = power_of<2>(epsc);

            return 0.008703312831206973 + 0.0032729575259871553*as + -0.0041554234514795205*epsb + 0.0041554234514795205*epsc + 0.008310846902959041*epsb*etaone + -0.008310846902959041*epsc*etaone + 0.008703312831206973*l1one*epsc2 + -0.0041554234514795205*l4one*epsc2;
        }

        double V1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06022799269853725 + 0.02982697398331807*as + -0.02875604015951835*epsb + -0.2785060105986231*chi2one*epsb + 0.8355180317958693*chi3pone*epsb + 0.02875604015951835*epsc + -0.2785060105986231*chi2one*epsc + 0.8355180317958693*chi3pone*epsc + 0.0575120803190367*epsb*etaone + -0.0575120803190367*epsc*etaone + 0.06648677522367233*epsb*etapone + -0.06648677522367233*epsc*etapone + 0.06962650264965578*xipone + 0.026183660207897242*as*xipone + -0.033243387611836164*epsb*xipone + 0.033243387611836164*epsc*xipone + 0.06648677522367233*epsb*etaone*xipone + -0.06648677522367233*epsc*etaone*xipone + 0.06022799269853725*l1one*epsc2 + 0.06962650264965578*l1pone*epsc2 + -0.02875604015951835*l4one*epsc2 + -0.033243387611836164*l4pone*epsc2 + 0.06962650264965578*l1one*xipone*epsc2 + -0.033243387611836164*l4one*xipone*epsc2;
        }

        double V1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);


            return 0.13065762648978307 + 0.017947557772235324*as + -0.06238288520246169*epsb + -2.4843077875504385*chi2one*epsb + -2.228048084788985*chi2pone*epsb + 7.4529233626513145*chi3pone*epsb + 3.342072127183477*chi3ppone*epsb + 0.06238288520246169*epsc + -2.4843077875504385*chi2one*epsc + -2.228048084788985*chi2pone*epsc + 7.4529233626513145*chi3pone*epsc + 3.342072127183477*chi3ppone*epsc + 0.12476577040492338*epsb*etaone + -0.12476577040492338*epsc*etaone + 0.5930701929996383*epsb*etapone + -0.5930701929996383*epsc*etapone + 0.2659471008946893*epsb*etappone + -0.2659471008946893*epsc*etappone + 0.6210769468876096*xipone + 0.29098311228233903*as*xipone + -0.29653509649981913*epsb*xipone + -2.228048084788985*chi2one*epsb*xipone + 6.684144254366954*chi3pone*epsb*xipone + 0.29653509649981913*epsc*xipone + -2.228048084788985*chi2one*epsc*xipone + 6.684144254366954*chi3pone*epsc*xipone + 0.5930701929996383*epsb*etaone*xipone + -0.5930701929996383*epsc*etaone*xipone + 0.5318942017893786*epsb*etapone*xipone + -0.5318942017893786*epsc*etapone*xipone + 0.2785060105986231*xippone + 0.10473464083158897*as*xippone + -0.13297355044734466*epsb*xippone + 0.13297355044734466*epsc*xippone + 0.2659471008946893*epsb*etaone*xippone + -0.2659471008946893*epsc*etaone*xippone + 0.13065762648978307*l1one*epsc2 + 0.481823941588298*l1pone*epsc2 + -0.06238288520246169*l4one*epsc2 + -0.2300483212761468*l4pone*epsc2 + 0.6210769468876096*l1one*xipone*epsc2 + 0.5570120211972462*l1pone*xipone*epsc2 + -0.29653509649981913*l4one*xipone*epsc2 + -0.2659471008946893*l4pone*xipone*epsc2 + 0.2785060105986231*l1one*xippone*epsc2 + -0.13297355044734466*l4one*xippone*epsc2;
        }

        double S1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.047334757401521106 + 0.011867083795070281*as + 0.047334757401521106*l1one*epsc2;
        }

        double S1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);


            return 0.284348884202279 + 0.21957539313187335*as + -0.3965605010079518*epsb + -1.5147122368486754*chi2one*epsb + 4.544136710546026*chi3pone*epsb + 0.3965605010079518*epsc + -1.5147122368486754*chi2one*epsc + 4.544136710546026*chi3pone*epsc + 0.7931210020159036*epsb*etaone + -0.7931210020159036*epsc*etaone + 0.37867805921216885*xipone + 0.09493667036056225*as*xipone + 0.284348884202279*l1one*epsc2 + 0.37867805921216885*l1pone*epsc2 + -0.3965605010079518*l4one*epsc2 + 0.37867805921216885*l1one*xipone*epsc2;
        }

        double S1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);


            return 0.6559703043190042 + 0.3344969447011169*as + -1.5890932982243728*epsb + -12.128588768170278*chi2one*epsb + -12.117697894789403*chi2pone*epsb + 36.38576630451084*chi3pone*epsb + 18.176546842184106*chi3ppone*epsb + 1.5890932982243728*epsc + -12.128588768170278*chi2one*epsc + -12.117697894789403*chi2pone*epsc + 36.38576630451084*chi3pone*epsc + 18.176546842184106*chi3ppone*epsc + 3.1781865964487457*epsb*etaone + -3.1781865964487457*epsc*etaone + 6.34496801612723*epsb*etapone + -6.34496801612723*epsc*etapone + 3.0321471920425696*xipone + 1.9464764857761112*as*xipone + -3.172484008063615*epsb*xipone + -12.117697894789403*chi2one*epsb*xipone + 36.35309368436821*chi3pone*epsb*xipone + 3.172484008063615*epsc*xipone + -12.117697894789403*chi2one*epsc*xipone + 36.35309368436821*chi3pone*epsc*xipone + 6.34496801612723*epsb*etaone*xipone + -6.34496801612723*epsc*etaone*xipone + 1.5147122368486754*xippone + 0.379746681442249*as*xippone + 0.6559703043190042*l1one*epsc2 + 2.2747910736182324*l1pone*epsc2 + -1.5890932982243728*l4one*epsc2 + -3.172484008063615*l4pone*epsc2 + 3.0321471920425696*l1one*xipone*epsc2 + 3.029424473697351*l1pone*xipone*epsc2 + -3.172484008063615*l4one*xipone*epsc2 + 1.5147122368486754*l1one*xippone*epsc2;
        }

        double fT_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.017287377854958173 + 0.02882855797466114*as + 0.017287377854958173*epsb + 0.017287377854958173*epsc + -0.03457475570991635*epsb*etaone + -0.03457475570991635*epsc*etaone + 0.017287377854958173*l1one*epsc2 + -0.017287377854958173*l4one*epsc2;
        }

        double fT_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.13580155950042527 + 0.2668270308398708*as + 0.13580155950042527*epsb + -0.5531960913586615*chi2one*epsb + 1.6595882740759849*chi3pone*epsb + 0.13580155950042527*epsc + -0.5531960913586615*chi2one*epsc + 1.6595882740759849*chi3pone*epsc + -0.27160311900085055*epsb*etaone + -0.27160311900085055*epsc*etaone + -0.2765980456793308*epsb*etapone + -0.2765980456793308*epsc*etapone + 0.1382990228396654*xipone + 0.23062846379728916*as*xipone + 0.1382990228396654*epsb*xipone + 0.1382990228396654*epsc*xipone + -0.2765980456793308*epsb*etaone*xipone + -0.2765980456793308*epsc*etaone*xipone + 0.13580155950042527*l1one*epsc2 + 0.1382990228396654*l1pone*epsc2 + -0.13580155950042527*l4one*epsc2 + -0.1382990228396654*l4pone*epsc2 + 0.1382990228396654*l1one*xipone*epsc2 + -0.1382990228396654*l4one*xipone*epsc2;
        }

        double fT_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.3875996470087848 + 0.8150376285991504*as + 0.3875996470087848*epsb + -5.452042086730931*chi2one*epsb + -4.425568730869292*chi2pone*epsb + 16.356126260192795*chi3pone*epsb + 6.638353096303939*chi3ppone*epsb + 0.3875996470087848*epsc + -5.452042086730931*chi2one*epsc + -4.425568730869292*chi2pone*epsc + 16.356126260192795*chi3pone*epsc + 6.638353096303939*chi3ppone*epsc + -0.7751992940175696*epsb*etaone + -0.7751992940175696*epsc*etaone + -2.7260210433654657*epsb*etapone + -2.7260210433654657*epsc*etapone + -1.106392182717323*epsb*etappone + -1.106392182717323*epsc*etappone + 1.3630105216827328*xipone + 2.595873174313545*as*xipone + 1.3630105216827328*epsb*xipone + -4.425568730869292*chi2one*epsb*xipone + 13.276706192607879*chi3pone*epsb*xipone + 1.3630105216827328*epsc*xipone + -4.425568730869292*chi2one*epsc*xipone + 13.276706192607879*chi3pone*epsc*xipone + -2.7260210433654657*epsb*etaone*xipone + -2.7260210433654657*epsc*etaone*xipone + -2.212784365434646*epsb*etapone*xipone + -2.212784365434646*epsc*etapone*xipone + 0.5531960913586615*xippone + 0.9225138551891566*as*xippone + 0.5531960913586615*epsb*xippone + 0.5531960913586615*epsc*xippone + -1.106392182717323*epsb*etaone*xippone + -1.106392182717323*epsc*etaone*xippone + 0.3875996470087848*l1one*epsc2 + 1.0864124760034022*l1pone*epsc2 + -0.3875996470087848*l4one*epsc2 + -1.0864124760034022*l4pone*epsc2 + 1.3630105216827328*l1one*xipone*epsc2 + 1.106392182717323*l1pone*xipone*epsc2 + -1.3630105216827328*l4one*xipone*epsc2 + -1.106392182717323*l4pone*xipone*epsc2 + 0.5531960913586615*l1one*xippone*epsc2 + -0.5531960913586615*l4one*xippone*epsc2;
        }
        // }}}

        // B -> D^* form factors
        // {{{
        double A1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.008363045281441904 + -0.003478702112200171*as + 0.008363045281441904*l2one*epsc2;
        }

        double A1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.0812300174749896 + 0.009087844216195124*as + 0.033452181125767616*epsb + -0.26761744900614093*chi2one*epsb + 0.8028523470184228*chi3pone*epsb + 0.033452181125767616*epsc + -0.26761744900614093*chi3pone*epsc + -0.06690436225153523*epsb*etaone + 0.06690436225153523*xipone + -0.02782961689760137*as*xipone + 0.0812300174749896*l2one*epsc2 + 0.06690436225153523*l2pone*epsc2 + -0.033452181125767616*l5one*epsc2 + 0.06690436225153523*l2one*xipone*epsc2;
        }

        double A1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.31721556277653973 + 0.1738424167676217*as + 0.25801570764842313*epsb + -3.134595457211949*chi2one*epsb + -2.1409395920491274*chi2pone*epsb + 9.403786371635848*chi3pone*epsb + 3.211409388073691*chi3ppone*epsb + 0.25801570764842313*epsc + -3.134595457211949*chi3pone*epsc + -1.0704697960245637*chi3ppone*epsc + -0.5160314152968463*epsb*etaone + -0.5352348980122819*epsb*etapone + 0.7836488643029873*xipone + 0.017043519934358188*as*xipone + 0.26761744900614093*epsb*xipone + -2.1409395920491274*chi2one*epsb*xipone + 6.422818776147382*chi3pone*epsb*xipone + 0.26761744900614093*epsc*xipone + -2.1409395920491274*chi3pone*epsc*xipone + -0.5352348980122819*epsb*etaone*xipone + 0.26761744900614093*xippone + -0.11131846759040548*as*xippone + 0.31721556277653973*l2one*epsc2 + 0.6498401397999168*l2pone*epsc2 + -0.25801570764842313*l5one*epsc2 + -0.26761744900614093*l5pone*epsc2 + 0.7836488643029873*l2one*xipone*epsc2 + 0.5352348980122819*l2pone*xipone*epsc2 + -0.26761744900614093*l5one*xipone*epsc2 + 0.26761744900614093*l2one*xippone*epsc2;
        }

        double A5_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.00560220278709002 + -0.0023302988340466895*as + 0.00560220278709002*l2one*epsc2;
        }

        double A5_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.04352386562353383 + 0.013706788896745767*as + -0.049962145016877284*epsb + -0.17927048918688063*chi2one*epsb + 0.5378114675606419*chi3pone*epsb + 0.04996214501687728*epsc + 0.17927048918688063*chi2one*epsc + -0.17927048918688063*chi3pone*epsc + 0.09992429003375457*epsb*etaone + 0.09992429003375457*epsc*etaone + 0.04481762229672016*xipone + -0.01864239067237352*as*xipone + 0.04352386562353383*l2one*epsc2 + 0.04481762229672016*l2pone*epsc2 + 0.04481762229672016*l3one*epsc2 + 0.049962145016877284*l5one*epsc2 + -0.09992429003375457*l6one*epsc2 + 0.04481762229672016*l2one*xipone*epsc2;
        }

        double A5_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1173127381910642 + 0.1179664643240928*as + -0.2882347553669587*epsb + -1.7513046783268438*chi2one*epsb + -1.434163913495045*chi2pone*epsb + 5.253914034980532*chi3pone*epsb + 2.1512458702425676*chi3ppone*epsb + 0.2882347553669586*epsc + 1.7513046783268438*chi2one*epsc + 1.434163913495045*chi2pone*epsc + -1.7513046783268438*chi3pone*epsc + -0.7170819567475225*chi3ppone*epsc + 0.5764695107339174*epsb*etaone + 0.5764695107339174*epsc*etaone + 0.7993943202700365*epsb*etapone + 0.7993943202700365*epsc*etapone + 0.43782616958171094*xipone + 0.07236952982921908*as*xipone + -0.39969716013501827*epsb*xipone + -1.434163913495045*chi2one*epsb*xipone + 4.302491740485135*chi3pone*epsb*xipone + 0.3996971601350182*epsc*xipone + 1.434163913495045*chi2one*epsc*xipone + -1.434163913495045*chi3pone*epsc*xipone + 0.7993943202700365*epsb*etaone*xipone + 0.7993943202700365*epsc*etaone*xipone + 0.17927048918688063*xippone + -0.07456956268949408*as*xippone + 0.1173127381910642*l2one*epsc2 + 0.34819092498827064*l2pone*epsc2 + 0.43782616958171094*l3one*epsc2 + 0.35854097837376125*l3pone*epsc2 + 0.2882347553669587*l5one*epsc2 + 0.39969716013501827*l5pone*epsc2 + -0.9761666708689356*l6one*epsc2 + -0.7993943202700365*l6pone*epsc2 + 0.43782616958171094*l2one*xipone*epsc2 + 0.35854097837376125*l2pone*xipone*epsc2 + 0.35854097837376125*l3one*xipone*epsc2 + 0.39969716013501827*l5one*xipone*epsc2 + -0.7993943202700365*l6one*xipone*epsc2 + 0.17927048918688063*l2one*xippone*epsc2;
        }

        double V4_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.0075025867890451845 + 0.006882664262345855*as + 0.0075025867890451845*epsb + 0.0075025867890451845*epsc + -0.015005173578090369*epsb*etaone + 0.0075025867890451845*l2one*epsc2 + -0.0075025867890451845*l5one*epsc2;
        }

        double V4_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.05963010511442657 + 0.06649203938674278*as + 0.05963010511442657*epsb + -0.2400827772494459*chi2one*epsb + 0.7202483317483377*chi3pone*epsb + 0.05963010511442657*epsc + -0.2400827772494459*chi3pone*epsc + -0.11926021022885314*epsb*etaone + -0.12004138862472295*epsb*etapone + 0.060020694312361476*xipone + 0.055061314098766835*as*xipone + 0.060020694312361476*epsb*xipone + 0.060020694312361476*epsc*xipone + -0.12004138862472295*epsb*etaone*xipone + 0.05963010511442657*l2one*epsc2 + 0.060020694312361476*l2pone*epsc2 + -0.05963010511442657*l5one*epsc2 + -0.060020694312361476*l5pone*epsc2 + 0.060020694312361476*l2one*xipone*epsc2 + -0.060020694312361476*l5one*xipone*epsc2;
        }

        double V4_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1543810939148951 + 0.1680554401420258*as + 0.1543810939148951*epsb + -2.388328918160542*chi2one*epsb + -1.9206622179955672*chi2pone*epsb + 7.164986754481626*chi3pone*epsb + 2.8809933269933508*chi3ppone*epsb + 0.1543810939148951*epsc + -2.388328918160542*chi3pone*epsc + -0.9603311089977836*chi3ppone*epsc + -0.3087621878297902*epsb*etaone + -1.194164459080271*epsb*etapone + -0.4801655544988918*epsb*etappone + 0.5970822295401355*xipone + 0.642058943291476*as*xipone + 0.5970822295401355*epsb*xipone + -1.9206622179955672*chi2one*epsb*xipone + 5.7619866539867015*chi3pone*epsb*xipone + 0.5970822295401355*epsc*xipone + -1.9206622179955672*chi3pone*epsc*xipone + -1.194164459080271*epsb*etaone*xipone + -0.9603311089977836*epsb*etapone*xipone + 0.2400827772494459*xippone + 0.22024525639506734*as*xippone + 0.2400827772494459*epsb*xippone + 0.2400827772494459*epsc*xippone + -0.4801655544988918*epsb*etaone*xippone + 0.1543810939148951*l2one*epsc2 + 0.47704084091541255*l2pone*epsc2 + -0.1543810939148951*l5one*epsc2 + -0.47704084091541255*l5pone*epsc2 + 0.5970822295401355*l2one*xipone*epsc2 + 0.4801655544988918*l2pone*xipone*epsc2 + -0.5970822295401355*l5one*xipone*epsc2 + -0.4801655544988918*l5pone*xipone*epsc2 + 0.2400827772494459*l2one*xippone*epsc2 + -0.2400827772494459*l5one*xippone*epsc2;
        }

        double P1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.03148856751270908 + -0.001974121999243681*as + -0.014123119863775352*epsb + 0.014123119863775352*epsc + 0.028246239727550703*epsb*etaone + 0.028246239727550703*epsc*etaone + 0.03148856751270908*l2one*epsc2 + 0.014123119863775352*l5one*epsc2 + -0.028246239727550703*l6one*epsc2;
        }

        double P1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.26142777033909836 + -0.04273192938942658*as + -0.11725448401958662*epsb + -1.0076341604066907*chi2one*epsb + 3.022902481220072*chi3pone*epsb + 0.11725448401958662*epsc + 1.0076341604066907*chi2one*epsc + -1.0076341604066907*chi3pone*epsc + 0.23450896803917323*epsb*etaone + 0.23450896803917323*epsc*etaone + 0.22596991782040562*epsb*etapone + 0.22596991782040562*epsc*etapone + 0.2519085401016727*xipone + -0.015792975993949448*as*xipone + -0.11298495891020281*epsb*xipone + 0.11298495891020281*epsc*xipone + 0.22596991782040562*epsb*etaone*xipone + 0.22596991782040562*epsc*etaone*xipone + 0.26142777033909836*l2one*epsc2 + 0.2519085401016727*l2pone*epsc2 + 0.2519085401016727*l3one*epsc2 + 0.11725448401958662*l5one*epsc2 + 0.11298495891020281*l5pone*epsc2 + -0.34749392694937603*l6one*epsc2 + -0.22596991782040562*l6pone*epsc2 + 0.2519085401016727*l2one*xipone*epsc2 + 0.11298495891020281*l5one*xipone*epsc2 + -0.22596991782040562*l6one*xipone*epsc2;
        }

        double P1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.7662176967652987 + -0.5869911628437433*as + -0.3436607387361962*epsb + -10.380956971664528*chi2one*epsb + -8.061073283253524*chi2pone*epsb + 31.142870914993587*chi3pone*epsb + 12.091609924880288*chi3ppone*epsb + 0.3436607387361962*epsc + 10.38095697166453*chi2one*epsc + 8.061073283253524*chi2pone*epsc + -10.380956971664528*chi3pone*epsc + -4.030536641626762*chi3ppone*epsc + 0.6873214774723924*epsb*etaone + 0.6873214774723924*epsc*etaone + 2.328011579954197*epsb*etapone + 2.328011579954197*epsc*etapone + 0.9038796712816225*epsb*etappone + 0.9038796712816225*epsc*etappone + 2.595239242916132*xipone + -0.3734413871033115*as*xipone + -1.1640057899770986*epsb*xipone + -8.061073283253524*chi2one*epsb*xipone + 24.183219849760576*chi3pone*epsb*xipone + 1.1640057899770986*epsc*xipone + 8.061073283253524*chi2one*epsc*xipone + -8.061073283253524*chi3pone*epsc*xipone + 2.328011579954197*epsb*etaone*xipone + 2.328011579954197*epsc*etaone*xipone + 1.807759342563245*epsb*etapone*xipone + 1.807759342563245*epsc*etapone*xipone + 1.0076341604066905*xippone + -0.06317190397579778*as*xippone + -0.45193983564081125*epsb*xippone + 0.45193983564081125*epsc*xippone + 0.9038796712816225*epsb*etaone*xippone + 0.9038796712816225*epsc*etaone*xippone + 0.7662176967652987*l2one*epsc2 + 2.091422162712787*l2pone*epsc2 + 2.5952392429161324*l3one*epsc2 + 2.015268320813381*l3pone*epsc2 + 0.3436607387361962*l5one*epsc2 + 0.9380358721566929*l5pone*epsc2 + -1.851327267449491*l6one*epsc2 + -2.7799514155950082*l6pone*epsc2 + 2.595239242916132*l2one*xipone*epsc2 + 2.015268320813381*l2pone*xipone*epsc2 + 2.015268320813381*l3one*xipone*epsc2 + 1.1640057899770986*l5one*xipone*epsc2 + 0.9038796712816225*l5pone*xipone*epsc2 + -3.23189125123582*l6one*xipone*epsc2 + -1.807759342563245*l6pone*xipone*epsc2 + 1.0076341604066905*l2one*xippone*epsc2 + 0.45193983564081125*l5one*xippone*epsc2 + -0.9038796712816225*l6one*xippone*epsc2;
        }

        double T1_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.008554226159037559 + 0.00487407755667564*as + -0.003836705538833072*epsb + 0.003836705538833072*epsc + 0.007673411077666144*epsb*etaone + 0.008554226159037559*l2one*epsc2 + -0.003836705538833072*l5one*epsc2;
        }

        double T1_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06798847109383505 + 0.04524757898962031*as + -0.03049390310389663*epsb + -0.2737352370892019*chi2one*epsb + 0.8212057112676057*chi3pone*epsb + 0.03049390310389663*epsc + -0.2737352370892019*chi3pone*epsc + 0.06098780620779326*epsb*etaone + 0.06138728862132917*epsb*etapone + 0.06843380927230047*xipone + 0.03899262045340512*as*xipone + -0.030693644310664583*epsb*xipone + 0.030693644310664583*epsc*xipone + 0.06138728862132917*epsb*etaone*xipone + 0.06798847109383505*l2one*epsc2 + 0.06843380927230047*l2pone*epsc2 + -0.03049390310389663*l5one*epsc2 + -0.030693644310664583*l5pone*epsc2 + 0.06843380927230047*l2one*xipone*epsc2 + -0.030693644310664583*l5one*xipone*epsc2;
        }

        double T1_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.17602072847139927 + 0.07214718628475085*as + -0.0789480768125533*epsb + -2.7231015491811252*chi2one*epsb + -2.189881896713615*chi2pone*epsb + 8.169304647543376*chi3pone*epsb + 3.2848228450704227*chi3ppone*epsb + 0.0789480768125533*epsc + -2.7231015491811252*chi3pone*epsc + -1.0949409483568076*chi3ppone*epsc + 0.1578961536251066*epsb*etaone + 0.6106770269050045*epsb*etapone + 0.24554915448531667*epsb*etappone + 0.6807753872952813*xipone + 0.4399658728237728*as*xipone + -0.30533851345250224*epsb*xipone + -2.189881896713615*chi2one*epsb*xipone + 6.569645690140845*chi3pone*epsb*xipone + 0.30533851345250224*epsc*xipone + -2.189881896713615*chi3pone*epsc*xipone + 0.6106770269050045*epsb*etaone*xipone + 0.49109830897063333*epsb*etapone*xipone + 0.2737352370892019*xippone + 0.15597048181362047*as*xippone + -0.12277457724265833*epsb*xippone + 0.12277457724265833*epsc*xippone + 0.24554915448531667*epsb*etaone*xippone + 0.17602072847139927*l2one*epsc2 + 0.5439077687506804*l2pone*epsc2 + -0.0789480768125533*l5one*epsc2 + -0.24395122483117307*l5pone*epsc2 + 0.6807753872952813*l2one*xipone*epsc2 + 0.5474704741784038*l2pone*xipone*epsc2 + -0.30533851345250224*l5one*xipone*epsc2 + -0.24554915448531667*l5pone*xipone*epsc2 + 0.2737352370892019*l2one*xippone*epsc2 + -0.12277457724265833*l5one*xippone*epsc2;
        }

        double T2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.0023576176034998757 + 0.0007880902517400929*as + 0.0023576176034998757*l2one*epsc2;
        }

        double T2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.020674096971735504 + 0.01926471690926969*as + -0.021025949448286316*epsb + -0.07544376331199602*chi2one*epsb + 0.22633128993598814*chi3pone*epsb + 0.021025949448286316*epsc + -0.07544376331199602*chi3pone*epsc + 0.04205189889657263*epsb*etaone + 0.018860940827999006*xipone + 0.006304722013920744*as*xipone + 0.020674096971735504*l2one*epsc2 + 0.018860940827999006*l2pone*epsc2 + -0.021025949448286316*l5one*epsc2 + 0.018860940827999006*l2one*xipone*epsc2;
        }

        double T2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06768609114155201 + 0.09651578481435866*as + -0.14232597356729948*epsb + -0.8124586297195282*chi2one*epsb + -0.6035501064959682*chi2pone*epsb + 2.437375889158585*chi3pone*epsb + 0.9053251597439524*chi3ppone*epsb + 0.14232597356729948*epsc + -0.8124586297195282*chi3pone*epsc + -0.3017750532479841*chi3ppone*epsc + 0.28465194713459896*epsb*etaone + 0.33641519117258106*epsb*etapone + 0.20311465742988205*xipone + 0.16672717930199904*as*xipone + -0.16820759558629053*epsb*xipone + -0.6035501064959682*chi2one*epsb*xipone + 1.810650319487905*chi3pone*epsb*xipone + 0.16820759558629053*epsc*xipone + -0.6035501064959682*chi3pone*epsc*xipone + 0.33641519117258106*epsb*etaone*xipone + 0.07544376331199602*xippone + 0.025218888055682977*as*xippone + 0.06768609114155201*l2one*epsc2 + 0.16539277577388406*l2pone*epsc2 + -0.14232597356729948*l5one*epsc2 + -0.16820759558629053*l5pone*epsc2 + 0.20311465742988205*l2one*xipone*epsc2 + 0.15088752662399205*l2pone*xipone*epsc2 + -0.16820759558629053*l5one*xipone*epsc2 + 0.07544376331199602*l2one*xippone*epsc2;
        }

        double T23_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.007038967893068947 + 0.002352943908547394*as + 0.007038967893068947*l2one*epsc2;
        }

        double T23_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06836929201240667 + 0.06431641970627533*as + 0.028155871572275844*epsb + -0.2252469725782063*chi2one*epsb + 0.6757409177346189*chi3pone*epsb + 0.02815587157227576*epsc + 0.22524697257820642*chi2one*epsc + -0.2252469725782063*chi3pone*epsc + -0.05631174314455169*epsb*etaone + 0.0563117431445516*epsc*etaone + 0.056311743144551576*xipone + 0.018823551268379153*as*xipone + 0.06836929201240667*l2one*epsc2 + 0.056311743144551576*l2pone*epsc2 + 0.056311743144551604*l3one*epsc2 + 0.028155871572275844*l5one*epsc2 + -0.056311743144551604*l6one*epsc2 + 0.056311743144551576*l2one*xipone*epsc2;
        }

        double T23_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.2669924753989712 + 0.40151283179319763*as + 0.21716542490507532*epsb + -2.638311289553426*chi2one*epsb + -1.8019757806256507*chi2pone*epsb + 7.914933868660278*chi3pone*epsb + 2.7029636709384763*chi3ppone*epsb + 0.21716542490507487*epsc + 2.6383112895534264*chi2one*epsc + 1.8019757806256507*chi2pone*epsc + -2.638311289553426*chi3pone*epsc + -0.9009878903128253*chi3ppone*epsc + -0.43433084981015063*epsb*etaone + 0.4343308498101503*epsc*etaone + -0.4504939451564132*epsb*etapone + 0.45049394515641267*epsc*etapone + 0.6595778223883565*xipone + 0.552178460186961*as*xipone + 0.2252469725782066*epsb*xipone + -1.8019757806256507*chi2one*epsb*xipone + 5.405927341876953*chi3pone*epsb*xipone + 0.22524697257820614*epsc*xipone + 1.8019757806256507*chi2one*epsc*xipone + -1.8019757806256507*chi3pone*epsc*xipone + -0.4504939451564132*epsb*etaone*xipone + 0.45049394515641267*epsc*etaone*xipone + 0.22524697257820633*xippone + 0.07529420507351661*as*xippone + 0.2669924753989712*l2one*epsc2 + 0.5469543360992533*l2pone*epsc2 + 0.6595778223883566*l3one*epsc2 + 0.45049394515641267*l3pone*epsc2 + 0.21716542490507532*l5one*epsc2 + 0.2252469725782066*l5pone*epsc2 + -0.6595778223883566*l6one*epsc2 + -0.45049394515641267*l6pone*epsc2 + 0.6595778223883565*l2one*xipone*epsc2 + 0.45049394515641267*l2pone*xipone*epsc2 + 0.45049394515641267*l3one*xipone*epsc2 + 0.2252469725782066*l5one*xipone*epsc2 + -0.45049394515641267*l6one*xipone*epsc2 + 0.22524697257820633*l2one*xippone*epsc2;
        }
        // }}}

        // B^* -> D form factors
        // {{{
        double P2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.043848623052024895 + -0.001635241305684815*as + -0.02108086828850025*epsb + 0.02108086828850025*epsc + -0.0421617365770005*epsb*etaone + -0.0421617365770005*epsc*etaone + 0.043848623052024895*l1one*epsc2 + -0.02108086828850025*l4one*epsc2;
        }

        double P2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.32960013132293325 + -0.04736745631879559*as + -0.1584600945860323*epsb + 1.4031559376647966*chi2one*epsb + -1.4031559376647966*chi3pone*epsb + 0.1584600945860323*epsc + -1.4031559376647966*chi2one*epsc + 4.209467812994389*chi3pone*epsc + -0.3169201891720646*epsb*etaone + -0.3169201891720646*epsc*etaone + -0.337293892616004*epsb*etapone + -0.337293892616004*epsc*etapone + 0.35078898441619916*xipone + -0.01308193044547852*as*xipone + -0.168646946308002*epsb*xipone + 0.168646946308002*epsc*xipone + -0.337293892616004*epsb*etaone*xipone + -0.337293892616004*epsc*etaone*xipone + 0.32960013132293325*l1one*epsc2 + 0.35078898441619916*l1pone*epsc2 + -0.1584600945860323*l4one*epsc2 + -0.168646946308002*l4pone*epsc2 + 0.35078898441619916*l1one*xipone*epsc2 + -0.168646946308002*l4one*xipone*epsc2;
        }

        double P2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.8529477601473066 + -0.7406079320469009*as + -0.4100671386488687*epsb + 13.353516077663457*chi2one*epsb + 11.225247501318373*chi2pone*epsb + -13.353516077663457*chi3pone*epsb + -5.6126237506591865*chi3ppone*epsb + 0.4100671386488687*epsc + -13.353516077663457*chi2one*epsc + -11.225247501318373*chi2pone*epsc + 40.06054823299037*chi3pone*epsc + 16.837871251977557*chi3ppone*epsc + -0.8201342772977374*epsb*etaone + -0.8201342772977374*epsc*etaone + -3.2099492986085245*epsb*etapone + -3.2099492986085245*epsc*etapone + -1.3491755704640158*epsb*etappone + -1.3491755704640158*epsc*etappone + 3.3383790194158642*xipone + -0.4051035114413218*as*xipone + -1.6049746493042623*epsb*xipone + 11.225247501318373*chi2one*epsb*xipone + -11.225247501318373*chi3pone*epsb*xipone + 1.6049746493042623*epsc*xipone + -11.225247501318373*chi2one*epsc*xipone + 33.675742503955114*chi3pone*epsc*xipone + -3.2099492986085245*epsb*etaone*xipone + -3.2099492986085245*epsc*etaone*xipone + -2.6983511409280316*epsb*etapone*xipone + -2.6983511409280316*epsc*etapone*xipone + 1.4031559376647966*xippone + -0.05232772178191408*as*xippone + -0.6745877852320079*epsb*xippone + 0.6745877852320079*epsc*xippone + -1.3491755704640158*epsb*etaone*xippone + -1.3491755704640158*epsc*etaone*xippone + 0.8529477601473066*l1one*epsc2 + 2.636801050583466*l1pone*epsc2 + -0.4100671386488687*l4one*epsc2 + -1.2676807566882584*l4pone*epsc2 + 3.3383790194158642*l1one*xipone*epsc2 + 2.8063118753295933*l1pone*xipone*epsc2 + -1.6049746493042623*l4one*xipone*epsc2 + -1.3491755704640158*l4pone*xipone*epsc2 + 1.4031559376647966*l1one*xippone*epsc2 + -0.6745877852320079*l4one*xippone*epsc2;
        }

        double V5_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.009099008915484865 + 0.008347177479748157*as + 0.009099008915484865*epsb + 0.009099008915484865*epsc + -0.01819801783096973*epsc*etaone + 0.009099008915484865*l1one*epsc2 + -0.009099008915484865*l4one*epsc2;
        }

        double V5_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.0669555721837609 + 0.07572072519390013*as + 0.0669555721837609*epsb + -0.2911682852955157*chi3pone*epsb + 0.0669555721837609*epsc + -0.2911682852955157*chi2one*epsc + 0.8735048558865469*chi3pone*epsc + -0.1339111443675218*epsc*etaone + -0.14558414264775785*epsc*etapone + 0.07279207132387892*xipone + 0.06677741983798526*as*xipone + 0.07279207132387892*epsb*xipone + 0.07279207132387892*epsc*xipone + -0.14558414264775785*epsc*etaone*xipone + 0.0669555721837609*l1one*epsc2 + 0.07279207132387892*l1pone*epsc2 + -0.0669555721837609*l4one*epsc2 + -0.07279207132387892*l4pone*epsc2 + 0.07279207132387892*l1one*xipone*epsc2 + -0.07279207132387892*l4one*xipone*epsc2;
        }

        double V5_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.15746000636048543 + 0.16807718349112652*as + 0.15746000636048543*epsb + -2.7249148804713803*chi3pone*epsb + -1.1646731411820628*chi3ppone*epsb + 0.15746000636048543*epsc + -2.7249148804713803*chi2one*epsc + -2.3293462823641256*chi2pone*epsc + 8.17474464141414*chi3pone*epsc + 3.4940194235461877*chi3ppone*epsc + -0.31492001272097087*epsc*etaone + -1.3624574402356902*epsc*etapone + -0.5823365705910314*epsc*etappone + 0.6812287201178451*xipone + 0.7393206412271716*as*xipone + 0.6812287201178451*epsb*xipone + -2.3293462823641256*chi3pone*epsb*xipone + 0.6812287201178451*epsc*xipone + -2.3293462823641256*chi2one*epsc*xipone + 6.988038847092375*chi3pone*epsc*xipone + -1.3624574402356902*epsc*etaone*xipone + -1.1646731411820628*epsc*etapone*xipone + 0.2911682852955157*xippone + 0.267109679351941*as*xippone + 0.2911682852955157*epsb*xippone + 0.2911682852955157*epsc*xippone + -0.5823365705910314*epsc*etaone*xippone + 0.15746000636048543*l1one*epsc2 + 0.5356445774700872*l1pone*epsc2 + -0.15746000636048543*l4one*epsc2 + -0.5356445774700872*l4pone*epsc2 + 0.6812287201178451*l1one*xipone*epsc2 + 0.5823365705910314*l1pone*xipone*epsc2 + -0.6812287201178451*l4one*xipone*epsc2 + -0.5823365705910314*l4pone*xipone*epsc2 + 0.2911682852955157*l1one*xippone*epsc2 + -0.2911682852955157*l4one*xippone*epsc2;
        }

        double A2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.013497592386244964 + -0.005614474340805168*as + 0.013497592386244964*l1one*epsc2;
        }

        double A2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.11544627665098015 + 0.021179445204968634*as + 0.053990369544979856*epsb + -0.43192295635983885*chi3pone*epsb + 0.053990369544979856*epsc + -0.43192295635983885*chi2one*epsc + 1.2957688690795166*chi3pone*epsc + -0.10798073908995971*epsc*etaone + 0.10798073908995971*xipone + -0.044915794726441347*as*xipone + 0.11544627665098015*l1one*epsc2 + 0.10798073908995971*l1pone*epsc2 + -0.053990369544979856*l4one*epsc2 + 0.10798073908995971*l1one*xipone*epsc2;
        }

        double A2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.39295357265360764 + 0.24981753639014034*as + 0.3538043675139609*epsb + -4.558126765551042*chi3pone*epsb + -1.7276918254393554*chi3ppone*epsb + 0.3538043675139609*epsc + -4.558126765551042*chi2one*epsc + -3.455383650878711*chi2pone*epsc + 13.674380296653126*chi3pone*epsc + 5.183075476318066*chi3ppone*epsc + -0.7076087350279218*epsc*etaone + -0.8638459127196777*epsc*etapone + 1.1395316913877604*xipone + 0.07960397218686636*as*xipone + 0.43192295635983885*epsb*xipone + -3.455383650878711*chi3pone*epsb*xipone + 0.43192295635983885*epsc*xipone + -3.455383650878711*chi2one*epsc*xipone + 10.366150952636133*chi3pone*epsc*xipone + -0.8638459127196777*epsc*etaone*xipone + 0.43192295635983885*xippone + -0.17966317890576539*as*xippone + 0.39295357265360764*l1one*epsc2 + 0.9235702132078412*l1pone*epsc2 + -0.3538043675139609*l4one*epsc2 + -0.43192295635983885*l4pone*epsc2 + 1.1395316913877604*l1one*xipone*epsc2 + 0.8638459127196777*l1pone*xipone*epsc2 + -0.43192295635983885*l4one*xipone*epsc2 + 0.43192295635983885*l1one*xippone*epsc2;
        }

        double A6_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.009779222979820786 + -0.004067777046606362*as + 0.009779222979820786*l1one*epsc2;
        }

        double A6_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06472590811759735 + 0.023998151689830795*as + -0.081363908984295*epsb + 0.3129351353542651*chi2one*epsb + -0.31293513535426515*chi3pone*epsb + 0.08136390898429502*epsc + -0.31293513535426515*chi2one*epsc + 0.9388054060627955*chi3pone*epsc + -0.16272781796859004*epsb*etaone + -0.16272781796859004*epsc*etaone + 0.07823378383856629*xipone + -0.032542216372850895*as*xipone + 0.06472590811759735*l1one*epsc2 + 0.07823378383856629*l1pone*epsc2 + -0.08136390898429502*l4one*epsc2 + 0.07823378383856629*l1one*xipone*epsc2;
        }

        double A6_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.14122128288513733 + 0.14055508337852513*as + -0.37579685908121613*epsb + 2.6970993304716457*chi2one*epsb + 2.5034810828341207*chi2pone*epsb + -2.6970993304716457*chi3pone*epsb + -1.2517405414170606*chi3ppone*epsb + 0.3757968590812163*epsc + -2.6970993304716457*chi2one*epsc + -2.503481082834121*chi2pone*epsc + 8.091297991414939*chi3pone*epsc + 3.755221624251182*chi3ppone*epsc + -0.7515937181624326*epsb*etaone + -0.7515937181624326*epsc*etaone + -1.3018225437487203*epsb*etapone + -1.3018225437487203*epsc*etapone + 0.6742748326179115*xipone + 0.1269007807729445*as*xipone + -0.65091127187436*epsb*xipone + 2.5034810828341207*chi2one*epsb*xipone + -2.503481082834121*chi3pone*epsb*xipone + 0.6509112718743602*epsc*xipone + -2.503481082834121*chi2one*epsc*xipone + 7.510443248502364*chi3pone*epsc*xipone + -1.3018225437487203*epsb*etaone*xipone + -1.3018225437487203*epsc*etaone*xipone + 0.31293513535426515*xippone + -0.13016886549140358*as*xippone + 0.14122128288513733*l1one*epsc2 + 0.5178072649407789*l1pone*epsc2 + -0.3757968590812163*l4one*epsc2 + -0.6509112718743602*l4pone*epsc2 + 0.6742748326179115*l1one*xipone*epsc2 + 0.6258702707085303*l1pone*xipone*epsc2 + -0.6509112718743602*l4one*xipone*epsc2 + 0.31293513535426515*l1one*xippone*epsc2;
        }

        double T1bar_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.010924441065250886 + 0.006409582853669373*as + -0.005252085178337256*epsb + 0.005252085178337256*epsc + -0.010504170356674511*epsc*etaone + 0.010924441065250886*l1one*epsc2 + -0.005252085178337256*l4one*epsc2;
        }

        double T1bar_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.08038811799237258 + 0.05563751173661732*as + -0.038647766096257516*epsb + -0.34958211408802836*chi3pone*epsb + 0.038647766096257516*epsc + -0.34958211408802836*chi2one*epsc + 1.048746342264085*chi3pone*epsc + -0.07729553219251503*epsc*etaone + -0.08403336285339609*epsc*etapone + 0.08739552852200709*xipone + 0.05127666282935498*as*xipone + -0.042016681426698045*epsb*xipone + 0.042016681426698045*epsc*xipone + -0.08403336285339609*epsc*etaone*xipone + 0.08038811799237258*l1one*epsc2 + 0.08739552852200709*l1pone*epsc2 + -0.038647766096257516*l4one*epsc2 + -0.042016681426698045*l4pone*epsc2 + 0.08739552852200709*l1one*xipone*epsc2 + -0.042016681426698045*l4one*xipone*epsc2;
        }

        double T1bar_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.18904944215317204 + 0.071218830706747*as + -0.09088829050155152*epsb + -3.27158400393198*chi3pone*epsb + -1.3983284563521134*chi3ppone*epsb + 0.09088829050155152*epsc + -3.27158400393198*chi2one*epsc + -2.796656912704227*chi2pone*epsc + 9.814752011795939*chi3pone*epsc + 4.19498536905634*chi3ppone*epsc + -0.18177658100310304*epsc*etaone + -0.7864309832469124*epsc*etapone + -0.3361334514135844*epsc*etappone + 0.817896000982995*xipone + 0.5476534195516486*as*xipone + -0.3932154916234562*epsb*xipone + -2.796656912704227*chi3pone*epsb*xipone + 0.3932154916234562*epsc*xipone + -2.796656912704227*chi2one*epsc*xipone + 8.38997073811268*chi3pone*epsc*xipone + -0.7864309832469124*epsc*etaone*xipone + -0.6722669028271688*epsc*etapone*xipone + 0.34958211408802836*xippone + 0.20510665131741992*as*xippone + -0.1680667257067922*epsb*xippone + 0.1680667257067922*epsc*xippone + -0.3361334514135844*epsc*etaone*xippone + 0.18904944215317204*l1one*epsc2 + 0.6431049439389808*l1pone*epsc2 + -0.09088829050155152*l4one*epsc2 + -0.3091821287700601*l4pone*epsc2 + 0.817896000982995*l1one*xipone*epsc2 + 0.6991642281760567*l1pone*xipone*epsc2 + -0.3932154916234562*l4one*xipone*epsc2 + -0.3361334514135844*l4pone*xipone*epsc2 + 0.34958211408802836*l1one*xippone*epsc2 + -0.1680667257067922*l4one*xippone*epsc2;
        }

        double T2bar_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return -0.004251589347148997 + -0.0014211957502845848*as + -0.004251589347148997*l1one*epsc2;
        }

        double T2bar_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return -0.03239165548585615 + -0.03177047341885325*as + 0.03537355978014203*epsb + 0.13605085910876794*chi3pone*epsb + -0.03537355978014203*epsc + 0.13605085910876794*chi2one*epsc + -0.4081525773263038*chi3pone*epsc + 0.07074711956028407*epsc*etaone + -0.034012714777191984*xipone + -0.011369566002276678*as*xipone + -0.03239165548585615*l1one*epsc2 + -0.034012714777191984*l1pone*epsc2 + 0.03537355978014203*l4one*epsc2 + -0.034012714777191984*l1one*xipone*epsc2;
        }

        double T2bar_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return -0.08953706088400629 + -0.12889154584954196*as + 0.19875401710877047*epsb + 1.308634693764933*chi3pone*epsb + 0.5442034364350717*chi3ppone*epsb + -0.19875401710877047*epsc + 1.308634693764933*chi2one*epsc + 1.0884068728701435*chi2pone*epsc + -3.9259040812947985*chi3pone*epsc + -1.6326103093052151*chi3ppone*epsc + 0.39750803421754094*epsc*etaone + 0.5659769564822726*epsc*etapone + -0.3271586734412332*xipone + -0.2769029193553795*as*xipone + 0.2829884782411363*epsb*xipone + 1.0884068728701435*chi3pone*epsb*xipone + -0.2829884782411363*epsc*xipone + 1.0884068728701435*chi2one*epsc*xipone + -3.2652206186104302*chi3pone*epsc*xipone + 0.5659769564822726*epsc*etaone*xipone + -0.13605085910876794*xippone + -0.045478264009106706*as*xippone + -0.08953706088400629*l1one*epsc2 + -0.2591332438868492*l1pone*epsc2 + 0.19875401710877047*l4one*epsc2 + 0.2829884782411363*l4pone*epsc2 + -0.3271586734412332*l1one*xipone*epsc2 + -0.27210171821753587*l1pone*xipone*epsc2 + 0.2829884782411363*l4one*xipone*epsc2 + -0.13605085910876794*l1one*xippone*epsc2;
        }

        double T23bar_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return -0.011736355765674576 + -0.003923158512284266*as + -0.011736355765674576*l1one*epsc2;
        }

        double T23bar_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return -0.10038224120467253 + -0.1026870142411883*as + -0.04694542306269827*epsb + -0.3755633845015864*chi2one*epsb + 0.37556338450158644*chi3pone*epsb + -0.046945423062698305*epsc + 0.37556338450158644*chi2one*epsc + -1.1266901535047595*chi3pone*epsc + -0.0938908461253966*epsb*etaone + 0.09389084612539661*epsc*etaone + -0.09389084612539661*xipone + -0.031385268098274126*as*xipone + -0.10038224120467253*l1one*epsc2 + -0.09389084612539661*l1pone*epsc2 + 0.046945423062698305*l4one*epsc2 + -0.09389084612539661*l1one*xipone*epsc2;
        }

        double T23bar_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return -0.34167893029244223 + -0.5546811615982133*as + -0.30763811869329344*epsb + -3.9633584875526924*chi2one*epsb + -3.004507076012691*chi2pone*epsb + 3.9633584875526937*chi3pone*epsb + 1.5022535380063458*chi3ppone*epsb + -0.30763811869329327*epsc + 3.9633584875526937*chi2one*epsc + 3.0045070760126915*chi2pone*epsc + -11.89007546265808*chi3pone*epsc + -4.506760614019037*chi3ppone*epsc + -0.6152762373865869*epsb*etaone + 0.6152762373865865*epsc*etaone + -0.7511267690031728*epsb*etapone + 0.7511267690031727*epsc*etapone + -0.9908396218881734*xipone + -0.8842666501260544*as*xipone + -0.37556338450158633*epsb*xipone + -3.004507076012691*chi2one*epsb*xipone + 3.0045070760126915*chi3pone*epsb*xipone + -0.37556338450158633*epsc*xipone + 3.0045070760126915*chi2one*epsc*xipone + -9.013521228038075*chi3pone*epsc*xipone + -0.7511267690031728*epsb*etaone*xipone + 0.7511267690031727*epsc*etaone*xipone + -0.37556338450158644*xippone + -0.12554107239309648*as*xippone + -0.34167893029244223*l1one*epsc2 + -0.8030579296373801*l1pone*epsc2 + 0.30763811869329327*l4one*epsc2 + 0.37556338450158633*l4pone*epsc2 + -0.9908396218881734*l1one*xipone*epsc2 + -0.7511267690031729*l1pone*xipone*epsc2 + 0.37556338450158633*l4one*xipone*epsc2 + -0.37556338450158644*l1one*xippone*epsc2;
        }
        // }}}

        // B^* -> D^* form factors
        // {{{
        double S2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.03876949172451258 + 0.009719724621896759*as + 0.03876949172451258*l2one*epsc2;
        }

        double S2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.2783004793038406 + 0.19603971918235066*as + -0.34313631589659094*epsb + -1.2406237351844025*chi3pone*epsb + 0.34313631589659094*epsc + -1.2406237351844025*chi3pone*epsc + 0.3101559337961006*xipone + 0.07775779697517408*as*xipone + 0.2783004793038406*l2one*epsc2 + 0.3101559337961006*l2pone*epsc2 + -0.34313631589659094*l5one*epsc2 + 0.3101559337961006*l2one*xipone*epsc2;
        }

        double S2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.7635413763985075 + 0.4936980746179209*as + -1.7768755017738411*epsb + -11.386862808091704*chi3pone*epsb + -4.96249494073761*chi3ppone*epsb + 1.7768755017738411*epsc + -11.386862808091704*chi3pone*epsc + -4.96249494073761*chi3ppone*epsc + 2.846715702022926*xipone + 1.7238333474091534*as*xipone + -2.7450905271727275*epsb*xipone + -9.92498988147522*chi3pone*epsb*xipone + 2.7450905271727275*epsc*xipone + -9.92498988147522*chi3pone*epsc*xipone + 1.2406237351844025*xippone + 0.31103118790069634*as*xippone + 0.7635413763985075*l2one*epsc2 + 2.226403834430725*l2pone*epsc2 + -1.7768755017738411*l5one*epsc2 + -2.7450905271727275*l5pone*epsc2 + 2.846715702022926*l2one*xipone*epsc2 + 2.481247470368805*l2pone*xipone*epsc2 + -2.7450905271727275*l5one*xipone*epsc2 + 1.2406237351844025*l2one*xippone*epsc2;
        }

        double S3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.027414170501558584 + 0.006872883191409052*as + 0.027414170501558584*l2one*epsc2;
        }

        double S3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.19678815612321215 + 0.1386210148157467*as + -0.2426340158418488*epsb + 0.8772534560498747*chi2one*epsb + -0.8772534560498748*chi3pone*epsb + 0.2426340158418488*epsc + 0.8772534560498747*chi2one*epsc + -0.8772534560498748*chi3pone*epsc + -0.4852680316836975*epsb*etaone + 0.4852680316836975*epsc*etaone + 0.2193133640124687*xipone + 0.05498306553127241*as*xipone + 0.19678815612321215*l2one*epsc2 + 0.2193133640124687*l2pone*epsc2 + 0.21931336401246868*l3one*epsc2 + 0.24263401584184868*l5one*epsc2 + -0.4852680316836975*l6one*epsc2 + 0.2193133640124687*l2one*xipone*epsc2;
        }

        double S3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.5399052849678948 + 0.34909725642107414*as + -1.2564407166285332*epsb + 8.051727908042537*chi2one*epsb + 7.018027648398997*chi2pone*epsb + -8.051727908042537*chi3pone*epsb + -3.509013824199499*chi3ppone*epsb + 1.2564407166285332*epsc + 8.051727908042537*chi2one*epsc + 7.018027648398997*chi2pone*epsc + -8.051727908042537*chi3pone*epsc + -3.509013824199499*chi3ppone*epsc + -2.5128814332570646*epsb*etaone + 2.5128814332570646*epsc*etaone + -3.88214425346958*epsb*etapone + 3.88214425346958*epsc*etapone + 2.012931977010634*xipone + 1.2189342495885183*as*xipone + -1.9410721267347903*epsb*xipone + 7.018027648398997*chi2one*epsb*xipone + -7.018027648398998*chi3pone*epsb*xipone + 1.9410721267347903*epsc*xipone + 7.018027648398997*chi2one*epsc*xipone + -7.018027648398998*chi3pone*epsc*xipone + -3.88214425346958*epsb*etaone*xipone + 3.88214425346958*epsc*etaone*xipone + 0.8772534560498747*xippone + 0.21993226212508965*as*xippone + 0.5399052849678948*l2one*epsc2 + 1.5743052489856972*l2pone*epsc2 + 2.012931977010634*l3one*epsc2 + 1.7545069120997492*l3pone*epsc2 + 1.2564407166285314*l5one*epsc2 + 1.9410721267347895*l5pone*epsc2 + -4.453953559991855*l6one*epsc2 + -3.88214425346958*l6pone*epsc2 + 2.012931977010634*l2one*xipone*epsc2 + 1.7545069120997494*l2pone*xipone*epsc2 + 1.7545069120997492*l3one*xipone*epsc2 + 1.9410721267347895*l5one*xipone*epsc2 + -3.88214425346958*l6one*xipone*epsc2 + 0.8772534560498747*l2one*xippone*epsc2;
        }

        double P3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.036840326316375654 + -0.0022102022896583733*as + -0.016649718028464256*epsb + 0.016649718028464256*epsc + 0.036840326316375654*l2one*epsc2 + -0.016649718028464256*l5one*epsc2;
        }

        double P3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.3236923482372954 + -0.050095453298419586*as + -0.1462904068720701*epsb + -1.1788904421240212*chi3pone*epsb + 0.1462904068720701*epsc + -1.1788904421240212*chi3pone*epsc + 0.2947226105310053*xipone + -0.017681618317266987*as*xipone + -0.13319774422771405*epsb*xipone + 0.13319774422771405*epsc*xipone + 0.3236923482372954*l2one*epsc2 + 0.2947226105310053*l2pone*epsc2 + -0.1462904068720701*l5one*epsc2 + -0.13319774422771405*l5pone*epsc2 + 0.2947226105310053*l2one*xipone*epsc2 + -0.13319774422771405*l5one*xipone*epsc2;
        }

        double P3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 1.0139125139266247 + -0.7049681479226149*as + -0.4582303999545686*epsb + -12.715936027841494*chi3pone*epsb + -4.715561768496084*chi3ppone*epsb + 0.4582303999545686*epsc + -12.715936027841494*chi3pone*epsc + -4.715561768496084*chi3ppone*epsc + 3.1789840069603734*xipone + -0.4361268630218906*as*xipone + -1.436718743431989*epsb*xipone + -9.431123536992168*chi3pone*epsb*xipone + 1.436718743431989*epsc*xipone + -9.431123536992168*chi3pone*epsc*xipone + 1.178890442124021*xippone + -0.07072647326906795*as*xippone + -0.5327909769108561*epsb*xippone + 0.5327909769108561*epsc*xippone + 1.0139125139266247*l2one*epsc2 + 2.589538785898363*l2pone*epsc2 + -0.4582303999545686*l5one*epsc2 + -1.1703232549765608*l5pone*epsc2 + 3.1789840069603734*l2one*xipone*epsc2 + 2.357780884248042*l2pone*xipone*epsc2 + -1.436718743431989*l5one*xipone*epsc2 + -1.0655819538217122*l5pone*xipone*epsc2 + 1.178890442124021*l2one*xippone*epsc2 + -0.5327909769108561*l5one*xippone*epsc2;
        }

        double V2_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.007335151340155248 + 0.0027093201917110336*as + -0.0033150683970844786*epsb + 0.0033150683970844786*epsc + 0.007335151340155248*l2one*epsc2 + -0.0033150683970844786*l5one*epsc2;
        }

        double V2_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.061094733790494654 + 0.02865790830563056*as + -0.02761132140633808*epsb + -0.23472484288496795*chi3pone*epsb + 0.02761132140633808*epsc + -0.23472484288496795*chi3pone*epsc + 0.05868121072124199*xipone + 0.02167456153368827*as*xipone + -0.02652054717667583*epsb*xipone + 0.02652054717667583*epsc*xipone + 0.061094733790494654*l2one*epsc2 + 0.05868121072124199*l2pone*epsc2 + -0.02761132140633808*l5one*epsc2 + -0.02652054717667583*l5pone*epsc2 + 0.05868121072124199*l2one*xipone*epsc2 + -0.02652054717667583*l5one*xipone*epsc2;
        }

        double V2_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.16740988852041994 + 0.0445707187184837*as + -0.07565968376894235*epsb + -2.424481167065765*chi3pone*epsb + -0.9388993715398718*chi3ppone*epsb + 0.07565968376894235*epsc + -2.424481167065765*chi3pone*epsc + -0.9388993715398718*chi3ppone*epsc + 0.6061202917664412*xipone + 0.272612389512421*as*xipone + -0.2739316656040563*epsb*xipone + -1.8777987430797436*chi3pone*epsb*xipone + 0.2739316656040563*epsc*xipone + -1.8777987430797436*chi3pone*epsc*xipone + 0.23472484288496795*xippone + 0.08669824613475308*as*xippone + -0.10608218870670331*epsb*xippone + 0.10608218870670331*epsc*xippone + 0.16740988852041994*l2one*epsc2 + 0.48875787032395723*l2pone*epsc2 + -0.07565968376894235*l5one*epsc2 + -0.2208905712507046*l5pone*epsc2 + 0.6061202917664412*l2one*xipone*epsc2 + 0.4694496857699359*l2pone*xipone*epsc2 + -0.2739316656040563*l5one*xipone*epsc2 + -0.21216437741340663*l5pone*xipone*epsc2 + 0.23472484288496795*l2one*xippone*epsc2 + -0.10608218870670331*l5one*xippone*epsc2;
        }

        double V3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.00518673525365337 + 0.0019157786799645093*as + -0.0023441073436756537*epsb + 0.0023441073436756537*epsc + -0.004688214687351307*epsb*etaone + 0.004688214687351307*epsc*etaone + 0.00518673525365337*l2one*epsc2 + 0.0023441073436756537*l5one*epsc2 + -0.004688214687351307*l6one*epsc2;
        }

        double V3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.04320050055804569 + 0.020264201297533653*as + -0.019524152603942933*epsb + 0.16597552811690783*chi2one*epsb + -0.16597552811690783*chi3pone*epsb + 0.019524152603942933*epsc + 0.16597552811690783*chi2one*epsc + -0.16597552811690783*chi3pone*epsc + -0.039048305207885874*epsb*etaone + 0.039048305207885874*epsc*etaone + -0.037505717498810466*epsb*etapone + 0.037505717498810466*epsc*etapone + 0.04149388202922696*xipone + 0.015326229439716075*as*xipone + -0.01875285874940523*epsb*xipone + 0.01875285874940523*epsc*xipone + -0.037505717498810466*epsb*etaone*xipone + 0.037505717498810466*epsc*etaone*xipone + 0.04320050055804569*l2one*epsc2 + 0.04149388202922696*l2pone*epsc2 + 0.04149388202922696*l3one*epsc2 + 0.01952415260394294*l5one*epsc2 + 0.018752858749405237*l5pone*epsc2 + -0.05780116395729111*l6one*epsc2 + -0.037505717498810466*l6pone*epsc2 + 0.04149388202922696*l2one*xipone*epsc2 + 0.018752858749405237*l5one*xipone*epsc2 + -0.037505717498810466*l6one*xipone*epsc2;
        }

        double V3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.11837666741047292 + 0.03151625744819802*as + -0.05349947545544891*epsb + 1.7143670740912775*chi2one*epsb + 1.3278042249352626*chi2pone*epsb + -1.7143670740912775*chi3pone*epsb + -0.6639021124676313*chi3ppone*epsb + 0.05349947545544891*epsc + 1.7143670740912775*chi2one*epsc + 1.3278042249352626*chi2pone*epsc + -1.7143670740912775*chi3pone*epsc + -0.6639021124676313*chi3ppone*epsc + -0.10699895091089782*epsb*etaone + 0.10699895091089782*epsc*etaone + -0.38739787666070796*epsb*etapone + 0.38739787666070796*epsc*etapone + -0.15002286999524184*epsb*etappone + 0.15002286999524184*epsc*etappone + 0.4285917685228194*xipone + 0.19276606925970147*as*xipone + -0.19369893833035395*epsb*xipone + 1.3278042249352626*chi2one*epsb*xipone + -1.3278042249352626*chi3pone*epsb*xipone + 0.19369893833035395*epsc*xipone + 1.3278042249352626*chi2one*epsc*xipone + -1.3278042249352626*chi3pone*epsc*xipone + -0.38739787666070796*epsb*etaone*xipone + 0.38739787666070796*epsc*etaone*xipone + -0.3000457399904837*epsb*etapone*xipone + 0.3000457399904837*epsc*etapone*xipone + 0.16597552811690783*xippone + 0.0613049177588643*as*xippone + -0.07501143499762092*epsb*xippone + 0.07501143499762092*epsc*xippone + -0.15002286999524184*epsb*etaone*xippone + 0.15002286999524184*epsc*etaone*xippone + 0.11837666741047292*l2one*epsc2 + 0.34560400446436546*l2pone*epsc2 + 0.4285917685228194*l3one*epsc2 + 0.33195105623381566*l3pone*epsc2 + 0.05349947545544891*l5one*epsc2 + 0.15619322083154355*l5pone*epsc2 + -0.3006978892412518*l6one*epsc2 + -0.46240931165832894*l6pone*epsc2 + 0.4285917685228194*l2one*xipone*epsc2 + 0.33195105623381566*l2pone*xipone*epsc2 + 0.33195105623381566*l3one*xipone*epsc2 + 0.19369893833035398*l5one*xipone*epsc2 + 0.15002286999524184*l5pone*xipone*epsc2 + -0.5374207466559499*l6one*xipone*epsc2 + -0.3000457399904837*l6pone*xipone*epsc2 + 0.16597552811690783*l2one*xippone*epsc2 + 0.07501143499762092*l5one*xippone*epsc2 + -0.15002286999524184*l6one*xippone*epsc2;
        }

        double V6_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.006543299374599415 + 0.006002640693093107*as + 0.006543299374599415*epsb + 0.006543299374599415*epsc + 0.01308659874919883*epsc*etaone + 0.006543299374599415*l2one*epsc2 + 0.006543299374599415*l5one*epsc2 + -0.01308659874919883*l6one*epsc2;
        }

        double V6_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.05449937088744518 + 0.060277889384701494*as + 0.05449937088744518*epsb + -0.20938557998718127*chi3pone*epsb + 0.05449937088744518*epsc + 0.20938557998718127*chi2one*epsc + -0.20938557998718127*chi3pone*epsc + 0.10899874177489036*epsc*etaone + 0.10469278999359063*epsc*etapone + 0.05234639499679532*xipone + 0.048021125544744865*as*xipone + 0.05234639499679532*epsb*xipone + 0.05234639499679532*epsc*xipone + 0.10469278999359063*epsc*etaone*xipone + 0.05449937088744518*l2one*epsc2 + 0.05234639499679532*l2pone*epsc2 + 0.05234639499679532*l3one*epsc2 + 0.05449937088744518*l5one*epsc2 + 0.05234639499679532*l5pone*epsc2 + -0.1613451367716857*l6one*epsc2 + -0.10469278999359063*l6pone*epsc2 + 0.05234639499679532*l2one*xipone*epsc2 + 0.05234639499679532*l5one*xipone*epsc2 + -0.10469278999359063*l6one*xipone*epsc2;
        }

        double V6_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1493374803135605 + 0.16396744598754961*as + 0.1493374803135605*epsb + -2.1627510283726084*chi3pone*epsb + -0.8375423199487251*chi3ppone*epsb + 0.1493374803135605*epsc + 2.1627510283726084*chi2one*epsc + 1.6750846398974502*chi2pone*epsc + -2.1627510283726084*chi3pone*epsc + -0.8375423199487251*chi3ppone*epsc + 0.298674960627121*epsc*etaone + 1.0813755141863042*epsc*etapone + 0.41877115997436254*epsc*etappone + 0.5406877570931521*xipone + 0.5782653661671017*as*xipone + 0.5406877570931521*epsb*xipone + -1.6750846398974502*chi3pone*epsb*xipone + 0.5406877570931521*epsc*xipone + 1.6750846398974502*chi2one*epsc*xipone + -1.6750846398974502*chi3pone*epsc*xipone + 1.0813755141863042*epsc*etaone*xipone + 0.8375423199487251*epsc*etapone*xipone + 0.20938557998718127*xippone + 0.19208450217897946*as*xippone + 0.20938557998718127*epsb*xippone + 0.20938557998718127*epsc*xippone + 0.41877115997436254*epsc*etaone*xippone + 0.1493374803135605*l2one*epsc2 + 0.4359949670995614*l2pone*epsc2 + 0.5406877570931521*l3one*epsc2 + 0.41877115997436254*l3pone*epsc2 + 0.1493374803135605*l5one*epsc2 + 0.4359949670995614*l5pone*epsc2 + -0.839362717720273*l6one*epsc2 + -1.2907610941734855*l6pone*epsc2 + 0.5406877570931521*l2one*xipone*epsc2 + 0.41877115997436254*l2pone*xipone*epsc2 + 0.41877115997436254*l3one*xipone*epsc2 + 0.5406877570931521*l5one*xipone*epsc2 + 0.41877115997436254*l5pone*xipone*epsc2 + -1.5001466741606668*l6one*xipone*epsc2 + -0.8375423199487251*l6pone*xipone*epsc2 + 0.20938557998718127*l2one*xippone*epsc2 + 0.20938557998718127*l5one*xippone*epsc2 + -0.41877115997436254*l6one*xippone*epsc2;
        }

        double V7_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.006543299374599415 + 0.006002640693093107*as + 0.006543299374599415*epsb + 0.006543299374599415*epsc + 0.01308659874919883*epsb*etaone + 0.006543299374599415*l2one*epsc2 + -0.006543299374599415*l5one*epsc2;
        }

        double V7_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.05449937088744518 + 0.060277889384701494*as + 0.05449937088744518*epsb + 0.20938557998718127*chi2one*epsb + -0.20938557998718127*chi3pone*epsb + 0.05449937088744518*epsc + -0.20938557998718127*chi3pone*epsc + 0.10899874177489036*epsb*etaone + 0.10469278999359063*epsb*etapone + 0.05234639499679532*xipone + 0.048021125544744865*as*xipone + 0.05234639499679532*epsb*xipone + 0.05234639499679532*epsc*xipone + 0.10469278999359063*epsb*etaone*xipone + 0.05449937088744518*l2one*epsc2 + 0.05234639499679532*l2pone*epsc2 + -0.05449937088744518*l5one*epsc2 + -0.05234639499679532*l5pone*epsc2 + 0.05234639499679532*l2one*xipone*epsc2 + -0.05234639499679532*l5one*xipone*epsc2;
        }

        double V7_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1493374803135605 + 0.16396744598754961*as + 0.1493374803135605*epsb + 2.1627510283726084*chi2one*epsb + 1.6750846398974502*chi2pone*epsb + -2.1627510283726084*chi3pone*epsb + -0.8375423199487251*chi3ppone*epsb + 0.1493374803135605*epsc + -2.1627510283726084*chi3pone*epsc + -0.8375423199487251*chi3ppone*epsc + 0.298674960627121*epsb*etaone + 1.0813755141863042*epsb*etapone + 0.41877115997436254*epsb*etappone + 0.5406877570931521*xipone + 0.5782653661671017*as*xipone + 0.5406877570931521*epsb*xipone + 1.6750846398974502*chi2one*epsb*xipone + -1.6750846398974502*chi3pone*epsb*xipone + 0.5406877570931521*epsc*xipone + -1.6750846398974502*chi3pone*epsc*xipone + 1.0813755141863042*epsb*etaone*xipone + 0.8375423199487251*epsb*etapone*xipone + 0.20938557998718127*xippone + 0.19208450217897946*as*xippone + 0.20938557998718127*epsb*xippone + 0.20938557998718127*epsc*xippone + 0.41877115997436254*epsb*etaone*xippone + 0.1493374803135605*l2one*epsc2 + 0.4359949670995614*l2pone*epsc2 + -0.1493374803135605*l5one*epsc2 + -0.4359949670995614*l5pone*epsc2 + 0.5406877570931521*l2one*xipone*epsc2 + 0.41877115997436254*l2pone*xipone*epsc2 + -0.5406877570931521*l5one*xipone*epsc2 + -0.41877115997436254*l5pone*xipone*epsc2 + 0.20938557998718127*l2one*xippone*epsc2 + -0.20938557998718127*l5one*xippone*epsc2;
        }

        double A3_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.006711208710683108 + -0.002791602237175056*as + 0.006711208710683108*l2one*epsc2;
        }

        double A3_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06881919420215553 + 0.005781487474467088*as + 0.02684483484273243*epsb + 0.21475867874185944*chi2one*epsb + -0.21475867874185944*chi3pone*epsb + 0.02684483484273243*epsc + -0.21475867874185944*chi3pone*epsc + 0.05368966968546486*epsb*etaone + 0.05368966968546486*xipone + -0.02233281789740045*as*xipone + 0.06881919420215553*l2one*epsc2 + 0.05368966968546486*l2pone*epsc2 + -0.02684483484273243*l5one*epsc2 + 0.05368966968546486*l2one*xipone*epsc2;
        }

        double A3_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.28459648348892935 + 0.1456399881376694*as + 0.22158710712315724*epsb + 2.6317315719526957*chi2one*epsb + 1.7180694299348755*chi2pone*epsb + -2.6317315719526957*chi3pone*epsb + -0.8590347149674378*chi3ppone*epsb + 0.22158710712315724*epsc + -2.6317315719526957*chi3pone*epsc + -0.8590347149674378*chi3ppone*epsc + 0.4431742142463145*epsb*etaone + 0.4295173574837189*epsb*etapone + 0.6579328929881739*xipone + 0.0015862640009358159*as*xipone + 0.21475867874185944*epsb*xipone + 1.7180694299348755*chi2one*epsb*xipone + -1.7180694299348755*chi3pone*epsb*xipone + 0.21475867874185944*epsc*xipone + -1.7180694299348755*chi3pone*epsc*xipone + 0.4295173574837189*epsb*etaone*xipone + 0.21475867874185944*xippone + -0.0893312715896018*as*xippone + 0.28459648348892935*l2one*epsc2 + 0.5505535536172442*l2pone*epsc2 + -0.22158710712315724*l5one*epsc2 + -0.21475867874185944*l5pone*epsc2 + 0.6579328929881739*l2one*xipone*epsc2 + 0.4295173574837189*l2pone*xipone*epsc2 + -0.21475867874185944*l5one*xipone*epsc2 + 0.21475867874185944*l2one*xippone*epsc2;
        }

        double A4_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.006711208710683108 + -0.002791602237175056*as + 0.006711208710683108*l2one*epsc2;
        }

        double A4_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.06881919420215553 + 0.005781487474467088*as + 0.02684483484273243*epsb + -0.21475867874185944*chi3pone*epsb + 0.02684483484273243*epsc + 0.21475867874185944*chi2one*epsc + -0.21475867874185944*chi3pone*epsc + 0.05368966968546486*epsc*etaone + 0.05368966968546486*xipone + -0.02233281789740045*as*xipone + 0.06881919420215553*l2one*epsc2 + 0.05368966968546486*l2pone*epsc2 + 0.05368966968546486*l3one*epsc2 + 0.02684483484273243*l5one*epsc2 + -0.05368966968546486*l6one*epsc2 + 0.05368966968546486*l2one*xipone*epsc2;
        }

        double A4_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.28459648348892935 + 0.1456399881376694*as + 0.22158710712315724*epsb + -2.6317315719526957*chi3pone*epsb + -0.8590347149674378*chi3ppone*epsb + 0.22158710712315724*epsc + 2.6317315719526957*chi2one*epsc + 1.7180694299348755*chi2pone*epsc + -2.6317315719526957*chi3pone*epsc + -0.8590347149674378*chi3ppone*epsc + 0.4431742142463145*epsc*etaone + 0.4295173574837189*epsc*etapone + 0.6579328929881739*xipone + 0.0015862640009358714*as*xipone + 0.21475867874185944*epsb*xipone + -1.7180694299348755*chi3pone*epsb*xipone + 0.21475867874185944*epsc*xipone + 1.7180694299348755*chi2one*epsc*xipone + -1.7180694299348755*chi3pone*epsc*xipone + 0.4295173574837189*epsc*etaone*xipone + 0.21475867874185944*xippone + -0.0893312715896018*as*xippone + 0.28459648348892935*l2one*epsc2 + 0.5505535536172442*l2pone*epsc2 + 0.6579328929881739*l3one*epsc2 + 0.4295173574837189*l3pone*epsc2 + 0.22158710712315724*l5one*epsc2 + 0.21475867874185944*l5pone*epsc2 + -0.6579328929881739*l6one*epsc2 + -0.4295173574837189*l6pone*epsc2 + 0.6579328929881739*l2one*xipone*epsc2 + 0.4295173574837189*l2pone*xipone*epsc2 + 0.4295173574837189*l3one*xipone*epsc2 + 0.21475867874185944*l5one*xipone*epsc2 + -0.4295173574837189*l6one*xipone*epsc2 + 0.21475867874185944*l2one*xippone*epsc2;
        }

        double A7_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.006412276517110341 + -0.0026672580517509013*as + 0.006412276517110341*l2one*epsc2;
        }

        double A7_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.05329514115021218 + 0.013900622461986156*as + -0.05675299940030692*epsb + -0.20519284854753092*chi3pone*epsb + 0.05675299940030692*epsc + -0.20519284854753092*chi3pone*epsc + 0.05129821213688273*xipone + -0.021338064414007214*as*xipone + 0.05329514115021218*l2one*epsc2 + 0.05129821213688273*l2pone*epsc2 + -0.05675299940030692*l5one*epsc2 + 0.05129821213688273*l2one*xipone*epsc2;
        }

        double A7_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1562776264382128 + 0.14297208702450237*as + -0.35819217355149036*epsb + -2.115830213901852*chi3pone*epsb + -0.8207713941901237*chi3ppone*epsb + 0.35819217355149036*epsc + -2.115830213901852*chi3pone*epsc + -0.8207713941901237*chi3ppone*epsc + 0.528957553475463*xipone + 0.06852885086787484*as*xipone + -0.4540239952024553*epsb*xipone + -1.6415427883802474*chi3pone*epsb*xipone + 0.4540239952024553*epsc*xipone + -1.6415427883802474*chi3pone*epsc*xipone + 0.20519284854753092*xippone + -0.08535225765602886*as*xippone + 0.1562776264382128*l2one*epsc2 + 0.4263611292016975*l2pone*epsc2 + -0.35819217355149036*l5one*epsc2 + -0.4540239952024553*l5pone*epsc2 + 0.528957553475463*l2one*xipone*epsc2 + 0.41038569709506184*l2pone*xipone*epsc2 + -0.4540239952024553*l5one*xipone*epsc2 + 0.20519284854753092*l2one*xippone*epsc2;
        }

        double T4_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.007954146388511356 + 0.002658864279089947*as + 0.007954146388511356*l2one*epsc2;
        }

        double T4_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.08156473276594307 + 0.07411808069728959*as + 0.03181658555404542*epsb + -0.2545326844323634*chi3pone*epsb + 0.03181658555404542*epsc + -0.2545326844323634*chi3pone*epsc + 0.06363317110809084*xipone + 0.021270914232719575*as*xipone + 0.08156473276594307*l2one*epsc2 + 0.06363317110809084*l2pone*epsc2 + -0.03181658555404542*l5one*epsc2 + 0.06363317110809084*l2one*xipone*epsc2;
        }

        double T4_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.3373046777286239 + 0.4909817757126442*as + 0.26262575995568144*epsb + -3.1191368173749052*chi3pone*epsb + -1.0181307377294535*chi3ppone*epsb + 0.26262575995568144*epsc + -3.1191368173749052*chi3pone*epsc + -1.0181307377294535*chi3ppone*epsc + 0.7797842043437263*xipone + 0.6354864740437558*as*xipone + 0.2545326844323634*epsb*xipone + -2.036261475458907*chi3pone*epsb*xipone + 0.2545326844323634*epsc*xipone + -2.036261475458907*chi3pone*epsc*xipone + 0.2545326844323634*xippone + 0.0850836569308783*as*xippone + 0.3373046777286239*l2one*epsc2 + 0.6525178621275446*l2pone*epsc2 + -0.26262575995568144*l5one*epsc2 + -0.2545326844323634*l5pone*epsc2 + 0.7797842043437263*l2one*xipone*epsc2 + 0.5090653688647268*l2pone*xipone*epsc2 + -0.2545326844323634*l5one*xipone*epsc2 + 0.2545326844323634*l2one*xippone*epsc2;
        }

        double T5_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.0018999627451430382 + 0.0006351081345396198*as + 0.0018999627451430382*l2one*epsc2;
        }

        double T5_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.017691356396422223 + 0.01580208901755895*as + -0.016815959861989967*epsb + -0.06079880784457723*chi3pone*epsb + 0.016815959861989967*epsc + 0.06079880784457723*chi2one*epsc + -0.06079880784457723*chi3pone*epsc + 0.033631919723979935*epsc*etaone + 0.015199701961144308*xipone + 0.0050808650763169575*as*xipone + 0.017691356396422223*l2one*epsc2 + 0.015199701961144308*l2pone*epsc2 + 0.015199701961144308*l3one*epsc2 + 0.016815959861989967*l5one*epsc2 + -0.033631919723979935*l6one*epsc2 + 0.015199701961144308*l2one*xipone*epsc2;
        }

        double T5_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.062096581417229 + 0.08515196368715783*as + -0.1229485920923181*epsb + -0.6877210203746656*chi3pone*epsb + -0.24319523137830895*chi3ppone*epsb + 0.1229485920923181*epsc + 0.6877210203746656*chi2one*epsc + 0.4863904627566179*chi2pone*epsc + -0.6877210203746656*chi3pone*epsc + -0.24319523137830895*chi3ppone*epsc + 0.2458971841846362*epsc*etaone + 0.2690553577918395*epsc*etapone + 0.1719302550936664*xipone + 0.1365784422931055*as*xipone + -0.13452767889591974*epsb*xipone + -0.4863904627566179*chi3pone*epsb*xipone + 0.13452767889591974*epsc*xipone + 0.4863904627566179*chi2one*epsc*xipone + -0.4863904627566179*chi3pone*epsc*xipone + 0.2690553577918395*epsc*etaone*xipone + 0.06079880784457724*xippone + 0.020323460305267833*as*xippone + 0.062096581417229*l2one*epsc2 + 0.1415308511713778*l2pone*epsc2 + 0.1719302550936664*l3one*epsc2 + 0.12159761568915448*l3pone*epsc2 + 0.1229485920923181*l5one*epsc2 + 0.13452767889591974*l5pone*epsc2 + -0.38042486308055595*l6one*epsc2 + -0.2690553577918395*l6pone*epsc2 + 0.1719302550936664*l2one*xipone*epsc2 + 0.12159761568915448*l2pone*xipone*epsc2 + 0.12159761568915448*l3one*xipone*epsc2 + 0.13452767889591974*l5one*xipone*epsc2 + -0.2690553577918395*l6one*xipone*epsc2 + 0.06079880784457724*l2one*xippone*epsc2;
        }

        double T6_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.0018999627451430382 + 0.0006351081345396198*as + 0.0018999627451430382*l2one*epsc2;
        }

        double T6_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.017691356396422223 + 0.01580208901755895*as + -0.016815959861989964*epsb + 0.060798807844577224*chi2one*epsb + -0.060798807844577224*chi3pone*epsb + 0.01681595986198997*epsc + -0.060798807844577224*chi3pone*epsc + -0.033631919723979935*epsb*etaone + 0.015199701961144306*xipone + 0.005080865076316958*as*xipone + 0.017691356396422223*l2one*epsc2 + 0.015199701961144306*l2pone*epsc2 + -0.01681595986198997*l5one*epsc2 + 0.015199701961144306*l2one*xipone*epsc2;
        }

        double T6_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.062096581417229 + 0.08515196368715786*as + -0.1229485920923181*epsb + 0.6877210203746656*chi2one*epsb + 0.4863904627566179*chi2pone*epsb + -0.6877210203746656*chi3pone*epsb + -0.24319523137830895*chi3ppone*epsb + 0.1229485920923181*epsc + -0.6877210203746656*chi3pone*epsc + -0.24319523137830895*chi3ppone*epsc + -0.2458971841846362*epsb*etaone + -0.2690553577918395*epsb*etapone + 0.1719302550936664*xipone + 0.13657844229310553*as*xipone + -0.1345276788959197*epsb*xipone + 0.4863904627566179*chi2one*epsb*xipone + -0.4863904627566179*chi3pone*epsb*xipone + 0.13452767889591974*epsc*xipone + -0.4863904627566179*chi3pone*epsc*xipone + -0.2690553577918395*epsb*etaone*xipone + 0.06079880784457724*xippone + 0.020323460305267833*as*xippone + 0.062096581417229*l2one*epsc2 + 0.1415308511713778*l2pone*epsc2 + -0.1229485920923181*l5one*epsc2 + -0.13452767889591974*l5pone*epsc2 + 0.1719302550936664*l2one*xipone*epsc2 + 0.12159761568915448*l2pone*xipone*epsc2 + -0.13452767889591974*l5one*xipone*epsc2 + 0.06079880784457724*l2one*xippone*epsc2;
        }

        double T7_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.014082253085597513 + 0.023483668425490135*as + 0.014082253085597513*epsb + 0.014082253085597513*epsc + 0.014082253085597513*l2one*epsc2 + -0.014082253085597513*l5one*epsc2;
        }

        double T7_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.13057035472750284 + 0.2506199645069021*as + 0.13057035472750284*epsb + -0.4506320987391204*chi3pone*epsb + 0.13057035472750284*epsc + -0.4506320987391204*chi3pone*epsc + 0.1126580246847801*xipone + 0.1878693474039211*as*xipone + 0.1126580246847801*epsb*xipone + 0.1126580246847801*epsc*xipone + 0.13057035472750284*l2one*epsc2 + 0.1126580246847801*l2pone*epsc2 + -0.13057035472750284*l5one*epsc2 + -0.1126580246847801*l5pone*epsc2 + 0.1126580246847801*l2one*xipone*epsc2 + -0.1126580246847801*l5one*xipone*epsc2;
        }

        double T7_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.44527688302044277 + 0.9265206779774454*as + 0.44527688302044277*epsb + -5.079515548758332*chi3pone*epsb + -1.8025283949564816*chi3ppone*epsb + 0.44527688302044277*epsc + -5.079515548758332*chi3pone*epsc + -1.8025283949564816*chi3ppone*epsc + 1.269878887189583*xipone + 2.380698410863059*as*xipone + 1.269878887189583*epsb*xipone + -3.6050567899129633*chi3pone*epsb*xipone + 1.269878887189583*epsc*xipone + -3.6050567899129633*chi3pone*epsc*xipone + 0.4506320987391204*xippone + 0.7514773896156844*as*xippone + 0.4506320987391204*epsb*xippone + 0.4506320987391204*epsc*xippone + 0.44527688302044277*l2one*epsc2 + 1.0445628378200227*l2pone*epsc2 + -0.44527688302044277*l5one*epsc2 + -1.0445628378200227*l5pone*epsc2 + 1.269878887189583*l2one*xipone*epsc2 + 0.9012641974782408*l2pone*xipone*epsc2 + -1.269878887189583*l5one*xipone*epsc2 + -0.9012641974782408*l5pone*xipone*epsc2 + 0.4506320987391204*l2one*xippone*epsc2 + -0.4506320987391204*l5one*xippone*epsc2;
        }

        double T8_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.007442866622238243 + 0.004254233697022471*as + -0.0033637495368397784*epsb + 0.0033637495368397784*epsc + -0.006727499073679558*epsb*etaone + 0.007442866622238243*l2one*epsc2 + -0.0033637495368397784*l5one*epsc2;
        }

        double T8_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.061991898167731534 + 0.041108351567110284*as + -0.028016788333473753*epsb + 0.23817173191162377*chi2one*epsb + -0.23817173191162377*chi3pone*epsb + 0.028016788333473753*epsc + -0.23817173191162377*chi3pone*epsc + -0.05603357666694753*epsb*etaone + -0.05381999258943646*epsb*etapone + 0.059542932977905944*xipone + 0.034033869576179765*as*xipone + -0.026909996294718228*epsb*xipone + 0.026909996294718228*epsc*xipone + -0.05381999258943646*epsb*etaone*xipone + 0.061991898167731534*l2one*epsc2 + 0.059542932977905944*l2pone*epsc2 + -0.028016788333473753*l5one*epsc2 + -0.026909996294718228*l5pone*epsc2 + 0.059542932977905944*l2one*xipone*epsc2 + -0.026909996294718228*l5one*xipone*epsc2;
        }

        double T8_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1698682704309257 + 0.07485642073680311*as + -0.07677073162624659*epsb + 2.4600842051906566*chi2one*epsb + 1.9053738552929902*chi2pone*epsb + -2.4600842051906566*chi3pone*epsb + -0.9526869276464951*chi3ppone*epsb + 0.07677073162624659*epsc + -2.4600842051906566*chi3pone*epsc + -0.9526869276464951*chi3ppone*epsc + -0.1535414632524932*epsb*etaone + -0.5559085985144531*epsb*etapone + -0.21527997035774585*epsb*etappone + 0.6150210512976642*xipone + 0.39693455168924174*as*xipone + -0.2779542992572265*epsb*xipone + 1.9053738552929902*chi2one*epsb*xipone + -1.9053738552929902*chi3pone*epsb*xipone + 0.2779542992572265*epsc*xipone + -1.9053738552929902*chi3pone*epsc*xipone + -0.5559085985144531*epsb*etaone*xipone + -0.4305599407154917*epsb*etapone*xipone + 0.23817173191162377*xippone + 0.13613547830471906*as*xippone + -0.10763998517887292*epsb*xippone + 0.10763998517887292*epsc*xippone + -0.21527997035774585*epsb*etaone*xippone + 0.1698682704309257*l2one*epsc2 + 0.49593518534185227*l2pone*epsc2 + -0.07677073162624659*l5one*epsc2 + -0.22413430666779002*l5pone*epsc2 + 0.6150210512976642*l2one*xipone*epsc2 + 0.47634346382324755*l2pone*xipone*epsc2 + -0.2779542992572265*l5one*xipone*epsc2 + -0.21527997035774585*l5pone*xipone*epsc2 + 0.23817173191162377*l2one*xippone*epsc2 + -0.10763998517887292*l5one*xippone*epsc2;
        }

        double T9_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.007442866622238245 + 0.00425423369702247*as + -0.0033637495368397793*epsb + 0.0033637495368397793*epsc + 0.006727499073679559*epsc*etaone + 0.007442866622238245*l2one*epsc2 + 0.0033637495368397793*l5one*epsc2 + -0.006727499073679559*l6one*epsc2;
        }

        double T9_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.061991898167731534 + 0.04110835156711028*as + -0.028016788333473763*epsb + -0.23817173191162383*chi3pone*epsb + 0.028016788333473763*epsc + 0.23817173191162383*chi2one*epsc + -0.23817173191162383*chi3pone*epsc + 0.05603357666694753*epsc*etaone + 0.05381999258943647*epsc*etapone + 0.05954293297790596*xipone + 0.03403386957617976*as*xipone + -0.026909996294718234*epsb*xipone + 0.026909996294718234*epsc*xipone + 0.05381999258943647*epsc*etaone*xipone + 0.061991898167731534*l2one*epsc2 + 0.05954293297790596*l2pone*epsc2 + 0.05954293297790596*l3one*epsc2 + 0.028016788333473763*l5one*epsc2 + 0.026909996294718234*l5pone*epsc2 + -0.08294357296166577*l6one*epsc2 + -0.05381999258943647*l6pone*epsc2 + 0.05954293297790596*l2one*xipone*epsc2 + 0.026909996294718234*l5one*xipone*epsc2 + -0.05381999258943647*l6one*xipone*epsc2;
        }

        double T9_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.1698682704309257 + 0.07485642073680306*as + -0.07677073162624659*epsb + -2.4600842051906566*chi3pone*epsb + -0.9526869276464951*chi3ppone*epsb + 0.07677073162624659*epsc + 2.4600842051906566*chi2one*epsc + 1.9053738552929902*chi2pone*epsc + -2.4600842051906566*chi3pone*epsc + -0.9526869276464951*chi3ppone*epsc + 0.15354146325249318*epsc*etaone + 0.5559085985144531*epsc*etapone + 0.21527997035774585*epsc*etappone + 0.6150210512976642*xipone + 0.39693455168924163*as*xipone + -0.27795429925722653*epsb*xipone + -1.9053738552929902*chi3pone*epsb*xipone + 0.27795429925722653*epsc*xipone + 1.9053738552929902*chi2one*epsc*xipone + -1.9053738552929902*chi3pone*epsc*xipone + 0.5559085985144531*epsc*etaone*xipone + 0.4305599407154917*epsc*etapone*xipone + 0.23817173191162377*xippone + 0.13613547830471903*as*xippone + -0.10763998517887292*epsb*xippone + 0.10763998517887292*epsc*xippone + 0.21527997035774585*epsc*etaone*xippone + 0.1698682704309257*l2one*epsc2 + 0.49593518534185227*l2pone*epsc2 + 0.6150210512976642*l3one*epsc2 + 0.47634346382324755*l3pone*epsc2 + 0.07677073162624659*l5one*epsc2 + 0.22413430666779008*l5pone*epsc2 + -0.43149576250971977*l6one*epsc2 + -0.663548583693326*l6pone*epsc2 + 0.6150210512976642*l2one*xipone*epsc2 + 0.47634346382324755*l2pone*xipone*epsc2 + 0.47634346382324755*l3one*xipone*epsc2 + 0.27795429925722653*l5one*xipone*epsc2 + 0.21527997035774585*l5pone*xipone*epsc2 + -0.7711885688721989*l6one*xipone*epsc2 + -0.4305599407154917*l6pone*xipone*epsc2 + 0.23817173191162377*l2one*xippone*epsc2 + 0.10763998517887292*l5one*xippone*epsc2 + -0.21527997035774585*l6one*xippone*epsc2;
        }

        double T10_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.009957656651211184 + 0.01660546119080049*as + 0.009957656651211184*epsb + 0.009957656651211184*epsc + 0.01991531330242237*epsb*etaone + 0.01991531330242237*epsc*etaone + 0.009957656651211184*l2one*epsc2 + 0.009957656651211184*l5one*epsc2 + -0.01991531330242237*l6one*epsc2;
        }

        double T10_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.09232718324975024 + 0.17721507640356235*as + 0.09232718324975024*epsb + 0.3186450128387579*chi2one*epsb + -0.3186450128387579*chi3pone*epsb + 0.09232718324975024*epsc + 0.3186450128387579*chi2one*epsc + -0.3186450128387579*chi3pone*epsc + 0.18465436649950048*epsb*etaone + 0.18465436649950048*epsc*etaone + 0.15932250641937895*epsb*etapone + 0.15932250641937895*epsc*etapone + 0.07966125320968948*xipone + 0.1328436895264039*as*xipone + 0.07966125320968948*epsb*xipone + 0.07966125320968948*epsc*xipone + 0.15932250641937895*epsb*etaone*xipone + 0.15932250641937895*epsc*etaone*xipone + 0.09232718324975024*l2one*epsc2 + 0.07966125320968948*l2pone*epsc2 + 0.07966125320968948*l3one*epsc2 + 0.09232718324975024*l5one*epsc2 + 0.07966125320968948*l5pone*epsc2 + -0.26431561970918993*l6one*epsc2 + -0.15932250641937895*l6pone*epsc2 + 0.07966125320968948*l2one*xipone*epsc2 + 0.07966125320968948*l5one*xipone*epsc2 + -0.15932250641937895*l6one*xipone*epsc2;
        }

        double T10_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b();
            const double epsc = _eps_c();
            const double epsc2 = power_of<2>(epsc);

            return 0.31485830348936417 + 0.6551490543074092*as + 0.31485830348936417*epsb + 3.591759889669524*chi2one*epsb + 2.549160102710063*chi2pone*epsb + -3.591759889669524*chi3pone*epsb + -1.2745800513550316*chi3ppone*epsb + 0.31485830348936417*epsc + 3.591759889669524*chi2one*epsc + 2.549160102710063*chi2pone*epsc + -3.591759889669524*chi3pone*epsc + -1.2745800513550316*chi3ppone*epsc + 0.6297166069787283*epsb*etaone + 0.6297166069787283*epsc*etaone + 1.795879944834762*epsb*etapone + 1.795879944834762*epsc*etapone + 0.6372900256775158*epsb*etappone + 0.6372900256775158*epsc*etappone + 0.897939972417381*xipone + 1.6834079902813064*as*xipone + 0.897939972417381*epsb*xipone + 2.549160102710063*chi2one*epsb*xipone + -2.549160102710063*chi3pone*epsb*xipone + 0.897939972417381*epsc*xipone + 2.549160102710063*chi2one*epsc*xipone + -2.549160102710063*chi3pone*epsc*xipone + 1.795879944834762*epsb*etaone*xipone + 1.795879944834762*epsc*etaone*xipone + 1.2745800513550316*epsb*etapone*xipone + 1.2745800513550316*epsc*etapone*xipone + 0.3186450128387579*xippone + 0.5313747581056155*as*xippone + 0.3186450128387579*epsb*xippone + 0.3186450128387579*epsc*xippone + 0.6372900256775158*epsb*etaone*xippone + 0.6372900256775158*epsc*etaone*xippone + 0.31485830348936417*l2one*epsc2 + 0.738617465998002*l2pone*epsc2 + 0.897939972417381*l3one*epsc2 + 0.6372900256775158*l3pone*epsc2 + 0.31485830348936417*l5one*epsc2 + 0.738617465998002*l5pone*epsc2 + -1.5276565793961092*l6one*epsc2 + -2.11452495767352*l6pone*epsc2 + 0.897939972417381*l2one*xipone*epsc2 + 0.6372900256775158*l2pone*xipone*epsc2 + 0.6372900256775158*l3one*xipone*epsc2 + 0.897939972417381*l5one*xipone*epsc2 + 0.6372900256775158*l5pone*xipone*epsc2 + -2.4331699705122776*l6one*xipone*epsc2 + -1.2745800513550316*l6pone*xipone*epsc2 + 0.3186450128387579*l2one*xippone*epsc2 + 0.3186450128387579*l5one*xippone*epsc2 + -0.6372900256775158*l6one*xippone*epsc2;
        }
        // }}}

        // B_s -> D_s form factors
        // {{{
        double V1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.004853927809740896 + 0.0018073507175946574*as + -0.0022489124823136626*epsb + 0.0022489124823136626*epsc + 0.004497824964627325*epsb*etasone + -0.004497824964627325*epsc*etasone + (0.004853927809740896*l1sone - 0.0022489124823136626*l4sone)*epsc2;
        }

        double V1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04085920220485458 + -0.018930806855320836*epsb + -0.15532568991170867*chi2sone*epsb + 0.465977069735126*chi3spone*epsb + 0.018930806855320836*epsc + -0.15532568991170867*chi2sone*epsc + 0.465977069735126*chi3spone*epsc + 0.03786161371064167*epsb*etasone + -0.03786161371064167*epsc*etasone + 0.0359825997170186*epsb*etaspone + -0.0359825997170186*epsc*etaspone + as*(0.01923252020198305 + 0.014458805740757258*xispone) + 0.03883142247792717*xispone + -0.0179912998585093*epsb*xispone + 0.0179912998585093*epsc*xispone + 0.0359825997170186*epsb*etasone*xispone + -0.0359825997170186*epsc*etasone*xispone + (0.03883142247792717*l1spone - 0.018930806855320836*l4sone - 0.0179912998585093*l4spone + l1sone*(0.04085920220485458 + 0.03883142247792717*xispone) - 0.0179912998585093*l4sone*xispone)*epsc2;
        }

        double V1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);


            return 0.11341891835977361 + -0.05254903476684066*epsb + -1.618145850378764*chi2sone*epsb + -1.2426055192936694*chi2spone*epsb + 4.854437551136292*chi3spone*epsb + 1.863908278940504*chi3sppone*epsb + 0.05254903476684066*epsc + -1.618145850378764*chi2sone*epsc + -1.2426055192936694*chi2spone*epsc + 4.854437551136292*chi3spone*epsc + 1.863908278940504*chi3sppone*epsc + 0.10509806953368132*epsb*etasone + -0.10509806953368132*epsc*etasone + 0.3748581091191706*epsb*etaspone + -0.3748581091191706*epsc*etaspone + 0.1439303988680744*epsb*etasppone + -0.1439303988680744*epsc*etasppone + 0.404536462594691*xispone + -0.1874290545595853*epsb*xispone + -1.2426055192936694*chi2sone*epsb*xispone + 3.727816557881008*chi3spone*epsb*xispone + 0.1874290545595853*epsc*xispone + -1.2426055192936694*chi2sone*epsc*xispone + 3.727816557881008*chi3spone*epsc*xispone + 0.3748581091191706*epsb*etasone*xispone + -0.3748581091191706*epsc*etasone*xispone + 0.2878607977361488*epsb*etaspone*xispone + -0.2878607977361488*epsc*etaspone*xispone + as*(0.03101685877692137 + 0.18277777309737894*xispone + 0.05783522296302903*xisppone) + 0.15532568991170867*xisppone + -0.0719651994340372*epsb*xisppone + 0.0719651994340372*epsc*xisppone + 0.1439303988680744*epsb*etasone*xisppone + -0.1439303988680744*epsc*etasone*xisppone + (0.32687361763883666*l1spone - 0.05254903476684066*l4sone - 0.1514464548425667*l4spone + 0.31065137982341734*l1spone*xispone - 0.1874290545595853*l4sone*xispone - 0.1439303988680744*l4spone*xispone + l1sone*(0.11341891835977361 + 0.404536462594691*xispone + 0.15532568991170867*xisppone) - 0.0719651994340372*l4sone*xisppone)*epsc2;
        }

        double S1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.027568784109983837 + 0.00691165359571278*as + 0.027568784109983837*l1sone*epsc2;
        }

        double S1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);


            return 0.1991754644782667 + 0.13815036774277073*as + -0.23801173042451831*epsb + -0.8822010915194828*chi2sone*epsb + 2.6466032745584482*chi3spone*epsb + 0.23801173042451831*epsc + -0.8822010915194828*chi2sone*epsc + 2.6466032745584482*chi3spone*epsc + 0.47602346084903663*epsb*etasone + -0.47602346084903663*epsc*etasone + 0.2205502728798707*xispone + 0.05529322876570224*as*xispone + 0.1991754644782667*l1sone*epsc2 + 0.2205502728798707*l1spone*epsc2 + -0.23801173042451831*l4sone*epsc2 + 0.2205502728798707*l1sone*xispone*epsc2;
        }

        double S1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);


            return 0.5502247632998395 + 0.35019441270807145*as + -1.2435335848826854*epsb + -8.138017046343501*chi2sone*epsb + -7.057608732155861*chi2spone*epsb + 24.414051139030498*chi3spone*epsb + 10.586413098233795*chi3sppone*epsb + 1.2435335848826854*epsc + -8.138017046343501*chi2sone*epsc + -7.057608732155861*chi2spone*epsc + 24.414051139030498*chi3spone*epsc + 10.586413098233795*chi3sppone*epsc + 2.487067169765371*epsb*etasone + -2.487067169765371*epsc*etasone + 3.808187686792293*epsb*etaspone + -3.808187686792293*epsc*etaspone + 2.0345042615858753*xispone + 1.2157893994735705*as*xispone + -1.9040938433961465*epsb*xispone + -7.057608732155861*chi2sone*epsb*xispone + 21.17282619646759*chi3spone*epsb*xispone + 1.9040938433961465*epsc*xispone + -7.057608732155861*chi2sone*epsc*xispone + 21.17282619646759*chi3spone*epsc*xispone + 3.808187686792293*epsb*etasone*xispone + -3.808187686792293*epsc*etasone*xispone + 0.8822010915194827*xisppone + 0.22117291506280895*as*xisppone + 0.5502247632998395*l1sone*epsc2 + 1.5934037158261338*l1spone*epsc2 + -1.2435335848826854*l4sone*epsc2 + -1.9040938433961465*l4spone*epsc2 + 2.0345042615858753*l1sone*xispone*epsc2 + 1.7644021830389653*l1spone*xispone*epsc2 + -1.9040938433961465*l4sone*xispone*epsc2 + 0.8822010915194827*l1sone*xisppone*epsc2;
        }

        double fTs_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.00966678746772539 + 0.016120405609235954*as + 0.00966678746772539*epsb + 0.00966678746772539*epsc + -0.01933357493545078*epsb*etasone + -0.01933357493545078*epsc*etasone + 0.00966678746772539*l1sone*epsc2 + -0.00966678746772539*l4sone*epsc2;
        }

        double fTs_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0904562210942245 + 0.1734159116203609*as + 0.0904562210942245*epsb + -0.30933719896721246*chi2sone*epsb + 0.9280115969016375*chi3spone*epsb + 0.0904562210942245*epsc + -0.30933719896721246*chi2sone*epsc + 0.9280115969016375*chi3spone*epsc + -0.180912442188449*epsb*etasone + -0.180912442188449*epsc*etasone + -0.15466859948360623*epsb*etaspone + -0.15466859948360623*epsc*etaspone + 0.07733429974180311*xispone + 0.12896324487388763*as*xispone + 0.07733429974180311*epsb*xispone + 0.07733429974180311*epsc*xispone + -0.15466859948360623*epsb*etasone*xispone + -0.15466859948360623*epsc*etasone*xispone + 0.07733429974180311*l1spone*epsc2 + 0.0904562210942245*l1sone*epsc2 + -0.07733429974180311*l4spone*epsc2 + -0.0904562210942245*l4sone*epsc2 + 0.07733429974180311*l1sone*xispone*epsc2 + -0.07733429974180311*l4sone*xispone*epsc2;
        }

        double fTs_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.3114246240345193 + 0.6475515597600012*as + 0.3114246240345193*epsb + -3.5132734729496087*chi2sone*epsb + -2.4746975917376997*chi2spone*epsb + 10.539820418848826*chi3spone*epsb + 3.71204638760655*chi3sppone*epsb + 0.3114246240345193*epsc + -3.5132734729496087*chi2sone*epsc + -2.4746975917376997*chi2spone*epsc + 10.539820418848826*chi3spone*epsc + 3.71204638760655*chi3sppone*epsc + -0.6228492480690386*epsb*etasone + -0.6228492480690386*epsc*etasone + -1.7566367364748043*epsb*etaspone + -1.7566367364748043*epsc*etaspone + -0.6186743979344249*epsb*etasppone + -0.6186743979344249*epsc*etasppone + 0.8783183682374022*xispone + 1.6452537827106624*as*xispone + 0.8783183682374022*epsb*xispone + -2.4746975917376997*chi2sone*epsb*xispone + 7.4240927752131*chi3spone*epsb*xispone + 0.8783183682374022*epsc*xispone + -2.4746975917376997*chi2sone*epsc*xispone + 7.4240927752131*chi3spone*epsc*xispone + -1.7566367364748043*epsb*etasone*xispone + -1.7566367364748043*epsc*etasone*xispone + -1.2373487958688498*epsb*etaspone*xispone + -1.2373487958688498*epsc*etaspone*xispone + 0.30933719896721246*xisppone + 0.5158529794955505*as*xisppone + 0.30933719896721246*epsb*xisppone + 0.30933719896721246*epsc*xisppone + -0.6186743979344249*epsb*etasone*xisppone + -0.6186743979344249*epsc*etasone*xisppone + 0.3114246240345193*l1sone*epsc2 + 0.723649768753796*l1spone*epsc2 + -0.3114246240345193*l4sone*epsc2 + -0.723649768753796*l4spone*epsc2 + 0.8783183682374022*l1sone*xispone*epsc2 + 0.6186743979344249*l1spone*xispone*epsc2 + -0.8783183682374022*l4sone*xispone*epsc2 + -0.6186743979344249*l4spone*xispone*epsc2 + 0.30933719896721246*l1sone*xisppone*epsc2 + -0.30933719896721246*l4sone*xisppone*epsc2;
        }
        // }}}

        // B_s -> D_s^* form factors
        // {{{
        double A1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0039029629021226243 + -0.0016234798992337527*as + 0.0039029629021226243*l2sone*epsc2;
        }

        double A1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.045389619854771975 + 0.0011297414860574153*as + 0.015611851608490495*epsb + -0.12489481286792396*chi2sone*epsb + 0.3746844386037719*chi3spone*epsb + 0.015611851608490495*epsc + -0.12489481286792396*chi3spone*epsc + -0.03122370321698099*epsb*etasone + 0.03122370321698099*xispone + -0.012987839193870022*as*xispone + 0.045389619854771975*l2sone*epsc2 + 0.03122370321698099*l2spone*epsc2 + -0.015611851608490495*l5sone*epsc2 + 0.03122370321698099*l2sone*xispone*epsc2;
        }

        double A1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.21375918823072088 + 0.09214542757833816*as + 0.15033477620210686*epsb + -1.7022574610885512*chi2sone*epsb + -0.9991585029433917*chi2spone*epsb + 5.106772383265654*chi3spone*epsb + 1.4987377544150877*chi3sppone*epsb + 0.15033477620210686*epsc + -1.7022574610885512*chi3spone*epsc + -0.49957925147169585*chi3sppone*epsc + -0.3006695524042137*epsb*etasone + -0.24978962573584793*epsb*etaspone + 0.4255643652721378*xispone + -0.01693774649928072*as*xispone + 0.12489481286792396*epsb*xispone + -0.9991585029433917*chi2sone*epsb*xispone + 2.9974755088301754*chi3spone*epsb*xispone + 0.12489481286792396*epsc*xispone + -0.9991585029433917*chi3spone*epsc*xispone + -0.24978962573584793*epsb*etasone*xispone + 0.12489481286792396*xisppone + -0.05195135677548009*as*xisppone + 0.21375918823072088*l2sone*epsc2 + 0.3631169588381758*l2spone*epsc2 + -0.15033477620210686*l5sone*epsc2 + -0.12489481286792396*l5spone*epsc2 + 0.4255643652721378*l2sone*xispone*epsc2 + 0.24978962573584793*l2spone*xispone*epsc2 + -0.12489481286792396*l5sone*xispone*epsc2 + 0.12489481286792396*l2sone*xisppone*epsc2;
        }

        double A5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.002527938212054363 + -0.0010515234135438955*as + 0.002527938212054363*l2sone*epsc2;
        }

        double A5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.024475416515784738 + -0.023236265467587044*epsb + -0.08089402278573964*chi2sone*epsb + 0.2426820683572189*chi3spone*epsb + 0.023236265467587044*epsc + 0.08089402278573964*chi2sone*epsc + -0.08089402278573964*chi3spone*epsc + 0.04647253093517409*epsb*etasone + 0.04647253093517409*epsc*etasone + as*(0.004718133656994875 - 0.008412187308351166*xispone) + 0.02022350569643491*xispone + (0.02022350569643491*l2spone + 0.02022350569643491*l3sone + 0.023236265467587044*l5sone - 0.04647253093517409*l6sone + l2sone*(0.024475416515784738 + 0.02022350569643491*xispone))*epsc2;
        }

        double A5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.08599311230540732 + -0.17850024445887427*epsb + -0.9450013740765909*chi2sone*epsb + -0.6471521822859172*chi2spone*epsb + 2.8350041222297726*chi3spone*epsb + 0.9707282734288758*chi3sppone*epsb + 0.1785002444588743*epsc + 0.9450013740765909*chi2sone*epsc + 0.6471521822859172*chi2spone*epsc + -0.9450013740765909*chi3spone*epsc + -0.3235760911429586*chi3sppone*epsc + 0.35700048891774855*epsb*etasone + 0.3570004889177486*epsc*etasone + 0.3717802474813927*epsb*etaspone + 0.3717802474813927*epsc*etaspone + 0.23625034351914773*xispone + -0.18589012374069636*epsb*xispone + -0.6471521822859172*chi2sone*epsb*xispone + 1.9414565468577516*chi3spone*epsb*xispone + 0.18589012374069636*epsc*xispone + 0.6471521822859172*chi2sone*epsc*xispone + -0.6471521822859172*chi3spone*epsc*xispone + 0.3717802474813927*epsb*etasone*xispone + 0.3717802474813927*epsc*etasone*xispone + as*(0.0719079842270774 + 0.02092069463925666*xispone - 0.03364874923340466*xisppone) + 0.08089402278573965*xisppone + (0.1958033321262779*l2spone + 0.23625034351914773*l3sone + 0.1617880455714793*l3spone + 0.17850024445887427*l5sone + 0.18589012374069636*l5spone - 0.5428906126584448*l6sone - 0.3717802474813927*l6spone + 0.1617880455714793*l2spone*xispone + 0.1617880455714793*l3sone*xispone + 0.18589012374069636*l5sone*xispone - 0.3717802474813927*l6sone*xispone + l2sone*(0.08599311230540732 + 0.23625034351914773*xispone + 0.08089402278573965*xisppone))*epsc2;
        }

        double V4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.004684646612240956 + 0.004297565736635373*as + 0.004684646612240956*epsb + 0.004684646612240956*epsc + -0.009369293224481911*epsb*etasone + 0.004684646612240956*l2sone*epsc2 + -0.004684646612240956*l5sone*epsc2;
        }

        double V4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04345453602144384 + 0.04722512441781111*as + 0.04345453602144384*epsb + -0.14990869159171058*chi2sone*epsb + 0.4497260747751317*chi3spone*epsb + 0.04345453602144384*epsc + -0.14990869159171058*chi3spone*epsc + -0.08690907204288768*epsb*etasone + -0.07495434579585529*epsb*etaspone + 0.037477172897927645*xispone + 0.034380525893082985*as*xispone + 0.037477172897927645*epsb*xispone + 0.037477172897927645*epsc*xispone + -0.07495434579585529*epsb*etasone*xispone + 0.04345453602144384*l2sone*epsc2 + 0.037477172897927645*l2spone*epsc2 + -0.04345453602144384*l5sone*epsc2 + -0.037477172897927645*l5spone*epsc2 + 0.037477172897927645*l2sone*xispone*epsc2 + -0.037477172897927645*l5sone*xispone*epsc2;
        }

        double V4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.13494178999195733 + 0.15007087537821776*as + 0.13494178999195733*epsb + -1.6903625358696244*chi2sone*epsb + -1.1992695327336846*chi2spone*epsb + 5.071087607608872*chi3spone*epsb + 1.7989042991005268*chi3sppone*epsb + 0.13494178999195733*epsc + -1.6903625358696244*chi3spone*epsc + -0.5996347663668423*chi3sppone*epsc + -0.26988357998391466*epsb*etasone + -0.8451812679348122*epsb*etaspone + -0.29981738318342116*epsb*etasppone + 0.4225906339674061*xispone + 0.44656204712865494*as*xispone + 0.4225906339674061*epsb*xispone + -1.1992695327336846*chi2sone*epsb*xispone + 3.5978085982010537*chi3spone*epsb*xispone + 0.4225906339674061*epsc*xispone + -1.1992695327336846*chi3spone*epsc*xispone + -0.8451812679348122*epsb*etasone*xispone + -0.5996347663668423*epsb*etaspone*xispone + 0.14990869159171058*xisppone + 0.13752210357233197*as*xisppone + 0.14990869159171058*epsb*xisppone + 0.14990869159171058*epsc*xisppone + -0.29981738318342116*epsb*etasone*xisppone + 0.13494178999195733*l2sone*epsc2 + 0.34763628817155073*l2spone*epsc2 + -0.13494178999195733*l5sone*epsc2 + -0.34763628817155073*l5spone*epsc2 + 0.4225906339674061*l2sone*xispone*epsc2 + 0.29981738318342116*l2spone*xispone*epsc2 + -0.4225906339674061*l5sone*xispone*epsc2 + -0.29981738318342116*l5spone*xispone*epsc2 + 0.14990869159171058*l2sone*xisppone*epsc2 + -0.14990869159171058*l5sone*xisppone*epsc2;
        }

        double P1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.017614604136229357 + -0.0012894515619258611*as + -0.007665367904889766*epsb + 0.007665367904889766*epsc + 0.015330735809779531*epsb*etasone + 0.015330735809779531*epsc*etasone + 0.017614604136229357*l2sone*epsc2 + 0.007665367904889766*l5sone*epsc2 + -0.015330735809779531*l6sone*epsc2;
        }

        double P1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.17385254509323594 + -0.027729362113336223*as + -0.07565561559229919*epsb + -0.5636673323593394*chi2sone*epsb + 1.6910019970780181*chi3spone*epsb + 0.07565561559229919*epsc + 0.5636673323593394*chi2sone*epsc + -0.5636673323593394*chi3spone*epsc + 0.15131123118459838*epsb*etasone + 0.15131123118459838*epsc*etasone + 0.12264588647823625*epsb*etaspone + 0.12264588647823625*epsc*etaspone + 0.14091683308983485*xispone + -0.01031561249540689*as*xispone + -0.061322943239118126*epsb*xispone + 0.061322943239118126*epsc*xispone + 0.12264588647823625*epsb*etasone*xispone + 0.12264588647823625*epsc*etasone*xispone + 0.17385254509323594*l2sone*epsc2 + 0.14091683308983485*l2spone*epsc2 + 0.14091683308983485*l3sone*epsc2 + 0.07565561559229919*l5sone*epsc2 + 0.061322943239118126*l5spone*epsc2 + -0.21263417442371652*l6sone*epsc2 + -0.12264588647823625*l6spone*epsc2 + 0.14091683308983485*l2sone*xispone*epsc2 + 0.061322943239118126*l5sone*xispone*epsc2 + -0.12264588647823625*l6sone*xispone*epsc2;
        }

        double P1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.6199338182961519 + -0.37274031827421483*as + -0.269777325517593*epsb + -6.69061610770223*chi2sone*epsb + -4.509338658874715*chi2spone*epsb + 20.071848323106686*chi3spone*epsb + 6.7640079883120725*chi3sppone*epsb + 0.269777325517593*epsc + 6.690616107702229*chi2sone*epsc + 4.509338658874715*chi2spone*epsc + -6.69061610770223*chi3spone*epsc + -2.2546693294373576*chi3sppone*epsc + 0.539554651035186*epsb*etasone + 0.539554651035186*epsc*etasone + 1.45578162243326*epsb*etaspone + 1.4557816224332596*epsc*etaspone + 0.49058354591294506*epsb*etasppone + 0.49058354591294506*epsc*etasppone + 1.6726540269255574*xispone + -0.24246612189750352*as*xispone + -0.72789081121663*epsb*xispone + -4.509338658874715*chi2sone*epsb*xispone + 13.528015976624145*chi3spone*epsb*xispone + 0.7278908112166298*epsc*xispone + 4.509338658874715*chi2sone*epsc*xispone + -4.509338658874715*chi3spone*epsc*xispone + 1.45578162243326*epsb*etasone*xispone + 1.4557816224332598*epsc*etasone*xispone + 0.9811670918258901*epsb*etaspone*xispone + 0.9811670918258901*epsc*etaspone*xispone + 0.5636673323593394*xisppone + -0.04126244998162756*as*xisppone + -0.24529177295647253*epsb*xisppone + 0.24529177295647253*epsc*xisppone + 0.49058354591294506*epsb*etasone*xisppone + 0.49058354591294506*epsc*etasone*xisppone + 0.6199338182961519*l2sone*epsc2 + 1.3908203607458876*l2spone*epsc2 + 1.6726540269255572*l3sone*epsc2 + 1.1273346647186788*l3spone*epsc2 + 0.269777325517593*l5sone*epsc2 + 0.6052449247383936*l5spone*epsc2 + -1.2674454622518159*l6sone*epsc2 + -1.7010733953897321*l6spone*epsc2 + 1.6726540269255574*l2sone*xispone*epsc2 + 1.1273346647186788*l2spone*xispone*epsc2 + 1.1273346647186788*l3sone*xispone*epsc2 + 0.72789081121663*l5sone*xispone*epsc2 + 0.49058354591294506*l5spone*xispone*epsc2 + -1.9463651683462044*l6sone*xispone*epsc2 + -0.9811670918258901*l6spone*xispone*epsc2 + 0.5636673323593394*l2sone*xisppone*epsc2 + 0.24529177295647253*l5sone*xisppone*epsc2 + -0.49058354591294506*l6sone*xisppone*epsc2;
        }

        double T1s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.005303371401576171 + 0.002984624789048728*as + -0.0023078743419353625*epsb + 0.0023078743419353625*epsc + 0.004615748683870725*epsb*etasone + 0.005303371401576171*l2sone*epsc2 + -0.0023078743419353625*l5sone*epsc2;
        }

        double T1s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.049193794682977326 + 0.031688272203984034*as + -0.021407721227583158*epsb + -0.1697078848504375*chi2sone*epsb + 0.5091236545513125*chi3spone*epsb + 0.021407721227583158*epsc + -0.1697078848504375*chi3spone*epsc + 0.042815442455166317*epsb*etasone + 0.0369259894709658*epsb*etaspone + 0.042426971212609375*xispone + 0.023876998312389822*as*xispone + -0.0184629947354829*epsb*xispone + 0.0184629947354829*epsc*xispone + 0.0369259894709658*epsb*etasone*xispone + 0.049193794682977326*l2sone*epsc2 + 0.042426971212609375*l2spone*epsc2 + -0.021407721227583158*l5sone*epsc2 + -0.0184629947354829*l5spone*epsc2 + 0.042426971212609375*l2sone*xispone*epsc2 + -0.0184629947354829*l5sone*xispone*epsc2;
        }

        double T1s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.152764229440672 + 0.07358933273759395*as + -0.0664785885798743*epsb + -1.9136171995561493*chi2sone*epsb + -1.3576630788034998*chi2spone*epsb + 5.740851598668448*chi3spone*epsb + 2.03649461820525*chi3sppone*epsb + 0.0664785885798743*epsc + -1.9136171995561493*chi3spone*epsc + -0.6788315394017499*chi3sppone*epsc + 0.1329571771597486*epsb*etasone + 0.41637551858326205*epsb*etaspone + 0.1477039578838632*epsb*etasppone + 0.4784042998890373*xispone + 0.30126017425665197*as*xispone + -0.20818775929163102*epsb*xispone + -1.3576630788034998*chi2sone*epsb*xispone + 4.0729892364105*chi3spone*epsb*xispone + 0.20818775929163102*epsc*xispone + -1.3576630788034998*chi3spone*epsc*xispone + 0.41637551858326205*epsb*etasone*xispone + 0.2954079157677264*epsb*etaspone*xispone + 0.16970788485043747*xisppone + 0.09550799324955929*as*xisppone + -0.0738519789419316*epsb*xisppone + 0.0738519789419316*epsc*xisppone + 0.1477039578838632*epsb*etasone*xisppone + 0.152764229440672*l2sone*epsc2 + 0.3935503574638186*l2spone*epsc2 + -0.0664785885798743*l5sone*epsc2 + -0.17126176982066524*l5spone*epsc2 + 0.4784042998890373*l2sone*xispone*epsc2 + 0.33941576970087495*l2spone*xispone*epsc2 + -0.20818775929163102*l5sone*xispone*epsc2 + -0.1477039578838632*l5spone*xispone*epsc2 + 0.16970788485043747*l2sone*xisppone*epsc2 + -0.0738519789419316*l5sone*xisppone*epsc2;
        }

        double T2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.00106406613984606 + 0.0003556896380373649*as + 0.00106406613984606*l2sone*epsc2;
        }

        double T2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.011366320304480485 + 0.009527969427653715*as + -0.009780667574323585*epsb + -0.034050116475073916*chi2sone*epsb + 0.10215034942522176*chi3spone*epsb + 0.009780667574323585*epsc + -0.034050116475073916*chi3spone*epsc + 0.01956133514864717*epsb*etasone + 0.008512529118768479*xispone + 0.002845517104298919*as*xispone + 0.011366320304480485*l2sone*epsc2 + 0.008512529118768479*l2spone*epsc2 + -0.009780667574323585*l5sone*epsc2 + 0.008512529118768479*l2sone*xispone*epsc2;
        }

        double T2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04649869228477457 + 0.061016848241349474*as + -0.08491544150876175*epsb + -0.43182248269352336*chi2sone*epsb + -0.2724009318005914*chi2spone*epsb + 1.2954674480805701*chi3spone*epsb + 0.4086013977008871*chi3sppone*epsb + 0.08491544150876175*epsc + -0.43182248269352336*chi3spone*epsc + -0.1362004659002957*chi3sppone*epsc + 0.1698308830175235*epsb*etasone + 0.15649068118917736*epsb*etaspone + 0.10795562067338084*xispone + 0.08191478962982757*as*xispone + -0.07824534059458868*epsb*xispone + -0.2724009318005914*chi2sone*epsb*xispone + 0.8172027954017742*chi3spone*epsb*xispone + 0.07824534059458868*epsc*xispone + -0.2724009318005914*chi3spone*epsc*xispone + 0.15649068118917736*epsb*etasone*xispone + 0.03405011647507392*xisppone + 0.011382068417195676*as*xisppone + 0.04649869228477457*l2sone*epsc2 + 0.09093056243584388*l2spone*epsc2 + -0.08491544150876175*l5sone*epsc2 + -0.07824534059458868*l5spone*epsc2 + 0.10795562067338084*l2sone*xispone*epsc2 + 0.06810023295014785*l2spone*xispone*epsc2 + -0.07824534059458868*l5sone*xispone*epsc2 + 0.03405011647507392*l2sone*xisppone*epsc2;
        }

        double T23s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0032856900136407955 + 0.0010983207226422545*as + 0.0032856900136407955*l2sone*epsc2;
        }

        double T23s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.03821102696074016 + 0.032126984281562174*as + 0.013142760054563206*epsb + -0.10514208043650546*chi2sone*epsb + 0.31542624130951635*chi3spone*epsb + 0.01314276005456317*epsc + 0.10514208043650554*chi2sone*epsc + -0.10514208043650546*chi3spone*epsc + -0.026285520109126412*epsb*etasone + 0.02628552010912638*epsc*etasone + 0.026285520109126364*xispone + 0.008786565781138036*as*xispone + 0.03821102696074016*l2sone*epsc2 + 0.026285520109126364*l2spone*epsc2 + 0.026285520109126385*l3sone*epsc2 + 0.013142760054563206*l5sone*epsc2 + -0.026285520109126385*l6sone*epsc2 + 0.026285520109126364*l2sone*xispone*epsc2;
        }

        double T23s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.1799521152793105 + 0.24300691099939695*as + 0.12655858773383452*epsb + -1.4330370236166963*chi2sone*epsb + -0.8411366434920439*chi2spone*epsb + 4.299111070850089*chi3spone*epsb + 1.2617049652380659*chi3sppone*epsb + 0.12655858773383424*epsc + 1.433037023616697*chi2sone*epsc + 0.8411366434920442*chi2spone*epsc + -1.4330370236166963*chi3spone*epsc + -0.42056832174602193*chi3sppone*epsc + -0.25311717546766904*epsb*etasone + 0.2531171754676687*epsc*etasone + -0.21028416087301116*epsb*etaspone + 0.21028416087301105*epsc*etaspone + 0.35825925590417407*xispone + 0.27458900581477363*as*xispone + 0.10514208043650558*epsb*xispone + -0.8411366434920439*chi2sone*epsb*xispone + 2.5234099304761317*chi3spone*epsb*xispone + 0.10514208043650543*epsc*xispone + 0.8411366434920442*chi2sone*epsc*xispone + -0.8411366434920439*chi3spone*epsc*xispone + -0.21028416087301116*epsb*etasone*xispone + 0.21028416087301105*epsc*etasone*xispone + 0.10514208043650548*xisppone + 0.03514626312455215*as*xisppone + 0.1799521152793105*l2sone*epsc2 + 0.3056882156859214*l2spone*epsc2 + 0.35825925590417423*l3sone*epsc2 + 0.21028416087301105*l3spone*epsc2 + 0.12655858773383452*l5sone*epsc2 + 0.10514208043650558*l5spone*epsc2 + -0.35825925590417423*l6sone*epsc2 + -0.21028416087301105*l6spone*epsc2 + 0.35825925590417407*l2sone*xispone*epsc2 + 0.21028416087301097*l2spone*xispone*epsc2 + 0.21028416087301105*l3sone*xispone*epsc2 + 0.10514208043650558*l5sone*xispone*epsc2 + -0.21028416087301105*l6sone*xispone*epsc2 + 0.10514208043650548*l2sone*xisppone*epsc2;
        }
        // }}}

        // B_s^* -> D_s form factors
        // {{{
        double P2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.03294794420017363 + -0.0015899447001142385*as + -0.01538157363810894*epsb + 0.01538157363810894*epsc + -0.03076314727621788*epsb*etasone + -0.03076314727621788*epsc*etasone + 0.03294794420017363*l1sone*epsc2 + -0.01538157363810894*l4sone*epsc2;
        }

        double P2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.27812193901117405 + -0.040298078998178166*as + -0.12983975750606846*epsb + 1.054334214405556*chi2sone*epsb + -1.054334214405556*chi3spone*epsb + 0.12983975750606844*epsc + -1.054334214405556*chi2sone*epsc + 3.1630026432166685*chi3spone*epsc + -0.2596795150121369*epsb*etasone + -0.2596795150121369*epsc*etasone + -0.24610517820974304*epsb*etaspone + -0.24610517820974304*epsc*etaspone + 0.263583553601389*xispone + -0.012719557600913908*as*xispone + -0.12305258910487152*epsb*xispone + 0.12305258910487152*epsc*xispone + -0.24610517820974304*epsb*etasone*xispone + -0.24610517820974304*epsc*etasone*xispone + 0.27812193901117405*l1sone*epsc2 + 0.263583553601389*l1spone*epsc2 + -0.12983975750606844*l4sone*epsc2 + -0.12305258910487152*l4spone*epsc2 + 0.263583553601389*l1sone*xispone*epsc2 + -0.12305258910487152*l4sone*xispone*epsc2;
        }

        double P2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.8002854055840344 + -0.600240705338798*as + -0.3736090125292201*epsb + 11.008570477168684*chi2sone*epsb + 8.43467371524445*chi2spone*epsb + -11.008570477168682*chi3spone*epsb + -4.217336857622223*chi3sppone*epsb + 0.37360901252921974*epsc + -11.008570477168682*chi2sone*epsc + -8.434673715244447*chi2spone*epsc + 33.02571143150605*chi3spone*epsc + 12.652010572866674*chi3sppone*epsc + -0.7472180250584398*epsb*etasone + -0.7472180250584395*epsc*etasone + -2.5696464765165814*epsb*etaspone + -2.569646476516581*epsc*etaspone + -0.9844207128389721*epsb*etasppone + -0.9844207128389721*epsc*etasppone + 2.7521426192921705*xispone + -0.3478237471872532*as*xispone + -1.2848232382582907*epsb*xispone + 8.43467371524445*chi2sone*epsb*xispone + -8.434673715244447*chi3spone*epsb*xispone + 1.2848232382582905*epsc*xispone + -8.434673715244447*chi2sone*epsc*xispone + 25.304021145733348*chi3spone*epsc*xispone + -2.5696464765165814*epsb*etasone*xispone + -2.569646476516581*epsc*etasone*xispone + -1.9688414256779443*epsb*etaspone*xispone + -1.9688414256779443*epsc*etaspone*xispone + 1.0543342144055559*xisppone + -0.050878230403655675*as*xisppone + -0.4922103564194861*epsb*xisppone + 0.4922103564194861*epsc*xisppone + -0.9844207128389721*epsb*etasone*xisppone + -0.9844207128389721*epsc*etasone*xisppone + 0.8002854055840344*l1sone*epsc2 + 2.224975512089393*l1spone*epsc2 + -0.37360901252921974*l4sone*epsc2 + -1.0387180600485475*l4spone*epsc2 + 2.7521426192921705*l1sone*xispone*epsc2 + 2.1086684288111117*l1spone*xispone*epsc2 + -1.2848232382582905*l4sone*xispone*epsc2 + -0.9844207128389721*l4spone*xispone*epsc2 + 1.0543342144055559*l1sone*xisppone*epsc2 + -0.4922103564194861*l4sone*xisppone*epsc2;
        }

        double V5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.005283133663761428 + 0.004846601268945023*as + 0.005283133663761428*epsb + 0.005283133663761428*epsc + -0.010566267327522857*epsc*etasone + 0.005283133663761428*l1sone*epsc2 + -0.005283133663761428*l4sone*epsc2;
        }

        double V5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.046629963052817115 + 0.05107859492650357*as + 0.046629963052817115*epsb + -0.1690602772403657*chi3spone*epsb + 0.046629963052817115*epsc + -0.1690602772403657*chi2sone*epsc + 0.5071808317210972*chi3spone*epsc + -0.09325992610563423*epsc*etasone + -0.08453013862018285*epsc*etaspone + 0.04226506931009143*xispone + 0.038772810151560186*as*xispone + 0.04226506931009143*epsb*xispone + 0.04226506931009143*epsc*xispone + -0.08453013862018285*epsc*etasone*xispone + 0.046629963052817115*l1sone*epsc2 + 0.04226506931009143*l1spone*epsc2 + -0.046629963052817115*l4sone*epsc2 + -0.04226506931009143*l4spone*epsc2 + 0.04226506931009143*l1sone*xispone*epsc2 + -0.04226506931009143*l4sone*xispone*epsc2;
        }

        double V5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.13687352617980045 + 0.15146661888422702*as + 0.13687352617980045*epsb + -1.8302793721708792*chi3spone*epsb + -0.6762411089614629*chi3sppone*epsb + 0.13687352617980045*epsc + -1.8302793721708792*chi2sone*epsc + -1.3524822179229259*chi2spone*epsc + 5.490838116512637*chi3spone*epsc + 2.0287233268843887*chi3sppone*epsc + -0.2737470523596009*epsc*etasone + -0.9151396860854396*epsc*etaspone + -0.33812055448073147*epsc*etasppone + 0.4575698430427198*xispone + 0.4861743797151489*as*xispone + 0.4575698430427198*epsb*xispone + -1.3524822179229259*chi3spone*epsb*xispone + 0.4575698430427198*epsc*xispone + -1.3524822179229259*chi2sone*epsc*xispone + 4.0574466537687774*chi3spone*epsc*xispone + -0.9151396860854396*epsc*etasone*xispone + -0.6762411089614629*epsc*etaspone*xispone + 0.16906027724036574*xisppone + 0.15509124060624077*as*xisppone + 0.16906027724036574*epsb*xisppone + 0.16906027724036574*epsc*xisppone + -0.33812055448073147*epsc*etasone*xisppone + 0.13687352617980045*l1sone*epsc2 + 0.3730397044225369*l1spone*epsc2 + -0.13687352617980045*l4sone*epsc2 + -0.3730397044225369*l4spone*epsc2 + 0.4575698430427198*l1sone*xispone*epsc2 + 0.33812055448073147*l1spone*xispone*epsc2 + -0.4575698430427198*l4sone*xispone*epsc2 + -0.33812055448073147*l4spone*xispone*epsc2 + 0.16906027724036574*l1sone*xisppone*epsc2 + -0.16906027724036574*l4sone*xisppone*epsc2;
        }

        double A2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.00503545648285162 + -0.0020945529302699613*as + 0.00503545648285162*l1sone*epsc2;
        }

        double A2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.05481409325271715 + 0.0030156952718291374*as + 0.02014182593140648*epsb + -0.16113460745125183*chi3spone*epsb + 0.02014182593140648*epsc + -0.16113460745125183*chi2sone*epsc + 0.4834038223537555*chi3spone*epsc + -0.04028365186281296*epsc*etasone + 0.04028365186281296*xispone + -0.016756423442159687*as*xispone + 0.05481409325271715*l1sone*epsc2 + 0.04028365186281296*l1spone*epsc2 + -0.02014182593140648*l4sone*epsc2 + 0.04028365186281296*l1sone*xispone*epsc2;
        }

        double A2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.24131157329136588 + 0.11401702767769033*as + 0.17897272114805565*epsb + -2.076320198989453*chi3spone*epsb + -0.6445384298050074*chi3sppone*epsb + 0.17897272114805565*epsc + -2.076320198989453*chi2sone*epsc + -1.2890768596100148*chi2spone*epsc + 6.228960596968358*chi3spone*epsc + 1.9336152894150223*chi3sppone*epsc + -0.3579454422961113*epsc*etasone + -0.3222692149025037*epsc*etaspone + 0.5190800497473632*xispone + -0.009387284709686276*as*xispone + 0.16113460745125185*epsb*xispone + -1.2890768596100148*chi3spone*epsb*xispone + 0.16113460745125185*epsc*xispone + -1.2890768596100148*chi2sone*epsc*xispone + 3.8672305788300445*chi3spone*epsc*xispone + -0.3222692149025037*epsc*etasone*xispone + 0.16113460745125185*xisppone + -0.06702569376863876*as*xisppone + 0.24131157329136588*l1sone*epsc2 + 0.43851274602173723*l1spone*epsc2 + -0.17897272114805565*l4sone*epsc2 + -0.16113460745125185*l4spone*epsc2 + 0.5190800497473632*l1sone*xispone*epsc2 + 0.3222692149025037*l1spone*xispone*epsc2 + -0.16113460745125185*l4sone*xispone*epsc2 + 0.16113460745125185*l1sone*xisppone*epsc2;
        }

        double A6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0035285573463436095 + -0.0014677418332536486*as + 0.0035285573463436095*l1sone*epsc2;
        }

        double A6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.031570021713464555 + 0.005931223319652669*as + -0.03023324226499239*epsb + 0.1129138350829955*chi2sone*epsb + -0.1129138350829955*chi3spone*epsb + 0.03023324226499239*epsc + -0.1129138350829955*chi2sone*epsc + 0.33874150524898655*chi3spone*epsc + -0.06046648452998478*epsb*etasone + -0.06046648452998478*epsc*etasone + 0.028228458770748876*xispone + -0.01174193466602919*as*xispone + 0.031570021713464555*l1sone*epsc2 + 0.028228458770748876*l1spone*epsc2 + -0.03023324226499239*l4sone*epsc2 + 0.028228458770748876*l1sone*xispone*epsc2;
        }

        double A6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.1012710007502295 + 0.07834959986071738*as + -0.21003049797287487*epsb + 1.2360683649968567*chi2sone*epsb + 0.903310680663964*chi2spone*epsb + -1.2360683649968567*chi3spone*epsb + -0.451655340331982*chi3sppone*epsb + 0.21003049797287487*epsc + -1.2360683649968567*chi2sone*epsc + -0.903310680663964*chi2spone*epsc + 3.708205094990571*chi3spone*epsc + 1.354966020995946*chi3sppone*epsc + -0.42006099594574975*epsb*etasone + -0.42006099594574975*epsc*etasone + -0.4837318762398783*epsb*etaspone + -0.4837318762398782*epsc*etaspone + 0.30901709124921417*xispone + 0.02396591722516295*as*xispone + -0.24186593811993914*epsb*xispone + 0.903310680663964*chi2sone*epsb*xispone + -0.903310680663964*chi3spone*epsb*xispone + 0.2418659381199391*epsc*xispone + -0.903310680663964*chi2sone*epsc*xispone + 2.709932041991892*chi3spone*epsc*xispone + -0.4837318762398783*epsb*etasone*xispone + -0.4837318762398782*epsc*etasone*xispone + 0.1129138350829955*xisppone + -0.046967738664116764*as*xisppone + 0.1012710007502295*l1sone*epsc2 + 0.25256017370771644*l1spone*epsc2 + -0.21003049797287487*l4sone*epsc2 + -0.2418659381199391*l4spone*epsc2 + 0.30901709124921417*l1sone*xispone*epsc2 + 0.225827670165991*l1spone*xispone*epsc2 + -0.2418659381199391*l4sone*xispone*epsc2 + 0.1129138350829955*l1sone*xisppone*epsc2;
        }

        double T1bars_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.006279445375704966 + 0.003638375183292366*as + -0.0029315258902368664*epsb + 0.0029315258902368664*epsc + -0.005863051780473733*epsc*etasone + 0.006279445375704966*l1sone*epsc2 + -0.0029315258902368664*l4sone*epsc2;
        }

        double T1bars_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.055423603583944527 + 0.03694318529992362*as + -0.02587421644993889*epsb + -0.20094225202255891*chi3spone*epsb + 0.02587421644993889*epsc + -0.20094225202255891*chi2sone*epsc + 0.6028267560676768*chi3spone*epsc + -0.05174843289987778*epsc*etasone + -0.046904414243789855*epsc*etaspone + 0.05023556300563973*xispone + 0.02910700146633893*as*xispone + -0.023452207121894927*epsb*xispone + 0.023452207121894927*epsc*xispone + -0.046904414243789855*epsc*etasone*xispone + 0.055423603583944527*l1sone*epsc2 + 0.05023556300563973*l1spone*epsc2 + -0.02587421644993889*l4sone*epsc2 + -0.023452207121894927*l4spone*epsc2 + 0.05023556300563973*l1sone*xispone*epsc2 + -0.023452207121894927*l4sone*xispone*epsc2;
        }

        double T1bars_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.16268561155695807 + 0.07824669591820689*as + -0.07594891805149258*epsb + -2.1754398187313426*chi3spone*epsb + -0.8037690080902358*chi3sppone*epsb + 0.07594891805149258*epsc + -2.1754398187313426*chi2sone*epsc + -1.6075380161804715*chi2spone*epsc + 6.526319456194027*chi3spone*epsc + 2.4113070242707066*chi3sppone*epsc + -0.15189783610298516*epsc*etasone + -0.507796291686602*epsc*etaspone + -0.18761765697515942*epsc*etasppone + 0.5438599546828357*xispone + 0.3537594853320668*as*xispone + -0.253898145843301*epsb*xispone + -1.6075380161804715*chi3spone*epsb*xispone + 0.253898145843301*epsc*xispone + -1.6075380161804715*chi2sone*epsc*xispone + 4.822614048541413*chi3spone*epsc*xispone + -0.507796291686602*epsc*etasone*xispone + -0.37523531395031884*epsc*etaspone*xispone + 0.20094225202255894*xisppone + 0.11642800586535572*as*xisppone + -0.09380882848757971*epsb*xisppone + 0.09380882848757971*epsc*xisppone + -0.18761765697515942*epsc*etasone*xisppone + 0.16268561155695807*l1sone*epsc2 + 0.4433888286715562*l1spone*epsc2 + -0.07594891805149258*l4sone*epsc2 + -0.20699373159951115*l4spone*epsc2 + 0.5438599546828357*l1sone*xispone*epsc2 + 0.4018845040451179*l1spone*xispone*epsc2 + -0.253898145843301*l4sone*xispone*epsc2 + -0.18761765697515942*l4spone*xispone*epsc2 + 0.20094225202255894*l1sone*xisppone*epsc2 + -0.09380882848757971*l4sone*xisppone*epsc2;
        }

        double T2bars_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return -0.0015316603889868848 + -0.0005119942350893066*as + -0.0015316603889868848*l1sone*epsc2;
        }

        double T2bars_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return -0.01523543419560388 + -0.012837094303612273*as + 0.013123510563295124*epsb + 0.049013132447580314*chi3spone*epsb + -0.013123510563295124*epsc + 0.049013132447580314*chi2sone*epsc + -0.14703939734274096*chi3spone*epsc + 0.026247021126590248*epsc*etasone + -0.012253283111895079*xispone + -0.004095953880714453*as*xispone + -0.01523543419560388*l1sone*epsc2 + -0.012253283111895079*l1spone*epsc2 + 0.013123510563295124*l4sone*epsc2 + -0.012253283111895079*l1sone*xispone*epsc2;
        }

        double T2bars_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return -0.05766303680766183 + -0.0742514207574954*as + 0.10429260961664973*epsb + 0.5855601591544847*chi3spone*epsb + 0.19605252979032126*chi3sppone*epsb + -0.10429260961664973*epsc + 0.5855601591544847*chi2sone*epsc + 0.3921050595806425*chi2spone*epsc + -1.7566804774634543*chi3spone*epsc + -0.5881575893709637*chi3sppone*epsc + 0.20858521923329945*epsc*etasone + 0.20997616901272195*epsc*etaspone + -0.14639003978862117*xispone + -0.1108886621903271*as*xispone + 0.10498808450636098*epsb*xispone + 0.3921050595806425*chi3spone*epsb*xispone + -0.10498808450636098*epsc*xispone + 0.3921050595806425*chi2sone*epsc*xispone + -1.1763151787419275*chi3spone*epsc*xispone + 0.20997616901272195*epsc*etasone*xispone + -0.049013132447580314*xisppone + -0.01638381552285781*as*xisppone + -0.05766303680766183*l1sone*epsc2 + -0.12188347356483102*l1spone*epsc2 + 0.10429260961664973*l4sone*epsc2 + 0.10498808450636098*l4spone*epsc2 + -0.14639003978862117*l1sone*xispone*epsc2 + -0.09802626489516063*l1spone*xispone*epsc2 + 0.10498808450636098*l4sone*xispone*epsc2 + -0.049013132447580314*l1sone*xisppone*epsc2;
        }

        double T23bars_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return -0.004371536851026707 + -0.0014612910814299973*as + -0.004371536851026707*l1sone*epsc2;
        }

        double T23bars_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return -0.047586912810369 + -0.04165717237084313*as + -0.01748614740410682*epsb + -0.13988917923285463*chi2sone*epsb + 0.13988917923285463*chi3spone*epsb + -0.01748614740410684*epsc + 0.13988917923285463*chi2sone*epsc + -0.41966753769856396*chi3spone*epsc + -0.03497229480821366*epsb*etasone + 0.03497229480821368*epsc*etasone + -0.03497229480821366*xispone + -0.011690328651439978*as*xispone + -0.047586912810369*l1sone*epsc2 + -0.03497229480821366*l1spone*epsc2 + 0.017486147404106835*l4sone*epsc2 + -0.03497229480821366*l1sone*xispone*epsc2;
        }

        double T23bars_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return -0.2094948966027084 + -0.2941557148897957*as + -0.1553753564332622*epsb + -1.8025595683975166*chi2sone*epsb + -1.1191134338628366*chi2spone*epsb + 1.802559568397517*chi3spone*epsb + 0.5595567169314184*chi3sppone*epsb + -0.1553753564332624*epsc + 1.802559568397517*chi2sone*epsc + 1.1191134338628368*chi2spone*epsc + -5.407678705192552*chi3spone*epsc + -1.6786701507942554*chi3sppone*epsc + -0.31075071286652456*epsb*etasone + 0.3107507128665248*epsc*etasone + -0.27977835846570914*epsb*etaspone + 0.27977835846570936*epsc*etaspone + -0.45063989209937927*xispone + -0.35663803626962504*as*xispone + -0.13988917923285446*epsb*xispone + -1.1191134338628366*chi2sone*epsb*xispone + 1.1191134338628368*chi3spone*epsb*xispone + -0.13988917923285468*epsc*xispone + 1.1191134338628368*chi2sone*epsc*xispone + -3.357340301588511*chi3spone*epsc*xispone + -0.27977835846570914*epsb*etasone*xispone + 0.27977835846570936*epsc*etasone*xispone + -0.1398891792328546*xisppone + -0.0467613146057599*as*xisppone + -0.2094948966027084*l1sone*epsc2 + -0.38069530248295197*l1spone*epsc2 + 0.1553753564332624*l4sone*epsc2 + 0.13988917923285468*l4spone*epsc2 + -0.45063989209937927*l1sone*xispone*epsc2 + -0.2797783584657092*l1spone*xispone*epsc2 + 0.13988917923285468*l4sone*xispone*epsc2 + -0.1398891792328546*l1sone*xisppone*epsc2;
        }
        // }}}

        // B_s^* -> D_s^* form factors
        // {{{
        double S2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.02706009741680484 + 0.006784122900199494*as + 0.02706009741680484*l2sone*epsc2;
        }

        double S2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.2189743421429717 + 0.14491175989743943*as + -0.24666697664657322*epsb + -0.8659231173377548*chi3spone*epsb + 0.24666697664657322*epsc + -0.8659231173377548*chi3spone*epsc + 0.2164807793344387*xispone + 0.05427298320159596*as*xispone + 0.2189743421429717*l2sone*epsc2 + 0.2164807793344387*l2spone*epsc2 + -0.24666697664657322*l5sone*epsc2 + 0.2164807793344387*l2sone*xispone*epsc2;
        }

        double S2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.6697186753961397 + 0.46925290063727393*as + -1.502731992354025*epsb + -8.739025183250604*chi3spone*epsb + -3.4636924693510194*chi3sppone*epsb + 1.502731992354025*epsc + -8.739025183250604*chi3spone*epsc + -3.4636924693510194*chi3sppone*epsc + 2.184756295812651*xispone + 1.2678400455827072*as*xispone + -1.9733358131725856*epsb*xispone + -6.927384938702039*chi3spone*epsb*xispone + 1.9733358131725856*epsc*xispone + -6.927384938702039*chi3spone*epsc*xispone + 0.8659231173377548*xisppone + 0.2170919328063838*as*xisppone + 0.6697186753961397*l2sone*epsc2 + 1.7517947371437732*l2spone*epsc2 + -1.502731992354025*l5sone*epsc2 + -1.9733358131725856*l5spone*epsc2 + 2.184756295812651*l2sone*xispone*epsc2 + 1.7318462346755097*l2spone*xispone*epsc2 + -1.9733358131725856*l5sone*xispone*epsc2 + 0.8659231173377548*l2sone*xisppone*epsc2;
        }

        double S3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.01913437838299128 + 0.00479709930713401*as + 0.01913437838299128*l2sone*epsc2;
        }

        double S3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.15483824223515844 + 0.1024680880971562*as + -0.17441989188157567*epsb + 0.6123001082557211*chi2sone*epsb + -0.6123001082557209*chi3spone*epsb + 0.17441989188157567*epsc + 0.6123001082557211*chi2sone*epsc + -0.6123001082557209*chi3spone*epsc + -0.34883978376315133*epsb*etasone + 0.34883978376315133*epsc*etasone + 0.15307502706393022*xispone + 0.03837679445707208*as*xispone + 0.15483824223515844*l2sone*epsc2 + 0.15307502706393022*l2spone*epsc2 + 0.15307502706393028*l3sone*epsc2 + 0.17441989188157567*l5sone*epsc2 + -0.34883978376315133*l6sone*epsc2 + 0.15307502706393022*l2sone*xispone*epsc2;
        }

        double S3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.4735626168598823 + 0.3318119081320737*as + -1.0625919820995022*epsb + 6.179423968036516*chi2sone*epsb + 4.898400866045769*chi2spone*epsb + -6.179423968036511*chi3spone*epsb + -2.449200433022884*chi3sppone*epsb + 1.0625919820995022*epsc + 6.179423968036516*chi2sone*epsc + 4.898400866045769*chi2spone*epsc + -6.179423968036511*chi3spone*epsc + -2.449200433022884*chi3sppone*epsc + -2.1251839641990045*epsb*etasone + 2.1251839641990045*epsc*etasone + -2.7907182701052107*epsb*etaspone + 2.7907182701052107*epsc*etaspone + 1.5448559920091278*xispone + 0.8964982936913937*as*xispone + -1.3953591350526053*epsb*xispone + 4.898400866045769*chi2sone*epsb*xispone + -4.898400866045768*chi3spone*epsb*xispone + 1.3953591350526053*epsc*xispone + 4.898400866045769*chi2sone*epsc*xispone + -4.898400866045768*chi3spone*epsc*xispone + -2.7907182701052107*epsb*etasone*xispone + 2.7907182701052107*epsc*etasone*xispone + 0.612300108255721*xisppone + 0.15350717782828832*as*xisppone + 0.4735626168598823*l2sone*epsc2 + 1.2387059378812677*l2spone*epsc2 + 1.544855992009129*l3sone*epsc2 + 1.2246002165114422*l3spone*epsc2 + 1.0625919820995022*l5sone*epsc2 + 1.3953591350526053*l5spone*epsc2 + -3.52054309925161*l6sone*epsc2 + -2.7907182701052107*l6spone*epsc2 + 1.5448559920091278*l2sone*xispone*epsc2 + 1.224600216511442*l2spone*xispone*epsc2 + 1.2246002165114422*l3sone*xispone*epsc2 + 1.3953591350526053*l5sone*xispone*epsc2 + -2.7907182701052107*l6sone*xispone*epsc2 + 0.612300108255721*l2sone*xisppone*epsc2;
        }

        double P3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.02181039634593078 + -0.0015340565139036527*as + -0.009570660132031263*epsb + 0.009570660132031263*epsc + 0.02181039634593078*l2sone*epsc2 + -0.009570660132031263*l5sone*epsc2;
        }

        double P3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.22407038485915534 + -0.03424640011175238*as + -0.09832473766761804*epsb + -0.697932683069785*chi3spone*epsb + 0.09832473766761804*epsc + -0.697932683069785*chi3spone*epsc + 0.17448317076744624*xispone + -0.012272452111229222*as*xispone + -0.07656528105625009*epsb*xispone + 0.07656528105625009*epsc*xispone + 0.22407038485915534*l2sone*epsc2 + 0.17448317076744624*l2spone*epsc2 + -0.09832473766761804*l5sone*epsc2 + -0.07656528105625009*l5spone*epsc2 + 0.17448317076744624*l2sone*xispone*epsc2 + -0.07656528105625009*l5sone*xispone*epsc2;
        }

        double P3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.8334754182339632 + -0.4704868585026707*as + -0.36573888111887276*epsb + -8.56611768163254*chi3spone*epsb + -2.79173073227914*chi3sppone*epsb + 0.36573888111887276*epsc + -8.56611768163254*chi3spone*epsc + -2.79173073227914*chi3sppone*epsc + 2.141529420408135*xispone + -0.2985161051164776*as*xispone + -0.9397284634534445*epsb*xispone + -5.58346146455828*chi3spone*epsb*xispone + 0.9397284634534445*epsc*xispone + -5.58346146455828*chi3spone*epsc*xispone + 0.697932683069785*xisppone + -0.04908980844491689*as*xisppone + -0.30626112422500035*epsb*xisppone + 0.30626112422500035*epsc*xisppone + 0.8334754182339632*l2sone*epsc2 + 1.7925630788732423*l2spone*epsc2 + -0.36573888111887276*l5sone*epsc2 + -0.7865979013409443*l5spone*epsc2 + 2.141529420408135*l2sone*xispone*epsc2 + 1.39586536613957*l2spone*xispone*epsc2 + -0.9397284634534445*l5sone*xispone*epsc2 + -0.6125222484500007*l5spone*xispone*epsc2 + 0.697932683069785*l2sone*xisppone*epsc2 + -0.30626112422500035*l5sone*xisppone*epsc2;
        }

        double V2s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.004609867550019866 + 0.0016868153974270927*as + -0.0020228644576260183*epsb + 0.0020228644576260183*epsc + (0.004609867550019866*l2sone - 0.0020228644576260183*l5sone)*epsc2;
        }

        double V2s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04453083033743383 + -0.01954065555696522*epsb + -0.1475157616006357*chi3spone*epsb + 0.01954065555696522*epsc + -0.1475157616006357*chi3spone*epsc + as*(0.020136736299265223 + 0.013494523179416741*xispone) + 0.036878940400158926*xispone + -0.01618291566100815*epsb*xispone + 0.01618291566100815*epsc*xispone + (0.036878940400158926*l2spone - 0.01954065555696522*l5sone - 0.01618291566100815*l5spone + l2sone*(0.04453083033743383 + 0.036878940400158926*xispone) - 0.01618291566100815*l5sone*xispone)*epsc2;
        }

        double V2s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.14479527883566448 + -0.06353788259869902*epsb + -1.720018093999154*chi3spone*epsb + -0.5900630464025428*chi3sppone*epsb + 0.06353788259869902*epsc + -1.720018093999154*chi3spone*epsc + -0.5900630464025428*chi3sppone*epsc + 0.4300045234997885*xispone + -0.18869107577773805*epsb*xispone + -1.1801260928050856*chi3spone*epsb*xispone + 0.18869107577773805*epsc*xispone + -1.1801260928050856*chi3spone*epsc*xispone + as*(0.04740742180967111 + 0.18808293675295523*xispone + 0.053978092717666966*xisppone) + 0.1475157616006357*xisppone + -0.06473166264403259*epsb*xisppone + 0.06473166264403259*epsc*xisppone + (0.3562466426994706*l2spone - 0.06353788259869902*l5sone - 0.15632524445572177*l5spone + 0.2950315232012714*l2spone*xispone - 0.18869107577773805*l5sone*xispone - 0.12946332528806517*l5spone*xispone + l2sone*(0.14479527883566448 + 0.4300045234997885*xispone + 0.1475157616006357*xisppone) - 0.06473166264403259*l5sone*xisppone)*epsc2;
        }

        double V3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0032596686049908637 + 0.0011927586061305786*as + -0.0014303811754086053*epsb + 0.0014303811754086053*epsc + -0.0028607623508172106*epsb*etasone + 0.0028607623508172106*epsc*etasone + (0.0032596686049908637*l2sone + 0.0014303811754086053*l5sone - 0.0028607623508172106*l6sone)*epsc2;
        }

        double V3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.03148805210346711 + -0.0138173300531607*epsb + 0.10430939535970762*chi2sone*epsb + -0.10430939535970764*chi3spone*epsb + 0.013817330053160699*epsc + 0.10430939535970761*chi2sone*epsc + -0.10430939535970764*chi3spone*epsc + -0.0276346601063214*epsb*etasone + 0.0276346601063214*epsc*etasone + -0.022886098806537684*epsb*etaspone + 0.022886098806537684*epsc*etaspone + as*(0.014238822788175745 + 0.009542068849044629*xispone) + 0.02607734883992691*xispone + -0.011443049403268842*epsb*xispone + 0.011443049403268842*epsc*xispone + -0.022886098806537684*epsb*etasone*xispone + 0.022886098806537684*epsc*etasone*xispone + (0.02607734883992691*l2spone + 0.026077348839926903*l3sone + 0.0138173300531607*l5sone + 0.011443049403268842*l5spone - 0.03907770950959024*l6sone - 0.022886098806537684*l6spone + l2sone*(0.03148805210346711 + 0.02607734883992691*xispone) + 0.011443049403268842*l5sone*xispone - 0.022886098806537684*l6sone*xispone)*epsc2;
        }

        double V3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.10238572354849537 + -0.044928067647774815*epsb + 1.2162364580303622*chi2sone*epsb + 0.8344751628776611*chi2spone*epsb + -1.2162364580303626*chi3spone*epsb + -0.41723758143883055*chi3sppone*epsb + 0.04492806764777479*epsc + 1.2162364580303622*chi2sone*epsc + 0.834475162877661*chi2spone*epsc + -1.2162364580303626*chi3spone*epsc + -0.41723758143883055*chi3sppone*epsc + -0.08985613529554964*epsb*etasone + 0.08985613529554964*epsc*etasone + -0.2668494784636466*epsb*etaspone + 0.2668494784636466*epsc*etaspone + -0.09154439522615074*epsb*etasppone + 0.09154439522615074*epsc*etasppone + 0.30405911450759066*xispone + -0.1334247392318233*epsb*xispone + 0.8344751628776611*chi2sone*epsb*xispone + -0.8344751628776611*chi3spone*epsb*xispone + 0.13342473923182327*epsc*xispone + 0.834475162877661*chi2sone*epsc*xispone + -0.8344751628776611*chi3spone*epsc*xispone + -0.2668494784636466*epsb*etasone*xispone + 0.2668494784636466*epsc*etasone*xispone + -0.18308879045230148*epsb*etaspone*xispone + 0.18308879045230148*epsc*etaspone*xispone + as*(0.03352210944018954 + 0.13299472000349521*xispone + 0.03816827539617852*xisppone) + 0.10430939535970764*xisppone + -0.04577219761307537*epsb*xisppone + 0.04577219761307537*epsc*xisppone + -0.09154439522615074*epsb*etasone*xisppone + 0.09154439522615074*epsc*etasone*xisppone + (0.2519044168277369*l2spone + 0.30405911450759054*l3sone + 0.20861879071941525*l3spone + 0.04492806764777482*l5sone + 0.11053864042528559*l5spone - 0.22328087452737289*l6sone - 0.3126216760767219*l6spone + 0.20861879071941528*l2spone*xispone + 0.20861879071941525*l3sone*xispone + 0.1334247392318233*l5sone*xispone + 0.09154439522615074*l5spone*xispone - 0.35839387368979736*l6sone*xispone - 0.18308879045230148*l6spone*xispone + l2sone*(0.10238572354849537 + 0.30405911450759066*xispone + 0.10430939535970764*xisppone) + 0.04577219761307537*l5sone*xisppone - 0.09154439522615074*l6sone*xisppone)*epsc2;
        }

        double V6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.00414233004658001 + 0.003800059460515357*as + 0.00414233004658001*epsb + 0.00414233004658001*epsc + 0.00828466009316002*epsc*etasone + 0.00414233004658001*l2sone*epsc2 + 0.00414233004658001*l5sone*epsc2 + -0.00828466009316002*l6sone*epsc2;
        }

        double V6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04001446777036228 + 0.04321714176889576*as + 0.04001446777036228*epsb + -0.13255456149056033*chi3spone*epsb + 0.04001446777036228*epsc + 0.13255456149056033*chi2sone*epsc + -0.13255456149056033*chi3spone*epsc + 0.08002893554072456*epsc*etasone + 0.06627728074528016*epsc*etaspone + 0.03313864037264008*xispone + 0.030400475684122855*as*xispone + 0.03313864037264008*epsb*xispone + 0.03313864037264008*epsc*xispone + 0.06627728074528016*epsc*etasone*xispone + 0.04001446777036228*l2sone*epsc2 + 0.03313864037264008*l2spone*epsc2 + 0.03313864037264008*l3sone*epsc2 + 0.04001446777036228*l5sone*epsc2 + 0.03313864037264008*l5spone*epsc2 + -0.11316757591336464*l6sone*epsc2 + -0.06627728074528016*l6spone*epsc2 + 0.03313864037264008*l2sone*xispone*epsc2 + 0.03313864037264008*l5sone*xispone*epsc2 + -0.06627728074528016*l6sone*xispone*epsc2;
        }

        double V6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.1301099928828365 + 0.14509523330458124*as + 0.1301099928828365*epsb + -1.5455720916327136*chi3spone*epsb + -0.5302182459622413*chi3sppone*epsb + 0.1301099928828365*epsc + 1.5455720916327136*chi2sone*epsc + 1.0604364919244826*chi2spone*epsc + -1.5455720916327136*chi3spone*epsc + -0.5302182459622413*chi3sppone*epsc + 0.260219985765673*epsc*etasone + 0.7727860458163568*epsc*etaspone + 0.26510912298112066*epsc*etasppone + 0.3863930229081784*xispone + 0.4065380855194118*as*xispone + 0.3863930229081784*epsb*xispone + -1.0604364919244826*chi3spone*epsb*xispone + 0.3863930229081784*epsc*xispone + 1.0604364919244826*chi2sone*epsc*xispone + -1.0604364919244826*chi3spone*epsc*xispone + 0.7727860458163568*epsc*etasone*xispone + 0.5302182459622413*epsc*etaspone*xispone + 0.13255456149056033*xisppone + 0.12160190273649144*as*xisppone + 0.13255456149056033*epsb*xisppone + 0.13255456149056033*epsc*xisppone + 0.26510912298112066*epsc*etasone*xisppone + 0.1301099928828365*l2sone*epsc2 + 0.32011574216289823*l2spone*epsc2 + 0.3863930229081784*l3sone*epsc2 + 0.26510912298112066*l3spone*epsc2 + 0.1301099928828365*l5sone*epsc2 + 0.32011574216289823*l5spone*epsc2 + -0.6466130086738514*l6sone*epsc2 + -0.9053406073069172*l6spone*epsc2 + 0.3863930229081784*l2sone*xispone*epsc2 + 0.26510912298112066*l2spone*xispone*epsc2 + 0.26510912298112066*l3sone*xispone*epsc2 + 0.3863930229081784*l5sone*xispone*epsc2 + 0.26510912298112066*l5spone*xispone*epsc2 + -1.0378951687974776*l6sone*xispone*epsc2 + -0.5302182459622413*l6spone*xispone*epsc2 + 0.13255456149056033*l2sone*xisppone*epsc2 + 0.13255456149056033*l5sone*xisppone*epsc2 + -0.26510912298112066*l6sone*xisppone*epsc2;
        }

        double V7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.00414233004658001 + 0.003800059460515357*as + 0.00414233004658001*epsb + 0.00414233004658001*epsc + 0.00828466009316002*epsb*etasone + 0.00414233004658001*l2sone*epsc2 + -0.00414233004658001*l5sone*epsc2;
        }

        double V7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04001446777036228 + 0.04321714176889576*as + 0.04001446777036228*epsb + 0.13255456149056033*chi2sone*epsb + -0.13255456149056033*chi3spone*epsb + 0.04001446777036228*epsc + -0.13255456149056033*chi3spone*epsc + 0.08002893554072456*epsb*etasone + 0.06627728074528016*epsb*etaspone + 0.03313864037264008*xispone + 0.030400475684122855*as*xispone + 0.03313864037264008*epsb*xispone + 0.03313864037264008*epsc*xispone + 0.06627728074528016*epsb*etasone*xispone + 0.04001446777036228*l2sone*epsc2 + 0.03313864037264008*l2spone*epsc2 + -0.04001446777036228*l5sone*epsc2 + -0.03313864037264008*l5spone*epsc2 + 0.03313864037264008*l2sone*xispone*epsc2 + -0.03313864037264008*l5sone*xispone*epsc2;
        }

        double V7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.1301099928828365 + 0.14509523330458124*as + 0.1301099928828365*epsb + 1.5455720916327136*chi2sone*epsb + 1.0604364919244826*chi2spone*epsb + -1.5455720916327136*chi3spone*epsb + -0.5302182459622413*chi3sppone*epsb + 0.1301099928828365*epsc + -1.5455720916327136*chi3spone*epsc + -0.5302182459622413*chi3sppone*epsc + 0.260219985765673*epsb*etasone + 0.7727860458163568*epsb*etaspone + 0.26510912298112066*epsb*etasppone + 0.3863930229081784*xispone + 0.4065380855194118*as*xispone + 0.3863930229081784*epsb*xispone + 1.0604364919244826*chi2sone*epsb*xispone + -1.0604364919244826*chi3spone*epsb*xispone + 0.3863930229081784*epsc*xispone + -1.0604364919244826*chi3spone*epsc*xispone + 0.7727860458163568*epsb*etasone*xispone + 0.5302182459622413*epsb*etaspone*xispone + 0.13255456149056033*xisppone + 0.12160190273649144*as*xisppone + 0.13255456149056033*epsb*xisppone + 0.13255456149056033*epsc*xisppone + 0.26510912298112066*epsb*etasone*xisppone + 0.1301099928828365*l2sone*epsc2 + 0.32011574216289823*l2spone*epsc2 + -0.1301099928828365*l5sone*epsc2 + -0.32011574216289823*l5spone*epsc2 + 0.3863930229081784*l2sone*xispone*epsc2 + 0.26510912298112066*l2spone*xispone*epsc2 + -0.3863930229081784*l5sone*xispone*epsc2 + -0.26510912298112066*l5spone*xispone*epsc2 + 0.13255456149056033*l2sone*xisppone*epsc2 + -0.13255456149056033*l5sone*xisppone*epsc2;
        }

        double A3s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0032646160303035363 + -0.0013579525700927507*as + 0.0032646160303035363*l2sone*epsc2;
        }

        double A3s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.03958930938988248 + 0.0002697117831198897*as + 0.013058464121214145*epsb + 0.10446771296971316*chi2sone*epsb + -0.10446771296971316*chi3spone*epsb + 0.013058464121214145*epsc + -0.10446771296971316*chi3spone*epsc + 0.02611692824242829*epsb*etasone + 0.02611692824242829*xispone + -0.010863620560742008*as*xispone + 0.03958930938988248*l2sone*epsc2 + 0.02611692824242829*l2spone*epsc2 + -0.013058464121214145*l5sone*epsc2 + 0.02611692824242829*l2sone*xispone*epsc2;
        }

        double A3s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.19459138473995258 + 0.07882797120349648*as + 0.1322403093171016*epsb + 1.4757933264156655*chi2sone*epsb + 0.8357417037577054*chi2spone*epsb + -1.4757933264156655*chi3spone*epsb + -0.4178708518788527*chi3sppone*epsb + 0.1322403093171016*epsc + -1.4757933264156655*chi3spone*epsc + -0.4178708518788527*chi3sppone*epsc + 0.2644806186342032*epsb*etasone + 0.20893542593942635*epsb*etaspone + 0.3689483316039164*xispone + -0.019569546856524914*as*xispone + 0.10446771296971318*epsb*xispone + 0.8357417037577054*chi2sone*epsb*xispone + -0.8357417037577054*chi3spone*epsb*xispone + 0.10446771296971318*epsc*xispone + -0.8357417037577054*chi3spone*epsc*xispone + 0.20893542593942635*epsb*etasone*xispone + 0.10446771296971318*xisppone + -0.04345448224296803*as*xisppone + 0.19459138473995258*l2sone*epsc2 + 0.3167144751190598*l2spone*epsc2 + -0.1322403093171016*l5sone*epsc2 + -0.10446771296971318*l5spone*epsc2 + 0.3689483316039164*l2sone*xispone*epsc2 + 0.20893542593942635*l2spone*xispone*epsc2 + -0.10446771296971318*l5sone*xispone*epsc2 + 0.10446771296971318*l2sone*xisppone*epsc2;
        }

        double A4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0032646160303035363 + -0.0013579525700927507*as + 0.0032646160303035363*l2sone*epsc2;
        }

        double A4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.03958930938988248 + 0.0002697117831198897*as + 0.013058464121214145*epsb + -0.10446771296971316*chi3spone*epsb + 0.013058464121214145*epsc + 0.10446771296971316*chi2sone*epsc + -0.10446771296971316*chi3spone*epsc + 0.02611692824242829*epsc*etasone + 0.02611692824242829*xispone + -0.010863620560742008*as*xispone + 0.03958930938988248*l2sone*epsc2 + 0.02611692824242829*l2spone*epsc2 + 0.02611692824242829*l3sone*epsc2 + 0.013058464121214145*l5sone*epsc2 + -0.02611692824242829*l6sone*epsc2 + 0.02611692824242829*l2sone*xispone*epsc2;
        }

        double A4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.19459138473995258 + 0.07882797120349648*as + 0.1322403093171016*epsb + -1.4757933264156655*chi3spone*epsb + -0.4178708518788527*chi3sppone*epsb + 0.1322403093171016*epsc + 1.4757933264156655*chi2sone*epsc + 0.8357417037577054*chi2spone*epsc + -1.4757933264156655*chi3spone*epsc + -0.4178708518788527*chi3sppone*epsc + 0.2644806186342032*epsc*etasone + 0.20893542593942635*epsc*etaspone + 0.3689483316039164*xispone + -0.019569546856524914*as*xispone + 0.10446771296971318*epsb*xispone + -0.8357417037577054*chi3spone*epsb*xispone + 0.10446771296971318*epsc*xispone + 0.8357417037577054*chi2sone*epsc*xispone + -0.8357417037577054*chi3spone*epsc*xispone + 0.20893542593942635*epsc*etasone*xispone + 0.10446771296971318*xisppone + -0.04345448224296803*as*xisppone + 0.19459138473995258*l2sone*epsc2 + 0.3167144751190598*l2spone*epsc2 + 0.3689483316039164*l3sone*epsc2 + 0.20893542593942635*l3spone*epsc2 + 0.1322403093171016*l5sone*epsc2 + 0.10446771296971318*l5spone*epsc2 + -0.3689483316039164*l6sone*epsc2 + -0.20893542593942635*l6spone*epsc2 + 0.3689483316039164*l2sone*xispone*epsc2 + 0.20893542593942635*l2spone*xispone*epsc2 + 0.20893542593942635*l3sone*xispone*epsc2 + 0.10446771296971318*l5sone*xispone*epsc2 + -0.20893542593942635*l6sone*xispone*epsc2 + 0.10446771296971318*l2sone*xisppone*epsc2;
        }

        double A7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0030181566001492937 + -0.0012554350876400194*as + 0.0030181566001492937*l2sone*epsc2;
        }

        double A7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.03072546294993736 + 0.0048262778721295235*as + -0.02751207995069488*epsb + -0.0965810112047774*chi3spone*epsb + 0.02751207995069488*epsc + -0.0965810112047774*chi3spone*epsc + 0.02414525280119435*xispone + -0.010043480701120155*as*xispone + 0.03072546294993736*l2sone*epsc2 + 0.02414525280119435*l2spone*epsc2 + -0.02751207995069488*l5sone*epsc2 + 0.02414525280119435*l2sone*xispone*epsc2;
        }

        double A7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.11437742691975071 + 0.08810044101490708*as + -0.2250545454775743*epsb + -1.1763768368075502*chi3spone*epsb + -0.3863240448191096*chi3sppone*epsb + 0.2250545454775743*epsc + -1.1763768368075502*chi3spone*epsc + -0.3863240448191096*chi3sppone*epsc + 0.29409420920188756*xispone + 0.01852326157479586*as*xispone + -0.22009663960555897*epsb*xispone + -0.7726480896382192*chi3spone*epsb*xispone + 0.22009663960555897*epsc*xispone + -0.7726480896382192*chi3spone*epsc*xispone + 0.0965810112047774*xisppone + -0.040173922804480615*as*xisppone + 0.11437742691975071*l2sone*epsc2 + 0.24580370359949882*l2spone*epsc2 + -0.2250545454775743*l5sone*epsc2 + -0.22009663960555897*l5spone*epsc2 + 0.29409420920188756*l2sone*xispone*epsc2 + 0.1931620224095548*l2spone*xispone*epsc2 + -0.22009663960555897*l5sone*xispone*epsc2 + 0.0965810112047774*l2sone*xisppone*epsc2;
        }

        double T4s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0038692335355298545 + 0.0012933841461525606*as + 0.0038692335355298545*l2sone*epsc2;
        }

        double T4s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.046921378231901255 + 0.03847593352404564*as + 0.015476934142119418*epsb + -0.12381547313695535*chi3spone*epsb + 0.015476934142119418*epsc + -0.12381547313695535*chi3spone*epsc + 0.030953868284238836*xispone + 0.010347073169220485*as*xispone + 0.046921378231901255*l2sone*epsc2 + 0.030953868284238836*l2spone*epsc2 + -0.015476934142119418*l5sone*epsc2 + 0.030953868284238836*l2sone*xispone*epsc2;
        }

        double T4s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.23063034199798754 + 0.30375562024325803*as + 0.1567316446433662*epsb + -1.749115049694751*chi3spone*epsb + -0.4952618925478214*chi3sppone*epsb + 0.1567316446433662*epsc + -1.749115049694751*chi3spone*epsc + -0.4952618925478214*chi3sppone*epsc + 0.43727876242368774*xispone + 0.32850161453080606*as*xispone + 0.12381547313695535*epsb*xispone + -0.9905237850956428*chi3spone*epsb*xispone + 0.12381547313695535*epsc*xispone + -0.9905237850956428*chi3spone*epsc*xispone + 0.12381547313695535*xisppone + 0.04138829267688194*as*xisppone + 0.23063034199798754*l2sone*epsc2 + 0.37537102585521004*l2spone*epsc2 + -0.1567316446433662*l5sone*epsc2 + -0.12381547313695535*l5spone*epsc2 + 0.43727876242368774*l2sone*xispone*epsc2 + 0.2476309462739107*l2spone*xispone*epsc2 + -0.12381547313695535*l5sone*xispone*epsc2 + 0.12381547313695535*l2sone*xisppone*epsc2;
        }

        double T5s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0008942822543584558 + 0.00029893530058386786*as + 0.0008942822543584558*l2sone*epsc2;
        }

        double T5s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.00999826190595043 + 0.008120804260731354*as + -0.0081518516564649*epsb + -0.028617032139470585*chi3spone*epsb + 0.008151851656464896*epsc + 0.028617032139470585*chi2sone*epsc + -0.028617032139470585*chi3spone*epsc + 0.0163037033129298*epsc*etasone + 0.007154258034867646*xispone + 0.002391482404670942*as*xispone + 0.00999826190595043*l2sone*epsc2 + 0.007154258034867646*l2spone*epsc2 + 0.007154258034867646*l3sone*epsc2 + 0.0081518516564649*l5sone*epsc2 + -0.0163037033129298*l6sone*epsc2 + 0.007154258034867646*l2sone*xispone*epsc2;
        }

        double T5s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04299410423573535 + 0.054628550455241724*as + -0.07483569644879517*epsb + -0.37717844526935496*chi3spone*epsb + -0.11446812855788234*chi3sppone*epsb + 0.07483569644879515*epsc + 0.37717844526935496*chi2sone*epsc + 0.22893625711576468*chi2spone*epsc + -0.37717844526935496*chi3spone*epsc + -0.11446812855788234*chi3sppone*epsc + 0.14967139289759035*epsc*etasone + 0.1304296265034384*epsc*etaspone + 0.09429461131733874*xispone + 0.0697493988951927*as*xispone + -0.06521481325171918*epsb*xispone + -0.22893625711576468*chi3spone*epsb*xispone + 0.0652148132517192*epsc*xispone + 0.22893625711576468*chi2sone*epsc*xispone + -0.22893625711576468*chi3spone*epsc*xispone + 0.1304296265034384*epsc*etasone*xispone + 0.028617032139470585*xisppone + 0.009565929618683768*as*xisppone + 0.04299410423573535*l2sone*epsc2 + 0.07998609524760344*l2spone*epsc2 + 0.09429461131733874*l3sone*epsc2 + 0.05723406427894117*l3spone*epsc2 + 0.07483569644879517*l5sone*epsc2 + 0.06521481325171918*l5spone*epsc2 + -0.21488620614930953*l6sone*epsc2 + -0.1304296265034384*l6spone*epsc2 + 0.09429461131733874*l2sone*xispone*epsc2 + 0.05723406427894117*l2spone*xispone*epsc2 + 0.05723406427894117*l3sone*xispone*epsc2 + 0.06521481325171918*l5sone*xispone*epsc2 + -0.1304296265034384*l6sone*xispone*epsc2 + 0.028617032139470585*l2sone*xisppone*epsc2;
        }

        double T6s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.0008942822543584558 + 0.00029893530058386786*as + 0.0008942822543584558*l2sone*epsc2;
        }

        double T6s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.009998261905950432 + 0.008120804260731355*as + -0.0081518516564649*epsb + 0.028617032139470585*chi2sone*epsb + -0.028617032139470585*chi3spone*epsb + 0.008151851656464901*epsc + -0.028617032139470585*chi3spone*epsc + -0.016303703312929792*epsb*etasone + 0.007154258034867646*xispone + 0.002391482404670943*as*xispone + 0.009998261905950432*l2sone*epsc2 + 0.007154258034867646*l2spone*epsc2 + -0.008151851656464901*l5sone*epsc2 + 0.007154258034867646*l2sone*xispone*epsc2;
        }

        double T6s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04299410423573535 + 0.05462855045524174*as + -0.07483569644879517*epsb + 0.37717844526935496*chi2sone*epsb + 0.22893625711576468*chi2spone*epsb + -0.37717844526935496*chi3spone*epsb + -0.11446812855788234*chi3sppone*epsb + 0.07483569644879517*epsc + -0.37717844526935496*chi3spone*epsc + -0.11446812855788234*chi3sppone*epsc + -0.1496713928975903*epsb*etasone + -0.13042962650343837*epsb*etaspone + 0.09429461131733874*xispone + 0.06974939889519273*as*xispone + -0.0652148132517192*epsb*xispone + 0.22893625711576468*chi2sone*epsb*xispone + -0.22893625711576468*chi3spone*epsb*xispone + 0.0652148132517192*epsc*xispone + -0.22893625711576468*chi3spone*epsc*xispone + -0.13042962650343837*epsb*etasone*xispone + 0.028617032139470585*xisppone + 0.009565929618683768*as*xisppone + 0.04299410423573535*l2sone*epsc2 + 0.07998609524760344*l2spone*epsc2 + -0.07483569644879517*l5sone*epsc2 + -0.0652148132517192*l5spone*epsc2 + 0.09429461131733874*l2sone*xispone*epsc2 + 0.05723406427894117*l2spone*xispone*epsc2 + -0.0652148132517192*l5sone*xispone*epsc2 + 0.028617032139470585*l2sone*xisppone*epsc2;
        }

        double T7s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.008880721929533939 + 0.014809556958285919*as + 0.008880721929533939*epsb + 0.008880721929533939*epsc + 0.008880721929533939*l2sone*epsc2 + -0.008880721929533939*l5sone*epsc2;
        }

        double T7s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.09419314822522039 + 0.17781229715364028*as + 0.09419314822522039*epsb + -0.28418310174508604*chi3spone*epsb + 0.09419314822522039*epsc + -0.28418310174508604*chi3spone*epsc + 0.07104577543627151*xispone + 0.11847645566628734*as*xispone + 0.07104577543627151*epsb*xispone + 0.07104577543627151*epsc*xispone + 0.09419314822522039*l2sone*epsc2 + 0.07104577543627151*l2spone*epsc2 + -0.09419314822522039*l5sone*epsc2 + -0.07104577543627151*l5spone*epsc2 + 0.07104577543627151*l2sone*xispone*epsc2 + -0.07104577543627151*l5sone*xispone*epsc2;
        }

        double T7s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.3685526663349023 + 0.7582917620602395*as + 0.3685526663349023*epsb + -3.5825469466972244*chi3spone*epsb + -1.1367324069803442*chi3sppone*epsb + 0.3685526663349023*epsc + -3.5825469466972244*chi3spone*epsc + -1.1367324069803442*chi3sppone*epsc + 0.8956367366743061*xispone + 1.6594512885616968*as*xispone + 0.8956367366743061*epsb*xispone + -2.2734648139606883*chi3spone*epsb*xispone + 0.8956367366743061*epsc*xispone + -2.2734648139606883*chi3spone*epsc*xispone + 0.28418310174508604*xisppone + 0.47390582266514936*as*xisppone + 0.28418310174508604*epsb*xisppone + 0.28418310174508604*epsc*xisppone + 0.3685526663349023*l2sone*epsc2 + 0.7535451858017631*l2spone*epsc2 + -0.3685526663349023*l5sone*epsc2 + -0.7535451858017631*l5spone*epsc2 + 0.8956367366743061*l2sone*xispone*epsc2 + 0.5683662034901721*l2spone*xispone*epsc2 + -0.8956367366743061*l5sone*xispone*epsc2 + -0.5683662034901721*l5spone*xispone*epsc2 + 0.28418310174508604*l2sone*xisppone*epsc2 + -0.28418310174508604*l5sone*xisppone*epsc2;
        }

        double T8s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.004677562565499274 + 0.002641374892133541*as + -0.0020525698318663577*epsb + 0.0020525698318663577*epsc + -0.004105139663732715*epsb*etasone + 0.004677562565499274*l2sone*epsc2 + -0.0020525698318663577*l5sone*epsc2;
        }

        double T8s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04518475698853479 + 0.029053820100103885*as + -0.01982760631337054*epsb + 0.14968200209597676*chi2sone*epsb + -0.14968200209597676*chi3spone*epsb + 0.019827606313370535*epsc + -0.14968200209597676*chi3spone*epsc + -0.03965521262674108*epsb*etasone + -0.03284111730986172*epsb*etaspone + 0.03742050052399419*xispone + 0.021130999137068326*as*xispone + -0.01642055865493086*epsb*xispone + 0.01642055865493086*epsc*xispone + -0.03284111730986172*epsb*etasone*xispone + 0.04518475698853479*l2sone*epsc2 + 0.03742050052399419*l2spone*epsc2 + -0.019827606313370535*l5sone*epsc2 + -0.01642055865493086*l5spone*epsc2 + 0.03742050052399419*l2sone*xispone*epsc2 + -0.01642055865493086*l5sone*xispone*epsc2;
        }

        double T8s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.1469215695665301 + 0.07347204786986515*as + -0.06447092414476885*epsb + 1.7452762278250669*chi2sone*epsb + 1.197456016767814*chi2spone*epsb + -1.7452762278250669*chi3spone*epsb + -0.598728008383907*chi3sppone*epsb + 0.06447092414476885*epsc + -1.7452762278250669*chi3spone*epsc + -0.598728008383907*chi3sppone*epsc + -0.12894184828953772*epsb*etasone + -0.38292393563365207*epsb*etaspone + -0.13136446923944692*epsb*etasppone + 0.4363190569562667*xispone + 0.27469255907496776*as*xispone + -0.19146196781682606*epsb*xispone + 1.197456016767814*chi2sone*epsb*xispone + -1.197456016767814*chi3spone*epsb*xispone + 0.19146196781682603*epsc*xispone + -1.197456016767814*chi3spone*epsc*xispone + -0.38292393563365207*epsb*etasone*xispone + -0.26272893847889384*epsb*etaspone*xispone + 0.14968200209597676*xisppone + 0.08452399654827332*as*xisppone + -0.06568223461972346*epsb*xisppone + 0.06568223461972346*epsc*xisppone + -0.13136446923944692*epsb*etasone*xisppone + 0.1469215695665301*l2sone*epsc2 + 0.3614780559082783*l2spone*epsc2 + -0.06447092414476885*l5sone*epsc2 + -0.15862085050696428*l5spone*epsc2 + 0.4363190569562667*l2sone*xispone*epsc2 + 0.2993640041919535*l2spone*xispone*epsc2 + -0.19146196781682603*l5sone*xispone*epsc2 + -0.13136446923944692*l5spone*xispone*epsc2 + 0.14968200209597676*l2sone*xisppone*epsc2 + -0.06568223461972346*l5sone*xisppone*epsc2;
        }

        double T9s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.004677562565499274 + 0.002641374892133541*as + -0.0020525698318663577*epsb + 0.0020525698318663577*epsc + 0.004105139663732715*epsc*etasone + 0.004677562565499274*l2sone*epsc2 + 0.0020525698318663577*l5sone*epsc2 + -0.004105139663732715*l6sone*epsc2;
        }

        double T9s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.04518475698853479 + 0.029053820100103885*as + -0.01982760631337054*epsb + -0.14968200209597676*chi3spone*epsb + 0.019827606313370542*epsc + 0.14968200209597676*chi2sone*epsc + -0.14968200209597676*chi3spone*epsc + 0.03965521262674108*epsc*etasone + 0.03284111730986172*epsc*etaspone + 0.03742050052399419*xispone + 0.021130999137068326*as*xispone + -0.01642055865493086*epsb*xispone + 0.01642055865493086*epsc*xispone + 0.03284111730986172*epsc*etasone*xispone + 0.04518475698853479*l2sone*epsc2 + 0.03742050052399419*l2spone*epsc2 + 0.03742050052399419*l3sone*epsc2 + 0.01982760631337054*l5sone*epsc2 + 0.01642055865493086*l5spone*epsc2 + -0.05607577128167193*l6sone*epsc2 + -0.03284111730986172*l6spone*epsc2 + 0.03742050052399419*l2sone*xispone*epsc2 + 0.01642055865493086*l5sone*xispone*epsc2 + -0.03284111730986172*l6sone*xispone*epsc2;
        }

        double T9s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.1469215695665301 + 0.07347204786986515*as + -0.06447092414476885*epsb + -1.7452762278250669*chi3spone*epsb + -0.598728008383907*chi3sppone*epsb + 0.06447092414476886*epsc + 1.7452762278250669*chi2sone*epsc + 1.197456016767814*chi2spone*epsc + -1.7452762278250669*chi3spone*epsc + -0.598728008383907*chi3sppone*epsc + 0.1289418482895377*epsc*etasone + 0.38292393563365207*epsc*etaspone + 0.13136446923944692*epsc*etasppone + 0.4363190569562667*xispone + 0.27469255907496776*as*xispone + -0.19146196781682606*epsb*xispone + -1.197456016767814*chi3spone*epsb*xispone + 0.19146196781682603*epsc*xispone + 1.197456016767814*chi2sone*epsc*xispone + -1.197456016767814*chi3spone*epsc*xispone + 0.38292393563365207*epsc*etasone*xispone + 0.26272893847889384*epsc*etaspone*xispone + 0.14968200209597676*xisppone + 0.08452399654827329*as*xisppone + -0.06568223461972346*epsb*xisppone + 0.06568223461972346*epsc*xisppone + 0.13136446923944692*epsc*etasone*xisppone + 0.1469215695665301*l2sone*epsc2 + 0.3614780559082783*l2spone*epsc2 + 0.4363190569562667*l3sone*epsc2 + 0.2993640041919535*l3spone*epsc2 + 0.06447092414476885*l5sone*epsc2 + 0.1586208505069643*l5spone*epsc2 + -0.32040381610636376*l6sone*epsc2 + -0.4486061702533755*l6spone*epsc2 + 0.4363190569562667*l2sone*xispone*epsc2 + 0.2993640041919535*l2spone*xispone*epsc2 + 0.2993640041919535*l3sone*xispone*epsc2 + 0.19146196781682606*l5sone*xispone*epsc2 + 0.13136446923944692*l5spone*xispone*epsc2 + -0.514288404873099*l6sone*xispone*epsc2 + -0.26272893847889384*l6spone*xispone*epsc2 + 0.14968200209597676*l2sone*xisppone*epsc2 + 0.06568223461972346*l5sone*xisppone*epsc2 + -0.13136446923944692*l6sone*xisppone*epsc2;
        }

        double T10s_a0() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.006279618698205529 + 0.010471938151572392*as + 0.006279618698205529*epsb + 0.006279618698205529*epsc + 0.012559237396411058*epsb*etasone + 0.012559237396411058*epsc*etasone + 0.006279618698205529*l2sone*epsc2 + 0.006279618698205529*l5sone*epsc2 + -0.012559237396411058*l6sone*epsc2;
        }

        double T10s_a1() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.06660461385136295 + 0.12573228109569648*as + 0.06660461385136295*epsb + 0.20094779834257692*chi2sone*epsb + -0.20094779834257692*chi3spone*epsb + 0.06660461385136295*epsc + 0.20094779834257692*chi2sone*epsc + -0.20094779834257692*chi3spone*epsc + 0.1332092277027259*epsb*etasone + 0.1332092277027259*epsc*etasone + 0.10047389917128846*epsb*etaspone + 0.10047389917128846*epsc*etaspone + 0.05023694958564423*xispone + 0.08377550521257913*as*xispone + 0.05023694958564423*epsb*xispone + 0.05023694958564423*epsc*xispone + 0.10047389917128846*epsb*etasone*xispone + 0.10047389917128846*epsc*etasone*xispone + 0.06660461385136295*l2sone*epsc2 + 0.05023694958564423*l2spone*epsc2 + 0.05023694958564423*l3sone*epsc2 + 0.06660461385136295*l5sone*epsc2 + 0.05023694958564423*l5spone*epsc2 + -0.18344617728837015*l6sone*epsc2 + -0.10047389917128846*l6spone*epsc2 + 0.05023694958564423*l2sone*xispone*epsc2 + 0.05023694958564423*l5sone*xispone*epsc2 + -0.10047389917128846*l6sone*xispone*epsc2;
        }

        double T10s_a2() const
        {
            const double as   = _alpha_s() / M_PI;
            const double epsb = _eps_b_s();
            const double epsc = _eps_c_s();
            const double epsc2 = power_of<2>(epsc);

            return 0.26060608958979237 + 0.5361932470706913*as + 0.26060608958979237*epsb + 2.533243239928768*chi2sone*epsb + 1.6075823867406152*chi2spone*epsb + -2.533243239928768*chi3spone*epsb + -0.8037911933703076*chi3sppone*epsb + 0.26060608958979237*epsc + 2.533243239928768*chi2sone*epsc + 1.6075823867406152*chi2spone*epsc + -2.533243239928768*chi3spone*epsc + -0.8037911933703076*chi3sppone*epsc + 0.5212121791795847*epsb*etasone + 0.5212121791795847*epsc*etasone + 1.266621619964384*epsb*etaspone + 1.266621619964384*epsc*etaspone + 0.4018955966851538*epsb*etasppone + 0.4018955966851538*epsc*etasppone + 0.633310809982192*xispone + 1.17340925919073*as*xispone + 0.633310809982192*epsb*xispone + 1.6075823867406152*chi2sone*epsb*xispone + -1.6075823867406152*chi3spone*epsb*xispone + 0.633310809982192*epsc*xispone + 1.6075823867406152*chi2sone*epsc*xispone + -1.6075823867406152*chi3spone*epsc*xispone + 1.266621619964384*epsb*etasone*xispone + 1.266621619964384*epsc*etasone*xispone + 0.8037911933703076*epsb*etaspone*xispone + 0.8037911933703076*epsc*etaspone*xispone + 0.2009477983425769*xisppone + 0.3351020208503165*as*xisppone + 0.2009477983425769*epsb*xisppone + 0.2009477983425769*epsc*xisppone + 0.4018955966851538*epsb*etasone*xisppone + 0.4018955966851538*epsc*etasone*xisppone + 0.26060608958979237*l2sone*epsc2 + 0.5328369108109036*l2spone*epsc2 + 0.633310809982192*l3sone*epsc2 + 0.4018955966851538*l3spone*epsc2 + 0.26060608958979237*l5sone*epsc2 + 0.5328369108109036*l5spone*epsc2 + -1.1545229891617768*l6sone*epsc2 + -1.467569418306961*l6spone*epsc2 + 0.633310809982192*l2sone*xispone*epsc2 + 0.4018955966851538*l2spone*xispone*epsc2 + 0.4018955966851538*l3sone*xispone*epsc2 + 0.633310809982192*l5sone*xispone*epsc2 + 0.4018955966851538*l5spone*xispone*epsc2 + -1.668517216649538*l6sone*xispone*epsc2 + -0.8037911933703076*l6spone*xispone*epsc2 + 0.2009477983425769*l2sone*xisppone*epsc2 + 0.2009477983425769*l5sone*xisppone*epsc2 + -0.4018955966851538*l6sone*xisppone*epsc2;
        }
        // }}}
    };

    const std::vector<OptionSpecification>
    Implementation<BGLCoefficients>::options
    {
        { "SU3F-limit-sslp"_ok, { "0"s, "1"s }, "0"s }
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
    double BGLCoefficients::fT_a0() const { return _imp->fT_a0(); }
    double BGLCoefficients::fT_a1() const { return _imp->fT_a1(); }
    double BGLCoefficients::fT_a2() const { return _imp->fT_a2(); }
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
    double BGLCoefficients::T1_a0() const { return _imp->T1_a0(); }
    double BGLCoefficients::T1_a1() const { return _imp->T1_a1(); }
    double BGLCoefficients::T1_a2() const { return _imp->T1_a2(); }
    double BGLCoefficients::T2_a0() const { return _imp->T2_a0(); }
    double BGLCoefficients::T2_a1() const { return _imp->T2_a1(); }
    double BGLCoefficients::T2_a2() const { return _imp->T2_a2(); }
    double BGLCoefficients::T23_a0() const { return _imp->T23_a0(); }
    double BGLCoefficients::T23_a1() const { return _imp->T23_a1(); }
    double BGLCoefficients::T23_a2() const { return _imp->T23_a2(); }
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
    double BGLCoefficients::T1bar_a0() const { return _imp->T1bar_a0(); }
    double BGLCoefficients::T1bar_a1() const { return _imp->T1bar_a1(); }
    double BGLCoefficients::T1bar_a2() const { return _imp->T1bar_a2(); }
    double BGLCoefficients::T2bar_a0() const { return _imp->T2bar_a0(); }
    double BGLCoefficients::T2bar_a1() const { return _imp->T2bar_a1(); }
    double BGLCoefficients::T2bar_a2() const { return _imp->T2bar_a2(); }
    double BGLCoefficients::T23bar_a0() const { return _imp->T23bar_a0(); }
    double BGLCoefficients::T23bar_a1() const { return _imp->T23bar_a1(); }
    double BGLCoefficients::T23bar_a2() const { return _imp->T23bar_a2(); }
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
    double BGLCoefficients::T4_a0() const { return _imp->T4_a0(); }
    double BGLCoefficients::T4_a1() const { return _imp->T4_a1(); }
    double BGLCoefficients::T4_a2() const { return _imp->T4_a2(); }
    double BGLCoefficients::T5_a0() const { return _imp->T5_a0(); }
    double BGLCoefficients::T5_a1() const { return _imp->T5_a1(); }
    double BGLCoefficients::T5_a2() const { return _imp->T5_a2(); }
    double BGLCoefficients::T6_a0() const { return _imp->T6_a0(); }
    double BGLCoefficients::T6_a1() const { return _imp->T6_a1(); }
    double BGLCoefficients::T6_a2() const { return _imp->T6_a2(); }
    double BGLCoefficients::T7_a0() const { return _imp->T7_a0(); }
    double BGLCoefficients::T7_a1() const { return _imp->T7_a1(); }
    double BGLCoefficients::T7_a2() const { return _imp->T7_a2(); }
    double BGLCoefficients::T8_a0() const { return _imp->T8_a0(); }
    double BGLCoefficients::T8_a1() const { return _imp->T8_a1(); }
    double BGLCoefficients::T8_a2() const { return _imp->T8_a2(); }
    double BGLCoefficients::T9_a0() const { return _imp->T9_a0(); }
    double BGLCoefficients::T9_a1() const { return _imp->T9_a1(); }
    double BGLCoefficients::T9_a2() const { return _imp->T9_a2(); }
    double BGLCoefficients::T10_a0() const { return _imp->T10_a0(); }
    double BGLCoefficients::T10_a1() const { return _imp->T10_a1(); }
    double BGLCoefficients::T10_a2() const { return _imp->T10_a2(); }
    // }}}

    // B_s -> D_s form factors
    // {{{
    double BGLCoefficients::V1s_a0() const { return _imp->V1s_a0(); }
    double BGLCoefficients::V1s_a1() const { return _imp->V1s_a1(); }
    double BGLCoefficients::V1s_a2() const { return _imp->V1s_a2(); }
    double BGLCoefficients::S1s_a0() const { return _imp->S1s_a0(); }
    double BGLCoefficients::S1s_a1() const { return _imp->S1s_a1(); }
    double BGLCoefficients::S1s_a2() const { return _imp->S1s_a2(); }
    double BGLCoefficients::fTs_a0() const { return _imp->fTs_a0(); }
    double BGLCoefficients::fTs_a1() const { return _imp->fTs_a1(); }
    double BGLCoefficients::fTs_a2() const { return _imp->fTs_a2(); }
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
    double BGLCoefficients::T1s_a0() const { return _imp->T1s_a0(); }
    double BGLCoefficients::T1s_a1() const { return _imp->T1s_a1(); }
    double BGLCoefficients::T1s_a2() const { return _imp->T1s_a2(); }
    double BGLCoefficients::T2s_a0() const { return _imp->T2s_a0(); }
    double BGLCoefficients::T2s_a1() const { return _imp->T2s_a1(); }
    double BGLCoefficients::T2s_a2() const { return _imp->T2s_a2(); }
    double BGLCoefficients::T23s_a0() const { return _imp->T23s_a0(); }
    double BGLCoefficients::T23s_a1() const { return _imp->T23s_a1(); }
    double BGLCoefficients::T23s_a2() const { return _imp->T23s_a2(); }
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
    double BGLCoefficients::T1bars_a0() const { return _imp->T1bars_a0(); }
    double BGLCoefficients::T1bars_a1() const { return _imp->T1bars_a1(); }
    double BGLCoefficients::T1bars_a2() const { return _imp->T1bars_a2(); }
    double BGLCoefficients::T2bars_a0() const { return _imp->T2bars_a0(); }
    double BGLCoefficients::T2bars_a1() const { return _imp->T2bars_a1(); }
    double BGLCoefficients::T2bars_a2() const { return _imp->T2bars_a2(); }
    double BGLCoefficients::T23bars_a0() const { return _imp->T23bars_a0(); }
    double BGLCoefficients::T23bars_a1() const { return _imp->T23bars_a1(); }
    double BGLCoefficients::T23bars_a2() const { return _imp->T23bars_a2(); }
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
    double BGLCoefficients::T4s_a0() const { return _imp->T4s_a0(); }
    double BGLCoefficients::T4s_a1() const { return _imp->T4s_a1(); }
    double BGLCoefficients::T4s_a2() const { return _imp->T4s_a2(); }
    double BGLCoefficients::T5s_a0() const { return _imp->T5s_a0(); }
    double BGLCoefficients::T5s_a1() const { return _imp->T5s_a1(); }
    double BGLCoefficients::T5s_a2() const { return _imp->T5s_a2(); }
    double BGLCoefficients::T6s_a0() const { return _imp->T6s_a0(); }
    double BGLCoefficients::T6s_a1() const { return _imp->T6s_a1(); }
    double BGLCoefficients::T6s_a2() const { return _imp->T6s_a2(); }
    double BGLCoefficients::T7s_a0() const { return _imp->T7s_a0(); }
    double BGLCoefficients::T7s_a1() const { return _imp->T7s_a1(); }
    double BGLCoefficients::T7s_a2() const { return _imp->T7s_a2(); }
    double BGLCoefficients::T8s_a0() const { return _imp->T8s_a0(); }
    double BGLCoefficients::T8s_a1() const { return _imp->T8s_a1(); }
    double BGLCoefficients::T8s_a2() const { return _imp->T8s_a2(); }
    double BGLCoefficients::T9s_a0() const { return _imp->T9s_a0(); }
    double BGLCoefficients::T9s_a1() const { return _imp->T9s_a1(); }
    double BGLCoefficients::T9s_a2() const { return _imp->T9s_a2(); }
    double BGLCoefficients::T10s_a0() const { return _imp->T10s_a0(); }
    double BGLCoefficients::T10s_a1() const { return _imp->T10s_a1(); }
    double BGLCoefficients::T10s_a2() const { return _imp->T10s_a2(); }
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
            opt_zorder_bound(o, "z-order-bound"_ok, { "1", "2" }, "2"),
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
                    result += power_of<2>(ff[i]) * nf; // to account for flavor symmetry
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
                    result += power_of<2>(ff[i]) * ns;
                }
            }

            return result;
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
                    result += power_of<2>(ff[i]) * nf; // to account for flavor symmetry
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
                    result += power_of<2>(ff[i]) * ns;
                }
            }

            return result;
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
                    result += power_of<2>(ff[i]) * nf; // to account for flavor symmetry
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
                    result += power_of<2>(ff[i]) * ns;
                }
            }

            return result;
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
                    result += power_of<2>(ff[i]) * nf; // to account for flavor symmetry
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
                    result += power_of<2>(ff[i]) * ns;
                }
            }

            return result;
        }

        double bound_1m_T() const
        {
            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=u,d
            const std::array<std::array<double, 3>, 7> bgl_coeffs_ud
            {{
                // B -> D f_T
                { bgl->fT_a0(),    bgl->fT_a1(),    bgl->fT_a2()    },
                // B -> D^* T_1
                { bgl->T1_a0(),    bgl->T1_a1(),    bgl->T1_a2()    },
                // B^* -> D Tbar_1
                { bgl->T1bar_a0(), bgl->T1bar_a1(), bgl->T1bar_a2() },
                // B^* -> D^* T_7
                { bgl->T7_a0(),    bgl->T7_a1(),    bgl->T7_a2()    },
                // B -> D^* T_8
                { bgl->T8_a0(),    bgl->T8_a1(),    bgl->T8_a2()    },
                // B -> D^* T_9
                { bgl->T9_a0(),    bgl->T9_a1(),    bgl->T9_a2()    },
                // B -> D^* T_10
                { bgl->T10_a0(),   bgl->T10_a1(),   bgl->T10_a2()   }
            }};

            double result = 0.0;
            for (const auto & ff : bgl_coeffs_ud)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += power_of<2>(ff[i]) * nf; // to account for flavor symmetry
                }
            }

            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=s
            const std::array<std::array<double, 3>, 7> bgl_coeffs_s
            {{
                // B_s -> D_s f_T
                { bgl->fTs_a0(),    bgl->fTs_a1(),    bgl->fTs_a2()    },
                // B_s -> D_s^* T_1
                { bgl->T1s_a0(),    bgl->T1s_a1(),    bgl->T1s_a2()    },
                // B_s^* -> D_s Tbar_1
                { bgl->T1bars_a0(), bgl->T1bars_a1(), bgl->T1bars_a2() },
                // B_s^* -> D_s^* T_7
                { bgl->T7s_a0(),    bgl->T7s_a1(),    bgl->T7s_a2()    },
                // B_s -> D_s^* T_8
                { bgl->T8s_a0(),    bgl->T8s_a1(),    bgl->T8s_a2()    },
                // B_s -> D_s^* T_9
                { bgl->T9s_a0(),    bgl->T9s_a1(),    bgl->T9s_a2()    },
                // B_s -> D_s^* T_10
                { bgl->T10s_a0(),   bgl->T10s_a1(),   bgl->T10s_a2()   }
            }};

            for (const auto & ff : bgl_coeffs_s)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += power_of<2>(ff[i]) * ns;
                }
            }

            return result;
        }

        double bound_1p_T() const
        {
            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=u,d
            const std::array<std::array<double, 3>, 7> bgl_coeffs_ud
            {{
                // B -> D^* T_2
                { bgl->T2_a0(), bgl->T2_a1(), bgl->T2_a2() },
                // B^* -> D Tbar_2
                { bgl->T2bar_a0(), bgl->T2bar_a1(), bgl->T2bar_a2() },
                // B -> D^* T_23
                { bgl->T23_a0(), bgl->T23_a1(), bgl->T23_a2() },
                // B^* -> D Tbar_23
                { bgl->T23bar_a0(), bgl->T23bar_a1(), bgl->T23bar_a2() },
                // B^* -> D^* T_4
                { bgl->T4_a0(), bgl->T4_a1(), bgl->T4_a2() },
                // B^* -> D^* T_5
                { bgl->T5_a0(), bgl->T5_a1(), bgl->T5_a2() },
                // B^* -> D^* T_6
                { bgl->T6_a0(), bgl->T6_a1(), bgl->T6_a2() }
            }};

            double result = 0.0;
            for (const auto & ff : bgl_coeffs_ud)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += power_of<2>(ff[i]) * nf; // to account for flavor symmetry
                }
            }


            // 7 rows of form factors with 3 columns (one column per z coefficient)
            // for spectator quark q=s
            const std::array<std::array<double, 3>, 7> bgl_coeffs_s
            {{
                // B_s -> D_s^* T_2
                { bgl->T2s_a0(), bgl->T2s_a1(), bgl->T2s_a2() },
                // B_s^* -> D_s Tbar_2
                { bgl->T2bars_a0(), bgl->T2bars_a1(), bgl->T2bars_a2() },
                // B_s -> D_s^* T_23
                { bgl->T23s_a0(), bgl->T23s_a1(), bgl->T23s_a2() },
                // B_s^* -> D_s Tbar_23
                { bgl->T23bars_a0(), bgl->T23bars_a1(), bgl->T23bars_a2() },
                // B_s^* -> D_s^* T_4
                { bgl->T4s_a0(), bgl->T4s_a1(), bgl->T4s_a2() },
                // B_s^* -> D_s^* T_5
                { bgl->T5s_a0(), bgl->T5s_a1(), bgl->T5s_a2() },
                // B_s^* -> D_s^* T_6
                { bgl->T6s_a0(), bgl->T6s_a1(), bgl->T6s_a2() }
            }};

            for (const auto & ff : bgl_coeffs_s)
            {
                for (unsigned i = 0 ; i <= zorder_bound ; ++i)
                {
                    result += power_of<2>(ff[i]) * ns;
                }
            }

            return result;
        }
        // }}}
    };

    const std::vector<OptionSpecification>
    Implementation<HQETUnitarityBounds>::options
    {
        { "SU3F-limit-sslp"_ok, { "0"s, "1"s }, "0"s },
        { "z-order-bound"_ok,   { "1"s, "2"s }, "2"s }
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
    HQETUnitarityBounds::bound_1p_T() const
    {
        return _imp->bound_1p_T();
    }

    double
    HQETUnitarityBounds::bound_1m_T() const
    {
        return _imp->bound_1m_T();
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
            cond_qq(-0.02/12), // (p["B->D^*::<qq>@BGL1997"], u),         // [BGL:1997] quark condensate
            cond_G2(0.02),     // (p["B->D^*::<alS/pi G^2>@BGL1997"], u), // [BGL:1997] gluon condensate
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
            // [BGL:1997] eq.(4.1) + (4.3) + (4.6, 4.9)
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
            opt_zorder_bound(o, "z-order-bound"_ok, { "1", "2" }, "2"),
            nf(p["B(*)->D(*)::n_f@BGL1997"], u)
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
                result += power_of<2>(_a_f_0[i]) * nf; // to account for flavor symmetry
            }

            return result;
        }

        double bound_0m() const
        {
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += power_of<2>(_a_F2[i]) * nf; // to account for flavor symmetry
            }

            return result;
        }

        double bound_1p() const
        {
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += power_of<2>(_a_f[i]) * nf; // to account for flavor symmetry
                result += power_of<2>(_a_F1[i]) * nf; // to account for flavor symmetry
            }

            return result;
        }

        double bound_1m() const
        {
            double result = 0.0;

            for (unsigned i = 0 ; i <= zorder_bound ; ++i)
            {
                result += power_of<2>(_a_f_p[i]) * nf; // to account for flavor symmetry
                result += power_of<2>(_a_g[i]) * nf; // to account for flavor symmetry
            }

            return result;
        }
        // }}}
    };

    const std::vector<OptionSpecification>
    Implementation<BGLUnitarityBounds>::options
    {
        { "z-order-bound"_ok, { "1"s, "2"s }, "2"s }
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
