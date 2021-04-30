/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Danny van Dyk
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

#ifndef MASTER_GUARD_EOS_FORM_FACTORS_MESONIC_HQET_HH
#define MASTER_GUARD_EOS_FORM_FACTORS_MESONIC_HQET_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/derivative.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <limits>

#include <iostream>

namespace eos
{
    using std::sqrt;

    /* HQET Form Factors, based on [BLPR2017] and [JS2018] */
    template <typename Process_, typename Transition_> class HQETFormFactors;

    class HQETFormFactorBase :
        public virtual ParameterUser
    {
        protected:
            std::shared_ptr<Model> _model;

            // spin avaraged mB mass
            UsedParameter _mBar;

            // parameter for modifying the z function
            UsedParameter _a;

            // option to determine the model for the leading-power IW function
            SwitchOption _opt_lp_model;
            std::function<double (const double &)> _xi;

            // option to determine if we use z^3 terms in the leading-power IW function
            SwitchOption _opt_lp_zorder;
            double _enable_lp_z3;
            double _enable_lp_z4;
            double _enable_lp_z5;

            // option to determine if we use z^2 terms in the subleading-power IW function
            SwitchOption _opt_slp_zorder;
            double _enable_slp_z2;

            // option to determine if we use z^2 terms in the subsubleading-power IW function
            SwitchOption _opt_sslp_zorder;
            double _enable_sslp_z1;
            double _enable_sslp_z2;

            // option to determine if we use the SU3_F-symmetry limit for the subsubleading-power IW functions
            SwitchOption _opt_sslp_limit;

            // parameters for the leading Isgur-Wise function xi
            UsedParameter _xipone, _xippone, _xipppone, _xippppone, _xipppppone;

            // parameters for the subleading Isgur-Wise function chi_2
            UsedParameter _chi2one, _chi2pone, _chi2ppone;

            // parameters for the subleading Isgur-Wise function chi_3
            UsedParameter _chi3pone, _chi3ppone;

            // parameters for the subleading Isgur-Wise function eta
            UsedParameter _etaone, _etapone, _etappone;

            // parameters for subsubleading 1/m_c corrections in h_+ (B->D), equal to delta_{h_+}
            UsedParameter _l1one, _l1pone, _l1ppone;

            // parameters for subsubleading 1/m_c corrections in h_A1 (B->D^*), equal to delta_{A_1}
            UsedParameter _l2one, _l2pone, _l2ppone;

            // parameters for subsubleading 1/m_c corrections
            UsedParameter _l3one, _l3pone, _l3ppone;
            UsedParameter _l4one, _l4pone, _l4ppone;
            UsedParameter _l5one, _l5pone, _l5ppone;
            UsedParameter _l6one, _l6pone, _l6ppone;

            std::string _sslp_prefix(const std::string & prefix)
            {
                if ("B(*)->D(*)" == prefix)
                    return prefix;

                if ("1" == _opt_sslp_limit.value())
                    return "B(*)->D(*)";

                return prefix;
            }

        public:
            HQETFormFactorBase(const Parameters & p, const Options & o, const std::string & prefix) :
                _model(Model::make("SM", p, o)),
                _mBar(p[prefix + "::mBar@HQET"], *this),
                _a(p[prefix + "::a@HQET"], *this),
                _opt_lp_model(o, "model-lp", { "power-series", "exponential" }, "power-series"),
                _opt_lp_zorder(o, "z-order-lp", { "2", "3", "4", "5" }, "3"),
                _enable_lp_z3(1.0 ? _opt_lp_zorder.value() >= "3" : 0.0),
                _enable_lp_z4(1.0 ? _opt_lp_zorder.value() >= "4" : 0.0),
                _enable_lp_z5(1.0 ? _opt_lp_zorder.value() >= "5" : 0.0),
                _opt_slp_zorder(o, "z-order-slp", { "1", "2" }, "2"),
                _enable_slp_z2(1.0 ? _opt_slp_zorder.value() >= "2" : 0.0),
                _opt_sslp_zorder(o, "z-order-sslp", { "0", "1", "2" }, "1"),
                _enable_sslp_z1(1.0 ? _opt_sslp_zorder.value() >= "1" : 0.0),
                _enable_sslp_z2(1.0 ? _opt_sslp_zorder.value() >= "2" : 0.0),
                _opt_sslp_limit(o, "SU3F-limit-sslp", { "0", "1" }, "0"),
                _xipone(p[prefix + "::xi'(1)@HQET"], *this),
                _xippone(p[prefix + "::xi''(1)@HQET"], *this),
                _xipppone(p[prefix + "::xi'''(1)@HQET"], *this),
                _xippppone(p[prefix + "::xi''''(1)@HQET"], *this),
                _xipppppone(p[prefix + "::xi'''''(1)@HQET"], *this),
                _chi2one(p[prefix + "::chi_2(1)@HQET"], *this),
                _chi2pone(p[prefix + "::chi_2'(1)@HQET"], *this),
                _chi2ppone(p[prefix + "::chi_2''(1)@HQET"], *this),
                _chi3pone(p[prefix + "::chi_3'(1)@HQET"], *this),
                _chi3ppone(p[prefix + "::chi_3''(1)@HQET"], *this),
                _etaone(p[prefix + "::eta(1)@HQET"], *this),
                _etapone(p[prefix + "::eta'(1)@HQET"], *this),
                _etappone(p[prefix + "::eta''(1)@HQET"], *this),
                _l1one(p[_sslp_prefix(prefix) + "::l_1(1)@HQET"], *this),
                _l1pone(p[_sslp_prefix(prefix) + "::l_1'(1)@HQET"], *this),
                _l1ppone(p[_sslp_prefix(prefix) + "::l_1''(1)@HQET"], *this),
                _l2one(p[_sslp_prefix(prefix) + "::l_2(1)@HQET"], *this),
                _l2pone(p[_sslp_prefix(prefix) + "::l_2'(1)@HQET"], *this),
                _l2ppone(p[_sslp_prefix(prefix) + "::l_2''(1)@HQET"], *this),
                _l3one(p[_sslp_prefix(prefix) + "::l_3(1)@HQET"], *this),
                _l3pone(p[_sslp_prefix(prefix) + "::l_3'(1)@HQET"], *this),
                _l3ppone(p[_sslp_prefix(prefix) + "::l_3''(1)@HQET"], *this),
                _l4one(p[_sslp_prefix(prefix) + "::l_4(1)@HQET"], *this),
                _l4pone(p[_sslp_prefix(prefix) + "::l_4'(1)@HQET"], *this),
                _l4ppone(p[_sslp_prefix(prefix) + "::l_4''(1)@HQET"], *this),
                _l5one(p[_sslp_prefix(prefix) + "::l_5(1)@HQET"], *this),
                _l5pone(p[_sslp_prefix(prefix) + "::l_5'(1)@HQET"], *this),
                _l5ppone(p[_sslp_prefix(prefix) + "::l_5''(1)@HQET"], *this),
                _l6one(p[_sslp_prefix(prefix) + "::l_6(1)@HQET"], *this),
                _l6pone(p[_sslp_prefix(prefix) + "::l_6'(1)@HQET"], *this),
                _l6ppone(p[_sslp_prefix(prefix) + "::l_6''(1)@HQET"], *this)
            {
                using std::placeholders::_1;

                if (_opt_lp_model.value() == "exponential")
                {
                    _xi = [=](const double & q2) -> double
                    {
                        return _xi_exponential(q2);
                    };
                }
                else
                {
                    _xi = [=](const double & q2) -> double
                    {
                        return _xi_power_series(q2);
                    };
                }
            }

            ~HQETFormFactorBase() = default;

        protected:
            /*
             * HQET parameters following [BLPR2017]
             */
            inline double _mu() const { return 2.31; } // mu^2 = m_b * m_c
            inline double _alpha_s() const { return 0.26; }
            inline double _m_b_1S() const { return 4.71; }
            inline double _m_b_pole() const { return _m_b_1S() * (1 + 2.0 / 9.0 * power_of<2>(_alpha_s())); }
            inline double _m_c_pole() const { return _m_b_pole() - 3.40; }
            inline double _lambda_1() const { return -0.30; }
            inline double _LambdaBar() const { return _mBar - _m_b_pole() + _lambda_1() / (2.0 * _m_b_1S()); }

            /*
             * Interface to Process_-specific kinematics.
             */
            virtual double _w(const double & q2) const = 0;
            virtual double _q2(const double & w) const = 0;

            /*
             * Isgur-Wise functions
             */
            double _zw(const double & w) const
            {
                return (std::sqrt(w + 1.0) - std::sqrt(2.0) * _a()) / (std::sqrt(w + 1.0) + std::sqrt(2.0) * _a());
            }

            double _z(const double & q2) const
            {
                const double w = _w(q2);

                return _zw(w);
            }
            // uses a power series ansatz
            double _xi_power_series(const double & q2) const
            {
                const double a = _a(), a2 = a * a, a3 = a * a2, a4 = a2 * a2, a5 = a3 * a2;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_z(q2) - z_0);
                const double z2   =  z *  z;
                const double z3   = z2 *  z * _enable_lp_z3;
                const double z4   = z2 * z2 * _enable_lp_z4;
                const double z5   = z3 * z2 * _enable_lp_z5;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2
                                  + (2.0 +       a) * pow(1.0 + a, 4) / (2.0 * a3) * z3
                                  + (5.0 + 3.0 * a) * pow(1.0 + a, 5) / (8.0 * a4) * z4
                                  + (3.0 + 2.0 * a) * pow(1.0 + a, 6) / (8.0 * a5) * z5;

                const double wm12 =   4.0                  * pow(1.0 + a, 4) / a2         * z2
                                  + ( 6.0 +  2.0 * a     ) * pow(1.0 + a, 5) / a3         * z3
                                  + (25.0 + 14.0 * a + a2) * pow(1.0 + a, 6) / (4.0 * a4) * z4
                                  + (11.0 +  8.0 * a + a2) * pow(1.0 + a, 7) / (2.0 * a5) * z5;

                const double wm13 =   8.0                  * pow(1.0 + a, 6) / a3         * z3
                                  + (18.0 +  6.0 * a     ) * pow(1.0 + a, 7) / a4         * z4
                                  + (51.0 + 30.0 * a + a2) * pow(1.0 + a, 8) / (2.0 * a5) * z5;

                const double wm14 =  16.0             * pow(1.0 + a, 8) / a4 * z4
                                  + (48.0 + 16.0 * a) * pow(1.0 + a, 9) / a5 * z5;

                const double wm15 = 32.0 * pow(1.0 + a, 5) / a5 * z5;

                return 1.0
                    + _xipone             * wm11
                    + _xippone    / 2.0   * wm12
                    + _xipppone   / 6.0   * wm13
                    + _xippppone  / 24.0  * wm14
                    + _xipppppone / 120.0 * wm15;
            }

            // uses an exponential ansatz and expands in (w-1) first, then in z*
            double _xi_exponential(const double & q2) const
            {
                const double a = _a(), a2 = a * a, a3 = a * a2, a4 = a2 * a2, a5 = a3 * a2;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_z(q2) - z_0);
                const double z2   =  z *  z;
                const double z3   = z2 *  z * _enable_lp_z3;
                const double z4   = z2 * z2 * _enable_lp_z4;
                const double z5   = z3 * z2 * _enable_lp_z5;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2
                                  + (2.0 +       a) * pow(1.0 + a, 4) / (2.0 * a3) * z3
                                  + (5.0 + 3.0 * a) * pow(1.0 + a, 5) / (8.0 * a4) * z4
                                  + (3.0 + 2.0 * a) * pow(1.0 + a, 6) / (8.0 * a5) * z5;

                const double wm12 =   4.0                  * pow(1.0 + a, 4) / a2         * z2
                                  + ( 6.0 +  2.0 * a     ) * pow(1.0 + a, 5) / a3         * z3
                                  + (25.0 + 14.0 * a + a2) * pow(1.0 + a, 6) / (4.0 * a4) * z4
                                  + (11.0 +  8.0 * a + a2) * pow(1.0 + a, 7) / (2.0 * a5) * z5;

                const double wm13 =   8.0                  * pow(1.0 + a, 6) / a3         * z3
                                  + (18.0 +  6.0 * a     ) * pow(1.0 + a, 7) / a4         * z4
                                  + (51.0 + 30.0 * a + a2) * pow(1.0 + a, 8) / (2.0 * a5) * z5;

                const double wm14 =  16.0             * pow(1.0 + a, 8) / a4 * z4
                                  + (48.0 + 16.0 * a) * pow(1.0 + a, 9) / a5 * z5;

                const double wm15 = 32.0 * pow(1.0 + a, 5) / a5 * z5;

                return (1.0
                    + _xipone              * wm11
                    - _xipone              * wm12
                    + _xipone * 2.0 /  3.0 * wm13
                    - _xipone       /  3.0 * wm14
                    + _xipone * 2.0 / 15.0 * wm15)
                    * (1.0 + _xippone      * wm11);
            }

            double _chi2(const double & q2) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_z(q2) - z_0);
                const double z2   =  z * z * _enable_slp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _chi2one + _chi2pone * wm11 + _chi2ppone / 2.0 * wm12;
            }

            double _chi3(const double & q2) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_z(q2) - z_0);
                const double z2   =  z * z * _enable_slp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return 0.0 + _chi3pone * wm11 + _chi3ppone / 2.0 * wm12;
            }

            double _eta(const double & q2) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_z(q2) - z_0);
                const double z2   =  z * z * _enable_slp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _etaone + _etapone * wm11 + _etappone / 2.0 * wm12;
            }

            /*
             * Auxilliary functions for the HQET Wilson coefficients
             *
             * We use a fixed scale mu = sqrt(m_b * m_c), with m_b = 4.2 and m_c = 1.27,
             * which yields mu = 2.31 GeV.
             */

            inline double _wz(const double & z) const
            {
                return 0.5 * (z + 1.0 / z);
            }

            inline double _wp(const double & w) const { return w + std::sqrt(w * w - 1.0); }
            inline double _wm(const double & w) const { return w - std::sqrt(w * w - 1.0); }

            double _r(const double & w) const
            {
                if (w < 1.0)
                    return std::numeric_limits<double>::quiet_NaN();

                if (w - 1.0 < 1.0e-5)
                    return 1.0 - (w - 1.0) / 3.0;

                return std::log(_wp(w)) / std::sqrt(w * w - 1.0);
            }

            inline double _Omega(const double & w, const double & z) const
            {
                if (w < 1.0)
                    return std::numeric_limits<double>::quiet_NaN();

                const double lnz = std::log(z);

                if (w - 1.0 < 1.0e-5)
                    return -1.0 - (1.0 + z) / (1.0 - z) * lnz;

                const double wm = _wm(w);
                const double wp = _wp(w);

                const complex<double> li2wmz = dilog(1.0 - wm * z);
                const complex<double> li2wpz = dilog(1.0 - wp * z);
                const complex<double> li2wm2 = dilog(1.0 - wm * wm);
                const complex<double> li2wp2 = dilog(1.0 - wp * wp);

                return w * real(2.0 * (li2wmz - li2wpz) + li2wp2 - li2wm2) / (2.0 * std::sqrt(w * w - 1.0))
                    - w * _r(w) * lnz + 1.0;
            }

            /* Power corrections */
            inline double _l1(const double & w) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
                const double z2   =  z * z * _enable_sslp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _l1one + _l1pone * wm11 + _l1ppone / 2.0 * wm12;
            }
            inline double _l2(const double & w) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
                const double z2   =  z * z * _enable_sslp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _l2one + _l2pone * wm11 + _l2ppone / 2.0 * wm12;
            }
            inline double _l3(const double & w) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
                const double z2   =  z * z * _enable_sslp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _l3one + _l3pone * wm11 + _l3ppone / 2.0 * wm12;
            }
            inline double _l4(const double & w) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
                const double z2   =  z * z * _enable_sslp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _l4one + _l4pone * wm11 + _l4ppone / 2.0 * wm12;
            }
            inline double _l5(const double & w) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
                const double z2   =  z * z * _enable_sslp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _l5one + _l5pone * wm11 + _l5ppone / 2.0 * wm12;
            }
            inline double _l6(const double & w) const
            {
                const double a = _a(), a2 = a * a;

                // expansion in z around z_0
                const double  z_0 = (1.0 - a) / (1.0 + a);
                const double  z   = (_zw(w) - z_0) * _enable_sslp_z1;
                const double z2   =  z * z * _enable_sslp_z2;

                const double wm11 =  2.0            * pow(1.0 + a, 2) / a          * z
                                  + (3.0 +       a) * pow(1.0 + a, 3) / (2.0 * a2) * z2;

                const double wm12 =   4.0           * pow(1.0 + a, 4) / a2         * z2;

                return _l6one + _l6pone * wm11 + _l6ppone / 2.0 * wm12;
            }

            /* Wilson Coefficients */

            inline double _CS(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = 2.0 * z * (w - wz) * _Omega(w, z);
                result -= (w - 1.0) * (z + 1.0) * (z + 1.0) * _r(w);
                result += (z2 - 1.0) * lnz;

                return result / (3.0 * z * (w - wz));
            }

            inline double _CP(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = 2.0 * z * (w - wz) * _Omega(w, z);
                result -= (w + 1.0) * (z - 1.0) * (z - 1.0) * _r(w);
                result += (z2 - 1.0) * lnz;

                return result / (3.0 * z * (w - wz));
            }

            inline double _CV1(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = 2.0 * (w + 1.0) * ((3.0 * w - 1.0) * z - z2 - 1.0) * _r(w);
                result += (12.0 * z * (wz - w) - (z2 - 1.0) * lnz);
                result += 4.0 * z * (w - wz) * _Omega(w, z);

                return result / (6.0 * z * (w - wz));
            }

            inline double _CV2(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = ((4.0 * w2 + 2.0 * w) * z2 - (2.0 * w2 + 5.0 * w - 1.0) * z - (1.0 + w) * z3 + 2.0) * _r(w);
                result += z * (2.0 * (z - 1.0) * (wz - w) + (z2 - (4.0 * w - 2.0) * z + (-2.0 * w + 3)) * lnz);

                return -1.0 * result / (6.0 * z2 * power_of<2>(w - wz));
            }

            inline double _CV3(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (-2.0 * z3 + (2.0 * w2 + 5.0 * w - 1.0) * z2 - (4.0 * w2 + 2.0 * w) * z + w + 1.0) * _r(w);
                result += 2.0 * z * (z - 1.0) * (wz - w) + ((-2.0 * w + 3.0) * z2 + (-4.0 * w + 2.0) * z + 1.0) * lnz;

                return +1.0 * result / (6.0 * z * power_of<2>(w - wz));
            }

            inline double _CA1(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = 2.0 * (w - 1.0) * ((3.0 * w + 1.0) * z - z2 - 1.0) * _r(w);
                result += (12.0 * z * (wz - w) - (z2 - 1.0) * lnz);
                result += 4.0 * z * (w - wz) * _Omega(w, z);

                return result / (6.0 * z * (w - wz));
            }

            inline double _CA2(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = ((4.0 * w2 - 2.0 * w) * z2 + (2.0 * w2 - 5.0 * w - 1.0) * z + (1.0 - w) * z3 + 2.0) * _r(w);
                result += z * (2.0 * (z + 1.0) * (wz - w) + (z2 - (4.0 * w + 2.0) * z + (2.0 * w + 3)) * lnz);

                return -1.0 * result / (6.0 * z2 * power_of<2>(w - wz));
            }

            inline double _CA3(const double & w, const double & z) const
            {
                const double z2  = z * z, z3 = z2 * z;
                const double w2  = w * w;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (2.0 * z3 + (2.0 * w2 - 5.0 * w - 1.0) * z2 + (4.0 * w2 - 2.0 * w) * z - w + 1.0) * _r(w);
                result += 2.0 * z * (z + 1.0) * (wz - w) - ((2.0 * w + 3.0) * z2 - (4.0 * w + 2.0) * z + 1.0) * lnz;

                return +1.0 * result / (6.0 * z * power_of<2>(w - wz));
            }

            inline double _CT1(const double & w, const double & z) const
            {
                const double z2  = z * z;
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (w - 1.0) * ((4.0 * w + 2.0) * z - z2 - 1.0) * _r(w);
                result += 6.0 * z * (wz - w) - (z2 - 1.0) * lnz;
                result += 2.0 * z * (w - wz) * _Omega(w, z);

                return +1.0 / (3.0 * z * (w - wz)) * result;
            }

            inline double _CT2(const double & w, const double & z) const
            {
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (1.0 - w * z) * _r(w) + z * lnz;

                return +2.0 / (3.0 * z * (w - wz)) * result;
            }

            inline double _CT3(const double & w, const double & z) const
            {
                const double wz  = _wz(z);
                const double lnz = std::log(z);

                double result = (w - z) * _r(w) + lnz;

                return +2.0 / (3.0 * (w - wz)) * result;
            }
    };

    class HQETIsgurWiseFunctionParameters :
        public virtual HQETFormFactorBase
    {
        private:
            /*
             * Kinematics
             */
            virtual double _w(const double & /*q2*/) const override
            {
                throw InternalError("Kinematic function _w() should not be used within HQETIsgurWiseFunctions");
                return 1.0;
            }

            virtual double _q2(const double & /*w*/) const override
            {
                throw InternalError("Kinematic function _q2() should not be used within HQETIsgurWiseFunctions");
                return 0.0;
            }

            std::string _prefix(const Options & o)
            {
                if (! o.has("q"))
                {
                    return "B(*)->D(*)";
                }

                std::string q = o.get("q", "u");

                switch (q[0])
                {
                    case 'u':
                    case 'd':
                        return "B(*)->D(*)";
                        break;

                    case 's':
                        return "B_s(*)->D_s(*)";
                        break;

                    default:
                        throw InvalidOptionValueError("q", q, "u,d,s");
                }
            }

        public:
            HQETIsgurWiseFunctionParameters(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o, _prefix(o))
            {
            }

            ~HQETIsgurWiseFunctionParameters() = default;

            double xipone()   const { return _xipone();   }
            double xippone()  const { return _xippone();  }
            double xipppone() const { return _xipppone();  }
            double chi2one()  const { return _chi2one();  }
            double chi2pone() const { return _chi2pone(); }
            double chi3pone() const { return _chi3pone(); }
            double etaone()   const { return _etaone();   }
            double etapone()  const { return _etapone();  }
    };

    template <typename Process_> class HQETFormFactors<Process_, PToP> :
        public virtual HQETFormFactorBase,
        public virtual FormFactors<PToP>
    {
        private:
            UsedParameter _m_B;
            UsedParameter _m_P;

            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override
            {
                const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
                const double m_P = this->_m_P(), m_P2 = power_of<2>(m_P);

                return (m_B2 + m_P2 - q2) / (2.0 * m_B * m_P);
            }

            virtual double _q2(const double & w) const override
            {
                const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
                const double m_P = this->_m_P(), m_P2 = power_of<2>(m_P);

                return m_B2 + m_P2 - 2.0 * m_B * m_P * w;
            }

            /* HQET form factors h_i */

            double _h_p(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;

                double result = 1.0 + as * (_CV1(w, z) + (w + 1.0) / 2.0 * (_CV2(w, z) + _CV3(w, z)));
                result += eps_c * (L1);
                result += eps_b * (L1);
                result += eps_c * eps_c * _l1(w);

                return result * xi;
            }

            double _h_m(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L4 = 2.0 * eta - 1.0;

                double result = (0.0 + as * (w + 1.0) / 2.0 * (_CV2(w, z) - _CV3(w, z)));
                result += eps_c * L4;
                result -= eps_b * L4;
                result += eps_c * eps_c * _l4(w);

                return result * xi;
            }

            double _h_S(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;

                double result = (1.0 + as * _CS(w, z));
                result += eps_c * (L1 - (w - 1.0) / (w + 1.0) * L4);
                result += eps_b * (L1 - (w - 1.0) / (w + 1.0) * L4);
                result += eps_c * eps_c * (_l1(w) - (w - 1.0) / (w + 1.0) * _l4(w));

                return result * xi;
            }

            double _h_T(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;

                double result = 1.0 + as * (_CT1(w, z) - _CT2(w, z) + _CT3(w, z));
                result += eps_c * (L1 - L4);
                result += eps_b * (L1 - L4);
                result += eps_c * eps_c * (_l1(w) - _l4(w));

                return result * xi;
            }

        public:
            HQETFormFactors(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o, Process_::hqe_prefix),
                _m_B(p[Process_::name_B], *static_cast<ParameterUser *>(this)),
                _m_P(p[Process_::name_P], *static_cast<ParameterUser *>(this))
            {
            }

            ~HQETFormFactors() = default;

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options)
            {
                return new HQETFormFactors<Process_, PToP>(parameters, options);
            }

            virtual double f_p(const double & q2) const
            {
                const double r = _m_P / _m_B;

                // cf. [FKKM2008], eq. (22)
                return 1.0 / (2.0 * sqrt(r)) * ((1.0 + r) * _h_p(q2) - (1.0 - r) * _h_m(q2));
            }

            double f_m(const double & q2) const
            {
                const double r = _m_P / _m_B;

                // cf. [FKKM2008], eq. (22)
                return 1.0 / (2.0 * sqrt(r)) * ((1.0 + r) * _h_m(q2) - (1.0 - r) * _h_p(q2));
            }

            virtual double f_0(const double & q2) const
            {
                // We do not use the relation between f_0 and the (scale-dependent) h_S.
                return f_p(q2) + q2 / (_m_B * _m_B - _m_P * _m_P) * f_m(q2);
            }

            virtual double f_t(const double & q2) const
            {
                const double r = _m_P / _m_B;

                // cf. [BJvD2019], eq. (A7)
                return (1.0 + r) / (2.0 * sqrt(r)) * _h_T(q2);
            }

            virtual double f_plus_T(const double & q2) const
            {
                return f_t(q2) * q2 / _m_B / (_m_B + _m_P);
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    const double m_b = _m_b_pole();
                    const double m_c = _m_c_pole();
                    const double z   = m_c / m_b;
                    const double wz  = _wz(z);

                    results.add(Diagnostics::Entry{ z,  "z = m_c_pole / m_b_pole" });
                    results.add(Diagnostics::Entry{ wz, "w_z"                     });
                }

                // Switches
                {
                    results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
                    results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
                }

                // z
                {
                    results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
                }

                // xi
                {
                    results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
                }

                // chi2
                {
                    results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
                }

                // chi3
                {
                    results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
                }

                // eta
                {
                    results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
                }

                // r(w)
                {
                    results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
                    results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
                    results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
                    results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
                    results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
                }

                // Omega(w, z = 0.25)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
                }

                // Omega(w, z = 0.20)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
                }

                // WCs at w = 1.2, z = 0.20
                {
                    results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
                }

                // WCs at w = 1.0, z = 0.25
                {
                    results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
                }

                // HQET definition of the form factors
                {
                    results.add(Diagnostics::Entry{ _h_p(_q2(1.4)), "h_+(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_m(_q2(1.4)), "h_-(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_T(_q2(1.4)), "h_T(w = 1.4)" });

                    results.add(Diagnostics::Entry{ _h_p(_q2(1.2)), "h_+(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_m(_q2(1.2)), "h_-(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_T(_q2(1.2)), "h_T(w = 1.2)" });

                    results.add(Diagnostics::Entry{ _h_p(_q2(1.0)), "h_+(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_m(_q2(1.0)), "h_-(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_T(_q2(1.0)), "h_T(w = 1.0)" });
                }

                return results;
            }
    };

    template <typename Process_> class HQETFormFactors<Process_, PToV> :
        public virtual HQETFormFactorBase,
        public virtual FormFactors<PToV>
    {
        private:
            UsedParameter _m_B;
            UsedParameter _m_V;

            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override
            {
                const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
                const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);

                return (m_B2 + m_V2 - q2) / (2.0 * m_B * m_V);
            }

            virtual double _q2(const double & w) const override
            {
                const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
                const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);

                return m_B2 + m_V2 - 2.0 * m_B * m_V * w;
            }

            /* HQET form factors h_i */

            double _h_a1(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (1.0 + as * _CA1(w, z));
                result += eps_c * (L2 - L5 * (w - 1.0) / (w + 1.0));
                result += eps_b * (L1 - L4 * (w - 1.0) / (w + 1.0));
                result += eps_c * eps_c * (_l2(w) - (w - 1.0) / (w + 1.0) * _l5(w));

                return result * xi;
            }

            double _h_a2(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 + as * _CA2(w, z));
                result += eps_c * (L3 + L6);
                result += eps_c * eps_c * (_l3(w) + _l6(w));

                return result * xi;
            }

            double _h_a3(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L3 = 4.0 * chi2;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (1.0 + as * (_CA1(w, z) +_CA3(w, z)));
                result += eps_c * (L2 - L3 + L6 - L5);
                result += eps_b * (L1 - L4);
                result += eps_c * eps_c * (_l2(w) - _l3(w) + _l6(w) - _l5(w));

                return result * xi;
            }

            double _h_v(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (1.0 + as * _CV1(w, z));
                result += eps_c * (L2 - L5);
                result += eps_b * (L1 - L4);
                result += eps_c * eps_c * (_l2(w) - _l5(w));

                return result * xi;
            }

            double _h_t1(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;

                double result = (1.0 + as * (_CT1(w, z) + (w - 1.0) / 2.0 * (_CT2(w, z) - _CT3(w, z))));
                result += eps_c * L2;
                result += eps_b * L1;
                result += eps_c * eps_c * _l2(w);

                return result * xi;
            }

            double _h_t2(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (0.0 + as * (w + 1.0) / 2.0 * (_CT2(w, z) + _CT3(w, z)));
                result += eps_c * L5;
                result -= eps_b * L4;
                result += eps_c * eps_c * _l5(w);

                return result * xi;
            }

            double _h_t3(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 + as * _CT2(w, z));
                result += eps_c * (L6 - L3);
                result += eps_c * eps_c * (_l6(w) - _l3(w));

                return result * xi;
            }

        public:
            HQETFormFactors(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o, Process_::hqe_prefix),
                _m_B(p[Process_::name_B], *static_cast<ParameterUser *>(this)),
                _m_V(p[Process_::name_V], *static_cast<ParameterUser *>(this))
            {
            }

            ~HQETFormFactors() = default;

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options)
            {
                return new HQETFormFactors(parameters, options);
            }

            virtual double v(const double & q2) const
            {
                const double r = _m_V / _m_B;

                // cf. [FKKM2008], eq. (22)
                return (1.0 + r) / 2.0 / sqrt(r) * _h_v(q2);
            }

            virtual double a_0(const double & q2) const
            {
                const double r = _m_V / _m_B;
                const double w = _w(q2);

                return 1.0 / (2.0 * sqrt(r)) * ((1.0 + w) * _h_a1(q2) + (r * w - 1.0) * _h_a2(q2) + (r - w) * _h_a3(q2));
                // cf. [FKKM2008], eq. (22)
                //const double a_30 = (1.0 + r * r - 2.0 * r * w) / (4.0 * r * sqrt(r)) * (r * _h_a2(q2) - _h_a3(q2));
                //return a_3(q2) - a_30;
            }

            virtual double a_1(const double & q2) const
            {
                const double r = _m_V / _m_B;
                const double w = _w(q2);

                // cf. [FKKM2008], eq. (22)
                return sqrt(r) * (1.0 + w) / (1.0 + r) * _h_a1(q2);
            }

            virtual double a_2(const double & q2) const
            {
                const double r = _m_V / _m_B;

                // cf. [FKKM2008], eq. (22)
                return (1.0 + r) / (2.0 * sqrt(r)) * (r * _h_a2(q2) + _h_a3(q2));
            }

            double a_3(const double & q2) const
            {
                const double r = _m_V / _m_B;

                // cf. [FKKM2008], below eq. (6)
                return ((1.0 + r) * a_1(q2) - (1.0 - r) * a_2(q2)) / (2.0 * r);
            }

            virtual double a_12(const double & q2) const
            {
                const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
                const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);
                const double lambda = eos::lambda(m_B2, m_V2, q2);

                double result = (m_B + m_V) * (m_B + m_V) * (m_B2 - m_V2 - q2) * a_1(q2) - lambda * a_2(q2);
                result /= 16.0 * m_B * m_V2 * (m_B + m_V);

                return result;
            }

            virtual double t_1(const double & q2) const
            {
                const double r = _m_V / _m_B;
                const double w = _w(q2);

                return -1.0 / (2.0 * sqrt(r)) * ((1.0 - r) * _h_t2(q2) - (1.0 + r) * _h_t1(q2));
            }

            virtual double t_2(const double & q2) const
            {
                const double r = _m_V / _m_B;
                const double w = _w(q2);

                return +1.0 / (2.0 * sqrt(r)) * (2.0 * r * (w + 1.0) / (1.0 + r) * _h_t1(q2) - 2.0 * r * (w - 1.0) / (1.0 - r) * _h_t2(q2));
            }

            virtual double t_3(const double & q2) const
            {
                const double r = _m_V / _m_B;

                return +1.0 / (2.0 * sqrt(r)) * ((1.0 - r) * _h_t1(q2) - (1.0 + r) * _h_t2(q2) + (1.0 - r * r) * _h_t3(q2));
            }

            virtual double t_23(const double & q2) const
            {
                const double m_B = this->_m_B(), m_B2 = power_of<2>(m_B);
                const double m_V = this->_m_V(), m_V2 = power_of<2>(m_V);
                const double lambda = eos::lambda(m_B2, m_V2, q2);

                return ((m_B2 - m_V2) * (m_B2 + 3.0 * m_V2 - q2) * t_2(q2) - lambda * t_3(q2)) / (8.0 * m_B * m_V2 * (m_B - m_V));
            }

            virtual double f_perp(const double &) const
            {
                return 0.;
            }

            virtual double f_para(const double &) const
            {
                return 0.;
            }

            virtual double f_long(const double &) const
            {
                return 0.;
            }

            virtual double f_perp_T(const double &) const
            {
                return 0.;
            }

            virtual double f_para_T(const double &) const
            {
                return 0.;
            }

            virtual double f_long_T(const double &) const
            {
                return 0.;
            }

            virtual double f_long_T_Normalized(const double &) const
            {
                return 0.;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    const double m_b = _m_b_pole();
                    const double m_c = _m_c_pole();
                    const double z   = m_c / m_b;
                    const double wz  = _wz(z);

                    results.add(Diagnostics::Entry{ z,  "z = m_c_pole / m_b_pole" });
                    results.add(Diagnostics::Entry{ wz, "w_z"                     });
                }

                // Switches
                {
                    results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
                    results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
                }

                // z
                {
                    results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
                }

                // xi
                {
                    results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
                }

                // chi2
                {
                    results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
                }

                // chi3
                {
                    results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
                }

                // eta
                {
                    results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
                }

                // r(w)
                {
                    results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
                    results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
                    results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
                    results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
                    results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
                }

                // Omega(w, z = 0.25)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
                }

                // Omega(w, z = 0.20)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
                }

                // WCs at w = 1.2, z = 0.20
                {
                    results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
                }

                // WCs at w = 1.0, z = 0.25
                {
                    results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
                }

                // HQET definition of the form factors
                {
                    results.add(Diagnostics::Entry{ _h_a1(_q2(1.4)), "h_A1(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_a2(_q2(1.4)), "h_A2(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_a3(_q2(1.4)), "h_A3(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_v (_q2(1.4)), "h_V (w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_t1(_q2(1.4)), "h_T1(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_t2(_q2(1.4)), "h_T2(w = 1.4)" });
                    results.add(Diagnostics::Entry{ _h_t3(_q2(1.4)), "h_T3(w = 1.4)" });

                    results.add(Diagnostics::Entry{ _h_a1(_q2(1.2)), "h_A1(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_a2(_q2(1.2)), "h_A2(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_a3(_q2(1.2)), "h_A3(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_v (_q2(1.2)), "h_V (w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_t1(_q2(1.2)), "h_T1(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_t2(_q2(1.2)), "h_T2(w = 1.2)" });
                    results.add(Diagnostics::Entry{ _h_t3(_q2(1.2)), "h_T3(w = 1.2)" });

                    results.add(Diagnostics::Entry{ _h_a1(_q2(1.0)), "h_A1(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_a2(_q2(1.0)), "h_A2(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_a3(_q2(1.0)), "h_A3(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_v (_q2(1.0)), "h_V (w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_t1(_q2(1.0)), "h_T1(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_t2(_q2(1.0)), "h_T2(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _h_t3(_q2(1.0)), "h_T3(w = 1.0)" });
                }

                return results;
            }
    };

    template <typename Process_> class HQETFormFactors<Process_, VToP> :
        public HQETFormFactorBase,
        public FormFactors<VToP>
    {
        private:
            UsedParameter _m_Bst;
            UsedParameter _m_P;

            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override
            {
                const double m_Bst = this->_m_Bst(), m_Bst2 = power_of<2>(m_Bst);
                const double m_P   = this->_m_P(),   m_P2   = power_of<2>(m_P);

                return (m_Bst2 + m_P2 - q2) / (2.0 * m_Bst * m_P);
            }

            virtual double _q2(const double & w) const override
            {
                const double m_Bst = this->_m_Bst(), m_Bst2 = power_of<2>(m_Bst);
                const double m_P   = this->_m_P(),   m_P2   = power_of<2>(m_P);

                return m_Bst2 + m_P2 - 2.0 * m_Bst * m_P * w;
            }

        public:
            HQETFormFactors(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o, Process_::hqe_prefix),
                _m_Bst(p[Process_::name_Bst], *static_cast<ParameterUser *>(this)),
                _m_P(p[Process_::name_P], *static_cast<ParameterUser *>(this))
            {
            }

            ~HQETFormFactors() = default;

            static FormFactors<VToP> * make(const Parameters & parameters, const Options & options)
            {
                return new HQETFormFactors<Process_, VToP>(parameters, options);
            }

            /* HQET form factors h_i */

            virtual double h_abar_1(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (1.0 + as * _CA1(w, z));
                result += eps_c * (L1 - L4 * (w - 1.0) / (w + 1.0));
                result += eps_b * (L2 - L5 * (w - 1.0) / (w + 1.0));
                result += eps_c * eps_c * (_l1(w) - _l4(w) * (w - 1.0) / (w + 1.0));

                return result * xi;
            }

            virtual double h_abar_2(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 - as * _CA3(w, z));
                result += eps_b * (L3 + L6);

                return result * xi;
            }

            virtual double h_abar_3(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L3 = 4.0 * chi2;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (1.0 + as * (_CA1(w, z) - _CA2(w, z)));
                result += eps_b * (L2 - L3 + L6 - L5);
                result += eps_c * (L1 - L4);
                result += eps_c * eps_c * (_l1(w) - _l4(w));

                return result * xi;
            }

            virtual double h_vbar(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L1 = -4.0 * (w - 1.0) * chi2 + 12.0 * chi3;
                const double L2 = -4.0 * chi3;
                const double L4 = 2.0 * eta - 1.0;
                const double L5 = -1.0;

                double result = (1.0 + as * _CV1(w, z));
                result += eps_b * (L2 - L5);
                result += eps_c * (L1 - L4);
                result += eps_c * eps_c * (_l1(w) - _l4(w));

                return result * xi;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    const double m_b = _m_b_pole();
                    const double m_c = _m_c_pole();
                    const double z   = m_c / m_b;
                    const double wz  = _wz(z);

                    results.add(Diagnostics::Entry{ z,  "z = m_c / m_b" });
                    results.add(Diagnostics::Entry{ wz, "w_z"           });
                }

                // Switches
                {
                    results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
                    results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
                }

                // z
                {
                    results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
                }

                // xi
                {
                    results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
                }

                // chi2
                {
                    results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
                }

                // chi3
                {
                    results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
                }

                // eta
                {
                    results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
                }

                // r(w)
                {
                    results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
                    results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
                    results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
                    results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
                    results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
                }

                // Omega(w, z = 0.25)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
                }

                // Omega(w, z = 0.20)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
                }

                // WCs at w = 1.2, z = 0.20
                {
                    results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
                }

                // WCs at w = 1.0, z = 0.25
                {
                    results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
                }

                // HQET definition of the form factors
                {
                    results.add(Diagnostics::Entry{ h_abar_1(_q2(1.4)), "h_Abar1(w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_abar_2(_q2(1.4)), "h_Abar2(w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_abar_3(_q2(1.4)), "h_Abar3(w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_vbar (_q2(1.4)),  "h_Vbar (w = 1.4)" });

                    results.add(Diagnostics::Entry{ h_abar_1(_q2(1.2)), "h_Abar1(w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_abar_2(_q2(1.2)), "h_Abar2(w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_abar_3(_q2(1.2)), "h_Abar3(w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_vbar (_q2(1.2)),  "h_Vbar (w = 1.2)" });

                    results.add(Diagnostics::Entry{ h_abar_1(_q2(1.0)), "h_Abar1(w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_abar_2(_q2(1.0)), "h_Abar2(w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_abar_3(_q2(1.0)), "h_Abar3(w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_vbar (_q2(1.0)),  "h_Vbar (w = 1.0)" });
                }

                return results;
            }
    };

    template <typename Process_> class HQETFormFactors<Process_, VToV> :
        public HQETFormFactorBase,
        public FormFactors<VToV>
    {
        private:
            /*
             * Kinematics
             */
            virtual double _w(const double & q2) const override
            {
                static constexpr double mV1 = Process_::mV1, mV12 = power_of<2>(Process_::mV1);
                static constexpr double mV2 = Process_::mV2, mV22 = power_of<2>(Process_::mV2);

                return (mV12 + mV22 - q2) / (2.0 * mV1 * mV2);
            }

            virtual double _q2(const double & w) const override
            {
                static constexpr double mV1 = Process_::mV1, mV12 = power_of<2>(Process_::mV1);
                static constexpr double mV2 = Process_::mV2, mV22 = power_of<2>(Process_::mV2);

                return mV12 + mV22 - 2.0 * mV1 * mV2 * w;
            }

        public:
            HQETFormFactors(const Parameters & p, const Options & o) :
                HQETFormFactorBase(p, o, Process_::hqe_prefix)
            {
            }

            ~HQETFormFactors() = default;

            static FormFactors<VToV> * make(const Parameters & parameters, const Options & options)
            {
                return new HQETFormFactors(parameters, options);
            }

            /* HQET form factors h_i */

            // vector current
            virtual double h_1(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L2 = -4.0 * chi3;

                double result = 1.0 + as * (_CV1(w, z) + (w + 1.0) / 2.0 * (_CV2(w, z) + _CV3(w, z)));
                result += eps_c * L2;
                result += eps_b * L2;
                result += eps_c * eps_c * _l2(w);

                return result * xi;
            }

            virtual double h_2(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L5 = -1.0;

                double result = as * (w + 1.0) / 2.0 * (_CV2(w, z) - _CV3(w, z));
                result += eps_c * L5;
                result -= eps_b * L5;
                result += eps_c * eps_c * _l5(w);

                return result * xi;
            }

            virtual double h_3(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L2 = -4.0 * chi3;
                const double L3 = 4.0 * chi2;
                const double L5 = -1.0;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (1.0 + as * _CV1(w, z));
                result += eps_c * (L2 + L5 + (w - 1.0) * L3 - (w + 1.0) * L6);
                result += eps_b * (L2 - L5);
                result += eps_c * eps_c * (_l2(w) + _l5(w) + (w - 1.0) * _l3(w) - (w + 1.0) * _l6(w));

                return result * xi;
            }

            virtual double h_4(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L2 = -4.0 * chi3;
                const double L3 = 4.0 * chi2;
                const double L5 = -1.0;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (1.0 + as * _CV1(w, z));
                result += eps_b * (L2 + L5 + (w - 1.0) * L3 - (w + 1.0) * L6);
                result += eps_c * (L2 - L5);
                result += eps_c * eps_c * (_l2(w) - _l5(w));

                return result * xi;
            }

            virtual double h_5(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 - as * _CV2(w, z));
                result += eps_c * (L3 - L6);
                result += eps_c * eps_c * (_l3(w) - _l6(w));

                return result * xi;
            }

            virtual double h_6(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 - as * _CV3(w, z));
                result += eps_b * (L3 - L6);

                return result * xi;
            }

            // axial current
            virtual double h_7(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double chi3 = _chi3(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L2 = -4.0 * chi3;

                double result = 1.0 + as * (_CA1(w, z) + (w - 1.0) / 2.0 * (_CA2(w, z) - _CA3(w, z)));
                result += eps_b * L2;
                result += eps_c * L2;
                result += eps_c * eps_c * _l2(w);

                return result * xi;
            }

            virtual double h_8(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);
                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L5 = -1.0;

                double result = as * (w + 1.0) / 2.0 * (_CA2(w, z) + _CA3(w, z));
                result += eps_c * L5;
                result -= eps_b * L5;
                result += eps_c * eps_c * _l5(w);

                return result * xi;
            }

            virtual double h_9(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_c = _LambdaBar() / (2.0 * m_c_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 - as * _CA2(w, z));
                result += eps_c * (L3 - L6);
                result += eps_c * eps_c * (_l3(w) - _l6(w));

                return result * xi;
            }

            virtual double h_10(const double & q2) const
            {
                const double m_b_pole = _m_b_pole();
                const double m_c_pole = _m_c_pole();

                const double w = this->_w(q2);
                const double z = m_c_pole / m_b_pole;

                const double as = _alpha_s() / M_PI;

                const double xi  = _xi(q2);
                const double eta = _eta(q2);
                const double chi2 = _chi2(q2);

                const double eps_b = _LambdaBar() / (2.0 * m_b_pole);

                // chi_1 is absorbed into def. of xi for LP and LV
                const double L3 = 4.0 * chi2;
                const double L6 = -2.0 * (1.0 + eta) / (w + 1.0);

                double result = (0.0 + as * _CA3(w, z));
                result += eps_b * (L3 - L6);

                return result * xi;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    const double m_b = _m_b_pole();
                    const double m_c = _m_c_pole();
                    const double z   = m_c / m_b;
                    const double wz  = _wz(z);

                    results.add(Diagnostics::Entry{ z,  "z = m_c / m_b" });
                    results.add(Diagnostics::Entry{ wz, "w_z"           });
                }

                // Switches
                {
                    results.add(Diagnostics::Entry{ _enable_lp_z3,  "enable LP  z^3 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z4,  "enable LP  z^4 terms" });
                    results.add(Diagnostics::Entry{ _enable_lp_z5,  "enable LP  z^5 terms" });
                    results.add(Diagnostics::Entry{ _enable_slp_z2, "enable SLP z^2 terms" });
                }

                // z
                {
                    results.add(Diagnostics::Entry{ _z(_q2(1.10)), "z(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.05)), "z(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _z(_q2(1.00)), "z(w = 1.00)" });
                }

                // xi
                {
                    results.add(Diagnostics::Entry{ _xi(_q2(2.10)), "xi(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.60)), "xi(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.10)), "xi(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.05)), "xi(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _xi(_q2(1.00)), "xi(w = 1.00)" });
                }

                // chi2
                {
                    results.add(Diagnostics::Entry{ _chi2(_q2(2.10)), "chi2(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.60)), "chi2(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.10)), "chi2(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.05)), "chi2(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi2(_q2(1.00)), "chi2(w = 1.00)" });
                }

                // chi3
                {
                    results.add(Diagnostics::Entry{ _chi3(_q2(2.10)), "chi3(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.60)), "chi3(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.10)), "chi3(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.05)), "chi3(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _chi3(_q2(1.00)), "chi3(w = 1.00)" });
                }

                // eta
                {
                    results.add(Diagnostics::Entry{ _eta(_q2(2.10)), "eta(w = 2.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.60)), "eta(w = 1.60)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.10)), "eta(w = 1.10)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.05)), "eta(w = 1.05)" });
                    results.add(Diagnostics::Entry{ _eta(_q2(1.00)), "eta(w = 1.00)" });
                }

                // r(w)
                {
                    results.add(Diagnostics::Entry{ _r(1.1),     "r(w = 1.1)"     });
                    results.add(Diagnostics::Entry{ _r(1.0007),  "r(w = 1.0007)"  });
                    results.add(Diagnostics::Entry{ _r(1.0001),  "r(w = 1.0001)"  });
                    results.add(Diagnostics::Entry{ _r(1.00005), "r(w = 1.00005)" });
                    results.add(Diagnostics::Entry{ _r(1.0),     "r(w = 1.0)"     });
                }

                // Omega(w, z = 0.25)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.25), "Omega(w = 1.1,     z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.25), "Omega(w = 1.0007,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.25), "Omega(w = 1.0001,  z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.25), "Omega(w = 1.00005, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.25), "Omega(w = 1.0,     z = 0.25)" });
                }

                // Omega(w, z = 0.20)
                {
                    results.add(Diagnostics::Entry{ _Omega(1.1,     0.20), "Omega(w = 1.1,     z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0007,  0.20), "Omega(w = 1.0007,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0001,  0.20), "Omega(w = 1.0001,  z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.00005, 0.20), "Omega(w = 1.00005, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _Omega(1.0,     0.20), "Omega(w = 1.0,     z = 0.20)" });
                }

                // WCs at w = 1.2, z = 0.20
                {
                    results.add(Diagnostics::Entry{ _CS( 1.2, 0.20), "C_{S  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CP( 1.2, 0.20), "C_{P  }(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV1(1.2, 0.20), "C_{V_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV2(1.2, 0.20), "C_{V_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CV3(1.2, 0.20), "C_{V_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA1(1.2, 0.20), "C_{A_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA2(1.2, 0.20), "C_{A_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CA3(1.2, 0.20), "C_{A_3}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT1(1.2, 0.20), "C_{T_1}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT2(1.2, 0.20), "C_{T_2}(w = 1.2, z = 0.20)" });
                    results.add(Diagnostics::Entry{ _CT3(1.2, 0.20), "C_{T_3}(w = 1.2, z = 0.20)" });
                }

                // WCs at w = 1.0, z = 0.25
                {
                    results.add(Diagnostics::Entry{ _CS( 1.0, 0.25), "C_{S  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CP( 1.0, 0.25), "C_{P  }(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV1(1.0, 0.25), "C_{V_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV2(1.0, 0.25), "C_{V_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CV3(1.0, 0.25), "C_{V_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA1(1.0, 0.25), "C_{A_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA2(1.0, 0.25), "C_{A_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CA3(1.0, 0.25), "C_{A_3}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT1(1.0, 0.25), "C_{T_1}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT2(1.0, 0.25), "C_{T_2}(w = 1.0, z = 0.25)" });
                    results.add(Diagnostics::Entry{ _CT3(1.0, 0.25), "C_{T_3}(w = 1.0, z = 0.25)" });
                }

                // HQET definition of the form factors
                {
                    results.add(Diagnostics::Entry{ h_1(_q2(1.4)),  "h_1 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_2(_q2(1.4)),  "h_2 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_3(_q2(1.4)),  "h_3 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_4(_q2(1.4)),  "h_4 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_5(_q2(1.4)),  "h_5 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_6(_q2(1.4)),  "h_6 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_7(_q2(1.4)),  "h_7 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_8(_q2(1.4)),  "h_8 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_9(_q2(1.4)),  "h_9 (w = 1.4)" });
                    results.add(Diagnostics::Entry{ h_10(_q2(1.4)), "h_10(w = 1.4)" });

                    results.add(Diagnostics::Entry{ h_1(_q2(1.2)),  "h_1 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_2(_q2(1.2)),  "h_2 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_3(_q2(1.2)),  "h_3 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_4(_q2(1.2)),  "h_4 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_5(_q2(1.2)),  "h_5 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_6(_q2(1.2)),  "h_6 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_7(_q2(1.2)),  "h_7 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_8(_q2(1.2)),  "h_8 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_9(_q2(1.2)),  "h_9 (w = 1.2)" });
                    results.add(Diagnostics::Entry{ h_10(_q2(1.2)), "h_10(w = 1.2)" });

                    results.add(Diagnostics::Entry{ h_1(_q2(1.0)),  "h_1 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_2(_q2(1.0)),  "h_2 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_3(_q2(1.0)),  "h_3 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_4(_q2(1.0)),  "h_4 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_5(_q2(1.0)),  "h_5 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_6(_q2(1.0)),  "h_6 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_7(_q2(1.0)),  "h_7 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_8(_q2(1.0)),  "h_8 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_9(_q2(1.0)),  "h_9 (w = 1.0)" });
                    results.add(Diagnostics::Entry{ h_10(_q2(1.0)), "h_10(w = 1.0)" });
                }

                return results;
            }
    };
}

#endif
