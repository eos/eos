/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Frederik Beaujean
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Keri Vos
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_IMPL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_IMPL_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/form-factors/analytic-b-to-pi.hh>
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/maths/derivative.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>

#include <array>
#include <limits>

#include <iostream> // <-- TODO: Remove!

namespace eos
{
    // Process    = B -> K^*, B -> D^*, B -> rho, B_s -> phi etc.
    // Transition = B -> V or B -> P

    /* Form Factors according to [KMPW2010] */
    template <typename Transition_> class KMPW2010FormFactors;

    template <> class KMPW2010FormFactors<PToV> :
        public FormFactors<PToV>
    {
        private:
            // fit parametrisation for P -> V according to [KMPW2010]
            UsedParameter
               _f0_V, _b1_V,
               _f0_A0, _b1_A0, _f0_A1, _b1_A1, _f0_A2, _b1_A2,
               _f0_T1, _b1_T1, _f0_T2, _b1_T2, _f0_T3, _b1_T3;

            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_Kstar, _m_Bs2_0m, _m_Bs2_1m, _m_Bs2_1p;

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

            static double ff_KMPW(const double & s, const double & f0, const double & b1, const double & m2)
            {
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                // cf. [KMPW2010], Eq. (8.8), p. 30
                return f0 / (1.0 - s / m2) * (1.0 + b1 * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

        // Remark: this is hardcoded B->K* FF's nuisance parameters (OK, because in [KMPW2010] only B->K^* calculated)
        public:
            KMPW2010FormFactors(const Parameters & p, const Options &) :
                _f0_V(p["B->K^*::F^V(0)@KMPW2010"],   *this),   _b1_V(p["B->K^*::b^V_1@KMPW2010"],    *this),
                _f0_A0(p["B->K^*::F^A0(0)@KMPW2010"], *this),   _b1_A0(p["B->K^*::b^A0_1@KMPW2010"],  *this),
                _f0_A1(p["B->K^*::F^A1(0)@KMPW2010"], *this),   _b1_A1(p["B->K^*::b^A1_1@KMPW2010"],  *this),
                _f0_A2(p["B->K^*::F^A2(0)@KMPW2010"], *this),   _b1_A2(p["B->K^*::b^A2_1@KMPW2010"],  *this),
                _f0_T1(p["B->K^*::F^T1(0)@KMPW2010"], *this),   _b1_T1(p["B->K^*::b^T1_1@KMPW2010"],  *this),
                _f0_T2(p["B->K^*::F^T2(0)@KMPW2010"], *this),   _b1_T2(p["B->K^*::b^T2_1@KMPW2010"],  *this),
                _f0_T3(p["B->K^*::F^T3(0)@KMPW2010"], *this),   _b1_T3(p["B->K^*::b^T3_1@KMPW2010"],  *this)
            {
            }

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options)
            {
                return new KMPW2010FormFactors(parameters, options);
            }

            virtual double v(const double & s) const
            {
                return ff_KMPW(s, _f0_V(), _b1_V(), _m_Bs2_1m);
            }

            virtual double a_0(const double & s) const
            {
                return ff_KMPW(s, _f0_A0(), _b1_A0(), _m_Bs2_0m);
            }

            virtual double a_1(const double & s) const
            {
                return ff_KMPW(s, _f0_A1(), _b1_A1(), _m_Bs2_1p);
            }

            virtual double a_2(const double & s) const
            {
                return ff_KMPW(s, _f0_A2(), _b1_A2(), _m_Bs2_1p);
            }

            virtual double a_12(const double & s) const
            {
                const double mB = BToKstar::mB, mB2 = mB * mB;
                const double mV = BToKstar::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB + mV) * (mB + mV) * (mB2 - mV2 - s) * this->a_1(s)
                    - lambda * this->a_2(s)) / (16.0 * mB * mV2 * (mB + mV));
            }

            virtual double t_1(const double & s) const
            {
                return ff_KMPW(s, _f0_T1(), _b1_T1(), _m_Bs2_1m);
            }

            virtual double t_2(const double & s) const
            {
                return ff_KMPW(s, _f0_T2(), _b1_T2(), _m_Bs2_1p);
            }

            virtual double t_3(const double & s) const
            {
                return ff_KMPW(s, _f0_T3(), _b1_T3(), _m_Bs2_1p);
            }

            virtual double t_23(const double & s) const
            {
                const double mB = BToKstar::mB, mB2 = mB * mB;
                const double mV = BToKstar::mV, mV2 = mV * mV;
                const double lambda = eos::lambda(mB2, mV2, s);

                return ((mB2 - mV2) * (mB2 + 3.0 * mV2 - s) * this->t_2(s)
                        - lambda * this->t_3(s)) / (8.0 * mB * mV2 * (mB - mV));
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
    };

    template <> class KMPW2010FormFactors<PToP> :
        public FormFactors<PToP>
    {
        private:
            // fit parametrisation for P -> P according to [KMPW2010]
            UsedParameter _b1_p, _b1_0, _b1_t;
            UsedParameter _f0_p, _f0_t;
            static const double _tau_p, _tau_m, _tau_0;
            static const double _m_B, _m_K, _m_Bs2;

            static double _calc_z(const double & s)
            {
                return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
            }

        public:
            KMPW2010FormFactors(const Parameters & p, const Options &) :
                _b1_p(p["B->K::b^p_1@KMPW2010"], *this),
                _b1_0(p["B->K::b^0_1@KMPW2010"], *this),
                _b1_t(p["B->K::b^t_1@KMPW2010"], *this),
                _f0_p(p["B->K::F^p(0)@KMPW2010"], *this),
                _f0_t(p["B->K::F^t(0)@KMPW2010"], *this)
            {
            }

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options)
            {
                return new KMPW2010FormFactors(parameters, options);
            }

            virtual double f_p(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_p() / (1 - s / _m_Bs2) * (1 + _b1_p() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double f_0(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_p() * (1 + _b1_0() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double f_t(const double & s) const
            {
                // cf. [KMPW2010], Eq. (8.8), p. 30
                const double zs = _calc_z(s), z0 = _calc_z(0.0);

                return _f0_t() / (1 - s / _m_Bs2) * (1 + _b1_t() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
            }

            virtual double f_plus_T(const double &) const
            {
                return 0.0;
            }
    };

    /* P -> PP Processes */

    struct BToPiPi {
        using Transition = PToPP;
        static constexpr const char * label = "B->pipi";
        static constexpr double mB  = 5.2795;
        static constexpr double mP1 = 0.13957;
        static constexpr double mP2 = 0.13957;

        // for pole and t_0 calculation in zhat
        static constexpr double mBst = 5.32465;

        // for pole calculation in z, depending on the current at hand
        static constexpr double mR2_1m = 5.32465;
        static constexpr double mR2_1p = 5.72590;
        static constexpr double mR2_0m = 5.27932;
    };

    template <typename Process_> class FvDV2018FormFactors :
        public FormFactors<PToPP>
    {
        private:
            UsedParameter _a_Fperp_0_0, _a_Fperp_0_1, _a_Fperp_0_2, _a_Fperp_0_3, _a_Fperp_1_0, _a_Fperp_1_1, _a_Fperp_1_2;
            UsedParameter _b_Fperp_0_0, _b_Fperp_0_1, _b_Fperp_0_2, _b_Fperp_0_3, _b_Fperp_1_0, _b_Fperp_1_1, _b_Fperp_1_2;
            UsedParameter _c_Fperp_0_0, _c_Fperp_0_1, _c_Fperp_0_2, _c_Fperp_0_3, _c_Fperp_1_0, _c_Fperp_1_1, _c_Fperp_1_2;

            UsedParameter _a_Fpara_0_0, _a_Fpara_0_1, _a_Fpara_0_2, _a_Fpara_0_3, _a_Fpara_1_0, _a_Fpara_1_1, _a_Fpara_1_2;
            UsedParameter _b_Fpara_0_0, _b_Fpara_0_1, _b_Fpara_0_2, _b_Fpara_0_3, _b_Fpara_1_0, _b_Fpara_1_1, _b_Fpara_1_2;
            UsedParameter _c_Fpara_0_0, _c_Fpara_0_1, _c_Fpara_0_2, _c_Fpara_0_3, _c_Fpara_1_0, _c_Fpara_1_1, _c_Fpara_1_2;

            UsedParameter _a_Flong_0_0, _a_Flong_0_1, _a_Flong_0_2, _a_Flong_0_3, _a_Flong_1_0, _a_Flong_1_1, _a_Flong_1_2;
            UsedParameter _b_Flong_0_0, _b_Flong_0_1, _b_Flong_0_2, _b_Flong_0_3, _b_Flong_1_0, _b_Flong_1_1, _b_Flong_1_2;
            UsedParameter _c_Flong_0_0, _c_Flong_0_1, _c_Flong_0_2, _c_Flong_0_3, _c_Flong_1_0, _c_Flong_1_1, _c_Flong_1_2;

            UsedParameter _a_Ftime_0_0, _a_Ftime_0_1, _a_Ftime_0_2, _a_Ftime_0_3, _a_Ftime_1_0, _a_Ftime_1_1, _a_Ftime_1_2;
            UsedParameter _b_Ftime_0_0, _b_Ftime_0_1, _b_Ftime_0_2, _b_Ftime_0_3, _b_Ftime_1_0, _b_Ftime_1_1, _b_Ftime_1_2;
            UsedParameter _c_Ftime_0_0, _c_Ftime_0_1, _c_Ftime_0_2, _c_Ftime_0_3, _c_Ftime_1_0, _c_Ftime_1_1, _c_Ftime_1_2;

            static double _calc_z(const double & t, const double & t_p, const double & t_0)
            {
                return (std::sqrt(t_p - t) - std::sqrt(t_p - t_0)) / (std::sqrt(t_p - t) + std::sqrt(t_p - t_0));
            }

            double _z(const double & t) const
            {
                static constexpr double mB  = Process_::mB;
                static constexpr double mP1 = Process_::mP1;
                static constexpr double mP2 = Process_::mP2;

                static constexpr double t_p = power_of<2>(mB + mP1 + mP2);
                static constexpr double t_0 = 0.0;

                return _calc_z(t, t_p, t_0);
            }

            double _zhat(const double & that) const
            {
                static constexpr double mB    = Process_::mB;
                static constexpr double mP2   = Process_::mP2;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                static constexpr double that_p = power_of<2>(mB + mP2);
                static const     double that_0 = that_p - std::sqrt(that_p * (that_p - mBst2));

                return _calc_z(that, that_p, that_0);
            }

            double _blaschke(const double & z, const double & zh) const
            {
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double zBst2  = _z(mBst2);
                const double zhBst2 = _zhat(mBst2);

                double result  = 1.0;
                result *= (1.0 - z  *  zBst2) / (z  -  zBst2);
                result *= (1.0 - zh * zhBst2) / (zh - zhBst2);

                return result;
            }

            double _blaschke_res_qhat2(const double & z) const
            {
                static constexpr double mB    = Process_::mB;
                static constexpr double mP2   = Process_::mP2;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                static constexpr double that_p = power_of<2>(mB + mP2);

                const double zBst2  = _z(mBst2);

                double result  = 4.0 * (mBst2 - that_p);
                result *= (1.0 - z  *  zBst2) / (z  -  zBst2);

                return result;
            }

          public:
            FvDV2018FormFactors(const Parameters & p, const Options &) :
                // perp
                _a_Fperp_0_0(p["B->pipi::a^Fperp_0_0@FvDV2018"], *this),
                _a_Fperp_0_1(p["B->pipi::a^Fperp_0_1@FvDV2018"], *this),
                _a_Fperp_0_2(p["B->pipi::a^Fperp_0_2@FvDV2018"], *this),
                _a_Fperp_0_3(p["B->pipi::a^Fperp_0_3@FvDV2018"], *this),
                _a_Fperp_1_0(p["B->pipi::a^Fperp_1_0@FvDV2018"], *this),
                _a_Fperp_1_1(p["B->pipi::a^Fperp_1_1@FvDV2018"], *this),
                _a_Fperp_1_2(p["B->pipi::a^Fperp_1_2@FvDV2018"], *this),
                _b_Fperp_0_0(p["B->pipi::b^Fperp_0_0@FvDV2018"], *this),
                _b_Fperp_0_1(p["B->pipi::b^Fperp_0_1@FvDV2018"], *this),
                _b_Fperp_0_2(p["B->pipi::b^Fperp_0_2@FvDV2018"], *this),
                _b_Fperp_0_3(p["B->pipi::b^Fperp_0_3@FvDV2018"], *this),
                _b_Fperp_1_0(p["B->pipi::b^Fperp_1_0@FvDV2018"], *this),
                _b_Fperp_1_1(p["B->pipi::b^Fperp_1_1@FvDV2018"], *this),
                _b_Fperp_1_2(p["B->pipi::b^Fperp_1_2@FvDV2018"], *this),
                _c_Fperp_0_0(p["B->pipi::c^Fperp_0_0@FvDV2018"], *this),
                _c_Fperp_0_1(p["B->pipi::c^Fperp_0_1@FvDV2018"], *this),
                _c_Fperp_0_2(p["B->pipi::c^Fperp_0_2@FvDV2018"], *this),
                _c_Fperp_0_3(p["B->pipi::c^Fperp_0_3@FvDV2018"], *this),
                _c_Fperp_1_0(p["B->pipi::c^Fperp_1_0@FvDV2018"], *this),
                _c_Fperp_1_1(p["B->pipi::c^Fperp_1_1@FvDV2018"], *this),
                _c_Fperp_1_2(p["B->pipi::c^Fperp_1_2@FvDV2018"], *this),
                // para
                _a_Fpara_0_0(p["B->pipi::a^Fpara_0_0@FvDV2018"], *this),
                _a_Fpara_0_1(p["B->pipi::a^Fpara_0_1@FvDV2018"], *this),
                _a_Fpara_0_2(p["B->pipi::a^Fpara_0_2@FvDV2018"], *this),
                _a_Fpara_0_3(p["B->pipi::a^Fpara_0_3@FvDV2018"], *this),
                _a_Fpara_1_0(p["B->pipi::a^Fpara_1_0@FvDV2018"], *this),
                _a_Fpara_1_1(p["B->pipi::a^Fpara_1_1@FvDV2018"], *this),
                _a_Fpara_1_2(p["B->pipi::a^Fpara_1_2@FvDV2018"], *this),
                _b_Fpara_0_0(p["B->pipi::b^Fpara_0_0@FvDV2018"], *this),
                _b_Fpara_0_1(p["B->pipi::b^Fpara_0_1@FvDV2018"], *this),
                _b_Fpara_0_2(p["B->pipi::b^Fpara_0_2@FvDV2018"], *this),
                _b_Fpara_0_3(p["B->pipi::b^Fpara_0_3@FvDV2018"], *this),
                _b_Fpara_1_0(p["B->pipi::b^Fpara_1_0@FvDV2018"], *this),
                _b_Fpara_1_1(p["B->pipi::b^Fpara_1_1@FvDV2018"], *this),
                _b_Fpara_1_2(p["B->pipi::b^Fpara_1_2@FvDV2018"], *this),
                _c_Fpara_0_0(p["B->pipi::c^Fpara_0_0@FvDV2018"], *this),
                _c_Fpara_0_1(p["B->pipi::c^Fpara_0_1@FvDV2018"], *this),
                _c_Fpara_0_2(p["B->pipi::c^Fpara_0_2@FvDV2018"], *this),
                _c_Fpara_0_3(p["B->pipi::c^Fpara_0_3@FvDV2018"], *this),
                _c_Fpara_1_0(p["B->pipi::c^Fpara_1_0@FvDV2018"], *this),
                _c_Fpara_1_1(p["B->pipi::c^Fpara_1_1@FvDV2018"], *this),
                _c_Fpara_1_2(p["B->pipi::c^Fpara_1_2@FvDV2018"], *this),
                // long
                _a_Flong_0_0(p["B->pipi::a^Flong_0_0@FvDV2018"], *this),
                _a_Flong_0_1(p["B->pipi::a^Flong_0_1@FvDV2018"], *this),
                _a_Flong_0_2(p["B->pipi::a^Flong_0_2@FvDV2018"], *this),
                _a_Flong_0_3(p["B->pipi::a^Flong_0_3@FvDV2018"], *this),
                _a_Flong_1_0(p["B->pipi::a^Flong_1_0@FvDV2018"], *this),
                _a_Flong_1_1(p["B->pipi::a^Flong_1_1@FvDV2018"], *this),
                _a_Flong_1_2(p["B->pipi::a^Flong_1_2@FvDV2018"], *this),
                _b_Flong_0_0(p["B->pipi::b^Flong_0_0@FvDV2018"], *this),
                _b_Flong_0_1(p["B->pipi::b^Flong_0_1@FvDV2018"], *this),
                _b_Flong_0_2(p["B->pipi::b^Flong_0_2@FvDV2018"], *this),
                _b_Flong_0_3(p["B->pipi::b^Flong_0_3@FvDV2018"], *this),
                _b_Flong_1_0(p["B->pipi::b^Flong_1_0@FvDV2018"], *this),
                _b_Flong_1_1(p["B->pipi::b^Flong_1_1@FvDV2018"], *this),
                _b_Flong_1_2(p["B->pipi::b^Flong_1_2@FvDV2018"], *this),
                _c_Flong_0_0(p["B->pipi::c^Flong_0_0@FvDV2018"], *this),
                _c_Flong_0_1(p["B->pipi::c^Flong_0_1@FvDV2018"], *this),
                _c_Flong_0_2(p["B->pipi::c^Flong_0_2@FvDV2018"], *this),
                _c_Flong_0_3(p["B->pipi::c^Flong_0_3@FvDV2018"], *this),
                _c_Flong_1_0(p["B->pipi::c^Flong_1_0@FvDV2018"], *this),
                _c_Flong_1_1(p["B->pipi::c^Flong_1_1@FvDV2018"], *this),
                _c_Flong_1_2(p["B->pipi::c^Flong_1_2@FvDV2018"], *this),
                // time
                _a_Ftime_0_0(p["B->pipi::a^Ftime_0_0@FvDV2018"], *this),
                _a_Ftime_0_1(p["B->pipi::a^Ftime_0_1@FvDV2018"], *this),
                _a_Ftime_0_2(p["B->pipi::a^Ftime_0_2@FvDV2018"], *this),
                _a_Ftime_0_3(p["B->pipi::a^Ftime_0_3@FvDV2018"], *this),
                _a_Ftime_1_0(p["B->pipi::a^Ftime_1_0@FvDV2018"], *this),
                _a_Ftime_1_1(p["B->pipi::a^Ftime_1_1@FvDV2018"], *this),
                _a_Ftime_1_2(p["B->pipi::a^Ftime_1_2@FvDV2018"], *this),
                _b_Ftime_0_0(p["B->pipi::b^Ftime_0_0@FvDV2018"], *this),
                _b_Ftime_0_1(p["B->pipi::b^Ftime_0_1@FvDV2018"], *this),
                _b_Ftime_0_2(p["B->pipi::b^Ftime_0_2@FvDV2018"], *this),
                _b_Ftime_0_3(p["B->pipi::b^Ftime_0_3@FvDV2018"], *this),
                _b_Ftime_1_0(p["B->pipi::b^Ftime_1_0@FvDV2018"], *this),
                _b_Ftime_1_1(p["B->pipi::b^Ftime_1_1@FvDV2018"], *this),
                _b_Ftime_1_2(p["B->pipi::b^Ftime_1_2@FvDV2018"], *this),
                _c_Ftime_0_0(p["B->pipi::c^Ftime_0_0@FvDV2018"], *this),
                _c_Ftime_0_1(p["B->pipi::c^Ftime_0_1@FvDV2018"], *this),
                _c_Ftime_0_2(p["B->pipi::c^Ftime_0_2@FvDV2018"], *this),
                _c_Ftime_0_3(p["B->pipi::c^Ftime_0_3@FvDV2018"], *this),
                _c_Ftime_1_0(p["B->pipi::c^Ftime_1_0@FvDV2018"], *this),
                _c_Ftime_1_1(p["B->pipi::c^Ftime_1_1@FvDV2018"], *this),
                _c_Ftime_1_2(p["B->pipi::c^Ftime_1_2@FvDV2018"], *this)
            {
            }

            static FormFactors<PToPP> * make(const Parameters & parameters, const Options & options)
            {
                return new FvDV2018FormFactors(parameters, options);
            }

            virtual complex<double> f_perp(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2);

                const double a = _a_Fperp_0_0 + _a_Fperp_1_0 * z + _a_Fperp_0_1 * zh + _a_Fperp_1_1 * z * zh + _a_Fperp_1_2 * z * zh * zh + _a_Fperp_0_2 * zh * zh + _a_Fperp_0_3 * zh * zh * zh;
                const double b = _b_Fperp_0_0 + _b_Fperp_1_0 * z + _b_Fperp_0_1 * zh + _b_Fperp_1_1 * z * zh + _b_Fperp_1_2 * z * zh * zh + _b_Fperp_0_2 * zh * zh + _b_Fperp_0_3 * zh * zh * zh;
                const double c = _c_Fperp_0_0 + _c_Fperp_1_0 * z + _c_Fperp_0_1 * zh + _c_Fperp_1_1 * z * zh + _c_Fperp_1_2 * z * zh * zh + _c_Fperp_0_2 * zh * zh + _c_Fperp_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * std::sqrt(lambda) / (mB * std::sqrt(k2)) };

                return result;
            }

            virtual double f_perp_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double lambda = eos::lambda(q2, k2, mB2);
                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Fperp_0_0 + _a_Fperp_1_0 * z + _a_Fperp_0_1 * zh + _a_Fperp_1_1 * z * zh + _a_Fperp_1_2 * z * zh * zh + _a_Fperp_0_2 * zh * zh + _a_Fperp_0_3 * zh * zh * zh;
                const double b = _b_Fperp_0_0 + _b_Fperp_1_0 * z + _b_Fperp_0_1 * zh + _b_Fperp_1_1 * z * zh + _b_Fperp_1_2 * z * zh * zh + _b_Fperp_0_2 * zh * zh + _b_Fperp_0_3 * zh * zh * zh;
                const double c = _c_Fperp_0_0 + _c_Fperp_1_0 * z + _c_Fperp_0_1 * zh + _c_Fperp_1_1 * z * zh + _c_Fperp_1_2 * z * zh * zh + _c_Fperp_0_2 * zh * zh + _c_Fperp_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * std::sqrt(lambda) / (mB * std::sqrt(k2));

                return result;
            }

            virtual complex<double> f_para(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2);

                const double a = _a_Fpara_0_0 + _a_Fpara_1_0 * z + _a_Fpara_0_1 * zh + _a_Fpara_1_1 * z * zh + _a_Fpara_1_2 * z * zh * zh + _a_Fpara_0_2 * zh * zh + _a_Fpara_0_3 * zh * zh * zh;
                const double b = _b_Fpara_0_0 + _b_Fpara_1_0 * z + _b_Fpara_0_1 * zh + _b_Fpara_1_1 * z * zh + _b_Fpara_1_2 * z * zh * zh + _b_Fpara_0_2 * zh * zh + _b_Fpara_0_3 * zh * zh * zh;
                const double c = _c_Fpara_0_0 + _c_Fpara_1_0 * z + _c_Fpara_0_1 * zh + _c_Fpara_1_1 * z * zh + _c_Fpara_1_2 * z * zh * zh + _c_Fpara_0_2 * zh * zh + _c_Fpara_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(k2) };

                return result;
            }

            virtual double f_para_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Fpara_0_0 + _a_Fpara_1_0 * z + _a_Fpara_0_1 * zh + _a_Fpara_1_1 * z * zh + _a_Fpara_1_2 * z * zh * zh + _a_Fpara_0_2 * zh * zh + _a_Fpara_0_3 * zh * zh * zh;
                const double b = _b_Fpara_0_0 + _b_Fpara_1_0 * z + _b_Fpara_0_1 * zh + _b_Fpara_1_1 * z * zh + _b_Fpara_1_2 * z * zh * zh + _b_Fpara_0_2 * zh * zh + _b_Fpara_0_3 * zh * zh * zh;
                const double c = _c_Fpara_0_0 + _c_Fpara_1_0 * z + _c_Fpara_0_1 * zh + _c_Fpara_1_1 * z * zh + _c_Fpara_1_2 * z * zh * zh + _c_Fpara_0_2 * zh * zh + _c_Fpara_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(k2);

                return result;
            }

            virtual complex<double> f_long(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2);

                const double a = _a_Flong_0_0 + _a_Flong_1_0 * z + _a_Flong_0_1 * zh + _a_Flong_1_1 * z * zh + _a_Flong_1_2 * z * zh * zh + _a_Flong_0_2 * zh * zh + _a_Flong_0_3 * zh * zh * zh;
                const double b = _b_Flong_0_0 + _b_Flong_1_0 * z + _b_Flong_0_1 * zh + _b_Flong_1_1 * z * zh + _b_Flong_1_2 * z * zh * zh + _b_Flong_0_2 * zh * zh + _b_Flong_0_3 * zh * zh * zh;
                const double c = _c_Flong_0_0 + _c_Flong_1_0 * z + _c_Flong_0_1 * zh + _c_Flong_1_1 * z * zh + _c_Flong_1_2 * z * zh * zh + _c_Flong_0_2 * zh * zh + _c_Flong_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(q2) * mB2 / std::sqrt(lambda) * mB2 / k2 };

                return result;
            }

            virtual double f_long_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double lambda = eos::lambda(q2, k2, mB2);
                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Flong_0_0 + _a_Flong_1_0 * z + _a_Flong_0_1 * zh + _a_Flong_1_1 * z * zh + _a_Flong_1_2 * z * zh * zh + _a_Flong_0_2 * zh * zh + _a_Flong_0_3 * zh * zh * zh;
                const double b = _b_Flong_0_0 + _b_Flong_1_0 * z + _b_Flong_0_1 * zh + _b_Flong_1_1 * z * zh + _b_Flong_1_2 * z * zh * zh + _b_Flong_0_2 * zh * zh + _b_Flong_0_3 * zh * zh * zh;
                const double c = _c_Flong_0_0 + _c_Flong_1_0 * z + _c_Flong_0_1 * zh + _c_Flong_1_1 * z * zh + _c_Flong_1_2 * z * zh * zh + _c_Flong_0_2 * zh * zh + _c_Flong_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB / std::sqrt(q2) * mB2 / std::sqrt(lambda) * mB2 / k2;

                return result;
            }

            virtual complex<double> f_time(const double & q2, const double & k2, const double & ctheta) const
            {
                static constexpr double mB  = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mP2 = Process_::mP2, mP22 = mP2 * mP2;

                const double lambda = eos::lambda(q2, k2, mB2);
                const double E2     = (mB2 + k2 - q2 - ctheta * std::sqrt(lambda)) / (4.0 * mB);
                const double qhat2  = mB2 + mP22 - 2.0 * mB * E2;

                const double z  = this->_z(q2);
                const double zh = this->_zhat(qhat2);

                const double a = _a_Ftime_0_0 + _a_Ftime_1_0 * z + _a_Ftime_0_1 * zh + _a_Ftime_1_1 * z * zh + _a_Ftime_1_2 * z * zh * zh + _a_Ftime_0_2 * zh * zh + _a_Ftime_0_3 * zh * zh * zh;
                const double b = _b_Ftime_0_0 + _b_Ftime_1_0 * z + _b_Ftime_0_1 * zh + _b_Ftime_1_1 * z * zh + _b_Ftime_1_2 * z * zh * zh + _b_Ftime_0_2 * zh * zh + _b_Ftime_0_3 * zh * zh * zh;
                const double c = _c_Ftime_0_0 + _c_Ftime_1_0 * z + _c_Ftime_0_1 * zh + _c_Ftime_1_1 * z * zh + _c_Ftime_1_2 * z * zh * zh + _c_Ftime_0_2 * zh * zh + _c_Ftime_0_3 * zh * zh * zh;

                const double blaschke = this->_blaschke(z, zh);

                complex<double> result{ 0.0, blaschke * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB * mB2 / std::sqrt(q2) / k2 };

                return result;
            }

            virtual double f_time_im_res_qhat2(const double & q2, const double & k2) const
            {
                static constexpr double mB    = Process_::mB,  mB2  = mB  * mB;
                static constexpr double mBst2 = power_of<2>(Process_::mBst);

                const double z  = this->_z(q2);
                const double zh = this->_z(mBst2);

                const double a = _a_Ftime_0_0 + _a_Ftime_1_0 * z + _a_Ftime_0_1 * zh + _a_Ftime_1_1 * z * zh + _a_Ftime_1_2 * z * zh * zh + _a_Ftime_0_2 * zh * zh + _a_Ftime_0_3 * zh * zh * zh;
                const double b = _b_Ftime_0_0 + _b_Ftime_1_0 * z + _b_Ftime_0_1 * zh + _b_Ftime_1_1 * z * zh + _b_Ftime_1_2 * z * zh * zh + _b_Ftime_0_2 * zh * zh + _b_Ftime_0_3 * zh * zh * zh;
                const double c = _c_Ftime_0_0 + _c_Ftime_1_0 * z + _c_Ftime_0_1 * zh + _c_Ftime_1_1 * z * zh + _c_Ftime_1_2 * z * zh * zh + _c_Ftime_0_2 * zh * zh + _c_Ftime_0_3 * zh * zh * zh;

                const double blaschke_res_qhat2 = this->_blaschke_res_qhat2(z);

                double result = blaschke_res_qhat2 * (a + b * (mB2 - k2) / mB2 + c * pow((mB2 - k2) / mB2, 2)) * mB * mB2 / std::sqrt(q2) / k2;

                return result;
            }
    };
}

#endif
