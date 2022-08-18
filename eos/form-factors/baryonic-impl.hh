/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2017 Elena Graverini
 * Copyright (c) 2017, 2018 Marzia Bordone
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_BARYONIC_IMPL_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/baryonic-processes.hh>
#include <eos/form-factors/hqet-b-to-c.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/utils/options.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    /*
     * J=1/2^+ -> J=1/2^+ transitions
     */

    /* Form Factors according to [MvD2016] for J=1/2^+ -> 1/2^+ transitions */
    template <typename Process_> class MvD2016FormFactors;

    template <typename Process_> class DKMR2017FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            UsedParameter _alpha_0_time_v, _alpha_1_time_v, _alpha_2_time_v;
            UsedParameter _alpha_0_time_a, _alpha_1_time_a, _alpha_2_time_a;

            UsedParameter _alpha_0_long_v, _alpha_1_long_v, _alpha_2_long_v;
            UsedParameter _alpha_0_long_a, _alpha_1_long_a, _alpha_2_long_a;
            UsedParameter _alpha_0_perp_v, _alpha_1_perp_v, _alpha_2_perp_v;
            UsedParameter                  _alpha_1_perp_a, _alpha_2_perp_a;

            UsedParameter _alpha_0_long_t,  _alpha_1_long_t,  _alpha_2_long_t;
            UsedParameter _alpha_0_long_t5, _alpha_1_long_t5, _alpha_2_long_t5;
            UsedParameter _alpha_0_perp_t,  _alpha_1_perp_t,  _alpha_2_perp_t;
            UsedParameter                   _alpha_1_perp_t5, _alpha_2_perp_t5;

            static constexpr double _z(const double & t, const double & tp, const double & t0)
        {
            return (std::sqrt(tp - t) - std::sqrt(tp - t0)) / (std::sqrt(tp - t) + std::sqrt(tp - t0));
        }

        public:
            DKMR2017FormFactors(const Parameters & p, const Options & o) :
                // time, V
                _alpha_0_time_v(p[stringify(Process_::label) + "::a_0_time^V@DKMR2017"], *this),
                _alpha_1_time_v(p[stringify(Process_::label) + "::a_1_time^V@DKMR2017"], *this),
                _alpha_2_time_v(p[stringify(Process_::label) + "::a_2_time^V@DKMR2017"], *this),
                // time, A
                _alpha_0_time_a(p[stringify(Process_::label) + "::a_0_time^A@DKMR2017"], *this),
                _alpha_1_time_a(p[stringify(Process_::label) + "::a_1_time^A@DKMR2017"], *this),
                _alpha_2_time_a(p[stringify(Process_::label) + "::a_2_time^A@DKMR2017"], *this),

                // long, V
                _alpha_0_long_v(p[stringify(Process_::label) + "::a_0_long^V@DKMR2017"], *this),
                _alpha_1_long_v(p[stringify(Process_::label) + "::a_1_long^V@DKMR2017"], *this),
                _alpha_2_long_v(p[stringify(Process_::label) + "::a_2_long^V@DKMR2017"], *this),
                // long, A
                _alpha_0_long_a(p[stringify(Process_::label) + "::a_0_long^A@DKMR2017"], *this),
                _alpha_1_long_a(p[stringify(Process_::label) + "::a_1_long^A@DKMR2017"], *this),
                _alpha_2_long_a(p[stringify(Process_::label) + "::a_2_long^A@DKMR2017"], *this),
                // perp, V
                _alpha_0_perp_v(p[stringify(Process_::label) + "::a_0_perp^V@DKMR2017"], *this),
                _alpha_1_perp_v(p[stringify(Process_::label) + "::a_1_perp^V@DKMR2017"], *this),
                _alpha_2_perp_v(p[stringify(Process_::label) + "::a_2_perp^V@DKMR2017"], *this),
                // perp, A
                _alpha_1_perp_a(p[stringify(Process_::label) + "::a_1_perp^A@DKMR2017"], *this),
                _alpha_2_perp_a(p[stringify(Process_::label) + "::a_2_perp^A@DKMR2017"], *this),

                // long, T
                _alpha_0_long_t(p[stringify(Process_::label) + "::a_0_long^T@DKMR2017"], *this),
                _alpha_1_long_t(p[stringify(Process_::label) + "::a_1_long^T@DKMR2017"], *this),
                _alpha_2_long_t(p[stringify(Process_::label) + "::a_2_long^T@DKMR2017"], *this),
                // long, T5
                _alpha_0_long_t5(p[stringify(Process_::label) + "::a_0_long^T5@DKMR2017"], *this),
                _alpha_1_long_t5(p[stringify(Process_::label) + "::a_1_long^T5@DKMR2017"], *this),
                _alpha_2_long_t5(p[stringify(Process_::label) + "::a_2_long^T5@DKMR2017"], *this),
                // perp, T
                _alpha_0_perp_t(p[stringify(Process_::label) + "::a_0_perp^T@DKMR2017"], *this),
                _alpha_1_perp_t(p[stringify(Process_::label) + "::a_1_perp^T@DKMR2017"], *this),
                _alpha_2_perp_t(p[stringify(Process_::label) + "::a_2_perp^T@DKMR2017"], *this),
                // perp, T5
                _alpha_1_perp_t5(p[stringify(Process_::label) + "::a_1_perp^T5@DKMR2017"], *this),
                _alpha_2_perp_t5(p[stringify(Process_::label) + "::a_2_perp^T5@DKMR2017"], *this)
        {
        }

        static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, const Options & options)
        {
            return new DKMR2017FormFactors(parameters, options);
        }

        // vector current
        virtual double f_time_v(const double & s) const
        {
            static const double mR2 = Process_::mR2_0p;

            const double z = _z(s, Process_::tp_0p, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_v() + _alpha_1_time_v() * z + _alpha_2_time_v() * z2);
        }

        virtual double f_long_v(const double & s) const
        {
            static const double mR2 = Process_::mR2_1m;

            const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_v() + _alpha_1_long_v() * z + _alpha_2_long_v() * z2);
        }

        virtual double f_perp_v(const double & s) const
        {
            static const double mR2 = Process_::mR2_1m;

            const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_v() + _alpha_1_perp_v() * z + _alpha_2_perp_v() * z2);
        }

        // axial vector current
        virtual double f_time_a(const double & s) const
        {
            static const double mR2 = Process_::mR2_0m;

            const double z = _z(s, Process_::tp_0m, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_a() + _alpha_1_time_a() * z + _alpha_2_time_a() * z2);
        }

        virtual double f_long_a(const double & s) const
        {
            static const double mR2 = Process_::mR2_1p;

            const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_long_a() * z + _alpha_2_long_a() * z2);
        }

        virtual double f_perp_a(const double & s) const
        {
            static const double mR2 = Process_::mR2_1p;

            const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

            // Using alpha_0_long_a instead of alpha_0_perp_a, in order to
            // fulfill relation eq. (7), [DM2016], p. 3.
            return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_perp_a() * z + _alpha_2_perp_a() * z2);
        }

        // tensor current
        virtual double f_long_t(const double & s) const
        {
            static const double mR2 = Process_::mR2_1m;

            const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t() + _alpha_1_long_t() * z + _alpha_2_long_t() * z2);
        }

        virtual double f_perp_t(const double & s) const
        {
            static const double mR2 = Process_::mR2_1m;

            const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_t() + _alpha_1_perp_t() * z + _alpha_2_perp_t() * z2);
        }

        // axial tensor current
        virtual double f_long_t5(const double & s) const
        {
            static const double mR2 = Process_::mR2_1p;

            const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

            return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_long_t5() * z + _alpha_2_long_t5() * z2);
        }

        virtual double f_perp_t5(const double & s) const
        {
            static const double mR2 = Process_::mR2_1p;

            const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

            // Using alpha_0_long_t5 instead of alpha_0_perp_t5, in order to
            // fulfill relation eq. (8), [DM2016], p. 3.
            return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_perp_t5() * z + _alpha_2_perp_t5() * z2);
        }
    };

    template <typename Process_> class DM2016FormFactors :
        public FormFactors<OneHalfPlusToOneHalfPlus>
    {
        private:
            UsedParameter _alpha_0_time_v, _alpha_1_time_v, _alpha_2_time_v;
            UsedParameter _alpha_0_time_a, _alpha_1_time_a, _alpha_2_time_a;

            UsedParameter _alpha_0_long_v, _alpha_1_long_v, _alpha_2_long_v;
            UsedParameter _alpha_0_long_a, _alpha_1_long_a, _alpha_2_long_a;
            UsedParameter _alpha_0_perp_v, _alpha_1_perp_v, _alpha_2_perp_v;
            UsedParameter                  _alpha_1_perp_a, _alpha_2_perp_a;

            UsedParameter _alpha_0_long_t,  _alpha_1_long_t,  _alpha_2_long_t;
            UsedParameter _alpha_0_long_t5, _alpha_1_long_t5, _alpha_2_long_t5;
            UsedParameter _alpha_0_perp_t,  _alpha_1_perp_t,  _alpha_2_perp_t;
            UsedParameter                   _alpha_1_perp_t5, _alpha_2_perp_t5;

            static constexpr double _z(const double & t, const double & tp, const double & t0)
            {
                return (std::sqrt(tp - t) - std::sqrt(tp - t0)) / (std::sqrt(tp - t) + std::sqrt(tp - t0));
            }

        public:
            DM2016FormFactors(const Parameters & p, const Options &) :
                // time, V
                _alpha_0_time_v(p[stringify(Process_::label) + "::a_0_time^V@DM2016"], *this),
                _alpha_1_time_v(p[stringify(Process_::label) + "::a_1_time^V@DM2016"], *this),
                _alpha_2_time_v(p[stringify(Process_::label) + "::a_2_time^V@DM2016"], *this),
                // time, A
                _alpha_0_time_a(p[stringify(Process_::label) + "::a_0_time^A@DM2016"], *this),
                _alpha_1_time_a(p[stringify(Process_::label) + "::a_1_time^A@DM2016"], *this),
                _alpha_2_time_a(p[stringify(Process_::label) + "::a_2_time^A@DM2016"], *this),

                // long, V
                _alpha_0_long_v(p[stringify(Process_::label) + "::a_0_long^V@DM2016"], *this),
                _alpha_1_long_v(p[stringify(Process_::label) + "::a_1_long^V@DM2016"], *this),
                _alpha_2_long_v(p[stringify(Process_::label) + "::a_2_long^V@DM2016"], *this),
                // long, A
                _alpha_0_long_a(p[stringify(Process_::label) + "::a_0_long^A@DM2016"], *this),
                _alpha_1_long_a(p[stringify(Process_::label) + "::a_1_long^A@DM2016"], *this),
                _alpha_2_long_a(p[stringify(Process_::label) + "::a_2_long^A@DM2016"], *this),
                // perp, V
                _alpha_0_perp_v(p[stringify(Process_::label) + "::a_0_perp^V@DM2016"], *this),
                _alpha_1_perp_v(p[stringify(Process_::label) + "::a_1_perp^V@DM2016"], *this),
                _alpha_2_perp_v(p[stringify(Process_::label) + "::a_2_perp^V@DM2016"], *this),
                // perp, A
                _alpha_1_perp_a(p[stringify(Process_::label) + "::a_1_perp^A@DM2016"], *this),
                _alpha_2_perp_a(p[stringify(Process_::label) + "::a_2_perp^A@DM2016"], *this),

                // long, T
                _alpha_0_long_t(p[stringify(Process_::label) + "::a_0_long^T@DM2016"], *this),
                _alpha_1_long_t(p[stringify(Process_::label) + "::a_1_long^T@DM2016"], *this),
                _alpha_2_long_t(p[stringify(Process_::label) + "::a_2_long^T@DM2016"], *this),
                // long, T5
                _alpha_0_long_t5(p[stringify(Process_::label) + "::a_0_long^T5@DM2016"], *this),
                _alpha_1_long_t5(p[stringify(Process_::label) + "::a_1_long^T5@DM2016"], *this),
                _alpha_2_long_t5(p[stringify(Process_::label) + "::a_2_long^T5@DM2016"], *this),
                // perp, T
                _alpha_0_perp_t(p[stringify(Process_::label) + "::a_0_perp^T@DM2016"], *this),
                _alpha_1_perp_t(p[stringify(Process_::label) + "::a_1_perp^T@DM2016"], *this),
                _alpha_2_perp_t(p[stringify(Process_::label) + "::a_2_perp^T@DM2016"], *this),
                // perp, T5
                _alpha_1_perp_t5(p[stringify(Process_::label) + "::a_1_perp^T5@DM2016"], *this),
                _alpha_2_perp_t5(p[stringify(Process_::label) + "::a_2_perp^T5@DM2016"], *this)
            {
            }

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, const Options & options)
            {
                return new DM2016FormFactors(parameters, options);
            }

            // vector current
            virtual double f_time_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_0p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_v() + _alpha_1_time_v() * z + _alpha_2_time_v() * z2);
            }

            virtual double f_long_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_v() + _alpha_1_long_v() * z + _alpha_2_long_v() * z2);
            }

            virtual double f_perp_v(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_v() + _alpha_1_perp_v() * z + _alpha_2_perp_v() * z2);
            }

            // axial vector current
            virtual double f_time_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_0m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_a() + _alpha_1_time_a() * z + _alpha_2_time_a() * z2);
            }

            virtual double f_long_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_long_a() * z + _alpha_2_long_a() * z2);
            }

            virtual double f_perp_a(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                // Using alpha_0_long_a instead of alpha_0_perp_a, in order to
                // fulfill relation eq. (7), [DM2016], p. 3.
                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_perp_a() * z + _alpha_2_perp_a() * z2);
            }

            // tensor current
            virtual double f_long_t(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t() + _alpha_1_long_t() * z + _alpha_2_long_t() * z2);
            }

            virtual double f_perp_t(const double & s) const
            {
                static const double mR2 = Process_::mR2_1m;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_t() + _alpha_1_perp_t() * z + _alpha_2_perp_t() * z2);
            }

            // axial tensor current
            virtual double f_long_t5(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_long_t5() * z + _alpha_2_long_t5() * z2);
            }

            virtual double f_perp_t5(const double & s) const
            {
                static const double mR2 = Process_::mR2_1p;

                const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

                // Using alpha_0_long_t5 instead of alpha_0_perp_t5, in order to
                // fulfill relation eq. (8), [DM2016], p. 3.
                return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_perp_t5() * z + _alpha_2_perp_t5() * z2);
            }
    };

    template <typename Transition_, typename Process_> class HQETFormFactors;

    /*
     * J=1/2^+ -> J=1/2^- transitions
     */

    template <typename Process_> class HQETFormFactors<OneHalfPlusToOneHalfMinus, Process_> :
        public FormFactors<OneHalfPlusToOneHalfMinus>
    {
        private:
            HQETBToC _b_to_c;

            UsedParameter _zeta_max, _rho, _delta_3b, _rho_3b;

            static constexpr double mLb = Process_::m1;
            static constexpr double mLcs = Process_::m2;
            static constexpr double mLb2 = mLb * mLb;
            static constexpr double mLcs2 = mLcs * mLcs;

            static constexpr double m_b_pole = 4.8;
            static constexpr double m_c_pole = 1.4;

            static constexpr double lambdabar = mLb - m_b_pole;
            static constexpr double lambdabarprime = mLcs - m_c_pole;

            static constexpr double s_max = (mLb - mLcs) * (mLb - mLcs);

            // auxiliary kinematics functions
            static constexpr double _s_plus(const double & s)
            {
                return power_of<2>(mLb + mLcs) - s;
            }
            static constexpr double _s_minus(const double & s)
            {
                return power_of<2>(mLb - mLcs) - s;
            }

            // parametrization of the Isgur-Wise functions
            double _z(const double & s) const
            {
                return _zeta_max * (1.0 + _rho * (s / s_max - 1.0));
            }

            double _z3b(const double & s) const
            {
                return _zeta_max * (_delta_3b + _rho_3b * (s / s_max - 1.0));
            }

            inline double omega(const double & s) const
            {
                return (mLb2 + mLcs2 - s) / (2.0 * mLb * mLcs);
            }

            inline double omegabar(const double & s) const
            {
                return omega(s) * (1.0 + lambdabar / m_b_pole + lambdabarprime / m_c_pole)
                    - (lambdabar / m_c_pole + lambdabarprime / m_b_pole);
            }

        public:
            HQETFormFactors(const Parameters & p) :
                _b_to_c(p, Options{ }),
                _zeta_max(p["Lambda_b->Lambda_c^*::zeta(q^2_max)@HQET"], *this),
                _rho(p["Lambda_b->Lambda_c^*::rho@HQET"], *this),
                _delta_3b(p["Lambda_b->Lambda_c^*::delta_3b@HQET"], *this),
                _rho_3b(p["Lambda_b->Lambda_c^*::rho_3b@HQET"], *this)
            {
                uses(_b_to_c);
            }

            static FormFactors<OneHalfPlusToOneHalfMinus> * make(const Parameters & parameters, const Options &)
            {
                return new HQETFormFactors(parameters);
            }

            // vector current
            virtual double f_time_v(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_vector(omegabar);
                const double C_2 = _b_to_c.c_2_vector(omegabar);
                const double C_3 = _b_to_c.c_3_vector(omegabar);

                // leading-power IWF
                double result = C_1 * sp;
                result += (mLb + mLcs) / (mLb - mLcs) * (mLb2 - mLcs2 + s) / (2.0 * mLb ) * (lambdabar      + C_2 * sp / (mLb + mLcs));
                result -= (mLb + mLcs) / (mLb - mLcs) * (mLb2 - mLcs2 - s) / (2.0 * mLcs) * (lambdabarprime - C_3 * sp / (mLb + mLcs));
                result *= _z(s);

                // next-to-leading-power IWF
                result -= 2.0 * (mLb + mLcs) * (mLb + mLcs) / (mLb - mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sm / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_long_v(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_vector(omegabar);
                const double C_2 = _b_to_c.c_2_vector(omegabar);
                const double C_3 = _b_to_c.c_3_vector(omegabar);

                // leading-power IWF
                double result = C_1 + sp * (C_2 * mLcs + C_3 * mLb) / (2.0 * mLb * mLcs * (mLb + mLcs));
                result *= sm;
                result += (mLb - mLcs) / (mLb + mLcs) * ((mLb2 - mLcs2 + s) / (2.0 * mLb) * lambdabar - (mLb2 - mLcs2 - s) / (2.0 * mLcs) * lambdabarprime);
                result *= _z(s);

                // next-to-leading-power IWF
                result -= 2.0 * (mLb - mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sp / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_perp_v(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_vector(omegabar);

                // leading-power IWF
                double result = C_1 * sm + (3.0 * mLb2 + mLcs2 - s) / (2.0 * mLb) * lambdabar - (mLb2 + 3.0 * mLcs2 - s) / (2.0 * mLcs) * lambdabarprime;
                result *= _z(s);

                // next-to-leading-power IWF
                result -= 2.0 * mLb * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sp / power_of<3>(mLb * mLcs));

                return result;
            }

            // axial vector current
            virtual double f_time_a(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_axialvector(omegabar);
                const double C_2 = _b_to_c.c_2_axialvector(omegabar);
                const double C_3 = _b_to_c.c_3_axialvector(omegabar);

                // leading-power IWF
                double result = C_1 * sm;
                result += (mLb - mLcs) / (mLb + mLcs) * (mLb2 - mLcs2 + s) / (2.0 * mLb ) * (lambdabar      - C_2 * sm / (mLb - mLcs));
                result -= (mLb - mLcs) / (mLb + mLcs) * (mLb2 - mLcs2 - s) / (2.0 * mLcs) * (lambdabarprime + C_3 * sm / (mLb - mLcs));
                result *= _z(s);

                // next-to-leading-power IWF
                result -= 2.0 * (mLb - mLcs) * (mLb - mLcs) / (mLb + mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sp / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_long_a(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_axialvector(omegabar);
                const double C_2 = _b_to_c.c_2_axialvector(omegabar);
                const double C_3 = _b_to_c.c_3_axialvector(omegabar);

                // leading-power IWF
                double result = C_1 - sm * (C_2 * mLcs + C_3 * mLb) / (2.0 * mLb * mLcs * (mLb - mLcs));
                result *= sp;
                result += (mLb + mLcs) / (mLb - mLcs) * ((mLb2 - mLcs2 + s) / (2.0 * mLb) * lambdabar - (mLb2 - mLcs2 - s) / (2.0 * mLcs) * lambdabarprime);
                result *= _z(s);

                // next-to-leading-power IWF
                result -= 2.0 * (mLb + mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sm / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_perp_a(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_axialvector(omegabar);

                // leading-power IWF
                double result = C_1 * sp + (3.0 * mLb2 + mLcs2 - s) / (2.0 * mLb) * lambdabar - (mLb2 + 3.0 * mLcs2 - s) / (2.0 * mLcs) * lambdabarprime;
                result *= _z(s);

                // next-to-leading-power IWF
                result -= 2.0 * mLb * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sm / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual Diagnostics diagnostics() const
            {
                Diagnostics results;

                // s = s_max
                {
                    const double s        = s_max;
                    const double omega    = this->omega(s);
                    const double omegabar = this->omegabar(s);
                    const double C_1_v    = _b_to_c.c_1_vector(omegabar);
                    const double C_2_v    = _b_to_c.c_2_vector(omegabar);
                    const double C_3_v    = _b_to_c.c_3_vector(omegabar);
                    const double C_1_a    = _b_to_c.c_1_axialvector(omegabar);
                    const double C_2_a    = _b_to_c.c_2_axialvector(omegabar);
                    const double C_3_a    = _b_to_c.c_3_axialvector(omegabar);

                    results.add({ s,              "s = s_max"       });
                    results.add({ s - 9.16430310, "s - s_max"       });
                    results.add({ omega,          "omega(s_max)"    });
                    results.add({ omegabar,       "omegabar(s_max)" });
                    results.add({ C_1_v,          "C_1_v(s_max)"    });
                    results.add({ C_2_v,          "C_2_v(s_max)"    });
                    results.add({ C_3_v,          "C_3_v(s_max)"    });
                    results.add({ C_1_a,          "C_1_a(s_max)"    });
                    results.add({ C_2_a,          "C_2_a(s_max)"    });
                    results.add({ C_3_a,          "C_3_a(s_max)"    });
                    results.add({ lambdabar,      "LambdaBar"       });
                    results.add({ lambdabarprime, "LambdaBar'"      });
                    results.add({ f_time_v(s),    "f_{time}"        });
                    results.add({ f_long_v(s),    "f_{long}"        });
                    results.add({ f_perp_v(s),    "f_{perp}"        });
                    results.add({ f_time_a(s),    "g_{time}"        });
                    results.add({ f_long_a(s),    "g_{long}"        });
                    results.add({ f_perp_a(s),    "g_{perp}"        });
                }

                // s = s_max - 3.0
                {
                    const double s        = s_max - 3.0;
                    const double omega    = this->omega(s);
                    const double omegabar = this->omegabar(s);
                    const double C_1_v    = _b_to_c.c_1_vector(omegabar);
                    const double C_2_v    = _b_to_c.c_2_vector(omegabar);
                    const double C_3_v    = _b_to_c.c_3_vector(omegabar);
                    const double C_1_a    = _b_to_c.c_1_axialvector(omegabar);
                    const double C_2_a    = _b_to_c.c_2_axialvector(omegabar);
                    const double C_3_a    = _b_to_c.c_3_axialvector(omegabar);

                    results.add({ s,              "s = s_max"       });
                    results.add({ s - 9.16430310, "s - s_max"       });
                    results.add({ omega,          "omega(s_max)"    });
                    results.add({ omegabar,       "omegabar(s_max)" });
                    results.add({ C_1_v,          "C_1_v(s_max)"    });
                    results.add({ C_2_v,          "C_2_v(s_max)"    });
                    results.add({ C_3_v,          "C_3_v(s_max)"    });
                    results.add({ C_1_a,          "C_1_a(s_max)"    });
                    results.add({ C_2_a,          "C_2_a(s_max)"    });
                    results.add({ C_3_a,          "C_3_a(s_max)"    });
                    results.add({ lambdabar,      "LambdaBar"       });
                    results.add({ lambdabarprime , "LambdaBar'"      });
                    results.add({ f_time_v(s),    "f_{time}"        });
                    results.add({ f_long_v(s),    "f_{long}"        });
                    results.add({ f_perp_v(s),    "f_{perp}"        });
                    results.add({ f_time_a(s),    "g_{time}"        });
                    results.add({ f_long_a(s),    "g_{long}"        });
                    results.add({ f_perp_a(s),    "g_{perp}"        });
                }

                return results;
            }
    };

    /*
     * J=1/2^+ -> J=3/2^- transitions
     */

    template <typename Process_> class HQETFormFactors<OneHalfPlusToThreeHalfMinus, Process_> :
        public FormFactors<OneHalfPlusToThreeHalfMinus>
    {
        private:
            HQETBToC _b_to_c;

            UsedParameter _zeta_max, _rho, _delta_3b, _rho_3b;

            static constexpr double mLb = Process_::m1;
            static constexpr double mLcs = Process_::m2;
            static constexpr double mLb2 = mLb * mLb;
            static constexpr double mLcs2 = mLcs * mLcs;

            static constexpr double m_b_pole = 4.8;
            static constexpr double m_c_pole = 1.4;

            static constexpr double lambdabar = mLb - m_b_pole;
            static constexpr double lambdabarprime = mLcs - m_c_pole;

            static constexpr double s_max = (mLb - mLcs) * (mLb - mLcs);

            // auxiliary kinematics functions
            static constexpr double _s_plus(const double & s)
            {
                return power_of<2>((mLb + mLcs)) - s;
            }
            static constexpr double _s_minus(const double & s)
            {
                return power_of<2>((mLb - mLcs)) - s;
            }

            // parametrization of the Isgur-Wise functions
            double _z(const double & s) const
            {
                return _zeta_max * (1.0 + _rho * (s / s_max - 1.0));
            }

            double _z3b(const double & s) const
            {
                return _zeta_max * (_delta_3b + _rho_3b * (s / s_max - 1.0));
            }

            inline double omega(const double & s) const
            {
                return (mLb2 + mLcs2 - s) / (2.0 * mLb * mLcs);
            }

            inline double omegabar(const double & s) const
            {
                return omega(s) * (1.0 + lambdabar / m_b_pole + lambdabarprime / m_c_pole)
                    - (lambdabar / m_c_pole + lambdabarprime / m_b_pole);
            }

        public:
            HQETFormFactors(const Parameters & p) :
                _b_to_c(p, Options{ }),
                _zeta_max(p["Lambda_b->Lambda_c^*::zeta(q^2_max)@HQET"], *this),
                _rho(p["Lambda_b->Lambda_c^*::rho@HQET"], *this),
                _delta_3b(p["Lambda_b->Lambda_c^*::delta_3b@HQET"], *this),
                _rho_3b(p["Lambda_b->Lambda_c^*::rho_3b@HQET"], *this)
            {
                uses(_b_to_c);
            }

            static FormFactors<OneHalfPlusToThreeHalfMinus> * make(const Parameters & parameters, const Options &)
            {
                return new HQETFormFactors(parameters);
            }

            // vector current
            virtual double f_time12_v(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_vector(omegabar);
                const double C_2 = _b_to_c.c_2_vector(omegabar);
                const double C_3 = _b_to_c.c_3_vector(omegabar);

                // leading-power IWF
                double result = C_1 * sp;
                result += (mLb + mLcs) / (mLb - mLcs) * (mLb2 - mLcs2 + s) / (2.0 * mLb ) * (lambdabar      + C_2 * sp / (mLb + mLcs));
                result -= (mLb + mLcs) / (mLb - mLcs) * (mLb2 - mLcs2 - s) / (2.0 * mLcs) * (lambdabarprime - C_3 * sp / (mLb + mLcs));
                result *= _z(s);

                // next-to-leading-power IWF
                result += (mLb + mLcs) * (mLb + mLcs) / (mLb - mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sm / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_long12_v(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_vector(omegabar);
                const double C_2 = _b_to_c.c_2_vector(omegabar);
                const double C_3 = _b_to_c.c_3_vector(omegabar);

                // leading-power IWF
                double result = C_1 + sp * (C_2 * mLcs + C_3 * mLb) / (2.0 * mLb * mLcs * (mLb + mLcs));
                result *= sm;
                result += (mLb - mLcs) / (mLb + mLcs) * ((mLb2 - mLcs2 + s) / (2.0 * mLb) * lambdabar - (mLb2 - mLcs2 - s) / (2.0 * mLcs) * lambdabarprime);
                result *= _z(s);

                // next-to-leading-power IWF
                result += (mLb - mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sp / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_perp12_v(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_vector(omegabar);

                // leading-power IWF
                double result = C_1 * sm + (3.0 * mLb2 + mLcs2 - s) / (2.0 * mLb) * lambdabar - (mLb2 + 3.0 * mLcs2 - s) / (2.0 * mLcs) * lambdabarprime;
                result *= _z(s);

                // next-to-leading-power IWF
                result += mLb * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sp / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_perp32_v(const double & s) const
            {
                const double sp   = _s_plus(s);

                // next-to-leading-power IWF
                double result = -0.5 * sqrt(sp / (mLcs * power_of<3>(mLb))) * _z3b(s);

                return result;
            }

            // axial vector current
            virtual double f_time12_a(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_axialvector(omegabar);
                const double C_2 = _b_to_c.c_2_axialvector(omegabar);
                const double C_3 = _b_to_c.c_3_axialvector(omegabar);

                // leading-power IWF
                double result = C_1 * sm;
                result += (mLb - mLcs) / (mLb + mLcs) * (mLb2 - mLcs2 + s) / (2.0 * mLb ) * (lambdabar      - C_2 * sm / (mLb - mLcs));
                result -= (mLb - mLcs) / (mLb + mLcs) * (mLb2 - mLcs2 - s) / (2.0 * mLcs) * (lambdabarprime + C_3 * sm / (mLb - mLcs));
                result *= _z(s);

                // next-to-leading-power IWF
                result += (mLb - mLcs) * (mLb - mLcs) / (mLb + mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sp / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_long12_a(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_axialvector(omegabar);
                const double C_2 = _b_to_c.c_2_axialvector(omegabar);
                const double C_3 = _b_to_c.c_3_axialvector(omegabar);

                // leading-power IWF
                double result = C_1 - sm * (C_2 * mLcs + C_3 * mLb) / (2.0 * mLb * mLcs * (mLb - mLcs));
                result *= sp;
                result += (mLb + mLcs) / (mLb - mLcs) * ((mLb2 - mLcs2 + s) / (2.0 * mLb) * lambdabar - (mLb2 - mLcs2 - s) / (2.0 * mLcs) * lambdabarprime);
                result *= _z(s);

                // next-to-leading-power IWF
                result += (mLb + mLcs) * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sm / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_perp12_a(const double & s) const
            {
                const double omegabar = this->omegabar(s);
                const double sm       = _s_minus(s);
                const double sp       = _s_plus(s);

                const double C_1 = _b_to_c.c_1_axialvector(omegabar);

                // leading-power IWF
                double result = C_1 * sp + (3.0 * mLb2 + mLcs2 - s) / (2.0 * mLb) * lambdabar - (mLb2 + 3.0 * mLcs2 - s) / (2.0 * mLcs) * lambdabarprime;
                result *= _z(s);

                // next-to-leading-power IWF
                result += mLb * _z3b(s);

                // normalisation
                result *= 0.5 * sqrt(sm / power_of<3>(mLb * mLcs));

                return result;
            }

            virtual double f_perp32_a(const double & s) const
            {
                const double sm   = _s_minus(s);

                // next-to-leading-power IWF
                double result = -0.5 * sqrt(sm / (mLcs * power_of<3>(mLb))) * _z3b(s);

                return result;
            }

            // tensor current
            virtual double f_long12_t(const double &) const { throw InternalError("HQETFormFactors::f_long12_t(): not implemented"); }
            virtual double f_perp12_t(const double &) const { throw InternalError("HQETFormFactors::f_perp12_t(): not implemented"); }
            virtual double f_perp32_t(const double &) const { throw InternalError("HQETFormFactors::f_perp32_t(): not implemented"); }
            virtual double f_long12_t5(const double &) const { throw InternalError("HQETFormFactors::f_long12_t5(): not implemented"); }
            virtual double f_perp12_t5(const double &) const { throw InternalError("HQETFormFactors::f_perp12_t5(): not implemented"); }
            virtual double f_perp32_t5(const double &) const { throw InternalError("HQETFormFactors::f_perp32_t5(): not implemented"); }


            virtual Diagnostics diagnostics() const
            {
                Diagnostics results;

                // s = s_max
                {
                    const double s        = s_max;
                    const double omega    = this->omega(s);
                    const double omegabar = this->omegabar(s);
                    const double C_1_v    = _b_to_c.c_1_vector(omegabar);
                    const double C_2_v    = _b_to_c.c_2_vector(omegabar);
                    const double C_3_v    = _b_to_c.c_3_vector(omegabar);
                    const double C_1_a    = _b_to_c.c_1_axialvector(omegabar);
                    const double C_2_a    = _b_to_c.c_2_axialvector(omegabar);
                    const double C_3_a    = _b_to_c.c_3_axialvector(omegabar);

                    results.add({ s,              "s = s_max"       });
                    results.add({ s - 8.94847396, "s - s_max"       });
                    results.add({ omega,          "omega(s_max)"    });
                    results.add({ omegabar,       "omegabar(s_max)" });
                    results.add({ C_1_v,          "C_1_v(s_max)"    });
                    results.add({ C_2_v,          "C_2_v(s_max)"    });
                    results.add({ C_3_v,          "C_3_v(s_max)"    });
                    results.add({ C_1_a,          "C_1_a(s_max)"    });
                    results.add({ C_2_a,          "C_2_a(s_max)"    });
                    results.add({ C_3_a,          "C_3_a(s_max)"    });
                    results.add({ lambdabar,      "LambdaBar"       });
                    results.add({ lambdabarprime, "LambdaBar'"      });
                    results.add({ f_time12_v(s),  "F_{1/2,time}"    });
                    results.add({ f_long12_v(s),  "F_{1/2,long}"    });
                    results.add({ f_perp12_v(s),  "F_{1/2,perp}"    });
                    results.add({ f_perp32_v(s),  "F_{3/2,perp}"    });
                    results.add({ f_time12_a(s),  "G_{1/2,time}"    });
                    results.add({ f_long12_a(s),  "G_{1/2,long}"    });
                    results.add({ f_perp12_a(s),  "G_{1/2,perp}"    });
                    results.add({ f_perp32_a(s),  "G_{3/2,perp}"    });
                }

                // s = s_max - 3.0
                {
                    const double s        = s_max - 3.0;
                    const double omega    = this->omega(s);
                    const double omegabar = this->omegabar(s);
                    const double C_1_v    = _b_to_c.c_1_vector(omegabar);
                    const double C_2_v    = _b_to_c.c_2_vector(omegabar);
                    const double C_3_v    = _b_to_c.c_3_vector(omegabar);
                    const double C_1_a    = _b_to_c.c_1_axialvector(omegabar);
                    const double C_2_a    = _b_to_c.c_2_axialvector(omegabar);
                    const double C_3_a    = _b_to_c.c_3_axialvector(omegabar);

                    results.add({ s,              "s = s_max"       });
                    results.add({ s - 8.94847396, "s - s_max"       });
                    results.add({ omega,          "omega(s_max)"    });
                    results.add({ omegabar,       "omegabar(s_max)" });
                    results.add({ C_1_v,          "C_1_v(s_max)"    });
                    results.add({ C_2_v,          "C_2_v(s_max)"    });
                    results.add({ C_3_v,          "C_3_v(s_max)"    });
                    results.add({ C_1_a,          "C_1_a(s_max)"    });
                    results.add({ C_2_a,          "C_2_a(s_max)"    });
                    results.add({ C_3_a,          "C_3_a(s_max)"    });
                    results.add({ lambdabar,      "LambdaBar"       });
                    results.add({ lambdabarprime, "LambdaBar'"      });
                    results.add({ f_time12_v(s),  "F_{1/2,time}"    });
                    results.add({ f_long12_v(s),  "F_{1/2,long}"    });
                    results.add({ f_perp12_v(s),  "F_{1/2,perp}"    });
                    results.add({ f_perp32_v(s),  "F_{3/2,perp}"    });
                    results.add({ f_time12_a(s),  "G_{1/2,time}"    });
                    results.add({ f_long12_a(s),  "G_{1/2,long}"    });
                    results.add({ f_perp12_a(s),  "G_{1/2,perp}"    });
                    results.add({ f_perp32_a(s),  "G_{3/2,perp}"    });
                }

                return results;
            }
    };
}

#endif
