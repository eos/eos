/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015, 2016 Danny van Dyk
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
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>

namespace eos
{
    /* Form Factors according to [MvD2016] for J=1/2^+ -> 1/2^+ transitions */
    template <typename Process_> class MvD2016FormFactors;

    /*
     * J=1/2^+ -> J=1/2^+ transitions
     */
    struct LambdaBToLambda {
        static constexpr const char * label = "Lambda_b->Lambda";
        // initial state mass
        static constexpr double m1 = 5.61951;
        // final state mass
        static constexpr double m2 = 1.115683;
        // semileptonic kinematic endpoint
        static constexpr double tm = (m1 - m2) * (m1 - m2);
        // pair production threshold: B + K
        static constexpr double tp = (5.279 + 0.494) * (5.279 + 0.494);
        // first resonances sorted by spin/parity
        static constexpr double mR2_0m = 5.367 * 5.367;
        static constexpr double mR2_0p = 5.711 * 5.711;
        static constexpr double mR2_1m = 5.416 * 5.416;
        static constexpr double mR2_1p = 5.750 * 5.750;
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
            DM2016FormFactors(const Parameters & p) :
                // time, V
                _alpha_0_time_v(p["Lambda_b->Lambda::a_0_time^V@DM2016"], *this),
                _alpha_1_time_v(p["Lambda_b->Lambda::a_1_time^V@DM2016"], *this),
                _alpha_2_time_v(p["Lambda_b->Lambda::a_2_time^V@DM2016"], *this),
                // time, A
                _alpha_0_time_a(p["Lambda_b->Lambda::a_0_time^A@DM2016"], *this),
                _alpha_1_time_a(p["Lambda_b->Lambda::a_1_time^A@DM2016"], *this),
                _alpha_2_time_a(p["Lambda_b->Lambda::a_2_time^A@DM2016"], *this),

                // long, V
                _alpha_0_long_v(p["Lambda_b->Lambda::a_0_long^V@DM2016"], *this),
                _alpha_1_long_v(p["Lambda_b->Lambda::a_1_long^V@DM2016"], *this),
                _alpha_2_long_v(p["Lambda_b->Lambda::a_2_long^V@DM2016"], *this),
                // long, A
                _alpha_0_long_a(p["Lambda_b->Lambda::a_0_long^A@DM2016"], *this),
                _alpha_1_long_a(p["Lambda_b->Lambda::a_1_long^A@DM2016"], *this),
                _alpha_2_long_a(p["Lambda_b->Lambda::a_2_long^A@DM2016"], *this),
                // perp, V
                _alpha_0_perp_v(p["Lambda_b->Lambda::a_0_perp^V@DM2016"], *this),
                _alpha_1_perp_v(p["Lambda_b->Lambda::a_1_perp^V@DM2016"], *this),
                _alpha_2_perp_v(p["Lambda_b->Lambda::a_2_perp^V@DM2016"], *this),
                // perp, A
                _alpha_1_perp_a(p["Lambda_b->Lambda::a_1_perp^A@DM2016"], *this),
                _alpha_2_perp_a(p["Lambda_b->Lambda::a_2_perp^A@DM2016"], *this),

                // long, T
                _alpha_0_long_t(p["Lambda_b->Lambda::a_0_long^T@DM2016"], *this),
                _alpha_1_long_t(p["Lambda_b->Lambda::a_1_long^T@DM2016"], *this),
                _alpha_2_long_t(p["Lambda_b->Lambda::a_2_long^T@DM2016"], *this),
                // long, T5
                _alpha_0_long_t5(p["Lambda_b->Lambda::a_0_long^T5@DM2016"], *this),
                _alpha_1_long_t5(p["Lambda_b->Lambda::a_1_long^T5@DM2016"], *this),
                _alpha_2_long_t5(p["Lambda_b->Lambda::a_2_long^T5@DM2016"], *this),
                // perp, T
                _alpha_0_perp_t(p["Lambda_b->Lambda::a_0_perp^T@DM2016"], *this),
                _alpha_1_perp_t(p["Lambda_b->Lambda::a_1_perp^T@DM2016"], *this),
                _alpha_2_perp_t(p["Lambda_b->Lambda::a_2_perp^T@DM2016"], *this),
                // perp, T5
                _alpha_1_perp_t5(p["Lambda_b->Lambda::a_1_perp^T5@DM2016"], *this),
                _alpha_2_perp_t5(p["Lambda_b->Lambda::a_2_perp^T5@DM2016"], *this)
            {
            }

            static FormFactors<OneHalfPlusToOneHalfPlus> * make(const Parameters & parameters, unsigned)
            {
                return new DM2016FormFactors(parameters);
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
}

#endif
