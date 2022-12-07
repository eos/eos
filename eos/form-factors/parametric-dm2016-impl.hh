/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014, 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2022 Méril Reboud
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DM2016_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DM2016_IMPL_HH 1

#include <eos/form-factors/parametric-dm2016.hh>

namespace eos
{
    template <typename Process_>
    DM2016FormFactors<Process_>::DM2016FormFactors(const Parameters & p, const Options &) :
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

    template <typename Process_>
    FormFactors<OneHalfPlusToOneHalfPlus> *
    DM2016FormFactors<Process_>::make(const Parameters & parameters, const Options & options)
    {
        return new DM2016FormFactors(parameters, options);
    }

    // vector current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_time_v(const double & s) const
    {
        static const double mR2 = Process_::mR2_0p;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_v() + _alpha_1_time_v() * z + _alpha_2_time_v() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_v(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_v() + _alpha_1_long_v() * z + _alpha_2_long_v() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_v(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_v() + _alpha_1_perp_v() * z + _alpha_2_perp_v() * z2);
    }

    // axial vector current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_time_a(const double & s) const
    {
        static const double mR2 = Process_::mR2_0m;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_a() + _alpha_1_time_a() * z + _alpha_2_time_a() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_a(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_long_a() * z + _alpha_2_long_a() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_a(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        // Using alpha_0_long_a instead of alpha_0_perp_a, in order to
        // fulfill relation eq. (7), [DM2016], p. 3.
        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_perp_a() * z + _alpha_2_perp_a() * z2);
    }

    // tensor current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_t(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t() + _alpha_1_long_t() * z + _alpha_2_long_t() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_t(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_t() + _alpha_1_perp_t() * z + _alpha_2_perp_t() * z2);
    }

    // axial tensor current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_t5(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_long_t5() * z + _alpha_2_long_t5() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_t5(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp, Process_::tm), z2 = z * z;

        // Using alpha_0_long_t5 instead of alpha_0_perp_t5, in order to
        // fulfill relation eq. (8), [DM2016], p. 3.
        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_perp_t5() * z + _alpha_2_perp_t5() * z2);
    }

    template <typename Process_>
    const std::set<ReferenceName> DM2016FormFactors<Process_>::references
    {
        "DM:2016A"_rn
    };

    template <typename Process_>
    const std::vector<OptionSpecification> DM2016FormFactors<Process_>::options
    {
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    DM2016FormFactors<Process_>::begin_options()
    {
        return options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    DM2016FormFactors<Process_>::end_options()
    {
        return options.cend();
    }
}

#endif
