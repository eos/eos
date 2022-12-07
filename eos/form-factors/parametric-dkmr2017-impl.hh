/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DKMR2017_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DKMR2017_IMPL_HH 1

#include <eos/form-factors/parametric-dkmr2017.hh>

namespace eos
{
    template <typename Process_>
    DKMR2017FormFactors<Process_>::DKMR2017FormFactors(const Parameters & p, const Options &) :
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

    template <typename Process_>
    FormFactors<OneHalfPlusToOneHalfPlus> *
    DKMR2017FormFactors<Process_>::make(const Parameters & parameters, const Options & options)
    {
        return new DKMR2017FormFactors(parameters, options);
    }

    // vector current
    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_time_v(const double & s) const
    {
        static const double mR2 = Process_::mR2_0p;

        const double z = _z(s, Process_::tp_0p, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_v() + _alpha_1_time_v() * z + _alpha_2_time_v() * z2);
    }

    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_long_v(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_v() + _alpha_1_long_v() * z + _alpha_2_long_v() * z2);
    }

    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_perp_v(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_v() + _alpha_1_perp_v() * z + _alpha_2_perp_v() * z2);
    }

    // axial vector current
    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_time_a(const double & s) const
    {
        static const double mR2 = Process_::mR2_0m;

        const double z = _z(s, Process_::tp_0m, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_a() + _alpha_1_time_a() * z + _alpha_2_time_a() * z2);
    }

    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_long_a(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_long_a() * z + _alpha_2_long_a() * z2);
    }

    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_perp_a(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

        // Using alpha_0_long_a instead of alpha_0_perp_a, in order to
        // fulfill relation eq. (7), [DM2016], p. 3.
        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_perp_a() * z + _alpha_2_perp_a() * z2);
    }

    // tensor current
    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_long_t(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t() + _alpha_1_long_t() * z + _alpha_2_long_t() * z2);
    }

    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_perp_t(const double & s) const
    {
        static const double mR2 = Process_::mR2_1m;

        const double z = _z(s, Process_::tp_1m, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_t() + _alpha_1_perp_t() * z + _alpha_2_perp_t() * z2);
    }

    // axial tensor current
    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_long_t5(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_long_t5() * z + _alpha_2_long_t5() * z2);
    }

    template <typename Process_>
    double
    DKMR2017FormFactors<Process_>::f_perp_t5(const double & s) const
    {
        static const double mR2 = Process_::mR2_1p;

        const double z = _z(s, Process_::tp_1p, Process_::tm), z2 = z * z;

        // Using alpha_0_long_t5 instead of alpha_0_perp_t5, in order to
        // fulfill relation eq. (8), [DM2016], p. 3.
        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_perp_t5() * z + _alpha_2_perp_t5() * z2);
    }

    template <typename Process_>
    const std::set<ReferenceName> DKMR2017FormFactors<Process_>::references
    {
        "DKMR:2017A"_rn
    };

    template <typename Process_>
    const std::vector<OptionSpecification> DKMR2017FormFactors<Process_>::options
    {
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    DKMR2017FormFactors<Process_>::begin_options()
    {
        return options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    DKMR2017FormFactors<Process_>::end_options()
    {
        return options.cend();
    }
}

#endif
