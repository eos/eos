/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2022 Danny van Dyk
 * Copyright (c) 2021-2022 Muslem Rahimi
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BMRvD2022_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BMRvD2022_IMPL_HH 1

#include <eos/form-factors/parametric-bmrvd2022.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>

#include <numeric>

#include <iostream>

namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BMRvD2022FormFactorTraits<Process_>::resonance_0m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BMRvD2022FormFactorTraits<Process_>::resonance_0p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,0@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BMRvD2022FormFactorTraits<Process_>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BMRvD2022FormFactorTraits<Process_>::resonance_1p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,1@BSZ2015" }
    };

    template <typename Process_>
    BMRvD2022FormFactors<Process_>::BMRvD2022FormFactors(const Parameters & p, const Options &) :
        _a_time_v{
            // a^(time,V)_0 replaced by equation of motion
            UsedParameter(p[_par_name("t", "V", 1)], *this),
            UsedParameter(p[_par_name("t", "V", 2)], *this),
            UsedParameter(p[_par_name("t", "V", 3)], *this),
            UsedParameter(p[_par_name("t", "V", 4)], *this)
        },
        _a_long_v{
            UsedParameter(p[_par_name("0", "V", 0)], *this),
            UsedParameter(p[_par_name("0", "V", 1)], *this),
            UsedParameter(p[_par_name("0", "V", 2)], *this),
            UsedParameter(p[_par_name("0", "V", 3)], *this),
            UsedParameter(p[_par_name("0", "V", 4)], *this)
        },
        _a_perp_v{
            UsedParameter(p[_par_name("perp", "V", 0)], *this),
            UsedParameter(p[_par_name("perp", "V", 1)], *this),
            UsedParameter(p[_par_name("perp", "V", 2)], *this),
            UsedParameter(p[_par_name("perp", "V", 3)], *this),
            UsedParameter(p[_par_name("perp", "V", 4)], *this)
        },
        _a_time_a{
            // a^(time,A)_0 replaced by equation of motion
            UsedParameter(p[_par_name("t", "A", 1)], *this),
            UsedParameter(p[_par_name("t", "A", 2)], *this),
            UsedParameter(p[_par_name("t", "A", 3)], *this),
            UsedParameter(p[_par_name("t", "A", 4)], *this)
        },
        _a_long_a{
            UsedParameter(p[_par_name("0", "A", 0)], *this),
            UsedParameter(p[_par_name("0", "A", 1)], *this),
            UsedParameter(p[_par_name("0", "A", 2)], *this),
            UsedParameter(p[_par_name("0", "A", 3)], *this),
            UsedParameter(p[_par_name("0", "A", 4)], *this)
        },
        _a_perp_a{
            // a^(perp,A)_0 replaced by equation of motion
            UsedParameter(p[_par_name("perp", "A", 1)], *this),
            UsedParameter(p[_par_name("perp", "A", 2)], *this),
            UsedParameter(p[_par_name("perp", "A", 3)], *this),
            UsedParameter(p[_par_name("perp", "A", 4)], *this)
        },
        _a_long_t{
            UsedParameter(p[_par_name("0", "T", 0)], *this),
            UsedParameter(p[_par_name("0", "T", 1)], *this),
            UsedParameter(p[_par_name("0", "T", 2)], *this),
            UsedParameter(p[_par_name("0", "T", 3)], *this),
            UsedParameter(p[_par_name("0", "T", 4)], *this)
        },
        _a_perp_t{
            // a^(perp,T)_0 replaced by equation of motion
            UsedParameter(p[_par_name("perp", "T", 1)], *this),
            UsedParameter(p[_par_name("perp", "T", 2)], *this),
            UsedParameter(p[_par_name("perp", "T", 3)], *this),
            UsedParameter(p[_par_name("perp", "T", 4)], *this)
        },
        _a_long_t5{
            // a^(long,T5)_0 replaced by equation of motion
            UsedParameter(p[_par_name("0", "T5", 1)], *this),
            UsedParameter(p[_par_name("0", "T5", 2)], *this),
            UsedParameter(p[_par_name("0", "T5", 3)], *this),
            UsedParameter(p[_par_name("0", "T5", 4)], *this)
        },
        _a_perp_t5{
            UsedParameter(p[_par_name("perp", "T5", 0)], *this),
            UsedParameter(p[_par_name("perp", "T5", 1)], *this),
            UsedParameter(p[_par_name("perp", "T5", 2)], *this),
            UsedParameter(p[_par_name("perp", "T5", 3)], *this),
            UsedParameter(p[_par_name("perp", "T5", 4)], *this)
        },
        _traits(p),
        _m_1(_traits.m_1),
        _m_2(_traits.m_2)
    {
    }

    template <typename Process_>
    FormFactors<OneHalfPlusToOneHalfPlus> *
    BMRvD2022FormFactors<Process_>::make(const Parameters & parameters, const Options & options)
    {
        return new BMRvD2022FormFactors(parameters, options);
    }

    template <typename Process_>
    QualifiedName
    BMRvD2022FormFactors<Process_>::_par_name(const std::string & pol, const std::string & current, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^(" + pol + "," + current + ")_" + stringify(idx) + "@BMRvD2022");
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_phi(const double & s, const double & chi, const double & s_p, const double & a,
                                         const double & b, const double & c, const double & d,
                                         const double & e, const double & f, const double & g) const
    {
        // general form:
        // phi = (mLb + mL)^a * (mLb - mL)^b / sqrt((16 + 8 * c) * d * pi^2 * chi)
        //     * (s_- / z(s, t_-))^(e / 4)
        //     * (s_+)^(f / 4)
        //     * (- z(s, 0) / s)^(g / 2)
        //     * sqrt(4.0 * (t_p - t_0)) * (1 + z(s, t_0))^(1/2) * (1 - z(s, t_0))^(-3/2)

        using std::abs;
        using std::pow;
        using std::sqrt;

        const double norm   = sqrt(4.0 * (s_p - _traits.t0)) * sqrt(1 + _traits.calc_z(s, s_p, _traits.t0)) * pow(1 - _traits.calc_z(s, s_p, _traits.t0), -3.0 / 2.0)
                            / sqrt((16.0 + 8.0 * c) * d * M_PI * M_PI * chi);
        const double base_a = _m_1 + _m_2;
        const double base_b = _m_1 - _m_2;
        const double base_e = abs(_traits.tm() - s) > 1.0e-7 ? (_traits.tm() - s) / _traits.calc_z(s, s_p, _traits.tm()) : 4.0 * (s_p - _traits.tm());
        const double base_f = power_of<2>(_m_1 + _m_2) - s;
        const double base_g = abs(s) > 1.0e-7 ? -1.0 * _traits.calc_z(s, s_p, 0.0) / s : 1.0 / (4.0 * s_p);

        return norm * pow(base_a, a) * pow(base_b, b) * pow(base_e, e / 4.0)
                * pow(base_f, f / 4.0) * pow(base_g, g / 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_time_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_0p, _traits.tp_v, 0.0, 1.0, 0.0, 1.0,       1.0, 3.0, 3.0 + 1.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m, _traits.tp_v, 1.0, 0.0, 1.0, 2.0,       3.0, 1.0, 3.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m, _traits.tp_v, 0.0, 0.0, 1.0, 1.0,       3.0, 1.0, 2.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_time_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_0m, _traits.tp_a, 1.0, 0.0, 1.0, 2.0 / 3.0, 3.0, 1.0, 3.0 + 1.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p, _traits.tp_a, 0.0, 1.0, 0.0, 3.0,       1.0, 3.0, 3.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p, _traits.tp_a, 0.0, 0.0, 1.0, 1.0,       1.0, 3.0, 2.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_t,  _traits.tp_v, 0.0, 0.0, 1.0, 2.0,       3.0, 1.0, 1.0 + 3.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_t,  _traits.tp_v, 1.0, 0.0, 1.0, 1.0,       3.0, 1.0, 2.0 + 3.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_t5, _traits.tp_a, 0.0, 0.0, 1.0, 2.0,       1.0, 3.0, 1.0 + 3.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_t5, _traits.tp_a, 0.0, 1.0, 1.0, 1.0,       1.0, 3.0, 2.0 + 3.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_time_v_0() const
    {
        const double x_time = _phi_time_v(0.0) * (power_of<2>(_traits.m_R_0p) < _traits.tp_v ? _traits.calc_z(0.0, _traits.tp_v, power_of<2>(_traits.m_R_0p)) : 1.0);
        const double x_long = _phi_long_v(0.0) * (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(0.0, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        std::array<double, 5> a;
        a[0] = x_time * _a_long_v[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_time * _a_long_v[i] - x_long * _a_time_v[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_v(_traits.calc_z(0.0, _traits.tp_v, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_time_a_0() const
    {
        const double x_time = _phi_time_a(0.0) * (power_of<2>(_traits.m_R_0m) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_0m)) : 1.0);
        const double x_long = _phi_long_a(0.0) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        std::array<double, 5> a;
        a[0] = x_time * _a_long_a[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_time * _a_long_a[i] - x_long * _a_time_a[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_a(_traits.calc_z(0.0, _traits.tp_a, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_perp_a_0() const
    {
        const double x_perp = _phi_perp_a(_traits.tm()) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double x_long = _phi_long_a(_traits.tm()) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        std::array<double, 5> a;
        a[0] = x_perp * _a_long_a[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp * _a_long_a[i] - x_long * _a_perp_a[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_a(_traits.calc_z(_traits.tm(), _traits.tp_a, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_perp_t_0() const
    {
        const double x_perp_t  = _phi_perp_t (0.0) * (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(0.0, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double x_perp_t5 = _phi_perp_t5(0.0) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const auto polynomials_perp_t  = _traits.orthonormal_polynomials_v(_traits.calc_z(0.0, _traits.tp_v, _traits.t0));
        const auto polynomials_perp_t5 = _traits.orthonormal_polynomials_a(_traits.calc_z(0.0, _traits.tp_a, _traits.t0));

        double a_perp_t_0 = x_perp_t * _a_perp_t5[0] * polynomials_perp_t5[0];
        for (unsigned i = 1 ; i < _a_perp_t5.size() ; ++i)
        {
            a_perp_t_0 += x_perp_t * _a_perp_t5[i] * polynomials_perp_t5[i] - x_perp_t5 * _a_perp_t[i - 1] * polynomials_perp_t[i];
        }

        return a_perp_t_0 / (polynomials_perp_t[0] * x_perp_t5);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_long_t5_0() const
    {
        const double x_long_t5 = _phi_long_t5(_traits.tm()) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double x_perp_t5 = _phi_perp_t5(_traits.tm()) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        std::array<double, 5> a;
        a[0] = x_long_t5 * _a_perp_t5[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_long_t5 * _a_perp_t5[i] - x_perp_t5 * _a_long_t5[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_a(_traits.calc_z(_traits.tm(), _traits.tp_a, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp_t5);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_time_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time_v_0();
        std::copy(_a_time_v.begin(), _a_time_v.end(), coefficients.begin() + 1);
        // resonances for 0^+
        const double blaschke     = (power_of<2>(_traits.m_R_0p) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_0p)) : 1.0);
        const double phi          = _phi_time_v(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_long_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_long_v.begin(), _a_long_v.end(), coefficients.begin());
        // resonances for 1^-
        const double blaschke     = (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_long_v(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_perp_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_perp_v.begin(), _a_perp_v.end(), coefficients.begin());
        // resonances for 1^-
        const double blaschke     = (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_perp_v(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_time_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time_a_0();
        std::copy(_a_time_a.begin(), _a_time_a.end(), coefficients.begin() + 1);
        // resonances for 0^-
        const double blaschke     = (power_of<2>(_traits.m_R_0m) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_0m)) : 1.0);
        const double phi          = _phi_time_a(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_long_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_long_a.begin(), _a_long_a.end(), coefficients.begin());
        // resonances for 1^+
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_long_a(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_perp_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp_a_0();
        std::copy(_a_perp_a.begin(), _a_perp_a.end(), coefficients.begin() + 1);
        // resonances for 1^+
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_perp_a(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_long_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_long_t.begin(), _a_long_t.end(), coefficients.begin());
        // resonances for T (1^- state)
        const double blaschke     = (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_long_t(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_perp_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp_t_0();
        std::copy(_a_perp_t.begin(), _a_perp_t.end(), coefficients.begin() + 1);
        // resonances for T (1^- state)
        const double blaschke     = (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_perp_t(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_long_t5(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_long_t5_0();
        std::copy(_a_long_t5.begin(), _a_long_t5.end(), coefficients.begin() + 1);
        // resonances for T5 (1^p state)
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_long_t5(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_perp_t5(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_perp_t5.begin(), _a_perp_t5.end(), coefficients.begin());
        // resonances for T5 (1^p state)
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_perp_t5(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_0p_v() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time_v_0();
        std::copy(_a_time_v.begin(), _a_time_v.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_0m_a() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time_a_0();
        std::copy(_a_time_a.begin(), _a_time_a.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_v_0() const
    {
        return std::inner_product(_a_long_v.begin(), _a_long_v.end(), _a_long_v.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_v_perp() const
    {
        // The perp_v form factor contribute equally to 1m_perp and 1m_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]
        return 0.5 * std::inner_product(_a_perp_v.begin(), _a_perp_v.end(), _a_perp_v.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_v_para() const
    {
        // The perp_v form factor contribute equally to 1m_perp and 1m_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]
        return 0.5 * std::inner_product(_a_perp_v.begin(), _a_perp_v.end(), _a_perp_v.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_v() const
    {
        // By convention, the sum is divided by 3 to follow the bound saturation < 1.0
        return (saturation_1m_v_0() + saturation_1m_v_perp() + saturation_1m_v_para()) / 3.0;
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_a_0() const
    {
        return std::inner_product(_a_long_a.begin(), _a_long_a.end(), _a_long_a.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_a_perp() const
    {
        // The perp_a form factor contribute equally to 1p_perp and 1p_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]
        std::array<double, 5> coefficients_perp;
        coefficients_perp[0] = _a_perp_a_0();
        std::copy(_a_perp_a.begin(), _a_perp_a.end(), coefficients_perp.begin() + 1);

        return 0.5 * std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_a_para() const
    {
        // The perp_a form factor contribute equally to 1p_perp and 1p_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]
        std::array<double, 5> coefficients_perp;
        coefficients_perp[0] = _a_perp_a_0();
        std::copy(_a_perp_a.begin(), _a_perp_a.end(), coefficients_perp.begin() + 1);

        return 0.5 * std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_a() const
    {
        // By convention, the sum is divided by 3 to follow the bound saturation < 1.0
        return (saturation_1p_a_0() + saturation_1p_a_perp() + saturation_1p_a_para()) / 3.0;
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_t_0() const
    {
        return std::inner_product(_a_long_t.begin(), _a_long_t.end(), _a_long_t.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_t_perp() const
    {
        // The perp_t form factor contribute equally to 1m_perp and 1m_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]
        std::array<double, 5> coefficients_perp;
        coefficients_perp[0] = _a_perp_t_0();
        std::copy(_a_perp_t.begin(), _a_perp_t.end(), coefficients_perp.begin() + 1);

        return 0.5 * std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_t_para() const
    {
        // The perp_t form factor contribute equally to 1m_perp and 1m_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]
        std::array<double, 5> coefficients_perp;
        coefficients_perp[0] = _a_perp_t_0();
        std::copy(_a_perp_t.begin(), _a_perp_t.end(), coefficients_perp.begin() + 1);

        return 0.5 * std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1m_t() const
    {
        // By convention, the sum is divided by 3 to follow the bound saturation < 1.0
        return (saturation_1m_t_0() + saturation_1m_t_perp() + saturation_1m_t_para()) / 3.0;
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_t5_0() const
    {
        std::array<double, 5> coefficients_long;
        coefficients_long[0] = _a_long_t5_0();
        std::copy(_a_long_t5.begin(), _a_long_t5.end(), coefficients_long.begin() + 1);

        return std::inner_product(coefficients_long.begin(), coefficients_long.end(), coefficients_long.begin(), 0.0);
    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_t5_perp() const
    {
        // The perp_t5 form factor contribute equally to 1p_perp and 1p_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]

        return 0.5 * std::inner_product(_a_perp_t5.begin(), _a_perp_t5.end(), _a_perp_t5.begin(), 0.0);

    }

    template<typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_t5_para() const
    {
        // The perp_t5 form factor contribute equally to 1p_perp and 1p_para.
        // The factor of 0.5 compensates the factor of 2.0 in the outer function of [BMRvD:2022A]

        return 0.5 * std::inner_product(_a_perp_t5.begin(), _a_perp_t5.end(), _a_perp_t5.begin(), 0.0);

    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::saturation_1p_t5() const
    {
        // By convention, the sum is divided by 3 to follow the bound saturation < 1.0
        return (saturation_1p_t5_0() + saturation_1p_t5_perp() + saturation_1p_t5_para()) / 3.0;
    }

    template <typename Process_>
    Diagnostics
    BMRvD2022FormFactors<Process_>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _traits.calc_z(0.0, _traits.tp_v, _traits.t0),  "z(q2 =  0)" });
        results.add({ _phi_time_v(0.0),               "phi(q2 =  0, f_time^V)"});
        results.add({ _phi_long_v(0.0),               "phi(q2 =  0, f_long^V)"});
        results.add({ _phi_perp_v(0.0),               "phi(q2 =  0, f_perp^V)"});
        results.add({ _phi_time_a(0.0),               "phi(q2 =  0, f_time^A)"});
        results.add({ _phi_long_a(0.0),               "phi(q2 =  0, f_long^A)"});
        results.add({ _phi_perp_a(0.0),               "phi(q2 =  0, f_perp^A)"});
        results.add({ _phi_long_t(0.0),               "phi(q2 =  0, f_long^T)"});
        results.add({ _phi_perp_t(0.0),               "phi(q2 =  0, f_perp^T)"});
        results.add({ _phi_long_t5(0.0),              "phi(q2 =  0, f_long^T5)"});
        results.add({ _phi_perp_t5(0.0),              "phi(q2 =  0, f_perp^T5)"});

        results.add({ _traits.calc_z(10.0, _traits.tp_v, _traits.t0), "z(q2 = 10)" });
        results.add({ _phi_time_v(10.0),              "phi(q2 = 10, f_time^V)"});
        results.add({ _phi_long_v(10.0),              "phi(q2 = 10, f_long^V)"});
        results.add({ _phi_perp_v(10.0),              "phi(q2 = 10, f_perp^V)"});
        results.add({ _phi_time_a(10.0),              "phi(q2 = 10, f_time^A)"});
        results.add({ _phi_long_a(10.0),              "phi(q2 = 10, f_long^A)"});
        results.add({ _phi_perp_a(10.0),              "phi(q2 = 10, f_perp^A)"});
        results.add({ _phi_long_t(10.0),              "phi(q2 = 10, f_long^T)"});
        results.add({ _phi_perp_t(10.0),              "phi(q2 = 10, f_perp^T)"});
        results.add({ _phi_long_t5(10.0),             "phi(q2 = 10, f_long^T5)"});
        results.add({ _phi_perp_t5(10.0),             "phi(q2 = 10, f_perp^T5)"});

        {
        const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials_v(0.0);
        results.add({ p0,                             "p_0(z = 0.0)" });
        results.add({ p1,                             "p_1(z = 0.0)" });
        results.add({ p2,                             "p_2(z = 0.0)" });
        results.add({ p3,                             "p_3(z = 0.0)" });
        results.add({ p4,                             "p_4(z = 0.0)" });
        results.add({ p5,                             "p_5(z = 0.0)" });
        }

        {
        const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials_v(_traits.calc_z(10.0, _traits.tp_v, _traits.t0));
        results.add({ p0,                             "p_0(z = z(q2 = 10))" });
        results.add({ p1,                             "p_1(z = z(q2 = 10))" });
        results.add({ p2,                             "p_2(z = z(q2 = 10))" });
        results.add({ p3,                             "p_3(z = z(q2 = 10))" });
        results.add({ p4,                             "p_4(z = z(q2 = 10))" });
        results.add({ p5,                             "p_5(z = z(q2 = 10))" });
        }

        return results;
    }

    template <typename Process_>
    const std::set<ReferenceName> BMRvD2022FormFactors<Process_>::references
    {
        "BMRvD:2022A"_rn
    };

    template <typename Process_>
    const std::vector<OptionSpecification> BMRvD2022FormFactors<Process_>::options
    {
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BMRvD2022FormFactors<Process_>::begin_options()
    {
        return options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BMRvD2022FormFactors<Process_>::end_options()
    {
        return options.cend();
    }
}

#endif
