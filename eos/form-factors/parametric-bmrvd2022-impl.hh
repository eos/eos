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
    BMRvD2022FormFactors<Process_>::BMRvD2022FormFactors(const Parameters & p, const Options &) :
        _m_1(Process_::m1),
        _m_2(Process_::m2),
        _t_0(Process_::t0),
        _t_m(Process_::tm),
        _t_p(Process_::tp),
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
        }
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
    BMRvD2022FormFactors<Process_>::_z(const double & t, const double & t_0) const
    {
        return (std::sqrt(_t_p - t) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - t) + std::sqrt(_t_p - t_0));
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_phi(const double & s, const double & chi, const double & a, const double & b,
                                         const double & c, const double & d,
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

        const double norm   = sqrt(4.0 * (_t_p - _t_0)) * sqrt(1 + _z(s, _t_0)) * pow(1 - _z(s, _t_0), -3.0 / 2.0)
                            / sqrt((16.0 + 8.0 * c) * d * M_PI * M_PI * chi);
        const double base_a = _m_1 + _m_2;
        const double base_b = _m_1 - _m_2;
        const double base_e = abs(_t_m - s) > 1.0e-7 ? (_t_m - s) / _z(s, _t_m) : 4.0 * (_t_p - _t_m);
        const double base_f = power_of<2>(_m_1 + _m_2) - s;
        const double base_g = abs(s) > 1.0e-7 ? -1.0 * _z(s, 0.0) / s : 1.0 / (4.0 * _t_p);

        return norm * pow(base_a, a) * pow(base_b, b) * pow(base_e, e / 4.0)
                * pow(base_f, f / 4.0) * pow(base_g, g / 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_time_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_0p, 0.0, 1.0, 0.0, 1.0,       1.0, 3.0, 3.0 + 1.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m, 1.0, 0.0, 1.0, 2.0,       3.0, 1.0, 3.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m, 0.0, 0.0, 1.0, 1.0,       3.0, 1.0, 2.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_time_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_0m, 1.0, 0.0, 1.0, 2.0 / 3.0, 3.0, 1.0, 3.0 + 1.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p, 0.0, 1.0, 0.0, 3.0,       1.0, 3.0, 3.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p, 0.0, 0.0, 1.0, 1.0,       1.0, 3.0, 2.0 + 2.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_t,  0.0, 0.0, 1.0, 2.0,       3.0, 1.0, 1.0 + 3.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_t,  1.0, 0.0, 1.0, 1.0,       3.0, 1.0, 2.0 + 3.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_long_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_t5, 0.0, 0.0, 1.0, 2.0,       1.0, 3.0, 1.0 + 3.0);
    }

    template <typename Process_>
    inline double
    BMRvD2022FormFactors<Process_>::_phi_perp_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_t5, 0.0, 1.0, 1.0, 1.0,       1.0, 3.0, 2.0 + 3.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_time_v_0() const
    {
        const double x_time = this->_z(0.0, Process_::mR2_0p) * this->_phi_time_v(0.0);
        const double x_long = this->_z(0.0, Process_::mR2_1m) * this->_phi_long_v(0.0);
        std::array<double, 5> a;
        a[0] = x_time * this->_a_long_v[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_time * this->_a_long_v[i] - x_long * this->_a_time_v[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, this->_t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_time_a_0() const
    {
        const double x_time = this->_z(0.0, Process_::mR2_0m) * this->_phi_time_a(0.0);
        const double x_long = this->_z(0.0, Process_::mR2_1p) * this->_phi_long_a(0.0);
        std::array<double, 5> a;
        a[0] = x_time * this->_a_long_a[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_time * this->_a_long_a[i] - x_long * this->_a_time_a[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, this->_t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_perp_a_0() const
    {
        const double x_perp = this->_z(Process_::tm, Process_::mR2_1p) * this->_phi_perp_a(Process_::tm);
        const double x_long = this->_z(Process_::tm, Process_::mR2_1p) * this->_phi_long_a(Process_::tm);
        std::array<double, 5> a;
        a[0] = x_perp * this->_a_long_a[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp * this->_a_long_a[i] - x_long * this->_a_perp_a[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(Process_::tm, this->_t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_perp_t_0() const
    {
        const double x_perp_t  = this->_z(0.0, Process_::mR2_1m) * this->_phi_perp_t (0.0);
        const double x_perp_t5 = this->_z(0.0, Process_::mR2_1p) * this->_phi_perp_t5(0.0);
        std::array<double, 5> a;
        a[0] = x_perp_t * this->_a_perp_t5[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp_t * this->_a_perp_t5[i] - x_perp_t5 * this->_a_perp_t[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, this->_t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp_t5);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::_a_long_t5_0() const
    {
        const double x_long_t5 = this->_z(Process_::tm, Process_::mR2_1p) * this->_phi_long_t5(Process_::tm);
        const double x_perp_t5 = this->_z(Process_::tm, Process_::mR2_1p) * this->_phi_perp_t5(Process_::tm);
        std::array<double, 5> a;
        a[0] = x_long_t5 * this->_a_perp_t5[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_long_t5 * this->_a_perp_t5[i] - x_perp_t5 * this->_a_long_t5[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(Process_::tm, this->_t_0));
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
        const double blaschke     = _z(q2, Process_::mR2_0p);
        const double phi          = _phi_time_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_long_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_perp_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_0m);
        const double phi          = _phi_time_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_long_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_perp_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_long_t(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_perp_t(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
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
        // no resonances for T5
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_long_t5(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::f_perp_t5(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_perp_t5.begin(), _a_perp_t5.end(), coefficients.begin());
        // no resonances for T5
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_perp_t5(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_0p() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time_v_0();
        std::copy(_a_time_v.begin(), _a_time_v.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_1m() const
    {
        std::array<double, 5> coefficients_long;
        std::copy(_a_long_v.begin(), _a_long_v.end(), coefficients_long.begin());

        std::array<double, 5> coefficients_perp;
        std::copy(_a_perp_v.begin(), _a_perp_v.end(), coefficients_perp.begin());

        return std::inner_product(coefficients_long.begin(), coefficients_long.end(), coefficients_long.begin(), 0.0)
                + std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_0m() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time_a_0();
        std::copy(_a_time_a.begin(), _a_time_a.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_1p() const
    {
        std::array<double, 5> coefficients_long;
        std::copy(_a_long_a.begin(), _a_long_a.end(), coefficients_long.begin());

        std::array<double, 5> coefficients_perp;
        coefficients_perp[0] = _a_perp_a_0();
        std::copy(_a_perp_a.begin(), _a_perp_a.end(), coefficients_perp.begin() + 1);

        return std::inner_product(coefficients_long.begin(), coefficients_long.end(), coefficients_long.begin(), 0.0)
                + std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_T() const
    {
        std::array<double, 5> coefficients_long;
        std::copy(_a_long_t.begin(), _a_long_t.end(), coefficients_long.begin());

        std::array<double, 5> coefficients_perp;
        coefficients_perp[0] = _a_perp_t_0();
        std::copy(_a_perp_t.begin(), _a_perp_t.end(), coefficients_perp.begin() + 1);

        return std::inner_product(coefficients_long.begin(), coefficients_long.end(), coefficients_long.begin(), 0.0)
                + std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_T5() const
    {
        std::array<double, 5> coefficients_long;
        coefficients_long[0] = _a_long_t5_0();
        std::copy(_a_long_t5.begin(), _a_long_t5.end(), coefficients_long.begin() + 1);

        std::array<double, 5> coefficients_perp;
        std::copy(_a_perp_t5.begin(), _a_perp_t5.end(), coefficients_perp.begin());

        return std::inner_product(coefficients_long.begin(), coefficients_long.end(), coefficients_long.begin(), 0.0)
                + std::inner_product(coefficients_perp.begin(), coefficients_perp.end(), coefficients_perp.begin(), 0.0);
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_0p_prior() const
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
            static const double sigma = 0.01; // 10% uncertainty
            return -power_of<2>((value - 1.0) / sigma) / 2.0;
        }
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_1m_prior() const
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
            static const double sigma = 0.01; // 10% uncertainty
            return -power_of<2>((value - 1.0) / sigma) / 2.0;
        }
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_0m_prior() const
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
            static const double sigma = 0.01; // 10% uncertainty
            return -power_of<2>((value - 1.0) / sigma) / 2.0;
        }
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_1p_prior() const
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
            static const double sigma = 0.01; // 10% uncertainty
            return -power_of<2>((value - 1.0) / sigma) / 2.0;
        }
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_T_prior() const
    {
        const double value = bound_T();

        if (value < 0.0)
        {
            throw InternalError("Contribution to T unitarity bound must be positive; found to be negative!");
        }
        else if ((0.0 <= value) && (value < 1.0))
        {
            return 0.0;
        }
        else
        {
            // add an r-fit like penalty
            static const double sigma = 0.01; // 10% uncertainty
            return -power_of<2>((value - 1.0) / sigma) / 2.0;
        }
    }

    template <typename Process_>
    double
    BMRvD2022FormFactors<Process_>::bound_T5_prior() const
    {
        const double value = bound_T5();

        if (value < 0.0)
        {
            throw InternalError("Contribution to T5 unitarity bound must be positive; found to be negative!");
        }
        else if ((0.0 <= value) && (value < 1.0))
        {
            return 0.0;
        }
        else
        {
            // add an r-fit like penalty
            static const double sigma = 0.01; // 10% uncertainty
            return -power_of<2>((value - 1.0) / sigma) / 2.0;
        }
    }

    template <typename Process_>
    Diagnostics
    BMRvD2022FormFactors<Process_>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _z(0.0, Process_::t0),  "z(q2 =  0)" });
        results.add({ _phi_time_v(0.0),       "phi(q2 =  0, f_time^V)"});
        results.add({ _phi_long_v(0.0),       "phi(q2 =  0, f_long^V)"});
        results.add({ _phi_perp_v(0.0),       "phi(q2 =  0, f_perp^V)"});
        results.add({ _phi_time_a(0.0),       "phi(q2 =  0, f_time^A)"});
        results.add({ _phi_long_a(0.0),       "phi(q2 =  0, f_long^A)"});
        results.add({ _phi_perp_a(0.0),       "phi(q2 =  0, f_perp^A)"});
        results.add({ _phi_long_t(0.0),       "phi(q2 =  0, f_long^T)"});
        results.add({ _phi_perp_t(0.0),       "phi(q2 =  0, f_perp^T)"});
        results.add({ _phi_long_t5(0.0),      "phi(q2 =  0, f_long^T5)"});
        results.add({ _phi_perp_t5(0.0),      "phi(q2 =  0, f_perp^T5)"});

        results.add({ _z(10.0, Process_::t0), "z(q2 = 10)" });
        results.add({ _phi_time_v(10.0),      "phi(q2 = 10, f_time^V)"});
        results.add({ _phi_long_v(10.0),      "phi(q2 = 10, f_long^V)"});
        results.add({ _phi_perp_v(10.0),      "phi(q2 = 10, f_perp^V)"});
        results.add({ _phi_time_a(10.0),      "phi(q2 = 10, f_time^A)"});
        results.add({ _phi_long_a(10.0),      "phi(q2 = 10, f_long^A)"});
        results.add({ _phi_perp_a(10.0),      "phi(q2 = 10, f_perp^A)"});
        results.add({ _phi_long_t(10.0),      "phi(q2 = 10, f_long^T)"});
        results.add({ _phi_perp_t(10.0),      "phi(q2 = 10, f_perp^T)"});
        results.add({ _phi_long_t5(10.0),     "phi(q2 = 10, f_long^T5)"});
        results.add({ _phi_perp_t5(10.0),     "phi(q2 = 10, f_perp^T5)"});

        {
        const auto & [p0, p1, p2, p3, p4, p5] = Process_::orthonormal_polynomials(0.0);
        results.add({ p0,                     "p_0(z = 0.0)" });
        results.add({ p1,                     "p_1(z = 0.0)" });
        results.add({ p2,                     "p_2(z = 0.0)" });
        results.add({ p3,                     "p_3(z = 0.0)" });
        results.add({ p4,                     "p_4(z = 0.0)" });
        results.add({ p5,                     "p_5(z = 0.0)" });
        }

        {
        const auto & [p0, p1, p2, p3, p4, p5] = Process_::orthonormal_polynomials(_z(10.0, Process_::t0));
        results.add({ p0,                     "p_0(z = z(q2 = 10))" });
        results.add({ p1,                     "p_1(z = z(q2 = 10))" });
        results.add({ p2,                     "p_2(z = z(q2 = 10))" });
        results.add({ p3,                     "p_3(z = z(q2 = 10))" });
        results.add({ p4,                     "p_4(z = z(q2 = 10))" });
        results.add({ p5,                     "p_5(z = z(q2 = 10))" });
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
