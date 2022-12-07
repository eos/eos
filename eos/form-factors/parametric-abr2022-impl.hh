/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_ABR2022_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_ABR2022_IMPL_HH 1

#include <eos/form-factors/parametric-abr2022.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>

#include <numeric>

#include <iostream>

namespace eos
{
    template <typename Process_>
    ABR2022FormFactors<Process_>::ABR2022FormFactors(const Parameters & p, const Options &) :
        _m_1(Process_::m1),
        _m_2(Process_::m2),
        _t_0(Process_::t0),
        _t_m(Process_::tm),
        _t_p(Process_::tp),
        _a_time12_v{
            UsedParameter(p[_par_name("t12", "V", 1)], *this),
            UsedParameter(p[_par_name("t12", "V", 2)], *this),
            UsedParameter(p[_par_name("t12", "V", 3)], *this),
            UsedParameter(p[_par_name("t12", "V", 4)], *this)
        },
        _a_long12_v{
            UsedParameter(p[_par_name("012", "V", 1)], *this),
            UsedParameter(p[_par_name("012", "V", 2)], *this),
            UsedParameter(p[_par_name("012", "V", 3)], *this),
            UsedParameter(p[_par_name("012", "V", 4)], *this)
        },
        _a_perp12_v{
            UsedParameter(p[_par_name("perp12", "V", 1)], *this),
            UsedParameter(p[_par_name("perp12", "V", 2)], *this),
            UsedParameter(p[_par_name("perp12", "V", 3)], *this),
            UsedParameter(p[_par_name("perp12", "V", 4)], *this)
        },
        _a_perp32_v{
            UsedParameter(p[_par_name("perp32", "V", 0)], *this),
            UsedParameter(p[_par_name("perp32", "V", 1)], *this),
            UsedParameter(p[_par_name("perp32", "V", 2)], *this),
            UsedParameter(p[_par_name("perp32", "V", 3)], *this),
            UsedParameter(p[_par_name("perp32", "V", 4)], *this)
        },
        _a_time12_a{
            UsedParameter(p[_par_name("t12", "A", 1)], *this),
            UsedParameter(p[_par_name("t12", "A", 2)], *this),
            UsedParameter(p[_par_name("t12", "A", 3)], *this),
            UsedParameter(p[_par_name("t12", "A", 4)], *this)
        },
        _a_long12_a{
            UsedParameter(p[_par_name("012", "A", 1)], *this),
            UsedParameter(p[_par_name("012", "A", 2)], *this),
            UsedParameter(p[_par_name("012", "A", 3)], *this),
            UsedParameter(p[_par_name("012", "A", 4)], *this)
        },
        _a_perp12_a{
            UsedParameter(p[_par_name("perp12", "A", 1)], *this),
            UsedParameter(p[_par_name("perp12", "A", 2)], *this),
            UsedParameter(p[_par_name("perp12", "A", 3)], *this),
            UsedParameter(p[_par_name("perp12", "A", 4)], *this)
        },
        _a_perp32_a{
            UsedParameter(p[_par_name("perp32", "A", 0)], *this),
            UsedParameter(p[_par_name("perp32", "A", 1)], *this),
            UsedParameter(p[_par_name("perp32", "A", 2)], *this),
            UsedParameter(p[_par_name("perp32", "A", 3)], *this),
            UsedParameter(p[_par_name("perp32", "A", 4)], *this)
        },
        _a_long12_t{
            UsedParameter(p[_par_name("012", "T", 1)], *this),
            UsedParameter(p[_par_name("012", "T", 2)], *this),
            UsedParameter(p[_par_name("012", "T", 3)], *this),
            UsedParameter(p[_par_name("012", "T", 4)], *this)
        },
        _a_perp12_t{
            UsedParameter(p[_par_name("perp12", "T", 1)], *this),
            UsedParameter(p[_par_name("perp12", "T", 2)], *this),
            UsedParameter(p[_par_name("perp12", "T", 3)], *this),
            UsedParameter(p[_par_name("perp12", "T", 4)], *this)
        },
        _a_perp32_t{
            UsedParameter(p[_par_name("perp32", "T", 0)], *this),
            UsedParameter(p[_par_name("perp32", "T", 1)], *this),
            UsedParameter(p[_par_name("perp32", "T", 2)], *this),
            UsedParameter(p[_par_name("perp32", "T", 3)], *this),
            UsedParameter(p[_par_name("perp32", "T", 4)], *this)
        },
        _a_long12_t5{
            UsedParameter(p[_par_name("012", "T5", 1)], *this),
            UsedParameter(p[_par_name("012", "T5", 2)], *this),
            UsedParameter(p[_par_name("012", "T5", 3)], *this),
            UsedParameter(p[_par_name("012", "T5", 4)], *this)
        },
        _a_perp12_t5{
            UsedParameter(p[_par_name("perp12", "T5", 1)], *this),
            UsedParameter(p[_par_name("perp12", "T5", 2)], *this),
            UsedParameter(p[_par_name("perp12", "T5", 3)], *this),
            UsedParameter(p[_par_name("perp12", "T5", 4)], *this)
        },
        _a_perp32_t5{
            UsedParameter(p[_par_name("perp32", "T5", 1)], *this),
            UsedParameter(p[_par_name("perp32", "T5", 2)], *this),
            UsedParameter(p[_par_name("perp32", "T5", 3)], *this),
            UsedParameter(p[_par_name("perp32", "T5", 4)], *this)
        }
    {
    }

    template <typename Process_>
    FormFactors<OneHalfPlusToThreeHalfMinus> *
    ABR2022FormFactors<Process_>::make(const Parameters & parameters, const Options & options)
    {
        return new ABR2022FormFactors(parameters, options);
    }

    template <typename Process_>
    QualifiedName
    ABR2022FormFactors<Process_>::_par_name(const std::string & pol, const std::string & current, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^(" + pol + "," + current + ")_" + stringify(idx) + "@ABR2022");
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_z(const double & t, const double & t_0) const
    {
        return (std::sqrt(_t_p - t) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - t) + std::sqrt(_t_p - t_0));
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_phi(const double & s, const double & chi, const double & A, const double & B,
                                       const double & d, const double & e, const double & f, const double & g,
                                       const double & n) const
    {
        using std::pow;
        using std::sqrt;

        const double z = _z(s, _t_0);
        const double t_pb = power_of<2>(_m_1 + _m_2);
        const double t_mb = power_of<2>(_m_1 - _m_2);

        const double norm   = pow(t_pb, 0.5 * A) * pow(t_mb, 0.5 * B) * pow(1 - z, n + g - 0.5 * (e + f + 3.0)) * sqrt(1 + z)
                            * sqrt(4.0 * (_t_p - _t_0)) / sqrt(16.0 * d * M_PI * M_PI * chi);
        const double phi1 = - 1.0 / (_t_0 * power_of<2>(1 + z) - 2.0 * _t_p * (1 + z*z) - 2.0 * (1 - z*z) * sqrt(_t_p) * sqrt(_t_p - _t_0));
        const double phi2 = - power_of<2>(1 - z) * t_mb - power_of<2>(1 + z) * _t_0 + 2.0 * _t_p * (1 + z*z) + 2.0 * (1 - z*z) * sqrt(_t_p - t_mb) * sqrt(_t_p - _t_0);
        const double phi3 = power_of<2>(1 - z) * t_pb - power_of<2>(1 + z) * _t_0 + 4.0 * _t_p * z;

        return norm * pow(phi1, 0.5 * (n + g)) * pow(phi2, e / 4.0) * pow(phi3, f / 4.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_time12_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_0p_v, 0.0, 1.0,  6.0, 3.0, 1.0, 3.0, 1.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_long12_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m_v, 1.0, 0.0, 18.0, 1.0, 3.0, 3.0, 2.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp12_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m_v, 0.0, 0.0,  9.0, 1.0, 3.0, 2.0, 2.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp32_v(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m_v, 0.0, 0.0,  3.0, 1.0, 3.0, 2.0, 2.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_time12_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_0m_a, 1.0, 0.0,  6.0, 1.0, 3.0, 3.0, 1.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_long12_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p_a, 0.0, 1.0, 18.0, 3.0, 1.0, 3.0, 2.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp12_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p_a, 0.0, 0.0,  9.0, 3.0, 1.0, 2.0, 2.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp32_a(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p_a, 0.0, 0.0,  3.0, 3.0, 1.0, 2.0, 2.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_long12_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m_t, 0.0, 0.0, 18.0, 1.0, 3.0, 1.0, 3.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp12_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m_t, 1.0, 0.0,  9.0, 1.0, 3.0, 2.0, 3.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp32_t(const double & q2) const
    {
        return _phi(q2, Process_::chi_1m_t, 1.0, 0.0,  3.0, 1.0, 3.0, 2.0, 3.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_long12_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p_t5, 0.0, 0.0, 18.0, 3.0, 1.0, 1.0, 3.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp12_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p_t5, 0.0, 1.0,  9.0, 3.0, 1.0, 2.0, 3.0);
    }

    template <typename Process_>
    inline double
    ABR2022FormFactors<Process_>::_phi_perp32_t5(const double & q2) const
    {
        return _phi(q2, Process_::chi_1p_t5, 0.0, 1.0,  3.0, 3.0, 1.0, 2.0, 3.0);
    }


    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_long12_v_0() const
    {
        const double x_long12_v = this->_phi_long12_v(_t_m) * 2.0 * (_m_1 - _m_2) / (_m_1 + _m_2);
        const double x_perp32_v = this->_phi_perp32_v(_t_m);
        std::array<double, 5> a;
        a[0] = x_long12_v * this->_a_perp32_v[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_long12_v * this->_a_perp32_v[i] - x_perp32_v * this->_a_long12_v[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp32_v);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_perp12_v_0() const
    {
        const double x_perp12_v = - this->_phi_perp12_v(_t_m);
        const double x_perp32_v =   this->_phi_perp32_v(_t_m);
        std::array<double, 5> a;
        a[0] = x_perp12_v * this->_a_perp32_v[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp12_v * this->_a_perp32_v[i] - x_perp32_v * this->_a_perp12_v[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp32_v);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_time12_a_0() const
    {
        std::array<double, 5> a;
        a[0] = 0;
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = - this->_a_time12_a[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / polynomials[0];
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_long12_t_0() const
    {
        const double x_long12_t = this->_phi_long12_t(_t_m) * 2.0 * (_m_1 + _m_2) / (_m_1 - _m_2);
        const double x_perp32_t = this->_phi_perp32_t(_t_m);
        std::array<double, 5> a;
        a[0] = x_long12_t * this->_a_perp32_t[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_long12_t * this->_a_perp32_t[i] - x_perp32_t * this->_a_long12_t[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp32_t);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_perp12_t_0() const
    {
        const double x_perp12_t = - this->_phi_perp12_t(_t_m);
        const double x_perp32_t =   this->_phi_perp32_t(_t_m);
        std::array<double, 5> a;
        a[0] = x_perp12_t * this->_a_perp32_t[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp12_t * this->_a_perp32_t[i] - x_perp32_t * this->_a_perp12_t[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp32_t);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_perp32_t5_0() const
    {
        const double x_perp32_t5 = - this->_z(0.0, Process_::mR2_1p) * this->_phi_perp32_t5(0.0) * power_of<2>((_m_1 + _m_2) / (_m_1 - _m_2));
        const double x_perp32_t  =   this->_z(0.0, Process_::mR2_1m) * this->_phi_perp32_t(0.0);
        std::array<double, 5> a;
        a[0] = x_perp32_t5 * this->_a_perp32_t[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp32_t5 * this->_a_perp32_t[i] - x_perp32_t * this->_a_perp32_t5[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp32_t);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_time12_v_0() const
    {
        const double x_time12_v = this->_z(0.0, Process_::mR2_0p) * this->_phi_time12_v(0.0) * power_of<2>((_m_1 + _m_2) / (_m_1 - _m_2));
        const double x_long12_v = this->_z(0.0, Process_::mR2_1m) * this->_phi_long12_v(0.0);
        std::array<double, 5> a;
        a[0] = x_time12_v * this->_a_long12_v_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_time12_v * this->_a_long12_v[i - 1] - x_long12_v * this->_a_time12_v[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_long12_v);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_long12_a_0() const
    {
        const double x_long12_a = this->_z(0.0, Process_::mR2_1p) * this->_phi_long12_a(0.0) * power_of<2>((_m_1 + _m_2) / (_m_1 - _m_2));
        const double x_time12_a = this->_z(0.0, Process_::mR2_0m) * this->_phi_time12_a(0.0);
        std::array<double, 5> a;
        a[0] = x_long12_a * this->_a_time12_a_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_long12_a * this->_a_time12_a[i - 1] - x_time12_a * this->_a_long12_a[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_time12_a);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_perp12_t5_0() const
    {
        const double x_perp12_t5 = this->_z(0.0, Process_::mR2_1p) * this->_phi_perp12_t5(0.0) * power_of<2>((_m_1 + _m_2) / (_m_1 - _m_2));
        const double x_perp12_t  = this->_z(0.0, Process_::mR2_1m) * this->_phi_perp12_t(0.0);
        std::array<double, 5> a;
        a[0] = x_perp12_t5 * this->_a_perp12_t_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp12_t5 * this->_a_perp12_t[i - 1] - x_perp12_t * this->_a_perp12_t5[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(0.0, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_perp12_t);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_perp12_a_0() const
    {
        const double x_perp12_a = this->_phi_perp12_a(_t_m);
        const double x_long12_a = this->_phi_long12_a(_t_m);
        const double x_perp32_a = this->_phi_perp32_a(_t_m);
        std::array<double, 5> a;
        a[0] = x_perp12_a * (this->_a_long12_a_0() / x_long12_a + this->_a_perp32_a[0] / x_perp32_a);
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_perp12_a * (this->_a_long12_a[i - 1] / x_long12_a + this->_a_perp32_a[i] / x_perp32_a) - this->_a_perp12_a[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / polynomials[0];
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::_a_long12_t5_0() const
    {
        const double x_long12_t5 = this->_phi_long12_t5(_t_m);
        const double x_perp12_t5 = this->_phi_perp12_t5(_t_m);
        const double x_perp32_t5 = this->_phi_perp32_t5(_t_m);
        std::array<double, 5> a;
        a[0] = x_long12_t5 * (this->_a_perp12_t5_0() / x_perp12_t5 - this->_a_perp32_t5_0() / x_perp32_t5);
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_long12_t5 * (this->_a_perp12_t5[i - 1] / x_perp12_t5 - this->_a_perp32_t5[i - 1] / x_perp32_t5) - this->_a_long12_t5[i - 1];
        }
        const auto polynomials = Process_::orthonormal_polynomials(_z(_t_m, _t_0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / polynomials[0];
    }


    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_time12_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time12_v_0();
        std::copy(_a_time12_v.begin(), _a_time12_v.end(), coefficients.begin() + 1);
        // resonances for 0^+
        const double blaschke     = _z(q2, Process_::mR2_0p);
        const double phi          = _phi_time12_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_long12_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_long12_v_0();
        std::copy(_a_long12_v.begin(), _a_long12_v.end(), coefficients.begin() + 1);
        // resonances for 1^-
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_long12_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp12_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp12_v_0();
        std::copy(_a_perp12_v.begin(), _a_perp12_v.end(), coefficients.begin() + 1);
        // resonances for 1^-
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_perp12_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp32_v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_perp32_v.begin(), _a_perp32_v.end(), coefficients.begin());
        // resonances for 1^-
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_perp32_v(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_time12_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time12_a_0();
        std::copy(_a_time12_a.begin(), _a_time12_a.end(), coefficients.begin() + 1);
        // resonances for 0^-
        const double blaschke     = _z(q2, Process_::mR2_0m);
        const double phi          = _phi_time12_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_long12_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_long12_a_0();
        std::copy(_a_long12_a.begin(), _a_long12_a.end(), coefficients.begin() + 1);
        // resonances for 1^+
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_long12_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp12_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp12_a_0();
        std::copy(_a_perp12_a.begin(), _a_perp12_a.end(), coefficients.begin() + 1);
        // resonances for 1^+
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_perp12_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp32_a(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_perp32_a.begin(), _a_perp32_a.end(), coefficients.begin());
        // resonances for 1^+
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_perp32_a(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_long12_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_long12_t_0();
        std::copy(_a_long12_t.begin(), _a_long12_t.end(), coefficients.begin() + 1);
        // resonances for T (1^- state)
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_long12_t(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp12_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp12_t_0();
        std::copy(_a_perp12_t.begin(), _a_perp12_t.end(), coefficients.begin() + 1);
        // resonances for T (1^- state)
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_perp12_t(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp32_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_perp32_t.begin(), _a_perp32_t.end(), coefficients.begin() );
        // resonances for T (1^- state)
        const double blaschke     = _z(q2, Process_::mR2_1m);
        const double phi          = _phi_perp32_t(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_long12_t5(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_long12_t5_0();
        std::copy(_a_long12_t5.begin(), _a_long12_t5.end(), coefficients.begin() + 1);
        // resonances for T5 (1^+ state)
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_long12_t5(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp12_t5(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp12_t5_0();
        std::copy(_a_perp12_t5.begin(), _a_perp12_t5.end(), coefficients.begin() + 1);
        // resonances for T5 (1^+ state)
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_perp12_t5(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::f_perp32_t5(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_perp32_t5_0();
        std::copy(_a_perp32_t5.begin(), _a_perp32_t5.end(), coefficients.begin() + 1);
        // resonances for T5 (1^+ state)
        const double blaschke     = _z(q2, Process_::mR2_1p);
        const double phi          = _phi_perp32_t5(q2);
        const double z            = _z(q2, _t_0);
        const auto   polynomials  = Process_::orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }


    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::saturation_0p_v() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time12_v_0();
        std::copy(_a_time12_v.begin(), _a_time12_v.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::saturation_1m_v() const
    {
        std::array<double, 5> coefficients_long12;
        coefficients_long12[0] = _a_long12_v_0();
        std::copy(_a_long12_v.begin(), _a_long12_v.end(), coefficients_long12.begin() + 1);

        std::array<double, 5> coefficients_perp12;
        coefficients_perp12[0] = _a_perp12_v_0();
        std::copy(_a_perp12_v.begin(), _a_perp12_v.end(), coefficients_perp12.begin() + 1);

        std::array<double, 5> coefficients_perp32;
        std::copy(_a_perp32_v.begin(), _a_perp32_v.end(), coefficients_perp32.begin());

        return std::inner_product(coefficients_long12.begin(), coefficients_long12.end(), coefficients_long12.begin(), 0.0)
             + std::inner_product(coefficients_perp12.begin(), coefficients_perp12.end(), coefficients_perp12.begin(), 0.0)
             + std::inner_product(coefficients_perp32.begin(), coefficients_perp32.end(), coefficients_perp32.begin(), 0.0);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::saturation_0m_a() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_time12_a_0();
        std::copy(_a_time12_a.begin(), _a_time12_a.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::saturation_1p_a() const
    {
        std::array<double, 5> coefficients_long12;
        coefficients_long12[0] = _a_long12_a_0();
        std::copy(_a_long12_a.begin(), _a_long12_a.end(), coefficients_long12.begin() + 1);

        std::array<double, 5> coefficients_perp12;
        coefficients_perp12[0] = _a_perp12_a_0();
        std::copy(_a_perp12_a.begin(), _a_perp12_a.end(), coefficients_perp12.begin() + 1);

        std::array<double, 5> coefficients_perp32;
        std::copy(_a_perp32_a.begin(), _a_perp32_a.end(), coefficients_perp32.begin());

        return std::inner_product(coefficients_long12.begin(), coefficients_long12.end(), coefficients_long12.begin(), 0.0)
             + std::inner_product(coefficients_perp12.begin(), coefficients_perp12.end(), coefficients_perp12.begin(), 0.0)
             + std::inner_product(coefficients_perp32.begin(), coefficients_perp32.end(), coefficients_perp32.begin(), 0.0);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::saturation_1m_t() const
    {
        std::array<double, 5> coefficients_long12;
        coefficients_long12[0] = _a_long12_t_0();
        std::copy(_a_long12_t.begin(), _a_long12_t.end(), coefficients_long12.begin() + 1);

        std::array<double, 5> coefficients_perp12;
        coefficients_perp12[0] = _a_perp12_t_0();
        std::copy(_a_perp12_t.begin(), _a_perp12_t.end(), coefficients_perp12.begin() + 1);

        std::array<double, 5> coefficients_perp32;
        std::copy(_a_perp32_t.begin(), _a_perp32_t.end(), coefficients_perp32.begin());

        return std::inner_product(coefficients_long12.begin(), coefficients_long12.end(), coefficients_long12.begin(), 0.0)
             + std::inner_product(coefficients_perp12.begin(), coefficients_perp12.end(), coefficients_perp12.begin(), 0.0)
             + std::inner_product(coefficients_perp32.begin(), coefficients_perp32.end(), coefficients_perp32.begin(), 0.0);
    }

    template <typename Process_>
    double
    ABR2022FormFactors<Process_>::saturation_1p_t5() const
    {
        std::array<double, 5> coefficients_long12;
        coefficients_long12[0] = _a_long12_t5_0();
        std::copy(_a_long12_t5.begin(), _a_long12_t5.end(), coefficients_long12.begin() + 1);

        std::array<double, 5> coefficients_perp12;
        coefficients_perp12[0] = _a_perp12_t5_0();
        std::copy(_a_perp12_t5.begin(), _a_perp12_t5.end(), coefficients_perp12.begin() + 1);

        std::array<double, 5> coefficients_perp32;
        coefficients_perp32[0] = _a_perp32_t5_0();
        std::copy(_a_perp32_t5.begin(), _a_perp32_t5.end(), coefficients_perp32.begin() + 1);

        return std::inner_product(coefficients_long12.begin(), coefficients_long12.end(), coefficients_long12.begin(), 0.0)
             + std::inner_product(coefficients_perp12.begin(), coefficients_perp12.end(), coefficients_perp12.begin(), 0.0)
             + std::inner_product(coefficients_perp32.begin(), coefficients_perp32.end(), coefficients_perp32.begin(), 0.0);
    }


    template <typename Process_>
    Diagnostics
    ABR2022FormFactors<Process_>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _z(0.0, Process_::t0),  "z(q2 =  0)" });
        results.add({ _z(10.0, Process_::t0), "z(q2 = 10)" });

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

        {
            results.add({ _phi_time12_v(1.0),      "phi_time12_v(z = z(q2 = 1))" });
            results.add({ _phi_long12_v(1.0),      "phi_long12_v(z = z(q2 = 1))" });
            results.add({ _phi_perp12_v(1.0),      "phi_perp12_v(z = z(q2 = 1))" });
            results.add({ _phi_perp32_v(1.0),      "phi_perp32_v(z = z(q2 = 1))" });
            results.add({ _phi_time12_a(1.0),      "phi_time12_a(z = z(q2 = 1))" });
            results.add({ _phi_long12_a(1.0),      "phi_long12_a(z = z(q2 = 1))" });
            results.add({ _phi_perp12_a(1.0),      "phi_perp12_a(z = z(q2 = 1))" });
            results.add({ _phi_perp32_a(1.0),      "phi_perp32_a(z = z(q2 = 1))" });
            results.add({ _phi_long12_t(1.0),      "phi_long12_t(z = z(q2 = 1))" });
            results.add({ _phi_perp12_t(1.0),      "phi_perp12_t(z = z(q2 = 1))" });
            results.add({ _phi_perp32_t(1.0),      "phi_perp32_t(z = z(q2 = 1))" });
            results.add({ _phi_long12_t5(1.0),     "phi_long12_t5(z = z(q2 = 1))" });
            results.add({ _phi_perp12_t5(1.0),     "phi_perp12_t5(z = z(q2 = 1))" });
            results.add({ _phi_perp32_t5(1.0),     "phi_perp32_t5(z = z(q2 = 1))" });
        }

        return results;
    }

    template <typename Process_>
    const std::set<ReferenceName> ABR2022FormFactors<Process_>::references
    {
        "ABR:2022A"_rn
    };

    template <typename Process_>
    const std::vector<OptionSpecification> ABR2022FormFactors<Process_>::options
    {
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    ABR2022FormFactors<Process_>::begin_options()
    {
        return options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    ABR2022FormFactors<Process_>::end_options()
    {
        return options.cend();
    }
}

#endif
