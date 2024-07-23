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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BFW2010_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BFW2010_IMPL_HH 1

#include <eos/form-factors/parametric-bfw2010.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>

#include <numeric>

namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BFW2010FormFactorTraits<Process_, PToV>::resonance_0m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BFW2010FormFactorTraits<Process_, PToV>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BFW2010FormFactorTraits<Process_, PToV>::resonance_1p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,1@BSZ2015" }
    };

    template<typename Process_>
    BFW2010FormFactors<Process_, PToV>::BFW2010FormFactors(const Parameters & p, const Options &) :
        _a_A0{  UsedParameter(p[_par_name("A0", 0)],  *this),
                UsedParameter(p[_par_name("A0", 1)],  *this),
                UsedParameter(p[_par_name("A0", 2)],  *this),
                UsedParameter(p[_par_name("A0", 3)],  *this),
                UsedParameter(p[_par_name("A0", 4)],  *this)
        },
        _a_V{   UsedParameter(p[_par_name("V", 0)],   *this),
                UsedParameter(p[_par_name("V", 1)],   *this),
                UsedParameter(p[_par_name("V", 2)],   *this),
                UsedParameter(p[_par_name("V", 3)],   *this),
                UsedParameter(p[_par_name("V", 4)],   *this)
        },
        _a_T1{  UsedParameter(p[_par_name("T1", 0)],  *this),
                UsedParameter(p[_par_name("T1", 1)],  *this),
                UsedParameter(p[_par_name("T1", 2)],  *this),
                UsedParameter(p[_par_name("T1", 3)],  *this),
                UsedParameter(p[_par_name("T1", 4)],  *this)
        },
        _a_A12{ UsedParameter(p[_par_name("A12", 1)], *this),
                UsedParameter(p[_par_name("A12", 2)], *this),
                UsedParameter(p[_par_name("A12", 3)], *this),
                UsedParameter(p[_par_name("A12", 4)], *this)
        },
        _a_T2{  UsedParameter(p[_par_name("T2", 1)],  *this),
                UsedParameter(p[_par_name("T2", 2)],  *this),
                UsedParameter(p[_par_name("T2", 3)],  *this),
                UsedParameter(p[_par_name("T2", 4)],  *this)
        },
        _a_A1{  UsedParameter(p[_par_name("A1", 1)],  *this),
                UsedParameter(p[_par_name("A1", 2)],  *this),
                UsedParameter(p[_par_name("A1", 3)],  *this),
                UsedParameter(p[_par_name("A1", 4)],  *this)
        },
        _a_T23{ UsedParameter(p[_par_name("T23", 1)], *this),
                UsedParameter(p[_par_name("T23", 2)], *this),
                UsedParameter(p[_par_name("T23", 3)], *this),
                UsedParameter(p[_par_name("T23", 4)], *this)
        },
        _traits(BFW2010FormFactorTraits<Process_, PToV>(p)),
        _mB(_traits.m_B),
        _mV(_traits.m_V)
    {
    }

    template<typename Process_>
    BFW2010FormFactors<Process_, PToV>::~BFW2010FormFactors() = default;

    template<typename Process_>
    FormFactors<PToV> *
    BFW2010FormFactors<Process_, PToV>::make(const Parameters & parameters, const Options & options)
    {
        return new BFW2010FormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    BFW2010FormFactors<Process_, PToV>::_par_name(const std::string & ff_name, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(idx) + "@BFW2010");
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::_phi(const double & t, const double & threshold_tp, const double & chi,
                                             const int & A, const unsigned B, const unsigned C, const unsigned k,
                                             const unsigned p, const unsigned n, const unsigned m) const
    {
        // [GvDV:2022B]
        const double z = _traits.calc_z(t, threshold_tp, _traits.t0),
            kinematic_tp = power_of<2>(_mB + _mV);
        const double norm = std::sqrt(Process_::eta * k * pow(kinematic_tp, A) * pow(_traits.tm(), B)
                                * pow(4 * _mB * _mV, C) / 96 / M_PI / M_PI / chi);

        // set Q^2 to 0
        const double invt = 1 / ( 2.0 * (std::sqrt(threshold_tp) * std::sqrt(threshold_tp - t) + threshold_tp) - t); // simplification of -_traits.calc_z(t, threshold_tp, 0) / t
        const double lambda_term = (kinematic_tp - t) * power_of<2>(std::sqrt(threshold_tp - t) + std::sqrt(threshold_tp - _traits.tm())); // simplification of lambda / z(t, threshold_tp, tm);
        const double sqrtjac = std::sqrt(4 * (1 + z) * (_traits.t0 - threshold_tp) / power_of<3>(z - 1)); // Abs[jacobian] = - jacobian

        return norm * sqrtjac * pow(lambda_term, 0.25 * m) * pow(invt, 0.5 * (p + n + 1.0));
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_v(const double & q2) const
    {
        return _phi(q2, _traits.tp_v, Process_::chi_1m_v, -1, 0, 0, 2, 1, 2, 3);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_a_0(const double & q2) const
    {
        return _phi(q2, _traits.tp_a, Process_::chi_0m_a, 0, 0, 0, 3, 2, 1, 3);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_a_1(const double & q2) const
    {
        return _phi(q2, _traits.tp_a, Process_::chi_1p_a, 1, 0, 0, 2, 1, 2, 1);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_a_12(const double & q2) const
    {
        return _phi(q2, _traits.tp_a, Process_::chi_1p_a, 0, 0, 2, 4, 2, 2, 1);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_t_1(const double & q2) const
    {
        return _phi(q2, _traits.tp_v, Process_::chi_1m_t, 0, 0, 0, 2, 1, 3, 3);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_t_2(const double & q2) const
    {
        return _phi(q2, _traits.tp_a, Process_::chi_1p_t5, 1, 1, 0, 2, 1, 3, 1);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToV>::_phi_t_23(const double & q2) const
    {
        return _phi(q2, _traits.tp_a, Process_::chi_1p_t5, -1, 0, 2, 1, 0, 3, 1);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::_a_A12_0() const
    {
        const double x_A12 = _phi_a_12(0.0) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0)
                             * (power_of<2>(_mB) - power_of<2>(_mV)) / 8.0 / _mB / _mV;
        const double x_A0  = _phi_a_0(0.0) * (power_of<2>(_traits.m_R_0m) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_0m)) : 1.0);
        std::array<double, 5> a;
        a[0] = x_A12 * this->_a_A0[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_A12 * this->_a_A0[i] - x_A0 * this->_a_A12[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_a(_traits.calc_z(0.0, _traits.tp_a, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_A0);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::_a_T2_0() const
    {
        const double x_T2 = _phi_t_2(0.0) * (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(0.0, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double x_T1 = _phi_t_1(0.0) * (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(0.0, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const auto polynomials_T2 = _traits.orthonormal_polynomials_a(_traits.calc_z(0.0, _traits.tp_a, _traits.t0));
        const auto polynomials_T1 = _traits.orthonormal_polynomials_v(_traits.calc_z(0.0, _traits.tp_v, _traits.t0));

        double a_T2_0 = x_T2 * this->_a_T1[0] * polynomials_T1[0];
        for (unsigned i = 1 ; i < _a_T1.size() ; ++i)
        {
            a_T2_0 += x_T2 * this->_a_T1[i] * polynomials_T1[i] - x_T1 * this->_a_T2[i - 1] * polynomials_T2[i];
        }

        return a_T2_0 / (polynomials_T2[0] * x_T1);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::_a_A1_0() const
    {
        const double x_A1  = _phi_a_1( _traits.tm()) * 16.0 * _mB * _mV * _mV / (_mB + _mV) / (_mB * _mB - _mV * _mV - _traits.tm());
        const double x_A12 = _phi_a_12(_traits.tm());
        std::array<double, 5> a;
        a[0] = x_A1 * this->_a_A12_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_A1 * this->_a_A12[i - 1] - x_A12 * this->_a_A1[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_a(_traits.calc_z(_traits.tm(), _traits.tp_a, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_A12);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::_a_T23_0() const
    {
        const double x_T23 = _phi_t_23(_traits.tm()) * (_mB + _mV) * (_mB * _mB + 3 * _mV * _mV - _traits.tm()) / 8.0 / _mB / _mV / _mV;
        const double x_T2  = _phi_t_2( _traits.tm());
        std::array<double, 5> a;
        a[0] = x_T23 * this->_a_T2_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_T23 * this->_a_T2[i - 1] - x_T2 * this->_a_T23[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials_a(_traits.calc_z(_traits.tm(), _traits.tp_a, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_T2);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_V.begin(), _a_V.end(), coefficients.begin());
        // resonances for 1^m
        const double blaschke     = (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_v(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::a_0(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_A0.begin(), _a_A0.end(), coefficients.begin());
        // resonances for 0^m
        const double blaschke     = (power_of<2>(_traits.m_R_0m) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_0m)) : 1.0);
        const double phi          = _phi_a_0(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::a_1(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_A1_0();
        std::copy(_a_A1.begin(), _a_A1.end(), coefficients.begin() + 1);
        // resonances for 1^p
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_a_1(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::a_12(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_A12_0();
        std::copy(_a_A12.begin(), _a_A12.end(), coefficients.begin() + 1);
        // resonances for 1^p
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_a_12(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::t_1(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_T1.begin(), _a_T1.end(), coefficients.begin());
        // resonances for T (1^m state)
        const double blaschke     = (power_of<2>(_traits.m_R_1m) < _traits.tp_v ? _traits.calc_z(q2, _traits.tp_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_t_1(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_v, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::t_2(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_T2_0();
        std::copy(_a_T2.begin(), _a_T2.end(), coefficients.begin() + 1);
        // resonances for T5 (1^p state)
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_t_2(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::t_23(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_T23_0();
        std::copy(_a_T23.begin(), _a_T23.end(), coefficients.begin() + 1);
        // resonances for T (1^p state)
        const double blaschke     = (power_of<2>(_traits.m_R_1p) < _traits.tp_a ? _traits.calc_z(q2, _traits.tp_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const double phi          = _phi_t_23(q2);
        const double z            = _traits.calc_z(q2, _traits.tp_a, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::a_2(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return (power_of<2>(_mB + _mV) * (power_of<2>(_mB) - power_of<2>(_mV) - s) * a_1(s)
                - 16.0 * _mB * power_of<2>(_mV) * (_mB + _mV) * a_12(s)) / lambda;
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::t_3(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return ((power_of<2>(_mB) - power_of<2>(_mV)) * (power_of<2>(_mB) + 3.0 * power_of<2>(_mV) - s) * t_2(s)
                - 8.0 * _mB * power_of<2>(_mV) * (_mB - _mV) * t_23(s)) / lambda;
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::f_perp(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return pow(2*lambda, 0.5) / _mB / (_mB + _mV) * v(s);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::f_para(const double & s) const
    {
        return pow(2, 0.5) * (_mB + _mV) / _mB * a_1(s);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::f_long(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return ((power_of<2>(_mB) - power_of<2>(_mV) - s) * pow(_mB + _mV, 2) * a_1(s) - lambda * a_2(s))
                / (2 * _mV * power_of<2>(_mB) * (_mB + _mV));
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::f_perp_T(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return pow(2 * lambda, 0.5) / power_of<2>(_mB) * t_1(s);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::f_para_T(const double & s) const
    {
        return pow(2, 0.5) * (power_of<2>(_mB) - power_of<2>(_mV)) / power_of<2>(_mB) * t_2(s);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::f_long_T(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return s * (power_of<2>(_mB) + 3*power_of<2>(_mV) - s) / (2 * pow(_mB, 3) * _mV) * t_2(s)
                - s * lambda / (2 * pow(_mB, 3) * _mV * (power_of<2>(_mB) - power_of<2>(_mV))) * t_3(s);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::saturation_0p_v() const
    {
        return 0.0;
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::saturation_0m_a() const
    {
        return std::inner_product(_a_A0.begin(), _a_A0.end(), _a_A0.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::saturation_1m_v() const
    {
        return std::inner_product(_a_V.begin(), _a_V.end(), _a_V.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::saturation_1p_a() const
    {
        std::array<double, 5> coefficients_A12;
        coefficients_A12[0] = _a_A12_0();
        std::copy(_a_A12.begin(), _a_A12.end(), coefficients_A12.begin() + 1);

        std::array<double, 5> coefficients_A1;
        coefficients_A1[0] = _a_A1_0();
        std::copy(_a_A1.begin(), _a_A1.end(), coefficients_A1.begin() + 1);


        return std::inner_product(coefficients_A12.begin(), coefficients_A12.end(), coefficients_A12.begin(), 0.0)
          + std::inner_product(coefficients_A1.begin(), coefficients_A1.end(), coefficients_A1.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::saturation_1m_t() const
    {
        return std::inner_product(_a_T1.begin(), _a_T1.end(), _a_T1.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToV>::saturation_1p_t5() const
    {
        std::array<double, 5> coefficients_T2;
        coefficients_T2[0] = _a_T2_0();
        std::copy(_a_T2.begin(), _a_T2.end(), coefficients_T2.begin() + 1);

        std::array<double, 5> coefficients_T23;
        coefficients_T23[0] = _a_T23_0();
        std::copy(_a_T23.begin(), _a_T23.end(), coefficients_T23.begin() + 1);

        return std::inner_product(coefficients_T2.begin(), coefficients_T2.end(), coefficients_T2.begin(), 0.0)
          + std::inner_product(coefficients_T23.begin(), coefficients_T23.end(), coefficients_T23.begin(), 0.0);
    }

    template <typename Process_>
    Diagnostics
    BFW2010FormFactors<Process_, PToV>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _traits.calc_z(0.0,  _traits.tp_a, _traits.t0), "z_a(q2 =  0)" });
        results.add({ _traits.calc_z(0.0,  _traits.tp_v, _traits.t0), "z_v(q2 =  0)" });
        results.add({ _traits.calc_z(10.0, _traits.tp_a, _traits.t0), "z_a(q2 = 10)" });
        results.add({ _traits.calc_z(10.0, _traits.tp_v, _traits.t0), "z_v(q2 = 10)" });

        {
            const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials_v(0.0);
            results.add({ p0,              "p_0(z = 0.0)" });
            results.add({ p1,              "p_1(z = 0.0)" });
            results.add({ p2,              "p_2(z = 0.0)" });
            results.add({ p3,              "p_3(z = 0.0)" });
            results.add({ p4,              "p_4(z = 0.0)" });
            results.add({ p5,              "p_5(z = 0.0)" });
        }

        {
            const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials_v(_traits.calc_z(10.0, _traits.tp_v, _traits.t0));
            results.add({ p0,              "p_0(z = z(q2 = 10))" });
            results.add({ p1,              "p_1(z = z(q2 = 10))" });
            results.add({ p2,              "p_2(z = z(q2 = 10))" });
            results.add({ p3,              "p_3(z = z(q2 = 10))" });
            results.add({ p4,              "p_4(z = z(q2 = 10))" });
            results.add({ p5,              "p_5(z = z(q2 = 10))" });
        }

        {
            results.add({ _phi_v(-2.0),     "phi_v(z = z(q2 = -2.0))" });
            results.add({ _phi_v( 1.0),     "phi_v(z = z(q2 =  1.0))" });
            results.add({ _phi_v( 4.0),     "phi_v(z = z(q2 =  4.0))" });

            results.add({ _phi_a_0(-2.0),   "phi_a_0(z = z(q2 = -2.0))" });
            results.add({ _phi_a_0( 1.0),   "phi_a_0(z = z(q2 =  1.0))" });
            results.add({ _phi_a_0( 4.0),   "phi_a_0(z = z(q2 =  4.0))" });

            results.add({ _phi_a_1(-2.0),   "phi_a_1(z = z(q2 = -2.0))" });
            results.add({ _phi_a_1( 1.0),   "phi_a_1(z = z(q2 =  1.0))" });
            results.add({ _phi_a_1( 4.0),   "phi_a_1(z = z(q2 =  4.0))" });

            results.add({ _phi_a_12(-2.0),  "phi_a_12(z = z(q2 = -2.0))" });
            results.add({ _phi_a_12( 1.0),  "phi_a_12(z = z(q2 =  1.0))" });
            results.add({ _phi_a_12( 4.0),  "phi_a_12(z = z(q2 =  4.0))" });

            results.add({ _phi_t_1(-2.0),   "phi_t_1(z = z(q2 = -2.0))" });
            results.add({ _phi_t_1( 1.0),   "phi_t_1(z = z(q2 =  1.0))" });
            results.add({ _phi_t_1( 4.0),   "phi_t_1(z = z(q2 =  4.0))" });

            results.add({ _phi_t_2(-2.0),   "phi_t_2(z = z(q2 = -2.0))" });
            results.add({ _phi_t_2( 1.0),   "phi_t_2(z = z(q2 =  1.0))" });
            results.add({ _phi_t_2( 4.0),   "phi_t_2(z = z(q2 =  4.0))" });

            results.add({ _phi_t_23(-2.0),  "phi_t_23(z = z(q2 = -2.0))" });
            results.add({ _phi_t_23( 1.0),  "phi_t_23(z = z(q2 =  1.0))" });
            results.add({ _phi_t_23( 4.0),  "phi_t_23(z = z(q2 =  4.0))" });
        }

        {
            results.add({ _a_A1_0(),      "a_A1_0"  });
            results.add({ _a_A12_0(),     "a_A12_0" });
            results.add({ _a_T2_0(),      "a_T2_0"  });
            results.add({ _a_T23_0(),     "a_T23_0" });
        }

        return results;
    }


    template<typename Process_>
    const std::set<ReferenceName> BFW2010FormFactors<Process_, PToV>::references
    {
        "BFW:2010A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> BFW2010FormFactors<Process_, PToV>::options
    {
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BFW2010FormFactors<Process_, PToV>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BFW2010FormFactors<Process_, PToV>::end_options()
    {
        return options.cend();
    }


    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BFW2010FormFactorTraits<Process_, PToP>::resonance_0p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,0@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BFW2010FormFactorTraits<Process_, PToP>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template<typename Process_>
    BFW2010FormFactors<Process_, PToP>::BFW2010FormFactors(const Parameters & p, const Options &) :
        _a_fp{ UsedParameter(p[_par_name("f+", 0)], *this),
               UsedParameter(p[_par_name("f+", 1)], *this),
               UsedParameter(p[_par_name("f+", 2)], *this),
               UsedParameter(p[_par_name("f+", 3)], *this),
               UsedParameter(p[_par_name("f+", 4)], *this)
        },
        _a_ft{ UsedParameter(p[_par_name("fT", 0)], *this),
               UsedParameter(p[_par_name("fT", 1)], *this),
               UsedParameter(p[_par_name("fT", 2)], *this),
               UsedParameter(p[_par_name("fT", 3)], *this),
               UsedParameter(p[_par_name("fT", 4)], *this)
        },
        _a_f0{ UsedParameter(p[_par_name("f0", 1)], *this),
               UsedParameter(p[_par_name("f0", 2)], *this),
               UsedParameter(p[_par_name("f0", 3)], *this),
               UsedParameter(p[_par_name("f0", 4)], *this)
        },
        _traits(BFW2010FormFactorTraits<Process_, PToP>(p)),
        _mB(_traits.m_B),
        _mP(_traits.m_P)
    {
    }

    template<typename Process_>
    BFW2010FormFactors<Process_, PToP>::~BFW2010FormFactors() = default;

    template<typename Process_>
    FormFactors<PToP> *
    BFW2010FormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new BFW2010FormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    BFW2010FormFactors<Process_, PToP>::_par_name(const std::string & ff_name, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(idx) + "@BFW2010");
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::_phi(const double & t, const double & threshold_tp, const double & chi,
                                             const int & A, const unsigned B, const unsigned C, const unsigned k,
                                             const unsigned p, const unsigned n, const unsigned m) const
    {
        // [GRvDV:2022B]
        const double z = _traits.calc_z(t, threshold_tp, _traits.t0),
            kinematic_tp = power_of<2>(_mB + _mP);
        const double norm = std::sqrt(Process_::eta * k * pow(kinematic_tp, A) * pow(_traits.tm(), B)
                                * pow(4 * _mB * _mP, C) / 96 / M_PI / M_PI / chi);

        // set Q^2 to 0
        const double invt = 1 / ( 2.0 * (std::sqrt(threshold_tp) * std::sqrt(threshold_tp - t) + threshold_tp) - t); // simplification of -_z(t, t_p, 0) / t
        const double lambda_term = (kinematic_tp - t) * power_of<2>(std::sqrt(threshold_tp - t) + std::sqrt(threshold_tp - _traits.tm())); // simplification of lambda / z(t, threshold_tp, tm);
        const double sqrtjac = std::sqrt(4 * (1 + z) * (_traits.t0 - threshold_tp) / power_of<3>(z - 1)); // Abs[jacobian] = - jacobian

        return norm * sqrtjac * pow(lambda_term, 0.25 * m) * pow(invt, 0.5 * (p + n + 1.0));
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToP>::_phi_f_p(const double & q2) const
    {
        return _phi(q2, _traits.tp, Process_::chi_1m_v, 0, 0, 0, 1, 2, 2, 3);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToP>::_phi_f_0(const double & q2) const
    {
        return _phi(q2, _traits.tp, Process_::chi_0p_v, 1, 1, 0, 3, 2, 1, 1);
    }

    template<typename Process_>
    inline double
    BFW2010FormFactors<Process_, PToP>::_phi_f_t(const double & q2) const
    {
        return _phi(q2, _traits.tp, Process_::chi_1m_t, -1, 0, 0, 1, 0, 3, 3);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::_a_f0_0() const
    {
        const double x_f0 = _phi_f_0(0.0) * (power_of<2>(_traits.m_R_0p) < _traits.tp ? _traits.calc_z(0.0, _traits.tp, power_of<2>(_traits.m_R_0p)) : 1.0);
        const double x_fp = _phi_f_p(0.0) * (power_of<2>(_traits.m_R_1m) < _traits.tp ? _traits.calc_z(0.0, _traits.tp, power_of<2>(_traits.m_R_1m)) : 1.0);
        std::array<double, 5> a;
        a[0] = x_f0 * this->_a_fp[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_f0 * this->_a_fp[i] - x_fp * this->_a_f0[i - 1];
        }
        const auto polynomials = _traits.orthonormal_polynomials(_traits.calc_z(0.0, _traits.tp, _traits.t0));
        return std::inner_product(a.begin(), a.end(), polynomials.begin(), 0.0) / (polynomials[0] * x_fp);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_p(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        // resonances for 1^m
        const double blaschke = (power_of<2>(_traits.m_R_1m) < _traits.tp ? _traits.calc_z(q2, _traits.tp, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_f_p(q2);
        const double z            = _traits.calc_z(q2, _traits.tp, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_0(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_f0_0();
        std::copy(_a_f0.begin(), _a_f0.end(), coefficients.begin() + 1);
        // resonances for 0^p
        const double blaschke = (power_of<2>(_traits.m_R_0p) < _traits.tp ? _traits.calc_z(q2, _traits.tp, power_of<2>(_traits.m_R_0p)) : 1.0);
        const double phi          = _phi_f_0(q2);
        const double z            = _traits.calc_z(q2, _traits.tp, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_ft.begin(), _a_ft.end(), coefficients.begin());
        // resonances for 1^m
        const double blaschke = (power_of<2>(_traits.m_R_1m) < _traits.tp ? _traits.calc_z(q2, _traits.tp, power_of<2>(_traits.m_R_1m)) : 1.0);
        const double phi          = _phi_f_t(q2);
        const double z            = _traits.calc_z(q2, _traits.tp, _traits.t0);
        const auto   polynomials  = _traits.orthonormal_polynomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_p_series(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.tp, 0.0), complex<double>(_traits.t0, 0.0));
        const auto polynomials       = _traits.orthonormal_polynomials(z);
        const complex<double> series = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), complex<double>(0.0, 0.0));

        return abs(series);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_p_series_prime(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.tp, 0.0), complex<double>(_traits.t0, 0.0));
        const auto polynomials_prime = _traits.orthonormal_polynomials_derivatives(z);
        const complex<double> series_prime = std::inner_product(coefficients.begin(), coefficients.end(), polynomials_prime.begin(), complex<double>(0.0, 0.0));

        return abs(series_prime);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_0_series(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.tp, 0.0), complex<double>(_traits.t0, 0.0));
        const auto polynomials       = _traits.orthonormal_polynomials(z);
        const complex<double> series = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), complex<double>(0.0, 0.0));

        return abs(series);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_0_series_prime(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.tp, 0.0), complex<double>(_traits.t0, 0.0));
        const auto polynomials_prime = _traits.orthonormal_polynomials_derivatives(z);
        const complex<double> series_prime = std::inner_product(coefficients.begin(), coefficients.end(), polynomials_prime.begin(), complex<double>(0.0, 0.0));

        return abs(series_prime);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_t_series(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.tp, 0.0), complex<double>(_traits.t0, 0.0));
        const auto polynomials       = _traits.orthonormal_polynomials(z);
        const complex<double> series = std::inner_product(coefficients.begin(), coefficients.end(), polynomials.begin(), complex<double>(0.0, 0.0));

        return abs(series);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_t_series_prime(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.tp, 0.0), complex<double>(_traits.t0, 0.0));
        const auto polynomials_prime = _traits.orthonormal_polynomials_derivatives(z);
        const complex<double> series_prime = std::inner_product(coefficients.begin(), coefficients.end(), polynomials_prime.begin(), complex<double>(0.0, 0.0));

        return abs(series_prime);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::f_plus_T(const double & q2) const
    {
        return f_t(q2) * q2 / _mB / (_mB + _mP);
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::saturation_0p_v() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_f0_0();
        std::copy(_a_f0.begin(), _a_f0.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::saturation_0m_a() const
    {
        return 0.;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::saturation_1m_v() const
    {
        return std::inner_product(_a_fp.begin(), _a_fp.end(), _a_fp.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::saturation_1p_a() const
    {
        return 0.;
    }

    template <typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::saturation_1m_t() const
    {
        return std::inner_product(_a_ft.begin(), _a_ft.end(), _a_ft.begin(), 0.0);
    }

    template<typename Process_>
    double
    BFW2010FormFactors<Process_, PToP>::saturation_1p_t5() const
    {
        return 0.;
    }

    template <typename Process_>
    Diagnostics
    BFW2010FormFactors<Process_, PToP>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _traits.calc_z(0.0,  _traits.tp, _traits.t0), "z(q2 =  0)" });
        results.add({ _traits.calc_z(10.0, _traits.tp, _traits.t0), "z(q2 = 10)" });

        {
            const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials(0.0);
            results.add({ p0,              "p_0(z = 0.0)" });
            results.add({ p1,              "p_1(z = 0.0)" });
        }

        {
            const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials(_traits.calc_z(10.0, _traits.tp,_traits.t0));
            results.add({ p0,              "p_0(z = z(q2 = 10))" });
            results.add({ p1,              "p_1(z = z(q2 = 10))" });
        }

        {
            results.add({ _phi_f_p(-2.0),   "phi_f_p(z = z(q2 = -2))" });
            results.add({ _phi_f_p( 1.0),   "phi_f_p(z = z(q2 =  1))" });
            results.add({ _phi_f_p( 4.0),   "phi_f_p(z = z(q2 =  4))" });

            results.add({ _phi_f_0(-2.0),   "phi_f_0(z = z(q2 = -2))" });
            results.add({ _phi_f_0( 1.0),   "phi_f_0(z = z(q2 =  1))" });
            results.add({ _phi_f_0( 4.0),   "phi_f_0(z = z(q2 =  4))" });

            results.add({ _phi_f_t(-2.0),   "phi_f_t(z = z(q2 = -2))" });
            results.add({ _phi_f_t( 1.0),   "phi_f_t(z = z(q2 =  1))" });
            results.add({ _phi_f_t( 4.0),   "phi_f_t(z = z(q2 =  4))" });
        }

        {
            results.add({ _a_f0_0(),  "_a_f0_0" });
        }

        return results;
    }

    template<typename Process_>
    const std::set<ReferenceName> BFW2010FormFactors<Process_, PToP>::references
    {
        "BFW:2010A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> BFW2010FormFactors<Process_, PToP>::options
    {
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BFW2010FormFactors<Process_, PToP>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BFW2010FormFactors<Process_, PToP>::end_options()
    {
        return options.cend();
    }
}

#endif
