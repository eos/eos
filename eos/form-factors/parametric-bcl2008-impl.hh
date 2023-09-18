/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BCL2008_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BCL2008_IMPL_HH 1

#include <eos/form-factors/parametric-bcl2008.hh>
#include <eos/utils/exception.hh>

namespace eos
{
    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 3u, false>::_z(const double & s) const
    {
        static const double m_B = Process_::m_B;
        static const double m_P = Process_::m_P;
        static const double tau_p = (m_B + m_P) * (m_B + m_P);
        static const double tau_0 = (m_B + m_P) * (std::sqrt(m_B) - std::sqrt(m_P)) * (std::sqrt(m_B) - std::sqrt(m_P));

        return (std::sqrt(tau_p - s) - std::sqrt(tau_p - tau_0))
            / (std::sqrt(tau_p - s) + std::sqrt(tau_p - tau_0));
    }

    template <typename Process_>
    BCL2008FormFactorBase<Process_, 3u, false>::BCL2008FormFactorBase(const Parameters & p, const Options &) :
        _f_plus_0(p[std::string(Process_::label) + "::f_+(0)@BCL2008"], *this),
        _b_plus_1(p[std::string(Process_::label) + "::b_+^1@BCL2008"],  *this),
        _b_plus_2(p[std::string(Process_::label) + "::b_+^2@BCL2008"],  *this),
        _b_zero_1(p[std::string(Process_::label) + "::b_0^1@BCL2008"],  *this),
        _b_zero_2(p[std::string(Process_::label) + "::b_0^2@BCL2008"],  *this),
        _b_zero_3(p[std::string(Process_::label) + "::b_0^3@BCL2008"],  *this)
    {
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 3u, false>::f_p(const double & s) const
    {
        const double z = _z(s), z2 = z * z, z3 = z * z2;
        const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

        return _f_plus_0 / (1.0 - s / Process_::mR2_1m) * (1.0 + _b_plus_1 * (zbar - z3bar / 3.0) + _b_plus_2 * (z2bar + 2.0 * z3bar / 3.0));
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 3u, false>::f_0(const double & s) const
    {
        const double z = _z(s), z2 = z * z, z3 = z * z2;
        const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

        // note that f_0(0) = f_+(0)!
        // for f_0(s) we do not have an equation of motion to express _b_zero_K in terms of the
        // other coefficients!
        return _f_plus_0 / (1.0 - s / Process_::mR2_0p) * (1.0 + _b_zero_1 * zbar + _b_zero_2 * z2bar + _b_zero_3 * z3bar);
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 3u, false>::f_t(const double &) const
    {
        throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

        return 0.0;
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 3u, false>::f_plus_T(const double &) const
    {
        throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

        return 0.0;
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 4u, false>::_z(const double & s) const
    {
        static const double m_B = Process_::m_B;
        static const double m_P = Process_::m_P;
        static const double tau_p = (m_B + m_P) * (m_B + m_P);
        static const double tau_0 = (m_B + m_P) * (std::sqrt(m_B) - std::sqrt(m_P)) * (std::sqrt(m_B) - std::sqrt(m_P));

        return (std::sqrt(tau_p - s) - std::sqrt(tau_p - tau_0))
            / (std::sqrt(tau_p - s) + std::sqrt(tau_p - tau_0));
    }

    template <typename Process_>
    BCL2008FormFactorBase<Process_, 4u, false>::BCL2008FormFactorBase(const Parameters & p, const Options &) :
        _f_plus_0(p[std::string(Process_::label) + "::f_+(0)@BCL2008"], *this),
        _b_plus_1(p[std::string(Process_::label) + "::b_+^1@BCL2008"],  *this),
        _b_plus_2(p[std::string(Process_::label) + "::b_+^2@BCL2008"],  *this),
        _b_plus_3(p[std::string(Process_::label) + "::b_+^3@BCL2008"],  *this),
        _b_zero_1(p[std::string(Process_::label) + "::b_0^1@BCL2008"],  *this),
        _b_zero_2(p[std::string(Process_::label) + "::b_0^2@BCL2008"],  *this),
        _b_zero_3(p[std::string(Process_::label) + "::b_0^3@BCL2008"],  *this),
        _b_zero_4(p[std::string(Process_::label) + "::b_0^4@BCL2008"],  *this)
    {
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 4u, false>::f_p(const double & s) const
    {
        const double z = _z(s), z2 = z * z, z3 = z * z2, z4 = z * z3;
        const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02, z04 = z0 * z03;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03, z4bar = z4 - z04;

        return _f_plus_0 / (1.0 - s / Process_::mR2_1m) * (1.0 + _b_plus_1 * (zbar + z4bar / 4.0) + _b_plus_2 * (z2bar - z4bar / 2.0) + _b_plus_3 * (z3bar + 3.0 * z4bar / 4.0));
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 4u, false>::f_0(const double & s) const
    {
        const double z = _z(s), z2 = z * z, z3 = z * z2, z4 = z * z3;
        const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02, z04 = z0 * z03;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03, z4bar = z4 - z04;

        // note that f_0(0) = f_+(0)!
        // for f_0(s) we do not have an equation of motion to express _b_zero_K in terms of the
        // other coefficients!
        return _f_plus_0 / (1.0 - s / Process_::mR2_0p) * (1.0 + _b_zero_1 * zbar + _b_zero_2 * z2bar + _b_zero_3 * z3bar + _b_zero_4 * z4bar);
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 4u, false>::f_t(const double &) const
    {
        throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

        return 0.0;
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 4u, false>::f_plus_T(const double &) const
    {
        throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

        return 0.0;
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 5u, false>::_z(const double & s) const
    {
        static const double m_B = Process_::m_B;
        static const double m_P = Process_::m_P;
        static const double tau_p = (m_B + m_P) * (m_B + m_P);
        static const double tau_0 = (m_B + m_P) * (std::sqrt(m_B) - std::sqrt(m_P)) * (std::sqrt(m_B) - std::sqrt(m_P));

        return (std::sqrt(tau_p - s) - std::sqrt(tau_p - tau_0))
            / (std::sqrt(tau_p - s) + std::sqrt(tau_p - tau_0));
    }

    template <typename Process_>
    BCL2008FormFactorBase<Process_, 5u, false>::BCL2008FormFactorBase(const Parameters & p, const Options &) :
        _f_plus_0(p[std::string(Process_::label) + "::f_+(0)@BCL2008"], *this),
        _b_plus_1(p[std::string(Process_::label) + "::b_+^1@BCL2008"],  *this),
        _b_plus_2(p[std::string(Process_::label) + "::b_+^2@BCL2008"],  *this),
        _b_plus_3(p[std::string(Process_::label) + "::b_+^3@BCL2008"],  *this),
        _b_plus_4(p[std::string(Process_::label) + "::b_+^4@BCL2008"],  *this),
        _b_zero_1(p[std::string(Process_::label) + "::b_0^1@BCL2008"],  *this),
        _b_zero_2(p[std::string(Process_::label) + "::b_0^2@BCL2008"],  *this),
        _b_zero_3(p[std::string(Process_::label) + "::b_0^3@BCL2008"],  *this),
        _b_zero_4(p[std::string(Process_::label) + "::b_0^4@BCL2008"],  *this),
        _b_zero_5(p[std::string(Process_::label) + "::b_0^5@BCL2008"],  *this)

    {
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 5u, false>::f_p(const double & s) const
    {
        const double z = _z(s), z2 = z * z, z3 = z * z2, z4 = z * z3, z5 = z * z4;
        const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02, z04 = z0 * z03, z05 = z0 * z04;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03, z4bar = z4 - z04, z5bar = z5 - z05;

        return _f_plus_0 / (1.0 - s / Process_::mR2_1m) * (1.0 + _b_plus_1 * (zbar - z5bar / 5.0) + _b_plus_2 * (z2bar + 2.0 * z5bar / 5.0) + _b_plus_3 * (z3bar - 3.0 * z5bar / 5.0) + _b_plus_4 * (z4bar + 4.0 * z5bar / 5.0));
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 5u, false>::f_0(const double & s) const
    {
        const double z = _z(s), z2 = z * z, z3 = z * z2, z4 = z * z3, z5 = z * z4;
        const double z0 = _z(0), z02 = z0 * z0, z03 = z0 * z02, z04 = z0 * z03, z05 = z0 * z04;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03, z4bar = z4 - z04, z5bar = z5 - z05;

        // note that f_0(0) = f_+(0)!
        // for f_0(s) we do not have an equation of motion to express _b_zero_K in terms of the
        // other coefficients!
        return _f_plus_0 / (1.0 - s / Process_::mR2_0p) * (1.0 + _b_zero_1 * zbar + _b_zero_2 * z2bar + _b_zero_3 * z3bar + _b_zero_4 * z4bar + _b_zero_5 * z5bar);
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 5u, false>::f_t(const double &) const
    {
        throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

        return 0.0;
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 5u, false>::f_plus_T(const double &) const
    {
        throw InternalError("This form factor parametrization has no inputs for tensor form factors.");

        return 0.0;
    }

    template <typename Process_>
    BCL2008FormFactorBase<Process_, 3u, true>::BCL2008FormFactorBase(const Parameters & p, const Options & o) :
        BCL2008FormFactorBase<Process_, 3u, false>(p, o),
        _f_t_0(p[std::string(Process_::label)    + "::f_T(0)@BCL2008"], *this),
        _b_t_1(p[std::string(Process_::label)    + "::b_T^1@BCL2008"],  *this),
        _b_t_2(p[std::string(Process_::label)    + "::b_T^2@BCL2008"],  *this)
    {
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 3u, true>::f_t(const double & s) const
    {
        const double z = this->_z(s), z2 = z * z, z3 = z * z2;
        const double z0 = this->_z(0), z02 = z0 * z0, z03 = z0 * z02;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03;

        return _f_t_0 / (1.0 - s / Process_::mR2_1m) * (1.0 + _b_t_1 * (zbar - z3bar / 3.0) + _b_t_2 * (z2bar + 2.0 * z3bar / 3.0));
    }

    template <typename Process_>
    BCL2008FormFactorBase<Process_, 4u, true>::BCL2008FormFactorBase(const Parameters & p, const Options & o) :
        BCL2008FormFactorBase<Process_, 4u, false>(p, o),
        _f_t_0(p[std::string(Process_::label)    + "::f_T(0)@BCL2008"], *this),
        _b_t_1(p[std::string(Process_::label)    + "::b_T^1@BCL2008"],  *this),
        _b_t_2(p[std::string(Process_::label)    + "::b_T^2@BCL2008"],  *this),
        _b_t_3(p[std::string(Process_::label)    + "::b_T^3@BCL2008"],  *this)

    {
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 4u, true>::f_t(const double & s) const
    {
        const double z = this->_z(s), z2 = z * z, z3 = z * z2, z4 = z * z3;
        const double z0 = this->_z(0), z02 = z0 * z0, z03 = z0 * z02, z04 = z0 * z03;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03, z4bar = z4 - z04;

        return _f_t_0 / (1.0 - s / Process_::mR2_1m) * (1.0 + _b_t_1 * (zbar + z4bar / 4.0) + _b_t_2 * (z2bar - z4bar / 2.0) + _b_t_3 * (z3bar + 3.0 * z4bar / 4.0));
    }

    template <typename Process_>
    BCL2008FormFactorBase<Process_, 5u, true>::BCL2008FormFactorBase(const Parameters & p, const Options & o) :
        BCL2008FormFactorBase<Process_, 5u, false>(p, o),
        _f_t_0(p[std::string(Process_::label)    + "::f_T(0)@BCL2008"], *this),
        _b_t_1(p[std::string(Process_::label)    + "::b_T^1@BCL2008"],  *this),
        _b_t_2(p[std::string(Process_::label)    + "::b_T^2@BCL2008"],  *this),
        _b_t_3(p[std::string(Process_::label)    + "::b_T^3@BCL2008"],  *this),
        _b_t_4(p[std::string(Process_::label)    + "::b_T^4@BCL2008"],  *this)

    {
    }

    template <typename Process_>
    double
    BCL2008FormFactorBase<Process_, 5u, true>::f_t(const double & s) const
    {
        const double z = this->_z(s), z2 = z * z, z3 = z * z2, z4 = z * z3, z5 = z * z4;
        const double z0 = this->_z(0), z02 = z0 * z0, z03 = z0 * z02, z04 = z0 * z03, z05 = z0 * z04;
        const double zbar = z - z0, z2bar = z2 - z02, z3bar = z3 - z03, z4bar = z4 - z04, z5bar = z5 - z05;

        return _f_t_0 / (1.0 - s / Process_::mR2_1m) * (1.0 + _b_t_1 * (zbar - z5bar / 5.0) + _b_t_2 * (z2bar + 2.0 * z5bar / 5.0) + _b_t_3 * (z3bar - 3.0 * z5bar / 5.0) + _b_t_4 * (z4bar + 4.0 * z5bar / 5.0));
    }

    template <typename Process_, unsigned K_>
    BCL2008FormFactors<Process_, K_>::BCL2008FormFactors(const Parameters & p, const Options & o) :
        BCL2008FormFactorBase<Process_, K_, Process_::uses_tensor_form_factors>(p, o)
    {
    }

    template <typename Process_, unsigned K_>
    FormFactors<PToP> *
    BCL2008FormFactors<Process_, K_>::make(const Parameters & parameters, const Options & options)
    {
        return new BCL2008FormFactors<Process_, K_>(parameters, options);
    }
}

#endif
