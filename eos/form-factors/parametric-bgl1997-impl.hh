/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
 * Copyright (c) 2020 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_IMPL_HH 1

#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/utils/kinematic.hh>
#include <eos/models/model.hh>
#include <eos/utils/power_of.hh>

#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    BGL1997FormFactorBase::BGL1997FormFactorBase(const Parameters &, const Options &, ParameterUser &, const double t_p, const double t_m) :
        _t_p(t_p),
        _t_m(t_m),
        _chi_1m( 5.131e-04), // TODO remove hard-coded numerical values
        _chi_0p( 6.204e-03),
        _chi_1p( 3.894e-04),
        _chi_0m(19.421e-03)
    {
    }

    BGL1997FormFactorBase::~BGL1997FormFactorBase() = default;

    double
    BGL1997FormFactorBase::_z(const double & t, const double & t_0) const
    {
        return (std::sqrt(_t_p - t) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - t) + std::sqrt(_t_p - t_0));
    }

    double
    BGL1997FormFactorBase::_phi(const double & s, const double & t_0, const unsigned & K, const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const
    {
        const double sq_tp    = std::sqrt(_t_p);
        const double sq_tp_t  = std::sqrt(_t_p - s);
        const double sq_tp_t0 = std::sqrt(_t_p - t_0);
        const double sq_tp_tm = std::sqrt(_t_p - _t_m);

        // [BGL:1997A] eq. (4.14) for OPE at Q^2 = -q^2 = 0
        // => generalization for q^2 != 0 possible, see eq.(4.15)
        return std::sqrt(1.0 / (K * M_PI * chi)) * (sq_tp_t + sq_tp_t0)
               * std::pow(sq_tp_t / sq_tp_t0, 1.0 / 2.0)
               * std::pow(_t_p - s, a / 4.0)
               * std::pow(sq_tp_t + sq_tp_tm, b / 2.0)
               * 1.0 / std::pow(sq_tp_t + sq_tp, c + 3.0);
    }


    // TODO hard-coded values of resonances from [BGS2017] table III with some changes (from Nico Gubernari)
    //  0^+ at 6.704, 7.122
    //  0^- at 6.275, 6.871, 7.250
    //  1^+ at 6.739, 6.750, 7.145, 7.150
    //  1^- at 6.329, 6.910, 7.020

    // TODO look on constraints on FF parametrization coming from kinematic constraints at q^2 = 0
    // 2 mV A_0 = (mB + mV) A_1 - (mB - mV) A_2
    // T_2 = T_3
    // f_+ = f_0
    std::string
    BGL1997FormFactors<BToDstar>::_par_name(const std::string & ff_name)
    {
        return std::string("B->D^*") + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    BGL1997FormFactors<BToDstar>::BGL1997FormFactors(const Parameters & p, const Options & o) :
        BGL1997FormFactorBase(p, o, *this, power_of<2>(BToDstar::mB + BToDstar::mV), power_of<2>(BToDstar::mB - BToDstar::mV)),
        _a_g{{   UsedParameter(p[_par_name("g_0")],  *this),
                 UsedParameter(p[_par_name("g_1")],  *this),
                 UsedParameter(p[_par_name("g_2")],  *this),
                 UsedParameter(p[_par_name("g_3")],  *this) }},
        _a_f{{   UsedParameter(p[_par_name("f_0")],  *this),
                 UsedParameter(p[_par_name("f_1")],  *this),
                 UsedParameter(p[_par_name("f_2")],  *this),
                 UsedParameter(p[_par_name("f_3")],  *this) }},
        _a_F1{{  UsedParameter(p[_par_name("F1_0")], *this),
                 UsedParameter(p[_par_name("F1_1")], *this),
                 UsedParameter(p[_par_name("F1_2")], *this),
                 UsedParameter(p[_par_name("F1_3")], *this) }},
        _a_F2{{  UsedParameter(p[_par_name("F2_0")], *this),
                 UsedParameter(p[_par_name("F2_1")], *this),
                 UsedParameter(p[_par_name("F2_2")], *this),
                 UsedParameter(p[_par_name("F2_3")], *this) }},
        _mB(BToDstar::mB),
        _mB2(power_of<2>(_mB)),
        _mV(BToDstar::mV),
        _mV2(power_of<2>(_mV)),
        // here optimal t_0 = sqrt(t_p) (sqrt(m_B) - sqrt(m_M))^2
        _t_0((_mB + _mV)* power_of<2>(std::sqrt(_mB) - std::sqrt(_mV)))
    {
    }

    BGL1997FormFactors<BToDstar>::~BGL1997FormFactors() = default;

    FormFactors<PToV> *
    BGL1997FormFactors<BToDstar>::make(const Parameters & parameters, const Options & options)
    {
        return new BGL1997FormFactors(parameters, options);
    }

    double
    BGL1997FormFactors<BToDstar>::g(const double & s) const
    {
        // resonances for 1^-
        const double blaschke = _z(s, 6.329 * 6.329) * _z(s, 6.910 * 6.910) * _z(s, 7.020 * 7.020);
        const double phi      = _phi(s, _t_0, 96, 3, 3, 1, _chi_1m);
        const double z        = _z(s, _t_0);
        const double series   = _a_g[0] + _a_g[1] * z + _a_g[2] * z * z + _a_g[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::f(const double & s) const
    {
        // resonances for 1^+
        const double blaschke = _z(s, 6.739 * 6.739) * _z(s, 6.750 * 6.750) * _z(s, 7.145 * 7.145) * _z(s, 7.150 * 7.150);
        const double phi      = _phi(s, _t_0, 24, 1, 1, 1, _chi_1p);
        const double z        = _z(s, _t_0);
        const double series   = _a_f[0] + _a_f[1] * z + _a_f[2] * z * z + _a_f[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::F1(const double & s) const
    {
        // resonances for 1^+
        const double blaschke = _z(s, 6.739 * 6.739) * _z(s, 6.750 * 6.750) * _z(s, 7.145 * 7.145) * _z(s, 7.150 * 7.150);
        const double phi      = _phi(s, _t_0, 48, 1, 1, 2, _chi_1p);
        const double z        = _z(s, _t_0);
        const double series   = _a_F1[0] + _a_F1[1] * z + _a_F1[2] * z * z + _a_F1[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::F2(const double & s) const
    {
        // resonances for 0^-
        const double blaschke = _z(s, 6.275 * 6.275) * _z(s, 6.871 * 6.871) * _z(s, 7.250 * 7.250);
        const double phi      = _phi(s, _t_0, 64, 3, 3, 1, _chi_0m);
        const double z        = _z(s, _t_0);
        const double series   = _a_F2[0] + _a_F2[1] * z + _a_F2[2] * z * z + _a_F2[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::v(const double & s) const
    {
        return (_mB + _mV) / 2.0 * g(s);
    }

    double
    BGL1997FormFactors<BToDstar>::a_0(const double & s) const
    {
        return F2(s) / 2.0;
    }

    double
    BGL1997FormFactors<BToDstar>::a_1(const double & s) const
    {
        return 1.0/(_mB + _mV) * f(s);
    }

    double
    BGL1997FormFactors<BToDstar>::a_2(const double & s) const
    {
        return (_mB + _mV) / eos::lambda(_mB2, _mV2, s) * ((_mB2 - _mV2 - s) * f(s) - 2.0 * _mV * F1(s));
    }

    double
    BGL1997FormFactors<BToDstar>::a_12(const double & s) const
    {
        return F1(s) / (8.0 * _mB * _mV);
    }

    double
    BGL1997FormFactors<BToDstar>::t_1(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::t_2(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::t_3(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::t_23(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_perp(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_para(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_long(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_perp_T(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_para_T(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_long_T(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::f_long_T_Normalized(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }



    std::string
    BGL1997FormFactors<BToD>::_par_name(const std::string & ff_name)
    {
        return std::string("B->D") + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    BGL1997FormFactors<BToD>::BGL1997FormFactors(const Parameters & p, const Options & o) :
        BGL1997FormFactorBase(p, o, *this, power_of<2>(BToD::m_B + BToD::m_P), power_of<2>(BToD::m_B - BToD::m_P)),
        _a_f_p{{ UsedParameter(p[_par_name("f+_0")], *this),
                 UsedParameter(p[_par_name("f+_1")], *this),
                 UsedParameter(p[_par_name("f+_2")], *this),
                 UsedParameter(p[_par_name("f+_3")], *this) }},
        _a_f_0{{ UsedParameter(p[_par_name("f0_0")], *this),
                 UsedParameter(p[_par_name("f0_1")], *this),
                 UsedParameter(p[_par_name("f0_2")], *this),
                 UsedParameter(p[_par_name("f0_3")], *this) }},
        _a_f_t{{ UsedParameter(p[_par_name("fT_0")], *this),
                 UsedParameter(p[_par_name("fT_1")], *this),
                 UsedParameter(p[_par_name("fT_2")], *this),
                 UsedParameter(p[_par_name("fT_3")], *this) }},
        _mB(BToD::m_B),
        _mB2(power_of<2>(_mB)),
        _mP(BToD::m_P),
        _mP2(power_of<2>(_mP)),
        // here optimal t_0 = sqrt(t_p) (sqrt(m_B) - sqrt(m_M))^2
        _t_0((_mB + _mP)* power_of<2>(std::sqrt(_mB) - std::sqrt(_mP)))
    {
    }

    BGL1997FormFactors<BToD>::~BGL1997FormFactors() = default;

    FormFactors<PToP> *
    BGL1997FormFactors<BToD>::make(const Parameters & parameters, const Options & options)
    {
        return new BGL1997FormFactors(parameters, options);
    }

    double
    BGL1997FormFactors<BToD>::f_p(const double & s) const
    {
        // resonances for 1^-
        const double blaschke = _z(s, 6.329 * 6.329) * _z(s, 6.910 * 6.910) * _z(s, 7.020 * 7.020);
        const double phi      = _phi(s, _t_0, 48, 3, 3, 2, _chi_1m);
        const double z        = _z(s, _t_0);
        const double series   = _a_f_p[0] + _a_f_p[1] * z + _a_f_p[2] * z * z + _a_f_p[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToD>::f_0(const double & s) const
    {
        // resonances for 0^+
        const double blaschke = _z(s, 6.704 * 6.704) * _z(s, 7.122 * 7.122);
        const double phi      = _phi(s, _t_0, 16, 1, 1, 1, _chi_0p);
        const double z        = _z(s, _t_0);
        const double series   = _a_f_0[0] + _a_f_0[1] * z + _a_f_0[2] * z * z + _a_f_0[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToD>::f_t(const double & /*s*/) const
    {
        return 0.0; //  TODO
    }

    double
    BGL1997FormFactors<BToD>::f_plus_T(const double & /*s*/) const
    {
        return 0.0; //  TODO
    }
}

#endif
