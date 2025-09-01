/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KMPW2010_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KMPW2010_IMPL_HH 1

#include <eos/form-factors/parametric-kmpw2010.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{
    double
    KMPW2010FormFactors<PToV>::_calc_z(const double & s)
    {
        return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
    }

    double
    KMPW2010FormFactors<PToV>::ff_KMPW(const double & s, const double & f0, const double & b1, const double & m2)
    {
        const double zs = _calc_z(s), z0 = _calc_z(0.0);

        // cf. [KMPW:2010A], Eq. (8.8), p. 30
        return f0 / (1.0 - s / m2) * (1.0 + b1 * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
    }

    KMPW2010FormFactors<PToV>::KMPW2010FormFactors(const Parameters & p, const Options &) :
        _f0_V(p["B->K^*::F^V(0)@KMPW2010"],   *this),   _b1_V(p["B->K^*::b^V_1@KMPW2010"],    *this),
        _f0_A0(p["B->K^*::F^A0(0)@KMPW2010"], *this),   _b1_A0(p["B->K^*::b^A0_1@KMPW2010"],  *this),
        _f0_A1(p["B->K^*::F^A1(0)@KMPW2010"], *this),   _b1_A1(p["B->K^*::b^A1_1@KMPW2010"],  *this),
        _f0_A2(p["B->K^*::F^A2(0)@KMPW2010"], *this),   _b1_A2(p["B->K^*::b^A2_1@KMPW2010"],  *this),
        _f0_T1(p["B->K^*::F^T1(0)@KMPW2010"], *this),   _b1_T1(p["B->K^*::b^T1_1@KMPW2010"],  *this),
        _f0_T2(p["B->K^*::F^T2(0)@KMPW2010"], *this),   _b1_T2(p["B->K^*::b^T2_1@KMPW2010"],  *this),
        _f0_T3(p["B->K^*::F^T3(0)@KMPW2010"], *this),   _b1_T3(p["B->K^*::b^T3_1@KMPW2010"],  *this)
    {
    }

    FormFactors<PToV> *
    KMPW2010FormFactors<PToV>::make(const Parameters & parameters, const Options & options)
    {
        return new KMPW2010FormFactors(parameters, options);
    }

    double
    KMPW2010FormFactors<PToV>::v(const double & s) const
    {
        return ff_KMPW(s, _f0_V(), _b1_V(), _m_Bs2_1m);
    }

    double
    KMPW2010FormFactors<PToV>::a_0(const double & s) const
    {
        return ff_KMPW(s, _f0_A0(), _b1_A0(), _m_Bs2_0m);
    }

    double
    KMPW2010FormFactors<PToV>::a_1(const double & s) const
    {
        return ff_KMPW(s, _f0_A1(), _b1_A1(), _m_Bs2_1p);
    }

    double
    KMPW2010FormFactors<PToV>::a_2(const double & s) const
    {
        return ff_KMPW(s, _f0_A2(), _b1_A2(), _m_Bs2_1p);
    }

    double
    KMPW2010FormFactors<PToV>::a_12(const double & s) const
    {
        const double mB = BToKstar::m_B, mB2 = mB * mB;
        const double mV = BToKstar::m_V, mV2 = mV * mV;
        const double lambda = eos::lambda(mB2, mV2, s);

        return ((mB + mV) * (mB + mV) * (mB2 - mV2 - s) * this->a_1(s)
            - lambda * this->a_2(s)) / (16.0 * mB * mV2 * (mB + mV));
    }

    double
    KMPW2010FormFactors<PToV>::t_1(const double & s) const
    {
        return ff_KMPW(s, _f0_T1(), _b1_T1(), _m_Bs2_1m);
    }

    double
    KMPW2010FormFactors<PToV>::t_2(const double & s) const
    {
        return ff_KMPW(s, _f0_T2(), _b1_T2(), _m_Bs2_1p);
    }

    double
    KMPW2010FormFactors<PToV>::t_3(const double & s) const
    {
        return ff_KMPW(s, _f0_T3(), _b1_T3(), _m_Bs2_1p);
    }

    double
    KMPW2010FormFactors<PToV>::t_23(const double & s) const
    {
        const double mB = BToKstar::m_B, mB2 = mB * mB;
        const double mV = BToKstar::m_V, mV2 = mV * mV;
        const double lambda = eos::lambda(mB2, mV2, s);

        return ((mB2 - mV2) * (mB2 + 3.0 * mV2 - s) * this->t_2(s)
                - lambda * this->t_3(s)) / (8.0 * mB * mV2 * (mB - mV));
    }

    double
    KMPW2010FormFactors<PToV>::f_perp(const double &) const
    {
        return 0.;
    }

    double
    KMPW2010FormFactors<PToV>::f_para(const double &) const
    {
        return 0.;
    }

    double
    KMPW2010FormFactors<PToV>::f_long(const double &) const
    {
        return 0.;
    }

    double
    KMPW2010FormFactors<PToV>::f_perp_T(const double &) const
    {
        return 0.;
    }

    double
    KMPW2010FormFactors<PToV>::f_para_T(const double &) const
    {
        return 0.;
    }

    double
    KMPW2010FormFactors<PToV>::f_long_T(const double &) const
    {
        return 0.;
    }


    double
    KMPW2010FormFactors<PToP>::_calc_z(const double & s)
    {
        return (std::sqrt(_tau_p - s) - std::sqrt(_tau_p - _tau_0)) / (std::sqrt(_tau_p - s) + std::sqrt(_tau_p - _tau_0));
    }

    KMPW2010FormFactors<PToP>::KMPW2010FormFactors(const Parameters & p, const Options &) :
        _b1_p(p["B->K::b^p_1@KMPW2010"], *this),
        _b1_0(p["B->K::b^0_1@KMPW2010"], *this),
        _b1_t(p["B->K::b^t_1@KMPW2010"], *this),
        _f0_p(p["B->K::F^p(0)@KMPW2010"], *this),
        _f0_t(p["B->K::F^t(0)@KMPW2010"], *this)
    {
    }

    FormFactors<PToP> *
    KMPW2010FormFactors<PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new KMPW2010FormFactors(parameters, options);
    }

    double
    KMPW2010FormFactors<PToP>::f_p(const double & s) const
    {
        // cf. [KMPW:2010A], Eq. (8.8), p. 30
        const double zs = _calc_z(s), z0 = _calc_z(0.0);

        return _f0_p() / (1 - s / _m_Bs2) * (1 + _b1_p() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
    }

    double
    KMPW2010FormFactors<PToP>::f_0(const double & s) const
    {
        // cf. [KMPW:2010A], Eq. (8.8), p. 30
        const double zs = _calc_z(s), z0 = _calc_z(0.0);

        return _f0_p() * (1 + _b1_0() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
    }

    double
    KMPW2010FormFactors<PToP>::f_t(const double & s) const
    {
        // cf. [KMPW:2010A], Eq. (8.8), p. 30
        const double zs = _calc_z(s), z0 = _calc_z(0.0);

        return _f0_t() / (1 - s / _m_Bs2) * (1 + _b1_t() * (zs - z0 + 0.5 * (zs * zs - z0 * z0)));
    }

    double
    KMPW2010FormFactors<PToP>::f_plus_T(const double &) const
    {
        return 0.0;
    }
}

#endif
