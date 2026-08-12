/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SSE_IMPL_ONEHALFPLUS_TO_ONEHALFPLUS_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SSE_IMPL_ONEHALFPLUS_TO_ONEHALFPLUS_HH 1

#include <eos/form-factors/parametric-sse.hh>
#include <eos/maths/power-of.hh>

namespace eos
{
    // 1/2^+ -> 1/2^+
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SSEFormFactorTraits<Process_, OneHalfPlusToOneHalfPlus>::resonance_0m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::up), "mass::D_u@HME" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SSEFormFactorTraits<Process_, OneHalfPlusToOneHalfPlus>::resonance_0p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,0@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,0@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,0@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,0@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d,0@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::up), "mass::D_u,0@HME" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SSEFormFactorTraits<Process_, OneHalfPlusToOneHalfPlus>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d^*@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::up), "mass::D_u^*@HME" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SSEFormFactorTraits<Process_, OneHalfPlusToOneHalfPlus>::resonance_1p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,1@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,1@HME" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,1@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,1@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d,1@HME" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::up), "mass::D_u,1@HME" }
    };

    template <typename Process_>
    std::string
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_par_name(const std::string & pol, const std::string & current, unsigned idx)
    {
        return std::string(Process_::label) + std::string("::alpha^(") + pol + "," + current + ")_" + std::to_string(idx) + std::string("@SSE");
    }

    template <typename Process_>
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::SSEFormFactors(const Parameters & p, const Options &) :
        _a_long_v{{ UsedParameter(p[_par_name("0", "V", 0)], *this),
                    UsedParameter(p[_par_name("0", "V", 1)], *this),
                    UsedParameter(p[_par_name("0", "V", 2)], *this) }},
        _a_perp_v{{ UsedParameter(p[_par_name("perp", "V", 0)], *this),
                    UsedParameter(p[_par_name("perp", "V", 1)], *this),
                    UsedParameter(p[_par_name("perp", "V", 2)], *this) }},
        _a_long_a{{ UsedParameter(p[_par_name("0", "A", 0)], *this),
                    UsedParameter(p[_par_name("0", "A", 1)], *this),
                    UsedParameter(p[_par_name("0", "A", 2)], *this) }},
        _a_long_t{{ UsedParameter(p[_par_name("0", "T", 0)], *this),
                    UsedParameter(p[_par_name("0", "T", 1)], *this),
                    UsedParameter(p[_par_name("0", "T", 2)], *this) }},
        _a_perp_t5{{ UsedParameter(p[_par_name("perp", "T5", 0)], *this),
                     UsedParameter(p[_par_name("perp", "T5", 1)], *this),
                     UsedParameter(p[_par_name("perp", "T5", 2)], *this) }},
        // leading coefficient replaced by an equation of motion
        _a_time_v{{ UsedParameter(p[_par_name("t", "V", 1)], *this),
                    UsedParameter(p[_par_name("t", "V", 2)], *this) }},
        _a_time_a{{ UsedParameter(p[_par_name("t", "A", 1)], *this),
                    UsedParameter(p[_par_name("t", "A", 2)], *this) }},
        _a_perp_a{{ UsedParameter(p[_par_name("perp", "A", 1)], *this),
                    UsedParameter(p[_par_name("perp", "A", 2)], *this) }},
        _a_perp_t{{ UsedParameter(p[_par_name("perp", "T", 1)], *this),
                    UsedParameter(p[_par_name("perp", "T", 2)], *this) }},
        _a_long_t5{{ UsedParameter(p[_par_name("0", "T5", 1)], *this),
                     UsedParameter(p[_par_name("0", "T5", 2)], *this) }},
        _traits(p),
        _m_1(_traits.m_1),
        _m_2(_traits.m_2)
    {
    }

    template <typename Process_>
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::~SSEFormFactors()
    {
    }

    template <typename Process_>
    FormFactors<OneHalfPlusToOneHalfPlus> *
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::make(const Parameters & parameters, const Options & options)
    {
        return new SSEFormFactors(parameters, options);
    }

    template <typename Process_>
    template <typename Parameter_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_calc_ff(const double & s, const double & m_R, const double & tp, const std::array<Parameter_, 3> & a) const
    {
        const double z = _traits.calc_z(s, tp);
        return 1.0 / (1.0 - s / power_of<2>(m_R)) * (a[0] + a[1] * z + a[2] * power_of<2>(z));
    }

    // Equation-of-motion relations that fix the leading coefficients. The z expansion is
    // not shifted (z(0) != 0), so the constant term absorbs the z(0)- (or z(t_-)-)dependent
    // difference of the sub-leading coefficients of the two form factors being related.

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_a_time_v_0() const
    {
        // f_time^V(0) = f_long^V(0); both expand in z(tp_v)
        const double z0 = _traits.calc_z(0.0, _traits.tp_v);
        return _a_long_v[0] + (_a_long_v[1] - _a_time_v[0]) * z0 + (_a_long_v[2] - _a_time_v[1]) * power_of<2>(z0);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_a_time_a_0() const
    {
        // f_time^A(0) = f_long^A(0); both expand in z(tp_a)
        const double z0 = _traits.calc_z(0.0, _traits.tp_a);
        return _a_long_a[0] + (_a_long_a[1] - _a_time_a[0]) * z0 + (_a_long_a[2] - _a_time_a[1]) * power_of<2>(z0);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_a_perp_a_0() const
    {
        // f_perp^A(t_-) = f_long^A(t_-); both share the 1^+ pole and expand in z(tp_a),
        // so the pole factors cancel.
        const double zm = _traits.calc_z(_traits.tm(), _traits.tp_a);
        return _a_long_a[0] + (_a_long_a[1] - _a_perp_a[0]) * zm + (_a_long_a[2] - _a_perp_a[1]) * power_of<2>(zm);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_a_perp_t_0() const
    {
        // f_perp^T(0) = f_perp^T5(0). f_perp^T expands in z(tp_v), f_perp^T5 in z(tp_a),
        // so the two are evaluated at different values of z(0).
        const double z0_v = _traits.calc_z(0.0, _traits.tp_v);
        const double z0_a = _traits.calc_z(0.0, _traits.tp_a);
        return _a_perp_t5[0] + _a_perp_t5[1] * z0_a + _a_perp_t5[2] * power_of<2>(z0_a)
            - _a_perp_t[0] * z0_v - _a_perp_t[1] * power_of<2>(z0_v);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::_a_long_t5_0() const
    {
        // f_long^T5(t_-) = f_perp^T5(t_-); both share the 1^+ pole and expand in z(tp_a),
        // so the pole factors cancel.
        const double zm = _traits.calc_z(_traits.tm(), _traits.tp_a);
        return _a_perp_t5[0] + (_a_perp_t5[1] - _a_long_t5[0]) * zm + (_a_perp_t5[2] - _a_long_t5[1]) * power_of<2>(zm);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_time_v(const double & s) const
    {
        const std::array<double, 3> a{{ _a_time_v_0(), _a_time_v[0], _a_time_v[1] }};
        return _calc_ff(s, _traits.m_R_0p, _traits.tp_v, a);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_long_v(const double & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _traits.tp_v, _a_long_v);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_perp_v(const double & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _traits.tp_v, _a_perp_v);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_time_a(const double & s) const
    {
        const std::array<double, 3> a{{ _a_time_a_0(), _a_time_a[0], _a_time_a[1] }};
        return _calc_ff(s, _traits.m_R_0m, _traits.tp_a, a);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_long_a(const double & s) const
    {
        return _calc_ff(s, _traits.m_R_1p, _traits.tp_a, _a_long_a);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_perp_a(const double & s) const
    {
        const std::array<double, 3> a{{ _a_perp_a_0(), _a_perp_a[0], _a_perp_a[1] }};
        return _calc_ff(s, _traits.m_R_1p, _traits.tp_a, a);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_long_t(const double & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _traits.tp_v, _a_long_t);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_perp_t(const double & s) const
    {
        const std::array<double, 3> a{{ _a_perp_t_0(), _a_perp_t[0], _a_perp_t[1] }};
        return _calc_ff(s, _traits.m_R_1m, _traits.tp_v, a);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_long_t5(const double & s) const
    {
        const std::array<double, 3> a{{ _a_long_t5_0(), _a_long_t5[0], _a_long_t5[1] }};
        return _calc_ff(s, _traits.m_R_1p, _traits.tp_a, a);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, OneHalfPlusToOneHalfPlus>::f_perp_t5(const double & s) const
    {
        return _calc_ff(s, _traits.m_R_1p, _traits.tp_a, _a_perp_t5);
    }
}

#endif
