/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SSE_IMPL_P_TO_P_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SSE_IMPL_P_TO_P_HH 1

#include <eos/form-factors/parametric-sse.hh>
#include <eos/maths/power-of.hh>

namespace eos
{
    // P -> P
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SSEFormFactorTraits<Process_, PToP>::resonance_0p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,0@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SSEFormFactorTraits<Process_, PToP>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template <typename Process_>
    template <typename Parameter_>
    complex<double>
    SSEFormFactors<Process_, PToP>::_calc_ff(const complex<double> & s, const double & m_R, const std::array<Parameter_, 3> & a) const
    {
        const complex<double> a_0(a[0]), a_1(a[1]), a_2(a[2]);

        const complex<double> z = _traits.calc_z(s);
        return 1.0 / (1.0 - s / power_of<2>(m_R)) *
                (a_0 + a_1 * z + a_2 * power_of<2>(z));
    }

    template <typename Process_>
    std::string
    SSEFormFactors<Process_, PToP>::_par_name(const std::string & ff_name)
    {
        return std::string(Process_::label) + std::string("::alpha^") + ff_name + std::string("@SSE");
    }

    template <typename Process_>
    SSEFormFactors<Process_, PToP>::SSEFormFactors(const Parameters & p, const Options &) :
        _a_fp{{ UsedParameter(p[_par_name("f+_0")], *this),
                UsedParameter(p[_par_name("f+_1")], *this),
                UsedParameter(p[_par_name("f+_2")], *this) }},
        _a_ft{{ UsedParameter(p[_par_name("fT_0")], *this),
                UsedParameter(p[_par_name("fT_1")], *this),
                UsedParameter(p[_par_name("fT_2")], *this) }},
        _a_fz{{ UsedParameter(p[_par_name("f0_1")], *this),
                UsedParameter(p[_par_name("f0_2")], *this) }},
        _traits(p),
        _mB(_traits.m_B),
        _mP(_traits.m_P)
    {
    }

    template <typename Process_>
    SSEFormFactors<Process_, PToP>::~SSEFormFactors()
    {
    }

    template <typename Process_>
    FormFactors<PToP> *
    SSEFormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new SSEFormFactors(parameters, options);
    }

    template <typename Process_>
    complex<double>
    SSEFormFactors<Process_, PToP>::f_p(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_fp);
    }

    template <typename Process_>
    complex<double>
    SSEFormFactors<Process_, PToP>::f_t(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_ft);
    }

    template <typename Process_>
    complex<double>
    SSEFormFactors<Process_, PToP>::f_0(const complex<double> & s) const
    {
        // Enforce the equation of motion f_0(0) = f_+(0), which removes the leading
        // coefficient of f_0 as a free parameter. The z expansion is not shifted
        // (z(0) != 0), so simply reusing the raw coefficient a^f+_0 does NOT enforce
        // f_0(0) = f_+(0). The constant term of f_0 must instead absorb the
        // z(0)-dependent difference of the sub-leading coefficients:
        //   c_0 = a^f+_0 + (a^f+_1 - a^f0_1) z0 + (a^f+_2 - a^f0_2) z0^2 .
        const double z0  = _traits.calc_z(0.0);
        const double c_0 = _a_fp[0]
            + (_a_fp[1] - _a_fz[0]) * z0
            + (_a_fp[2] - _a_fz[1]) * power_of<2>(z0);

        std::array<double, 3> values
        {{
            c_0,
            _a_fz[1 - 1],
            _a_fz[2 - 1],
        }};

        return _calc_ff(s, _traits.m_R_0p, values);
    }

    template <typename Process_>
    complex<double>
    SSEFormFactors<Process_, PToP>::f_plus_T(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_ft) * s / _mB() / (_mB + _mP);
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, PToP>::f_p(const double & s) const
    {
        return real(f_p(complex<double>(s)));
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, PToP>::f_t(const double & s) const
    {
        return real(f_t(complex<double>(s)));
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, PToP>::f_0(const double & s) const
    {
        return real(f_0(complex<double>(s)));
    }

    template <typename Process_>
    double
    SSEFormFactors<Process_, PToP>::f_plus_T(const double & s) const
    {
        return real(f_plus_T(complex<double>(s)));
    }
}

#endif
