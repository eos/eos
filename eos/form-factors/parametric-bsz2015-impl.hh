/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Frederik Beaujean
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BSZ2015_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BSZ2015_IMPL_HH 1

#include <eos/form-factors/parametric-bsz2015.hh>
#include <eos/maths/power-of.hh>

namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BSZ2015FormFactorTraits<Process_, PToV>::resonance_0m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BSZ2015FormFactorTraits<Process_, PToV>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BSZ2015FormFactorTraits<Process_, PToV>::resonance_1p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::down), "mass::D_d,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,1@BSZ2015" }
    };

    template <typename Process_>
    template <typename Parameter_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::_calc_ff(const complex<double> & s, const double & m_R, const std::array<Parameter_, 3> & a) const
    {
        const complex<double> a_0(a[0]), a_1(a[1]), a_2(a[2]);

        const complex<double> diff_z = _traits.calc_z(s) - _traits.calc_z(0.0);
        return 1.0 / (1.0 - s / power_of<2>(m_R)) *
                (a_0 + a_1 * diff_z + a_2 * power_of<2>(diff_z));
    }

    template <typename Process_>
    std::string
    BSZ2015FormFactors<Process_, PToV>::_par_name(const std::string & ff_name)
    {
        return std::string(Process_::label) + std::string("::alpha^") + ff_name + std::string("@BSZ2015");
    }

    template <typename Process_>
    BSZ2015FormFactors<Process_, PToV>::BSZ2015FormFactors(const Parameters & p, const Options &) :
        _a_A0{{  UsedParameter(p[_par_name("A0_0")],  *this),
                    UsedParameter(p[_par_name("A0_1")],  *this),
                    UsedParameter(p[_par_name("A0_2")],  *this) }},
        _a_A1{{  UsedParameter(p[_par_name("A1_0")],  *this),
                    UsedParameter(p[_par_name("A1_1")],  *this),
                    UsedParameter(p[_par_name("A1_2")],  *this) }},
        _a_V{{   UsedParameter(p[_par_name("V_0")],   *this),
                    UsedParameter(p[_par_name("V_1")],   *this),
                    UsedParameter(p[_par_name("V_2")],   *this) }},
        _a_T1{{  UsedParameter(p[_par_name("T1_0")],  *this),
                    UsedParameter(p[_par_name("T1_1")],  *this),
                    UsedParameter(p[_par_name("T1_2")],  *this) }},
        _a_T23{{ UsedParameter(p[_par_name("T23_0")], *this),
                    UsedParameter(p[_par_name("T23_1")], *this),
                    UsedParameter(p[_par_name("T23_2")], *this) }},
        _a_A12{{ UsedParameter(p[_par_name("A12_1")], *this),
                    UsedParameter(p[_par_name("A12_2")], *this) }},
        _a_T2{{  UsedParameter(p[_par_name("T2_1")],  *this),
                    UsedParameter(p[_par_name("T2_2")],  *this) }},
        _traits(p),
        _mB(_traits.m_B),
        _mV(_traits.m_V)
    {
    }

    template <typename Process_>
    BSZ2015FormFactors<Process_, PToV>::~BSZ2015FormFactors()
    {
    }

    template <typename Process_>
    FormFactors<PToV> *
    BSZ2015FormFactors<Process_, PToV>::make(const Parameters & parameters, const Options & options)
    {
        return new BSZ2015FormFactors(parameters, options);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::v(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_V);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::a_0(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_0m, _a_A0);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::a_1(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1p, _a_A1);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::a_12(const complex<double> & s) const
    {
        // use constraint (B.6) in [BSZ2015] to remove A_12(0)
        std::array<double, 3> values
        {{
            (power_of<2>(_mB) - power_of<2>(_mV)) / (8.0 * _mB * _mV) * _a_A0[0],
            _a_A12[1 - 1],
            _a_A12[2 - 1],
        }};

        return _calc_ff(s, _traits.m_R_1p, values);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::a_2(const complex<double> & s) const
    {
        const complex<double> lambda = eos::lambda(complex<double>(power_of<2>(_mB), 0.0), complex<double>(power_of<2>(_mV), 0.0), s);

        return (power_of<2>(_mB + _mV) * (power_of<2>(_mB) - power_of<2>(_mV) - s) * a_1(s)
                - 16.0 * _mB * power_of<2>(_mV) * (_mB + _mV) * a_12(s)) / lambda;
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::t_1(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_T1);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::t_2(const complex<double> & s) const
    {
        // use constraint T_1(0) = T_2(0) to replace T_2(0)
        std::array<double, 3> values
        {{
            _a_T1[0],
            _a_T2[1 - 1],
            _a_T2[2 - 1],
        }};
        return _calc_ff(s, _traits.m_R_1p, values);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToV>::t_23(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1p, _a_T23);
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::v(const double & s) const
    {
        return real(v(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::a_0(const double & s) const
    {
        return real(a_0(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::a_1(const double & s) const
    {
        return real(a_1(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::a_12(const double & s) const
    {
        return real(a_12(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::a_2(const double & s) const
    {
        return real(a_2(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::t_1(const double & s) const
    {
        return real(t_1(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::t_2(const double & s) const
    {
        return real(t_2(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::t_23(const double & s) const
    {
        return real(t_23(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::t_3(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return ((power_of<2>(_mB) - power_of<2>(_mV)) * (power_of<2>(_mB) + 3.0 * power_of<2>(_mV) - s) * t_2(s)
                - 8.0 * _mB * power_of<2>(_mV) * (_mB - _mV) * t_23(s)) / lambda;
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::f_perp(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return pow(2 * lambda, 0.5) / _mB / (_mB + _mV) * v(s);
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::f_para(const double & s) const
    {
        return pow(2, 0.5) * (_mB + _mV) / _mB * a_1(s);
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::f_long(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return ((power_of<2>(_mB) - power_of<2>(_mV) - s) * power_of<2>(_mB + _mV) * a_1(s) - lambda * a_2(s))
                / (2 * _mV * power_of<2>(_mB) * (_mB + _mV));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::f_perp_T(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return pow(2*lambda, 0.5) / power_of<2>(_mB) * t_1(s);
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::f_para_T(const double & s) const
    {
        return pow(2, 0.5) * (power_of<2>(_mB) - power_of<2>(_mV)) / power_of<2>(_mB) * t_2(s);
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToV>::f_long_T(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return s * (power_of<2>(_mB) + 3 * power_of<2>(_mV) - s) / (2 * power_of<3>(_mB) * _mV) * t_2(s)
                - s * lambda / (2 * power_of<3>(_mB) * _mV * (power_of<2>(_mB) - power_of<2>(_mV))) * t_3(s);
    }


    // P -> P
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    BSZ2015FormFactorTraits<Process_, PToP>::resonance_0p_names
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
    BSZ2015FormFactorTraits<Process_, PToP>::resonance_1m_names
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
    BSZ2015FormFactors<Process_, PToP>::_calc_ff(const complex<double> & s, const double & m_R, const std::array<Parameter_, 3> & a) const
    {
        const complex<double> a_0(a[0]), a_1(a[1]), a_2(a[2]);

        const complex<double> diff_z = _traits.calc_z(s) - _traits.calc_z(0.0);
        return 1.0 / (1.0 - s / power_of<2>(m_R)) *
                (a_0 + a_1 * diff_z + a_2 * power_of<2>(diff_z));
    }

    template <typename Process_>
    std::string
    BSZ2015FormFactors<Process_, PToP>::_par_name(const std::string & ff_name)
    {
        return std::string(Process_::label) + std::string("::alpha^") + ff_name + std::string("@BSZ2015");
    }

    template <typename Process_>
    BSZ2015FormFactors<Process_, PToP>::BSZ2015FormFactors(const Parameters & p, const Options &) :
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
    BSZ2015FormFactors<Process_, PToP>::~BSZ2015FormFactors()
    {
    }

    template <typename Process_>
    FormFactors<PToP> *
    BSZ2015FormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new BSZ2015FormFactors(parameters, options);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToP>::f_p(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_fp);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToP>::f_t(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_ft);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToP>::f_0(const complex<double> & s) const
    {
        // use equation of motion to replace f_0(0) by f_+(0)
        std::array<double, 3> values
        {{
            _a_fp[0],
            _a_fz[1 - 1],
            _a_fz[2 - 1],
        }};

        return _calc_ff(s, _traits.m_R_0p, values);
    }

    template <typename Process_>
    complex<double>
    BSZ2015FormFactors<Process_, PToP>::f_plus_T(const complex<double> & s) const
    {
        return _calc_ff(s, _traits.m_R_1m, _a_ft) * s / _mB() / (_mB + _mP);
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToP>::f_p(const double & s) const
    {
        return real(f_p(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToP>::f_t(const double & s) const
    {
        return real(f_t(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToP>::f_0(const double & s) const
    {
        return real(f_0(complex<double>(s)));
    }

    template <typename Process_>
    double
    BSZ2015FormFactors<Process_, PToP>::f_plus_T(const double & s) const
    {
        return real(f_plus_T(complex<double>(s)));
    }
}

#endif
