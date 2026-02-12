/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020-2025 Danny van Dyk
 * Copyright (c) 2020      Nico Gubernari
 * Copyright (c) 2020      Christoph Bobeth
 * Copyright (c) 2025      Maximilian Hoverath
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
#include <eos/maths/power-of.hh>

#include <gsl/gsl_sf_dilog.h>

#include <numeric>

namespace eos
{
    using namespace std::literals::string_literals;

    template<typename Process_>
    std::string BGL1997FormFactors<Process_, PToV>::_par_name(const std::string & ff_name)
    {
        return std::string(Process_::label) + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    template<typename Process_>
    BGL1997FormFactors<Process_, PToV>::BGL1997FormFactors(const Parameters & p, const Options & o) :
        _a_g{{ UsedParameter(p[_par_name("g_0")], *this),
               UsedParameter(p[_par_name("g_1")], *this),
               UsedParameter(p[_par_name("g_2")], *this),
               UsedParameter(p[_par_name("g_3")], *this)
        }},
        _a_f{{ UsedParameter(p[_par_name("f_0")], *this),
               UsedParameter(p[_par_name("f_1")], *this),
               UsedParameter(p[_par_name("f_2")], *this),
               UsedParameter(p[_par_name("f_3")], *this)
        }},
        _a_F1{{  /* F1_0 parameter determined by identity F1(t_-) = (mB - mV) * f(t_-) */
                 UsedParameter(p[_par_name("F1_1")], *this),
                 UsedParameter(p[_par_name("F1_2")], *this),
                 UsedParameter(p[_par_name("F1_3")], *this)
        }},
        _a_F2{{  /* F2_0 parameter determined by identity between F2 and F1 at q2 = 0 */
                 UsedParameter(p[_par_name("F2_1")], *this),
                 UsedParameter(p[_par_name("F2_2")], *this),
                 UsedParameter(p[_par_name("F2_3")], *this)
        }},
        _a_T1{{  UsedParameter(p[_par_name("T1_0")], *this),
                 UsedParameter(p[_par_name("T1_1")], *this),
                 UsedParameter(p[_par_name("T1_2")], *this),
                 UsedParameter(p[_par_name("T1_3")], *this)
        }},
        _a_T2{{ /* T2_0 parameter determined by identity T1(0) = T2(0) */
                 UsedParameter(p[_par_name("T2_1")], *this),
                 UsedParameter(p[_par_name("T2_2")], *this),
                 UsedParameter(p[_par_name("T2_3")], *this)
        }},
        _a_T23{{ /* T23_0 parameter determined by identity between T2 and T23 at q2 = t_- */
                 UsedParameter(p[_par_name("T23_1")], *this),
                 UsedParameter(p[_par_name("T23_2")], *this),
                 UsedParameter(p[_par_name("T23_3")], *this)
        }},
        _traits(BGL1997FormFactorTraits<Process_, PToV>(p, o, _options)),
        _mB(_traits.m_B),
        _mV(_traits.m_V),
        t_0(_traits.t_0)
    {
    }

    template<typename Process_>
    BGL1997FormFactors<Process_, PToV>::~BGL1997FormFactors() = default;

    template<typename Process_>
    FormFactors<PToV> *
    BGL1997FormFactors<Process_, PToV>::make(const Parameters & parameters, const Options & options)
    {
        return new BGL1997FormFactors(parameters, options);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::_phi(const double & s, const double & s_0, const double & K, const unsigned & a, const unsigned & b,                                                 const unsigned & c, const double & chi) const
    {
        const double sq_tp    = std::sqrt(_traits.tp());
        const double sq_tp_t  = std::sqrt(_traits.tp() - s);
        const double sq_tp_t0 = std::sqrt(_traits.tp() - s_0);
        const double sq_tp_tm = std::sqrt(_traits.tp() - _traits.tm());

        // [BGL:1997A] eq. (4.14) for OPE at Q^2 = -q^2 = 0
        // => generalization for q^2 != 0 possible, see eq.(4.15)
        return std::sqrt(1.0 / (K * M_PI * chi)) * (sq_tp_t + sq_tp_t0)
               * std::pow(sq_tp_t / sq_tp_t0, 1.0 / 2.0)
               * std::pow(_traits.tp() - s, a / 4.0)
               * std::pow(sq_tp_t + sq_tp_tm, b / 2.0)
               * 1.0 / std::pow(sq_tp_t + sq_tp, c + 3.0);
    }

    // Functions
    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::g(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 96, 3, 3, 1, _traits.chi_1m);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = _a_g[0] + _a_g[1] * z + _a_g[2] * z * z + _a_g[3] * z * z * z;
        const double blaschke = _traits.blaschke_1m(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 24, 1, 1, 1, _traits.chi_1p);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = _a_f[0] + _a_f[1] * z + _a_f[2] * z * z + _a_f[3] * z * z * z;
        const double blaschke = _traits.blaschke_1p(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_F1_0() const
    {
        const double z    = _traits._z(_traits.tm(), _traits.t_0, _traits.tp());
        const double x_f  = _traits.blaschke_1p(_traits.tm()) * _phi(_traits.tm(), _traits.t_0, 24.0, 1, 1, 1, _traits.chi_1p);
        const double x_F1 = _traits.blaschke_1p(_traits.tm()) * _phi(_traits.tm(), _traits.t_0, 48.0, 1, 1, 2, _traits.chi_1p) * (_mB - _mV);
        std::array<double, 4> an, zn;
        zn[0] = 1.0;
        an[0]  = x_F1 * this->_a_f[0] * zn[0];
        for (unsigned i = 1 ; i < an.size() ; ++i)
        {
            an[i] = x_F1 * this->_a_f[i] - x_f * this->_a_F1[i - 1];
            zn[i] = z * zn[i - 1];
        }
        return std::inner_product(an.begin(), an.end(), zn.begin(), 0.0) / (zn[0] * x_f);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::F1(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 48, 1, 1, 2, _traits.chi_1p);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = a_F1_0() + _a_F1[0] * z + _a_F1[1] * z * z + _a_F1[2] * z * z * z;
        const double blaschke = _traits.blaschke_1p(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_F2_0() const
    {
        const double r    = _mV / _mB;
        const double wmax = (power_of<2>(_mB) + power_of<2>(_mV)) / (2.0 * _mB * _mV);

        const double z = _traits._z(0.0, _traits.t_0, _traits.tp());
        const double x_F1 = _traits.blaschke_1p(0.0) * _phi(0.0, _traits.t_0, 48.0, 1, 1, 2, _traits.chi_1p);
        const double x_F2 = _traits.blaschke_0m(0.0) * _phi(0.0, _traits.t_0, 64.0, 3, 3, 1, _traits.chi_0m) * (1.0 + r) / ((1.0 - r) * (1.0 + wmax) * r * power_of<2>(_mB));
        std::array<double, 4> an, zn;
        zn[0] = 1.0;
        an[0]  = x_F2 * this->a_F1_0() * zn[0]; // a_F1[0] is the linear coefficient; we need the constant part
        for (unsigned i = 1 ; i < an.size() ; ++i)
        {
            an[i] = x_F2 * this->_a_F1[i - 1] - x_F1 * this->_a_F2[i - 1];
            zn[i] = z * zn[i - 1];
        }
        return std::inner_product(an.begin(), an.end(), zn.begin(), 0.0) / (zn[0] * x_F1);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::F2(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 64, 3, 3, 1, _traits.chi_0m);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = a_F2_0() + _a_F2[0] * z + _a_F2[1] * z * z + _a_F2[2] * z * z * z;
        const double blaschke = _traits.blaschke_0m(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::v(const double & s) const
    {
        return (_mB + _mV) / 2.0 * g(s);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_0(const double & s) const
    {
        return F2(s) / 2.0;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_1(const double & s) const
    {
        return 1.0/(_mB + _mV) * f(s);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_2(const double & s) const
    {
        return (_mB + _mV) / eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s) * ((power_of<2>(_mB) - power_of<2>(_mV) - s) * f(s) - 2.0 * _mV * F1(s));
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_12(const double & s) const
    {
        return F1(s) / (8.0 * _mB * _mV);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::t_1(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 24.0, 3, 3, 2, _traits.chi_T_1m);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = _a_T1[0] + _a_T1[1] * z + _a_T1[2] * z * z + _a_T1[3] * z * z * z;
        const double blaschke = _traits.blaschke_1m(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_T2_0() const
    {
        const double z = _traits._z(0.0, _traits.t_0, _traits.tp());
        const double x_T2 = _traits.blaschke_1p(0.0) * _phi(0.0, _traits.t_0, 24.0 / (_traits.tp() * _traits.tm()), 1, 1, 2, _traits.chi_T_1p);
        const double x_T1 = _traits.blaschke_1m(0.0) * _phi(0.0, _traits.t_0, 24.0, 3, 3, 2, _traits.chi_T_1m);
        std::array<double, 4> an, zn;
        zn[0] = 1.0;
        an[0]  = x_T2 * this->_a_T1[0] * zn[0];
        for (unsigned i = 1 ; i < an.size() ; ++i)
        {
            an[i] = x_T2 * this->_a_T1[i] - x_T1 * this->_a_T2[i - 1];
            zn[i] = z * zn[i - 1];
        }
        return std::inner_product(an.begin(), an.end(), zn.begin(), 0.0) / (zn[0] * x_T1);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::t_2(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 24.0 / (_traits.tp() * _traits.tm()), 1, 1, 2, _traits.chi_T_1p);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = a_T2_0() + _a_T2[0] * z + _a_T2[1] * z * z + _a_T2[2] * z * z * z;
        const double blaschke = _traits.blaschke_1p(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::t_3(const double & s) const
    {
        return (
                (power_of<2>(_mB) - power_of<2>(_mV)) * (power_of<2>(_mB) + 3.0 * power_of<2>(_mV) - s) * this->t_2(s)
                - 8.0 * _mB * power_of<2>(_mV) * (_mB - _mV) * this->t_23(s)
        ) / eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::a_T23_0() const
    {
        const double z = _traits._z(_traits.tm(), _traits.t_0, _traits.tp());
        const double x_T2 = _traits.blaschke_1p(_traits.tm()) * _phi(_traits.tm(), _traits.t_0, 24.0 / (_traits.tp() * _traits.tm()), 1, 1, 2, _traits.chi_T_1p);
        const double x_T23 = _traits.blaschke_1p(_traits.tm()) * _phi(_traits.tm(), _traits.t_0, 3.0 * _traits.tp() / (power_of<2>(_mB) * power_of<2>(_mV)), 1, 1, 1, _traits.chi_T_1p) / (8.0 * _mB * power_of<2>(_mV)) * ((_mB + _mV) * (power_of<2>(_mB) + 3.0 * power_of<2>(_mV) - _traits.tm()));
        std::array<double, 4> an, zn;
        zn[0] = 1.0;
        an[0]  = x_T23 * this->a_T2_0() * zn[0]; // a_T2[0] is the linear coefficient; we need the constant part
        for (unsigned i = 1 ; i < an.size() ; ++i)
        {
            an[i] = x_T23 * this->_a_T2[i - 1] - x_T2 * this->_a_T23[i - 1];
            zn[i] = z * zn[i - 1];
        }
        return std::inner_product(an.begin(), an.end(), zn.begin(), 0.0) / (zn[0] * x_T2);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::t_23(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 3.0 * _traits.tp() / (power_of<2>(_mB) * power_of<2>(_mV)), 1.0, 1.0, 1.0, _traits.chi_T_1p);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = a_T23_0() + _a_T23[0] * z + _a_T23[1] * z * z + _a_T23[2] * z * z * z;
        const double blaschke = _traits.blaschke_1p(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f_perp(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f_para(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f_long(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f_perp_T(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f_para_T(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToV>::f_long_T(const double & /*s*/) const
    {
        return 0.0;  //  TODO
    }

    template<typename Process_>
    const std::set<ReferenceName> BGL1997FormFactors<Process_, PToV>::references
    {
        "BGL:1997A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> BGL1997FormFactors<Process_, PToV>::_options
    {
        { "n-bound-states-1m"_ok, { "1"s, "2"s, "3"s, "4"s }, "3"s },
        { "n-bound-states-1p"_ok, { "1"s, "2"s, "3"s, "4"s }, "4"s },
        { "n-bound-states-0m"_ok, { "1"s, "2"s, "3"s       }, "3"s },
        { "n-bound-states-0p"_ok, { "1"s, "2"s             }, "2"s }
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BGL1997FormFactors<Process_, PToV>::begin_options()
    {
        return _options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BGL1997FormFactors<Process_, PToV>::end_options()
    {
        return _options.cend();
    }


    template<typename Process_>
    std::string BGL1997FormFactors<Process_, PToP>::_par_name(const std::string & ff_name)
    {
        return std::string(Process_::label) + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    template<typename Process_>
    BGL1997FormFactors<Process_, PToP>::BGL1997FormFactors(const Parameters & p, const Options & o) :
        _a_f_p{{ UsedParameter(p[_par_name("f+_0")], *this),
                 UsedParameter(p[_par_name("f+_1")], *this),
                 UsedParameter(p[_par_name("f+_2")], *this),
                 UsedParameter(p[_par_name("f+_3")], *this)
        }},
        _a_f_0{{ UsedParameter(p[_par_name("f0_0")], *this),
                 UsedParameter(p[_par_name("f0_1")], *this),
                 UsedParameter(p[_par_name("f0_2")], *this),
                 UsedParameter(p[_par_name("f0_3")], *this)
        }},
        _a_f_t{{ UsedParameter(p[_par_name("fT_0")], *this),
                 UsedParameter(p[_par_name("fT_1")], *this),
                 UsedParameter(p[_par_name("fT_2")], *this),
                 UsedParameter(p[_par_name("fT_3")], *this) }},
        _traits(BGL1997FormFactorTraits<Process_, PToP>(p, o, _options)),
        _mB(_traits.m_B),
        _mP(_traits.m_P),
        t_0(_traits.t_0)
    {
    }

    template<typename Process_>
    BGL1997FormFactors<Process_, PToP>::~BGL1997FormFactors() = default;

    template<typename Process_>
    FormFactors<PToP> *
    BGL1997FormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new BGL1997FormFactors(parameters, options);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToP>::_phi(const double & s, const double & s_0, const double & K, const unsigned & a, const unsigned & b,                                                 const unsigned & c, const double & chi) const
    {
        const double sq_tp    = std::sqrt(_traits.tp());
        const double sq_tp_t  = std::sqrt(_traits.tp() - s);
        const double sq_tp_t0 = std::sqrt(_traits.tp() - s_0);
        const double sq_tp_tm = std::sqrt(_traits.tp() - _traits.tm());

        // [BGL:1997A] eq. (4.14) for OPE at Q^2 = -q^2 = 0
        // => generalization for q^2 != 0 possible, see eq.(4.15)
        return std::sqrt(1.0 / (K * M_PI * chi)) * (sq_tp_t + sq_tp_t0)
               * std::pow(sq_tp_t / sq_tp_t0, 1.0 / 2.0)
               * std::pow(_traits.tp() - s, a / 4.0)
               * std::pow(sq_tp_t + sq_tp_tm, b / 2.0)
               * 1.0 / std::pow(sq_tp_t + sq_tp, c + 3.0);
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToP>::f_p(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 48, 3, 3, 2, _traits.chi_1m);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = _a_f_p[0] + _a_f_p[1] * z + _a_f_p[2] * z * z + _a_f_p[3] * z * z * z;
        const double blaschke = _traits.blaschke_1m(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToP>::f_0(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 16, 1, 1, 1, _traits.chi_0p);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = _a_f_0[0] + _a_f_0[1] * z + _a_f_0[2] * z * z + _a_f_0[3] * z * z * z;
        const double blaschke = _traits.blaschke_0p(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToP>::f_t(const double & s) const
    {
        const double phi      = _phi(s, _traits.t_0, 48.0 * _traits.tp(), 3, 3, 1, _traits.chi_T_1m);
        const double z        = _traits._z(s, _traits.t_0, _traits.tp());
        const double series   = _a_f_t[0] + _a_f_t[1] * z + _a_f_t[2] * z * z + _a_f_t[3] * z * z * z;
        const double blaschke = _traits.blaschke_1m(s);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double BGL1997FormFactors<Process_, PToP>::f_plus_T(const double & /*s*/) const
    {
        return 0.0; //  TODO
    }

    template<typename Process_>
    const std::set<ReferenceName> BGL1997FormFactors<Process_, PToP>::references
    {
        "BGL:1997A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> BGL1997FormFactors<Process_, PToP>::_options
    {
        { "n-bound-states-1m"_ok, { "1"s, "2"s, "3"s, "4"s }, "3"s },
        { "n-bound-states-0p"_ok, { "1"s, "2"s             }, "2"s }
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BGL1997FormFactors<Process_, PToP>::begin_options()
    {
        return _options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    BGL1997FormFactors<Process_, PToP>::end_options()
    {
        return _options.cend();
    }
}

#endif
