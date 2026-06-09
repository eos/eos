/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Nico Gubernari
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_G2026_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_G2026_IMPL_HH 1

#include <eos/form-factors/parametric-g2026.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>

#include <numeric>

namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>>
    G2026FormFactorTraits<Process_, PToV>::pole_A0_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), {
            "mass::B_u,A^0[1]@G2026",
            "mass::B_u,A^0[2]@G2026",
            "mass::B_u,A^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), {
            "mass::B_d,A^0[1]@G2026",
            "mass::B_d,A^0[2]@G2026",
            "mass::B_d,A^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), {
            "mass::B_s,A^0[1]@G2026",
            "mass::B_s,A^0[2]@G2026",
            "mass::B_s,A^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), {
            "mass::B_c,A^0[1]@G2026",
            "mass::B_c,A^0[2]@G2026",
            "mass::B_c,A^0[3]@G2026"
        }},
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>>
    G2026FormFactorTraits<Process_, PToV>::pole_V1_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), {
            "mass::B_u,V^1[1]@G2026",
            "mass::B_u,V^1[2]@G2026",
            "mass::B_u,V^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), {
            "mass::B_d,V^1[1]@G2026",
            "mass::B_d,V^1[2]@G2026",
            "mass::B_d,V^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), {
            "mass::B_s,V^1[1]@G2026",
            "mass::B_s,V^1[2]@G2026",
            "mass::B_s,V^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), {
            "mass::B_c,V^1[1]@G2026",
            "mass::B_c,V^1[2]@G2026",
            "mass::B_c,V^1[3]@G2026"
        }},
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>>
    G2026FormFactorTraits<Process_, PToV>::pole_A1_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), {
            "mass::B_u,A^1[1]@G2026",
            "mass::B_u,A^1[2]@G2026",
            "mass::B_u,A^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), {
            "mass::B_d,A^1[1]@G2026",
            "mass::B_d,A^1[2]@G2026",
            "mass::B_d,A^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), {
            "mass::B_s,A^1[1]@G2026",
            "mass::B_s,A^1[2]@G2026",
            "mass::B_s,A^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), {
            "mass::B_c,A^1[1]@G2026",
            "mass::B_c,A^1[2]@G2026",
            "mass::B_c,A^1[3]@G2026"
        }},
    };


    template<typename Process_>
    G2026FormFactors<Process_, PToV>::G2026FormFactors(const Parameters & p, const Options &) :
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
        _traits(G2026FormFactorTraits<Process_, PToV>(p)),
        _mB(_traits.m_B),
        _mV(_traits.m_V)
    {
        this->uses(_traits);
    }

    template<typename Process_>
    G2026FormFactors<Process_, PToV>::~G2026FormFactors() = default;

    template<typename Process_>
    FormFactors<PToV> *
    G2026FormFactors<Process_, PToV>::make(const Parameters & parameters, const Options & options)
    {
        return new G2026FormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    G2026FormFactors<Process_, PToV>::_par_name(const std::string & ff_name, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(idx) + "@G2026");
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::_phi(
        const double & s,
        const double & sG,
        const double & s0,
        const std::vector<double> & Mres,
        const double & tchi,
        const unsigned Kn,
        const int Ksp,
        const int Ksm,
        const int Kspm,
        const unsigned a,
        const unsigned b,
        const unsigned c,
        const unsigned d
    ) const
    {
        // Common kinematic invariants
        const double sp = power_of<2>(_mB + _mV);     // s_+
        const double sm = _traits.sm();               // s_-

        const double sqrt_sG_s   = std::sqrt(sG - s);
        const double sqrt_sG_s0  = std::sqrt(sG - s0);
        const double sqrt_sG_sm  = std::sqrt(sG - sm);
        const double sqrt_sG     = std::sqrt(sG);
        const double sqrt_sG_Q2  = std::sqrt(sG + _traits.Q2);

        // Overall normalization
        const double norm = std::sqrt(
            _traits.eta / (
                Kn * pow(sp, Ksp) * pow(sm, Ksm) * pow(sp - sm, Kspm) * M_PI * tchi
            )
        );

        // Individual factors from the new definition
        const double factor_const =
            std::pow((sG - s) / (sG - s0), 0.25) *
            (sqrt_sG_s + sqrt_sG_s0);

        const double factor_a =
            std::pow(sp - s, 0.25 * a);

        const double factor_b =
            std::pow(sqrt_sG_s + sqrt_sG_sm, 0.5 * b);

        const double factor_c =
            std::pow(sqrt_sG_s + sqrt_sG, -(c + 3.0));

        const double factor_d =
            std::pow((sqrt_sG_s + sqrt_sG) / (sqrt_sG_s + sqrt_sG_Q2), d);

        double factor_Mres = 1.0;
        for (const auto & m : Mres)
        {
            const double Mres2 = power_of<2>(m);
            const bool include_factor_Mres = (Mres2 > sG) && (Mres2 < sp);
            if (include_factor_Mres)
                factor_Mres *= (Mres2 - s) / power_of<2>(sqrt_sG_s + sqrt_sG_Q2);
        }

        return norm
            * factor_const
            * factor_a
            * factor_b
            * factor_c
            * factor_d
            * factor_Mres;
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_v(const double & q2) const
    {
        return _phi(
            q2, _traits.sV, _traits.s0V, _traits.m_R_V1, _traits.tchi_V1,
            24 /*Kn*/, 1 /*Ksp*/, 0 /*Ksm*/, 0 /*Kspm*/,
            3 /*a*/, 3 /*b*/, 1 /*c*/, 3 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_a_0(const double & q2) const
    {
        return _phi(
            q2, _traits.sA, _traits.s0A, _traits.m_R_A0, _traits.tchi_A0,
            16 /*Kn*/, 0 /*Ksp*/, 0 /*Ksm*/, 0 /*Kspm*/,
            3 /*a*/, 3 /*b*/, 1 /*c*/, 2 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_a_1(const double & q2) const
    {
        return _phi(
            q2, _traits.sA, _traits.s0A, _traits.m_R_A1, _traits.tchi_A1,
            24 /*Kn*/, -1 /*Ksp*/, 0 /*Ksm*/, 0 /*Kspm*/,
            1 /*a*/, 1 /*b*/, 1 /*c*/, 3 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_a_12(const double & q2) const
    {
        return _phi(
            q2, _traits.sA, _traits.s0A, _traits.m_R_A1, _traits.tchi_A1,
            12 /*Kn*/, 0 /*Ksp*/, 0 /*Ksm*/, -2 /*Kspm*/,
            1 /*a*/, 1 /*b*/, 2 /*c*/, 3 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_t_1(const double & q2) const
    {
        return _phi(
            q2, _traits.sV, _traits.s0V, _traits.m_R_V1, _traits.tchi_T1,
            24 /*Kn*/, 0 /*Ksp*/, 0 /*Ksm*/, 0 /*Kspm*/,
            3 /*a*/, 3 /*b*/, 2 /*c*/, 4 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_t_2(const double & q2) const
    {
        return _phi(
            q2, _traits.sA, _traits.s0A, _traits.m_R_A1, _traits.tchi_AT1,
            24 /*Kn*/, -1 /*Ksp*/, -1 /*Ksm*/, 0 /*Kspm*/,
            1 /*a*/, 1 /*b*/, 2 /*c*/, 4 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToV>::_phi_t_23(const double & q2) const
    {
        return _phi(
            q2, _traits.sA, _traits.s0A, _traits.m_R_A1, _traits.tchi_AT1,
            48 /*Kn*/, 1 /*Ksp*/, 0 /*Ksm*/, -2 /*Kspm*/,
            1 /*a*/, 1 /*b*/, 1 /*c*/, 4 /*d*/
        );
    }

    // Extra coefficient a_A12_0 calculated from the endpoint relation A_12(0) = (mB^2 - mV^2) / (8 mB mV) A_0(0)
    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::_a_A12_0() const
    {
        const double x_A12 = _phi_a_12(0.0) * _traits.blaschke_product(0.0, _traits.sA, _traits.m_R_A1)
                             * (power_of<2>(_mB) - power_of<2>(_mV)) / 8.0 / _mB / _mV;
        const double x_A0  = _phi_a_0(0.0) * _traits.blaschke_product(0.0, _traits.sA, _traits.m_R_A0);
        std::array<double, 5> a;
        a[0] = x_A12 * this->_a_A0[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_A12 * this->_a_A0[i] - x_A0 * this->_a_A12[i - 1];
        }
        const auto monomials = _traits.monomials_a(_traits.calc_z(0.0, _traits.sA, _traits.s0A));
        return std::inner_product(a.begin(), a.end(), monomials.begin(), 0.0) / (monomials[0] * x_A0);
    }

    // Extra coefficient a_A1_0 calculated from the endpoint relation A_12(sm) = (mB^2 - mV^2) / (8 mB mV) A_1(sm)
    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::_a_A1_0() const
    {
        const double x_A1  = _phi_a_1( _traits.sm()) * 16.0 * _mB * _mV * _mV / (_mB + _mV) / (_mB * _mB - _mV * _mV - _traits.sm());
        const double x_A12 = _phi_a_12(_traits.sm());
        std::array<double, 5> a;
        a[0] = x_A1 * this->_a_A12_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_A1 * this->_a_A12[i - 1] - x_A12 * this->_a_A1[i - 1];
        }
        const auto monomials = _traits.monomials_a(_traits.calc_z(_traits.sm(), _traits.sA, _traits.s0A));
        return std::inner_product(a.begin(), a.end(), monomials.begin(), 0.0) / (monomials[0] * x_A12);
    }

    // Extra coefficient a_T2_0 calculated from the endpoint relation T_1(0) = T_2(0)
    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::_a_T2_0() const
    {
        const double x_T2 = _phi_t_2(0.0) * _traits.blaschke_product(0.0, _traits.sA, _traits.m_R_A1);
        const double x_T1 = _phi_t_1(0.0) * _traits.blaschke_product(0.0, _traits.sV, _traits.m_R_V1);
        const auto monomials_T2 = _traits.monomials_a(_traits.calc_z(0.0, _traits.sA, _traits.s0A));
        const auto monomials_T1 = _traits.monomials_v(_traits.calc_z(0.0, _traits.sV, _traits.s0V));

        double a_T2_0 = x_T2 * this->_a_T1[0] * monomials_T1[0];
        for (unsigned i = 1 ; i < _a_T1.size() ; ++i)
        {
            a_T2_0 += x_T2 * this->_a_T1[i] * monomials_T1[i] - x_T1 * this->_a_T2[i - 1] * monomials_T2[i];
        }

        return a_T2_0 / (monomials_T2[0] * x_T1);
    }

    // Extra coefficient a_T23_0 calculated from the endpoint relation T_23(sm) = (mB + mV)^2 / (4 mB mV) T_2(sm)
    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::_a_T23_0() const
    {
        const double x_T23 = _phi_t_23(_traits.sm()) * (_mB + _mV) * (_mB * _mB + 3 * _mV * _mV - _traits.sm()) / 8.0 / _mB / _mV / _mV;
        const double x_T2  = _phi_t_2( _traits.sm());
        std::array<double, 5> a;
        a[0] = x_T23 * this->_a_T2_0();
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_T23 * this->_a_T2[i - 1] - x_T2 * this->_a_T23[i - 1];
        }
        const auto monomials = _traits.monomials_a(_traits.calc_z(_traits.sm(), _traits.sA, _traits.s0A));
        return std::inner_product(a.begin(), a.end(), monomials.begin(), 0.0) / (monomials[0] * x_T2);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::v(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_V.begin(), _a_V.end(), coefficients.begin());
        // resonances for 1^m
        const double blaschke     = _traits.blaschke_product(q2, _traits.sV, _traits.m_R_V1);
        const double phi          = _phi_v(q2);
        const double z            = _traits.calc_z(q2, _traits.sV, _traits.s0V);
        const auto   monomials  = _traits.monomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::a_0(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_A0.begin(), _a_A0.end(), coefficients.begin());
        // resonances for 0^m
        const double blaschke     = _traits.blaschke_product(q2, _traits.sA, _traits.m_R_A0);
        const double phi          = _phi_a_0(q2);
        const double z            = _traits.calc_z(q2, _traits.sA, _traits.s0A);
        const auto   monomials  = _traits.monomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::a_1(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_A1_0();
        std::copy(_a_A1.begin(), _a_A1.end(), coefficients.begin() + 1);
        // resonances for 1^p
        const double blaschke     = _traits.blaschke_product(q2, _traits.sA, _traits.m_R_A1);
        const double phi          = _phi_a_1(q2);
        const double z            = _traits.calc_z(q2, _traits.sA, _traits.s0A);
        const auto   monomials  = _traits.monomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::a_12(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_A12_0();
        std::copy(_a_A12.begin(), _a_A12.end(), coefficients.begin() + 1);
        // resonances for 1^p
        const double blaschke     = _traits.blaschke_product(q2, _traits.sA, _traits.m_R_A1);
        const double phi          = _phi_a_12(q2);
        const double z            = _traits.calc_z(q2, _traits.sA, _traits.s0A);
        const auto   monomials  = _traits.monomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::t_1(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_T1.begin(), _a_T1.end(), coefficients.begin());
        // resonances for T (1^m state)
        const double blaschke     = _traits.blaschke_product(q2, _traits.sV, _traits.m_R_V1);
        const double phi          = _phi_t_1(q2);
        const double z            = _traits.calc_z(q2, _traits.sV, _traits.s0V);
        const auto   monomials  = _traits.monomials_v(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::t_2(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_T2_0();
        std::copy(_a_T2.begin(), _a_T2.end(), coefficients.begin() + 1);
        // resonances for T5 (1^p state)
        const double blaschke     = _traits.blaschke_product(q2, _traits.sA, _traits.m_R_A1);
        const double phi          = _phi_t_2(q2);
        const double z            = _traits.calc_z(q2, _traits.sA, _traits.s0A);
        const auto   monomials  = _traits.monomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToV>::t_23(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_T23_0();
        std::copy(_a_T23.begin(), _a_T23.end(), coefficients.begin() + 1);
        // resonances for T (1^p state)
        const double blaschke     = _traits.blaschke_product(q2, _traits.sA, _traits.m_R_A1);
        const double phi          = _phi_t_23(q2);
        const double z            = _traits.calc_z(q2, _traits.sA, _traits.s0A);
        const auto   monomials  = _traits.monomials_a(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::a_2(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return (power_of<2>(_mB + _mV) * (power_of<2>(_mB) - power_of<2>(_mV) - s) * a_1(s)
                - 16.0 * _mB * power_of<2>(_mV) * (_mB + _mV) * a_12(s)) / lambda;
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::t_3(const double & s) const
    {
        const double lambda = eos::lambda(power_of<2>(_mB), power_of<2>(_mV), s);

        return ((power_of<2>(_mB) - power_of<2>(_mV)) * (power_of<2>(_mB) + 3.0 * power_of<2>(_mV) - s) * t_2(s)
                - 8.0 * _mB * power_of<2>(_mV) * (_mB - _mV) * t_23(s)) / lambda;
    }

    // Unused but needed to satisfy the FormFactors interface
    template<typename Process_> double G2026FormFactors<Process_, PToV>::f_perp(const double &) const { return 0.0; }
    template<typename Process_> double G2026FormFactors<Process_, PToV>::f_para(const double &) const { return 0.0; }
    template<typename Process_> double G2026FormFactors<Process_, PToV>::f_long(const double &) const { return 0.0; }
    template<typename Process_> double G2026FormFactors<Process_, PToV>::f_perp_T(const double &) const { return 0.0; }
    template<typename Process_> double G2026FormFactors<Process_, PToV>::f_para_T(const double &) const { return 0.0; }
    template<typename Process_> double G2026FormFactors<Process_, PToV>::f_long_T(const double &) const { return 0.0; }

    // EOS uses channel labels for the dispersive-bound saturations.
    // In [G:2026A], these same channels are denoted V0, A0, V1, A1, T1, and AT1.
    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::saturation_0p_v() const // [G:2026A]: saturation_V0
    {
        return 0.0;
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::saturation_0m_a() const // [G:2026A]: saturation_A0
    {
        return std::inner_product(_a_A0.begin(), _a_A0.end(), _a_A0.begin(), 0.0);
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::saturation_1m_v() const // [G:2026A]: saturation_V1
    {
        return std::inner_product(_a_V.begin(), _a_V.end(), _a_V.begin(), 0.0);
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::saturation_1p_a() const // [G:2026A]: saturation_A1
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
    G2026FormFactors<Process_, PToV>::saturation_1m_t() const // [G:2026A]: saturation_T1
    {
        return std::inner_product(_a_T1.begin(), _a_T1.end(), _a_T1.begin(), 0.0);
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToV>::saturation_1p_t5() const // [G:2026A]: saturation_AT1
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
    G2026FormFactors<Process_, PToV>::diagnostics() const
     {
        Diagnostics results;

        results.add({ _traits.calc_z(0.0,  _traits.sV, _traits.s0V), "z(q2 =  0, sV)" });
        results.add({ _traits.calc_z(11.0, _traits.sV, _traits.s0V), "z(q2 = 11, sV)" });
        results.add({ _traits.calc_z(18.0, _traits.sA, _traits.s0A), "z(q2 = 18, sA)" });
        results.add({ _traits.calc_z(-3.0, _traits.sA, _traits.s0A), "z(q2 = -3, sA)" });

        {
            const auto & [z0, z1, z2, z3, z4, z5] = _traits.monomials_v(_traits.calc_z(0.0, _traits.sV, _traits.s0V));
            results.add({ z0,              "z_0(z = z(q2 = 0, sV))" });
            results.add({ z1,              "z_1(z = z(q2 = 0, sV))" });
            results.add({ z2,              "z_2(z = z(q2 = 0, sV))" });
            results.add({ z3,              "z_3(z = z(q2 = 0, sV))" });
            results.add({ z4,              "z_4(z = z(q2 = 0, sV))" });
            results.add({ z5,              "z_5(z = z(q2 = 0, sV))" });
        }

        {
            const auto & [z0, z1, z2, z3, z4, z5] = _traits.monomials_v(_traits.calc_z(19.0, _traits.sV, _traits.s0V));
            results.add({ z0,              "z_0(z = z(q2 = 19, sV))" });
            results.add({ z1,              "z_1(z = z(q2 = 19, sV))" });
            results.add({ z2,              "z_2(z = z(q2 = 19, sV))" });
            results.add({ z3,              "z_3(z = z(q2 = 19, sV))" });
            results.add({ z4,              "z_4(z = z(q2 = 19, sV))" });
            results.add({ z5,              "z_5(z = z(q2 = 19, sV))" });
        }

        {
            const auto & [z0, z1, z2, z3, z4, z5] = _traits.monomials_a(_traits.calc_z(-3.0, _traits.sA, _traits.s0A));
            results.add({ z0,              "z_0(z = z(q2 = -3, sA))" });
            results.add({ z1,              "z_1(z = z(q2 = -3, sA))" });
            results.add({ z2,              "z_2(z = z(q2 = -3, sA))" });
            results.add({ z3,              "z_3(z = z(q2 = -3, sA))" });
            results.add({ z4,              "z_4(z = z(q2 = -3, sA))" });
            results.add({ z5,              "z_5(z = z(q2 = -3, sA))" });
        }

        {
            results.add({ _phi_v(-2.0),   "phi_v(q2 = -2)" });
            results.add({ _phi_v( 1.0),   "phi_v(q2 =  1)" });
            results.add({ _phi_v( 4.0),   "phi_v(q2 =  4)" });

            results.add({ _phi_a_0(-2.0), "phi_a_0(q2 = -2)" });
            results.add({ _phi_a_0( 1.0), "phi_a_0(q2 =  1)" });
            results.add({ _phi_a_0( 4.0), "phi_a_0(q2 =  4)" });

            results.add({ _phi_a_1(-2.0), "phi_a_1(q2 = -2)" });
            results.add({ _phi_a_1( 1.0), "phi_a_1(q2 =  1)" });
            results.add({ _phi_a_1( 4.0), "phi_a_1(q2 =  4)" });

            results.add({ _phi_a_12(-2.0), "phi_a_12(q2 = -2)" });
            results.add({ _phi_a_12( 1.0), "phi_a_12(q2 =  1)" });
            results.add({ _phi_a_12( 4.0), "phi_a_12(q2 =  4)" });

            results.add({ _phi_t_1(-2.0), "phi_t_1(q2 = -2)" });
            results.add({ _phi_t_1( 1.0), "phi_t_1(q2 =  1)" });
            results.add({ _phi_t_1( 4.0), "phi_t_1(q2 =  4)" });

            results.add({ _phi_t_2(-2.0), "phi_t_2(q2 = -2)" });
            results.add({ _phi_t_2( 1.0), "phi_t_2(q2 =  1)" });
            results.add({ _phi_t_2( 4.0), "phi_t_2(q2 =  4)" });

            results.add({ _phi_t_23(-2.0), "phi_t_23(q2 = -2)" });
            results.add({ _phi_t_23( 1.0), "phi_t_23(q2 =  1)" });
            results.add({ _phi_t_23( 4.0), "phi_t_23(q2 =  4)" });
        }

        {
            results.add({ _a_A12_0(),     "a_A12_0" });
            results.add({ _a_A1_0(),      "a_A1_0"  });
            results.add({ _a_T2_0(),      "a_T2_0"  });
            results.add({ _a_T23_0(),     "a_T23_0" });
        }

        return results;
    }


    template<typename Process_>
    const std::set<ReferenceName> G2026FormFactors<Process_, PToV>::references
    {
        "G:2026A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> G2026FormFactors<Process_, PToV>::options
    {
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    G2026FormFactors<Process_, PToV>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    G2026FormFactors<Process_, PToV>::end_options()
    {
        return options.cend();
    }


    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>>
    G2026FormFactorTraits<Process_, PToP>::pole_V0_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), {
            "mass::B_u,V^0[1]@G2026",
            "mass::B_u,V^0[2]@G2026",
            "mass::B_u,V^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), {
            "mass::B_d,V^0[1]@G2026",
            "mass::B_d,V^0[2]@G2026",
            "mass::B_d,V^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), {
            "mass::B_s,V^0[1]@G2026",
            "mass::B_s,V^0[2]@G2026",
            "mass::B_s,V^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), {
            "mass::B_c,V^0[1]@G2026",
            "mass::B_c,V^0[2]@G2026",
            "mass::B_c,V^0[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), {
            "mass::D_s,V^0[1]@G2026",
            "mass::D_s,V^0[2]@G2026",
            "mass::D_s,V^0[3]@G2026"
        }}
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::vector<std::string>>
    G2026FormFactorTraits<Process_, PToP>::pole_V1_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), {
            "mass::B_u,V^1[1]@G2026",
            "mass::B_u,V^1[2]@G2026",
            "mass::B_u,V^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), {
            "mass::B_d,V^1[1]@G2026",
            "mass::B_d,V^1[2]@G2026",
            "mass::B_d,V^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), {
            "mass::B_s,V^1[1]@G2026",
            "mass::B_s,V^1[2]@G2026",
            "mass::B_s,V^1[3]@G2026"
        }},
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), {
            "mass::B_c,V^1[1]@G2026",
            "mass::B_c,V^1[2]@G2026",
            "mass::B_c,V^1[3]@G2026"
        }},
    };


    template<typename Process_>
    G2026FormFactors<Process_, PToP>::G2026FormFactors(const Parameters & p, const Options &) :
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
        _traits(G2026FormFactorTraits<Process_, PToP>(p)),
        _mB(_traits.m_B),
        _mP(_traits.m_P)
    {
        this->uses(_traits);
    }

    template<typename Process_>
    G2026FormFactors<Process_, PToP>::~G2026FormFactors() = default;

    template<typename Process_>
    FormFactors<PToP> *
    G2026FormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new G2026FormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    G2026FormFactors<Process_, PToP>::_par_name(const std::string & ff_name, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(idx) + "@G2026");
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToP>::_phi(
        const double & s,
        const double & sG,
        const std::vector<double> & Mres,
        const double & tchi,
        const unsigned Kn,
        const int Ksp,
        const int Ksm,
        const int Kspm,
        const unsigned a,
        const unsigned b,
        const unsigned c,
        const unsigned d
    ) const
    {
        // Common kinematic invariants
        const double sp = power_of<2>(_mB + _mP);     // s_+
        const double sm = _traits.sm();               // s_-

        const double sqrt_sG_s   = std::sqrt(sG - s);
        const double sqrt_sG_s0  = std::sqrt(sG - _traits.s0V);
        const double sqrt_sG_sm  = std::sqrt(sG - sm);
        const double sqrt_sG     = std::sqrt(sG);
        const double sqrt_sG_Q2  = std::sqrt(sG + _traits.Q2);

        // Overall normalization
        const double norm = std::sqrt(
            _traits.eta / (
                Kn * pow(sp, Ksp) * pow(sm, Ksm) * pow(sp - sm, Kspm) * M_PI * tchi
            )
        );

        // Individual factors from the new definition
        const double factor_const =
            std::pow((sG - s) / (sG - _traits.s0V), 0.25) *
            (sqrt_sG_s + sqrt_sG_s0);

        const double factor_a =
            std::pow(sp - s, 0.25 * a);

        const double factor_b =
            std::pow(sqrt_sG_s + sqrt_sG_sm, 0.5 * b);

        const double factor_c =
            std::pow(sqrt_sG_s + sqrt_sG, -(c + 3.0));

        const double factor_d =
            std::pow((sqrt_sG_s + sqrt_sG) / (sqrt_sG_s + sqrt_sG_Q2), d);

        double factor_Mres = 1.0;
        for (const auto & m : Mres)
        {
            const double Mres2 = power_of<2>(m);
            const bool include_factor_Mres = (Mres2 > sG) && (Mres2 < sp);
            if (include_factor_Mres)
                factor_Mres *= (Mres2 - s) / power_of<2>(sqrt_sG_s + sqrt_sG_Q2);
        }

        return norm
            * factor_const
            * factor_a
            * factor_b
            * factor_c
            * factor_d
            * factor_Mres;
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToP>::_phi_f_p(const double & q2) const
    {
        return _phi(
            q2, _traits.sV, _traits.m_R_V1, _traits.tchi_V1,
            48 /*Kn*/, 0 /*Ksp*/, 0 /*Ksm*/, 0 /*Kspm*/,
            3 /*a*/, 3 /*b*/, 2 /*c*/, 3 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToP>::_phi_f_0(const double & q2) const
    {
        return _phi(
            q2, _traits.sV, _traits.m_R_V0, _traits.tchi_V0,
            16 /*Kn*/, -1 /*Ksp*/, -1 /*Ksm*/, 0 /*Kspm*/,
            1 /*a*/, 1 /*b*/, 1 /*c*/, 2 /*d*/
        );
    }

    template<typename Process_>
    inline double
    G2026FormFactors<Process_, PToP>::_phi_f_t(const double & q2) const
    {
        return _phi(
            q2, _traits.sV, _traits.m_R_V1, _traits.tchi_T1,
            48 /*Kn*/, 1 /*Ksp*/, 0 /*Ksm*/, 0 /*Kspm*/,
            3 /*a*/, 3 /*b*/, 1 /*c*/, 4 /*d*/
        );
    }

    // Extra coefficient a_f0_0 calculated from the endpoint relation f_0(0) = f+(0)
    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::_a_f0_0() const
    {
        const double x_f0 = _phi_f_0(0.0) * _traits.blaschke_product(0.0, _traits.sV, _traits.m_R_V0);
        const double x_fp = _phi_f_p(0.0) * _traits.blaschke_product(0.0, _traits.sV, _traits.m_R_V1);
        std::array<double, 5> a;
        a[0] = x_f0 * this->_a_fp[0];
        for (unsigned i = 1 ; i < a.size() ; ++i)
        {
            a[i] = x_f0 * this->_a_fp[i] - x_fp * this->_a_f0[i - 1];
        }
        const auto monomials = _traits.monomials(_traits.calc_z(0.0, _traits.sV, _traits.s0V));
        return std::inner_product(a.begin(), a.end(), monomials.begin(), 0.0) / (monomials[0] * x_fp);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_p(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        // resonances for 1^m
        const double blaschke = _traits.blaschke_product(q2, _traits.sV, _traits.m_R_V1);
        const double phi          = _phi_f_p(q2);
        const double z            = _traits.calc_z(q2, _traits.sV, _traits.s0V);
        const auto   monomials  = _traits.monomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_0(const double & q2) const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_f0_0();
        std::copy(_a_f0.begin(), _a_f0.end(), coefficients.begin() + 1);
        // resonances for 0^p
        const double blaschke = _traits.blaschke_product(q2, _traits.sV, _traits.m_R_V0);
        const double phi          = _phi_f_0(q2);
        const double z            = _traits.calc_z(q2, _traits.sV, _traits.s0V);
        const auto   monomials  = _traits.monomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_t(const double & q2) const
    {
        std::array<double, 5> coefficients;
        std::copy(_a_ft.begin(), _a_ft.end(), coefficients.begin());
        // resonances for 1^m
        const double blaschke = _traits.blaschke_product(q2, _traits.sV, _traits.m_R_V1);
        const double phi          = _phi_f_t(q2);
        const double z            = _traits.calc_z(q2, _traits.sV, _traits.s0V);
        const auto   monomials  = _traits.monomials(z);
        const double series       = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), 0.0);

        return series / phi / blaschke;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_p_series(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.sV, 0.0), complex<double>(_traits.s0V, 0.0));
        const auto monomials       = _traits.monomials(z);
        const complex<double> series = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), complex<double>(0.0, 0.0));

        return abs(series);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_p_series_prime(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.sV, 0.0), complex<double>(_traits.s0V, 0.0));
        const auto monomials_prime = _traits.monomials_derivatives(z);
        const complex<double> series_prime = std::inner_product(coefficients.begin(), coefficients.end(), monomials_prime.begin(), complex<double>(0.0, 0.0));

        return abs(series_prime);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_0_series(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.sV, 0.0), complex<double>(_traits.s0V, 0.0));
        const auto monomials       = _traits.monomials(z);
        const complex<double> series = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), complex<double>(0.0, 0.0));

        return abs(series);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_0_series_prime(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.sV, 0.0), complex<double>(_traits.s0V, 0.0));
        const auto monomials_prime = _traits.monomials_derivatives(z);
        const complex<double> series_prime = std::inner_product(coefficients.begin(), coefficients.end(), monomials_prime.begin(), complex<double>(0.0, 0.0));

        return abs(series_prime);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_t_series(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.sV, 0.0), complex<double>(_traits.s0V, 0.0));
        const auto monomials       = _traits.monomials(z);
        const complex<double> series = std::inner_product(coefficients.begin(), coefficients.end(), monomials.begin(), complex<double>(0.0, 0.0));

        return abs(series);
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::f_t_series_prime(const double & q2) const
    {
        std::array<complex<double>, 5> coefficients;
        std::copy(_a_fp.begin(), _a_fp.end(), coefficients.begin());
        const complex<double> z      = this->_traits.calc_z(complex<double>(q2, 0.0), complex<double>(_traits.sV, 0.0), complex<double>(_traits.s0V, 0.0));
        const auto monomials_prime = _traits.monomials_derivatives(z);
        const complex<double> series_prime = std::inner_product(coefficients.begin(), coefficients.end(), monomials_prime.begin(), complex<double>(0.0, 0.0));

        return abs(series_prime);
    }

    // Unused but needed to satisfy the FormFactors interface
    template<typename Process_> double G2026FormFactors<Process_, PToP>::f_plus_T(const double &) const { return 0.0; }

    // EOS uses channel labels for the dispersive-bound saturations.
    // In [G:2026A], these same channels are denoted V0, A0, V1, A1, T1, and AT1.
    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::saturation_0p_v() const // [G:2026A]: saturation_V0
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_f0_0();
        std::copy(_a_f0.begin(), _a_f0.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToP>::saturation_0m_a() const // [G:2026A]: saturation_A0
    {
        return 0.;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::saturation_1m_v() const // [G:2026A]: saturation_V1
    {
        return std::inner_product(_a_fp.begin(), _a_fp.end(), _a_fp.begin(), 0.0);
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToP>::saturation_1p_a() const // [G:2026A]: saturation_A1
    {
        return 0.;
    }

    template <typename Process_>
    double
    G2026FormFactors<Process_, PToP>::saturation_1m_t() const // [G:2026A]: saturation_T1
    {
        return std::inner_product(_a_ft.begin(), _a_ft.end(), _a_ft.begin(), 0.0);
    }

    template<typename Process_>
    double
    G2026FormFactors<Process_, PToP>::saturation_1p_t5() const // [G:2026A]: saturation_AT1
    {
        return 0.;
    }

    template <typename Process_>
    Diagnostics
    G2026FormFactors<Process_, PToP>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _traits.calc_z(0.0,  _traits.sV, _traits.s0V), "z(q2 =  0)" });
        results.add({ _traits.calc_z(11.0, _traits.sV, _traits.s0V), "z(q2 = 11)" });
        results.add({ _traits.calc_z(18.0, _traits.sV, _traits.s0V), "z(q2 = 18)" });
        results.add({ _traits.calc_z(-3.0, _traits.sV, _traits.s0V), "z(q2 = -3)" });

        {
            const auto & [z0, z1, z2, z3, z4, z5] = _traits.monomials(_traits.calc_z(0.0, _traits.sV, _traits.s0V));
            results.add({ z0,              "z_0(z = z(q2 = 0))" });
            results.add({ z1,              "z_1(z = z(q2 = 0))" });
            results.add({ z2,              "z_2(z = z(q2 = 0))" });
            results.add({ z3,              "z_3(z = z(q2 = 0))" });
            results.add({ z4,              "z_4(z = z(q2 = 0))" });
            results.add({ z5,              "z_5(z = z(q2 = 0))" });
        }

        {
            const auto & [z0, z1, z2, z3, z4, z5] = _traits.monomials(_traits.calc_z(19.0, _traits.sV, _traits.s0V));
            results.add({ z0,              "z_0(z = z(q2 = 19))" });
            results.add({ z1,              "z_1(z = z(q2 = 19))" });
            results.add({ z2,              "z_2(z = z(q2 = 19))" });
            results.add({ z3,              "z_3(z = z(q2 = 19))" });
            results.add({ z4,              "z_4(z = z(q2 = 19))" });
            results.add({ z5,              "z_5(z = z(q2 = 19))" });
        }

        {
            const auto & [z0, z1, z2, z3, z4, z5] = _traits.monomials(_traits.calc_z(-3.0, _traits.sV, _traits.s0V));
            results.add({ z0,              "z_0(z = z(q2 = -3))" });
            results.add({ z1,              "z_1(z = z(q2 = -3))" });
            results.add({ z2,              "z_2(z = z(q2 = -3))" });
            results.add({ z3,              "z_3(z = z(q2 = -3))" });
            results.add({ z4,              "z_4(z = z(q2 = -3))" });
            results.add({ z5,              "z_5(z = z(q2 = -3))" });
        }

        {
            results.add({ _phi_f_p(-2.0),   "phi_f_p(q2 = -2)" });
            results.add({ _phi_f_p( 1.0),   "phi_f_p(q2 =  1)" });
            results.add({ _phi_f_p( 4.0),   "phi_f_p(q2 =  4)" });

            results.add({ _phi_f_0(-2.0),   "phi_f_0(q2 = -2)" });
            results.add({ _phi_f_0( 1.0),   "phi_f_0(q2 =  1)" });
            results.add({ _phi_f_0( 4.0),   "phi_f_0(q2 =  4)" });

            results.add({ _phi_f_t(-2.0),   "phi_f_t(q2 = -2)" });
            results.add({ _phi_f_t( 1.0),   "phi_f_t(q2 =  1)" });
            results.add({ _phi_f_t( 4.0),   "phi_f_t(q2 =  4)" });
        }

        {
            results.add({ _a_f0_0(),  "_a_f0_0" });
        }

        return results;
    }

    template<typename Process_>
    const std::set<ReferenceName> G2026FormFactors<Process_, PToP>::references
    {
        "G:2026A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> G2026FormFactors<Process_, PToP>::options
    {
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    G2026FormFactors<Process_, PToP>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    G2026FormFactors<Process_, PToP>::end_options()
    {
        return options.cend();
    }
}

#endif
