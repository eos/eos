/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_HKVT2025_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_HKVT2025_IMPL_HH 1

#include <eos/form-factors/parametric-hkvt2025.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>

#include <boost/math/special_functions/legendre.hpp>

#include <numeric>
#include <iostream>
namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    HKVT2025FormFactorTraits<Process_, PToPP>::resonance_0m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    HKVT2025FormFactorTraits<Process_, PToPP>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    HKVT2025FormFactorTraits<Process_, PToPP>::resonance_1p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,1@BSZ2015" }
    };

    template<typename Process_>
    HKVT2025FormFactors<Process_, PToPP>::HKVT2025FormFactors(const Parameters & p, const Options &) :
        _traits(HKVT2025FormFactorTraits<Process_, PToPP>(p)),
        _mP1(_traits.m_P1),
        _mP2(_traits.m_P2),
        _mP3(_traits.m_P3),
        _numWaves(Process_::num_waves)
    {
        for (unsigned l = 0 ; l < _numWaves ; l++)
        {
            _a_g_1.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("g_1", l,  0)],  *this), UsedParameter(p[_par_name("g_1", l,  1)],  *this),
                                                                                  UsedParameter(p[_par_name("g_1", l,  2)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("g_1", l,  5)],  *this), UsedParameter(p[_par_name("g_1", l,  6)],  *this),
                                                                                  UsedParameter(p[_par_name("g_1", l,  7)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("g_1", l, 10)],  *this), UsedParameter(p[_par_name("g_1", l, 11)],  *this),
                                                                                  UsedParameter(p[_par_name("g_1", l, 12)],  *this) } } });

            _a_g_2.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("g_2", l,  0)],  *this), UsedParameter(p[_par_name("g_2", l,  1)],  *this),
                                                                                  UsedParameter(p[_par_name("g_2", l,  2)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("g_2", l,  5)],  *this), UsedParameter(p[_par_name("g_2", l,  6)],  *this),
                                                                                  UsedParameter(p[_par_name("g_2", l,  7)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("g_2", l, 10)],  *this), UsedParameter(p[_par_name("g_2", l, 11)],  *this),
                                                                                  UsedParameter(p[_par_name("g_2", l, 12)],  *this) } } });

            _a_f_1.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("f_1", l,  0)],  *this), UsedParameter(p[_par_name("f_1", l,  1)],  *this),
                                                                                  UsedParameter(p[_par_name("f_1", l,  2)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("f_1", l,  5)],  *this), UsedParameter(p[_par_name("f_1", l,  6)],  *this),
                                                                                  UsedParameter(p[_par_name("f_1", l,  7)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("f_1", l, 10)],  *this), UsedParameter(p[_par_name("f_1", l, 11)],  *this),
                                                                                  UsedParameter(p[_par_name("f_1", l, 12)],  *this) } } });

            _a_f_2.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("f_2", l,  0)],  *this), UsedParameter(p[_par_name("f_2", l,  1)],  *this),
                                                                                  UsedParameter(p[_par_name("f_2", l,  2)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("f_2", l,  5)],  *this), UsedParameter(p[_par_name("f_2", l,  6)],  *this),
                                                                                  UsedParameter(p[_par_name("f_2", l,  7)],  *this) } ,
                                                                                { UsedParameter(p[_par_name("f_2", l, 10)],  *this), UsedParameter(p[_par_name("f_2", l, 11)],  *this),
                                                                                  UsedParameter(p[_par_name("f_2", l, 12)],  *this) } } });

            _a_F1_1.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("F1_1", l,  0)],  *this), UsedParameter(p[_par_name("F1_1", l,  1)],  *this),
                                                                                   UsedParameter(p[_par_name("F1_1", l,  2)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F1_1", l,  5)],  *this), UsedParameter(p[_par_name("F1_1", l,  6)],  *this),
                                                                                   UsedParameter(p[_par_name("F1_1", l,  7)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F1_1", l, 10)],  *this), UsedParameter(p[_par_name("F1_1", l, 11)],  *this),
                                                                                   UsedParameter(p[_par_name("F1_1", l, 12)],  *this) } } });

            _a_F1_2.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("F1_2", l,  0)],  *this), UsedParameter(p[_par_name("F1_2", l,  1)],  *this),
                                                                                   UsedParameter(p[_par_name("F1_2", l,  2)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F1_2", l,  5)],  *this), UsedParameter(p[_par_name("F1_2", l,  6)],  *this),
                                                                                   UsedParameter(p[_par_name("F1_2", l,  7)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F1_2", l, 10)],  *this), UsedParameter(p[_par_name("F1_2", l, 11)],  *this),
                                                                                   UsedParameter(p[_par_name("F1_2", l, 12)],  *this) } } });

            _a_F2_1.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("F2_1", l,  0)],  *this), UsedParameter(p[_par_name("F2_1", l,  1)],  *this),
                                                                                   UsedParameter(p[_par_name("F2_1", l,  2)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F2_1", l,  5)],  *this), UsedParameter(p[_par_name("F2_1", l,  6)],  *this),
                                                                                   UsedParameter(p[_par_name("F2_1", l,  7)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F2_1", l, 10)],  *this), UsedParameter(p[_par_name("F2_1", l, 11)],  *this),
                                                                                   UsedParameter(p[_par_name("F2_1", l, 12)],  *this) } } });

            _a_F2_2.emplace_back(std::array<std::array<UsedParameter, 3>, 3> { { { UsedParameter(p[_par_name("F2_2", l,  0)],  *this), UsedParameter(p[_par_name("F2_2", l,  1)],  *this),
                                                                                   UsedParameter(p[_par_name("F2_2", l,  2)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F2_2", l,  5)],  *this), UsedParameter(p[_par_name("F2_2", l,  6)],  *this),
                                                                                   UsedParameter(p[_par_name("F2_2", l,  7)],  *this) } ,
                                                                                 { UsedParameter(p[_par_name("F2_2", l, 10)],  *this), UsedParameter(p[_par_name("F2_2", l, 11)],  *this),
                                                                                   UsedParameter(p[_par_name("F2_2", l, 12)],  *this) } } });
        }
    }

    template<typename Process_>
    HKVT2025FormFactors<Process_, PToPP>::~HKVT2025FormFactors() = default;

    template<typename Process_>
    FormFactors<PToPP> *
    HKVT2025FormFactors<Process_, PToPP>::make(const Parameters & parameters, const Options & options)
    {
        return new HKVT2025FormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    HKVT2025FormFactors<Process_, PToPP>::_par_name(const std::string & ff_name, unsigned wave, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(wave) + "_" + stringify(idx) + "@HKvT2025");
    }

    template<typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::_phi(const int & l, const double & s, const double & q2, const double & q2p, const double & chi,
                        const unsigned a, const unsigned b, const double N) const
    {
        // Note: The factor 1 / 2 / pi at the end originates from the change of q^2 -> z and the 1/pi in front of the dispersive integral
        const double norm = std::sqrt(N / 256 / M_PI / M_PI / M_PI / chi / (2 * l + 1) / 2 / M_PI);

        const double z = _traits.calc_z(q2, q2p, _traits.q20);

        // kinematic_q2p, kinematic_q2m depend on s
        const double kinematic_q2p = power_of<2>(_mP1 + std::sqrt(s));
        const double kinematic_q2m = power_of<2>(_mP1 - std::sqrt(s));

        // set Q^2 to 0
        const double q2_term = 1 / (2 * (q2p + std::sqrt(q2p) * std::sqrt(q2p - q2)) - q2);
        const double lambda_q2_term = (kinematic_q2p - q2) * power_of<2>(std::sqrt(q2p - q2) + std::sqrt(q2p - kinematic_q2m));
        const double sqrtjac_q2 = std::sqrt(4 * (1 + z) * (_traits.q20 - q2p) / power_of<3>(z - 1));

        return norm * sqrtjac_q2 * pow(lambda_q2_term, (2 * l + 1 - 2 * a) * 0.25) * pow(q2_term, b * 0.5);
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_g(const double & q2, const double & s, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, s, q2, _traits.q2p_v, Process_::chi_1m_v, 0, 4, l * (l + 1) / 48.0);
        }
        else
        {
            return 1.0;
        }
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_f(const double & q2, const double & s, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, s, q2, _traits.q2p_a, Process_::chi_1p_a, 1, 4, l * (l + 1) / 3.0);
        }
        else
        {
            return 1.0;
        }
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_F1(const double & q2, const double & s, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, s, q2, _traits.q2p_a, Process_::chi_1p_a, 1, 5, 1.0 / 12.0);
        }
        else
        {
            return _phi(0, s, q2, _traits.q2p_a, Process_::chi_1p_a, -1, 5, 1.0 / 12.0);
        }
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_F2(const double & q2, const double & s, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, s, q2, _traits.q2p_a, Process_::chi_0m_a, 0, 4, 1.0);
        }
        else
        {
            return _phi(0, s, q2, _traits.q2p_a, Process_::chi_0m_a, 0, 4, 1.0);
        }
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::v_perp(const double & q2, const double & s, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_1m) < _traits.q2p_v ? _traits.calc_z(q2, _traits.q2p_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const double z             = _traits.calc_z(q2, _traits.q2p_v, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_v(z, s);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);
        const double kinpref = std::pow(_traits.kappa(s, q2), l - 1) / std::sqrt(s);

        complex<double> res = 0.0;

        const double phi = _phi_g(q2, s, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_g_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_g_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        return kinpref * res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::a_par(const double & q2, const double & s, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_1p) < _traits.q2p_a ? _traits.calc_z(q2, _traits.q2p_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const double z             = _traits.calc_z(q2, _traits.q2p_a, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_a(z, s);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);
        const double kinpref = std::pow(_traits.kappa(s, q2), l - 1) / std::sqrt(s);

        complex<double> res = 0.0;

        const double phi = _phi_f(q2, s, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_f_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_f_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        return kinpref * res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::a_0(const double & q2, const double & s, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_1p) < _traits.q2p_a ? _traits.calc_z(q2, _traits.q2p_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const double z             = _traits.calc_z(q2, _traits.q2p_a, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_a(z, s);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);
        // Note, we absorb one power of sqrt(lambda_q3) in the tensorial structure and one by reducing the power of kappa
        const double kinpref = std::pow(_traits.kappa(s, q2), l - 1) * std::sqrt((s - power_of<2>(_mP2 + _mP3))*(s - power_of<2>(_mP2 - _mP3))) / s; // * (_mP1 * _mP1 - s) / s;

        complex<double> res = 0.0;

        const double phi = _phi_F1(q2, s, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F1_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F1_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        if (l != 0)
        {
            res *= kinpref;
        }

        return res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::a_t(const double & q2, const double & s, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_0m) < _traits.q2p_a ? _traits.calc_z(q2, _traits.q2p_a, power_of<2>(_traits.m_R_0m)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const double z             = _traits.calc_z(q2, _traits.q2p_a, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_a(z, s);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);
        const double kinpref = std::pow(_traits.kappa(s, q2), l);

        complex<double> res = 0.0;

        const double phi = _phi_F2(q2, s, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F2_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F2_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        if (l == 0)
        {
            res *= (_mP1 * _mP1 - s) / s;
        }
        else
        {
            res *= kinpref;
        }

        return res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::unitarity_integrand_0m(const double & s, const unsigned & l, const bool & iso) const
    {
        double res = 0.0;
        double kinpref = std::pow(std::sqrt((s - power_of<2>(_mP2 + _mP3))*(s - power_of<2>(_mP2 - _mP3))) / s, 2 * l + 1);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const std::array<complex<double>, 3> y_polynomials  = _traits.threshold_improved_polynomials(y, l);

        if (l == 0)
        {
            kinpref *= power_of<2>((_mP1 * _mP1 - s) / s);
        }

        for (unsigned i = 0 ; i < 3 ; i++)
        {
            if (!iso)
            {
                res += std::norm(double(_a_F2_1[l][0][i]) * y_polynomials[0] + double(_a_F2_1[l][1][i]) * y_polynomials[1] + double(_a_F2_1[l][2][i]) * y_polynomials[2]);
            }
            else
            {
                res += std::norm(double(_a_F2_2[l][0][i]) * y_polynomials[0] + double(_a_F2_2[l][1][i]) * y_polynomials[1] + double(_a_F2_2[l][2][i]) * y_polynomials[2]);
            }
        }

        return res * kinpref;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::unitarity_integrand_1p(const double & s, const unsigned & l, const bool & iso) const
    {
        double res_f = 0.0;
        double res_F1 = 0.0;
        double kinpref_f = s * std::pow(std::sqrt((s - power_of<2>(_mP2 + _mP3))*(s - power_of<2>(_mP2 - _mP3))) / s, 2 * l + 1);
        double kinpref_F1 = std::pow(std::sqrt((s - power_of<2>(_mP2 + _mP3))*(s - power_of<2>(_mP2 - _mP3))) / s, 2 * l + 1);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const std::array<complex<double>, 3> y_polynomials  = _traits.threshold_improved_polynomials(y, l);

        for (unsigned i = 0 ; i < 3 ; i++)
        {
            if (!iso)
            {
                res_f += std::norm(double(_a_f_1[l][0][i]) * y_polynomials[0] + double(_a_f_1[l][1][i]) * y_polynomials[1] + double(_a_f_1[l][2][i]) * y_polynomials[2]);
                res_F1 += std::norm(double(_a_F1_1[l][0][i]) * y_polynomials[0] + double(_a_F1_1[l][1][i]) * y_polynomials[1] + double(_a_F1_1[l][2][i]) * y_polynomials[2]);
            }
            else
            {
                res_f += std::norm(double(_a_f_2[l][0][i]) * y_polynomials[0] + double(_a_f_2[l][1][i]) * y_polynomials[1] + double(_a_f_2[l][2][i]) * y_polynomials[2]);
                res_F1 += std::norm(double(_a_F1_2[l][0][i]) * y_polynomials[0] + double(_a_F1_2[l][1][i]) * y_polynomials[1] + double(_a_F1_2[l][2][i]) * y_polynomials[2]);
            }
        }

        return res_f * kinpref_f + res_F1 * kinpref_F1;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::unitarity_integrand_1m(const double & s, const unsigned & l, const bool & iso) const
    {
        double res = 0.0;
        double kinpref = s * std::pow(std::sqrt((s - power_of<2>(_mP2 + _mP3))*(s - power_of<2>(_mP2 - _mP3))) / s, 2 * l + 1);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(s), complex<double>(_traits.sin1), complex<double>(_traits.s0)) : _traits.calc_y(complex<double>(s), complex<double>(_traits.sin0), complex<double>(_traits.s0));
        const std::array<complex<double>, 3> y_polynomials  = _traits.threshold_improved_polynomials(y, l);

        for (unsigned i = 0 ; i < 3 ; i++)
        {
            if (!iso)
            {
                res += std::norm(double(_a_g_1[l][0][i]) * y_polynomials[0] + double(_a_g_1[l][1][i]) * y_polynomials[1] + double(_a_g_1[l][2][i]) * y_polynomials[2]);
            }
            else
            {
                res += std::norm(double(_a_g_2[l][0][i]) * y_polynomials[0] + double(_a_g_2[l][1][i]) * y_polynomials[1] + double(_a_g_2[l][2][i]) * y_polynomials[2]);
            }
        }

        return res * kinpref;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_perp(const double & /* q2 */, const double & /* k2 */, const double & /* z */) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_para(const double & /* q2 */, const double & /* k2 */, const double & /* z */) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_long(const double & /* q2 */, const double & /* k2 */, const double & /* z */) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_time(const double & /* q2 */, const double & /* k2 */, const double & /* z */) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::f_perp_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::f_para_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::f_long_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::f_time_im_res_qhat2(const double & /*q2*/, const double & /*k2*/) const
    {
        throw InternalError("Not yet implemented");

        return 0.0;
    }

    template <typename Process_>
    Diagnostics
    HKVT2025FormFactors<Process_, PToPP>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _traits.calc_y(4 * 0.135 * 0.135,  _traits.sin1, _traits.s0), "y(s = 4*0.135^2)" });
        results.add({ _traits.calc_y(0.1,                _traits.sin1, _traits.s0), "y(s = 0.1)"       });

        results.add({ _traits.calc_z(0.0,  _traits.q2p_a, _traits.q20), "z_a(q2 =  0)" });
        results.add({ _traits.calc_z(0.0,  _traits.q2p_v, _traits.q20), "z_v(q2 =  0)" });
        results.add({ _traits.calc_z(10.0, _traits.q2p_a, _traits.q20), "z_a(q2 = 10)" });
        results.add({ _traits.calc_z(10.0, _traits.q2p_v, _traits.q20), "z_v(q2 = 10)" });

        {
            const auto & [p0, p1, p2] = _traits.orthonormal_polynomials_v(0.0, 0.1);
            results.add({ p0,              "p_0(z = 0.0, s = 0.1)" });
            results.add({ p1,              "p_1(z = 0.0, s = 0.1)" });
            results.add({ p2,              "p_2(z = 0.0, s = 0.1)" });
        }

        {
            const auto & [p0, p1, p2] = _traits.orthonormal_polynomials_v(_traits.calc_z(10.0, _traits.q2p_v, _traits.q20), 0.1);
            results.add({ p0,              "p_0(z = z(q2 = 10, s = 0.1))" });
            results.add({ p1,              "p_1(z = z(q2 = 10, s = 0.1))" });
            results.add({ p2,              "p_2(z = z(q2 = 10, s = 0.1))" });
        }

        {
            const auto & [p0, p1, p2] = _traits.threshold_improved_polynomials(_traits.calc_y(0.5, _traits.sin1, _traits.s0), 0);
            results.add({ p0,              "p_0(y = y(s = 0.5), 0)" });
            results.add({ p1,              "p_1(y = y(s = 0.5), 0)" });
            results.add({ p2,              "p_2(y = y(s = 0.5), 0)" });
        }

        {
            const auto & [p0, p1, p2] = _traits.threshold_improved_polynomials(_traits.calc_y(0.5, _traits.sin1, _traits.s0), 1);
            results.add({ p0,              "p_0(y = y(s = 0.5), 1)" });
            results.add({ p1,              "p_1(y = y(s = 0.5), 1)" });
            results.add({ p2,              "p_2(y = y(s = 0.5), 1)" });
        }

        {
            const auto & [p0, p1, p2] = _traits.threshold_improved_polynomials(_traits.calc_y(0.5, _traits.sin1, _traits.s0), 2);
            results.add({ p0,              "p_0(y = y(s = 0.5), 2)" });
            results.add({ p1,              "p_1(y = y(s = 0.5), 2)" });
            results.add({ p2,              "p_2(y = y(s = 0.5), 2)" });
        }

        {
            results.add({ _phi_g(-2.0, 0.1, 1),     "phi_g(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)" });
            results.add({ _phi_g( 1.0, 0.2, 1),     "phi_g(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)" });
            results.add({ _phi_g( 4.0, 0.1, 2),     "phi_g(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)" });

            results.add({ _phi_f(-2.0, 0.1, 1),     "phi_f(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)" });
            results.add({ _phi_f( 1.0, 0.2, 1),     "phi_f(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)" });
            results.add({ _phi_f( 4.0, 0.1, 2),     "phi_f(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)" });

            results.add({ _phi_F1(-2.0, 0.1, 1),    "phi_F1(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)" });
            results.add({ _phi_F1( 1.0, 0.2, 1),    "phi_F1(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)" });
            results.add({ _phi_F1( 4.0, 0.1, 2),    "phi_F1(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)" });
            results.add({ _phi_F1( 4.0, 0.1, 0),    "phi_F1(z = z(q2 =  4.0), y = y(s = 0.1), l = 0)" });

            results.add({ _phi_F2(-2.0, 0.1, 1),    "phi_F2(z = z(q2 = -2.0), y = y(s = 0.1), l = 1)" });
            results.add({ _phi_F2( 1.0, 0.2, 1),    "phi_F2(z = z(q2 =  1.0), y = y(s = 0.2), l = 1)" });
            results.add({ _phi_F2( 4.0, 0.1, 2),    "phi_F2(z = z(q2 =  4.0), y = y(s = 0.1), l = 2)" });
            results.add({ _phi_F2( 4.0, 0.1, 0),    "phi_F2(z = z(q2 =  4.0), y = y(s = 0.1), l = 0)" });
        }

        return results;
    }


    template<typename Process_>
    const std::set<ReferenceName> HKVT2025FormFactors<Process_, PToPP>::references
    {
        "HKvT:2025A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> HKVT2025FormFactors<Process_, PToPP>::options
    {
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HKVT2025FormFactors<Process_, PToPP>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HKVT2025FormFactors<Process_, PToPP>::end_options()
    {
        return options.cend();
    }
}
#endif
