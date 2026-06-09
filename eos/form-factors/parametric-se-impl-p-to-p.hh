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


#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SE_IMPL_P_TO_P_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_SE_IMPL_P_TO_P_HH 1

#include <eos/form-factors/parametric-se.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/diagnostics.hh>

#include <numeric>

namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SEFormFactorTraits<Process_, PToP>::resonance_0p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,0@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,0@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    SEFormFactorTraits<Process_, PToP>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template<typename Process_>
    SEFormFactors<Process_, PToP>::SEFormFactors(const Parameters & p, const Options &) :
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
        _traits(SEFormFactorTraits<Process_, PToP>(p)),
        _mB(_traits.m_B),
        _mP(_traits.m_P)
    {
    }

    template<typename Process_>
    SEFormFactors<Process_, PToP>::~SEFormFactors() = default;

    template<typename Process_>
    FormFactors<PToP> *
    SEFormFactors<Process_, PToP>::make(const Parameters & parameters, const Options & options)
    {
        return new SEFormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    SEFormFactors<Process_, PToP>::_par_name(const std::string & ff_name, unsigned idx) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(idx) + "@SE");
    }

    template<typename Process_>
    double
    SEFormFactors<Process_, PToP>::_phi(const double & t, const double & threshold_tp, const double & chi,
                                             const int & A, const unsigned B, const unsigned C, const unsigned k,
                                             const unsigned p, const unsigned n, const unsigned m) const
    {
        // [GRvDV:2023A]
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
    SEFormFactors<Process_, PToP>::_phi_f_p(const double & q2) const
    {
        return _phi(q2, _traits.tp, Process_::chi_1m_v, 0, 0, 0, 1, 2, 2, 3);
    }

    template<typename Process_>
    inline double
    SEFormFactors<Process_, PToP>::_phi_f_0(const double & q2) const
    {
        return _phi(q2, _traits.tp, Process_::chi_0p_v, 1, 1, 0, 3, 2, 1, 1);
    }

    template<typename Process_>
    inline double
    SEFormFactors<Process_, PToP>::_phi_f_t(const double & q2) const
    {
        return _phi(q2, _traits.tp, Process_::chi_1m_t, -1, 0, 0, 1, 0, 3, 3);
    }

    template <typename Process_>
    double
    SEFormFactors<Process_, PToP>::_a_f0_0() const
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
    SEFormFactors<Process_, PToP>::f_p(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_0(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_t(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_p_series(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_p_series_prime(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_0_series(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_0_series_prime(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_t_series(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_t_series_prime(const double & q2) const
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
    SEFormFactors<Process_, PToP>::f_plus_T(const double & q2) const
    {
        return f_t(q2) * q2 / _mB / (_mB + _mP);
    }

    template <typename Process_>
    double
    SEFormFactors<Process_, PToP>::saturation_0p_v() const
    {
        std::array<double, 5> coefficients;
        coefficients[0] = _a_f0_0();
        std::copy(_a_f0.begin(), _a_f0.end(), coefficients.begin() + 1);

        return std::inner_product(coefficients.begin(), coefficients.end(), coefficients.begin(), 0.0);
    }

    template<typename Process_>
    double
    SEFormFactors<Process_, PToP>::saturation_0m_a() const
    {
        return 0.;
    }

    template <typename Process_>
    double
    SEFormFactors<Process_, PToP>::saturation_1m_v() const
    {
        return std::inner_product(_a_fp.begin(), _a_fp.end(), _a_fp.begin(), 0.0);
    }

    template<typename Process_>
    double
    SEFormFactors<Process_, PToP>::saturation_1p_a() const
    {
        return 0.;
    }

    template <typename Process_>
    double
    SEFormFactors<Process_, PToP>::saturation_1m_t() const
    {
        return std::inner_product(_a_ft.begin(), _a_ft.end(), _a_ft.begin(), 0.0);
    }

    template<typename Process_>
    double
    SEFormFactors<Process_, PToP>::saturation_1p_t5() const
    {
        return 0.;
    }

    template <typename Process_>
    Diagnostics
    SEFormFactors<Process_, PToP>::diagnostics() const
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
            const auto & [p0, p1, p2, p3, p4, p5] = _traits.orthonormal_polynomials(_traits.calc_z(10.0, _traits.tp, _traits.t0));
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
    const std::set<ReferenceName> SEFormFactors<Process_, PToP>::references
    {
        "BFW:2010A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> SEFormFactors<Process_, PToP>::options
    {
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    SEFormFactors<Process_, PToP>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    SEFormFactors<Process_, PToP>::end_options()
    {
        return options.cend();
    }
}

#endif
