/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2017 Danny van Dyk
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

#include <eos/form-factors/parametric-bfvd2014.hh>

namespace eos
{
    BFvD2014FormFactors::BFvD2014FormFactors(const Parameters & p, const Options &) :
        _f_long_v(p["Lambda_b->Lambda::f_0^V(0)@BFvD2014"], *this),
        _b_1_long_v(p["Lambda_b->Lambda::b_1_0^V@BFvD2014"], *this),
        _f_long_a(p["Lambda_b->Lambda::f_0^A(0)@BFvD2014"], *this),
        _b_1_long_a(p["Lambda_b->Lambda::b_1_0^A@BFvD2014"], *this),
        _f_perp_v(p["Lambda_b->Lambda::f_perp^V(0)@BFvD2014"], *this),
        _b_1_perp_v(p["Lambda_b->Lambda::b_1_perp^V@BFvD2014"], *this),
        _f_perp_a(p["Lambda_b->Lambda::f_perp^A(0)@BFvD2014"], *this),
        _b_1_perp_a(p["Lambda_b->Lambda::b_1_perp^A@BFvD2014"], *this),
        _m_lambda_b(p["mass::Lambda_b"], *this),
        _m_lambda(p["mass::Lambda"], *this)
    {
    }

    FormFactors<OneHalfPlusToOneHalfPlus> *
    BFvD2014FormFactors::make(const Parameters & parameters, const Options & options)
    {
        return new BFvD2014FormFactors(parameters, options);
    }

    double
    BFvD2014FormFactors::f_long_v(const double & s) const
    {
        const double tp = power_of<2>(_m_lambda_b + _m_lambda);
        const double zt = _z(s, tp, 12.0), z0 = _z(0.0, tp, 12.0);

        return _f_long_v() / (1.0 - s / mv2) * (1.0 + _b_1_long_v() * (zt - z0));
    }

    double
    BFvD2014FormFactors::f_perp_v(const double & s) const
    {
        const double tp = power_of<2>(_m_lambda_b + _m_lambda);
        const double zt = _z(s, tp, 12.0), z0 = _z(0, tp, 12.0);

        return _f_perp_v() / (1.0 - s / mv2) * (1.0 + _b_1_perp_v() * (zt - z0));
    }

    double
    BFvD2014FormFactors::f_long_a(const double & s) const
    {
        const double tp = power_of<2>(_m_lambda_b + _m_lambda);
        const double zt = _z(s, tp, 12.0), z0 = _z(0, tp, 12.0);

        return _f_long_a() / (1.0 - s / ma2) * (1.0 + _b_1_long_a() * (zt - z0));
    }

    double
    BFvD2014FormFactors::f_perp_a(const double & s) const
    {
        const double tp = power_of<2>(_m_lambda_b + _m_lambda);
        const double zt = _z(s, tp, 12.0), z0 = _z(0, tp, 12.0);

        return _f_perp_a() / (1.0 - s / ma2) * (1.0 + _b_1_perp_a() * (zt - z0));
    }

    const std::set<ReferenceName> BFvD2014FormFactors::references
    {
        "DKMR:2017A"_rn
    };

    const std::vector<OptionSpecification> BFvD2014FormFactors::options
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BFvD2014FormFactors::begin_options()
    {
        return options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BFvD2014FormFactors::end_options()
    {
        return options.cend();
    }
}
