/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_WILSON_COEFFICIENTS_HH
#define EOS_GUARD_SRC_UTILS_WILSON_COEFFICIENTS_HH 1

#include <src/utils/parameters.hh>
#include <src/utils/qcd.hh>

#include <array>
#include <cmath>

namespace eos
{
    template <typename Tag_> struct WilsonCoefficients;

    struct BToS {};

    template <> struct WilsonCoefficients<BToS>
    {
        /* Order: c1..c6, cq3..cq6, c2b, c7..c10 */
        std::array<double, 15> _coefficients;

        double _alpha_s;

        // Unknown Basis
        inline double c1() const { return _coefficients[0]; }
        inline double c2() const { return _coefficients[1]; }
        inline double c3() const { return _coefficients[2]; }
        inline double c4() const { return _coefficients[3]; }
        inline double c5() const { return _coefficients[4]; }
        inline double c6() const { return _coefficients[5]; }

        inline double cq3() const { return _coefficients[6]; }
        inline double cq4() const { return _coefficients[7]; }
        inline double cq5() const { return _coefficients[8]; }
        inline double cq6() const { return _coefficients[9]; }

        inline double c2b() const { return _coefficients[10]; }

        inline double c7() const { return 4.0 * M_PI / _alpha_s * _coefficients[11]; }
        inline double c8() const { return 4.0 * M_PI / _alpha_s * _coefficients[12]; }
        inline double c9() const { return 4.0 * M_PI / _alpha_s * _coefficients[13]; }
        inline double c10() const { return 4.0 * M_PI / _alpha_s * _coefficients[14]; }
    };

    WilsonCoefficients<BToS> evolve(const std::array<double, 15> & wc_qcd_0,
            const std::array<double, 15> & wc_qcd_1,
            const std::array<double, 15> & wc_qcd_2,
            const double & alpha_s_0, const double & alpha_s,
            const double & nf, const QCD::BetaFunction & beta);

    void calculate_wilson_coefficients(const double & mu, Parameters & parameters);
}

#endif
