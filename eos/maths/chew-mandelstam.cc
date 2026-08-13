/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
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

#include <eos/maths/chew-mandelstam.hh>
#include <eos/utils/exception.hh>

#include <cmath>

namespace eos
{
    namespace chew_mandelstam
    {
        namespace impl
        {
            // std::atan has branch points at +-i. s_wave evaluates it at
            // z = s / sqrt(s * (mp^2 - s)), which tends to exactly +i as
            // mp -> 0 (physically negligible channel masses, e.g. the e^+e^-
            // channel) whenever the Mandelstam variable s is (numerically) real.
            // Right at that limit, z's real part -- mathematically zero -- is
            // instead the result of several chained sqrt/division roundings, so
            // its sign (and hence which side of atan's branch cut is selected)
            // is not portable across standard library implementations. Pin it
            // to zero, but only when s itself is real: when the caller supplies
            // a genuinely complex s (e.g. to move onto a different Riemann
            // sheet), z's small real part is physically meaningful, not noise,
            // and must not be touched.
            complex<double>
            atan_near_branch_point(const complex<double> & z, const complex<double> & s)
            {
                if ((std::abs(s.imag()) < 1.0e-10) && (std::abs(z.real()) < 1.0e-9 * std::abs(z.imag())))
                {
                    return std::atan(complex<double>(0.0, z.imag()));
                }

                return std::atan(z);
            }
        } // namespace impl

        complex<double>
        s_wave(const complex<double> & S, const double & m1, const double & m2)
        {
            if (m1 != m2)
            {
                throw InternalError("chew_mandelstam::s_wave: only equal masses (m1 == m2) are currently implemented");
            }

            static const double pi = M_PI;

            const double          mp = m1 + m2;
            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s  = S + complex<double>(0.0, 1e-15);

            return -1.0 / 8.0 / pi / pi * std::sqrt(mp * mp - s) * impl::atan_near_branch_point(s / std::sqrt(s * (mp * mp - s)), s) / std::sqrt(s);
        }
    } // namespace chew_mandelstam
} // namespace eos
