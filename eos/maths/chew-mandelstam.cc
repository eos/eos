/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2026 Simon Mutke
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

            // The derivative of s_wave with respect to s.
            complex<double>
            s_wave_prime(const complex<double> & S, const double & m)
            {
                static const double pi = M_PI;

                // Adapt s to match Mathematica's behaviour on the branch cut
                const complex<double> s = S + complex<double>(0.0, 1e-15);

                const complex<double> sqkallen = std::sqrt(s * (4.0 * m * m - s));

                return -1.0 / 16.0 / pi / pi / s * (1.0 - 4.0 * m * m / sqkallen * atan_near_branch_point(s / sqkallen, s));
            }

            complex<double>
            s_wave_prime(const complex<double> & S, const double & m1, const double & m2)
            {
                if (m1 == m2)
                {
                    return s_wave_prime(S, m1);
                }

                static const double pi = M_PI;

                const double mp = m1 + m2;
                const double mm = m1 - m2;

                // Adapt s to match Mathematica's behaviour on the branch cut
                const complex<double> s = S + complex<double>(0.0, 1e-15);

                const complex<double> sqkallen = std::sqrt((s - mp * mp) * (s - mm * mm));

                return -1.0 / 16.0 / pi / pi / s
                       * (1.0 + (sqkallen / s + (m1 * m1 + m2 * m2 - s) / sqkallen) * std::log((m1 * m1 + m2 * m2 - s + sqkallen) / 2.0 / m1 / m2)
                          - mp * mm / s * std::log(m1 / m2));
            }

            // One pole of the partial-fraction decomposition of n_L(s)^2, regular at s = pole.
            complex<double>
            L_wave_pole(const complex<double> & S, const complex<double> & pole, const double & m1, const double & m2)
            {
                const double mp = m1 + m2;

                // The difference quotient degenerates at the pole; use its limit there.
                if (std::abs(S - pole) < 1e-7)
                {
                    return s_wave(pole, m1, m2) / (mp * mp - pole) + s_wave_prime(pole, m1, m2);
                }

                return (s_wave(S, m1, m2) + (S - mp * mp) / (mp * mp - pole) * s_wave(pole, m1, m2)) / (S - pole);
            }
        } // namespace impl

        complex<double>
        s_wave(const complex<double> & S, const double & m)
        {
            static const double pi = M_PI;

            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s = S + complex<double>(0.0, 1e-15);

            return -1.0 / 8.0 / pi / pi * std::sqrt(4.0 * m * m - s) * impl::atan_near_branch_point(s / std::sqrt(s * (4.0 * m * m - s)), s) / std::sqrt(s);
        }

        complex<double>
        s_wave(const complex<double> & S, const double & m1, const double & m2)
        {
            if (m1 == m2)
            {
                return s_wave(S, m1);
            }

            static const double pi = M_PI;

            const double mp = m1 + m2;
            const double mm = m1 - m2;

            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s = S + complex<double>(0.0, 1e-15);

            const complex<double> sqkallen = std::sqrt((s - mp * mp) * (s - mm * mm));

            return 1.0 / 16.0 / pi / pi * (sqkallen / s * std::log((m1 * m1 + m2 * m2 - s + sqkallen) / 2.0 / m1 / m2) - mp * mm / s * (1.0 - s / mp / mp) * std::log(m1 / m2));
        }

        complex<double>
        p_wave(const complex<double> & S, const double & m, const double & q0)
        {
            const double          mp    = 2.0 * m;
            const complex<double> delta = mp * mp - 4.0 * q0 * q0;

            // Partial fractions leave a single difference quotient of s_wave, singular only at s = delta.
            if (std::abs(S - delta) < 1e-7)
            {
                return -4.0 * q0 * q0 * impl::s_wave_prime(delta, m);
            }

            return (S - mp * mp) / (S - delta) * (s_wave(S, m) - s_wave(delta, m));
        }

        complex<double>
        p_wave(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            if (m1 == m2)
            {
                return p_wave(S, m1, q0);
            }

            const double          mp      = m1 + m2;
            const double          mm      = m1 - m2;
            const double          q0sq    = q0 * q0;
            const complex<double> a       = m1 * m1 + m2 * m2 - 2.0 * q0sq;
            const complex<double> b       = std::sqrt(a * a - mp * mp * mm * mm);
            const complex<double> s1plus  = a + b;
            const complex<double> s1minus = a - b;

            return s_wave(S, m1, m2) + 2.0 * q0sq / b * (s1minus * impl::L_wave_pole(S, s1minus, m1, m2) - s1plus * impl::L_wave_pole(S, s1plus, m1, m2));
        }

        complex<double>
        d_wave(const complex<double> & S, const double & m, const double & q0)
        {
            return d_wave(S, m, m, q0);
        }

        complex<double>
        d_wave(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            const double          mp      = m1 + m2;
            const double          mm      = m1 - m2;
            const double          q0sq    = q0 * q0;
            const complex<double> u1      = 3.0 / 2.0 * complex<double>(-1, std::sqrt(3.0));
            const complex<double> u2      = -3.0 / 2.0 * complex<double>(1, std::sqrt(3.0));
            const complex<double> a1      = m1 * m1 + m2 * m2 + 2.0 * q0sq * u1;
            const complex<double> b1      = std::sqrt(a1 * a1 - mp * mp * mm * mm);
            const complex<double> a2      = m1 * m1 + m2 * m2 + 2.0 * q0sq * u2;
            const complex<double> b2      = std::sqrt(a2 * a2 - mp * mp * mm * mm);
            const complex<double> s1plus  = a1 + b1;
            const complex<double> s1minus = a1 - b1;
            const complex<double> s2plus  = a2 + b2;
            const complex<double> s2minus = a2 - b2;

            return s_wave(S, m1, m2)
                   + 2.0 * q0sq
                             * (u1 * u1 / (2.0 * u1 + 3.0) / b1 * (s1plus * impl::L_wave_pole(S, s1plus, m1, m2) - s1minus * impl::L_wave_pole(S, s1minus, m1, m2))
                                + u2 * u2 / (2.0 * u2 + 3.0) / b2 * (s2plus * impl::L_wave_pole(S, s2plus, m1, m2) - s2minus * impl::L_wave_pole(S, s2minus, m1, m2)));
        }

    } // namespace chew_mandelstam
} // namespace eos
