/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2026 Simon Mutke
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

            // Limit of s_wave(s, m1, m2) for m2 -> 0. The logarithms of the general expression diverge
            // individually as log(m2), and only their sum remains finite.
            complex<double>
            s_wave_massless(const complex<double> & S, const double & m)
            {
                static const double pi = M_PI;

                // Adapt s to match Mathematica's behaviour on the branch cut
                const complex<double> s = S + complex<double>(0.0, 1e-15);

                return 1.0 / 16.0 / pi / pi * (s - m * m) / s * std::log(m * m / (m * m - s));
            }

            complex<double>
            s_wave_massless_prime(const complex<double> & S, const double & m)
            {
                static const double pi = M_PI;

                // Adapt s to match Mathematica's behaviour on the branch cut
                const complex<double> s = S + complex<double>(0.0, 1e-15);

                return 1.0 / 16.0 / pi / pi / s * (m * m / s * std::log(m * m / (m * m - s)) - 1.0);
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

                if ((m1 == 0.0) || (m2 == 0.0))
                {
                    return s_wave_massless_prime(S, m1 + m2);
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

            // The pair of poles s_+- = a +- b of the partial-fraction decomposition of n_L(s)^2, contributing
            // (s_+ L(s_+) - s_- L(s_-)) / b. The pair is even in b and regular as b -> 0, where the poles coincide
            // and the difference quotient cancels catastrophically. For |b| below 1e-4 (|a| + (m1 + m2)^2), evaluate
            // the pair at that value of |b| instead, which approximates the limit to relative accuracy O(1e-8).
            complex<double>
            pole_pair(const complex<double> & S, const complex<double> & a, const complex<double> & b, const double & m1, const double & m2)
            {
                const double          h = 1e-4 * (std::abs(a) + (m1 + m2) * (m1 + m2));
                const complex<double> c = (std::abs(b) < h) ? complex<double>(h, 0.0) : b;

                return ((a + c) * L_wave_pole(S, a + c, m1, m2) - (a - c) * L_wave_pole(S, a - c, m1, m2)) / c;
            }

            // 2 i rho(s) n_l(s)^2 with i rho(s) = -sqrt(-lambda(s)) / (16 pi s), where n_l(s)^2 is rational in
            // u = lambda(s) / (4 s q0^2). The principal square root of -lambda = (s_+ - s) (s - s_-) is negative
            // real only for real s above threshold s_+ or below pseudothreshold s_-, which places the branch
            // cuts there and nowhere else. On the cut above threshold, the small positive imaginary part
            // added to s selects the upper lip, where i rho(s) = +i |rho(s)|.
            template <unsigned l_>
            complex<double>
            l_wave_disc(const complex<double> & S, const double & m1, const double & m2, const double & q0)
            {
                static const double pi = M_PI;

                const double mp = m1 + m2;
                const double mm = m1 - m2;

                const complex<double> s       = S + complex<double>(0.0, 1e-15);
                const complex<double> mlambda = (mp * mp - s) * (s - mm * mm);
                const complex<double> disc    = -2.0 * std::sqrt(mlambda) / 16.0 / pi / s;

                const complex<double> u = -mlambda / (4.0 * s * q0 * q0);

                if constexpr (l_ == 0)
                {
                    return disc;
                }
                else if constexpr (l_ == 1)
                {
                    return disc * u / (u + 1.0);
                }
                else if constexpr (l_ == 2)
                {
                    return disc * u * u / (u * u + 3.0 * u + 9.0);
                }
                else
                {
                    static_assert(l_ == 3);

                    return disc * u * u * u / (u * u * u + 6.0 * u * u + 45.0 * u + 225.0);
                }
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

            if ((m1 == 0.0) || (m2 == 0.0))
            {
                return impl::s_wave_massless(S, m1 + m2);
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

            const double          mp   = m1 + m2;
            const double          mm   = m1 - m2;
            const double          q0sq = q0 * q0;
            const complex<double> a    = m1 * m1 + m2 * m2 - 2.0 * q0sq;
            const complex<double> b    = std::sqrt(a * a - mp * mp * mm * mm);

            return s_wave(S, m1, m2) - 2.0 * q0sq * impl::pole_pair(S, a, b, m1, m2);
        }

        complex<double>
        d_wave(const complex<double> & S, const double & m, const double & q0)
        {
            return d_wave(S, m, m, q0);
        }

        complex<double>
        d_wave(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            const double          mp   = m1 + m2;
            const double          mm   = m1 - m2;
            const double          q0sq = q0 * q0;
            const complex<double> u1   = 3.0 / 2.0 * complex<double>(-1, std::sqrt(3.0));
            const complex<double> u2   = -3.0 / 2.0 * complex<double>(1, std::sqrt(3.0));
            const complex<double> a1   = m1 * m1 + m2 * m2 + 2.0 * q0sq * u1;
            const complex<double> b1   = std::sqrt(a1 * a1 - mp * mp * mm * mm);
            const complex<double> a2   = m1 * m1 + m2 * m2 + 2.0 * q0sq * u2;
            const complex<double> b2   = std::sqrt(a2 * a2 - mp * mp * mm * mm);

            return s_wave(S, m1, m2)
                   + 2.0 * q0sq * (u1 * u1 / (2.0 * u1 + 3.0) * impl::pole_pair(S, a1, b1, m1, m2) + u2 * u2 / (2.0 * u2 + 3.0) * impl::pole_pair(S, a2, b2, m1, m2));
        }

        complex<double>
        f_wave(const complex<double> & S, const double & m, const double & q0)
        {
            return f_wave(S, m, m, q0);
        }

        complex<double>
        f_wave(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            const double          mp   = m1 + m2;
            const double          mm   = m1 - m2;
            const double          q0sq = q0 * q0;
            const double          acbr = std::cbrt((75.0 * std::sqrt(5.0) - 151.0) / 2.0);
            const double          u1   = -2.0 + acbr - 11.0 / acbr;
            const complex<double> u2   = -2.0 + acbr * complex<double>(-1, std::sqrt(3.0)) / 2.0 + 5.5 / acbr * complex<double>(1, std::sqrt(3.0));
            const complex<double> u3   = -2.0 - acbr * complex<double>(1, std::sqrt(3.0)) / 2.0 - 5.5 / acbr * complex<double>(-1, std::sqrt(3.0));
            const complex<double> a1   = m1 * m1 + m2 * m2 + 2.0 * q0sq * u1;
            const complex<double> b1   = std::sqrt(a1 * a1 - mp * mp * mm * mm);
            const complex<double> a2   = m1 * m1 + m2 * m2 + 2.0 * q0sq * u2;
            const complex<double> b2   = std::sqrt(a2 * a2 - mp * mp * mm * mm);
            const complex<double> a3   = m1 * m1 + m2 * m2 + 2.0 * q0sq * u3;
            const complex<double> b3   = std::sqrt(a3 * a3 - mp * mp * mm * mm);

            return s_wave(S, m1, m2)
                   + 2.0 * q0sq
                             * (u1 * u1 * u1 / (3.0 * u1 * u1 + 12.0 * u1 + 45.0) * impl::pole_pair(S, a1, b1, m1, m2)
                                + u2 * u2 * u2 / (3.0 * u2 * u2 + 12.0 * u2 + 45.0) * impl::pole_pair(S, a2, b2, m1, m2)
                                + u3 * u3 * u3 / (3.0 * u3 * u3 + 12.0 * u3 + 45.0) * impl::pole_pair(S, a3, b3, m1, m2));
        }

        complex<double>
        s_wave_disc(const complex<double> & S, const double & m)
        {
            return impl::l_wave_disc<0>(S, m, m, 1.0);
        }

        complex<double>
        s_wave_disc(const complex<double> & S, const double & m1, const double & m2)
        {
            return impl::l_wave_disc<0>(S, m1, m2, 1.0);
        }

        complex<double>
        p_wave_disc(const complex<double> & S, const double & m, const double & q0)
        {
            return impl::l_wave_disc<1>(S, m, m, q0);
        }

        complex<double>
        p_wave_disc(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            return impl::l_wave_disc<1>(S, m1, m2, q0);
        }

        complex<double>
        d_wave_disc(const complex<double> & S, const double & m, const double & q0)
        {
            return impl::l_wave_disc<2>(S, m, m, q0);
        }

        complex<double>
        d_wave_disc(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            return impl::l_wave_disc<2>(S, m1, m2, q0);
        }

        complex<double>
        f_wave_disc(const complex<double> & S, const double & m, const double & q0)
        {
            return impl::l_wave_disc<3>(S, m, m, q0);
        }

        complex<double>
        f_wave_disc(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            return impl::l_wave_disc<3>(S, m1, m2, q0);
        }
    } // namespace chew_mandelstam
} // namespace eos
