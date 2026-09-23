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
#include <eos/maths/power-of.hh>
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
            static const double pi = M_PI;

            const double          mp    = 2.0 * m;
            const complex<double> delta = mp * mp - 4.0 * q0 * q0;

            // Fix the behavior near delta by Taylor expanding to first order
            if (std::abs(S - delta) < 1e-5)
            {
                return s_wave(S, m) - s_wave(delta, m) - 4.0 * q0 * q0 * (s_wave(delta + 1e-8, m) - s_wave(delta - 1e-8, m)) / 2e-8;
            }

            // Adapt s to match Mathematica's behaviour on the branch cut
            const complex<double> s = S + complex<double>(0.0, 1e-15);

            // Squared Blatt-Weisskopf form factor for l = 1, cf. PDG's resonance review, eq. (50.26):
            // F(z)^2 = 1 / (z^2 + 1), with z = sqrt(s - mp^2) / (2 q0)
            const complex<double> Fsq = 1.0 / ((s - mp * mp) / (4.0 * q0 * q0) + 1.0);

            complex<double> leading_term;
            // Fix the behavior near threshold by Taylor expanding to second order
            if (std::abs(s - mp * mp) < 1e-7)
            {
                leading_term = Fsq * (mp * mp - s) / 16.0 / mp / mp / pi / pi * (-2.0 * (mp * mp - s) + mp * pi * std::sqrt(mp * mp - s));
            }
            else
            {
                leading_term = Fsq * power_of<3>(std::sqrt(mp * mp - s)) * impl::atan_near_branch_point(s / std::sqrt(s * (mp * mp - s)), s) / 8.0 / pi / pi / std::sqrt(s);
            }

            const complex<double> loop_correction = -power_of<3>(q0) * (mp * mp - s) * std::atan(std::sqrt(delta) / 2.0 / q0) / pi / pi / std::sqrt(delta) / (s - delta);

            return (leading_term + loop_correction) / 4.0 / q0 / q0;
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
            const complex<double> zsq     = (S - mp * mp) * (S - mm * mm) / 4.0 / q0sq / S;

            // Fix the behavior near s1plus by Taylor expanding to first order
            if (std::abs(S - s1plus) < 1e-5)
            {
                return s_wave(S, m1, m2) * (1.0 + 2.0 * q0sq / b * s1minus / (S - s1minus))
                       + 2.0 * q0sq / b
                                 * (s1minus / (mp * mp - s1minus) / (S - s1minus) * (S - mp * mp) * s_wave(s1minus, m1, m2)
                                    - s1plus * (s_wave(s1plus, m1, m2) / (mp * mp - s1plus) + (s_wave(s1plus + 1e-8, m1, m2) - s_wave(s1plus - 1e-8, m1, m2)) / 2e-8));
            }

            // Fix the behavior near s1minus by Taylor expanding to first order
            if (std::abs(S - s1minus) < 1e-5)
            {
                return s_wave(S, m1, m2) * (1.0 - 2.0 * q0sq / b * s1plus / (S - s1plus))
                       - 2.0 * q0sq / b
                                 * (s1plus / (mp * mp - s1plus) / (S - s1plus) * (S - mp * mp) * s_wave(s1plus, m1, m2)
                                    - s1minus * (s_wave(s1minus, m1, m2) / (mp * mp - s1minus) + (s_wave(s1minus + 1e-8, m1, m2) - s_wave(s1minus - 1e-8, m1, m2)) / 2e-8));
            }

            return s_wave(S, m1, m2) * zsq / (1.0 + zsq)
                   + 2.0 * q0sq / b * (S - mp * mp)
                             * (s_wave(s1minus, m1, m2) * s1minus / (mp * mp - s1minus) / (S - s1minus) - s_wave(s1plus, m1, m2) * s1plus / (mp * mp - s1plus) / (S - s1plus));
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
            const complex<double> zsq     = (S - mp * mp) * (S - mm * mm) / 4.0 / q0sq / S;

            return s_wave(S, m1, m2) * zsq * zsq / (9.0 + 3.0 * zsq + zsq * zsq)
                   - 2.0 * q0sq * (S - mp * mp)
                             * (u1 * u1 / (2.0 * u1 + 3.0) / b1
                                        * (s_wave(s1minus, m1, m2) * s1minus / (mp * mp - s1minus) / (S - s1minus)
                                           - s_wave(s1plus, m1, m2) * s1plus / (mp * mp - s1plus) / (S - s1plus))
                                + u2 * u2 / (2.0 * u2 + 3.0) / b2
                                          * (s_wave(s2minus, m1, m2) * s2minus / (mp * mp - s2minus) / (S - s2minus)
                                             - s_wave(s2plus, m1, m2) * s2plus / (mp * mp - s2plus) / (S - s2plus)));
        }

        complex<double>
        f_wave(const complex<double> & S, const double & m, const double & q0)
        {
            return f_wave(S, m, m, q0);
        }

        complex<double>
        f_wave(const complex<double> & S, const double & m1, const double & m2, const double & q0)
        {
            const double          mp      = m1 + m2;
            const double          mm      = m1 - m2;
            const double          q0sq    = q0 * q0;
            const double          acbr    = std::cbrt((75.0 * std::sqrt(5.0) - 151.0) / 2.0);
            const double          u1      = -2.0 + acbr - 11.0 / acbr;
            const complex<double> u2      = -2.0 + acbr * complex<double>(-1, std::sqrt(3.0)) / 2.0 + 5.5 / acbr * complex<double>(1, std::sqrt(3.0));
            const complex<double> u3      = -2.0 - acbr * complex<double>(1, std::sqrt(3.0)) / 2.0 - 5.5 / acbr * complex<double>(-1, std::sqrt(3.0));
            const complex<double> a1      = m1 * m1 + m2 * m2 + 2.0 * q0sq * u1;
            const complex<double> b1      = std::sqrt(a1 * a1 - mp * mp * mm * mm);
            const complex<double> a2      = m1 * m1 + m2 * m2 + 2.0 * q0sq * u2;
            const complex<double> b2      = std::sqrt(a2 * a2 - mp * mp * mm * mm);
            const complex<double> a3      = m1 * m1 + m2 * m2 + 2.0 * q0sq * u3;
            const complex<double> b3      = std::sqrt(a3 * a3 - mp * mp * mm * mm);
            const complex<double> s1plus  = a1 + b1;
            const complex<double> s1minus = a1 - b1;
            const complex<double> s2plus  = a2 + b2;
            const complex<double> s2minus = a2 - b2;
            const complex<double> s3plus  = a3 + b3;
            const complex<double> s3minus = a3 - b3;
            const complex<double> zsq     = (S - mp * mp) * (S - mm * mm) / 4.0 / q0sq / S;

            return s_wave(S, m1, m2) * zsq * zsq * zsq / (225.0 + 45.0 * zsq + 6.0 * zsq * zsq + zsq * zsq * zsq)
                   - 2.0 * q0sq * (S - mp * mp)
                             * (u1 * u1 * u1 / (3.0 * u1 * u1 + 12.0 * u1 + 45.0) / b1
                                        * (s_wave(s1minus, m1, m2) * s1minus / (mp * mp - s1minus) / (S - s1minus)
                                           - s_wave(s1plus, m1, m2) * s1plus / (mp * mp - s1plus) / (S - s1plus))
                                + u2 * u2 * u2 / (3.0 * u2 * u2 + 12.0 * u2 + 45.0) / b2
                                          * (s_wave(s2minus, m1, m2) * s2minus / (mp * mp - s2minus) / (S - s2minus)
                                             - s_wave(s2plus, m1, m2) * s2plus / (mp * mp - s2plus) / (S - s2plus))
                                + u3 * u3 * u3 / (3.0 * u3 * u3 + 12.0 * u3 + 45.0) / b3
                                          * (s_wave(s3minus, m1, m2) * s3minus / (mp * mp - s3minus) / (S - s3minus)
                                             - s_wave(s3plus, m1, m2) * s3plus / (mp * mp - s3plus) / (S - s3plus)));
        }
    } // namespace chew_mandelstam
} // namespace eos
