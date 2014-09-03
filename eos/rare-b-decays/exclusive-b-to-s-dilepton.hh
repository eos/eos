/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Christoph Bobeth
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH 1

#include <eos/utils/power_of.hh>

#include <array>

namespace eos
{
    namespace btovll
    {
        struct Amplitudes
        {
            complex<double> a_long_right, a_long_left;
            complex<double> a_perp_right, a_perp_left;
            complex<double> a_par_right, a_par_left;
            complex<double> a_timelike;
            complex<double> a_scalar;
            complex<double> a_par_perp, a_t_long;
            complex<double> a_t_perp, a_long_perp;
            complex<double> a_t_par, a_long_par;
        };

        struct AngularCoefficients
        {
            double j1s, j1c;
            double j2s, j2c;
            double j3;
            double j4;
            double j5;
            double j6s, j6c;
            double j7;
            double j8;
            double j9;
        };

        inline AngularCoefficients array_to_angular_coefficients(const std::array<double, 12> & arr)
        {
            AngularCoefficients a_c = { arr[0], arr[1], arr[2], arr[3], arr[4],  arr[5],
                arr[6], arr[7], arr[8], arr[9], arr[10], arr[11] };

            return a_c;
        }

        inline double decay_width(const AngularCoefficients & a_c)
        {
            // cf. [BHvD2010], p. 6, eq. (2.7)
            return 2.0 * a_c.j1s + a_c.j1c - 1.0 / 3.0 * (2.0 * a_c.j2s + a_c.j2c);
        }

        inline std::array<double, 12> angular_coefficients_array(const Amplitudes & A, const double & s, const double & m_l)
        {
            // cf. [BHvD2010], p. 26, eqs. (A1)-(A11)
            // cf. [BHvD2012], app B, eqs. (B1)-(B12)
            std::array<double, 12> result;

            double z = 4.0 * power_of<2>(m_l) / s;
            double y = m_l / std::sqrt(s);
            double beta2 = 1.0 - z;
            double beta = std::sqrt(beta2);

            // j1s
            result[0] = 3.0 / 4.0 * (
                  (2.0 + beta2) / 4.0 * (norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_par_left) + norm(A.a_par_right))
                  + z * real(A.a_perp_left * conj(A.a_perp_right) + A.a_par_left * conj(A.a_par_right))
                  + 4.0 * beta2 * (norm(A.a_long_perp) + norm(A.a_long_par))
                  + 4.0 * (4.0 - 3.0 * beta2) * (norm(A.a_t_perp) + norm(A.a_t_par))
                  + 8.0 * std::sqrt(2.0) * y * real(
                       (A.a_par_left + A.a_par_right)   * conj(A.a_t_par)
                     + (A.a_perp_left + A.a_perp_right) * conj(A.a_t_perp)
                  )
               );
            // j1c
            result[1] = 3.0 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  + z * (norm(A.a_timelike) + 2.0 * real(A.a_long_left * conj(A.a_long_right)))
                  + beta2 * norm(A.a_scalar)
                  + 8.0 * (2.0 - beta2) * norm(A.a_t_long)
                  + 8.0 * beta2 * norm(A.a_par_perp)
                  + 16.0 * y * real((A.a_long_left + A.a_long_right) * conj(A.a_t_long))
               );
            // j2s
            result[2] = 3.0 * beta2 / 16.0 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_par_left) + norm(A.a_par_right)
                  - 16.0 * (norm(A.a_t_perp) + norm(A.a_t_par) + norm(A.a_long_perp) + norm(A.a_long_par))
               );
            // j2c
            result[3] = -3.0 * beta2 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  - 8.0 * (norm(A.a_t_long) + norm(A.a_par_perp))
               );
            // j3
            result[4] = 3.0 / 8.0 * beta2 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) - norm(A.a_par_left) - norm(A.a_par_right)
                  + 16.0 * (norm(A.a_t_par) - norm(A.a_t_perp) + norm(A.a_long_par) - norm(A.a_long_perp))
               );
            // j4
            result[5] = 3.0 / (4.0 * std::sqrt(2.0)) * beta2 * real(
                  A.a_long_left * conj(A.a_par_left) + A.a_long_right * conj(A.a_par_right)
                  - 8.0 * std::sqrt(2.0) * (A.a_t_long * conj(A.a_t_par) + A.a_par_perp * conj(A.a_long_par))
               );
            // j5
            result[6] = 3.0 * std::sqrt(2.0) / 4.0 * beta * real(
                  A.a_long_left * conj(A.a_perp_left) - A.a_long_right * conj(A.a_perp_right)
                  - 2.0 * std::sqrt(2.0) * A.a_t_par * conj(A.a_scalar)
                  - y * (
                     (A.a_par_left + A.a_par_right) * conj(A.a_scalar)
                     + 4.0 * std::sqrt(2.0) * A.a_long_par * conj(A.a_timelike)
                     - 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_t_perp)
                     - 4.0 * (A.a_perp_left - A.a_perp_right) * conj(A.a_t_long)
                  )
               );
            // j6s
            result[7] = 3.0 / 2.0 * beta * real(
                  A.a_par_left * conj(A.a_perp_left) - A.a_par_right * conj(A.a_perp_right)
                  + 4.0 * std::sqrt(2.0) * y * (
                     (A.a_perp_left - A.a_perp_right) * conj(A.a_t_par)
                     + (A.a_par_left - A.a_par_right) * conj(A.a_t_perp)
                  )
               );
            // j6c
            result[8] = 3.0 * beta * real(
                  2.0 * A.a_t_long * conj(A.a_scalar)
                  + y * (
                     (A.a_long_left + A.a_long_right) * conj(A.a_scalar)
                     + 4.0 * A.a_par_perp * conj(A.a_timelike)
                  )
               );
            // j7
            result[9] = 3.0 * std::sqrt(2.0) / 4.0 * beta * imag(
                  A.a_long_left * conj(A.a_par_left) - A.a_long_right * conj(A.a_par_right)
                  + 2.0 * std::sqrt(2.0) * A.a_t_perp * conj(A.a_scalar)
                  + y * (
                     (A.a_perp_left + A.a_perp_right) * conj(A.a_scalar)
                     + 4.0 * std::sqrt(2.0) * A.a_long_perp * conj(A.a_timelike)
                     + 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_t_par)
                     - 4.0 * (A.a_par_left - A.a_par_right) * conj(A.a_t_long)
                  )
               );
            // j8
            result[10] = 3.0 / 4.0 / std::sqrt(2.0) * beta2 * imag(
                  A.a_long_left * conj(A.a_perp_left) + A.a_long_right * conj(A.a_perp_right)
               );
            // j9
            result[11] = 3.0 / 4.0 * beta2 * imag(
                  conj(A.a_par_left) * A.a_perp_left + conj(A.a_par_right) * A.a_perp_right
               );

            return result;
        }
    }
}

#endif
