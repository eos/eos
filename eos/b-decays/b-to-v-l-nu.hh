/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2013, 2015 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_V_L_NU_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_V_L_NU_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>

#include <array>

namespace eos
{
    namespace btovlnu
    {
        struct Amplitudes
        {
            complex<double> a_long_left;
            complex<double> a_perp_left;
            complex<double> a_para_left;
            complex<double> a_time_left;

            complex<double> a_paraperp;
            complex<double> a_longpara;
            complex<double> a_timeperp;
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

            AngularCoefficients(const std::array<double, 12> & a) :
                j1s(a[0]),
                j1c(a[1]),
                j2s(a[2]),
                j2c(a[3]),
                j3(a[4]),
                j4(a[5]),
                j5(a[6]),
                j6s(a[7]),
                j6c(a[8]),
                j7(a[9]),
                j8(a[10]),
                j9(a[11])
            {
            }
        };

        inline double decay_width(const AngularCoefficients & a_c)
        {
            // cf. [BHvD:2010A], p. 6, eq. (2.7)
            return 2.0 * a_c.j1s + a_c.j1c - 1.0 / 3.0 * (2.0 * a_c.j2s + a_c.j2c);
        }

        inline std::array<double, 12> angular_coefficients_array(const Amplitudes & a)
        {
            // cf. [BHvD:2010A], p. 26, eqs. (A1)-(A11)
            std::array<double, 12> result;

            // j1s
            result[0] = 3.0 / 16.0 * (3.0 * (norm(a.a_perp_left) + norm(a.a_para_left))
                    + 16.0 * (norm(a.a_longpara) + norm(a.a_timeperp)));
            // j1c
            result[1] = 3.0 / 4.0 * (norm(a.a_long_left) + 2.0 * norm(a.a_time_left)
                    + 8.0 * norm(a.a_paraperp));
            // j2s
            result[2] = 3.0 / 16.0 * (norm(a.a_perp_left) + norm(a.a_para_left)
                    - 16.0 * norm(a.a_longpara) - 16.0 * norm(a.a_timeperp));
            // j2c
            result[3] = -3.0 / 4.0 * (norm(a.a_long_left) - 8.0 * norm(a.a_paraperp));
            // j3
            result[4] = 3.0 / 8.0 * (norm(a.a_perp_left) - norm(a.a_para_left)
                    + 16.0 * norm(a.a_longpara) - 16.0 * norm(a.a_timeperp));
            // j4
            result[5] = 3.0 / (4.0 * sqrt(2.0)) * real(a.a_long_left * conj(a.a_para_left)
                    - 8.0 * sqrt(2.0) * a.a_paraperp * conj(a.a_longpara)
                );
            // j5
            result[6] = 3.0 * sqrt(2.0) / 4.0 * real(a.a_long_left * conj(a.a_perp_left)
                    + 2.0 * sqrt(2.0) * a.a_longpara * conj(a.a_time_left));
            // j6s
            result[7] = 3.0 / 2.0 * real(a.a_para_left * conj(a.a_perp_left));
            // j6c
            result[8] = -6.0 * real(a.a_paraperp * conj(a.a_time_left));
            // j7
            result[9] = 3.0 * sqrt(2.0) / 4.0 * imag(a.a_long_left * conj(a.a_para_left)
                    - 2.0 * sqrt(2.0) * a.a_timeperp * a.a_time_left);
            // j8
            result[10] = 3.0 / 4.0 / sqrt(2.0) * imag(a.a_long_left * conj(a.a_perp_left));
            // j9
            result[11] = 3.0 / 4.0 * imag(a.a_perp_left * conj(a.a_para_left));

            return result;
        }
    }
}

#endif
