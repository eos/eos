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

#ifndef EOS_GUARD_EOS_MATHS_CHEW_MANDELSTAM_HH
#define EOS_GUARD_EOS_MATHS_CHEW_MANDELSTAM_HH 1

#include <eos/maths/complex.hh>

namespace eos
{
    namespace chew_mandelstam
    {
        /*!
         * The S-wave Chew-Mandelstam function for a two-body channel of two
         * particles with masses m1 and m2, evaluated at the Mandelstam
         * variable s.
         *
         * The Chew-Mandelstam function is the once-subtracted dispersive
         * integral of the two-body phase space
         *
         *   rho(s) = sqrt((s - (m1 + m2)^2) * s) / (16 pi s),
         *
         * and is the S-wave (l_orbital = 0) analytic continuation of i * rho(s)
         * used e.g. within the K-matrix formalism to unitarize scattering
         * amplitudes below and above threshold, and onto other Riemann sheets.
         *
         * Only equal masses (m1 == m2) are currently supported, and an
         * InternalError is raised otherwise.
         *
         * @param s     Mandelstam variable at which to evaluate the function.
         * @param m1    Mass of the first particle in the two-body channel.
         * @param m2    Mass of the second particle in the two-body channel.
         */
        complex<double> s_wave(const complex<double> & s, const double & m1, const double & m2);

        /*!
         * The P-wave Chew-Mandelstam function for a two-body channel of two
         * particles with masses m1 and m2, evaluated at the Mandelstam
         * variable s.
         *
         * This is the P-wave (l_orbital = 1) analytic continuation of
         * i * rho(s), including the squared Blatt-Weisskopf form factor
         * for l = 1 (cf. the PDG's resonance review, eq. (50.26)) and the
         * loop-correction term that keeps the amplitude finite at the
         * pseudo-threshold s = mp^2 - 4 q0^2.
         *
         * Only equal masses (m1 == m2) are currently supported, and an
         * InternalError is raised otherwise.
         *
         * @param s     Mandelstam variable at which to evaluate the function.
         * @param m1    Mass of the first particle in the two-body channel.
         * @param m2    Mass of the second particle in the two-body channel.
         * @param q0    Effective momentum scale entering the Blatt-Weisskopf form factor.
         */
        complex<double> p_wave(const complex<double> & s, const double & m1, const double & m2, const double & q0);
    } // namespace chew_mandelstam
} // namespace eos

#endif
