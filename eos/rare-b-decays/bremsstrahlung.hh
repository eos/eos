/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2013 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_BREMSSTRAHLUNG_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_BREMSSTRAHLUNG_HH 1

#include <eos/maths/complex.hh>

namespace eos
{
    struct Bremsstrahlung
    {
        // cf. [AAGW:2001A], Eqs. (30), (31), p. 12
        static complex<double> G_m1(const double & t);
        static complex<double> G_0(const double & t);

        // cf. [AAGW:2001A], Eqs. (28), (29), p. 11
        static complex<double> Deltai_23(const double & s_hat, const double & w, const double & z);
        static complex<double> Deltai_27(const double & s_hat, const double & w, const double & z);

        // cf. [AAGW:2001A], Eqs. (23)-(26), p. 10
        static complex<double> tau_22(const double & s_hat, const double & w, const double & z);
        static complex<double> tau_27(const double & s_hat, const double & w, const double & z);
        static complex<double> tau_28(const double & s_hat, const double & w, const double & z);
        static complex<double> tau_29(const double & s_hat, const double & w, const double & z);
        // cf. [AAGW:2001A], Eqs. (15)-(17), p. 8
        static double tau_78(const double & s_hat);
        static double tau_88(const double & s_hat);
        static double tau_89(const double & s_hat);

        // Integrals of tau_2x from w = s_hat to w = 1, cf [AAGW:2001A], Eq. (22)
        static complex<double> itau_22(const double & s_hat, const double & z);
        static complex<double> itau_27(const double & s_hat, const double & z);
        static complex<double> itau_28(const double & s_hat, const double & z);
        static complex<double> itau_29(const double & s_hat, const double & z);
    };
}

#endif
