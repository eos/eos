/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_MATHS_ANGULAR_INTEGRALS_HH
#define EOS_GUARD_EOS_MATHS_ANGULAR_INTEGRALS_HH 1

namespace eos
{
    // Wigner 3j Symbols
    double wigner_3j(const int & j1, const int & m1, const int & j2, const int & m2, const int & j3, const int & m3);

    // Integral from -1 to 1 over P_l^m(x)
    double legendre_integral(const int & l, const int & m);

    // Integral from -1 to 1 over P_l1^m1(x) * P_l2^m2(x)
    double two_legendre_integral(const int & l1, const int & m1, const int & l2, const int & m2);

    // Integral from -1 to 1 over P_l1^m1(x) * P_l2^m2(x) * P_l3^m3(x)
    double three_legendre_integral(const int & l1, const int & m1, const int & l2, const int & m2, const int & l3, const int & m3);
}

#endif
