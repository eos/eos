/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024-2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_MATHS_OUTER_FUNCTION_HH
#define EOS_GUARD_EOS_MATHS_OUTER_FUNCTION_HH 1

#include <eos/maths/complex.hh>

#include <functional>

namespace eos
{
    /* Computes the outer function of a given function f numerically using the integral representation */
    complex<double> outer(const std::function<complex<double> (const complex<double> &)> & f, complex<double> z, double relative_precision);

}

#endif
