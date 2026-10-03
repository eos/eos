/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/* Faddeeva's scaled complex error function
 *      w(z) = exp(-z^2) * erfc(-i*z),
 * nameless function (7.1.3) of Abramowitz & Stegun (1964), also known as the
 * plasma dispersion function.
 *
 * Translated to C++ from the library libcerf,
 *   https://jugit.fz-juelich.de/mlz/lib/cerf
 *
 * Copyright:
 *   (C) 2012      Massachusetts Institute of Technology
 *   (C) 2013-2026 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, 2025
 */

#ifndef EOS_GUARD_EOS_MATHS_FADDEEVA_HH
#define EOS_GUARD_EOS_MATHS_FADDEEVA_HH 1

#include <eos/maths/complex.hh>

namespace eos
{
    //! Returns Faddeeva's scaled complex error function w(z) = exp(-z^2) erfc(-iz).
    complex<double> faddeeva_w(const complex<double> & z);
} // namespace eos

#endif
