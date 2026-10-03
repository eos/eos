/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Tables for Faddeeva's scaled complex error function, see faddeeva.hh.
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

#ifndef EOS_GUARD_EOS_MATHS_FADDEEVA_CONST_HH
#define EOS_GUARD_EOS_MATHS_FADDEEVA_CONST_HH 1

namespace eos
{
    namespace faddeeva_impl
    {
        namespace imwofx
        {
            extern const double coeffs0[144 * 8];
            extern const double coeffs1[144 * 4];
        } // namespace imwofx

        namespace dimwofx
        {
            inline constexpr double x_min = 2.5; // begin of domain
            inline constexpr double x_max = 12;  // end of domain

            extern const double coeffs0[72 * 8];
            extern const double coeffs1[72 * 4];
        } // namespace dimwofx

        namespace erfcxtab
        {
            extern const double coeffs0[208 * 8];
            extern const double coeffs1[208 * 2];
        } // namespace erfcxtab

        namespace derfcxtab
        {
            inline constexpr double x_min = 0.5; // begin of domain
            inline constexpr double x_max = 12;  // end of domain

            extern const double coeffs0[144 * 8];
            extern const double coeffs1[144 * 4];
        } // namespace derfcxtab

        namespace anchor
        {
            inline constexpr int n_zeros = 10;
            inline constexpr int degree  = 15;

            //! One zero: the lattice tile that holds it, the squared reach, and the
            //! Taylor coefficients at each of the four corners, real part first.
            struct Anchor
            {
                    double xm, xp, ym, yp;
                    double tau2;
                    double a[4][2 * (degree + 1)];
            };

            extern const Anchor anchors[n_zeros];
        } // namespace anchor

        namespace taylor
        {
            inline constexpr int    NTay     = 23; // Taylor truncation order
            inline constexpr double inverseA = 7;  // 1/a, inverse of lattice constant a
            inline constexpr int    nXcover  = 50; // number of tiles per axis

            extern const signed char Tiles[5000];
            extern const double      TaylorCoeffs[2 * 24 * 54];
        } // namespace taylor
    } // namespace faddeeva_impl
} // namespace eos

#endif
