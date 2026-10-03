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

#include <eos/maths/faddeeva-const.hh>
#include <eos/maths/faddeeva.hh>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>

#if defined(__FAST_MATH__)
#  error "faddeeva.cc relies on IEEE 754 semantics, which -ffast-math breaks"
#endif
#if FLT_EVAL_METHOD != 0
#  error "faddeeva.cc requires double arithmetic without extended precision"
#endif

namespace eos
{
    namespace faddeeva_impl
    {
        namespace
        {
            /******************************************************************************/
            /*  Exact products                                                            */
            /******************************************************************************/

            /* The exponent Re(-z^2) = (y-x)(x+y), Im(-z^2) = -2xy rounded to one double costs exp(-z^2)
               a relative error of order eps |z|^2, so it is carried as a double-word (Dekker, Numer.
               Math. 18, 224 (1971); Joldes, Muller & Popescu, ACM Trans. Math. Softw. 44, 15 (2017)). */

            //! Returns x*x - h exactly, with h = x*x rounded, so that h + result = x^2.

            double
            sqr_res(const double x, const double h)
            {
#ifdef FP_FAST_FMA
                return std::fma(x, x, -h);
#else
                const double c  = 134217729. * x; // 2^27 + 1: Veltkamp's splitting constant
                const double hx = c - (c - x);    // 26 leading bits of x
                const double tx = x - hx;         // the remaining bits, exactly
                return ((hx * hx - h) + 2. * hx * tx) + tx * tx;
#endif
            }

            //! Returns a*b - p exactly, with p = a*b rounded, so that p + result = a*b.

            double
            prod_res(const double a, const double b, const double p)
            {
#ifdef FP_FAST_FMA
                return std::fma(a, b, -p);
#else
                const double ca = 134217729. * a;
                const double ha = ca - (ca - a);
                const double ta = a - ha;
                const double cb = 134217729. * b;
                const double hb = cb - (cb - b);
                const double tb = b - hb;
                return (((ha * hb - p) + ha * tb) + ta * hb) + ta * tb;
#endif
            }

            /******************************************************************************/
            /*  Real arguments                                                            */
            /******************************************************************************/

            //! Returns exp(x^2), with the exponent carried as a double-word.

            double
            exp_x2(const double x)
            {
                const double h = x * x;
                if (! (h < 750.))
                {
                    return std::exp(h); // overflow, or NaN; the correction would be void
                }
                return std::exp(h) * (1. + sqr_res(x, h));
            }

            //! Returns exp(-x^2), with the exponent carried as a double-word.

            double
            exp_mx2(const double x)
            {
                const double h = x * x;
                if (! (h < 750.))
                {
                    return std::exp(-h); // underflow, or NaN; the correction would be void
                }
                return std::exp(-h) * (1. - sqr_res(x, h));
            }

            /******************************************************************************/
            /*  Complex arguments                                                         */
            /******************************************************************************/

            //! Returns |exp(-z^2)| = exp(y^2-x^2), with the exponent carried as a double-word.

            double
            exp_mz2_mod(const double x, const double y)
            {
                const double th = x * x;
                const double uh = y * y;
                if (! (th <= std::numeric_limits<double>::max() && uh <= std::numeric_limits<double>::max()))
                {
                    // A square overflows, or an argument is NaN. Fall back to the form that
                    // is safe against overflow, though inaccurate; the result is 0 or inf.
                    return std::exp((y - x) * (x + y));
                }
                // 2Sum of u_h and -t_h, plus the two residuals
                const double d  = uh - th;
                const double bb = d - uh;
                const double e  = ((uh - (d - bb)) + (-th - bb)) + (sqr_res(y, uh) - sqr_res(x, th));

                // Renormalize: where the two squares cancel, the residual e is not small
                // compared to their difference d, and the linear correction below would
                // leave a term e^2/2. One more 2Sum restores |s_l| <= eps |s_h| / 2.
                const double sh = d + e;
                const double cc = sh - d;
                const double sl = (d - (sh - cc)) + (e - cc);

                if (sh > 709.7)
                {
                    if (sh > 1100.)
                    {
                        return std::exp(sh); // inf; also keeps sh - k from turning into a NaN
                    }
                    // e^k split off so the correction is applied before the final overflow
                    const double k = std::ceil(sh - 700.);
                    return std::exp(sh - k) * (1. + sl) * std::exp(k);
                }
                if (! (sh > -750.))
                {
                    return std::exp(sh); // underflow, or NaN; the correction would be void
                }
                return std::exp(sh) * (1. + sl);
            }

            //! Returns exp(-z^2)/|exp(-z^2)| = exp(-2ixy), with the phase carried as a double-word.

            complex<double>
            exp_mz2_arg(const double x, const double y)
            {
                const double a  = -2. * x; // exact unless it overflows, which the next test catches
                const double bh = a * y;
                if (! (std::fabs(bh) <= std::numeric_limits<double>::max()))
                {
                    return complex<double>(std::cos(bh), std::sin(bh)); // the result is NaN in any case
                }

                const double bl = prod_res(a, y, bh);
                const double cb = std::cos(bh);
                const double sb = std::sin(bh);
                if (std::fabs(bl) < 1e-8)
                {
                    // The quadratic term of the correction stays below eps/2.
                    return complex<double>(cb - bl * sb, sb + bl * cb);
                }
                // |b_h| > 9e7: the quadratic term of the linear correction would exceed eps.
                const double cl = std::cos(bl);
                const double sl = std::sin(bl);
                return complex<double>(cb * cl - sb * sl, sb * cl + cb * sl);
            }

            //! Returns exp(-z^2) from the library function, in the corners where the
            //! double-word machinery has nothing left to contribute.

            complex<double>
            exp_mz2_lim(const double x, const double y)
            {
                const double a = (y - x) * (x + y); // Re(-z^2) = y^2 - x^2
                const double b = -2. * x * y;       // Im(-z^2)
                if (a == -std::numeric_limits<double>::infinity() && (std::isinf(b) || std::isnan(b)))
                {
                    // The modulus has underflowed to zero, whatever the phase, so C99
                    // Annex G prescribes zero here.
                    return complex<double>(0., 0.);
                }
                return std::exp(complex<double>(a, b));
            }

            //! Returns exp(-z^2) for z = x+iy, with the exponent carried as a double-word.

            complex<double>
            exp_mz2(const double x, const double y)
            {
                const double th = x * x;
                const double uh = y * y;
                if (! (th <= std::numeric_limits<double>::max() && uh <= std::numeric_limits<double>::max()))
                {
                    // A square overflows, or an argument is NaN.
                    return exp_mz2_lim(x, y);
                }

                const double m = exp_mz2_mod(x, y);
                if (! (m > 0. && m <= std::numeric_limits<double>::max()))
                {
                    // The exponential over- or underflows. Nothing is left to correct, and
                    // the limits are delicate, so leave them to the library function;
                    // m > 0 also excludes a NaN.
                    return exp_mz2_lim(x, y);
                }

                // A positive finite real factor, so that the product cannot turn a zero of
                // the phase into a NaN.
                return m * exp_mz2_arg(x, y);
            }

            /******************************************************************************/
            /*  Double-word arithmetic                                                    */
            /******************************************************************************/

            /* The operations of Joldes, Muller and Popescu, ACM Trans. Math. Softw. 44, 15
               (2017), over the exact products above. w_of_z uses them for the phase that
               places a point relative to a zero; the path is rare, so nothing here is
               written for speed. */

            struct DoubleWord
            {
                    double h, l;

                    //! Returns h+l, rounded to a double.
                    double
                    value() const
                    {
                        return h + l;
                    }
            };

            struct ComplexDoubleWord
            {
                    DoubleWord re, im;
            };

            //! Returns a+b exactly.
            DoubleWord
            two_sum(const double a, const double b)
            {
                const double s  = a + b;
                const double bb = s - a;
                return { s, (a - (s - bb)) + (b - bb) };
            }

            //! Returns a+b exactly, given that |a| >= |b|, in three operations instead of six.
            DoubleWord
            quick_two_sum(const double a, const double b)
            {
                const double s = a + b;
                return { s, b - (s - a) };
            }

            //! Returns a*b exactly.
            DoubleWord
            two_prod(const double a, const double b)
            {
                const double p = a * b;
                return { p, prod_res(a, b, p) };
            }

            DoubleWord
            operator- (const DoubleWord x)
            {
                return { -x.h, -x.l };
            }

            DoubleWord
            operator+ (const DoubleWord x, const DoubleWord y)
            {
                const DoubleWord s = two_sum(x.h, y.h);
                const DoubleWord t = two_sum(x.l, y.l);
                const DoubleWord v = quick_two_sum(s.h, s.l + t.h);
                return quick_two_sum(v.h, t.l + v.l);
            }

            DoubleWord
            operator- (const DoubleWord x, const DoubleWord y)
            {
                return x + (-y);
            }

            DoubleWord
            operator* (const DoubleWord x, const DoubleWord y)
            {
                const DoubleWord c = two_prod(x.h, y.h);
                return quick_two_sum(c.h, c.l + (x.h * y.l + x.l * y.h + x.l * y.l));
            }

            //! Returns x*a. Exact in the second factor, so a must be exactly a double;
            //! for the reciprocal of a series index divide by the index instead.
            DoubleWord
            operator* (const DoubleWord x, const double a)
            {
                const DoubleWord c = two_prod(x.h, a);
                return quick_two_sum(c.h, c.l + x.l * a);
            }

            DoubleWord
            operator/ (const DoubleWord x, const double a)
            {
                const double     th = x.h / a;
                const DoubleWord p  = two_prod(th, a);
                return quick_two_sum(th, ((x.h - p.h) - p.l + x.l) / a);
            }

            DoubleWord
            operator/ (const DoubleWord x, const DoubleWord y)
            {
                const double     th = x.h / y.h;
                const DoubleWord p  = x - y * th;
                return quick_two_sum(th, (p.h + p.l) / y.h);
            }

            //! Returns sqrt(x) for positive x.
            DoubleWord
            sqrt(const DoubleWord x)
            {
                const double     s = std::sqrt(x.h);
                const DoubleWord r = x - two_prod(s, s);
                return quick_two_sum(s, (r.h + r.l) / (2 * s));
            }

            //! A term below this no longer moves a double-word, and ends a series.
            const double dw_tiny = 1.9259299443872359e-34; // 2^-112

            //! ln 2, as a double-word.
            const DoubleWord        dw_ln2     = { 6.9314718055994531e-01, 2.3190468138462996e-17 };
            //! pi, as a double-word.
            const DoubleWord        dw_pi      = { 3.1415926535897931e+00, 1.2246467991473532e-16 };
            //! pi/2, as a double-word.
            const DoubleWord        dw_pi_half = { 1.5707963267948966e+00, 6.1232339957367660e-17 };
            //! 1, as a complex double-word.
            const ComplexDoubleWord cdw_one    = {
                DoubleWord{ 1., 0. },
                DoubleWord{ 0., 0. }
            };

            //! Returns ln(x) for positive x. Table-free: the exponent comes out exactly,
            //! and repeated square roots leave a mantissa that the atanh series settles in
            //! a few terms. Nothing therefore bounds the argument.
            DoubleWord
            log(const DoubleWord x)
            {
                int e;
                std::frexp(x.h, &e);
                e            -= 1;                      // so that x / 2^e lies in [1,2)
                DoubleWord y  = x * std::ldexp(1., -e); // exact: a power of two
                for (int i = 0; i < 6; ++i)
                {
                    y = sqrt(y);
                }
                // ln y = 2 atanh(u), u = (y-1)/(y+1), and six square roots have left
                // |u| below 0.011, so eight terms of the series suffice.
                const DoubleWord one = { 1., 0. };
                const DoubleWord u   = (y - one) / (y + one);
                const DoubleWord u2  = u * u;
                DoubleWord       t   = u;
                DoubleWord       s   = u;
                for (int k = 3; k < 60; k += 2)
                {
                    t                     = t * u2;
                    const DoubleWord term = t / k;
                    if (std::fabs(term.value()) < dw_tiny)
                    {
                        break;
                    }
                    s = s + term;
                }
                return dw_ln2 * e + s * 128.; // 2 * 2^6
            }

            //! Returns arctan(t) for |t| <= 1.
            DoubleWord
            atan(DoubleWord t)
            {
                const DoubleWord one = { 1., 0. };
                int              j   = 0;
                // t -> t/(1+sqrt(1+t^2)) halves the argument for one square root
                while (std::fabs(t.h) > 0.02 && j < 12)
                {
                    t = t / (one + sqrt(one + t * t));
                    ++j;
                }
                const DoubleWord t2 = t * t;
                DoubleWord       p  = t;
                DoubleWord       s  = t;
                for (int k = 3; k < 60; k += 2)
                {
                    p                     = p * t2;
                    const DoubleWord term = p / ((k % 4 == 1) ? k : -k);
                    if (std::fabs(term.value()) < dw_tiny)
                    {
                        break;
                    }
                    s = s + term;
                }
                return s * std::ldexp(1., j);
            }

            //! Returns the argument of x+iy, in the correct quadrant, table-free.
            DoubleWord
            dw_atan2(const double y, const double x)
            {
                const DoubleWord ax = { std::fabs(x), 0. };
                const DoubleWord ay = { std::fabs(y), 0. };
                DoubleWord       a  = ay.h <= ax.h ? atan(ay / ax) : dw_pi_half - atan(ax / ay);
                if (x < 0)
                {
                    a = dw_pi - a;
                }
                return y < 0 ? -a : a;
            }

            ComplexDoubleWord
            operator+ (const ComplexDoubleWord a, const ComplexDoubleWord b)
            {
                return { a.re + b.re, a.im + b.im };
            }

            ComplexDoubleWord
            operator- (const ComplexDoubleWord a, const ComplexDoubleWord b)
            {
                return { a.re - b.re, a.im - b.im };
            }

            ComplexDoubleWord
            operator* (const ComplexDoubleWord a, const ComplexDoubleWord b)
            {
                return { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
            }

            ComplexDoubleWord
            operator* (const ComplexDoubleWord a, const double s)
            {
                return { a.re * s, a.im * s };
            }

            ComplexDoubleWord
            operator/ (const ComplexDoubleWord a, const double s)
            {
                return { a.re / s, a.im / s };
            }

            ComplexDoubleWord
            operator/ (const ComplexDoubleWord a, const ComplexDoubleWord b)
            {
                const DoubleWord q = b.re * b.re + b.im * b.im;
                return { (a.re * b.re + a.im * b.im) / q, (a.im * b.re - a.re * b.im) / q };
            }

            //! Returns |a|, as an ordinary double: only its magnitude is ever wanted.
            double
            abs(const ComplexDoubleWord a)
            {
                const double re = a.re.value();
                const double im = a.im.value();
                return std::sqrt(re * re + im * im);
            }

            //! Returns h+l of both parts, rounded to a complex double.
            complex<double>
            to_complex(const ComplexDoubleWord a)
            {
                return complex<double>(a.re.value(), a.im.value());
            }

            //! Returns the piecewise polynomial approximation tabulated in c0 and c1, with 2^M
            //! subdomains per octave, the first octave starting at 2^(j0-1) and the domain at its
            //! subdomain l0. Code taken from https://jugit.fz-juelich.de/mlz/app/ppapp, see
            //! Wuttke & Kleinsorge, ACM Trans. Math. Softw. 52, 13 (2026), doi:10.1145/3805698.
            template <int M, int j0, int l0, int q_stride, int q_terms>
            double
            cheb_approximant(const double * c0, const double * c1, const double x)
            {
                constexpr int loff = (j0 + 1) * (1 << M) + l0;

                int          je;
                const double xm = std::frexp(x, &je);

                const int    ip  = static_cast<int>((1 << (M + 1)) * xm); // index in octave + 2^M
                const int    lij = je * (1 << M) + ip - loff;             // index in lookup table
                const double t   = (1 << (M + 2)) * xm - (1 + 2 * ip);

                const double * const P = c0 + lij * 8;
                const double * const Q = c1 + lij * q_stride;
                double               r = P[0];
                for (int k = 1; k < 8; ++k)
                {
                    r = r * t + P[k];
                }
                for (int k = 0; k < q_terms; ++k)
                {
                    r = r * t + Q[k];
                }
                return r;
            }

            //! Gamma(n+1/2)/pi, the coefficients of the asymptotic expansion of w.
            const double asymptotic[23] = {
                5.6418958354775628e-01, 2.8209479177387814e-01, 4.2314218766081724e-01, // n = 0..2
                1.0578554691520430e+00, 3.7024941420321507e+00, 1.6661223639144676e+01, // n = 3..5
                9.1636730015295726e+01, 5.9563874509942218e+02, 4.4672905882456671e+03, // n = 6..8
                3.7971970000088164e+04, 3.6073371500083758e+05, 3.7877040075087948e+06, // n = 9..11
                4.3558596086351141e+07, 5.4448245107938921e+08, 7.3505130895717545e+09, // n = 12..14
                1.0658243979879044e+11, 1.6520278168812520e+12, 2.7258458978540656e+13, // n = 15..17
                4.7702303212446150e+14, 8.8249260943025370e+15, 1.7208605883889949e+17, // n = 18..20
                3.5277642061974395e+18, 7.5846930433244938e+19,                         // n = 21..22
            };

            //! Returns sum_{k=k0}^{k1} c_k t^(k-k0), where c_k = asymptotic[k], with signs
            //! alternating from + at k0 if requested.
            double
            asymptotic_series(const int k0, const int k1, const double t, const bool alternating)
            {
                const auto c = [&](const int k) { return alternating && (k - k0) % 2 ? -asymptotic[k] : asymptotic[k]; };
                double     p = c(k1);
                for (int k = k1 - 1; k >= k0; --k)
                {
                    p = p * t + c(k);
                }
                return p;
            }

            //! Returns sum_{k<n} c[k] t^(n-1-k), leading coefficient first.
            double
            horner(const double * c, const int n, const double t)
            {
                double p = c[0];
                for (int k = 1; k < n; ++k)
                {
                    p = p * t + c[k];
                }
                return p;
            }

            //! Maclaurin coefficients of Im w(x), from x^25 down to x^1 in steps of x^2.
            const double im_w_maclaurin[13] = {
                5.8461000084165966602290712e-10, -7.30762501052074563638866034e-9, 8.40376876209885782941868884e-8, -8.82395720020380130481012927e-7, 8.38275934019361123956e-6,
                -7.1253454391645686483238e-5,    0.00053440090793734269229,        -0.0034736059015927275001,       0.019104832458760001251,          -0.085971746064420005629,
                0.30090111122547001970,          -0.75225277806367504925,          1.1283791670955125739,
            };

            //! Maclaurin coefficients of erfcx(x), from x^14 down to x^0.
            const double erfcx_maclaurin[15] = {
                1.9841269841269841e-04,  -5.3440090793734269e-04, 1.3888888888888889e-03,  -3.4736059015927274e-03, 8.3333333333333332e-03,
                -1.9104832458760001e-02, 4.1666666666666664e-02,  -8.5971746064419999e-02, 1.6666666666666666e-01,  -3.0090111122547003e-01,
                5.0000000000000000e-01,  -7.5225277806367508e-01, 1.0000000000000000e+00,  -1.1283791670955126e+00, 1.0000000000000000e+00,
            };

            /******************************************************************************/
            /*  im_w_of_x, the scaled Dawson integral 2*dawson(x)/sqrt(pi)                */
            /******************************************************************************/

            double
            im_w_of_x(double x)
            {
                const double ax = std::fabs(x); // very fast

                if (ax < .51)
                {
                    const int n = ax < .003 ? 4 : (ax < .083 ? 7 : (ax < .272 ? 10 : 13));
                    return horner(im_w_maclaurin + 13 - n, n, x * x) * x;
                }

                if (ax < 12.)
                {
                    return std::copysign(cheb_approximant<5, 0, 0, 4, 3>(imwofx::coeffs0, imwofx::coeffs1, ax), x);
                }

                if (ax < 6.9e7)
                {
                    const double r = 1 / x;
                    const int    n = ax < 23.2 ? 10 : (ax < 150 ? 6 : 3);
                    return asymptotic_series(0, n, r * r, false) * r;
                }

                // 1-term expansion, important to avoid overflow
                return std::numbers::inv_sqrtpi / x;
            } // im_w_of_x

            /******************************************************************************/
            /*  d_im_w_of_x, the derivative of im_w_of_x                                  */
            /******************************************************************************/

            /* Written out it is 2/sqrt(pi) - 2 x Im w(x), a difference smaller than either
               of its terms by 1/(2 x^2); forming it that way costs 2 x^2 times the
               relative error of im_w_of_x, some 200 eps at x = 7. Hence this dedicated
               implementation. */

            double
            d_im_w_of_x(double x)
            {
                const double ax = std::fabs(x); // the function is even

                if (ax < dimwofx::x_min)
                {
                    return 2 * (std::numbers::inv_sqrtpi - ax * im_w_of_x(ax));
                }

                if (ax < dimwofx::x_max)
                {
                    return -cheb_approximant<5, 2, 8, 4, 3>(dimwofx::coeffs0, dimwofx::coeffs1, ax);
                }

                // Asymptotic expansion, coefficients 2 Gamma(n+1/2)/pi for n >= 1
                const double r2 = 1 / (ax * ax);
                if (ax < 6.9e7)
                {
                    const int n = ax < 23.2 ? 12 : (ax < 150. ? 8 : 4);
                    return -r2 * (2 * asymptotic_series(1, n, r2, false));
                }
                // 1-term expansion, important to avoid underflow; r2 itself underflows for x > 1.34e154
                return -(std::numbers::inv_sqrtpi / ax) / ax;
            } // d_im_w_of_x

            /******************************************************************************/
            /*  erfcx, the scaled complementary error function exp(x^2) erfc(x)           */
            /******************************************************************************/

            double
            erfcx(double x)
            {
                const double ax = std::fabs(x);

                if (ax < .125)
                {
                    return horner(erfcx_maclaurin, 15, x);
                }

                if (x < 0)
                {
                    if (x < -26.7)
                    {
                        return std::numeric_limits<double>::infinity();
                    }
                    if (x < -6.1)
                    {
                        return 2 * exp_x2(x);
                    }
                    return 2 * exp_x2(x) - cheb_approximant<5, -2, 0, 2, 2>(erfcxtab::coeffs0, erfcxtab::coeffs1, -x);
                }

                if (x < 12)
                {
                    return cheb_approximant<5, -2, 0, 2, 2>(erfcxtab::coeffs0, erfcxtab::coeffs1, x);
                }

                const double r = 1 / x;
                if (x < 6.9e7)
                {
                    const int n = x < 23.2 ? 10 : (x < 150 ? 6 : 3);
                    return asymptotic_series(0, n, r * r, true) * r;
                }
                // 1-term expansion, important to avoid overflow
                return std::numbers::inv_sqrtpi * r;
            } // erfcx

            /******************************************************************************/
            /*  d_erfcx, the derivative of erfcx                                          */
            /******************************************************************************/

            /* Written out it is -2 (1/sqrt(pi) - y erfcx(y)), a difference smaller than
               either of its terms by 1/(2 y^2); forming it that way costs 2 y^2 times the
               relative error of erfcx, some 7000 eps at y = 100. Hence this dedicated
               implementation. */

            double
            d_erfcx(double x)
            {
                if (x < derfcxtab::x_min)
                {
                    return 2 * x * erfcx(x) - 2 * std::numbers::inv_sqrtpi;
                }

                if (x < derfcxtab::x_max)
                {
                    return -cheb_approximant<5, 0, 0, 4, 3>(derfcxtab::coeffs0, derfcxtab::coeffs1, x);
                }

                // Asymptotic expansion, coefficients 2 Gamma(n+1/2)/pi for n >= 1, alternating
                const double r2 = 1 / (x * x);
                if (x < 6.9e7)
                {
                    const int n = x < 23.2 ? 12 : (x < 150. ? 8 : 4);
                    return -r2 * (2 * asymptotic_series(1, n, r2, true));
                }
                // 1-term expansion, important to avoid underflow; r2 itself underflows for x > 1.34e154
                return -(std::numbers::inv_sqrtpi / x) / x;
            } // d_erfcx

            /******************************************************************************/
            /*  asymptotic expansion, summed by a real three-term recurrence              */
            /******************************************************************************/

            /* w(z) = i sum_{n<N} a_n / z^(2n+1) with a_n = Gamma(n+1/2)/pi, for Im z >= 0 and
               |z| >= 7; xs + i*ya is the evaluation point and s = xs^2 + ya^2.

               With eta = 1/z = mu - i*nu, the sum is i*eta*P(eta^2), P(t) = sum a_n t^n. Since
               the a_n are real and eta^2 is a root of the real quadratic t^2 - rho*t - sigma,
               with rho = 2 Re eta^2 and sigma = -|eta^2|^2, synthetic division of P by that
               quadratic (Knuth, TAOCP 2, sect. 4.6.4; the device of Goertzel and Clenshaw)
               evaluates P(eta^2) = gamma + beta_1 eta^2 through the real recurrence for the
               beta_n, and the prefactor i*eta is absorbed as well:
               Re w = nu (a_0 + g+), Im w = mu (a_0 + g-), g+- = (rho +- q) beta_1 + sigma beta_2.
               No complex product ever occurs, the reconstruction is free of cancellation, and
               Im w is strictly proportional to mu, so that its relative accuracy survives
               x -> 0.

               The recurrence is stable because consecutive terms of the sum shrink by
               (n+1/2) q <= 43/98 < 1/2 for q = 1/s <= 1/49 and N <= 23, which keeps the
               beta_n within a factor 2.6 of the a_n; a smaller |z| or a larger N would need
               a new analysis. Two details of the evaluation order carry the error bound:
               q is formed as mu^2 + nu^2, consistent with the computed mu and nu, and not as
               1/s; and a_0 is added last, to the small correction g+-, so that only two
               roundings are incurred at the magnitude of the result. */

            complex<double>
            w_asymptotic(const double xs, const double ya, const double s, const int N)
            {
                const double mu    = xs / s; // eta = 1/z = mu - i*nu
                const double nu    = ya / s;
                const double U     = mu * mu;
                const double V     = nu * nu;
                const double q     = U + V;       // |eta|^2, from U and V for consistency, not 1/s
                const double rho   = 2 * (U - V); // 2 Re eta^2
                const double sigma = -(q * q);    // -|eta^2|^2

                // beta_n = a_n + sigma beta_{n+2} + rho beta_{n+1}, n = N-2..1, beta_{N-1} = a_{N-1}
                double b2 = 0;                 // beta_{n+2}
                double b1 = asymptotic[N - 1]; // beta_{n+1}
                for (int n = N - 2; n >= 1; --n)
                {
                    const double b0 = asymptotic[n] + sigma * b2 + rho * b1;
                    b2              = b1;
                    b1              = b0;
                }
                // now b1 = beta_1, b2 = beta_2
                const double sb2 = sigma * b2;
                const double gp  = sb2 + (rho + q) * b1;
                const double gm  = sb2 + (rho - q) * b1;
                return complex<double>(nu * (asymptotic[0] + gp), mu * (asymptotic[0] + gm));
            }

            /******************************************************************************/
            /*  Re w close to the real axis                                               */
            /******************************************************************************/

            /* Where |y| << |x|, Re w is smaller than |w| by a factor of about y/x. The
               tilewise Taylor approximation bounds the relative error of |w|, so that of
               Re w is amplified by x/y. Here Re w is obtained from an expansion about the
               real axis instead, to fifth order in the small component; Im w, and with it
               |w|, is left to the tiles.

               With w(z) = exp(-z^2) + 2i/sqrt(pi) F(z), where F is Dawson's integral,
               entire and real on the axis, the Taylor coefficients F^(n)(x)/n! are real.
               Scaling them as g_n = 2 F^(n)(x) / (n! sqrt(pi)), so that g_0 = Im w(x),
               only odd n enter Im F:
                 Re w = exp(y^2-x^2) cos 2xy - sum_m (-1)^m g_{2m+1} y^{2m+1}.
               Successive terms shrink by (y/x)^2, so five of them suffice for |y| < 1e-2|x|.
               Run forwards, the recurrence below amplifies rounding errors by about
               (2x^2)^n/n!, but they reach Re w only through the higher terms, suppressed by
               (y/x)^(n-1); Re w stays within a few ulps.

               The coefficients follow from F'(z) = 1 - 2 z F(z):
                 g_{n+1} = (-2 x g_n - 2 g_{n-1}) / (n+1),
                 g_0 = Im w(x),  g_1 = 2/sqrt(pi) - 2 x Im w(x).
               Note that g_1 is the derivative of Im w. It comes from d_im_w_of_x rather
               than from that difference, which would cancel to relative size 1/(2x^2). */

            double
            re_w_near_axis(const double x, const double y)
            {
                double g[10];
                g[0] = im_w_of_x(x);
                g[1] = d_im_w_of_x(x);
                for (int n = 1; n < 9; ++n)
                {
                    g[n + 1] = (-2. * x * g[n] - 2. * g[n - 1]) / (n + 1);
                }

                // exp(y^2-x^2), with x*x carried as an unevaluated sum of two doubles: its
                // single rounding would otherwise cost eps x^2 / 2, that is 25 eps at x = 7
                const double xh = x * x;
                const double xl = sqr_res(x, xh);
                const double e2 = std::exp(y * y - xh) * (1. - xl);

                const double y2 = y * y;
                double       so = 0.; // sum_m (-1)^{m+1} g_{2m+1} (y^2)^m
                for (int m = 4; m >= 0; --m)
                {
                    so = so * y2 + ((m % 2) ? g[2 * m + 1] : -g[2 * m + 1]);
                }

                return e2 * std::cos(2. * x * y) + y * so;
            }

            /******************************************************************************/
            /*  At the innermost zeros: an anchored polynomial                            */
            /******************************************************************************/

            /* The phase form below reaches only from the eleventh zero outward. The ten innermost zeros
               use one Taylor polynomial per sector of the lattice tile holding the zero, anchored at its
               corner, so the relative error stays bounded up to the zero (Wuttke, "Guaranteed relative
               accuracy around a real or complex zero without extended precision"). */

            //! Returns w(x+iy) if that point lies within the reach of one of the ten
            //! innermost zeros. Expects y < 0, x > 0.
            std::optional<complex<double>>
            w_at_inner_zero(const double x, const double y)
            {
                // Im(z^2) = 2xy places z along the string: at the n-th zero it is
                // -2 pi (n - 1/8), up to a term of order |z|^-2
                const double nd = std::round(0.125 - x * y * std::numbers::inv_pi);
                if (! (nd >= 1. && nd <= anchor::n_zeros))
                {
                    return std::nullopt;
                }
                const anchor::Anchor * Z = &anchor::anchors[static_cast<int>(nd) - 1];

                /* One comparison per coordinate picks the sector, and with it the corner it
                   is anchored at. No machine number lies strictly inside the tile, so the
                   comparisons split the lattice cleanly, and the reach is below half the
                   smallest anchor coordinate, so both reductions are exact by Sterbenz. */
                const int    sx = x > Z->xm;
                const int    sy = y > Z->ym;
                const double tx = x - (sx ? Z->xp : Z->xm);
                const double ty = y - (sy ? Z->yp : Z->ym);
                if (! (tx * tx + ty * ty <= Z->tau2))
                {
                    return std::nullopt;
                }

                const double * a  = Z->a[2 * sx + sy];
                double         pr = a[2 * anchor::degree];
                double         pi = a[2 * anchor::degree + 1];
                for (int k = anchor::degree - 1; k >= 0; --k)
                {
                    const double qr = pr * tx - pi * ty;
                    const double qi = pr * ty + pi * tx;
                    pr              = qr + a[2 * k];
                    pi              = qi + a[2 * k + 1];
                }
                return complex<double>(pr, pi);
            }

            /******************************************************************************/
            /*  Near a zero of w: the zero carried in a phase                             */
            /******************************************************************************/

            /* Below the real axis w is 2 exp(-z^2) - w(-z), and at a zero the two terms
               cancel completely, so that no relative accuracy survives. Substituting the
               asymptotic expansion for w(-z) instead of evaluating it turns the reflection
               into

                   w(z) = -A (exp(phi) - 1),        A = -w(-z),
                   phi  = -z^2 + Log(2i sqrt(pi) z) - Log S(z^2) - 2 pi i n,

               with S the series of the expansion and n the index of the nearest zero, at
               which phi vanishes exactly. The cancellation is thereby no longer between two
               computed function values but between -z^2 and the exactly known 2 pi i n,
               which a double-word settles. The amplitude A is what w_asymptotic returned
               anyway, and no table of zeros is needed: the branch of the logarithm is
               immaterial, being absorbed into n. */

            const DoubleWord dw_ln2sqpi = { 1.26551212348464537e+00, 2.83234437198169119e-17 };

            /* 2 pi as a triple-word, good to 2^-164. Two doubles would not do: the phase
               is reduced from a quantity of size |z|^2, so a carrier of 106 bits leaves an
               absolute 2^-106|z|^2, which is the whole of the second term of the error
               bound. Three leave nothing that matters. */
            const double two_pi_0 = 6.28318530717958623e+00;
            const double two_pi_1 = 2.44929359829470641e-16;
            const double two_pi_2 = -5.98953961943667933e-33;

            //! Returns p - 2 pi n, for a p of the size of 2 pi n, with no loss.
            DoubleWord
            reduce_two_pi(const DoubleWord p, const double n)
            {
                const DoubleWord q = two_prod(two_pi_0, n);
                /* p.h and q.h agree to within a few units out of |z|^2, so their difference
                   is exact by Sterbenz's lemma, and the cancellation is over in one
                   subtraction. What is left is small, and the remaining terms are added to
                   it without further loss. */
                return DoubleWord{ p.h - q.h, 0. } + DoubleWord{ p.l, 0. } - DoubleWord{ q.l, 0. } - two_prod(two_pi_1, n) - DoubleWord{ two_pi_2 * n, 0. };
            }

            //! Returns S(zeta), the series of the asymptotic expansion, as a double-word.
            ComplexDoubleWord
            phase_S(const ComplexDoubleWord zeta)
            {
                const ComplexDoubleWord inv  = cdw_one / (zeta * 2.);
                ComplexDoubleWord       term = inv;
                ComplexDoubleWord       s    = cdw_one + term;
                double                  prev = abs(term);
                for (int m = 2; m < 400; ++m)
                {
                    term           = term * inv * (2 * m - 1);
                    const double a = abs(term);
                    if (a < dw_tiny || a > prev) // negligible, or the series turning
                    {
                        break;
                    }
                    s    = s + term;
                    prev = a;
                }
                return s;
            }

            //! Returns Log(1+u) for a small u, by its series.
            ComplexDoubleWord
            phase_log1p(const ComplexDoubleWord u)
            {
                ComplexDoubleWord p = u;
                ComplexDoubleWord r = u;
                for (int k = 2; k < 200; ++k)
                {
                    p                         = p * u;
                    // divide by the index, never multiply by its reciprocal: 1/3 is not a
                    // double, and the rounding would cost half of the double-word
                    const ComplexDoubleWord t = p / ((k % 2) ? k : -k);
                    if (abs(t) < dw_tiny)
                    {
                        break;
                    }
                    r = r + t;
                }
                return r;
            }

            //! Returns exp(phi)-1 for a small phi, by its series.
            ComplexDoubleWord
            phase_expm1(const ComplexDoubleWord phi)
            {
                ComplexDoubleWord term  = phi;
                ComplexDoubleWord r     = phi;
                const double      scale = abs(phi);
                for (int k = 2; k < 60; ++k)
                {
                    term = term * phi / k;
                    if (abs(term) < dw_tiny * scale)
                    {
                        break;
                    }
                    r = r + term;
                }
                return r;
            }

            //! Returns w(x+iy) if that point lies within the disc |phi| < 1 about a zero
            //! of index at least 11. Expects y < 0, x > 0 and |z| >= 7, and A = -w(-z).
            std::optional<complex<double>>
            w_near_zero(const double x, const double y, const complex<double> A)
            {
                /* The caller has already bracketed Re(phi), to within the half ln 2 that
                   the exponent of |z|^2 leaves open, which is tight enough that no further
                   filter earns its keep here. */

                const DoubleWord        x2      = two_prod(x, x);
                const DoubleWord        y2      = two_prod(y, y);
                const ComplexDoubleWord zeta    = { x2 - y2, two_prod(x, y) * 2. };
                const ComplexDoubleWord log_s   = phase_log1p(phase_S(zeta) - cdw_one);
                // Log(2i sqrt(pi) z) = ln(2 sqrt(pi)) + ln|z| + i (pi/2 + arg z)
                const DoubleWord        log_abs = dw_ln2sqpi + log(x2 + y2) / 2.;
                const DoubleWord        arg_z   = dw_pi_half + dw_atan2(y, x);

                /* Only Im(-z^2) is of the size of 2 pi n; the rest of Im(Phi) is of order
                   one. Keeping the two apart until the reduction is what lets the
                   reduction be exact: were they added first, the sum could carry no more
                   than 106 bits, and the small one would be truncated at 2^-106|z|^2. */
                const DoubleWord big_im   = -zeta.im;
                const DoubleWord small_im = arg_z - log_s.im;
                const double     n        = std::round((big_im.value() + small_im.value()) / (2 * std::numbers::pi));
                if (n < 11.)
                {
                    return std::nullopt; // the expansion cannot reach the innermost zeros
                }
                const ComplexDoubleWord phi = { -zeta.re + log_abs - log_s.re, reduce_two_pi(big_im, n) + small_im };
                if (! (abs(phi) < 1.))
                {
                    return std::nullopt; // outside the disc, the reflection is accurate enough
                }

                return -A * to_complex(phase_expm1(phi));
            }

            /******************************************************************************/
            /*  w_of_z, Faddeeva's scaled complex error function                          */
            /******************************************************************************/

            //! Maclaurin coefficients of w(z), from z^19 down to z^1.
            const complex<double> maclaurin[19] = {
                {                       0, -8.8239572002038009e-07 },
                { -2.7557319223985893e-06,                       0 },
                {                       0, +8.3827593401936105e-06 },
                { +2.4801587301587302e-05,                       0 },
                {                       0, -7.1253454391645692e-05 },
                { -1.9841269841269841e-04,                       0 },
                {                       0, +5.3440090793734269e-04 },
                { +1.3888888888888889e-03,                       0 },
                {                       0, -3.4736059015927274e-03 },
                { -8.3333333333333332e-03,                       0 },
                {                       0, +1.9104832458760001e-02 },
                { +4.1666666666666664e-02,                       0 },
                {                       0, -8.5971746064419999e-02 },
                { -1.6666666666666666e-01,                       0 },
                {                       0, +3.0090111122547003e-01 },
                { +5.0000000000000000e-01,                       0 },
                {                       0, -7.5225277806367508e-01 },
                { -1.0000000000000000e+00,                       0 },
                {                       0, +1.1283791670955126e+00 },
            };

            complex<double>
            w_of_z(const complex<double> z)
            {
                const double x  = std::real(z);
                const double xa = std::fabs(x);
                const double y  = std::imag(z);
                const double ya = std::fabs(y);
                const double z2 = xa * xa + y * y;

                // ------------------------------------------------------------------------------
                // Case: the innermost zeros
                // ------------------------------------------------------------------------------

                // The ten innermost zeros (2.4 < |z| < 7.9) use anchored polynomials, the third quadrant by
                // w(-z*) = w(z)*. Their discs lie within 5.4 < |z|^2 < 62.6 and 1.73 < x^2-y^2 < 3.73.

                if (y < 0 && z2 > 5. && z2 < 63.)
                {
                    const double h = 2 * xa * xa - z2; // = x^2 - y^2
                    if (h > 1.7 && h < 3.8)
                    {
                        if (const auto wz = w_at_inner_zero(xa, y))
                        {
                            return x >= 0 ? *wz : std::conj(*wz);
                        }
                    }
                }

                // ------------------------------------------------------------------------------
                // Case |y| << |x|
                // ------------------------------------------------------------------------------

                // Here |Re w| << |Im w|, so Re w comes from the real axis: exp(-x^2) - y d/dx Im w(x).

                if (ya < 1e-9 * xa)
                {
                    const double wi = im_w_of_x(x);
                    const double e2 = xa > 27. ? 0. : exp_mx2(xa);
                    if (ya == 0)
                    {
                        return complex<double>(e2, wi); // also works for x=+-inf
                    }
                    // beyond the asymptotic range of d_im_w_of_x, -y d/dx Im w = y / (sqrt(pi) x^2) without the subnormal 1 / x^2
                    const double re = xa < 6.9e7 ? e2 - y * d_im_w_of_x(x) : (y / x) * (std::numbers::inv_sqrtpi / x);
                    return complex<double>(re, e2 == 0 ? wi : wi - 2 * x * (y * e2)); // also for x = +-inf
                }

                // ------------------------------------------------------------------------------
                // Case |x| << |y|
                // ------------------------------------------------------------------------------

                // Here |Im w| << |Re w|, so Im w comes from the imaginary axis: -x erfcx'(y). For y < 0, w is close to
                // 2 exp(-z^2), whose neglected relative term 2 x^2 y^2 must stay below the last bit.

                if (xa < 1e-9 * ya && (y > 0 || xa == 0 || xa * ya < 1e-8))
                {
                    const double wr = erfcx(y);
                    if (xa == 0)
                    {
                        return complex<double>(wr, 0); // also works for y=+inf
                    }
                    // likewise, -x erfcx'(y) = x / (sqrt(pi) y^2) for large y
                    return complex<double>(wr, y < 6.9e7 ? -x * d_erfcx(y) : (x / y) * (std::numbers::inv_sqrtpi / y));
                }

                // ------------------------------------------------------------------------------
                // Case |z| -> 0: Maclaurin series
                // ------------------------------------------------------------------------------

                if (z2 < .053)
                {
                    const int       degree = z2 < 4e-7 ? 4 : (z2 < .00689 ? 13 : 19);
                    complex<double> ret    = maclaurin[19 - degree];
                    for (int k = 20 - degree; k < 19; ++k)
                    {
                        ret = ret * z + maclaurin[k];
                    }
                    return ret * z + 1.;
                }

                // ------------------------------------------------------------------------------
                // Case |z| -> infty: Asymptotic expansion
                // ------------------------------------------------------------------------------

                if (z2 >= 49)
                {
                    complex<double> ret = 0.;                                   // return value
                    const double    xs  = y < 0 ? -std::real(z) : std::real(z); // compute for -z if y < 0

                    if (z2 > 4.6e16)
                    {
                        // Scale to prevent overflow.
                        if (xa > ya)
                        {
                            const double yax   = ya / xs;
                            const double denom = std::numbers::inv_sqrtpi / (xs + yax * ya);
                            ret                = complex<double>(denom * yax, denom);
                        }
                        else if (std::isinf(ya))
                        {
                            // xa cannot be NaN here, since then z2 were NaN, contradicting z2 >= 49
                            return (y < 0 ? complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()) : complex<double>(0, 0));
                        }
                        else
                        {
                            const double xya   = xs / ya;
                            const double denom = std::numbers::inv_sqrtpi / (xya * xs + ya);
                            ret                = complex<double>(denom, denom * xya);
                        }
                    }
                    else if (z2 > 28000)
                    {
                        ret = w_asymptotic(xs, ya, z2, 4);
                    }
                    else if (z2 > 120)
                    {
                        ret = w_asymptotic(xs, ya, z2, 12);
                    }
                    else
                    {
                        ret = w_asymptotic(xs, ya, z2, 23);
                    }
                    if (y < 0)
                    {
                        // Near a zero the reflection cancels and the phase form takes over. Re(phi) is about
                        // y^2-x^2 + ln(2 sqrt(pi)|z|), and the exponent of |z|^2 pins the logarithm to within
                        // ln(2)/2, so one frexp and one comparison filter for the discs.
                        int zexp;
                        std::frexp(z2, &zexp);                                                            // z2 in [2^(zexp-1), 2^zexp)
                        const double h     = z2 - 2 * xa * xa;                                            // = y^2 - x^2
                        const double rephi = h + 1.2655121234846454 + 0.34657359027997264 * (zexp - 0.5); // ln2/2
                        if (std::fabs(rephi) < 1.25 && z2 <= 4.6e16)
                        {
                            const complex<double> amplitude = x >= 0 ? -ret : -std::conj(ret);
                            if (const auto wz = w_near_zero(xa, y, amplitude))
                            {
                                return x >= 0 ? *wz : std::conj(*wz);
                            }
                        }
                        // Use w(z) = 2.0*exp(-z*z) - w(-z),
                        // but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
                        return 2.0 * exp_mz2(xs, ya) - ret;
                    }
                    else
                    {
                        // The expansion lacks exp(-z^2), which is Re w on the real axis; it matters only where it is not
                        // negligible against Re w ~ y / (sqrt(pi) x^2), that is for small y and |x| < 27.
                        if (ya < 1e-3 && xa < 27)
                        {
                            ret += exp_mz2(x, y);
                        }
                        return ret;
                    }
                }

                // ------------------------------------------------------------------------------
                // Taylor around some z
                // ------------------------------------------------------------------------------

                if (std::isnan(xa))
                {
                    return complex<double>(xa, xa);
                }

                if (std::isnan(y))
                {
                    return complex<double>(y, y);
                }

                const int iTile = 2 * (static_cast<int>(taylor::inverseA * xa) * taylor::nXcover + static_cast<int>(taylor::inverseA * ya));
                const int kP    = taylor::Tiles[iTile];
                assert(kP >= 0); // tiles marked -2 lie in the Maclaurin or asymptotic domain
                const double *        T = &taylor::TaylorCoeffs[2 * (taylor::NTay + 1) * kP];
                const complex<double> dz(xa - T[0], ya - T[1]);
                const int             Nk = taylor::Tiles[iTile + 1];
                complex<double>       ret(T[2 * Nk], T[2 * Nk + 1]);
                for (int k = Nk - 1; k >= 1; --k)
                {
                    ret = ret * dz + complex<double>(T[2 * k], T[2 * k + 1]);
                }

                // Close to the real axis, recompute the small component; see above.
                if (ya < 1e-2 * xa)
                {
                    ret = complex<double>(re_w_near_axis(xa, ya), std::imag(ret));
                }

                if (y < 0)
                {
                    if (x < 0)
                    {
                        return 2.0 * exp_mz2(x, y) - ret;
                    }
                    return 2.0 * exp_mz2(x, y) - std::conj(ret);
                }
                if (x < 0)
                {
                    ret = std::conj(ret);
                }
                return ret;

            } // w_of_z
        } // namespace
    } // namespace faddeeva_impl

    complex<double>
    faddeeva_w(const complex<double> & z)
    {
        return faddeeva_impl::w_of_z(z);
    }
} // namespace eos
