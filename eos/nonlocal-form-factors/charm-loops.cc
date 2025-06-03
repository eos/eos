/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023      Méril Reboud
 * Copyright (c) 2022      Philip Lüghausen
 * Copyright (c) 2010-2025 Danny van Dyk
 * Copyright (c) 2010      Christoph Bobeth
 * Copyright (c) 2010-2011 Christian Wacker
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

#include <eos/maths/power-of.hh>
#include <eos/nonlocal-form-factors/charm-loops.hh>
#include <eos/nonlocal-form-factors/charm-loops-impl.hh>
#include <eos/nonlocal-form-factors/long-distance.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <complex>
#include <cstring>
#include <vector>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    using std::atan;
    using std::complex;
    using std::log;
    using std::pow;
    using std::sqrt;

    complex<double>
    ShortDistanceLowRecoil::c7eff(const double & s, const double & mu, const double & alpha_s, const double & m_b_PS, bool use_nlo,
                const WilsonCoefficients<BToS> & wc)
    {
        // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
        complex<double> result = wc.c7();
        // LO contribution
        result += -1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
        if (use_nlo)
        {
            complex<double> nlo = -1.0 * (
                  wc.c1() * CharmLoops::F17_massless(mu, s, m_b_PS)
                + wc.c2() * CharmLoops::F27_massless(mu, s, m_b_PS)
                + wc.c8() * CharmLoops::F87_massless(mu, s, m_b_PS));
            result += (alpha_s / (4.0 * M_PI)) * nlo;
        }

        return result;
    }

    complex<double>
    ShortDistanceLowRecoil::c8eff(const WilsonCoefficients<BToS> & wc)
    {
        // cf. [BFS2001], below Eq. (26), p. 8
        complex<double> lo = wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

        return wc.c8() + lo;
    }

    complex<double>
    ShortDistanceLowRecoil::c9eff(const double & s, const double & mu, const double & alpha_s, const double & m_b_PS, const double & m_c_MSbar,
                bool use_nlo, bool ccbar_resonance, const complex<double> & lambda_hat_u,
                const WilsonCoefficients<BToS> & wc)
    {
        // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
        // Substitute pole mass by PS mass
        complex<double> c = -2.0 / 27.0 * (8.0 * wc.c1() + 6.0 * wc.c2() - 6.0 * wc.c3() - 8.0 * wc.c4() - 12.0 * wc.c5() - 160.0 * wc.c6());
        complex<double> c_0 = -2.0 / 27.0 * (48.0 * wc.c1() + 36.0 * wc.c2() + 198.0 * wc.c3() - 24.0 * wc.c4() + 1872.0 * wc.c5() - 384.0 * wc.c6());
        complex<double> c_b = +2.0 / 27.0 * (126.0 * wc.c3() + 24.0 * wc.c4() + 1368.0 * wc.c5() + 384.0 * wc.c6());
        complex<double> G0 = -3.0 / 8.0 * ((ccbar_resonance ? LongDistance::g_had_ccbar(s, m_c_MSbar) : CharmLoops::h(mu, s)) + 4.0 / 9.0);
        complex<double> Gb = -3.0 / 8.0 * (CharmLoops::h(mu, s, m_b_PS) + 4.0 / 9.0);

        complex<double> lo = c_b * Gb + c_0 * G0 + c;
        complex<double> nlo_alpha_s = -1.0 * (wc.c1() * CharmLoops::F19_massless(mu, s, m_b_PS)
                                            + wc.c2() * CharmLoops::F29_massless(mu, s, m_b_PS)
                                            + wc.c8() * CharmLoops::F89_massless(s, m_b_PS));
        complex<double> nlo_mc = m_c_MSbar * m_c_MSbar / s * 8.0 *
            ((4.0/9.0 * wc.c1() + 1.0/3.0 * wc.c2()) * (1.0 + lambda_hat_u) + 2.0 * wc.c3() + 20.0 * wc.c5());

        complex<double> result = wc.c9() + lo;
        if ((! ccbar_resonance) && (use_nlo))
            result += (alpha_s / (4.0 * M_PI)) * nlo_alpha_s + nlo_mc;

        return result;
    }

    complex<double>
    CharmLoops::h(const double & mu, const double & s)
    {
        // cf. [BFS2001], Eq. (11), p. 4 in the limit m_q -> 0
        return 4.0 / 9.0 * complex<double>(2.0 / 3.0 + 2.0 * std::log(mu) - std::log(s), M_PI);
    }

    complex<double>
    CharmLoops::h(const double & mu, const double & s, const double & m_q)
    {
        if (m_q < 1e-4)
            return h(mu, s);

        const double z = 4.0 * m_q * m_q / s;
        // treat s smaller than dielectron threshold as zero
        if ((std::abs(s) < 1e-6) || (std::abs(z) < 1e-10))
            return complex<double>(-4.0 / 9.0 * (1.0 + 2.0 * std::log(m_q / mu)), 0.0);

        const double sqrt1z = std::sqrt(std::abs(z - 1.0));
        const double a = 2.0 * std::log(m_q / mu) - 2.0 / 3.0 - z;
        const double b = (2.0 + z) * sqrt1z;
        double rc, ic;

        if ((s > 0) && (z > 1.0))
        {
            ic = 0.0;
            rc = std::atan(1.0 / sqrt1z);
        }
        else if ((s > 0) && (z > 0))
        {
            ic = -M_PI / 2.0;
            rc = std::log((1.0 + sqrt1z) / std::sqrt(z));
        }
        else if (s < 0) // we use [KMPW2010], Eq. (12), p. 7
        {
            ic = 0.0;
            // note that our prefactor b varies from eq. (12) by a factor of 2.
            // therefore, c = 0.5 * log(...)
            //
            rc = -0.5 * std::log((sqrt1z - 1.0) / (sqrt1z + 1.0));
        }
        else
        {
            throw InternalError("CharmLoops::h not prepared for its arguments");
        }

        // cf. [BFS2001], Eq. (11), p. 4
        return -4.0 / 9.0 * (a + b * complex<double>(rc, ic));
    }

    complex<double>
    CharmLoops::A(const double & mu, const double & s, const double & m_b)
    {
        /* in the limit s -> 0 all terms vanish, except for the mu-dependent log term. */
        if (std::abs(s) <= 1e-6) // treat s smaller than dielectron threshold as zero
        {
            return -104.0 / 243.0 * 2.0 * log(m_b / mu);
        }

        /* cf. [S2004], Eq. (29), p. 8
         * We have three different cases for the evaluation of the formula depending on s_hat.
         * 1. s_hat < 1: We can use the formula without modification
         * 2. s_hat > 1: 1 - s_hat is negative. We have to take care, because in this regime
         *    the logarithm and the dilogarithm have a branch cut. The real part ist continuous.
         *    And for small epsilon > 0 is Im(dilog(s +- i*epsilon)) = +-pi*ln(s) and
         *    Im(log(1 - s -+ i*epsilon)) = -+pi. So regardless of the epsilon chosen the final result
         *    does not change. If we use compatible definitions for the branch cuts (as in clib and gsl),
         *    we don't need to specify epsilon.
         * 3. s_hat = 1: The fomula cannot be used since denom = 0, instead we use a taylor approximation.
         */

        double s_hat = s / m_b / m_b, s_hat2 = s_hat * s_hat, denom = 1 - s_hat;
        double ln = log(s_hat), ln2 = ln * ln;

        double z = 4.0 / s_hat, sqrt1z = sqrt(z - 1.0);

        double a = -104.0 / 243.0 * 2.0 * log(m_b / mu);

        if (std::abs(s_hat - 1.0) < 1e-2)
        {
            // We use taylor approximation with a maximum error of 4*10^-8, as the exact expression is numerical instable.
            static const std::complex<double> c0((997.0 + 18.0 * std::sqrt(3.0) * M_PI) / 1458.0, 64.0 / 243.0 * M_PI);
            static const std::complex<double> c1((215.0 +  9.0 * std::sqrt(3.0) * M_PI) / 1215.0, -1.0 / 27.0 * M_PI);
            static const std::complex<double> c2(( 95.0 + 12.0 * std::sqrt(3.0) * M_PI) / 2430.0, -7.0 / 405.0 * M_PI);

            return a + c0 + c1 * denom + c2 * power_of<2>(denom);
        }

        complex<double> ln1s = std::log(std::complex<double> (1.0 - s_hat, 0.0));

        // the dilogarithm function expects the radius and angle of z, but returns the real and imaginary part
        gsl_sf_result li_2s_real, li_2s_imag;
        gsl_sf_complex_dilog_e(s_hat, 0, &li_2s_real, &li_2s_imag);
        complex<double> li_2s(li_2s_real.val, li_2s_imag.val);

        complex<double> b = +4.0 * s_hat / 27.0 / denom * (li_2s + ln * ln1s);

        complex<double> c = +1.0 / 729.0 / power_of<2>(denom) * complex<double>(
                6.0 * s_hat * (29.0 - 47.0 * s_hat) * ln + 785.0 - 1600.0 * s_hat + 833.0 * s_hat2,
                6.0 * M_PI * (20.0 - 49.0 * s_hat + 47.0 * s_hat2));

        // identity: arccot(x) = pi / 2 - arctan(x)
        complex<double> d = -2.0 / 243.0 / power_of<3>(denom) * complex<double>(
                2.0 * sqrt1z * (-4.0 + 9.0 * s_hat - 15.0 * s_hat2 + 4.0 * s_hat * s_hat2) * (M_PI / 2.0 - atan(sqrt1z))
                + 9.0 * s_hat * s_hat2 * ln2,
                18.0 * M_PI * s_hat * (1.0 - 2.0 * s_hat) * ln);

        double e = +2.0 * s_hat / 243.0 / power_of<4>(denom) * (36.0 * power_of<2>(M_PI / 2.0 - atan(sqrt1z)) +
                M_PI * M_PI * (-4.0 + 9.0 * s_hat - 9.0 * s_hat2 + 3.0 * s_hat * s_hat2));

        return a + b + c + d + e;
    }

    complex<double>
    CharmLoops::B(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (30), p. 8
        // See remarks in CharmLoops::A
        double s_hat = s / m_b / m_b, s_hat2 = s_hat * s_hat, denom = 1 - s_hat;
        double ln = log(s_hat), ln2 = ln * ln;
        double z = 4.0 / s_hat, sqrt1z = sqrt(z - 1.0);
        double lnmu = 2.0 * log(m_b / mu);

        complex<double> x1(0.5, 0.5 * sqrt1z), x2(0.5, -0.5 * sqrt1z), x3(0.5, 0.5 / sqrt1z), x4(0.5, -0.5 / sqrt1z);
        complex<double> lx1 = log(x1), lx2 = log(x2), lx3 = log(x3), lx4 = log(x4);

        complex<double> a = 8.0 / 243.0 / s_hat * (complex<double>(4.0 - 34.0 * s_hat, -17.0 * M_PI * s_hat) * lnmu
                + 8.0 * s_hat * lnmu * lnmu + 17.0 * s_hat * ln * lnmu);

        if (std::abs(s_hat - 1.0) < 1e-2)
        {
            /* We use a taylor approximation around s_hat = 1 with a maximum error of 5e-7, as the exact expression is
             * numerical instable. The terms containing masses are implemented exactly.
             * The taylor polynomial reads
             *                1 / 2187 * (287 - 6 * (-366 * i + sqrt(3)) * pi + (-342 + 8 * i * sqrt(3)) * pi^2 +
             *                + (9 - 3 * i * sqrt(3)) * Polygamma(1, 1/6) + (9 + 3 * i * sqrt(3)) * Polygamma(1, 1/3) +
             *                + 3 * i * (3 * i + sqrt(3)) * Polygamma(1, 2/3) + (-9 - 3 * i * sqrt(3)) * Polygamma(1, 5/6))
             * + (s - 1)   * (1 / 32805 * (2 * pi (-7155 * i + 573 * sqrt(3) + (-540 - 80 * i * sqrt(3)) * pi) +
             *                60 * (-423 + i (3 * i + sqrt(3)) * Polygamma(1, 1/6) + (-3 - I sqrt(3)) * Polygamma(1, 1/3) +
             *                (3 - i sqrt(3)) Polygamma(1, 2/3) + (3 + i * sqrt(3)) * Polygamma(1, 5/6)))
             * + (s - 1)^2 * (-1 / 196830 * (-1)^(5/6) (2 * pi * (-16587 + 6999 * i * sqrt(3) + 20 * (123 * i + 67 * sqrt(3)) * pi) +
             *                15 * (2925 * (i + sqrt(3)) + 56 * sqrt(3) * Polygamma(1, 1/6) + 28 * (3 * i + sqrt(3)) * Polygamma(1, 1/3) -
             *                56 * sqrt(3) * Polygamma(1, 2/3) - 28 * (3 * i + sqrt(3)) * Polygamma(1, 5/6)))
             * But we use the numerical values.
             */
            return std::complex<double>(-1.2534705628994441, 3.1545210184193809)
                    + std::complex<double>(-1.1399966466176837, -1.3704066719362884) * (s_hat - 1.0)
                    + std::complex<double>(0.77575942579740349, 0.59987612809286587) * power_of<2>(s_hat - 1.0)
                    + a - 16.0 / 243.0 * (2.0 + s_hat) / s_hat * sqrt1z * lnmu * (M_PI / 2.0 - atan(sqrt1z));
        }

        // calculate Li_2(-x_2 / x_1)
        gsl_sf_result li_2x2x1_real, li_2x2x1_imag;
        gsl_sf_complex_dilog_e(1.0, arg(-x2 / x1), &li_2x2x1_real, &li_2x2x1_imag);
        complex<double> li_2x2x1(li_2x2x1_real.val, li_2x2x1_imag.val);

        complex<double> ln1s = std::log(std::complex<double> (1.0 - s_hat, 0.0));

        gsl_sf_result li_2s_real, li_2s_imag;
        gsl_sf_complex_dilog_e(s_hat, 0, &li_2s_real, &li_2s_imag);
        complex<double> li_2s(li_2s_real.val, li_2s_imag.val);

        complex<double> b = (2.0 + s_hat) * sqrt1z / 729.0 / s_hat * (
                    -48.0 * lnmu * (M_PI / 2.0 - atan(sqrt1z)) - 18.0 * M_PI * 2.0 * log(sqrt1z) - 12.0 * M_PI * (2.0 * lx1 + lx3 + lx4)
                    + complex<double>(0.0, 1.0) * (
                        3.0 * power_of<2>(log(z - 1.0)) - 24.0 * li_2x2x1 - 5.0 * M_PI * M_PI
                        + 6.0 * (-9.0 * power_of<2>(lx1) + power_of<2>(lx2) - 2.0 * power_of<2>(lx4) + 6.0 * lx1 * lx2
                        - 4.0 * lx1 * lx3 + 8.0 * lx1 * lx4)));

        complex<double> c = -2.0 / 243.0 / s_hat / denom * (
                4.0 * s_hat * (-8.0 + 17.0 * s_hat) * (li_2s + ln * ln1s)
                + 3.0 * (2.0 + s_hat) * (3.0 - s_hat) * power_of<2>(lx2 - lx1)
                + 12.0 * M_PI * (-6.0 - s_hat + s_hat2) * (M_PI / 2.0 - atan(sqrt1z)));

        complex<double> d = 2.0 / (2187.0 * s_hat * power_of<2>(denom)) * complex<double>(
                -18.0 * s_hat * (120.0 - 211.0 * s_hat + 73.0 * s_hat2) * ln
                -288.0 - 8.0 * s_hat + 934.0 * s_hat2 - 692.0 * s_hat * s_hat2,
                18.0 * M_PI * s_hat * (82.0 - 173.0 * s_hat + 73.0 * s_hat2));

        complex<double> e = -4.0 / (243.0 * s_hat * power_of<3>(denom)) * std::complex<double>(
                -2.0 * sqrt1z * (4.0 - 3.0 * s_hat - 18.0 * s_hat2 + 16.0 * s_hat * s_hat2 - 5.0 * s_hat2 * s_hat2) * (M_PI / 2 - atan(sqrt1z))
                -9.0 * s_hat * s_hat2 * ln2,
                2.0 * M_PI * s_hat * (8.0 - 33.0 * s_hat + 51.0 * s_hat2 - 17.0 * s_hat * s_hat2) * ln);

        complex<double> f = 2.0 / (729.0 * s_hat * power_of<4>(denom)) * (72.0 * (3.0 - 8.0 * s_hat + 2.0 * s_hat2)
                * power_of<2>(M_PI / 2.0 - atan(sqrt1z))
                - M_PI * M_PI * (54.0 - 53.0 * s_hat - 286.0 * s_hat2 + 612.0 * s_hat * s_hat2 - 446 * s_hat2 * s_hat2 + 113.0 * s_hat2 * s_hat2 * s_hat));

        return a + b + c + d + e + f;
    }

    complex<double>
    CharmLoops::C(const double & mu, const double & s)
    {
        static const double zeta3 = 1.20206;

        // cf. [S2004], Eq. (31), p. 9
        return complex<double>(16.0 / 81.0 * log(mu * mu / s) + 428.0 / 243.0 - 64.0 / 27.0 * zeta3,
                16.0 / 81.0 * M_PI);
    }

    /* Two-Loop functions for massless quarks from[S2004], suitable for up-quark loops */
    complex<double>
    CharmLoops::F17_massless(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (22), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return -A(mu, s, m_b);
    }

    complex<double>
    CharmLoops::F19_massless(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (24), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return -B(mu, s, m_b) - 4.0 * C(mu, s);
    }

    complex<double>
    CharmLoops::F27_massless(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (23), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return 6.0 * A(mu, s, m_b);
    }

    complex<double>
    CharmLoops::F29_massless(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (25), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return 6.0 * B(mu, s, m_b) - 3.0 * C(mu, s);
    }

    /* Two-Loop functions for charm-quark loops */
    // cf. [AAGW2001], Eq. (56), p. 20
    complex<double>
    CharmLoops::F17_massive(const double & mu, const double & s, const double & m_b, const double & m_c)
    {
        // cf. [ABGW2001], Appendix B, pp. 34-38
        static double kap1700[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-1.14266, -0.517135}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-2.20356, 1.59186}, {-5.21743, 1.86168}, {0.592593, 3.72337}, {0.395062, 0}, {0, 0}},
            {{1.86366, -3.06235}, {-4.66347, 0}, {0, 3.72337}, {0.395062, 0}, {0, 0}},
            {{-1.21131, 2.89595}, {2.99588, -2.48225}, {-4.14815, 0}, {0, 0}, {0, 0}}
        };

        static double kap1710[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-2.07503, 1.39626}, {-0.444444, 0.930842}, {0, 0}, {0, 0}, {0, 0}},
            {{-25.9259, 5.78065}, {-3.40101, 13.0318}, {-4.4917, 3.72337}, {0.395062, 0}, {-0.395062, 0}},
            {{11.4229, -15.2375}, {-34.0806, 11.1701}, {10.3704, 18.6168}, {2.37037, 0}, {0, 0}},
            {{11.7509, 15.6984}, {18.9564, -24.8225}, {-14.6173, 0}, {0, 0}, {0, 0}}
        };

        static double kap1711[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0164609, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{1.03704, 0.930842}, {0.592593, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-4.66347, 0}, {0, 7.44674}, {2.37037, 0}, {0, 0}, {0, 0}},
            {{6.73754, 1.86168}, {1.18519, -7.44674}, {-2.37037, 0}, {0, 0}, {0, 0}}
        };

        static double kap1720[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.00555556, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-19.4691, 1.59019}, {-11.6779, 0.930842}, {-2.96296, 0}, {-0.395062, 0}, {0, 0}},
            {{-90.4953, 14.7788}, {14.9329, 22.3402}, {-24.438, 3.72337}, {1.18519, 0}, {-1.18519, 0}},
            {{23.8816, -32.8021}, {-82.7915, 39.0954}, {32.2963, 44.6804}, {5.92593, 0}, {0, 0}},
            {{38.1415, 34.8683}, {38.6436, -80.673}, {-41.5802, 0}, {0, 0}, {0, 0}}
        };

        static double kap1721[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0164609, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{2.37037, 1.86168}, {1.18519, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-13.9904, 3.72337}, {2.37037, 22.3402}, {7.11111, 0}, {0, 0}, {0, 0}},
            {{27.5428, 3.72337}, {2.37037, -29.787}, {-9.48148, 0}, {0, 0}, {0, 0}}
        };

        static double kap1730[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.00010778, 0.00258567}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.946811, -0.0258567}, {0.488889, 0}, {0.0987654, 0}, {0, 0}, {0, 0}},
            {{-41.9952, 1.63673}, {-30.2091, 0.930842}, {-6.22222, 0}, {-1.18519, 0}, {0, 0}},
            {{-189.354, 25.8196}, {42.6566, 31.0281}, {-57.765, 3.72337}, {2.76543, 0}, {-2.37037, 0}},
            {{45.1784, -52.4207}, {-145.181, 88.7403}, {70.9136, 81.9141}, {11.0617, 0}, {0, 0}},
            {{77.3602, 54.2499}, {58.4491, -184.927}, {-96.0988, 0}, {0, 0}, {0, 0}}
        };

        static double kap1731[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0164609, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{3.85185, 2.79253}, {1.77778, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-27.3882, 13.0318}, {8.2963, 44.6804}, {14.2222, 0}, {0, 0}, {0, 0}},
            {{69.4495, 1.86168}, {1.18519, -74.4674}, {-23.7037, 0}, {0, 0}, {0, 0}}
        };

        double m_c_hat = m_c / m_b, z = power_of<2>(m_c_hat);
        double s_hat = s / power_of<2>(m_b);

        complex<double> log_s_hat = { std::log(std::abs(s_hat)), 0.0 };
        if ((0.0 < s_hat) && (s_hat <= 0.45))
        {
            log_s_hat.imag(0.0);
        }
        else if ((-0.45 <= s_hat) && (s_hat <= -0.00))
        {
            log_s_hat.imag(+M_PI);
        }
        else
        {
            throw InternalError("CharmLoop::F17_massive used outside its domain of validity, s_hat = " + stringify(s_hat));
        }

        const double rho17[4] = {
            1.94955 * power_of<3>(m_c_hat), 11.6973 * m_c_hat, 70.1839 * m_c_hat, -3.8991 / m_c_hat + 159.863 * m_c_hat
        };

        // real part
        complex<double> r = -208.0 / 243.0 * log(mu / m_b);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 4 ; m++)
                r = r + kap1700[l][m][0] * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap1710[l][m][0] * s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1711[l][m][0] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 2 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap1720[l][m][0] * s_hat * s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1721[l][m][0] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap1730[l][m][0] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1731[l][m][0] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 0 ; l < 4; l++)
            r = r + rho17[l] * pow(s_hat, l);

        // imaginary part
        complex<double> i = 0.0;

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1700[l][m][1] * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1710[l][m][1] * s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1711[l][m][1] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1720[l][m][1] * s_hat * s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1721[l][m][1] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1730[l][m][1] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_c_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1731[l][m][1] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_c_hat), m);

        return r + complex<double>(0.0, 1.0) * i;
    }

    namespace impl
    {
        complex<double> f27_0(const double & mu, const double & m_b, const double & m_q)
        {
            static long double kap2700[7][5][2] = {
                {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
                {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
                {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
                {{6.85597, 3.10281}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
                {{13.2214, -9.55118}, {31.3046, -11.1701}, {-3.55556, -22.3402}, {-2.37037, 0}, {0, 0}},
                {{-11.182, 18.3741}, {27.9808, 0}, {0, -22.3402}, {-2.37037, 0}, {0, 0}},
                {{7.26787, -17.3757}, {-17.9753, 14.8935}, {24.8889, 0}, {0, 0}, {0, 0}}
            };

            const double m_q_hat = m_q / m_b, z = power_of<2>(m_q_hat);

            const long double rho27 = -11.6973 * power_of<3>(m_q_hat);

            // real part
            double r = 416.0 / 81.0 * log(mu / m_b);

            for (int l = 3 ; l < 7 ; l++)
                for (int m = 0 ; m < 4 ; m++)
                    r = r + kap2700[l][m][0] * pow(z, l-3) * pow(log(m_q_hat), m);

            r = r + rho27 ;

            // imaginary part
            double i = 0.0;

            for (int l = 3 ; l < 7 ; l++)
                for (int m = 0 ; m < 3 ; m++)
                    i = i + kap2700[l][m][1] * pow(z, l-3) * pow(log(m_q_hat), m);

            return complex<double>(r, i);
        }
    }

    // cf. [AAGW2001], Eq. (56), p. 20
    complex<double>
    CharmLoops::F27_massive(const double & mu, const double & s, const double & m_b, const double & m_q)
    {
        // cf. [ABGW2001], Appendix B, pp. 34-38
        static double kap2700[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{6.85597, 3.10281}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{13.2214, -9.55118}, {31.3046, -11.1701}, {-3.55556, -22.3402}, {-2.37037, 0}, {0, 0}},
            {{-11.182, 18.3741}, {27.9808, 0}, {0, -22.3402}, {-2.37037, 0}, {0, 0}},
            {{7.26787, -17.3757}, {-17.9753, 14.8935}, {24.8889, 0}, {0, 0}, {0, 0}}
        };

        static double kap2710[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{12.4502, -8.37758}, {2.66667, -5.58505}, {0, 0}, {0, 0}, {0, 0}},
            {{155.555, -34.6839}, {20.4061, -78.1908}, {26.9502, -22.3402}, {-2.37037, 0}, {2.37037, 0}},
            {{-68.5374, 91.4251}, {204.484, -67.0206}, {-62.2222, -111.701}, {-14.2222, 0}, {0, 0}},
            {{-70.5057, -94.1903}, {-113.738, 148.935}, {87.7037, 0}, {0, 0}, {0, 0}}
        };

        static double kap2711[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.0987654, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-6.22222, -5.58505}, {-3.55556, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{27.9808, 0}, {0, -44.6804}, {-14.2222, 0}, {0, 0}, {0, 0}},
            {{-40.4253, -11.1701}, {-7.11111, 44.6804}, {14.2222, 0}, {0, 0}, {0, 0}}
        };

        static double kap2720[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0333333, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{116.815, -9.54113}, {70.0677, -5.58505}, {17.7778, 0}, {2.37037, 0}, {0, 0}},
            {{542.972, -88.6728}, {-89.5971, -134.041}, {146.628, -22.3402}, {-7.11111, 0}, {7.11111, 0}},
            {{-143.29, 196.813}, {496.749, -234.572}, {-193.778, -268.083}, {-35.5556, 0}, {0, 0}},
            {{-228.849, -209.21}, {-231.862, 484.038}, {249.481, 0}, {0, 0}, {0, 0}}
        };

        static double kap2721[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.0987654, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-14.2222, -11.1701}, {-7.11111, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{83.9424, -22.3402}, {-14.2222, -134.041}, {-42.6667, 0}, {0, 0}, {0, 0}},
            {{-165.257, -22.3402}, {-14.2222, 178.722}, {56.8889, 0}, {0, 0}, {0, 0}}
        };

        static double kap2730[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.000646678, -0.015514}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-5.68087, 0.15514}, {-2.93333, 0}, {-0.592593, 0}, {0, 0}, {0, 0}},
            {{251.971, -9.82039}, {181.255, -5.58505}, {37.3333, 0}, {7.11111, 0}, {0, 0}},
            {{1136.13, -154.918}, {-255.94, -186.168}, {346.59, -22.3402}, {-16.5926, 0}, {14.2222, 0}},
            {{-271.07, 314.524}, {871.089, -532.442}, {-425.481, -491.485}, {-66.3704, 0}, {0, 0}},
            {{-464.161, -325.499}, {-350.695, 1109.56}, {576.593, 0}, {0, 0}, {0, 0}}
        };

        static double kap2731[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.0987654, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-23.1111, -16.7552}, {-10.6667, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{164.329, -78.1908}, {-49.7778, -268.083}, {-85.3333, 0}, {0, 0}, {0, 0}},
            {{-416.697, -11.1701}, {-7.11111, 446.804}, {142.222, 0}, {0, 0}, {0, 0}}
        };

        double m_q_hat = m_q / m_b, z = power_of<2>(m_q_hat);
        double s_hat = s / m_b / m_b;

        const double rho27[4] = {
            -11.6973 * power_of<3>(m_q_hat), -70.1839 * m_q_hat, -421.103 * m_q_hat, 23.3946 / m_q_hat - 959.179 * m_q_hat
        };

        if (s_hat == 0)
        {
            return impl::f27_0(mu, m_b, m_q);
        }

        complex<double> log_s_hat = { std::log(std::abs(s_hat)), 0.0 };
        if ((0.0 < s_hat) && (s_hat <= 0.45))
        {
            log_s_hat.imag(0.0);
        }
        else if ((-0.45 <= s_hat) && (s_hat <= -0.00))
        {
            log_s_hat.imag(+M_PI);
        }
        else
        {
            throw InternalError("CharmLoop::F27_massive used outside its domain of validity, s_hat = " + stringify(s_hat));
        }

        // real part
        complex<double> r = 416.0 / 81.0 * log(mu / m_b);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 4 ; m++)
                r = r + kap2700[l][m][0] * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap2710[l][m][0] * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2711[l][m][0] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 2 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap2720[l][m][0] * s_hat * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2721[l][m][0] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap2730[l][m][0] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2731[l][m][0] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0 ; l < 4; l++)
            r = r + rho27[l] * pow(s_hat, l);

        // imaginary part
        complex<double> i = 0.0;

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2700[l][m][1] * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2710[l][m][1] * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2711[l][m][1] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2720[l][m][1] * s_hat * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2721[l][m][1] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2730[l][m][1] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2731[l][m][1] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        return r + complex<double>(0.0, 1.0) * i;
    }

    // cf. [AAGW2001], Eq. (54), p. 19
    complex<double>
    CharmLoops::F19_massive(const double & mu, const double & s, const double & m_b, const double & m_q)
    {
        // F19(s) diverges for s -> 0. However, s * F19(s) -> 0 for s -> 0.
        if (abs(s) < 1e-6) // allow for s = 1e-6, corresponding roughly to the dielectron threshold
            throw InternalError("CharmLoops::F19_massive: F19 diverges for s -> 0. Check that F19 enters via 's * F19(s)' and replace by zero.");

        // cf. [ABGW2001], Appendix B, pp. 34-38
        static double kap1900[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-4.61812, 3.67166}, {5.62963, 1.86168}, {0, 0}, {0, 0}, {0, 0}},
            {{14.4621, -16.2155}, {9.59321, -11.1701}, {-1.18519, -7.44674}, {-0.790123, 0}, {0, 0}},
            {{-16.0864, 26.7517}, {54.2439, -14.8935}, {-15.4074, -29.787}, {-3.95062, 0}, {0, 0}},
            {{-14.73, -23.6892}, {-28.5761, 34.7514}, {20.1481, 0}, {0, 0}, {0, 0}}
        };

        static double kap1901[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0493827, -0.103427}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.592593, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{4.95977, -1.86168}, {-1.18519, -7.44674}, {-2.37037, 0}, {0, 0}, {0, 0}},
            {{-9.20287, -1.65483}, {-1.0535, 9.92898}, {3.16049, 0}, {0, 0}, {0, 0}}
        };

        static double kap1910[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-2.48507, -0.186168}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{4.47441, -0.310281}, {1.48148, -1.86168}, {0, 0}, {0, 0}, {0, 0}},
            {{71.3855, -30.7987}, {8.47677, -33.5103}, {12.5389, -7.44674}, {-0.790123, 0}, {0.790123, 0}},
            {{-18.1301, 66.1439}, {149.596, -67.0206}, {-49.1852, -81.9141}, {-11.0617, 0}, {0, 0}},
            {{-72.89, -63.7828}, {-68.135, 134.041}, {63.6049, 0}, {0, 0}, {0, 0}}
        };

        static double kap1911[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-2.66667, -1.86168}, {-1.18519, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{18.6539, -7.44674}, {-4.74074, -29.787}, {-9.48148, 0}, {0, 0}, {0, 0}},
            {{-41.6104, -3.72337}, {-2.37037, 44.6804}, {14.2222, 0}, {0, 0}, {0, 0}}
        };

        static double kap1920[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.403158, -0.0199466}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0613169, 0.0620562}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{37.1282, -1.36524}, {22.0621, -1.86168}, {5.33333, 0}, {0.790123, 0}, {0, 0}},
            {{212.74, -52.2081}, {-21.9215, -52.1272}, {57.1724, -7.44674}, {-2.37037, 0}, {2.37037, 0}},
            {{-44.6829, 108.713}, {272.015, -163.828}, {-119.111, -156.382}, {-21.3333, 0}, {0, 0}},
            {{-137.203, -106.832}, {-99.437, 330.139}, {168.889, 0}, {0, 0}, {0, 0}}
        };

        static double kap1921[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.0164609, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-5.33333, -3.72337}, {-2.37037, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{40.786, -22.3402}, {-14.2222, -67.0206}, {-21.3333, 0}, {0, 0}, {0, 0}},
            {{-111.356, 0}, {0, 119.148}, {37.9259, 0}, {0, 0}, {0, 0}}
        };

        static double kap1930[7][5][2] = {
            {{-0.0759415, -0.00295505}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.00480894, 0.00369382}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-1.81002, 0.0871741}, {-0.919459, 0}, {-0.197531, 0}, {0, 0}, {0, 0}},
            {{79.7475, -1.72206}, {57.3171, -1.86168}, {11.2593, 0}, {2.37037, 0}, {0, 0}},
            {{425.579, -76.6479}, {-68.8016, -69.5029}, {129.357, -7.44674}, {-5.53086, 0}, {4.74074, 0}},
            {{-87.8946, 148.481}, {417.612, -311.522}, {-227.16, -253.189}, {-34.7654, 0}, {0, 0}},
            {{-279.268, -135.118}, {-146.853, 652.831}, {331.259, 0}, {0, 0}, {0, 0}}
        };

        static double kap1931[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.0219479, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-8.2963, -5.58505}, {-3.55556, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{70.2698, -49.6449}, {-31.6049, -119.148}, {-37.9259, 0}, {0, 0}, {0, 0}},
            {{-231.893, 18.6168}, {11.8519, 248.225}, {79.0123, 0}, {0, 0}, {0, 0}}
        };

        double m_q_hat = m_q / m_b, z = power_of<2>(m_q_hat);
        double s_hat = s / m_b / m_b;

        complex<double> log_s_hat = { std::log(std::abs(s_hat)), 0.0 };
        if ((0.000 <= s_hat) && (s_hat <= 0.45))
        {
            log_s_hat.imag(0.0);
        }
        else if ((-0.45 <= s_hat) && (s_hat <= -0.000))
        {
            log_s_hat.imag(+M_PI);
        }
        else
        {
            throw InternalError("CharmLoop::F19_massive used outside its domain of validity, s_hat = " + stringify(s_hat));
        }

        const double rho19[4] = {
            3.8991 * power_of<3>(m_q_hat), -23.3946 * m_q_hat, -140.368 * m_q_hat, 7.79821 / m_q_hat - 319.726 * m_q_hat
        };

        // real part
        complex<double> r = (-1424.0 / 729.0 + 64.0 / 27.0 * log(m_q_hat)) * log(mu/m_b)
            - 16.0 / 243.0 * log(mu/m_b) * log_s_hat
            + (16.0 / 1215.0 - 32.0 / 135.0 /power_of<2>(m_q_hat)) * log(mu/m_b) * s_hat
            + (4.0 / 2835.0 - 8.0 / 315.0 /power_of<4>(m_q_hat)) * log(mu/m_b) * s_hat * s_hat
            + (16.0 / 76545.0 - 32.0 /8505.0 / power_of<6>(m_q_hat)) * log(mu/m_b) * power_of<3>(s_hat)
            - 256.0 / 243.0 * power_of<2>(log(mu/m_b));

        for (int l = 3  ; l < 7 ; l++)
            for (int m = 0  ; m < 4  ; m++)
                r = r + kap1900[l][m][0] * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3  ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1901[l][m][0] * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 2  ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap1910[l][m][0] * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4  ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1911[l][m][0] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 1  ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap1920[l][m][0] * s_hat * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3  ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1921[l][m][0] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0  ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap1930[l][m][0] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3  ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap1931[l][m][0] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0  ; l < 4; l++)
            r = r + rho19[l] * pow(s_hat, l);

        // imaginary part
        complex<double> i = 16.0 / 243.0 * M_PI * log(mu/m_b);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1900[l][m][1] * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1901[l][m][1] * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 2 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1910[l][m][1] * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1911[l][m][1] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1920[l][m][1] * s_hat * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1921[l][m][1] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap1930[l][m][1] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap1931[l][m][1] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        return r + complex<double>(0.0, 1.0) * i;
    }

    // cf. [AAGW2001], Eq. (54), p. 19
    complex<double>
    CharmLoops::F29_massive(const double & mu, const double & s, const double & m_b, const double & m_q)
    {
        // F29(s) diverges for s -> 0. However, s * F29(s) -> 0 for s -> 0.
        if (abs(s) < 1e-6) // allow for s = 1e-6, corresponding roughly to the dielectron threshold
            throw InternalError("CharmLoops::F29_massive: F29 diverges for s -> 0. Check that F29 enters via 's * F29(s)' and replace by zero.");

        // cf. [ABGW2001], Appendix B, pp. 34-38
        static double kap2900[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-24.2913, -22.0299}, {-23.1111, -11.1701}, {0, 0}, {0, 0}, {0, 0}},
            {{-86.7723, 97.2931}, {-57.5593, 67.0206}, {7.11111, 44.6804}, {4.74074, 0}, {0, 0}},
            {{96.5187, -160.51}, {-325.463, 89.3609}, {92.4444, 178.722}, {23.7037, 0}, {0, 0}},
            {{88.3801, 142.135}, {171.457, -208.509}, {-120.889, 0}, {0, 0}, {0, 0}}
        };

        static double kap2901[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.296296, 0.620562}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{3.55556, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-29.7586, 11.1701}, {7.11111, 44.6804}, {14.2222, 0}, {0, 0}, {0, 0}},
            {{55.2172, 9.92898}, {6.32099, -59.5739}, {-18.963, 0}, {0, 0}, {0, 0}}
        };

        static double kap2910[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.8462, 1.11701}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-26.8464, 1.86168}, {-8.88889, 11.1701}, {0, 0}, {0, 0}, {0, 0}},
            {{-428.313, 184.792}, {-50.8606, 201.062}, {-75.2337, 44.6804}, {4.74074, 0}, {-4.74074, 0}},
            {{108.781, -396.864}, {-897.575, 402.124}, {295.111, 491.485}, {66.3704, 0}, {0, 0}},
            {{437.34, 382.697}, {408.81, -804.248}, {-381.63, 0}, {0, 0}, {0, 0}}
        };

        static double kap2911[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{16., 11.1701}, {7.11111, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-111.923, 44.6804}, {28.4444, 178.722}, {56.8889, 0}, {0, 0}, {0, 0}},
            {{249.663, 22.3402}, {14.2222, -268.083}, {-85.3333, 0}, {0, 0}, {0, 0}}
        };

        static double kap2920[7][5][2] = {{{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0132191, 0.11968}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.367901, -0.372337}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-222.769, 8.19141}, {-132.372, 11.1701}, {-32., 0}, {-4.74074, 0}, {0, 0}},
            {{-1276.44, 313.249}, {131.529, 312.763}, {-343.034, 44.6804}, {14.2222, 0}, {-14.2222, 0}},
            {{268.098, -652.279}, {-1632.09, 982.969}, {714.667, 938.289}, {128., 0}, {0, 0}},
            {{823.218, 640.989}, {596.622, -1980.83}, {-1013.33, 0}, {0, 0}, {0, 0}}
        };

        static double kap2921[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.0987654, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{32., 22.3402}, {14.2222, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-244.716, 134.041}, {85.3333, 402.124}, {128., 0}, {0, 0}, {0, 0}},
            {{668.137, 0}, {0, -714.887}, {-227.556, 0}, {0, 0}, {0, 0}}
        };

        static double kap2930[7][5][2] = {
            {{-0.0142243, 0.0177303}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0.0288536, -0.0221629}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{10.8601, -0.523045}, {5.51675, 0}, {1.18519, 0}, {0, 0}, {0, 0}},
            {{-478.485, 10.3323}, {-343.902, 11.1701}, {-67.5556, 0}, {-14.2222, 0}, {0, 0}},
            {{-2553.47, 459.887}, {412.809, 417.017}, {-776.143, 44.6804}, {33.1852, 0}, {-28.4444, 0}},
            {{527.368, -890.889}, {-2505.67, 1869.13}, {1362.96, 1519.13}, {208.593, 0}, {0, 0}},
            {{1675.61, 810.709}, {881.117, -3916.98}, {-1987.56, 0}, {0, 0}, {0, 0}}
        };

        static double kap2931[7][5][2] = {
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-0.131687, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{49.7778, 33.5103}, {21.3333, 0}, {0, 0}, {0, 0}, {0, 0}},
            {{-421.619, 297.87}, {189.63, 714.887}, {227.556, 0}, {0, 0}, {0, 0}},
            {{1391.36, -111.701}, {-71.1111, -1489.35}, {-474.074, 0}, {0, 0}, {0, 0}}
        };

        double m_q_hat = m_q / m_b, z = power_of<2>(m_q_hat);
        double s_hat = s / m_b / m_b;

        complex<double> log_s_hat = { std::log(std::abs(s_hat)), 0.0 };
        if ((0.000 <= s_hat) && (s_hat <= 0.45))
        {
            log_s_hat.imag(0.0);
        }
        else if ((-0.45 <= s_hat) && (s_hat <= -0.000))
        {
            log_s_hat.imag(+M_PI);
        }
        else
        {
            throw InternalError("CharmLoop::F29_massive used outside its domain of validity, s_hat = " + stringify(s_hat));
        }

        const double rho29[4] = {
            -23.3946 * power_of<3>(m_q_hat), 140.368 * m_q_hat, 842.206 * m_q_hat, -46.7892 / m_q_hat + 1918.36 * m_q_hat
        };

        // real part
        complex<double> r = (256.0 / 243.0 - 128.0 / 9.0 * log(m_q_hat)) * log(mu / m_b)
            + 32.0 / 81.0 * log(mu / m_b) * log_s_hat
            + (-32.0 / 405.0 + 64.0 / 45 / power_of<2>(m_q_hat)) * log(mu / m_b) * s_hat
            + (-8.0 / 945.0 + 16.0 / 105 / power_of<4>(m_q_hat)) * log(mu / m_b) * s_hat * s_hat
            + (-32.0 / 25515.0 + 64.0 / 2835 / power_of<6>(m_q_hat)) * log(mu / m_b) * power_of<3>(s_hat)
            + 512.0 / 81.0 * power_of<2>(log(mu / m_b));

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 4 ; m++)
                r = r + kap2900[l][m][0] * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2901[l][m][0] * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 2 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap2910[l][m][0] * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2911[l][m][0] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap2920[l][m][0] * s_hat * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2921[l][m][0] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0 ; l < 7 ; l++)
            for (int m = 0 ; m < 5 ; m++)
                r = r + kap2930[l][m][0] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                r = r + kap2931[l][m][0] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0 ; l < 4; l++)
            r = r + rho29[l] * pow(s_hat, l);

        // imaginary part
        complex<double> i = - 32.0 / 81.0 * M_PI * log(mu/m_b);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2900[l][m][1] * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 3 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2901[l][m][1] * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 2 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2910[l][m][1] * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2911[l][m][1] * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 1 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2920[l][m][1] * s_hat * s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2921[l][m][1] * s_hat * s_hat * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 0 ; l < 7 ; l++)
            for (int m = 0 ; m < 3 ; m++)
                i = i + kap2930[l][m][1] * power_of<3>(s_hat) * pow(z, l-3) * pow(log(m_q_hat), m);

        for (int l = 4 ; l < 7 ; l++)
            for (int m = 0 ; m < 2 ; m++)
                i = i + kap2931[l][m][1] * power_of<3>(s_hat) * log_s_hat * pow(z, l-3) * pow(log(m_q_hat), m);

        return r + complex<double>(0.0, 1.0) * i;
    }

    // cf. [AAGW2001], eqs. (48) and (49), p. 18
    complex<double>
    CharmLoops::delta_F29_massive(const double & mu, const double & s, const double & m_q)
    {
        const double x = s / (4.0 * m_q * m_q);

        return 64.0 / 945.0 * (2.0 / 3.0 + log(mu / m_q)) * (105.0 + 84.0 * x + 72.0 * x * x + 64.0 * x * x * x);
    }

    // cf. [BFS2001], Eq. (82), p. 30
    complex<double>
    CharmLoops::F87_massless(const double & mu, const double & s, const double & m_q)
    {
        if (abs(s) < 1e-6) // allow for s = 1e-6, roughly corresponding to the dielectron threshold
            return -4.0 / 9.0 * (complex<double>(8.0 * std::log(mu / m_q) + 11.0, 2.0 * M_PI) + 4.0 * C0(0.0, m_q));

        // Loop-Functions are calculated for the pole mass!
        double s_hat = s / (m_q * m_q);
        double s_hat2 = s_hat * s_hat;
        double denom = (1.0 - s_hat);
        double denom2 = denom * denom;

        complex<double> a = complex<double>(-32.0 * std::log(mu / m_q), -8.0 * M_PI);

        if (std::abs(s_hat - 1.0) < 1e-2)
        {
            // We use a taylor approximation with a maximum error of 1e-9, as the exact expression ist numerical instable.
            static const double c0 = (-67.0 + 6.0 * std::sqrt(3.0) * M_PI) / 27.0;
            static const double c1 = -1.0 + 58.0 * M_PI / (135.0 * std::sqrt(3.0));
            static const double c2 = 4.0 * (-180.0 + 23.0 * std::sqrt(3.0) * M_PI) / 1215.0;
            static const double c3 = -74.0 / 45.0 + 4436.0 * M_PI / (5103.0 * std::sqrt(3.0));

            return a / 9.0 + c0 + c1 * denom + c2 * power_of<2>(denom) + c3 * power_of<3>(denom);
        }

        complex<double> b(-8.0 * s_hat / denom * std::log(s_hat) - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2);
        complex<double> c((4.0 / power_of<3>(denom)) * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_q) - (4.0 + 2.0 * s_hat) * C0(s, m_q)));

        return (1.0 / 9.0) * (a + b + c);
    }

    // cf. [BFS2001], Eq. (83), p. 30
    complex<double>
    CharmLoops::F89_massless(const double & s, const double & m_q)
    {
        // F89(s) diverges for s -> 0. However, s * F89(s) -> 0 for s -> 0.
        if (abs(s) < 1e-6) // allow for s = 1e-6, roughly corresponding to the dielectron threshold
            throw InternalError("CharmLoops::F89_massless: F89 diverges for s -> 0. Check that F89 enters via 's * F89(s)' and replace by zero.");

        // Loop-Functions are calculated for the pole mass!
        double s_hat = s / (m_q * m_q);
        double denom = (1.0 - s_hat);
        double denom2 = denom * denom;

        if (std::abs(s_hat - 1.0) < 1e-2)
        {
            // We use a taylor approximation with a maximum error of 1e-9, as the exact expression ist numerical instable
            static const double c0 = 2.0 * (2.0 * sqrt(3.0) * M_PI - 37.0) / 27.0;
            static const double c1 = 28.0 * M_PI / (45.0 * sqrt(3.0)) - 2.0;
            static const double c2 = 8.0 * (11.0 * sqrt(3.0) * M_PI - 90.0) / 405.0;
            static const double c3 = 4.0 * (17577.0 - 4790.0 * sqrt(3.0) * M_PI) / 76545.0;

            return c0 + c1 * denom + c2 * power_of<2>(denom) + c3 * power_of<3>(denom);;
        }

        double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
        complex<double> b = (-8.0 * (4.0 - s_hat) / power_of<3>(denom)) * ((1.0 + s_hat) * B0(s, m_q) - 2.0 * C0(s, m_q));

        return (1.0 / 9.0) * (a + b);
    }

    // cf. [BFS2001], Eq. (29), p. 8
    complex<double>
    CharmLoops::B0(const double & s, const double & m_q)
    {
        if ((0.0 == m_q) && (0.0 == s))
        {
            throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q == 0 & s == 0");
        }

        if (0.0 == s)
            return complex<double>(-2.0, 0.0);

        double z = 4.0 * m_q * m_q / s;
        complex<double> result;

        if (z > 1.0) // s > 4 m_q^2
        {
            result = complex<double>(-2.0 * std::sqrt(z - 1.0) * std::atan(1.0 / std::sqrt(z - 1.0)), 0.0);
        }
        else if (z > 0) // s > 0 && s <= 4 m_q^2
        {
            result = complex<double>(std::sqrt(1.0 - z) * std::log((1.0 - std::sqrt(1 - z)) / (1.0 + std::sqrt(1.0 - z))),
                    std::sqrt(1.0 - z) * M_PI);
        }
        else if (z < 0) // s < 0
        {
            result = complex<double>(sqrt(1.0 - z) * log((sqrt(1.0 - z) - 1.0) / (sqrt(1.0 - z) + 1.0)), 0.0);
        }

        return result;
    }

    // cf. [BFS2001], Eq. (84), p. 30
    complex<double>
    CharmLoops::C0(const double & s, const double & m_q)
    {
        double s_hat = s / (m_q * m_q);

        if (s_hat < 0.0)
            throw InternalError("CharmLoops::C0: s < 0 is unphysical");

        if (s_hat > 2.0)
            throw InternalError("CharmLoops::C0: support for s > 2.0 * m_q^2 is not implemented, here s / m_q^2 = " + stringify(s_hat));

        if (s_hat < 0.01)
        {
            // the following approximation via linear interpolation yields a difference < 5e-7
            static const double a = -M_PI * M_PI / 6.0;
            static const double b = -0.145134;

            return complex<double>(a + s_hat * b, 0.0);
        }

        if ((0.99 <= s_hat) && (s_hat <= 1.01))
        {
            // the following quadratic approximation yields a difference < 1e-8
            static const double a = -M_PI / std::sqrt(3.0);
            static const double b = (-9.0 + std::sqrt(3.0) * M_PI) / 18.0;
            static const double c = (9.0 - 2.0 * std::sqrt(3.0) * M_PI) / 54.0;

            return complex<double>(a + (s_hat - 1.0) * b + power_of<2>(s_hat - 1.0) * c, 0.0);
        }

        double A = sqrt(s_hat * (4.0 - s_hat));
        double at1 = atan(A / (2.0 - s_hat));
        double at2 = atan(A / s_hat);
        double log1 = log(2.0 - s_hat);

        complex<double> arg;
        gsl_sf_result res_re, res_im;

        arg = 0.5 * complex<double>(2.0 - s_hat, -A);
        gsl_sf_complex_dilog_e(abs(arg), atan2(imag(arg), real(arg)), &res_re, &res_im);
        complex<double> Li_1(res_re.val, res_im.val);

        arg = 0.5 * complex<double>(2.0 - s_hat, +A);
        gsl_sf_complex_dilog_e(abs(arg), atan2(imag(arg), real(arg)), &res_re, &res_im);
        complex<double> Li_2(res_re.val, res_im.val);

        arg = 0.5 * complex<double>(1.0, -A / (2.0 - s_hat));
        gsl_sf_complex_dilog_e(abs(arg), atan2(imag(arg), real(arg)), &res_re, &res_im);
        complex<double> Li_3(res_re.val, res_im.val);

        arg = 0.5 * complex<double>(1.0, +A / (2.0 - s_hat));
        gsl_sf_complex_dilog_e(abs(arg), atan2(imag(arg), real(arg)), &res_re, &res_im);
        complex<double> Li_4(res_re.val, res_im.val);

        return 1.0 / (1.0 - s_hat) * (2.0 * at1 * (at1 - at2) + log1 * log1 - Li_1 - Li_2 + Li_3 + Li_4);
    }

    complex<double>
    CharmLoops::F17_massive_Qsb(const double & s)
    {
        const static thread_local CharmLoopsInterpolation F17_massive_Qsb_Interpolation = {
                {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.7, 2.4, 3.1,
                 3.8, 4.5, 5.2, 5.9, 6.6, 7.3, 8., 8.7, 9.4, 10.1, 10.8, 11.5, 12.2,
                 12.9, 13.6, 14.3, 15.},
                {-0.0717552, -0.071871, -0.0719185, -0.0719426, -0.0719525,
                -0.0719526, -0.0719453, -0.0719321, -0.0719141, -0.0718923,
                -0.071867, -0.0716271, -0.0713227, -0.0709846, -0.0706271,
                -0.0702578, -0.0698812, -0.0695003, -0.0691169, -0.0687322,
                -0.0683473, -0.0679627, -0.0675788, -0.067196, -0.0668144,
                -0.0664342, -0.0660555, -0.0656784, -0.0653028, -0.0649287,
                -0.0645562},
                {-0.00421001, -0.00468545, -0.00504867, -0.00536975, -0.00566363,
                -0.00593744, -0.00619537, -0.0064402, -0.00667389, -0.00689793,
                -0.00711346, -0.00844384, -0.00955976, -0.0105336, -0.0114032,
                -0.0121916, -0.0129146, -0.013583, -0.0142054, -0.0147879,
                -0.0153357, -0.0158529, -0.0163428, -0.0168083, -0.0172516,
                -0.0176749, -0.0180799, -0.018468, -0.0188407, -0.019199, -0.0195441}};

        if (s < 0.0)
        {
            Log::instance()->message("[CharmLoop::F17_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 0 GeV^2.";
            return F17_massive_Qsb_Interpolation(0.);
        }
        else if (s > 15.0)
        {
            Log::instance()->message("[CharmLoop::F17_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 15 GeV^2.";
            return F17_massive_Qsb_Interpolation(15.);
        }
        else return F17_massive_Qsb_Interpolation(s);

    }

    complex<double>
    CharmLoops::F19_massive_Qsb(const double & s)
    {
        const static thread_local CharmLoopsInterpolation F19_massive_Qsb_Interpolation = {
                {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.7, 2.4, 3.1,
                 3.8, 4.5, 5.2, 5.9, 6.6, 7.3, 8., 8.7, 9.4, 10.1, 10.8, 11.5, 12.2,
                 12.9, 13.6, 14.3, 15.},
                {1.51695, 0.420975, 0.367386, 0.33627, 0.314336, 0.297424, 0.283681,
                0.272122, 0.262159, 0.253412, 0.245623, 0.206987, 0.182514, 0.164744,
                0.15088, 0.139568, 0.130052, 0.121864, 0.114698, 0.108341, 0.102639,
                0.0974782, 0.0927718, 0.0884517, 0.0844637, 0.0807641, 0.0773173,
                0.0740936, 0.0710682, 0.0682201, 0.0655313},
                {0.742634, 0.37831, 0.351979, 0.336317, 0.325048, 0.316199, 0.308888,
                0.302645, 0.297185, 0.292328, 0.287947, 0.265299, 0.249968, 0.238219,
                0.228621, 0.220467, 0.213353, 0.207026, 0.201318, 0.19611, 0.191314,
                0.186865, 0.182712, 0.178816, 0.175143, 0.171667, 0.168367, 0.165223,
                0.162222, 0.159348, 0.156591}};

        if (s < 0.0)
        {
            Log::instance()->message("[CharmLoop::F19_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 0 GeV^2.";
            return F19_massive_Qsb_Interpolation(0.);
        }
        else if (s > 15.0)
        {
            Log::instance()->message("[CharmLoop::F19_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 15 GeV^2.";
            return F19_massive_Qsb_Interpolation(15.);
        }
        else return F19_massive_Qsb_Interpolation(s);

    }

    complex<double>
    CharmLoops::F27_massive_Qsb(const double & s)
    {
        const static thread_local CharmLoopsInterpolation F27_massive_Qsb_Interpolation = {
                {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.7, 2.4, 3.1,
                 3.8, 4.5, 5.2, 5.9, 6.6, 7.3, 8., 8.7, 9.4, 10.1, 10.8, 11.5, 12.2,
                 12.9, 13.6, 14.3, 15.},
                {0.430531, 0.431226, 0.431511, 0.431655, 0.431715, 0.431716,
                0.431672, 0.431592, 0.431485, 0.431354, 0.431202, 0.429763, 0.427936,
                0.425908, 0.423762, 0.421547, 0.419287, 0.417002, 0.414701, 0.412393,
                0.410084, 0.407776, 0.405473, 0.403176, 0.400886, 0.398605, 0.396333,
                0.39407, 0.391817, 0.389572, 0.387337},
                {0.0252601, 0.0281127, 0.030292, 0.0322185, 0.0339818, 0.0356246,
                0.0371722, 0.0386412, 0.0400433, 0.0413876, 0.0426807, 0.050663,
                0.0573586, 0.0632015, 0.068419, 0.0731498, 0.0774874, 0.0814983,
                0.0852322, 0.0887275, 0.0920144, 0.0951175, 0.098057, 0.10085,
                0.10351, 0.106049, 0.108479, 0.110808, 0.113044, 0.115194, 0.117265}};

        if (s < 0.0)
        {
            Log::instance()->message("[CharmLoop::F27_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 0 GeV^2.";
            return F27_massive_Qsb_Interpolation(0.);
        }
        else if (s > 15.0)
        {
            Log::instance()->message("[CharmLoop::F27_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 15 GeV^2.";
            return F27_massive_Qsb_Interpolation(15.);
        }
        else return F27_massive_Qsb_Interpolation(s);

    }

    complex<double>
    CharmLoops::F29_massive_Qsb(const double & s)
    {
        const static thread_local CharmLoopsInterpolation F29_massive_Qsb_Interpolation = {
                {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.7, 2.4, 3.1,
                 3.8, 4.5, 5.2, 5.9, 6.6, 7.3, 8., 8.7, 9.4, 10.1, 10.8, 11.5, 12.2,
                 12.9, 13.6, 14.3, 15.},
                {-9.10167, -2.52585, -2.20432, -2.01762, -1.88602, -1.78454,
                -1.70209, -1.63273, -1.57295, -1.52047, -1.47374, -1.24192, -1.09508,
                -0.988462, -0.905279, -0.837411, -0.780312, -0.731185, -0.688188,
                -0.650044, -0.615832, -0.584869, -0.556631, -0.53071, -0.506782,
                -0.484585, -0.463904, -0.444562, -0.426409, -0.40932, -0.393188},
                {-4.4558, -2.26986, -2.11188, -2.0179, -1.95029, -1.89719, -1.85333,
                -1.81587, -1.78311, -1.75397, -1.72768, -1.59179, -1.49981, -1.42931,
                -1.37173, -1.3228, -1.28012, -1.24216, -1.20791, -1.17666, -1.14788,
                -1.12119, -1.09627, -1.07289, -1.05086, -1.03, -1.0102, -0.99134,
                -0.973329, -0.956087, -0.939546}};

        if (s < 0.0)
        {
            Log::instance()->message("[CharmLoop::F29_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 0 GeV^2.";
            return F29_massive_Qsb_Interpolation(0.);
        }
        else if (s > 15.0)
        {
            Log::instance()->message("[CharmLoop::F29_massive_Qsb]", ll_error)
                    << "This function is evaluated outside its domain of validity, at s = " <<  stringify(s) << " GeV^2. Returning the value at s = 15 GeV^2.";
            return F29_massive_Qsb_Interpolation(15.);
        }
        else return F29_massive_Qsb_Interpolation(s);
    }

    complex<double>
    agv_2019a::delta_c7_Qc(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo)
    {
        const agv_2019a::CharmLoopsParameters clp = {/*muhat =*/ mu / m_b, /*s =*/ s, /*z =*/ (m_c * m_c) / (m_b * m_b), /*feynepsilonhat*/ 1e-12};

        // cf. [AGV:2019A] Eq. (2.11), p. 6, and Eq. (2.21), p. 7
        complex<double> result;

        // LO contribution
        result += 0.0;

        if (use_nlo)
        {
            // NLO contribution
            complex<double> nlo = -1.0 * (
                wc.c1() * agv_2019a::F17_Qc(clp)
                + wc.c2() * agv_2019a::F27_Qc(clp));

            result += (alpha_s / (4.0 * M_PI)) * nlo;
        }

        return result;
    }

    complex<double>
    agv_2019a::delta_c7(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo)
    {
        const agv_2019a::CharmLoopsParameters clp = {/*muhat =*/ mu / m_b, /*s =*/ s, /*z =*/ (m_c * m_c) / (m_b * m_b), /*feynepsilonhat*/ 1e-12};

        // cf. [AGV:2019A] Eq. (2.11), p. 6, and Eq. (2.21), p. 7
        complex<double> result;

        // LO contribution
        result += 0.0;

        if (use_nlo)
        {
            // NLO contribution
            complex<double> nlo = -1.0 * (
                wc.c1() * (agv_2019a::F17_Qc(clp) + agv_2019a::F17_Qsb(clp))
                + wc.c2() * (agv_2019a::F27_Qc(clp) + agv_2019a::F27_Qsb(clp)));

            result += (alpha_s / (4.0 * M_PI)) * nlo;
        }

        return result;
    }

    complex<double>
    agv_2019a::delta_c9_Qc(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo)
    {
        const agv_2019a::CharmLoopsParameters clp(mu / m_b, s, (m_c * m_c) / (m_b * m_b), 1e-12);

        // cf. [AGV:2019A] Eq. (2.11), p. 6, and Eq. (2.21), p. 7
        complex<double> result;

        // LO contribution cf. [AGV:2019A] p. 31
        result += wc.c1() * agv_2019a::f190(clp) + wc.c2() * agv_2019a::f290(clp);

        if (use_nlo)
        {
            complex<double> nlo = -1.0 * (
                wc.c1() * agv_2019a::F19_Qc(clp)
                + wc.c2() * agv_2019a::F29_Qc(clp));

            result += (alpha_s / (4.0 * M_PI)) * nlo;
        }

        return result;
    }

    complex<double>
    agv_2019a::delta_c9(const complex<double> & s, const double & mu, const double & alpha_s, const double & m_c, const double & m_b, const WilsonCoefficients<BToS> & wc, bool use_nlo)
    {
        const agv_2019a::CharmLoopsParameters clp(mu / m_b, s, (m_c * m_c) / (m_b * m_b), 1e-12);

        // cf. [AGV:2019A] Eq. (2.11), p. 6, and Eq. (2.21), p. 7
        complex<double> result;

        // LO contribution cf. [AGV:2019A] p. 31
        result += wc.c1() * agv_2019a::f190(clp) + wc.c2() * agv_2019a::f290(clp);

        if (use_nlo)
        {
            complex<double> nlo = -1.0 * (
                wc.c1() * (agv_2019a::F19_Qc(clp) + agv_2019a::F19_Qsb(clp))
                + wc.c2() * (agv_2019a::F29_Qc(clp) + agv_2019a::F29_Qsb(clp)));

            result += (alpha_s / (4.0 * M_PI)) * nlo;
        }

        return result;
    }


    /*
     * Adapter that exports the charn loops function as observables
     */
    template <>
    struct Implementation<CharmLoopsAdapter>
    {
        RestrictedOption opt_contribution;

        UsedParameter m_b;
        UsedParameter m_c;
        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        static const std::map<std::string, std::tuple<double, double, double, double, double, double, double, double, double>> contribution_map;

        double flag_0 = 0.0,
               flag_a = 0.0,
               flag_b = 0.0,
               flag_c = 0.0,
               flag_d = 0.0,
               flag_e = 0.0,
               flag_ctQc = 0.0,
               flag_ctQs = 0.0,
               flag_ctQb = 0.0;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_contribution(o, options, "contribution"_ok),
            m_b(p["mass::b(MSbar)"], u),
            m_c(p["mass::c"], u),
            mu(p["sb::mu"], u)
        {
            auto i = contribution_map.find(opt_contribution.value());
            if (i == contribution_map.end())
                throw InternalError("Unknown charm loops option: " + opt_contribution.value());

            std::tie(flag_0, flag_a, flag_b, flag_c, flag_d, flag_e, flag_ctQc, flag_ctQs, flag_ctQb) = i->second;
        }

        inline agv_2019a::CharmLoopsParameters clp(const complex<double> & s) const
        {
            return agv_2019a::CharmLoopsParameters(mu / m_b, s / m_b() / m_b(), (m_c * m_c) / (m_b * m_b), 1e-12);
        }

        complex<double> F17(const complex<double> & s) const
        {
            auto params = clp(s);
            return flag_a * agv_2019a::f17a(params)
                 + flag_b * agv_2019a::f17b(params)
                 + flag_c * agv_2019a::f17c(params)
                 + flag_d * agv_2019a::f17d(params)
                 + flag_e * agv_2019a::f17e(params)
                 + flag_ctQc * agv_2019a::f17ctQc(params)
                 + flag_ctQs * agv_2019a::f17ctQs(params)
                 + flag_ctQb * agv_2019a::f17ctQb(params);
        }
        complex<double> F19(const complex<double> & s) const
        {
            auto params = clp(s);
            return flag_0 * agv_2019a::f190(params)
                 + flag_a * agv_2019a::f19a(params)
                 + flag_b * agv_2019a::f19b(params)
                 + flag_c * agv_2019a::f19c(params)
                 + flag_d * agv_2019a::f19d(params)
                 + flag_e * agv_2019a::f19e(params)
                 + flag_ctQc * agv_2019a::f19ctQc(params)
                 + flag_ctQs * agv_2019a::f19ctQs(params)
                 + flag_ctQb * agv_2019a::f19ctQb(params);
        }
        complex<double> F27(const complex<double> & s) const
        {
            auto params = clp(s);
            return flag_a * agv_2019a::f27a(params)
                 + flag_b * agv_2019a::f27b(params)
                 + flag_c * agv_2019a::f27c(params)
                 + flag_d * agv_2019a::f27d(params)
                 + flag_e * agv_2019a::f27e(params)
                 + flag_ctQc * agv_2019a::f27ctQc(params)
                 + flag_ctQs * agv_2019a::f27ctQs(params)
                 + flag_ctQb * agv_2019a::f27ctQb(params);
        }
        complex<double> F29(const complex<double> & s) const
        {
            auto params = clp(s);
            return flag_0 * agv_2019a::f290(params)
                 + flag_a * agv_2019a::f29a(params)
                 + flag_b * agv_2019a::f29b(params)
                 + flag_c * agv_2019a::f29c(params)
                 + flag_d * agv_2019a::f29d(params)
                 + flag_e * agv_2019a::f29e(params)
                 + flag_ctQc * agv_2019a::f29ctQc(params)
                 + flag_ctQs * agv_2019a::f29ctQs(params)
                 + flag_ctQb * agv_2019a::f29ctQb(params);
        }
    };

    const std::vector<OptionSpecification>
    Implementation<CharmLoopsAdapter>::options
    {
        { "contribution"_ok, { "0", "Qc", "Qsb", "a", "b", "c", "d", "e", "ctQc", "ctQs", "ctQb", "all" }, "all" },
    };

    const std::map<std::string, std::tuple<double, double, double, double, double, double, double, double, double>>
    Implementation<CharmLoopsAdapter>::Implementation::contribution_map
    {
        std::make_pair("0",     std::make_tuple(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        std::make_pair("Qc",    std::make_tuple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0)),
        std::make_pair("Qsb",   std::make_tuple(0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)),
        std::make_pair("a",     std::make_tuple(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        std::make_pair("b",     std::make_tuple(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        std::make_pair("c",     std::make_tuple(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        std::make_pair("d",     std::make_tuple(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)),
        std::make_pair("e",     std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)),
        std::make_pair("ctQc",  std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)),
        std::make_pair("ctQs",  std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)),
        std::make_pair("ctQsb", std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
        std::make_pair("all",   std::make_tuple(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
    };

    CharmLoopsAdapter::CharmLoopsAdapter(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<CharmLoopsAdapter>(new Implementation<CharmLoopsAdapter>(parameters, options, *this))
    {
    }

    CharmLoopsAdapter::~CharmLoopsAdapter()
    {
    }

    double
    CharmLoopsAdapter::real_F17(const double & re_q2, const double & im_q2) const
    {
        return real(_imp->F17(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::imag_F17(const double & re_q2, const double & im_q2) const
    {
        return imag(_imp->F17(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::real_F27(const double & re_q2, const double & im_q2) const
    {
        return real(_imp->F27(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::imag_F27(const double & re_q2, const double & im_q2) const
    {
        return imag(_imp->F27(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::real_F19(const double & re_q2, const double & im_q2) const
    {
        return real(_imp->F19(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::imag_F19(const double & re_q2, const double & im_q2) const
    {
        return imag(_imp->F19(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::real_F29(const double & re_q2, const double & im_q2) const
    {
        return real(_imp->F29(complex<double>(re_q2, im_q2)));
    }

    double
    CharmLoopsAdapter::imag_F29(const double & re_q2, const double & im_q2) const
    {
        return imag(_imp->F29(complex<double>(re_q2, im_q2)));
    }

    const std::set<ReferenceName>
    CharmLoopsAdapter::references
    {
        "AGV:2019A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    CharmLoopsAdapter::begin_options()
    {
        return Implementation<CharmLoopsAdapter>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    CharmLoopsAdapter::end_options()
    {
        return Implementation<CharmLoopsAdapter>::options.cend();
    }

    namespace agv_2019a
    {
        complex<double> F17_Qc(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f17c(clp) + f17d(clp) + f17e(clp) + f17ctQc(clp);

            return result;
        }

        complex<double> F27_Qc(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f27c(clp) + f27d(clp) + f27e(clp) + f27ctQc(clp);

            return result;
        }

        complex<double> F19_Qc(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f19c(clp) + f19d(clp) + f19e(clp) + f19ctQc(clp);

            return result;
        }

        complex<double> F29_Qc(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f29c(clp) + f29d(clp) + f29e(clp) + f29ctQc(clp);

            return result;
        }

        complex<double> F17_Qsb(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f17a(clp) + f17b(clp) + f17ctQs(clp) + f17ctQb(clp);

            return result;
        }

        complex<double> F27_Qsb(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f27a(clp) + f27b(clp) + f27ctQs(clp) + f27ctQb(clp);

            return result;
        }

        complex<double> F19_Qsb(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f19a(clp) + f19b(clp) + f19ctQs(clp) + f19ctQb(clp);

            return result;
        }

        complex<double> F29_Qsb(const CharmLoopsParameters & clp)
        {
            const complex<double> result = f29a(clp) + f29b(clp) + f29ctQs(clp) + f29ctQb(clp);

            return result;
        }
    }
}
