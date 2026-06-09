/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2016, 2017 Danny van Dyk
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
#include <eos/maths/polylog.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/rare-b-decays/qcdf-integrals-impl.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <limits>

namespace eos
{
    namespace impl
    {
        /* s = 0, cases for B->V gamma */

        // J1
        inline complex<double> j1_szero_bottom(const double & mh, const double & a1, const double & a2)
        {
            static const complex<double> i(0.0, 1.0);
            static const double ln2 = std::log(2.0);
            static const double pi = M_PI, pi2 = M_PI * M_PI, pi3 = pi2 * pi;
            static const double zeta3 = 1.2020569031595942854; // Apery's constant

            double mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh4 * mh2;

            double radix = std::sqrt(4.0 * mh2 - 1.0);
            double atan = std::atan(radix);
            complex<double> zm = complex<double>(0.5, -0.5 * radix);
            complex<double> lnzm = std::log(zm);
            complex<double> dilogzm = dilog(zm);
            complex<double> trilog1 = trilog((radix + i) / (radix - i));

            // Asymptotic result
            complex<double> asymp = 4.0 * (3.0 - 12.0 * mh2 * (1.0 - radix * pi) + 6.0 * mh2 * pi2 - 12.0 * mh4 * pi2 + 2.0 * i * mh2 * pi2 * pi
                    + 6.0 * mh2 * pi2 * std::log(4.0 * mh2) - 12.0 * i * mh2 * pi * lnzm * lnzm + 24.0 * mh2 * zeta3);
            asymp -= 96.0 * mh2 * trilog1;
            asymp += atan * (4.0 * (48.0 * mh4 * pi - 24.0 * mh2 * (radix + pi) + 24.0 * i * mh2 * lnzm * lnzm) + 192.0 * i * mh2 * dilogzm);
            asymp -= 96.0 * mh2 * i * pi * dilogzm;
            asymp += 4.0 * atan * atan * 24.0 * mh2 * (1.0 - 2.0 * mh2 - pi * i - 2.0 * lnzm);
            asymp -= 64.0 * i * mh2 * power_of<3>(atan);
            asymp -= 48.0 * ln2 * pi2 * mh2;

            // First gegenbauer moment
            complex<double> gb1 = 12.0 * (1.0 + 36.0 * mh6 * pi2 - 36.0 * mh4 * (-1.0 + radix * pi + pi2)
                    + mh2 * (-33.0 + 30.0 * radix * pi + 12.0 * pi2 + 2.0 * i * pi3 + 24.0 * zeta3));
            gb1 -= 288.0 * mh2 * trilog1;
            gb1 += atan * (144.0 * mh2 * (-5.0 * radix - 4.0 * pi - 12.0 * mh4 * pi + 6.0 * mh2 * (radix + 2.0 * pi))
                + 288.0 * i * mh2 * lnzm * lnzm + 576.0 * i * mh2 * dilogzm);
            gb1 -= 288.0 * i * mh2 * pi * dilogzm;
            gb1 += atan * atan * (288.0 * mh2 * (2.0 - 6.0 * mh2 + 6.0 * mh4 - i * pi) - 576.0 * mh2 * lnzm);
            gb1 -= 144.0 * i * mh2 * pi * lnzm * lnzm;
            gb1 += 72.0 * mh2 * pi2 * std::log(4.0 * mh2);
            gb1 -= 192.0 * i * mh2 * power_of<3>(atan);
            gb1 -= 144.0 * ln2 * mh2 * pi2;

            // Second gegenbauer moment
            complex<double> gb2 = 8.0 / 3.0 * mh2 * (-533.0 + 438.0 * radix * pi + 144.0 * pi2 + 18.0 * i * pi3 + 216.0 * zeta3);
            gb2 -= 576.0 * mh2 * trilog1;
            gb2 += atan * (-32.0 * mh2 * (73.0 * radix + 48.0 * pi - 600 * mh6 * pi + 60 * mh4 * (5.0 * radix + 9 * pi)
                        - 4 * mh2 * (55.0 * radix + 54.0 * pi))
                + 576.0 * i * mh2 * lnzm * lnzm + 1152.0 * i * mh2 * dilogzm);
            gb2 -= 576.0 * i * mh2 * pi * dilogzm;
            gb2 += atan * atan * (-192.0 * mh2 * (-8.0 + 36.0 * mh2 - 90.0 * mh4 + 100.0 * mh6 + 3.0 * pi * i) - 1152.0 * mh2 * lnzm);
            gb2 -= 288.0 * i * mh2 * pi * lnzm * lnzm;
            gb2 += 144.0 * mh2 * pi2 * std::log(4.0 * mh2);
            gb2 -= 384.0 * i * mh2 * power_of<3>(atan);
            gb2 += 12.0 - 288.0 * ln2 * mh2 * pi2 - 4800.0 * mh4 * mh4 * pi2 + 480.0 * mh6 * (-10.0 + 10.0 * radix * pi + 9.0 * pi2)
                - 16.0 * mh4 * (-245.0 + 220.0 * radix * pi + 108.0 * pi2);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j1_szero_charm(const double & mh, const double & a1, const double & a2)
        {
            static const complex<double> i(0.0, 1.0);
            static const double ln2 = std::log(2.0);
            static const double pi = M_PI, pi2 = M_PI * M_PI;
            static const double zeta3 = 1.2020569031595942854; // Apery's constant

            double mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh4 * mh2, lnmh = std::log(mh);
            double radix = std::sqrt(1.0 - 4.0 * mh2);
            double atanh = std::atanh(radix);
            complex<double> zm = 0.5 * (1.0 - radix);
            complex<double> lnzm = std::log(zm);
            complex<double> dilogzm = dilog(zm);
            complex<double> trilog1 = trilog((-1.0 + radix) / (1.0 + radix));

            complex<double> asymp = -4.0 * (-3.0 + 2.0 * mh2 * (6.0 - 3.0 * (1.0 + 2.0 * lnmh - 2.0 * mh2) * pi2
                    - 12.0 * radix * atanh - 4.0 * atanh * (atanh * (-3.0 + 6.0 * mh2 + 2.0 * atanh) + 6.0 * (dilogzm + atanh * lnzm) + 3.0 * lnzm * lnzm)
                    + 6.0 * i * pi * (1.0 - 4.0 * lnmh * ln2 - 4.0 * mh2 + 2.0 * (atanh + lnzm - ln2 * ln2 + power_of<2>(atanh + lnzm)
                            + 2.0 * ln2 * (2.0 * mh2 + atanh + lnzm)))
                    + 12.0 * trilog1 - 12.0 * zeta3));
            complex<double> gb1 = -12.0 * (-1.0 + mh2 * (33.0 - 6.0 * i * (-5.0 + 4.0 * ln2 * (2.0 * lnmh + ln2)) * pi
                        - 60.0 * radix * atanh + 4.0 * (-3.0 * (1.0 + lnmh) * pi2
                            + 9.0 * mh2 * (-1.0 + pi * (-4.0 * i + 8.0 * i * ln2 + pi) + 2.0 * (radix - 2.0 * atanh) * atanh)
                            + mh4 * (-3.0 * pi * (-8.0 * i + 32.0 * i * ln2 + 3.0 * pi) + 36.0 * atanh * atanh)
                            + 6.0 * i * pi * (atanh + lnzm) * (2.0 + 2.0 * ln2 + atanh + lnzm)
                            - 2.0 * atanh * (2.0 * (-3.0 + atanh) * atanh + 6.0 * (dilogzm + atanh * lnzm) + 3.0 * lnzm * lnzm))
                        + 24.0 * trilog1 - 24.0 * zeta3));
            complex<double> gb2 = 4.0 / 3.0 * (9.0 + 2.0 * mh2 * (-533.0 + 6.0 * i * (-73.0 + 36.0 * ln2 * (2.0 * lnmh + ln2)) * pi +876.0 * radix * atanh
                        + 6.0 * (30.0 * mh4 * (-10.0 + 24.0 * i * (4.0 * ln2 - 1.0) * pi + 9.0 * pi2 + 4.0 * (5.0 * radix - 9.0 * atanh) * atanh)
                            + 20.0 * mh6 * (-32.0 * i * (6.0 * ln2 - 1.0) * pi - 15.0 * pi2 + 60.0 * atanh * atanh)
                            + mh2 * (245.0 - 108.0 * pi * (-4.0 * i + 8.0 * i * ln2 + pi) + 8.0 * atanh * (-55.0 * radix + 54.0 * atanh))
                            + 6.0 * ((4.0 + 3.0 * lnmh) * pi2 - 2.0 * i * pi * (atanh + lnzm) * (8.0 + 6.0 * ln2 + 3.0 * (atanh + lnzm))
                                + 2.0 * atanh * (2.0 * (atanh - 4.0) * atanh + 6.0 * (dilogzm + atanh * lnzm) + 3.0 * lnzm * lnzm)))
                            - 216.0 * trilog1 + 216.0 * zeta3));

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J4
        inline complex<double> j4_szero_bottom(const double & m_b, const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            double mh = m_b / m_B, mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh2 * mh4, mh8 = mh4 * mh4;
            double radix = std::sqrt(4.0 * mh * mh - 1.0);
            double atan = std::atan(1.0 / radix);

            double asymp = 2.0 / 9.0 * (3.0 + 32.0 * mh2 + 48.0 * mh4 - 4 * radix * (1.0 + 2.0 * mh2 + 24.0 * mh4) * atan
                    + 48.0 * mh4 * (4.0 * mh2 - 3.0) * atan * atan - 4.0 * std::log(m_b / mu));
            double gb1 = 1.0 / 3.0 * (1.0 + 24.0 * mh2 + 252.0 * mh4 - 432.0 * mh6 + 432.0 * mh4 * (2.0 * mh2 - 1.0) * radix * atan
                    - 288.0 * (mh4 - 4.0 * mh6 + 6.0 * mh8) * atan * atan);
            double gb2 = 2.0 / 15.0 * (1.0 + 60.0 * mh2 + 2140.0 * mh4 - 9600.0 * mh6 + 14400.0 * mh8
                    - 240.0 * mh4 * radix * (13.0 - 70.0 * mh2 + 120.0 * mh4) * atan
                    + 1440.0 * mh4 * (-1.0 + 8.0 * mh2 - 30.0 * mh4 + 40.0 * mh6) * atan * atan);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j4_szero_charm(const double & m_c, const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            static const complex<double> i(0.0, 1.0);
            static const double pi = M_PI;

            double mh = m_c / m_B, mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh2 * mh4;
            double radix = std::sqrt(1.0 - 4.0 * mh2);
            double ln = std::log((1.0 + radix) / 2.0 / mh);
            complex<double> ln2 = (2.0 * ln - i * pi), ln22 = ln2 * ln2;

            complex<double> asymp = 2.0 / 9.0 * (3.0 - 4.0 * std::log(m_c / mu)
                    + 4.0 * mh2 * (8.0 + 3.0 * mh2 * (4.0 - (-3.0 + 4.0 * mh2) * ln22))
                    + 2.0 * i * radix * (1.0 + 2.0 * mh2 + 24.0 * mh4) * (2.0 * i * ln + pi));
            complex<double> gb1 = 1.0 / 3.0 * (1.0 + 12.0 * mh2 * (2.0 + 3.0 * mh2 * (7.0 - 12.0 * mh2 + 6.0 * radix * (-1.0 + 2.0 * mh2) * ln2
                    - 2.0 * (1.0 - 4.0 * mh2 + 6.0 * mh4) * (-4.0 * ln * ln + 4.0 * i * ln * pi + pi * pi))));
            complex<double> gb2 = 2.0 / 15.0 * (1.0 + 20.0 * mh2 * (3.0 + mh2 * (107.0 - 480.0 * mh2
                            + 6.0 * (i * radix * (13.0 - 70.0 * mh2 + 120.0 * mh4) * (2.0 * i * ln + pi)
                                + 3.0 * ((4.0 * ln * ln - 4.0 * i * ln * pi - pi * pi) * (1.0 - 8.0 * mh2 + 30.0 * mh4 - 40.0 * mh6) + 40.0 * mh4)))));

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j4_szero_massless(const double & mB, const double & mu, const double & a1, const double & a2)
        {
            return 2.0 / 3.0 + 4.0 / 9.0 * complex<double>(-2.0 * std::log(mB / mu), M_PI) + a1 / 3.0 + 2.0 * a2 / 15.0;
        }

        /* J5 */
        inline complex<double> j5_szero_bottom(const double & m_b, const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            double mh = m_b / m_B, mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh2 * mh4;
            double radix = std::sqrt(4.0 * mh * mh - 1.0);
            double atan = std::atan(1.0 / radix);
            double lnmbmu = std::log(m_b / mu);

            double asymp = 2.0 / 9.0 * (13.0 - 156.0 * mh2 + 12.0 * radix * (10.0 * mh2 - 1.0) * atan + 144.0 * mh4 * atan * atan - 12.0 * lnmbmu);
            double gb1 = -2.0 / 3.0 * (-7.0 + 220.0 * mh2 + 96.0 * mh4 - 4.0 * radix * (-1.0 + 34.0 * mh2 + 48.0 * mh4) * atan
                    + 48.0 * mh4 * (-9.0 + 8.0 * mh2) * atan * atan + 4.0 * lnmbmu);
            double gb2 = 1.0 / 3.0 * (17.0 - 1064.0 * mh2 - 1740.0 * mh4 + 2160.0 * mh6
                    - 8.0 * radix * (1.0 - 70.0 * mh2 - 390.0 * mh4 + 540.0 * mh6) * atan
                    + 192.0 * mh4 * (18.0 - 40.0 * mh2 + 45.0 * mh4) * atan * atan - 8.0 * lnmbmu);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j5_szero_charm(const double & m_c, const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            static const complex<double> i(0.0, 1.0);
            static const double pi = M_PI, pi2 = pi * pi;

            double mh = m_c / m_B, mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh2 * mh4;
            double radix = std::sqrt(1.0 - 4.0 * mh * mh);
            double ln = std::log((1.0 + radix) / (2.0 * mh)), lnmqmu = std::log(m_c / mu);
            double atanh = std::atanh(radix);

            complex<double> asymp = 2.0 / 9.0 * (13.0 - 12.0 * lnmqmu - 156.0 * mh2 + 6.0 * radix * (10.0 * mh2 - 1.0) * complex<double>(2.0 * ln, -pi)
                    + 36.0 * mh4 * complex<double>(pi2 - 4.0 * ln * ln, 4.0 * pi * atanh));
            complex<double> gb1 = +2.0 / 3.0 * (7.0 - 4.0 * lnmqmu + 2.0 * radix * (-1.0 + 34.0 * mh2 + 48.0 * mh4) * complex<double>(2.0 * ln, -pi)
                    + 4.0 * mh2 * (-55.0 + 3.0 * mh2 * (4.0 * ln * ln * (8.0 * mh2 - 9.0) - 8.0 + (9.0 - 8.0 * mh2) * pi2))
                    - 48.0 * i * mh4 * (8.0 * mh2 - 9.0) * pi * atanh);
            complex<double> gb2 = 1.0 / 3.0 * (17.0 - 8.0 * lnmqmu - 4.0 * radix * (1.0 - 70.0 * mh2 - 390.0 * mh4 + 540.0 * mh6) * complex<double>(2.0 * ln, -pi)
                    - 4.0 * mh2 * (266.0 + 3.0 * mh2 * (145.0 + 16.0 * ln * ln * (18.0 - 40.0 * mh2 + 45.0 * mh4) - 72.0 * pi2
                            - 20.0 * mh2 * (9.0 + (-8.0 + 9.0 * mh2) * pi2)))
                    + 192.0 * i * mh4 * (18.0 - 40.0 * mh2 + 45.0 * mh4) * pi * atanh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j5_szero_massless(const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            double lnmu = std::log(m_B / mu);

            return 1.0 / 9.0 * (complex<double>(26.0 - 24.0 * lnmu, 12.0 * M_PI)
                    + a1 * complex<double>(42.0 - 24.0 * lnmu, 12.0 * M_PI)
                    + a2 * complex<double>(51.0 - 24.0 * lnmu, 12.0 * M_PI));
        }

        // J6
        inline complex<double> j6_szero_bottom(const double & m_b, const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            double mh = m_b / m_B, mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh2 * mh4;
            double radix = std::sqrt(4.0 * mh * mh - 1.0);
            double atan = std::atan(1.0 / radix);

            complex<double> asymp = -2.0 / 9.0 * (-5.0 + 94.0 * mh2 + 24.0 * mh4 - 4.0 * radix * (-1.0 + 16.0 * mh2 + 12.0 * mh4) * atan
                    + 48.0 * mh4 * (-3.0 + 2.0 * mh2) * atan * atan + 4.0 * std::log(m_b / mu));
            complex<double> gb1 = 1.0 / 3.0 - 92.0 / 3.0 * mh2 - 44.0 * mh4 + 48.0 * mh6
                    + 16.0 * mh2 * atan * (radix * (1.0 + 5.0 * mh2 - 6.0 * mh4) + 6.0 * (mh2 - 2.0 * mh4 + 2.0 * mh6) * atan);
            complex<double> gb2 = 2.0 / 15.0 - 2.0 / 3.0 * mh2 * (55.0 + 236.0 * mh2 - 660.0 * mh4 + 720.0 * mh6)
                    + 16.0 * mh2 * atan * (radix * (1.0 + 16.0 * mh2 - 50.0 * mh4 + 60.0 * mh6)
                            + 12.0 * mh2 * (1.0 - 4.0 * mh2 + 10.0 * mh4 - 10.0 * mh6) * atan);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j6_szero_charm(const double & m_c, const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            static const complex<double> i(0.0, 1.0);
            static const double pi = M_PI;

            double mh = m_c / m_B, mh2 = mh * mh, mh4 = mh2 * mh2, mh6 = mh2 * mh4, mh8 = mh4 * mh4;
            double radix = std::sqrt(1.0 - 4.0 * mh * mh);
            double ln = std::log((1.0 + radix) / 2.0 / mh);
            complex<double> ln2 = (2.0 * ln - i * pi), ln22 = ln2 * ln2;

            complex<double> asymp = -2.0 / 9.0 * (-5.0 + 4.0 * std::log(m_c / mu) + 94.0 * mh2 + 12.0 * mh4 * (2.0 - (-3.0 + 2.0 * mh2) * ln22)
                    + 2.0 * i * radix * (-1.0 + 16.0 * mh2 + 12.0 * mh4) * (2.0 * i * ln + pi));
            complex<double> gb1 = 1.0 / 3.0 * (1.0 + 4.0 * mh2
                    * (-23.0 + 3.0 * mh2 * (-11.0 + 12.0 * mh2 - 6.0 * (1.0 - 2.0 * mh2 + 2.0 * mh4) * ln22)
                        + 6.0 * i * radix * (-1.0 - 5.0 * mh2 + 6.0 * mh4) * (2.0 * i * ln + pi)));
            complex<double> gb2 = 2.0 / 15.0 + 2.0 / 3.0 * mh2 * (-55.0 + 4.0 * mh2 * (-59.0 - 18.0 * ln22) + 12.0 * mh4 * (55.0 + 24.0 * ln22)
                    + 12.0 * radix * (1.0 + 16.0 * mh2 - 50.0 * mh4 + 60.0 * mh6) * ln2 + 720.0 * mh8 * ln22
                    + 720.0 * mh6 * (-1.0 - 4.0 * ln * ln + 4.0 * i * ln * pi + pi * pi));

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j6_szero_massless(const double & m_B, const double & mu, const double & a1, const double & a2)
        {
            return 1.0 / 45.0 * complex<double>(50.0 + 15.0 * a1 + 6.0 * a2 - 40.0 * std::log(m_B / mu), 20.0 * M_PI);
        }

        // J7
        // integration up to a cut-off x ~= Lambda / m_B
        inline double j7_szero(const double & x, const double & a1, const double & a2)
        {
            double lnx = std::log(x);

            return -6.0 * (1.0 - x + lnx) - 6.0 * a1 * (-3.0 * (1.0 - x) * (x - 2.0) + 3.0 * lnx)
                - 6.0 * a2 * (2.0 * (1.0 - x) * (8.0 - 10.0 * x + 5.0 * x * x) + 6.0 * lnx);
        }

        // Massive case: bottom quarks
        struct DileptonIntegralsBottom
        {
            double sh, sh2, sh3, sh4, lnsh;
            double mh, mh2, mh3, mh4, mh6, mh8, mh10, mh12, lnmh;
            double lnmqmu;
            double rho, rho2, rho3, rho4, rho5, rho6, rho7, lnrho, lnrhom1;
            double radixrho, radix4mh2;
            double lnradixrho, lndeltarho4mh2;
            double atanrho, atan4mh2, atanh4mh2rho;
            complex<double> aminus, aplus;
            complex<double> bminus, bplus;
            complex<double> lnam, lnbm, lnradices;
            complex<double> dilogam2;
            complex<double> dilogapbm;
            complex<double> dilogambm;
            complex<double> trilogam2;
            complex<double> trilogapbm;
            complex<double> trilogambm;

            DileptonIntegralsBottom(const double & sh, const double & mh, const double & mB, const double & mu) :
                sh(sh), sh2(sh * sh), sh3(sh2 * sh), sh4(sh2 * sh2),
                lnsh(std::log(sh)),
                mh(mh), mh2(mh * mh), mh3(mh2 * mh), mh4(mh2 * mh2), mh6(mh4 * mh2), mh8(mh4 * mh4), mh10(mh8 * mh2), mh12(mh8 * mh4),
                lnmh(std::log(mh)), lnmqmu(2.0 * std::log(mh * mB / mu)),
                rho(4.0 * mh * mh / sh), rho2(rho * rho), rho3(rho2 * rho), rho4(rho2 * rho2), rho5(rho3 * rho2), rho6(rho3 * rho3), rho7(rho4 * rho3),
                lnrho(std::log(rho)), lnrhom1(std::log(rho - 1.0)),
                radixrho(std::sqrt(rho - 1.0)), radix4mh2(std::sqrt(4.0 * mh2 - 1.0)),
                lnradixrho(0.5 * lnrhom1), lndeltarho4mh2(std::log(rho - 4.0 * mh2)),
                atanrho(std::atan(radixrho)), atan4mh2(std::atan(radix4mh2)), atanh4mh2rho(std::atanh(radix4mh2 / radixrho)),
                aminus(0.5 * complex<double>(1.0, -radixrho)), aplus(1.0 - aminus),
                bminus(0.5 * complex<double>(1.0, -radix4mh2)), bplus(1.0 - bminus),
                lnam(std::log(aminus)), lnbm(std::log(bminus)), lnradices(std::log(radixrho - radix4mh2)),
                dilogam2(dilog(power_of<2>(aminus / aplus))),
                dilogapbm(dilog((aplus * bminus) / (aminus * bplus))),
                dilogambm(dilog((aminus * bminus) / (aplus * bplus))),
                trilogam2(trilog(power_of<2>(aminus / aplus))),
                trilogapbm(trilog((aplus * bminus) / (aminus * bplus))),
                trilogambm(trilog((aminus * bminus) / (aplus * bplus)))
            {
            }

            complex<double> j1(const double & a1, const double & a2) const;
            complex<double> j2(const double & a1, const double & a2) const;
            complex<double> j3(const double & a1, const double & a2) const;
            complex<double> j4(const double & a1, const double & a2) const;
            complex<double> j5(const double & a1, const double & a2) const;
            complex<double> j6(const double & a1, const double & a2) const;
        };

        // J1
        complex<double>
        DileptonIntegralsBottom::j1(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi, pi3 = pi2 * pi;
            static const double ln2 = std::log(2.0);
            static const double zeta3 = 1.2020569031595942854;

            // Asymptotic part
            complex<double> asymp = dilogambm*((complex<double>(0.0,24.0)*atan4mh2*mh2*rho)/( 4.0 *mh2 - rho) - (complex<double>(0.0,12.0)*mh2*pi*rho)/( 4.0 *mh2 - rho)) +
                dilogam2*((complex<double>(0.0,-24.0)*atanrho*mh2*rho)/( 4.0 *mh2 - rho) + (complex<double>(0.0,12.0)*mh2*pi*rho)/( 4.0 *mh2 - rho)) +
                (trilogambm*( 48.0 *mh4*rho -  12.0 *mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) + (trilogapbm*( 48.0 *mh4*rho -  12.0 *mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                (complex<double>(0.0,80.0)*power_of<3>(atan4mh2)*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) -
                (complex<double>(0.0,80.0)*power_of<3>(atanrho)*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                ( 12.0 *trilogam2*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                power_of<2>(atanrho)*((complex<double>(0.0,-48.0)*atan4mh2*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                        ( 48.0 *atanh4mh2rho*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) -
                        ( 24.0 *(- 8.0 *ln2*mh4*rho -  8.0 *lnam*mh4*rho -  8.0 *lnbm*mh4*rho +  8.0 *lnradices*mh4*rho + complex<double>(0.0,12.0)*mh4*pi*rho + mh2*rho2 +
                              2.0 *ln2*mh2*rho2 +  2.0 *lnam*mh2*rho2 +  2.0 *lnbm*mh2*rho2 -  2.0 *lnradices*mh2*rho2 -  2.0 *mh4*rho2 - complex<double>(0.0,3.0)*mh2*pi*rho2))/
                        power_of<2>(- 4.0 *mh2 + rho)) + power_of<2>(atan4mh2)*((- 48.0 *atanh4mh2rho*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                        ( 24.0 *(- 8.0 *ln2*mh4*rho -  8.0 *lnam*mh4*rho -  8.0 *lnbm*mh4*rho +  8.0 *lnradices*mh4*rho + complex<double>(0.0,12.0)*mh4*pi*rho + mh2*rho2 +
                              2.0 *ln2*mh2*rho2 +  2.0 *lnam*mh2*rho2 +  2.0 *lnbm*mh2*rho2 -  2.0 *lnradices*mh2*rho2 -  2.0 *mh4*rho2 - complex<double>(0.0,3.0)*mh2*pi*rho2))/
                        power_of<2>(- 4.0 *mh2 + rho)) + dilogapbm*((complex<double>(0.0,-24.0)*atan4mh2*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                        (complex<double>(0.0,12.0)*(- 4.0 *mh4*pi*rho + mh2*pi*rho2))/power_of<2>(- 4.0 *mh2 + rho)) +
                        atan4mh2*(( 48.0 *atanh4mh2rho*(- 4.0 *mh4*pi*rho + mh2*pi*rho2))/power_of<2>(- 4.0 *mh2 + rho) -
                                ( 24.0 *(- 8.0 *ln2*mh4*pi*rho -  8.0 *lnam*mh4*pi*rho -  8.0 *lnbm*mh4*pi*rho +  8.0 *lnradices*mh4*pi*rho + radix4mh2*mh2*rho2 + mh2*pi*rho2 +
                                      2.0 *ln2*mh2*pi*rho2 +  2.0 *lnam*mh2*pi*rho2 +  2.0 *lnbm*mh2*pi*rho2 -  2.0 *lnradices*mh2*pi*rho2 -  2.0 *mh4*pi*rho2))/
                                power_of<2>(- 4.0 *mh2 + rho)) + (- 8.0 *mh4*(-6.0 + (-6.0 +  6.0 *zeta3 +  6.0 *radixrho*pi - complex<double>(0.0,1.0)*pi3)*rho) +  3.0 *rho2 -
                                 2.0 *mh2*( 12.0 *rho + (6.0 -  6.0 *zeta3 -  6.0 *radix4mh2*pi + complex<double>(0.0,1.0)*pi3)*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                                atanrho*((complex<double>(0.0,48.0)*power_of<2>(atan4mh2)*(- 4.0 *mh4*rho + mh2*rho2))/power_of<2>(- 4.0 *mh2 + rho) -
                                        ( 48.0 *atanh4mh2rho*(- 4.0 *mh4*pi*rho + mh2*pi*rho2))/power_of<2>(- 4.0 *mh2 + rho) +
                                        ( 4.0 *( 24.0 *radixrho*mh4*rho -  48.0 *lnam*mh4*pi*rho -  48.0 *lnbm*mh4*pi*rho +  48.0 *lnradices*mh4*pi*rho - complex<double>(0.0,4.0)*mh4*pi2*rho +
                                             6.0 *mh2*pi*rho2 +  12.0 *lnam*mh2*pi*rho2 +  12.0 *lnbm*mh2*pi*rho2 -  12.0 *lnradices*mh2*pi*rho2 -  12.0 *mh4*pi*rho2 +
                                            complex<double>(0.0,1.0)*mh2*pi2*rho2 +  12.0 *ln2*pi*(- 4.0 *mh4*rho + mh2*rho2)))/power_of<2>(- 4.0 *mh2 + rho));
            // End of asymptotic part

            // 1st Gegenbauer moment
            complex<double> gb1 =  dilogambm*((complex<double>(0.0,72.0)*atan4mh2*mh2*rho)/( 4.0 *mh2 - rho) - (complex<double>(0.0,36.0)*mh2*pi*rho)/( 4.0 *mh2 - rho)) +
                dilogam2*((complex<double>(0.0,-72.0)*atanrho*mh2*rho)/( 4.0 *mh2 - rho) + (complex<double>(0.0,36.0)*mh2*pi*rho)/( 4.0 *mh2 - rho)) -
                (complex<double>(0.0,240.0)*power_of<3>(atan4mh2)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) +
                (complex<double>(0.0,240.0)*power_of<3>(atanrho)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) -
                ( 36.0 *trilogam2*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) +
                ( 36.0 *trilogambm*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) +
                ( 36.0 *trilogapbm*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) +
                dilogapbm*((complex<double>(0.0,72.0)*atan4mh2*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) -
                        (complex<double>(0.0,36.0)*pi*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho)) +
                power_of<2>(atanrho)*((complex<double>(0.0,144.0)*atan4mh2*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) -
                        ( 144.0 *atanh4mh2rho*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) +
                        ( 36.0 *( 64.0 *ln2*mh4*mh6*rho +  64.0 *lnam*mh4*mh6*rho - complex<double>(0.0,96.0)*power_of<2>(mh2)*mh6*pi*rho -  8.0 *power_of<2>(mh2)*mh4*rho2 -
                              32.0 *ln2*power_of<2>(mh4)*rho2 -  32.0 *lnam*power_of<2>(mh4)*rho2 +  16.0 *power_of<2>(mh2)*mh6*rho2 + complex<double>(0.0,48.0)*power_of<2>(mh2)*mh4*pi*rho2 +
                              2.0 *power_of<3>(mh2)*rho3 +  4.0 *ln2*mh2*mh4*rho3 +  4.0 *lnam*mh2*mh4*rho3 -  12.0 *power_of<2>(mh2)*mh4*rho3 +  2.0 *mh6*rho3 - mh2*mh6*rho3 -
                             complex<double>(0.0,1.0)*radix4mh2*mh2*mh6*rho3 +  14.0 *power_of<2>(mh2)*mh6*rho3 + mh8*rho3 + complex<double>(0.0,1.0)*radix4mh2*mh8*rho3 -
                              2.0 *mh2*mh8*rho3 - complex<double>(0.0,6.0)*power_of<3>(mh2)*pi*rho3 +  4.0 *lnbm*power_of<2>(mh2)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3) -
                              4.0 *lnradices*power_of<2>(mh2)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3)))/(mh4*power_of<3>( 4.0 *mh2 - rho))) +
                power_of<2>(atan4mh2)*(( 144.0 *atanh4mh2rho*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) -
                        ( 36.0 *( 64.0 *ln2*mh4*mh6*rho +  64.0 *lnam*mh4*mh6*rho - complex<double>(0.0,96.0)*power_of<2>(mh2)*mh6*pi*rho -  8.0 *power_of<2>(mh2)*mh4*rho2 -
                              32.0 *ln2*power_of<2>(mh4)*rho2 -  32.0 *lnam*power_of<2>(mh4)*rho2 +  16.0 *power_of<2>(mh2)*mh6*rho2 + complex<double>(0.0,48.0)*power_of<2>(mh2)*mh4*pi*rho2 +
                              6.0 *power_of<3>(mh2)*rho3 +  4.0 *ln2*mh2*mh4*rho3 +  4.0 *lnam*mh2*mh4*rho3 -  12.0 *power_of<2>(mh2)*mh4*rho3 -  2.0 *mh6*rho3 + mh2*mh6*rho3 +
                             complex<double>(0.0,1.0)*radix4mh2*mh2*mh6*rho3 +  10.0 *power_of<2>(mh2)*mh6*rho3 - mh8*rho3 - complex<double>(0.0,1.0)*radix4mh2*mh8*rho3 +
                              2.0 *mh2*mh8*rho3 - complex<double>(0.0,6.0)*power_of<3>(mh2)*pi*rho3 +  4.0 *lnbm*power_of<2>(mh2)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3) -
                              4.0 *lnradices*power_of<2>(mh2)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3)))/(mh4*power_of<3>( 4.0 *mh2 - rho))) +
                ( 3.0 *(power_of<2>(complex<double>(0.0,-1.0) + radix4mh2)*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*
                    ( 4.0 *(1.0 - complex<double>(0.0,1.0)*radix4mh2 -  2.0 *mh2 +  12.0 *power_of<2>(lnam)*mh8 -  24.0 *lnam*lnbm*mh8 +  12.0 *power_of<2>(lnbm)*mh8)*rho3 +
                      2.0 *power_of<2>(complex<double>(0.0,1.0) + radix4mh2)*mh4*(- 48.0 *rho +
                          16.0 *(-12.0 +  6.0 *zeta3 +  3.0 *radix4mh2*pi +  9.0 *radixrho*pi - complex<double>(0.0,1.0)*pi3)*rho2 +
                          3.0 *(-9.0 +  2.0 *radix4mh2*(complex<double>(0.0,-2.0) + complex<double>(0.0,3.0)*lnam - complex<double>(0.0,3.0)*lnbm +  6.0 *pi))*rho3) +
                     power_of<2>(complex<double>(0.0,1.0) + radix4mh2)*mh2*( 24.0 *rho2 +
                         (69.0 +  12.0 *power_of<2>(lnam) +  6.0 *lnbm +  12.0 *power_of<2>(lnbm) -  6.0 *lnam*(1.0 +  4.0 *lnbm - complex<double>(0.0,1.0)*radix4mh2) -
                          complex<double>(0.0,3.0)*radix4mh2 - complex<double>(0.0,6.0)*lnbm*radix4mh2 -  24.0 *zeta3 -  60.0 *radix4mh2*pi + complex<double>(0.0,4.0)*pi3)*rho3)) -
                     16.0 *mh6*(complex<double>(0.0,-3.0)*(complex<double>(0.0,-2.0) - radix4mh2 +  2.0 *radixrho - complex<double>(0.0,2.0)*radix4mh2*radixrho +
                            radix4mh2*power_of<2>(radixrho) +  2.0 *lnbm*
                            ( 2.0 *radixrho + radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) + complex<double>(0.0,1.0)*(-2.0 + rho)) +
                            power_of<2>(lnam)*( 8.0 *radixrho + complex<double>(0.0,4.0)*(-2.0 + rho)) + power_of<2>(lnbm)*( 8.0 *radixrho + complex<double>(0.0,4.0)*(-2.0 + rho)) -
                             2.0 *lnam*(complex<double>(0.0,-2.0) +  2.0 *radixrho + radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) +
                                lnbm*( 8.0 *radixrho + complex<double>(0.0,4.0)*(-2.0 + rho)) + complex<double>(0.0,1.0)*rho) + complex<double>(0.0,1.0)*rho)*rho3 -
                        complex<double>(0.0,6.0)*mh2*(complex<double>(0.0,-6.0) -  4.0 *radix4mh2 +  6.0 *radixrho - complex<double>(0.0,8.0)*radix4mh2*radixrho +
                             4.0 *radix4mh2*power_of<2>(radixrho) +  6.0 *lnbm*radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) +
                            lnam*(- 6.0 *radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) +
                                lnbm*( 4.0 *radixrho -  2.0 *radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) + complex<double>(0.0,2.0)*(-2.0 + rho))) +
                            power_of<2>(lnam)*(- 2.0 *radixrho + radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) - complex<double>(0.0,1.0)*(-2.0 + rho)) +
                            power_of<2>(lnbm)*(- 2.0 *radixrho + radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho) - complex<double>(0.0,1.0)*(-2.0 + rho)) + complex<double>(0.0,3.0)*rho)
                        *rho3 +  4.0 *power_of<2>(mh2)*( 4.0 *(-33.0 +  24.0 *zeta3 - complex<double>(0.0,72.0)*pi - complex<double>(0.0,4.0)*pi3)*rho2 +
                             12.0 *power_of<3>(radixrho)*(complex<double>(0.0,-1.0)*rho +  12.0 *pi*rho - complex<double>(0.0,2.0)*rho2 +  6.0 *pi*rho2) +
                            rho*(208.0 -  192.0 *zeta3 + complex<double>(0.0,32.0)*pi3 - complex<double>(0.0,144.0)*pi*(-2.0 + rho2) -  102.0 *rho2 +  3.0 *power_of<2>(lnam)*rho3 -
                                 6.0 *lnam*lnbm*rho3 +  3.0 *power_of<2>(lnbm)*rho3) -
                            complex<double>(0.0,2.0)*radixrho*(-32.0 +  2.0 *(-57.0 +  48.0 *zeta3 - complex<double>(0.0,36.0)*pi - complex<double>(0.0,8.0)*pi3)*rho +
                                (-66.0 - complex<double>(0.0,36.0)*pi)*rho2 -  12.0 *rho3 +  3.0 *power_of<2>(lnam)*rho3 -  6.0 *lnam*lnbm*rho3 +  3.0 *power_of<2>(lnbm)*rho3) +
                             2.0 *(32.0 + (78.0 + complex<double>(0.0,72.0)*pi)*rho2 -  3.0 *(-5.0 + power_of<2>(lnam) -  2.0 *lnam*lnbm + power_of<2>(lnbm))*rho3)))))/
                            (32.*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*mh4*power_of<3>( 4.0 *mh2 - rho)) +
                            atan4mh2*((- 144.0 *atanh4mh2rho*pi*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) -
                                    ( 9.0 *(- 2048.0 *ln2*mh4*mh6*pi*rho -  2048.0 *lnam*mh4*mh6*pi*rho -  16.0 *radix4mh2*mh4*rho2 +  16.0 *power_of<5>(radix4mh2)*mh4*rho2 +
                                         128.0 *radix4mh2*mh2*mh4*rho2 +  64.0 *lnam*mh4*pi*rho2 +  128.0 *lnam*power_of<2>(radix4mh2)*mh4*pi*rho2 +
                                         64.0 *lnam*power_of<4>(radix4mh2)*mh4*pi*rho2 +  256.0 *power_of<2>(mh2)*mh4*pi*rho2 +  1024.0 *ln2*power_of<2>(mh4)*pi*rho2 -
                                         512.0 *power_of<2>(mh2)*mh6*pi*rho2 + complex<double>(0.0,4.0)*lnam*mh2*rho3 +  11.0 *radix4mh2*mh2*rho3 +
                                        complex<double>(0.0,8.0)*lnam*power_of<2>(radix4mh2)*mh2*rho3 + complex<double>(0.0,4.0)*lnam*power_of<4>(radix4mh2)*mh2*rho3 -
                                         11.0 *power_of<5>(radix4mh2)*mh2*rho3 -  88.0 *radix4mh2*power_of<2>(mh2)*rho3 - complex<double>(0.0,16.0)*power_of<3>(mh2)*rho3 -  6.0 *radix4mh2*mh4*rho3 +
                                         6.0 *power_of<5>(radix4mh2)*mh4*rho3 +  48.0 *radix4mh2*mh2*mh4*rho3 + complex<double>(0.0,16.0)*mh6*rho3 - complex<double>(0.0,64.0)*lnam*mh6*rho3 +
                                         16.0 *radix4mh2*mh6*rho3 + complex<double>(0.0,32.0)*lnam*mh2*mh6*rho3 +  96.0 *radix4mh2*mh2*mh6*rho3 -  32.0 *lnam*radix4mh2*mh2*mh6*rho3 -
                                        complex<double>(0.0,64.0)*lnam*power_of<2>(mh2)*mh6*rho3 - complex<double>(0.0,16.0)*lnam*mh8*rho3 +  32.0 *lnam*radix4mh2*mh8*rho3 +
                                        complex<double>(0.0,16.0)*lnam*power_of<2>(radix4mh2)*mh8*rho3 -  8.0 *lnam*mh2*pi*rho3 -  16.0 *lnam*power_of<2>(radix4mh2)*mh2*pi*rho3 -
                                         8.0 *lnam*power_of<4>(radix4mh2)*mh2*pi*rho3 -  128.0 *power_of<3>(mh2)*pi*rho3 -  128.0 *ln2*mh2*mh4*pi*rho3 +  384.0 *power_of<2>(mh2)*mh4*pi*rho3 -
                                         384.0 *power_of<2>(mh2)*mh6*pi*rho3 +  128.0 *lnradices*power_of<2>(mh2)*pi*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3) -
                                         32.0 *lnbm*(mh2*(-((complex<double>(0.0,-1.0) + radix4mh2)*mh6) + complex<double>(0.0,2.0)*mh8)*rho3 +
                                            (complex<double>(0.0,-2.0)*mh6 + (complex<double>(0.0,-1.0) + radix4mh2)*mh8)*rho3 +  2.0 *power_of<3>(mh2)*(complex<double>(0.0,1.0) +  2.0 *pi)*rho3 +
                                            power_of<2>(mh2)*( 64.0 *mh6*pi*rho -  32.0 *mh4*pi*rho2 - complex<double>(0.0,2.0)*mh6*rho3))))/(8.*mh4*power_of<3>( 4.0 *mh2 - rho))) +
                            atanrho*((- 36.0 *atan4mh2*( 4.0 *mh6 +  4.0 *power_of<2>(mh2)*mh6 -  2.0 *mh2*( 2.0 *mh4 + mh6 + complex<double>(0.0,1.0)*radix4mh2*mh6) -
                                            power_of<2>(complex<double>(0.0,-1.0) + radix4mh2)*mh8)*rho3)/(mh4*power_of<3>( 4.0 *mh2 - rho)) -
                                    (complex<double>(0.0,144.0)*power_of<2>(atan4mh2)*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) +
                                    ( 144.0 *atanh4mh2rho*pi*( 16.0 *mh6*rho -  8.0 *mh4*rho2 + mh2*rho3))/power_of<3>( 4.0 *mh2 - rho) -
                                    ( 3.0 *( 6144.0 *lnbm*power_of<2>(mh2)*mh6*pi*rho -  6144.0 *lnradices*power_of<2>(mh2)*mh6*pi*rho +  6144.0 *ln2*mh4*mh6*pi*rho +
                                        complex<double>(0.0,512.0)*power_of<2>(mh2)*mh6*pi2*rho -  768.0 *power_of<2>(mh2)*mh4*pi*rho2 -  3072.0 *lnbm*power_of<2>(mh2)*mh4*pi*rho2 +
                                         3072.0 *lnradices*power_of<2>(mh2)*mh4*pi*rho2 -  3072.0 *ln2*power_of<2>(mh4)*pi*rho2 +  1536.0 *power_of<2>(mh2)*mh6*pi*rho2 -
                                        complex<double>(0.0,256.0)*power_of<2>(mh2)*mh4*pi2*rho2 +  2304.0 *radixrho*mh4*(mh4*rho2 - mh6*( 2.0 *rho + rho2)) -  3.0 *radix4mh2*mh2*rho3 +
                                         3.0 *power_of<5>(radix4mh2)*mh2*rho3 +  24.0 *radix4mh2*power_of<2>(mh2)*rho3 + complex<double>(0.0,48.0)*power_of<3>(mh2)*rho3 +
                                        complex<double>(0.0,192.0)*lnbm*power_of<3>(mh2)*rho3 -  18.0 *radix4mh2*mh4*rho3 +  18.0 *power_of<5>(radix4mh2)*mh4*rho3 +
                                         144.0 *radix4mh2*mh2*mh4*rho3 - complex<double>(0.0,48.0)*mh6*rho3 - complex<double>(0.0,192.0)*lnbm*mh6*rho3 -  48.0 *radix4mh2*mh6*rho3 +
                                        complex<double>(0.0,96.0)*lnbm*mh2*mh6*rho3 -  288.0 *radix4mh2*mh2*mh6*rho3 -  96.0 *lnbm*radix4mh2*mh2*mh6*rho3 -
                                        complex<double>(0.0,192.0)*lnbm*power_of<2>(mh2)*mh6*rho3 - complex<double>(0.0,96.0)*lnbm*mh8*rho3 +  96.0 *lnbm*radix4mh2*mh8*rho3 +
                                        complex<double>(0.0,192.0)*lnbm*mh2*mh8*rho3 +  384.0 *power_of<3>(mh2)*pi*rho3 +  384.0 *lnbm*power_of<3>(mh2)*pi*rho3 -
                                         384.0 *lnradices*power_of<3>(mh2)*pi*rho3 +  384.0 *ln2*mh2*mh4*pi*rho3 -  1152.0 *power_of<2>(mh2)*mh4*pi*rho3 +
                                         1152.0 *power_of<2>(mh2)*mh6*pi*rho3 + complex<double>(0.0,32.0)*power_of<3>(mh2)*pi2*rho3 +
                                         96.0 *lnam*(mh2*((complex<double>(0.0,-1.0) + radix4mh2)*mh6 - complex<double>(0.0,2.0)*mh8)*rho3 +
                                            complex<double>(0.0,1.0)*( 2.0 *mh6 + mh8 + complex<double>(0.0,1.0)*radix4mh2*mh8)*rho3 +  2.0 *power_of<3>(mh2)*(complex<double>(0.0,-1.0) +  2.0 *pi)*rho3 +
                                            power_of<2>(mh2)*( 64.0 *mh6*pi*rho -  32.0 *mh4*pi*rho2 + complex<double>(0.0,2.0)*mh6*rho3))))/(8.*mh4*power_of<3>( 4.0 *mh2 - rho)));
            // End of 1st Gegenbauer moment

            // 2nd Gegenbauer moment
            complex<double> gb2 = dilogambm*((complex<double>(0.0,144.0)*atan4mh2*mh2*rho)/( 4.0 *mh2 - rho) - (complex<double>(0.0,72.0)*mh2*pi*rho)/( 4.0 *mh2 - rho)) +
                dilogam2*((complex<double>(0.0,-144.0)*atanrho*mh2*rho)/( 4.0 *mh2 - rho) + (complex<double>(0.0,72.0)*mh2*pi*rho)/( 4.0 *mh2 - rho)) +
                (complex<double>(0.0,480.0)*power_of<3>(atan4mh2)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) -
                (complex<double>(0.0,480.0)*power_of<3>(atanrho)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                ( 72.0 *trilogam2*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) -
                ( 72.0 *trilogambm*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) -
                ( 72.0 *trilogapbm*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                dilogapbm*((complex<double>(0.0,-144.0)*atan4mh2*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                        (complex<double>(0.0,72.0)*pi*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho)) +
                power_of<2>(atan4mh2)*((- 288.0 *atanh4mh2rho*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                        ( 24.0 *(- 768.0 *ln2*mh4*mh8*rho -  768.0 *lnam*mh4*mh8*rho + complex<double>(0.0,1152.0)*power_of<2>(mh2)*mh8*pi*rho +  46.0 *power_of<3>(mh2)*rho4 +
                              12.0 *ln2*mh2*mh4*rho4 +  12.0 *lnam*mh2*mh4*rho4 -  72.0 *power_of<2>(mh2)*mh4*rho4 -  30.0 *mh6*rho4 +
                              15.0 *mh2*mh6*rho4 + complex<double>(0.0,15.0)*radix4mh2*mh2*mh6*rho4 +  150.0 *power_of<2>(mh2)*mh6*rho4 -
                              15.0 *mh8*rho4 - complex<double>(0.0,15.0)*radix4mh2*mh8*rho4 +  30.0 *mh2*mh8*rho4 -  200.0 *power_of<2>(mh2)*mh8*rho4 -
                             complex<double>(0.0,18.0)*power_of<3>(mh2)*pi*rho4 +  96.0 *power_of<2>(mh2)*mh6*rho2 +  576.0 *ln2*mh4*mh6*rho2 +  576.0 *lnam*mh4*mh6*rho2 -
                              192.0 *power_of<2>(mh2)*mh8*rho2 - complex<double>(0.0,864.0)*power_of<2>(mh2)*mh6*pi*rho2 -  48.0 *power_of<2>(mh2)*mh4*rho3 -  144.0 *ln2*power_of<2>(mh4)*rho3 -
                              144.0 *lnam*power_of<2>(mh4)*rho3 +  96.0 *power_of<2>(mh2)*mh6*rho3 + complex<double>(0.0,216.0)*power_of<2>(mh2)*mh4*pi*rho3 +
                              12.0 *lnbm*power_of<2>(mh2)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3) -
                              12.0 *lnradices*power_of<2>(mh2)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3)))/(mh4*power_of<4>(- 4.0 *mh2 + rho)))\
                + power_of<2>(atanrho)*((complex<double>(0.0,-288.0)*atan4mh2*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                        ( 288.0 *atanh4mh2rho*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                        ( 24.0 *( 768.0 *ln2*mh4*mh8*rho +  768.0 *lnam*mh4*mh8*rho - complex<double>(0.0,1152.0)*power_of<2>(mh2)*mh8*pi*rho +  14.0 *power_of<3>(mh2)*rho4 -
                              12.0 *ln2*mh2*mh4*rho4 -  12.0 *lnam*mh2*mh4*rho4 +  72.0 *power_of<2>(mh2)*mh4*rho4 -  30.0 *mh6*rho4 +
                              15.0 *mh2*mh6*rho4 + complex<double>(0.0,15.0)*radix4mh2*mh2*mh6*rho4 -  210.0 *power_of<2>(mh2)*mh6*rho4 -
                              15.0 *mh8*rho4 - complex<double>(0.0,15.0)*radix4mh2*mh8*rho4 +  30.0 *mh2*mh8*rho4 +  200.0 *power_of<2>(mh2)*mh8*rho4 +
                             complex<double>(0.0,18.0)*power_of<3>(mh2)*pi*rho4 -  96.0 *power_of<2>(mh2)*mh6*rho2 -  576.0 *ln2*mh4*mh6*rho2 -  576.0 *lnam*mh4*mh6*rho2 +
                              192.0 *power_of<2>(mh2)*mh8*rho2 + complex<double>(0.0,864.0)*power_of<2>(mh2)*mh6*pi*rho2 +  48.0 *power_of<2>(mh2)*mh4*rho3 +  144.0 *ln2*power_of<2>(mh4)*rho3 +
                              144.0 *lnam*power_of<2>(mh4)*rho3 -  96.0 *power_of<2>(mh2)*mh6*rho3 - complex<double>(0.0,216.0)*power_of<2>(mh2)*mh4*pi*rho3 -
                              12.0 *lnbm*power_of<2>(mh2)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3) +
                              12.0 *lnradices*power_of<2>(mh2)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3)))/(mh4*power_of<4>(- 4.0 *mh2 + rho)))\
                + atanrho*(( 360.0 *atan4mh2*( 4.0 *mh6 +  4.0 *power_of<2>(mh2)*mh6 -  2.0 *mh2*( 2.0 *mh4 + mh6 + complex<double>(0.0,1.0)*radix4mh2*mh6) -
                                power_of<2>(complex<double>(0.0,-1.0) + radix4mh2)*mh8)*rho4)/(mh4*power_of<4>(- 4.0 *mh2 + rho)) +
                        (complex<double>(0.0,288.0)*power_of<2>(atan4mh2)*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                        ( 288.0 *atanh4mh2rho*pi*( 64.0 *mh8*rho - mh2*rho4 -  48.0 *mh6*rho2 +  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                        (- 1179648.0 *lnbm*power_of<2>(mh2)*mh8*pi*rho +  1179648.0 *lnradices*power_of<2>(mh2)*mh8*pi*rho -  1179648.0 *ln2*mh4*mh8*pi*rho -
                         complex<double>(0.0,98304.0)*power_of<2>(mh2)*mh8*pi2*rho - complex<double>(0.0,200.0)*mh2*rho4 + complex<double>(0.0,705.0)*lnbm*mh2*rho4 -
                          720.0 *radix4mh2*mh2*rho4 +  570.0 *lnbm*radix4mh2*mh2*rho4 +  720.0 *power_of<5>(radix4mh2)*mh2*rho4 -
                          1170.0 *lnbm*power_of<5>(radix4mh2)*mh2*rho4 - complex<double>(0.0,200.0)*power_of<6>(radix4mh2)*mh2*rho4 +
                         complex<double>(0.0,930.0)*lnbm*power_of<6>(radix4mh2)*mh2*rho4 -  810.0 *lnbm*power_of<7>(radix4mh2)*mh2*rho4 +
                         complex<double>(0.0,225.0)*lnbm*power_of<8>(radix4mh2)*mh2*rho4 -  210.0 *lnbm*power_of<9>(radix4mh2)*mh2*rho4 +
                         complex<double>(0.0,2400.0)*power_of<2>(mh2)*rho4 - complex<double>(0.0,7560.0)*lnbm*power_of<2>(mh2)*rho4 +
                          5760.0 *radix4mh2*power_of<2>(mh2)*rho4 -  3000.0 *lnbm*radix4mh2*power_of<2>(mh2)*rho4 +
                         complex<double>(0.0,1920.0)*power_of<3>(mh2)*rho4 + complex<double>(0.0,69120.0)*lnbm*power_of<3>(mh2)*rho4 - complex<double>(0.0,1200.0)*mh4*rho4 +
                         complex<double>(0.0,1260.0)*lnbm*mh4*rho4 -  4320.0 *radix4mh2*mh4*rho4 -  2280.0 *lnbm*radix4mh2*mh4*rho4 +
                          4320.0 *power_of<5>(radix4mh2)*mh4*rho4 +  3480.0 *lnbm*power_of<5>(radix4mh2)*mh4*rho4 -
                         complex<double>(0.0,1200.0)*power_of<6>(radix4mh2)*mh4*rho4 + complex<double>(0.0,1860.0)*lnbm*power_of<6>(radix4mh2)*mh4*rho4 +
                          1200.0 *lnbm*power_of<7>(radix4mh2)*mh4*rho4 + complex<double>(0.0,600.0)*lnbm*power_of<8>(radix4mh2)*mh4*rho4 +
                         complex<double>(0.0,14400.0)*mh2*mh4*rho4 - complex<double>(0.0,12720.0)*lnbm*mh2*mh4*rho4 +  34560.0 *radix4mh2*mh2*mh4*rho4 +
                          13440.0 *lnbm*radix4mh2*mh2*mh4*rho4 - complex<double>(0.0,44800.0)*power_of<2>(mh2)*mh4*rho4 +
                         complex<double>(0.0,29760.0)*lnbm*power_of<2>(mh2)*mh4*rho4 - complex<double>(0.0,10440.0)*mh6*rho4 -
                         complex<double>(0.0,46080.0)*lnbm*mh6*rho4 -  11520.0 *radix4mh2*mh6*rho4 +  1440.0 *lnbm*radix4mh2*mh6*rho4 -
                          1440.0 *lnbm*power_of<5>(radix4mh2)*mh6*rho4 + complex<double>(0.0,1080.0)*power_of<6>(radix4mh2)*mh6*rho4 -
                         complex<double>(0.0,12960.0)*mh2*mh6*rho4 + complex<double>(0.0,23040.0)*lnbm*mh2*mh6*rho4 -  69120.0 *radix4mh2*mh2*mh6*rho4 -
                          34560.0 *lnbm*radix4mh2*mh2*mh6*rho4 + complex<double>(0.0,128640.0)*power_of<2>(mh2)*mh6*rho4 -
                         complex<double>(0.0,69120.0)*lnbm*power_of<2>(mh2)*mh6*rho4 - complex<double>(0.0,23040.0)*lnbm*mh8*rho4 +
                          23040.0 *lnbm*radix4mh2*mh8*rho4 + complex<double>(0.0,46080.0)*lnbm*mh2*mh8*rho4 -
                         complex<double>(0.0,69120.0)*power_of<2>(mh2)*mh8*rho4 - complex<double>(0.0,153600.0)*lnbm*power_of<2>(mh2)*mh8*rho4 +
                          24576.0 *power_of<3>(mh2)*pi*rho4 +  18432.0 *lnbm*power_of<3>(mh2)*pi*rho4 -  18432.0 *lnradices*power_of<3>(mh2)*pi*rho4 +
                          18432.0 *ln2*mh2*mh4*pi*rho4 -  110592.0 *power_of<2>(mh2)*mh4*pi*rho4 +  276480.0 *power_of<2>(mh2)*mh6*pi*rho4 -
                          307200.0 *power_of<2>(mh2)*mh8*pi*rho4 + complex<double>(0.0,1536.0)*power_of<3>(mh2)*pi2*rho4 +  147456.0 *power_of<2>(mh2)*mh6*pi*rho2 +
                          884736.0 *lnbm*power_of<2>(mh2)*mh6*pi*rho2 -  884736.0 *lnradices*power_of<2>(mh2)*mh6*pi*rho2 +  884736.0 *ln2*mh4*mh6*pi*rho2 -
                          294912.0 *power_of<2>(mh2)*mh8*pi*rho2 + complex<double>(0.0,73728.0)*power_of<2>(mh2)*mh6*pi2*rho2 -  73728.0 *power_of<2>(mh2)*mh4*pi*rho3 -
                          221184.0 *lnbm*power_of<2>(mh2)*mh4*pi*rho3 +  221184.0 *lnradices*power_of<2>(mh2)*mh4*pi*rho3 -  221184.0 *ln2*power_of<2>(mh4)*pi*rho3 +
                          147456.0 *power_of<2>(mh2)*mh6*pi*rho3 - complex<double>(0.0,18432.0)*power_of<2>(mh2)*mh4*pi2*rho3 +
                          4096.0 *radixrho*mh4*( 2.0 *mh8*( 112.0 *rho +  50.0 *rho2 +  75.0 *rho3) +  27.0 *(- 6.0 *mh6*rho2 +  2.0 *mh4*rho3 -  5.0 *mh6*rho3)) +
                          4608.0 *lnam*( 5.0 *mh2*((complex<double>(0.0,-1.0) + radix4mh2)*mh6 - complex<double>(0.0,2.0)*mh8)*rho4 +
                                 complex<double>(0.0,5.0)*( 2.0 *mh6 + mh8 + complex<double>(0.0,1.0)*radix4mh2*mh8)*rho4 +  2.0 *power_of<3>(mh2)*(complex<double>(0.0,-5.0) +  2.0 *pi)*rho4 -
                                  2.0 *power_of<2>(mh2)*( 128.0 *mh8*pi*rho - complex<double>(0.0,5.0)*mh6*rho4 -  96.0 *mh6*pi*rho2 +  24.0 *mh4*pi*rho3)))/
                                 (64.*mh4*power_of<4>(- 4.0 *mh2 + rho))) + ( 9.0 *mh6*
                                     (- 4096.0 *power_of<2>(complex<double>(0.0,-1.0) + radixrho)*mh4*rho +
                                       5.0 *rho4*( 96.0 *power_of<2>(lnam)*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*(2.0 + (-1.0 - complex<double>(0.0,1.0)*radix4mh2)*mh2 +  2.0 *power_of<2>(mh2)) +
                                           12.0 *power_of<2>(lnbm)*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*
                                          (16.0 + complex<double>(0.0,1.0)*radix4mh2 - complex<double>(0.0,1.0)*power_of<5>(radix4mh2) -  8.0 *mh2 - complex<double>(0.0,16.0)*radix4mh2*mh2 +  32.0 *power_of<2>(mh2)) -
                                          lnbm*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*(-87.0 +  9.0 *power_of<6>(radix4mh2) -  108.0 *mh2 +  1072.0 *power_of<2>(mh2) +
                                              complex<double>(0.0,96.0)*radix4mh2*(1.0 +  6.0 *mh2)) + lnam*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*
                                          (-87.0 +  9.0 *power_of<6>(radix4mh2) -  108.0 *mh2 +  1072.0 *power_of<2>(mh2) + complex<double>(0.0,96.0)*radix4mh2*(1.0 +  6.0 *mh2) +
                                            12.0 *lnbm*(complex<double>(0.0,1.0)*power_of<5>(radix4mh2) + complex<double>(0.0,1.0)*radix4mh2*(-1.0 +  24.0 *mh2) +  16.0 *(-2.0 + mh2 -  3.0 *power_of<2>(mh2)))) +
                                          power_of<2>(complex<double>(0.0,1.0) + radix4mh2)*(-61.0 + complex<double>(0.0,22.0)*radixrho + power_of<2>(complex<double>(0.0,-1.0) + radixrho)*power_of<2>(1.0 -  4.0 *mh2) +
                                               2.0 *radix4mh2*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*(-1.0 +  4.0 *mh2)*(complex<double>(0.0,-11.0) +  40.0 *pi) -  59.0 *(-1.0 + rho) -
                                               2.0 *(-1.0 +  4.0 *mh2)*(-162.0 +  10.0 *radixrho*(complex<double>(0.0,-15.0) +  16.0 *pi) + complex<double>(0.0,80.0)*pi*(-2.0 + rho) +  51.0 *rho) +
                                               2.0 *radix4mh2*( 2.0 *radixrho*(71.0 + complex<double>(0.0,40.0)*pi) -  40.0 *pi*(-2.0 + rho) + complex<double>(0.0,1.0)*(-166.0 +  23.0 *rho)))) +
                                       384.0 *power_of<2>(complex<double>(0.0,-1.0) + radixrho)*mh4*( 2.0 *(-29.0 +  24.0 *zeta3 +  8.0 *radix4mh2*pi + radixrho*(complex<double>(0.0,-5.0) +  36.0 *pi) -
                                              complex<double>(0.0,4.0)*pi3)*rho2 +  5.0 *(-9.0 +  4.0 *radixrho*(complex<double>(0.0,-1.0) +  3.0 *pi))*rho3)) +
                                     power_of<2>(complex<double>(0.0,-1.0) + radixrho)*mh4*(( 8.0 *(18.0 +  10.0 *mh4*
                                                 ( 12.0 *lnbm*(-5.0 + complex<double>(0.0,27.0)*radix4mh2 +  30.0 *mh2) +
                                                   9.0 *power_of<2>(lnbm)*(1.0 +  18.0 *mh2 -  80.0 *power_of<2>(mh2) + complex<double>(0.0,1.0)*radix4mh2*(-1.0 +  40.0 *mh2)) +
                                                   3.0 *lnam*( 4.0 *(5.0 - complex<double>(0.0,27.0)*radix4mh2 -  30.0 *mh2) +
                                                       3.0 *lnbm*(-1.0 + complex<double>(0.0,1.0)*radix4mh2 -  18.0 *mh2 - complex<double>(0.0,40.0)*radix4mh2*mh2 +  80.0 *power_of<2>(mh2))) +
                                                   2.0 *(216.0 -  81.0 *mh2 + complex<double>(0.0,1.0)*radix4mh2*(101.0 +  90.0 *mh2 + complex<double>(0.0,264.0)*pi)))) +
                                             mh2*(- 8640.0 *power_of<2>(lnam) +  15.0 *lnbm*(-261.0 + complex<double>(0.0,288.0)*radix4mh2 +  27.0 *power_of<6>(radix4mh2) -  4.0 *mh2 +
                                                      1296.0 *power_of<2>(mh2)) - complex<double>(0.0,180.0)*power_of<2>(lnbm)*
                                                 ( 7.0 *power_of<5>(radix4mh2) + radix4mh2*(-7.0 +  52.0 *mh2) - complex<double>(0.0,4.0)*(12.0 + mh2 +  30.0 *power_of<2>(mh2))) +
                                                  60.0 *lnam*( 8.0 *(9.0 - complex<double>(0.0,9.0)*radix4mh2 -  10.0 *mh2) +
                                                      3.0 *lnbm*(complex<double>(0.0,7.0)*power_of<5>(radix4mh2) + complex<double>(0.0,1.0)*radix4mh2*(-7.0 +  52.0 *mh2) +  4.0 *(24.0 + mh2 +  30.0 *power_of<2>(mh2)))) +
                                                  16.0 *(-1201.0 +  216.0 *zeta3 -  30.0 *mh2 + radix4mh2*(complex<double>(0.0,135.0) + complex<double>(0.0,70.0)*mh2 +  876.0 *pi) - complex<double>(0.0,36.0)*pi3)))*
                                         rho4 -  2304.0 *(mh2*rho3 -  3.0 *mh4*( 2.0 *rho2 + (16.0 -  6.0 *zeta3 -  4.0 *radix4mh2*pi -  12.0 *radixrho*pi + complex<double>(0.0,1.0)*pi3)*rho3))) -
                                      4.0 *power_of<2>(complex<double>(0.0,-1.0) + radix4mh2)*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*mh8*
                                     (576.0 +  180.0 *( 3.0 *power_of<2>(lnam) -  11.0 *lnam*lnbm +  8.0 *power_of<2>(lnbm))*rho4 +  1890.0 *rho2 - complex<double>(0.0,620.0)*radixrho*rho2 -
                                      complex<double>(0.0,900.0)*power_of<3>(radixrho)*rho2 -  2400.0 *radixrho*pi*rho2 +  810.0 *(-1.0 + rho)*rho2 +  2745.0 *rho3 +
                                      complex<double>(0.0,900.0)*radixrho*rho3 -  3600.0 *radixrho*pi*rho3 -  45.0 *(-1.0 + rho)*rho3 -
                                      rho*( 8.0 *(-512.0 + complex<double>(0.0,35.0)*power_of<3>(radixrho) +  432.0 *zeta3 +  7.0 *radixrho*(complex<double>(0.0,5.0) +  96.0 *pi) - complex<double>(0.0,72.0)*pi3 -  15.0 *rho) +
                                           405.0 *lnam*rho3) + complex<double>(0.0,2.0)*radix4mh2*( 900.0 *(lnam - lnbm)*lnbm*rho4 +
                                               10.0 *(complex<double>(0.0,90.0)*power_of<3>(radixrho) + radixrho*(complex<double>(0.0,62.0) +  240.0 *pi) -  27.0 *(4.0 +  3.0 *rho))*rho2 +
                                              rho*( 8.0 *(-512.0 + complex<double>(0.0,35.0)*power_of<3>(radixrho) +  432.0 *zeta3 +  7.0 *radixrho*(complex<double>(0.0,5.0) +  96.0 *pi) - complex<double>(0.0,72.0)*pi3 -
                                                       15.0 *rho) +  405.0 *lnam*rho3) +  9.0 *(-64.0 +  5.0 *(-62.0 +  20.0 *radixrho*(complex<double>(0.0,-1.0) +  4.0 *pi) + rho)*rho3)) +
                                      (-1.0 +  4.0 *mh2)*( 900.0 *(lnam - lnbm)*lnbm*rho4 +
                                           10.0 *(complex<double>(0.0,90.0)*power_of<3>(radixrho) + radixrho*(complex<double>(0.0,62.0) +  240.0 *pi) -  27.0 *(4.0 +  3.0 *rho))*rho2 +
                                          rho*( 8.0 *(-512.0 + complex<double>(0.0,35.0)*power_of<3>(radixrho) +  432.0 *zeta3 +  7.0 *radixrho*(complex<double>(0.0,5.0) +  96.0 *pi) - complex<double>(0.0,72.0)*pi3 -
                                                   15.0 *rho) +  405.0 *lnam*rho3) +  9.0 *(-64.0 +  5.0 *(-62.0 +  20.0 *radixrho*(complex<double>(0.0,-1.0) +  4.0 *pi) + rho)*rho3))))/
                                                  (48.*power_of<2>(complex<double>(0.0,-1.0) + radixrho)*mh4*power_of<4>(- 4.0 *mh2 + rho)) +
                                                  atan4mh2*(( 288.0 *atanh4mh2rho*pi*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3))/power_of<4>(- 4.0 *mh2 + rho) +
                                                          ( 8.0 *(complex<double>(0.0,25.0)*mh2*rho4 +  382.0 *radix4mh2*mh2*rho4 -  382.0 *power_of<5>(radix4mh2)*mh2*rho4 +
                                                              complex<double>(0.0,25.0)*power_of<6>(radix4mh2)*mh2*rho4 - complex<double>(0.0,300.0)*power_of<2>(mh2)*rho4 -
                                                               3056.0 *radix4mh2*power_of<2>(mh2)*rho4 - complex<double>(0.0,240.0)*power_of<3>(mh2)*rho4 + complex<double>(0.0,150.0)*mh4*rho4 -
                                                               340.0 *radix4mh2*mh4*rho4 +  340.0 *power_of<5>(radix4mh2)*mh4*rho4 +
                                                              complex<double>(0.0,150.0)*power_of<6>(radix4mh2)*mh4*rho4 - complex<double>(0.0,1800.0)*mh2*mh4*rho4 +
                                                               2720.0 *radix4mh2*mh2*mh4*rho4 + complex<double>(0.0,5600.0)*power_of<2>(mh2)*mh4*rho4 + complex<double>(0.0,1305.0)*mh6*rho4 +
                                                               2640.0 *radix4mh2*mh6*rho4 -  1200.0 *power_of<5>(radix4mh2)*mh6*rho4 -
                                                              complex<double>(0.0,135.0)*power_of<6>(radix4mh2)*mh6*rho4 + complex<double>(0.0,1620.0)*mh2*mh6*rho4 -
                                                               960.0 *radix4mh2*mh2*mh6*rho4 - complex<double>(0.0,16080.0)*power_of<2>(mh2)*mh6*rho4 +
                                                              complex<double>(0.0,8640.0)*power_of<2>(mh2)*mh8*rho4 -  3072.0 *power_of<3>(mh2)*pi*rho4 +  13824.0 *power_of<2>(mh2)*mh4*pi*rho4 -
                                                               34560.0 *power_of<2>(mh2)*mh6*pi*rho4 +  38400.0 *power_of<2>(mh2)*mh8*pi*rho4 +  1152.0 *radix4mh2*mh6*rho2 -
                                                               1152.0 *power_of<5>(radix4mh2)*mh6*rho2 -  9216.0 *radix4mh2*mh2*mh6*rho2 -  18432.0 *power_of<2>(mh2)*mh6*pi*rho2 +
                                                               36864.0 *power_of<2>(mh2)*mh8*pi*rho2 -  576.0 *radix4mh2*mh4*rho3 +  576.0 *power_of<5>(radix4mh2)*mh4*rho3 +  4608.0 *radix4mh2*mh2*mh4*rho3 +
                                                               9216.0 *power_of<2>(mh2)*mh4*pi*rho3 -  18432.0 *power_of<2>(mh2)*mh6*pi*rho3 +
                                                               2304.0 *lnradices*power_of<2>(mh2)*pi*(- 64.0 *mh8*rho + mh2*rho4 +  48.0 *mh6*rho2 -  12.0 *mh4*rho3) +
                                                               2304.0 *ln2*mh4*pi*( 64.0 *mh8*rho - mh2*rho4 -  48.0 *mh6*rho2 +  12.0 *mh4*rho3)) +
                                                            2304.0 *lnam*(power_of<2>(complex<double>(0.0,-1.0) + radix4mh2)*mh8*rho*( 32.0 *power_of<2>(complex<double>(0.0,1.0) + radix4mh2)*pi + complex<double>(0.0,5.0)*rho3) -
                                                                2.0 *(complex<double>(0.0,10.0)*mh6*rho4 + mh2*( 5.0 *(complex<double>(0.0,-1.0) + radix4mh2)*mh6 +  2.0 *mh4*(complex<double>(0.0,-5.0) +  2.0 *pi))*rho4 +
                                                                    2.0 *power_of<2>(mh2)*(complex<double>(0.0,5.0)*mh6*rho4 +  96.0 *mh6*pi*rho2 -  24.0 *mh4*pi*rho3))) -
                                                            3.0 *lnbm*(- 5.0 *mh2*( 78.0 *power_of<5>(radix4mh2) - complex<double>(0.0,62.0)*power_of<6>(radix4mh2) +  54.0 *power_of<7>(radix4mh2) -
                                                                   complex<double>(0.0,15.0)*power_of<8>(radix4mh2) +  14.0 *power_of<9>(radix4mh2) + radix4mh2*(-38.0 -  896.0 *mh4 +  2304.0 *mh6) +
                                                                   complex<double>(0.0,1.0)*(-47.0 +  848.0 *mh4 -  1536.0 *mh6 -  3072.0 *mh8))*rho4 +
                                                                20.0 *((complex<double>(0.0,21.0) -  38.0 *radix4mh2 +  58.0 *power_of<5>(radix4mh2) + complex<double>(0.0,31.0)*power_of<6>(radix4mh2) +  20.0 *power_of<7>(radix4mh2) +
                                                                       complex<double>(0.0,10.0)*power_of<8>(radix4mh2))*mh4 -
                                                                    24.0 *((complex<double>(0.0,32.0) - radix4mh2 + power_of<5>(radix4mh2))*mh6 -  16.0 *(complex<double>(0.0,-1.0) + radix4mh2)*mh8))*rho4 +
                                                                1536.0 *power_of<3>(mh2)*(complex<double>(0.0,15.0) +  4.0 *pi)*rho4 -
                                                                8.0 *power_of<2>(mh2)*( 5.0 *( 25.0 *radix4mh2 + complex<double>(0.0,1.0)*(63.0 -  248.0 *mh4 +  576.0 *mh6))*rho4 +
                                                                    256.0 *mh8*( 192.0 *pi*rho + complex<double>(0.0,25.0)*rho4) +  9216.0 *(- 4.0 *mh6*pi*rho2 + mh4*pi*rho3))))
                / (64.*mh4*power_of<4>(- 4.0 *mh2 + rho)));
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J2
        complex<double>
        DileptonIntegralsBottom::j2(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            static const double ln2 = std::log(2.0);

            // Begin of asymptotic part
            complex<double> asymp = complex<double>(0,6)*dilogam2*radixrho + complex<double>(0,6)*dilogapbm*radixrho + 12.0 *atanh4mh2rho*radixrho*pi - 3.0 *lnrhom1*radixrho*pi -
                (complex<double>(0,6)*dilogambm*(-1 + rho))/radixrho - (12.0 *power_of<2>(atan4mh2)*(2.0 *mh2*(-2 + rho) + rho))/(4.0 *mh2 - rho) +
                (12.0 *power_of<2>(atanrho)*(2.0 *mh2*(radixrho*(-2 + rho) - complex<double>(0,4)*(-1 + rho)) + (radixrho + complex<double>(0,2)*(-1 + rho))*rho))/
                (radixrho*(4.0 *mh2 - rho)) + atan4mh2*(-24*atanh4mh2rho*radixrho +
                        (12.0 *(2.0 *mh2*pi*(radixrho*(-2 + rho) - complex<double>(0,2)*(-1 + rho)) + radixrho*(radix4mh2 + pi + complex<double>(0,1)*radixrho*pi)*rho))/
                        (radixrho*(4.0 *mh2 - rho))) + (48.0 *lnmh*mh2*pi + 24*lnradixrho*mh2*pi - 24*lnrho*mh2*pi - complex<double>(0,24)*mh2*power_of<2>(pi) +
                            complex<double>(0,28)*mh2*pi2 - 6*pi*rho - 12.0 *lnmh*pi*rho - 6*lnradixrho*pi*rho + 6*lnrho*pi*rho - 48*lnmh*mh2*pi*rho -
                            24*lnradixrho*mh2*pi*rho + 24*lnrho*mh2*pi*rho - 6.0 * radixrho * radix4mh2 * pi * rho + complex<double>(0,6)*power_of<2>(pi)*rho +
                            complex<double>(0,24)*mh2*power_of<2>(pi)*rho - complex<double>(0,7)*pi2*rho - complex<double>(0,28)*mh2*pi2*rho + 6*radixrho*(-4*mh2 + rho) +
                            24*ln2*power_of<2>(radixrho)*pi*(-4*mh2 + rho) + 6*lndeltarho4mh2*pi*(4.0 *mh2*(-1 + rho) + rho - rho2) + 6*pi*rho2 +
                            12.0 *lnmh*pi*rho2 + 6*lnradixrho*pi*rho2 - 6*lnrho*pi*rho2 - complex<double>(0,6)*power_of<2>(pi)*rho2 + complex<double>(0,7)*pi2*rho2)/
                        (radixrho*(4.0 *mh2 - rho)) + atanrho*(complex<double>(0,24)*atan4mh2*radixrho + 6*lnrhom1*radixrho +
                                (12.0 *((1.0 + 4*ln2 - lndeltarho4mh2 + 2*lnmh + lnradixrho - lnrho + complex<double>(0,1)*pi - radixrho*pi)*rho +
                                        2*mh2*(-2 + 2*lnrho - 2*power_of<2>(radixrho) + 8*ln2*power_of<2>(radixrho) - 2*lndeltarho4mh2*power_of<2>(radixrho) -
                                            complex<double>(0,2)*pi + 2*radixrho*pi + 4*lnmh*(-1 + rho) + 2*lnradixrho*(-1 + rho) + 2*rho - 2*lnrho*rho +
                                            complex<double>(0,2)*pi*rho - radixrho*pi*rho) +
                                        (-1 - 4*ln2 + lndeltarho4mh2 - 2*lnmh - lnradixrho + lnrho - complex<double>(0,1)*pi)*rho2))/(radixrho*(4.0 *mh2 - rho)));
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (complex<double>(0,18)*dilogam2*radixrho*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) -
                (complex<double>(0,18)*dilogambm*radixrho*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) +
                (complex<double>(0,18)*dilogapbm*radixrho*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) +
                (36.0 *atanh4mh2rho*radixrho*pi*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) -
                (9.0 *lnrhom1*radixrho*pi*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) +
                (36.0 *power_of<2>(atan4mh2)*(rho*(rho + mh2*(-8 + 6*rho)) - 4*mh4*(-4.0 + 2*rho + rho2)))/power_of<2>(-4*mh2 + rho) -
                (36.0 *power_of<2>(atanrho)*(rho*((radixrho + complex<double>(0,2)*(-1 + rho))*rho + 2*mh2*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-4.0 + 3*rho))) -
                                             4*mh4*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-4.0 + 2*rho + rho2))))/(radixrho*power_of<2>(-4*mh2 + rho)) +
                atan4mh2*((-72*atanh4mh2rho*radixrho*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) +
                        (9.0*(-8*power_of<2>(mh2)*pi*rho*(-4.0 - complex<double>(0,4)*radixrho + 3*rho) - 3*radix4mh2*(1.0 + power_of<2>(radix4mh2))*rho2 +
                              2*mh2*(complex<double>(0,-2)*(complex<double>(0,-1) + radixrho)*pi*power_of<2>(rho) + power_of<3>(radix4mh2)*(2.0 *rho - rho2) +
                                  radix4mh2*(2.0 *rho + rho2) + 8*mh4*pi*(-4.0 - complex<double>(0,4)*radixrho + 2*rho + rho2))))/(mh2*power_of<2>(-4*mh2 + rho)))\
                + (3.0*(48.0 *power_of<3>(radixrho)*mh2*(2.0 *mh4*(-2 + rho) + 3*mh2*rho) +
                            12.0 *power_of<2>(radixrho)*mh2*pi*rho*(-4*ln2*(16.0 *mh4 + rho*(-8*mh2 + rho)) + lndeltarho4mh2*(16.0 *mh4 - 8*mh2*rho + rho2)) +
                            3*radixrho*(48.0 *power_of<2>(mh2)*rho - 3*(1.0 + power_of<2>(radix4mh2))*rho*rho2 +
                                mh2*(-64*mh4 + 2*rho*(2.0 *(1.0 + power_of<2>(radix4mh2))*rho - power_of<2>(radix4mh2)*rho2))) -
                            2*mh2*rho*(96.0 *lnrho*mh4*pi + complex<double>(0,96)*mh4*power_of<2>(pi) - complex<double>(0,112)*mh4*pi2 + 24*mh2*pi*rho - 48*lnrho*mh2*pi*rho +
                                24*mh2*radixrho * radix4mh2 *pi*rho - 48*mh4*pi*rho - 96*lnrho*mh4*pi*rho -
                                complex<double>(0,48)*mh2*power_of<2>(pi)*rho - complex<double>(0,96)*mh4*power_of<2>(pi)*rho + complex<double>(0,56)*mh2*pi2*rho + complex<double>(0,112)*mh4*pi2*rho -
                                12.0 *pi*rho2 + 6*lnrho*pi*rho2 - 24*mh2*pi*rho2 + 48*lnrho*mh2*pi*rho2 - 12.0 * radixrho * radix4mh2 *pi*rho2 -
                                12.0 *mh2*radixrho * radix4mh2 * pi*rho2 + 48*mh4*pi*rho2 + complex<double>(0,6)*power_of<2>(pi)*rho2 +
                                complex<double>(0,48)*mh2*power_of<2>(pi)*rho2 - complex<double>(0,7)*pi2*rho2 - complex<double>(0,56)*mh2*pi2*rho2 + 18*pi*rho*rho2 - 6*pi*rho3 -
                                6*lnrho*pi*rho3 - complex<double>(0,6)*power_of<2>(pi)*rho3 + complex<double>(0,7)*pi2*rho3 +
                                12.0 *lnmh*pi*(16.0 *mh4*(-1 + rho) + 8*mh2*rho - rho2 - 8*mh2*rho2 + rho3) +
                                6*lnradixrho*pi*(16.0 *mh4*(-1 + rho) + 8*mh2*rho - rho2 - 8*mh2*rho2 + rho3))))/
                (2.*radixrho*mh2*rho*power_of<2>(-4*mh2 + rho)) +
                atanrho*((complex<double>(0,72)*atan4mh2*radixrho*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) +
                        (18.0 *lnrhom1*radixrho*(16.0 *mh4 + rho*(-8*mh2 + rho)))/power_of<2>(-4*mh2 + rho) +
                        (36.0*(4.0 *power_of<4>(radixrho)*(2.0 *mh4*(-2 + rho) + 3*mh2*rho) +
                               4*power_of<2>(radixrho)*(3.0 *mh2*rho + 2*mh4*(-2.0 +
                                       (3.0 + 8*ln2 - 2*lndeltarho4mh2 + 4*lnmh + 2*lnradixrho - 2*lnrho + complex<double>(0,2)*pi)*rho)) +
                               radixrho*pi*rho*(rho*(rho + mh2*(-8 + 6*rho)) - 4*mh4*(-4.0 + 2*rho + rho2)) +
                               rho*(8.0 *mh2*(2.0 + 4*ln2 - lndeltarho4mh2 + 2*lnmh + lnradixrho - lnrho + complex<double>(0,1)*pi)*(rho - rho2) +
                                   (-2 - 4*ln2 + lndeltarho4mh2 - 2*lnmh - lnradixrho + lnrho - complex<double>(0,1)*pi + 3*rho)*rho2 +
                                   (-1 + 4*ln2 - lndeltarho4mh2 + 2*lnmh + lnradixrho - lnrho + complex<double>(0,1)*pi)*rho3)))/
                        (radixrho*rho*power_of<2>(-4*mh2 + rho)));
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (complex<double>(0,36)*dilogam2*radixrho*(-64*mh6 + rho*(48*mh4 + rho*(-12*mh2 + rho))))/power_of<3>(-4*mh2 + rho) -
                (complex<double>(0,36)*dilogambm*radixrho*(-64*mh6 + rho*(48*mh4 + rho*(-12*mh2 + rho))))/power_of<3>(rho - 4.0 * mh2) +
                (complex<double>(0,36)*dilogapbm*radixrho*(-64*mh6 + rho*(48*mh4 + rho*(-12*mh2 + rho))))/power_of<3>(rho - 4.0 * mh2) +
                (72*atanh4mh2rho*radixrho*pi*(-64*mh6 + rho*(48*mh4 + rho*(-12*mh2 + rho))))/power_of<3>(rho - 4.0 * mh2) -
                (72*power_of<2>(atan4mh2)*(rho*(rho*(12*(rho - 1.0)*mh2 + rho) - 4.0 *mh4*(-12 + 4.0 *rho + 5.0 *rho2)) + 4.0 *mh6*(-16 + 8.0 *rho + 5.0 *rho3)))/
                power_of<3>(4.0 * mh2 - rho) + (72*power_of<2>(atanrho)*(rho*(rho*
                                (complex<double>(0,-24)*radixrho*mh2 + 12.0 *(rho - 1.0)*mh2 + rho + complex<double>(0,2)*radixrho*rho) +
                                4.0 *mh4*(12.0 + complex<double>(0,24)*radixrho - 4.0 *rho - 5.0 *rho2)) + 4.0 *mh6*(-16.0 - complex<double>(0,32)*radixrho + 8.0 *rho + 5.0 *rho3)))/
                power_of<3>(4.0 * mh2 - rho) + atan4mh2*((-144*atanh4mh2rho*radixrho*(-64*mh6 + rho*(48*mh4 + rho*(-12*mh2 + rho))))/
                        power_of<3>(rho - 4.0 * mh2) + (3.0 *(-96*mh8*pi*rho*(-12.0 - complex<double>(0,12)*radixrho + 4.0 *rho + 5.0 *rho2) -
                                3.0 *radix4mh2*(mh2*(4*(1 + power_of<2>(4.0 * mh2 - 1.0))*rho2 + 8.0 *(4.0 * mh2 - 1.0)*rho2 -
                                        5.0 *(-1.0 + power_of<2>(4.0 * mh2 - 1.0))*rho*rho2) - 3.0 *power_of<2>(1.0 + (4.0 * mh2 - 1.0))*rho3) +
                                mh4*(24.0*(complex<double>(0,-1) + radixrho)*pi*(12*radixrho*mh2 + complex<double>(0,1)*rho)*rho2 +
                                    8.0 *radix4mh2 * (4.0 * mh2 - 1.0)*(6*rho - 5.0 *rho3) + 3.0 *radix4mh2 * power_of<2>(4.0 * mh2 - 1.0)*(8*rho - 5.0 *rho3) + 3.0 *radix4mh2*(8*rho + 5.0 *rho3) +
                                    96.0 *mh6*pi*(-16.0 - complex<double>(0,16)*radixrho + 8.0 *rho + 5.0 *rho3))))/(mh4*power_of<3>(4.0 * mh2 - rho))) +
                atanrho*((complex<double>(0,144)*atan4mh2*radixrho*(-64*mh6 + rho*(48*mh4 + rho*(-12*mh2 + rho))))/power_of<3>(rho - 4.0 * mh2) +
                        (24.0*(-36.0*radixrho * (rho - 1.0)*mh2*pi*rho2*rho2 +
                               24.0 *power_of<3>(rho - 1.0)*(rho*(mh4*(4.0 - 5.0 *rho) - 3.0*mh2*rho) + mh6*(-8.0 + 5.0 *rho2)) +
                               16.0 *power_of<2>(rho - 1.0)*(3.0*rho*(4*mh4 - 3.0*mh2*rho) + 4.0 *mh6*(-6.0 + 5.0 *rho2)) +
                               8.0 *(rho - 1.0)*(3.0*rho*(-3.0*mh2*rho + mh4*(4 + 5.0 *rho)) +
                                   mh6*(-24.0 + (49 + 96.0 *ln2 - 24.0 *lndeltarho4mh2 + 48.0 *lnmh + 48.0 *lnradixrho - 24.0 *lnrho + complex<double>(0,24)*pi)*rho2)) -
                               3.0*radixrho*pi*rho2*(rho*(rho2 - 4.0 *mh4*(-12 + 4.0 *rho + 5.0 *rho2)) + 4.0 *mh6*(-16.0 + 8.0 *rho + 5.0 *rho3)) +
                               rho2*(-((5.0 + 12.0 *ln2 - 3.0*lndeltarho4mh2 + 6.0 *lnmh + 6.0 *lnradixrho - 3.0*lnrho + complex<double>(0,3)*pi)*rho4) +
                                   48.0 *mh4*(8 + 12.0 *ln2 - 3.0*lndeltarho4mh2 + 6.0 *lnmh + 6.0 *lnradixrho - 3.0*lnrho + complex<double>(0,3)*pi)*(rho - rho2) -
                                   96.0 *mh2*rho2 - 144.0 *ln2*mh2*rho2 + 36.0 *lndeltarho4mh2*mh2*rho2 - 72.0 *lnmh*mh2*rho2 - 72.0 *lnradixrho*mh2*rho2 +
                                   36.0 *lnrho*mh2*rho2 - complex<double>(0,36)*mh2*pi*rho2 - 18.0 *rho4 + 23.0*rho3 + 12.0 *ln2*rho3 - 3.0*lndeltarho4mh2*rho3 +
                                   6.0 *lnmh*rho3 + 6.0 *lnradixrho*rho3 - 3.0*lnrho*rho3 + 72.0 *mh2*rho3 + 144.0 *ln2*mh2*rho3 - 36.0 *lndeltarho4mh2*mh2*rho3 +
                                   72.0 *lnmh*mh2*rho3 + 72.0 *lnradixrho*mh2*rho3 - 36.0 *lnrho*mh2*rho3 + complex<double>(0,3)*pi*rho3 + complex<double>(0,36)*mh2*pi*rho3 +
                                   3.0*rho*(-5*rho2 + 8.0 *mh2*rho2 + 5.0 *rho3))))/(radixrho*power_of<3>(4.0 * mh2 - rho)*rho2)) +
                (-576*radixrho * power_of<2>(rho - 1.0)*mh4*(2*rho*(3.0*mh2*rho + mh4*(-4 + 5.0 *rho)) + mh6*(16 - 15.0 *rho2)) -
                 384.0 *(rho - 1.0)*mh4*rho2*(96*ln2*mh6*pi - 24.0 *lndeltarho4mh2*mh6*pi + 48.0 *lnmh*mh6*pi + 48.0 *lnradixrho*mh6*pi -
                     24.0 *lnrho*mh6*pi - complex<double>(0,24)*mh6*pi2 + complex<double>(0,28)*mh6*pi2 + 42.0 *mh4*pi*rho + 10.0 *mh6*pi*rho +
                     complex<double>(0,18)*mh4*pi2*rho - complex<double>(0,21)*mh4*pi2*rho - 9.0 *mh2*pi*rho2 - 15.0 *mh4*pi*rho2 + 15.0 *mh6*pi*rho2) -
                 192.0 *radixrho * (rho - 1.0)*mh4*(6*rho*(6*mh2*rho + mh4*(-8 + 5.0 *rho)) + mh6*(96 - 80.0 *rho2 + 15.0 *rho3)) +
                 24.0 *mh4*rho2*(-96*(8 + 12.0 *ln2 - 3.0*lndeltarho4mh2 + 6.0 *lnmh + 6.0 *lnradixrho - 3.0*lnrho)*mh4*pi*(rho - rho2) -
                     complex<double>(0,6)*pi2*(rho4 + 12.0 *mh2*(rho2 - rho3) - rho3) +
                     complex<double>(0,7)*pi2*(rho4 + 12.0 *mh2*(rho2 - rho3) - rho3) +
                     2.0 *pi*((5 + 12.0 *ln2 - 3.0 *lndeltarho4mh2 + 6.0 *lnmh + 6.0 *lnradixrho - 3.0 *lnrho)*rho4 + 18.0 *rho4 - 23.0 *rho3 -
                         12.0 *ln2*rho3 + 3.0 *lndeltarho4mh2*rho3 - 6.0 *lnmh*rho3 - 6.0 *lnradixrho*rho3 + 3.0 *lnrho*rho3 -
                         3.0 *rho*((-5 + 8.0 *mh2)*rho2 + 5.0 *rho3) + 12.0 *mh2*((8 + 12.0 *ln2 - 3.0 *lndeltarho4mh2 + 6.0 *lnmh + 6.0 *lnradixrho - 3.0 *lnrho)*rho2 +
                             3.0 *(-2 - 4.0 *ln2 + lndeltarho4mh2 - 2.0 *lnmh - 2.0 *lnradixrho + lnrho)*rho3))) +
                 radixrho*(-18*(1 + (4.0 * mh2 - 1.0))*rho2*
                     (mh2*rho*(4*(1 + (4.0 * mh2 - 1.0))*rho - 5.0 *(4.0 * mh2 - 1.0)*rho2) - 3.0 *(1 + (4.0 * mh2 - 1.0))*rho3) -
                     288.0 *mh8*(8*rho*(-2 + radix4mh2*pi*rho2) - 5.0 *radix4mh2*pi*rho5) +
                     mh4*(144*rho3 + 288.0 *(4.0 * mh2 - 1.0)*rho3 + 144.0 *power_of<2>(4.0 * mh2 - 1.0)*rho3 +
                         720.0 *radix4mh2*pi*rho5 - 25.0 *rho5 - 375.0 *(4.0 * mh2 - 1.0)*rho5 -
                         90.0 *power_of<2>(4.0 * mh2 - 1.0)*rho5 + 135.0 *(4.0 * mh2 - 1.0)*rho5 - 45.0 *power_of<2>(4.0 * mh2 - 1.0)*rho5 -
                         1104.0 *radix4mh2*pi*rho5 + 64.0 *mh6*(-144 + 25.0 *rho2 + 15.0 *rho3) +
                         12.0 *mh2*(-120*radix4mh2*pi*rho5 + 4.0 *radix4mh2*pi*rho2*(24*rho2 + 5.0 *rho3) +
                             rho2*(-288 - 5.0 *rho3 + 15.0 *(4.0 * mh2 - 1.0)*rho3)))))/
                             (4.*radixrho*mh4*power_of<3>(4.0 * mh2 - rho)*rho2);
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J3
        complex<double>
        DileptonIntegralsBottom::j3(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            static const double ln2 = std::log(2.0);

            // Begin of asymptotic part
            complex<double> asymp =(complex<double>(0,-24)*dilogam2*(-1 + rho)*(-4*mh4 + mh2*rho))/(radixrho*(4*mh2 - rho)*rho) +
   (complex<double>(0,24)*dilogambm*(-1 + rho)*(-4*mh4 + mh2*rho))/(radixrho*(4*mh2 - rho)*rho) -
   (complex<double>(0,24)*dilogapbm*(-1 + rho)*(-4*mh4 + mh2*rho))/(radixrho*(4*mh2 - rho)*rho) -
   (48*atanh4mh2rho*pi*(-1 + rho)*(-4*mh4 + mh2*rho))/(radixrho*(4*mh2 - rho)*rho) +
   (24*power_of<2>(atan4mh2)*(mh2*(-2 + rho)*rho - mh4*(-8 + 4*rho + rho2)))/((4*mh2 - rho)*rho) -
   (24*power_of<2>(atanrho)*(mh2*(radixrho*(-2 + rho) - complex<double>(0,4)*(-1 + rho))*rho -
        mh4*(complex<double>(0,-16)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2))))/(radixrho*(4*mh2 - rho)*rho) +
   atan4mh2*((96*atanh4mh2rho*(-1 + rho)*(-4*mh4 + mh2*rho))/(radixrho*(4*mh2 - rho)*rho) +
      (3.0*(-8*mh4*pi*(radixrho*(-2 + rho) - complex<double>(0,2)*(-1 + rho))*rho -
           radix4mh2 * (4.0 * mh2 - 1.0)*radixrho*(rho2 + mh2*(-4*rho + rho2)) +
           8*mh2*mh4*pi*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2)) + radix4mh2*radixrho*(-rho2 + mh2*(4*rho + rho2))))/
       (radixrho*mh2*(4*mh2 - rho)*rho)) + atanrho*
    ((complex<double>(0,-96)*atan4mh2*(-1 + rho)*(-4*mh4 + mh2*rho))/(radixrho*(4*mh2 - rho)*rho) +
      (6.0*(8*power_of<2>(rho - 1.0)*(mh4*(-4 + rho) + mh2*rho) - 8*(rho - 1.0)*(-(mh2*rho) + mh4*(4 + rho)) -
           4*radixrho*pi*rho*(-(mh2*(-2 + rho)*rho) + mh4*(-8 + 4*rho + rho2)) +
           rho*(16*mh4*(3 + 8*ln2 - 2*lndeltarho4mh2 + 4*lnmh + 4*lnradixrho - 2*lnrho + complex<double>(0,2)*pi)*(-1 + rho) +
              8*mh2*(2 + 4*ln2 - lndeltarho4mh2 + 2*lnmh + 2*lnradixrho - lnrho + complex<double>(0,1)*pi)*(rho - rho2) - rho2 +
              2*rho*rho2 - rho3)))/(radixrho*(4*mh2 - rho)*rho2)) +
   (mh2*(-48*(rho - 1.0)*pi*rho*(mh4*(-6 + rho) + mh2*rho) + 48*radixrho * (rho - 1.0)*(mh4*(-4 + rho) + mh2*rho) -
        3*radixrho*(64*mh4 + rho*((4 + (4.0 * mh2 - 1.0) - 2*radix4mh2*pi)*rho2 +
              4*mh2*(-4 - 4*rho + 4*radix4mh2*pi*rho - radix4mh2*pi*rho2))) -
        2.0*rho*(16*mh4*(3*(3 + 8*ln2 - 2*lndeltarho4mh2 + 4*lnmh + 4*lnradixrho - 2*lnrho)*pi - complex<double>(0,6)*pi2 +
              complex<double>(0,7)*pi2)*(-1 + rho) + 4*mh2*(complex<double>(0,6)*pi2*(-1 + rho)*rho - complex<double>(0,7)*pi2*(-1 + rho)*rho +
              6*(2 + 4*ln2 - lndeltarho4mh2 + 2*lnmh + 2*lnradixrho - lnrho)*pi*(rho - rho2)) + 3*pi*((-1 + 2*rho)*rho2 - rho3))))/
    (2.*radixrho*mh2*(4*mh2 - rho)*rho2);
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (complex<double>(0,72)*dilogam2*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) -
   (complex<double>(0,72)*dilogambm*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) +
   (complex<double>(0,72)*dilogapbm*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) +
   (144*atanh4mh2rho*pi*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) -
   (72*power_of<2>(atan4mh2)*(rho*(mh2*(-2 + rho)*rho + mh4*(16 - 8*rho - 3*rho2)) + 4*mh6*(-8 + 4*rho + rho2 + rho3)))/
    (rho*power_of<2>(-4*mh2 + rho)) + (72*power_of<2>(atanrho)*
      (rho*(mh2*(radixrho*(-2 + rho) - complex<double>(0,4)*(-1 + rho))*rho + mh4*(complex<double>(0,32)*(-1 + rho) + radixrho*(16 - 8*rho - 3*rho2))) +
        4.0*mh6*(complex<double>(0,-16)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2 + rho3))))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) +
   atanrho*((complex<double>(0,288)*atan4mh2*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) -
      (6.0*(-16*power_of<4>(radixrho)*(-3*rho*(-8*mh4 + mh2*rho) + 16*mh6*(-3 + rho2)) -
           24*power_of<6>(radixrho)*(-(rho*(mh2*rho + mh4*(-8 + 3*rho))) + 4*mh6*(-4 + rho + rho2)) +
           24*power_of<2>(radixrho)*(rho*(mh2*rho - mh4*(8 + 3*rho)) + 4*mh6*(4 + rho + rho2)) +
           rho2*(-64*mh6*(13 + 24*ln2 - 6*lndeltarho4mh2 + 12*lnmh + 12*lnradixrho - 6*lnrho + complex<double>(0,6)*pi)*(-1 + rho) +
              4*power_of<4>(rho) - 48*mh4*(9 + 16*ln2 - 4*lndeltarho4mh2 + 8*lnmh + 8*lnradixrho - 4*lnrho + complex<double>(0,4)*pi)*
               (rho - rho2) + 60*mh2*rho2 + 96*ln2*mh2*rho2 - 24*lndeltarho4mh2*mh2*rho2 + 48*lnmh*mh2*rho2 +
              48*lnradixrho*mh2*rho2 - 24*lnrho*mh2*rho2 + complex<double>(0,24)*mh2*pi*rho2 + 9*rho*rho2 - 48*mh2*rho*rho2 + 6*power_of<2>(rho2) -
              10*rho3 - 12*mh2*rho3 - 96*ln2*mh2*rho3 + 24*lndeltarho4mh2*mh2*rho3 - 48*lnmh*mh2*rho3 - 48*lnradixrho*mh2*rho3 +
              24*lnrho*mh2*rho3 - complex<double>(0,24)*mh2*pi*rho3 - 9*rho*rho3) +
           12.0*radixrho*pi*rho2*(rho*(mh2*(-2 + rho)*rho + mh4*(16 - 8*rho - 3*rho2)) + 4*mh6*(-8 + 4*rho + rho2 + rho3))))/
       (radixrho*rho3*power_of<2>(-4*mh2 + rho))) +
   atan4mh2*((-288*atanh4mh2rho*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)))/(radixrho*rho*power_of<2>(-4*mh2 + rho)) +
      (3.0*(96*mh4*pi*rho*(mh2*(radixrho*(-2 + rho) - complex<double>(0,2)*(-1 + rho))*rho +
              mh4*(complex<double>(0,16)*(-1 + rho) + radixrho*(16 - 8*rho - 3*rho2))) +
           3*power_of<5>(radix4mh2)*radixrho*(mh2*rho*(-8*rho + 3*rho2) + 4*mh4*(4*rho - rho2 - rho3) + rho3) -
           2*power_of<3>(radix4mh2)*radixrho*(24*mh2*rho2 - 3*rho3 + 16*mh4*(-3*rho + rho3)) +
           384*mh4*mh6*pi*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2 + rho3)) +
           3.0*radix4mh2*radixrho*(-(mh2*rho*(8*rho + 3*rho2)) + rho3 + 4*mh4*(4*rho + rho2 + rho3))))/
       (4.*radixrho*mh4*rho*power_of<2>(-4*mh2 + rho))) +
   (-288*power_of<5>(radixrho)*mh4*(rho*(mh2*rho + mh4*(-8 + 3*rho)) - 4*mh6*(-4 + rho + rho2)) -
      96*power_of<3>(radixrho)*mh4*(3*rho*(2*mh2*rho + mh4*(-16 + 3*rho)) - 4*mh6*(-24 + 3*rho + 11*rho2)) -
      12.0*mh4*rho2*(-768*lnmh*mh6*pi - 768*lnradixrho*mh6*pi + 384*lnrho*mh6*pi + complex<double>(0,384)*mh6*power_of<2>(pi) -
         complex<double>(0,448)*mh6*pi2 + 96*mh4*pi*rho + 384*lnmh*mh4*pi*rho + 384*lnradixrho*mh4*pi*rho - 192*lnrho*mh4*pi*rho -
         160*mh6*pi*rho + 768*lnmh*mh6*pi*rho + 768*lnradixrho*mh6*pi*rho - 384*lnrho*mh6*pi*rho - complex<double>(0,192)*mh4*power_of<2>(pi)*rho -
         complex<double>(0,384)*mh6*power_of<2>(pi)*rho + complex<double>(0,224)*mh4*pi2*rho + complex<double>(0,448)*mh6*pi2*rho - 4*pi*power_of<4>(rho) +
         96*ln2*pi*(-1 + rho)*(16*mh6 + rho*(-8*mh4 + mh2*rho)) - 36*mh2*pi*rho2 - 48*lnmh*mh2*pi*rho2 - 48*lnradixrho*mh2*pi*rho2 +
         24*lnrho*mh2*pi*rho2 - 24*mh4*pi*rho2 - 384*lnmh*mh4*pi*rho2 - 384*lnradixrho*mh4*pi*rho2 + 192*lnrho*mh4*pi*rho2 +
         64*mh6*pi*rho2 + complex<double>(0,24)*mh2*power_of<2>(pi)*rho2 + complex<double>(0,192)*mh4*power_of<2>(pi)*rho2 - complex<double>(0,28)*mh2*pi2*rho2 -
         complex<double>(0,224)*mh4*pi2*rho2 - 9*pi*rho*rho2 + 48*mh2*pi*rho*rho2 - 72*mh4*pi*rho*rho2 - 6*pi*power_of<2>(rho2) -
         24*lndeltarho4mh2*pi*(-1 + rho)*(16*mh6 - 8*mh4*rho + mh2*rho2) + 10*pi*rho3 - 12*mh2*pi*rho3 + 48*lnmh*mh2*pi*rho3 +
         48*lnradixrho*mh2*pi*rho3 - 24*lnrho*mh2*pi*rho3 + 96*mh6*pi*rho3 - complex<double>(0,24)*mh2*power_of<2>(pi)*rho3 +
         complex<double>(0,28)*mh2*pi2*rho3 + 9*pi*rho*rho3) + radixrho*
       (-288*power_of<2>(mh4)*rho*(-8 + radix4mh2*pi*rho*(4*rho - rho2 - rho3)) +
         18*power_of<2>(mh2)*power_of<2>(rho)*(-32*mh2*power_of<2>(rho) + 3*power_of<2>(radix4mh2)*rho*rho2 + 4*rho3) +
         mh4*(512*mh6*(4*rho2 - 3*(3 + rho3)) - power_of<2>(rho)*
             (18*power_of<2>(radix4mh2)*rho2 + 18*power_of<4>(radix4mh2)*rho2 -
               36*rho*(2 + 4*power_of<2>(radix4mh2) + 2*power_of<4>(radix4mh2) + 3*radix4mh2*pi*rho2) + 5*rho3 + 48*power_of<2>(radix4mh2)*rho3 +
               27*power_of<4>(radix4mh2)*rho3 + 120*radix4mh2*pi*rho3 -
               12*mh2*(-24 - rho3 + 3*power_of<2>(radix4mh2)*rho3 + 2*radix4mh2*pi*(-9*(-2 + rho)*rho2 + 2*rho3))))))/
    (4.*radixrho*mh4*rho3*power_of<2>(-4*mh2 + rho));
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (complex<double>(0,144)*dilogam2*(-1 + rho)*(-64*mh8 + rho*(48*mh6 + rho*(-12*mh4 + mh2*rho))))/(radixrho*rho*power_of<3>(rho - 4.0 * mh2)) -
   (complex<double>(0,144)*dilogambm*(-1 + rho)*(-64*mh8 + rho*(48*mh6 + rho*(-12*mh4 + mh2*rho))))/(radixrho*rho*power_of<3>(rho - 4.0 * mh2)) +
   (complex<double>(0,144)*dilogapbm*(-1 + rho)*(-64*mh8 + rho*(48*mh6 + rho*(-12*mh4 + mh2*rho))))/(radixrho*rho*power_of<3>(rho - 4.0 * mh2)) +
   (288*atanh4mh2rho*pi*(-1 + rho)*(-64*mh8 + rho*(48*mh6 + rho*(-12*mh4 + mh2*rho))))/(radixrho*rho*power_of<3>(rho - 4.0 * mh2)) +
   (144*power_of<2>(atan4mh2)*(-(mh8*(64*rho + 25*rho4 + 16*(-8 + rho2))) +
        rho*(rho*(mh2*(-2 + rho)*rho - 6*mh4*(-4 + 2*rho + rho2)) + 4*mh6*(-24 + 12*rho + 2*rho2 + 5*rho3))))/(power_of<3>(4.0 * mh2 - rho)*rho) -
   (144*power_of<2>(atanrho)*(-(mh8*(complex<double>(0,-256)*(-1 + rho) + radixrho*(64*rho + 25*rho4 + 16*(-8 + rho2)))) +
        rho*(rho*(mh2*(radixrho*(-2 + rho) - complex<double>(0,4)*(-1 + rho))*rho -
              6.0*mh4*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-4 + 2*rho + rho2))) +
           4*mh6*(complex<double>(0,-48)*(-1 + rho) + radixrho*(-24 + 12*rho + 2*rho2 + 5*rho3)))))/(radixrho*power_of<3>(4.0 * mh2 - rho)*rho) +
   atan4mh2*((-576*atanh4mh2rho*(-1 + rho)*(-64*mh8 + rho*(48*mh6 + rho*(-12*mh4 + mh2*rho))))/(radixrho*rho*power_of<3>(rho - 4.0 * mh2)) +
      (3.0*(mh6*(384*mh8*pi*(complex<double>(0,-128)*(-1 + rho) + radixrho*(64*rho + 25*rho4 + 16*(-8 + rho2))) +
              radix4mh2 * power_of<2>(4.0 * mh2 - 1.0)*radixrho*(576*rho - 275*rho4 - 48*rho2) -
              3*radix4mh2 * power_of<3>(4.0 * mh2 - 1.0)*radixrho*(-64*rho + 25*rho4 + 16*rho2) +
              3*radix4mh2*radixrho*(64*rho + 25*rho4 + 16*rho2) +
              radix4mh2 * (4.0 * mh2 - 1.0)*radixrho*(576*rho - 365*rho4 + 48*rho2) -
              384*pi*rho2*(mh2*(radixrho*(-2 + rho) - complex<double>(0,2)*(-1 + rho))*rho -
                 6*mh4*(complex<double>(0,-4)*(-1 + rho) + radixrho*(-4 + 2*rho + rho2)))) -
           1536*mh12*pi*rho*(complex<double>(0,-24)*(-1 + rho) + radixrho*(-24 + 12*rho + 2*rho2 + 5*rho3)) -
           radix4mh2*(1 + (4.0 * mh2 - 1.0))*radixrho*
            (3*(1 + (4.0 * mh2 - 1.0))*((1 + (4.0 * mh2 - 1.0))*rho4 - 12*(1 + (4.0 * mh2 - 1.0))*mh2*rho*rho2 +
                 6*(-1 + (4.0 * mh2 - 1.0))*mh2*rho4) +
              4*mh4*rho*(36*(16.0 * mh4)*rho - 6*(-1 + power_of<2>(4.0 * mh2 - 1.0))*rho2 -
                 5*(-3 + 8*(4.0 * mh2 - 1.0) + 3*power_of<2>(4.0 * mh2 - 1.0))*rho3))))/(8.*radixrho*mh6*power_of<3>(4.0 * mh2 - rho)*rho)) +
   atanrho*((complex<double>(0,576)*atan4mh2*(-1 + rho)*(-64*mh8 + rho*(48*mh6 + rho*(-12*mh4 + mh2*rho))))/
       (radixrho*rho*power_of<3>(rho - 4.0 * mh2)) + (6.0*(rho3*
            (768*mh8*(11 + 16*ln2 - 4*lndeltarho4mh2 + 8*lnmh + 8*lnradixrho - 4*lnrho + complex<double>(0,4)*pi)*(-1 + rho) +
              3*rho4 - 48*mh2*rho4 - 192*ln2*mh2*rho4 + 48*lndeltarho4mh2*mh2*rho4 -
              96*lnmh*mh2*rho4 - 96*lnradixrho*mh2*rho4 + 48*lnrho*mh2*rho4 -
              complex<double>(0,48)*mh2*pi*rho4 + 25*rho5 +
              256*mh6*(25 + 36*ln2 - 9*lndeltarho4mh2 + 18*lnmh + 18*lnradixrho - 9*lnrho + complex<double>(0,9)*pi)*(rho - rho2) -
              1632*mh4*rho2 - 2304*ln2*mh4*rho2 + 576*lndeltarho4mh2*mh4*rho2 - 1152*lnmh*mh4*rho2 - 1152*lnradixrho*mh4*rho2 +
              576*lnrho*mh4*rho2 - complex<double>(0,576)*mh4*pi*rho2 - 48*mh2*rho*rho2 + 576*mh4*rho*rho2 + 36*rho4 -
              144*mh2*rho4 + 192*mh2*rho3 + 192*ln2*mh2*rho3 - 48*lndeltarho4mh2*mh2*rho3 + 96*lnmh*mh2*rho3 +
              96*lnradixrho*mh2*rho3 - 48*lnrho*mh2*rho3 + 1056*mh4*rho3 + 2304*ln2*mh4*rho3 - 576*lndeltarho4mh2*mh4*rho3 +
              1152*lnmh*mh4*rho3 + 1152*lnradixrho*mh4*rho3 - 576*lnrho*mh4*rho3 + complex<double>(0,48)*mh2*pi*rho3 +
              complex<double>(0,576)*mh4*pi*rho3 - 40*rho*rho3 + 48*mh2*rho*rho3 - 24*rho2*rho3) +
           48*power_of<4>(rho - 1.0)*(rho*(rho*(6*mh4*(-2 + rho) + mh2*rho) - 4*mh6*(-12 + 2*rho + 5*rho2)) + mh8*(-64 + 16*rho + 25*rho3)) -
           48*(-1 + rho)*(-(rho*(rho*(mh2*rho - 6*mh4*(2 + rho)) + 4*mh6*(12 + 2*rho + 5*rho2))) + mh8*(64 + 16*rho + 25*rho3)) +
           16*power_of<3>(rho - 1.0)*(rho*(9*rho*(2*mh4*(-6 + rho) + mh2*rho) - 4*mh6*(-108 + 6*rho + 55*rho2)) +
              mh8*(-576 + 48*rho + 275*rho3)) + 16*power_of<2>(rho - 1.0)*
            (rho*(9*rho*(mh2*rho - 2*mh4*(6 + rho)) + 4*mh6*(108 + 6*rho - 25*rho2)) + mh8*(-576 - 48*rho + 365*rho3)) -
           24*radixrho*pi*rho3*(mh8*(64*rho + 25*rho4 + 16*(-8 + rho2)) -
              rho*(rho*(mh2*(-2 + rho)*rho - 6*mh4*(-4 + 2*rho + rho2)) + 4*mh6*(-24 + 12*rho + 2*rho2 + 5*rho3)))))/
       (radixrho*power_of<3>(4.0 * mh2 - rho)*rho4)) +
   (-9216*radixrho * power_of<3>(rho - 1.0)*mh6*(rho*(-(rho*(6*mh4*(-2 + rho) + mh2*rho)) + 4*mh6*(-12 + 2*rho + 5*rho2)) + mh8*(64 - 16*rho - 25*rho3)) +
      1920*power_of<2>(rho - 1.0)*mh6*(-72*mh6 + mh8*(-14 + 3*rho))*(complex<double>(0,-1)*(rho - 1.0) + complex<double>(0.0, rho - 1.0))*rho3 -
      512*radixrho * (rho - 1.0)*mh6*(2*rho*(-27*rho*(2*mh4*(-6 + rho) + mh2*rho) +
            2*mh6*(-648 + 36*rho + 400*rho2 - 45*(rho - 1.0)*rho2 + complex<double>(0,45)*complex<double>(0.0, rho - 1.0)*rho2 - 30*rho3)) +
         mh8*(3456 - 288*rho + 15*rho4*(40 + 3*(rho - 1.0) - complex<double>(0,3)*complex<double>(0.0, rho - 1.0)) - 3050*rho3 -
            240*(rho - 1.0)*rho3 + complex<double>(0,240)*complex<double>(0.0, rho - 1.0)*rho3)) -
      1536*radixrho * power_of<2>(rho - 1.0)*mh6*(-3*mh8*(-384 + 64*rho + 5.0*(36 + (rho - 1.0) - complex<double>(0,1)*complex<double>(0.0, rho - 1.0))*rho3) -
         2*rho*(9*rho*(4*mh4*(-3 + rho) + mh2*rho) - mh6*
             (48*rho + 5.0*(62 + 3*(rho - 1.0) - complex<double>(0,3)*complex<double>(0.0, rho - 1.0))*rho2 - 6*(72 + 5*rho3)))) +
      192*rho3*(complex<double>(0,60)*mh6*(-1 + rho)*(8*mh6 + mh8*(26 + 3*rho))*((rho - 1.0) + complex<double>(0,1)*complex<double>(0.0, rho - 1.0)) +
         8*mh4*mh4*(complex<double>(0,-6)*pi2*(rho4 - rho3) + complex<double>(0,7)*pi2*(rho4 - rho3) +
            6.0*pi*((1 + 4*ln2 - lndeltarho4mh2 + 2*lnmh + 2*lnradixrho - lnrho)*rho4 + 3*rho4 +
               rho*(rho2 - rho3) + (-4 - 4*ln2 + lndeltarho4mh2 - 2*lnmh - 2*lnradixrho + lnrho +
                  complex<double>(0,1)*complex<double>(0.0, rho - 1.0))*rho3)) +
         mh2*mh4*(-3*pi*rho4 - 25*pi*rho5 + 1632*mh4*pi*rho2 + 2304*ln2*mh4*pi*rho2 - 576*lndeltarho4mh2*mh4*pi*rho2 +
            1152*lnmh*mh4*pi*rho2 + 1152*lnradixrho*mh4*pi*rho2 - 576*lnrho*mh4*pi*rho2 - complex<double>(0,576)*mh4*pi2*rho2 +
            complex<double>(0,672)*mh4*pi2*rho2 - 576*mh4*pi*rho*rho2 - complex<double>(0,1152)*mh4*pi*complex<double>(0.0, rho - 1.0)*rho2 - 36*pi*rho4 +
            16*mh6*(complex<double>(0,-144)*pi2*(-1 + rho)*rho +
               3.0*(complex<double>(0,-5)*(rho - 1.0) + complex<double>(0,56)*pi2*(-1 + rho)*rho + 5.0*complex<double>(0.0, rho - 1.0)) +
               4.0*pi*(4.0*(25.0 + 36*ln2 - 9*lndeltarho4mh2 + 18*lnmh + 18*lnradixrho - 9*lnrho - complex<double>(0,4)*complex<double>(0.0, rho - 1.0))*
                   rho2 + rho*(-100 - 144*ln2 + 36*lndeltarho4mh2 - 72*lnmh - 72*lnradixrho + 36*lnrho +
                     complex<double>(0,88)*complex<double>(0.0, rho - 1.0) - complex<double>(0,15)*complex<double>(0.0, rho - 1.0)*rho2))) - 1056*mh4*pi*rho3 -
            2304*ln2*mh4*pi*rho3 + 576*lndeltarho4mh2*mh4*pi*rho3 - 1152*lnmh*mh4*pi*rho3 - 1152*lnradixrho*mh4*pi*rho3 +
            576*lnrho*mh4*pi*rho3 + complex<double>(0,576)*mh4*pi2*rho3 - complex<double>(0,672)*mh4*pi2*rho3 + 40*pi*rho*rho3 +
            complex<double>(0,288)*mh4*pi*complex<double>(0.0, rho - 1.0)*rho3 + 24*pi*rho2*rho3 +
            2*mh8*(complex<double>(0,1536)*pi2*(-1 + rho) - complex<double>(0,1792)*pi2*(-1 + rho) +
               5.0*(34 + 3*rho)*(complex<double>(0,-1)*(rho - 1.0) + complex<double>(0.0, rho - 1.0)) +
               8*pi*(528 + 384*lnmh + 384*lnradixrho - 192*lnrho - 768*ln2*(-1 + rho) + 192*lndeltarho4mh2*(-1 + rho) - 528*rho -
                  384*lnmh*rho - 384*lnradixrho*rho + 192*lnrho*rho - complex<double>(0,528)*complex<double>(0.0, rho - 1.0) +
                  complex<double>(0,88)*rho*complex<double>(0.0, rho - 1.0) + complex<double>(0,50)*complex<double>(0.0, rho - 1.0)*rho2 +
                  complex<double>(0,75)*complex<double>(0.0, rho - 1.0)*rho3)))) +
      radixrho*(27648*mh8*rho4*rho2 -
         1152*mh6*(96*mh4*rho5 + 2*rho7 + 3*(4.0 * mh2 - 1.0)*rho3*rho4) +
         1024.0*mh12*rho*(432.0 + 5.0*(-10 + 27*(rho - 1.0) - complex<double>(0,27)*complex<double>(0.0, rho - 1.0))*rho2 - 30*rho3) +
         6*radix4mh2*mh4*mh4*(5.0*(complex<double>(0,-3) - 12*radix4mh2 + complex<double>(0,18)*(4.0 * mh2 - 1.0) + 12*radix4mh2 * (4.0 * mh2 - 1.0) -
               complex<double>(0,3)*power_of<2>(4.0 * mh2 - 1.0) + 64*pi)*rho4 + 256*pi*rho*(6*rho2 - 5*rho3) + 2304*pi*(rho4 - 2*rho3))*
          rho3 + 180*radix4mh2*(complex<double>(0,-1) + radix4mh2)*power_of<4>(complex<double>(0,1) + radix4mh2)*mh4*rho*rho6 +
         mh2*mh4*(800*rho4*rho3 - complex<double>(0,990)*radix4mh2 * (4.0 * mh2 - 1.0)*rho4*rho3 +
            complex<double>(0,105)*radix4mh2 * power_of<2>(4.0 * mh2 - 1.0)*rho4*rho3 - 90*power_of<3>(4.0 * mh2 - 1.0)*rho4*rho3 +
            48*power_of<2>(4.0 * mh2 - 1.0)*rho4*(24*rho2 + 65*rho3) + 6*(4.0 * mh2 - 1.0)*rho4*(192*rho2 + 1795*rho3) +
            3*radix4mh2*rho3*((complex<double>(0,155) + 64*(-3 + 50*mh4)*pi)*rho4 + 768*pi*(40*mh4 - 3*rho2)*rho2 -
               512*pi*rho*(-5*rho3 + 6*mh4*(2*rho2 + 5*rho3)))) +
         mh6*(512*mh8*(-1152.0 + 15*rho4*(8 + 3*(rho - 1.0) - complex<double>(0,3)*complex<double>(0.0, rho - 1.0)) -
               5.0*(-40 + 93*(rho - 1.0) - complex<double>(0,93)*complex<double>(0.0, rho - 1.0))*rho3) +
            2304*mh4*(-48*rho2 - 64*radix4mh2*mh2*pi*rho*rho3 + 25*radix4mh2*mh2*pi*rho4*rho3 +
               16*radix4mh2*mh2*pi*rho2*rho3) + rho3*
             (-5.0*(80.0 - complex<double>(0,51)*radix4mh2 + 1034*(4.0 * mh2 - 1.0) + complex<double>(0,234)*radix4mh2 * (4.0 * mh2 - 1.0) +
                  744*power_of<2>(4.0 * mh2 - 1.0) + complex<double>(0,21)*radix4mh2 * power_of<2>(4.0 * mh2 - 1.0) + 198*power_of<3>(4.0 * mh2 - 1.0))*rho4 -
               576*(4.0 * mh2 - 1.0) * 16.0 * mh4*rho2 +
               6*mh2*(1536.0 + 5.0*(-32.0 + complex<double>(0,3)*radix4mh2 + 148*(4.0 * mh2 - 1.0) - complex<double>(0,18)*radix4mh2 * (4.0 * mh2 - 1.0) +
                     12.0*power_of<2>(4.0 * mh2 - 1.0) + complex<double>(0,3)*radix4mh2 * power_of<2>(4.0 * mh2 - 1.0))*rho4 +
                  64*rho*(24 + 24*power_of<2>(4.0 * mh2 - 1.0) + (4.0 * mh2 - 1.0)*(48 - 15*rho3) + 5*rho3))))))/
    (64.*radixrho*mh6*power_of<3>(4.0 * mh2 - rho)*rho4);
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J4
        complex<double>
        DileptonIntegralsBottom::j4(const double & a1, const double & a2) const
        {
            static const double pi = M_PI;
            const double acotrho = pi / 2.0 - atanrho, acot4mh2 = pi / 2.0 - atan4mh2;

            // Begin of the asymptotic part
            double asymp = (32*power_of<2>(acotrho)*mh4*(4*mh2*(-3 + rho) - 3*rho)*rho2)/(3.*power_of<3>(4.0 * mh2 - rho)) +
   (32*power_of<2>(acot4mh2)*mh4*mh6*(-4*mh2*(-3 + rho) + 3*rho)*rho2)/(3.*mh6*power_of<3>(4.0 * mh2 - rho)) -
   (64*acotrho*radixrho*mh4*(3*rho*(-2 + 5*rho) + 4*mh2*(2 + rho + 3*rho2)))/(9.*power_of<3>(4.0 * mh2 - rho)) -
   (4*lnmqmu*mh4*(16*mh6*(4*mh2 - 3*rho) + rho*(12*mh4*rho - mh2*rho2)))/(9.*mh6*power_of<3>(4.0 * mh2 - rho)) +
   (8*acot4mh2*radix4mh2*mh4*rho*(6*(16*mh6 + mh4*(-4 + rho))*rho - 4*mh4*(-3 + rho)*rho + mh2*(rho2 + 24*mh4*(rho + rho2))))/
    (9.*mh6*power_of<3>(4.0 * mh2 - rho)) + (2*mh4*(16*mh6*(rho*(-13 + 15*rho) + 4*mh2*(3 + 8*rho + 3*rho2)) +
        rho*(-3*mh2*rho2 - 4*mh4*(-13*rho + 60*mh2*rho + 8*rho2 + 12*mh2*rho2))))/(9.*mh6*power_of<3>(4.0 * mh2 - rho));
            // End of the asymptotic part

            // Begin of the 1st Gegenbauer moment
            double gb1 = (96*power_of<2>(acot4mh2)*mh4*((4*mh2*(-4 + rho) - rho)*rho + 2*mh4*(-8 + 8*rho - 3*rho2))*rho2)/power_of<4>(rho - 4.0 * mh2) +
   (96*power_of<2>(acotrho)*mh4*rho2*(rho*(-4*mh2*(-4 + rho) + rho) + mh4*(16 - 16*rho + 6*rho2)))/power_of<4>(rho - 4.0 * mh2) +
   (64*acotrho*radixrho*mh4*((-2 + 24*mh2)*rho2 + 36*mh4*rho2 + 12*mh2*rho*rho2 + 5*rho3 - 18*mh4*rho3))/power_of<4>(rho - 4.0 * mh2) -
   (4*lnmqmu*mh4*rho*(16*mh8*rho - 8*mh6*rho*(4*mh2 + rho) + mh4*rho*(16*mh4 + rho*(16*mh2 + rho)) - 2*mh4*(4*mh2 + rho)*rho2 + mh4*rho3))/
    (mh8*power_of<4>(rho - 4.0 * mh2)) - (mh4*(24*mh4*rho*(mh2*(4 - 3*rho) + rho)*rho2 +
        2*mh4*rho2*(7*rho*(-rho + 4*mh2*(-4 + 3*rho)) + 2*mh4*(-56 + 168*rho + 9*rho2)) +
        32*mh8*(rho*(rho*(-13 + 15*rho) + 4*mh2*(-4 + 39*rho + 9*rho2)) + 2*mh4*(4 + 24*rho + 63*rho2 - 27*rho3)) -
        8*mh6*rho*(144*mh4*rho2 + 5*(-4 + 3*rho)*rho2 + 4*mh2*rho*(-20 + 60*rho + 9*rho2) + 6*mh4*(40*rho - 9*rho3)) - 11*mh4*rho*rho3))/
    (3.*mh8*power_of<4>(rho - 4.0 * mh2)) - (8*acot4mh2*radix4mh2*mh4*rho*
      (8*(-1 + 4*mh2)*mh6*rho*(4*mh2 + rho) + 2*mh8*(48*mh2*rho2 + 4*rho*(2 + 4*mh2 + 3*rho2) - 3*rho3) +
        mh4*(16*mh4*rho2 - 2*rho*rho2 + 2*rho3 - 36*mh6*rho3 + 2*mh2*(4*rho2 - 2*rho*rho2 + 2*rho3))))/(mh8*power_of<4>(rho - 4.0 * mh2));
            // End of the 1st Gegenbauer moment

            // Begin of the 2nd Gegenbauer moment
            double gb2 =(-192*power_of<2>(acotrho)*mh4*rho2*(rho*(rho*(36*mh2 + rho - 8*mh2*rho) + mh4*(144 - 96*rho + 30*rho2)) -
        8*mh6*(-8 + 16*rho - 15*rho2 + 5*rho3)))/power_of<5>(4.0 * mh2 - rho) -
   (192*power_of<2>(acot4mh2)*mh4*rho2*(rho*(rho*(-36*mh2 - rho + 8*mh2*rho) + 6*mh4*(-24 + 16*rho - 5*rho2)) +
        8*mh6*(-8 + 16*rho - 15*rho2 + 5*rho3)))/power_of<5>(4.0 * mh2 - rho) +
   (8*lnmqmu*mh4*rho*(16*mh10*rho*(4*mh2 + rho) - 8*mh8*rho*(16*mh4 + rho*(12*mh2 + rho)) - 4*mh2*mh4*(16*mh4 + rho*(12*mh2 + rho))*rho2 +
        mh6*rho*(64*mh6 + rho*(144*mh4 + 36*mh2*rho + rho2)) + 5*mh2*mh4*(4*mh2 + rho)*rho3 - 2*mh6*rho4))/(mh10*power_of<5>(4.0 * mh2 - rho))
     - (128*acotrho*radixrho*mh4*(228*mh4*rho*rho2 + 6*mh4*rho*(16*rho - 15*rho3) - 2*rho3 +
        4*mh2*(52*mh4*rho2 + rho2*(-2 + 6*rho2) + 17*rho3 - 70*mh4*rho3) + 5*rho4 + 120*mh6*rho4))/power_of<5>(4.0 * mh2 - rho) +
   (2*mh4*(-40*mh2*mh4*rho*rho2*(3*rho*(-2*rho + 3*mh2*(-8 + 5*rho)) + 4*mh4*(-24 + 45*rho + rho2)) +
        10*mh6*rho2*(rho*(7*rho*(-36*mh2 - rho + 24*mh2*rho) + 18*mh4*(-56 + 112*rho + 5*rho2)) +
           8*mh6*(-56 + 336*rho + 45*rho2 - 15*rho3)) + 55*mh2*mh4*rho*(-5*rho + 4*mh2*(-5 + 3*rho))*rho3 + 104*mh6*rho*rho4 -
        40*mh8*rho*(rho*(864*mh4*rho2 + 5*(-4 + 3*rho)*rho2 + 12*mh2*rho*(-20 + 45*rho + 6*rho2)) +
           2*mh4*(1080*rho2 + 36*mh2*(16*rho2 - 15*rho3) - 5*rho*(32 + 27*rho3)) + 120*mh6*(8*rho + 3*rho4)) +
        32*mh10*(5*rho*(rho*(rho*(-13 + 15*rho) + 12*mh2*(-7 + 31*rho + 6*rho2)) + 6*mh4*(-4 + 144*rho + 129*rho2 - 45*rho3)) +
           8*mh6*(4 + 60*rho + 535*rho2 - 600*rho3 + 225*rho4))))/(15.*mh10*power_of<5>(4.0 * mh2 - rho)) +
   (16*acot4mh2*radix4mh2*mh4*rho*(768*mh2*mh12*rho2 + 8*(-1 + 4*mh2)*mh8*(12*mh2 + rho)*rho2 +
        2*mh2*mh8*(8*rho2*(-4 + 4*mh2*(-2 + 9*rho) + 3*rho2) - 15*(4*mh2 + rho)*rho3) +
        mh4*(128*(-1 + 4*mh2)*mh8*rho - mh2*(1 + 2*mh2)*(48*mh2*rho*rho2 + 4*rho4 - 20*mh2*rho3 - 5*rho*rho3)) +
        8*mh12*(8*(1 + 2*mh2)*rho + 5*rho4) + mh6*
         (72*mh4*rho3 + 2*mh2*(18*rho3 - rho4) - rho4 + 240*mh8*rho4 +
           4*mh4*(36*(1 + 2*mh2)*rho2 - 45*mh2*rho*rho3 + 2*(-90*mh4*rho3 + rho4)))))/(mh10*power_of<5>(4.0 * mh2 - rho));
            // End of the 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J5
        complex<double>
        DileptonIntegralsBottom::j5(const double & a1, const double & a2) const
        {
            static const double pi = M_PI;
            const double acotrho = pi / 2.0 - atanrho, acot4mh2 = pi / 2.0 - atan4mh2;
            const double ln4mh2 = 2.0 * (std::log(2.0) + lnmh);

            double asymp = (2*mh2*rho*(208*mh2*mh4 - 96*lnmqmu*mh2*mh4 - 192*acotrho*radixrho*mh2*mh4 + 96*power_of<2>(acot4mh2)*mh4*rho -
       96*power_of<2>(acotrho)*mh4*rho + 256*acot4mh2*radix4mh2*mh4*rho - 256*acotrho*radixrho*mh4*rho - 368*mh2*mh4*rho +
       128*acot4mh2*radix4mh2*mh2*mh4*rho + 480*acotrho*radixrho*mh2*mh4*rho - 384*acot4mh2*radix4mh2*mh6*rho - 13*mh2*rho2 +
       6*lnmqmu*mh2*rho2 + 12*acot4mh2*radix4mh2*mh2*rho2 + 92*mh4*rho2 - 120*acot4mh2*radix4mh2*mh4*rho2 +
       64*acotrho*radixrho*mh4*rho2 - 144*power_of<2>(acot4mh2)*mh2*mh4*rho2 + 144*power_of<2>(acotrho)*mh2*mh4*rho2 +
       48*(1 + lnmqmu)*mh4*rho*ln4mh2 - 48*(1 + lnmqmu)*mh4*rho*lnrho))/(9.*mh4*power_of<3>(4.0 * mh2 - rho));

            double gb1 = (2*mh2*rho*(-448*mh4*mh6 + 128*lnmqmu*mh4*mh6 + 256*acotrho*radixrho*mh4*mh6 - 944*mh2*mh6*rho - 384*power_of<2>(acot4mh2)*mh2*mh6*rho +
       384*power_of<2>(acotrho)*mh2*mh6*rho + 480*lnmqmu*mh2*mh6*rho - 640*acot4mh2*radix4mh2*mh2*mh6*rho +
       1600*acotrho*radixrho*mh2*mh6*rho + 2496*mh4*mh6*rho - 512*acot4mh2*radix4mh2*mh4*mh6*rho -
       2176*acotrho*radixrho*mh4*mh6*rho + 320*mh8*rho - 192*lnmqmu*mh8*rho - 384*acot4mh2*radix4mh2*mh8*rho +
       1536*acot4mh2*radix4mh2*mh2*mh8*rho + 236*mh2*mh4*rho2 - 120*lnmqmu*mh2*mh4*rho2 - 240*acot4mh2*radix4mh2*mh2*mh4*rho2 -
       80*mh6*rho2 - 96*power_of<2>(acot4mh2)*mh6*rho2 + 96*power_of<2>(acotrho)*mh6*rho2 + 48*lnmqmu*mh6*rho2 -
       160*acot4mh2*radix4mh2*mh6*rho2 + 256*acotrho*radixrho*mh6*rho2 + 720*mh2*mh6*rho2 -
       128*acot4mh2*radix4mh2*mh2*mh6*rho2 - 1696*acotrho*radixrho*mh2*mh6*rho2 + 384*mh4*mh6*rho2 +
       1728*power_of<2>(acot4mh2)*mh4*mh6*rho2 - 1728*power_of<2>(acotrho)*mh4*mh6*rho2 - 768*acotrho*radixrho*mh4*mh6*rho2 - 720*mh8*rho2 +
       1824*acot4mh2*radix4mh2*mh8*rho2 - 21*mh4*rho*rho2 + 18*lnmqmu*mh4*rho*rho2 + 36*acot4mh2*radix4mh2*mh4*rho*rho2 +
       56*mh2*mh4*rho*rho2 + 72*acot4mh2*radix4mh2*mh2*mh4*rho*rho2 + 432*power_of<2>(acot4mh2)*mh2*mh6*rho*rho2 -
       432*power_of<2>(acotrho)*mh2*mh6*rho*rho2 - 384*power_of<2>(acot4mh2)*mh4*mh6*rho*rho2 + 384*power_of<2>(acotrho)*mh4*mh6*rho*rho2 -
       96*mh8*rho*rho2 + 28*mh4*rho3 - 20*lnmqmu*mh4*rho3 - 40*acot4mh2*radix4mh2*mh4*rho3 - 180*mh2*mh4*rho3 +
       64*acot4mh2*radix4mh2*mh2*mh4*rho3 - 32*mh6*rho3 - 64*acotrho*radixrho*mh6*rho3 + 192*acot4mh2*radix4mh2*mh8*rho3 -
       48*(1 + lnmqmu)*mh6*rho*(4*mh2 + rho)*ln4mh2 + 48*(1 + lnmqmu)*mh6*rho*(4*mh2 + rho)*lnrho))/
   (3.*mh6*power_of<4>(rho - 4.0 * mh2));

            double gb2 = (mh2*rho*(4352*mh6*mh8 - 1024*lnmqmu*mh6*mh8 - 2048*acotrho*radixrho*mh6*mh8 - 5120*mh12*rho +
       3072*lnmqmu*mh12*rho + 6144*acot4mh2*radix4mh2*mh12*rho - 24576*acot4mh2*radix4mh2*mh2*mh12*rho +
       24064*mh4*mh8*rho + 6144*power_of<2>(acot4mh2)*mh4*mh8*rho - 6144*power_of<2>(acotrho)*mh4*mh8*rho - 10240*lnmqmu*mh4*mh8*rho +
       10240*acot4mh2*radix4mh2*mh4*mh8*rho - 30720*acotrho*radixrho*mh4*mh8*rho - 8192*mh2*mh4*mh8*rho +
       8192*acot4mh2*radix4mh2*mh2*mh4*mh8*rho - 43520*mh6*mh8*rho + 35840*acotrho*radixrho*mh6*mh8*rho -
       9216*acot4mh2*radix4mh2*mh12*rho2 - 11520*mh4*mh6*rho2 + 6912*lnmqmu*mh4*mh6*rho2 +
       13824*acot4mh2*radix4mh2*mh4*mh6*rho2 - 55296*acot4mh2*radix4mh2*mh2*mh4*mh6*rho2 + 23040*mh12*rho2 +
       11520*mh2*mh8*rho2 + 4608*power_of<2>(acot4mh2)*mh2*mh8*rho2 - 4608*power_of<2>(acotrho)*mh2*mh8*rho2 - 6912*lnmqmu*mh2*mh8*rho2 +
       3072*acot4mh2*radix4mh2*mh2*mh8*rho2 - 16896*acotrho*radixrho*mh2*mh8*rho2 - 55296*mh4*mh8*rho2 +
       6144*acot4mh2*radix4mh2*mh4*mh8*rho2 + 70656*acotrho*radixrho*mh4*mh8*rho2 + 49920*acotrho*radixrho*mh2*mh4*mh8*rho2 -
       27840*mh6*mh8*rho2 - 55296*power_of<2>(acot4mh2)*mh6*mh8*rho2 + 55296*power_of<2>(acotrho)*mh6*mh8*rho2 -
       6912*acot4mh2*radix4mh2*mh4*mh6*rho*rho2 + 7680*mh12*rho*rho2 + 2016*mh8*rho*rho2 - 1728*lnmqmu*mh8*rho*rho2 -
       3456*acot4mh2*radix4mh2*mh8*rho*rho2 - 4480*mh2*mh8*rho*rho2 - 7680*mh4*mh8*rho*rho2 -
       41472*power_of<2>(acot4mh2)*mh4*mh8*rho*rho2 + 41472*power_of<2>(acotrho)*mh4*mh8*rho*rho2 +
       15360*acotrho*radixrho*mh4*mh8*rho*rho2 + 30720*power_of<2>(acot4mh2)*mh6*mh8*rho*rho2 -
       30720*power_of<2>(acotrho)*mh6*mh8*rho*rho2 + 168*mh2*mh4*rho4 - 144*lnmqmu*mh2*mh4*rho4 -
       288*acot4mh2*radix4mh2*mh2*mh4*rho4 + 1920*mh4*mh6*rho4 - 1120*mh8*rho4 -
       576*acot4mh2*radix4mh2*mh8*rho4 - 180*mh2*mh8*rho4 - 3456*power_of<2>(acot4mh2)*mh2*mh8*rho4 +
       3456*power_of<2>(acotrho)*mh2*mh8*rho4 - 8640*power_of<2>(acot4mh2)*mh6*mh8*rho4 +
       8640*power_of<2>(acotrho)*mh6*mh8*rho4 - 2880*mh2*mh6*rho3 + 1728*lnmqmu*mh2*mh6*rho3 +
       3456*acot4mh2*radix4mh2*mh2*mh6*rho3 + 17280*mh4*mh6*rho3 - 13824*acot4mh2*radix4mh2*mh4*mh6*rho3 - 320*mh8*rho3 +
       384*power_of<2>(acot4mh2)*mh8*rho3 - 384*power_of<2>(acotrho)*mh8*rho3 + 448*lnmqmu*mh8*rho3 + 1920*acot4mh2*radix4mh2*mh8*rho3 -
       1024*acotrho*radixrho*mh8*rho3 - 4736*mh2*mh8*rho3 + 3072*acot4mh2*radix4mh2*mh2*mh8*rho3 +
       14592*acotrho*radixrho*mh2*mh8*rho3 - 15360*acot4mh2*radix4mh2*mh4*mh8*rho3 + 8640*mh6*mh8*rho3 -
       17280*acotrho*radixrho*mh6*mh8*rho3 - 160*mh2*mh4*rho*rho3 + 160*lnmqmu*mh2*mh4*rho*rho3 +
       320*acot4mh2*radix4mh2*mh2*mh4*rho*rho3 - 2160*mh12*rho*rho3 + 360*mh8*rho*rho3 +
       640*acot4mh2*radix4mh2*mh8*rho*rho3 - 3840*acot4mh2*radix4mh2*mh2*mh8*rho*rho3 + 7680*power_of<2>(acot4mh2)*mh4*mh8*rho*rho3 -
       7680*power_of<2>(acotrho)*mh4*mh8*rho*rho3 - 25*mh6*rho4 - 12*lnmqmu*mh6*rho4 - 24*acot4mh2*radix4mh2*mh6*rho4 +
       1440*mh2*mh6*rho4 - 624*acot4mh2*radix4mh2*mh2*mh6*rho4 + 720*acot4mh2*radix4mh2*mh4*mh6*rho4 +
       4320*acot4mh2*radix4mh2*mh12*rho4 + 128*mh8*rho4 + 256*acotrho*radixrho*mh8*rho4 +
       192*(1 + lnmqmu)*mh8*rho*(16*mh4 + rho*(12*mh2 + rho))*ln4mh2 - 192*(1 + lnmqmu)*mh8*rho*(16*mh4 + 12*mh2*rho + rho2)*lnrho
       ))/(3.*mh8*power_of<5>(4.0 * mh2 - rho));

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J6
        complex<double>
        DileptonIntegralsBottom::j6(const double & a1, const double & a2) const
        {
            static const double pi = M_PI;
            const double acotrho = pi / 2.0 - atanrho, acot4mh2 = pi / 2.0 - atan4mh2;

            double asymp = (32*acotrho*radixrho*mh2*(3*rho*((-1 + rho)*rho + 2*mh2*(2 + rho)) - 4*mh4*(4 + 2*rho - 3*rho2)))/(9.*power_of<3>(4.0 * mh2 - rho)) +
   (32*power_of<2>(acotrho)*(3 - 2*mh2)*mh4*rho*rho2)/(3.*power_of<3>(4.0 * mh2 - rho)) -
   (8*acot4mh2*radix4mh2*mh2*rho*(4*mh2*mh4*rho2 + mh4*(-1 + 12*mh4)*rho2 + 12*mh6*rho2))/(9.*mh6*power_of<3>(4.0 * mh2 - rho)) +
   (4*lnmqmu*mh2*(-4*mh6*(16*mh4 + 3*rho*(-4*mh2 + rho)) + mh4*rho*rho2))/(9.*mh6*power_of<3>(4.0 * mh2 - rho)) -
   (2*mh2*(-8*mh6*(rho*(5*rho - 2*mh2*(7 + 15*rho)) + 4*mh4*(4 + 7*rho - 3*rho2)) +
        rho*(-24*mh8*rho2 + mh4*((3 - 60*mh2)*rho2 + 2*(rho2 + 7*mh2*rho2)))))/(9.*mh6*power_of<3>(4.0 * mh2 - rho)) +
   (32*power_of<2>(acot4mh2)*(-3 + 2*mh2)*mh4*mh6*rho3)/(3.*mh6*power_of<3>(4.0 * mh2 - rho));

            double gb1 = (32*power_of<2>(acotrho)*mh4*(mh4*(8 - 6*rho) + 6*mh2*(-2 + rho) - 3*rho)*rho*rho2)/power_of<4>(rho - 4.0 * mh2) -
   (32*power_of<2>(acot4mh2)*mh4*(mh4*(8 - 6*rho) + 6*mh2*(-2 + rho) - 3*rho)*rho3)/power_of<4>(rho - 4.0 * mh2) +
   (4*lnmqmu*mh2*rho*(4*mh8*rho2 + 3*mh2*mh4*(4*mh2 + rho)*rho2 - mh6*(12*mh2 + rho)*rho2 - mh2*mh4*(4*mh2 + 3*rho)*rho2 + mh6*rho3))/
    (3.*mh8*power_of<4>(rho - 4.0 * mh2)) + (8*acot4mh2*radix4mh2*mh2*rho*
      (16*mh4*mh6*rho2 - mh6*(rho3 + (-1 + 6*mh4 + 36*mh6)*rho3) + mh2*(36*mh8*rho*rho2 + 2*mh6*(-6*rho2 + 3*rho3)) +
        2*mh4*(24*mh6*rho2 + 24*mh8*rho2 + mh4*(4*rho2 - 3*rho*rho2 + 3*rho3))))/(3.*mh8*power_of<4>(rho - 4.0 * mh2)) +
   (mh2*(12*mh2*mh4*rho*(mh2*(4 - 6*rho) + 3*rho)*rho2 -
        4*mh6*rho*(144*mh6*rho2 - 108*mh6*rho*rho2 + 30*mh2*(-2 + 3*rho)*rho2 + 36*mh4*(10 + 3*rho)*rho2 - 5*rho3) - 11*mh6*rho*rho3 +
        6*mh2*mh4*(-7*rho + 14*mh2*(-2 + 3*rho) + mh4*(56 + 6*rho))*rho3 -
        16*mh8*(-(rho*(rho*(-18*mh2 - 5*rho + 90*mh2*rho) + 12*mh4*(2 + 9*rho + 9*rho2))) + 4*mh6*(4 + 12*rho - 27*rho2 + 27*rho3))))/
    (9.*mh8*power_of<4>(rho - 4.0 * mh2)) - (32*acotrho*radixrho*mh2*
      (24*mh2*mh4*rho2 + 36*mh4*rho*rho2 + (-1 + 18*mh2)*rho3 - 36*mh6*rho3 + rho4))/(3.*power_of<4>(rho - 4.0 * mh2));

            double gb2 = (-64*power_of<2>(acotrho)*mh4*rho*rho2*(3*(4*mh2*(-3 + rho) - rho)*rho + mh6*(32 - 60*rho + 30*rho2) + mh4*(64*rho - 6*(8 + 5*rho2))))/
    power_of<5>(4.0 * mh2 - rho) + (64*power_of<2>(acot4mh2)*mh4*(3*(4*mh2*(-3 + rho) - rho)*rho + mh6*(32 - 60*rho + 30*rho2) +
        mh4*(64*rho - 6*(8 + 5*rho2)))*rho3)/power_of<5>(4.0 * mh2 - rho) -
   (4*lnmqmu*mh2*rho*(4*mh10*(8*mh2 + rho)*rho2 + 6*mh2*mh6*(16*mh4 + rho*(12*mh2 + rho))*rho2 - mh8*(96*mh4 + rho*(32*mh2 + rho))*rho2 -
        4*mh8*(8*mh4 + rho*(16*mh2 + 3*rho))*rho2 + 10*mh2*mh6*(2*mh2 + rho)*rho3 - 3*mh8*rho4))/(3.*mh10*power_of<5>(4.0 * mh2 - rho)) +
   (32*acotrho*radixrho*mh2*(rho5 + 96*mh8*rho2 + 144*mh4*rho2*(rho + rho2) - 360*mh6*rho*rho3 - rho4 + 360*mh8*rho4 +
        mh2*(528*mh4*rho*rho2 - 8*rho3 - 480*mh6*rho3 + 44*rho4)))/(3.*power_of<5>(4.0 * mh2 - rho)) -
   (2*mh2*(-120*mh8*rho*rho2*(rho*(-3*rho + mh2*(-16 + 15*rho)) + mh4*(-8 + 30*rho + rho2)) +
        55*mh2*mh6*rho*(-5*rho + mh2*(-10 + 9*rho))*rho3 +
        30*mh2*mh6*(7*(12*mh2*(-1 + rho) - rho)*rho + mh6*(224 + 60*rho - 30*rho2) + 2*mh4*(-56 + 224*rho + 15*rho2))*rho3 -
        8*mh10*(-5*rho*(rho*(rho*(-76*mh2 - 5*rho + 180*mh2*rho) + 48*mh4*(-1 + 24*rho + 9*rho2)) +
              8*mh6*(4 + 48*rho + 243*rho2 - 135*rho3)) + 8*mh8*(8 + 60*rho - 390*rho2 + 1125*rho3 - 675*rho4)) + 78*mh8*rho*rho4 -
        10*mh8*rho*(1152*mh8*rho2 - 2160*mh8*rho*rho2 + 72*mh6*(40 + 32*rho - 15*rho2)*rho2 + 1080*mh8*rho4 +
           48*mh4*rho2*(-10 + 45*rho + 9*rho2) - 160*mh2*rho3 - 5*rho4 + 180*mh2*rho4)))/(45.*mh10*power_of<5>(4.0 * mh2 - rho)) -
   (8*acot4mh2*radix4mh2*mh2*rho*(184*mh12*rho3 +
        2*mh2*(2*mh4*(16*mh8*(5*rho2 + 12*rho*rho2) + 3*mh6*rho*(8*rho - 5*rho3)) + mh6*rho*(8*rho3 - 180*mh6*rho3) +
           mh8*(-32*rho*rho2 - 12*rho4 - 16*rho3 - rho4)) +
        mh8*(384*mh8*rho2 - 12*rho4 - rho4 + 3*(-1 + 20*mh6 + 120*mh8)*rho4 + 4*mh4*(-8*rho2 + 12*rho2*(-2 + 3*rho2) + 3*rho4)) +
        4*mh4*(6*mh4*mh6*(8*rho2 - 5*rho3) - 180*mh12*rho3 + mh6*(23*rho3 + 5*rho*rho3 + 3*rho4))))/
    (3.*mh10*power_of<5>(4.0 * mh2 - rho));

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // Massive case: charm quarks
        struct DileptonIntegralsCharm
        {
            double sh, sh2, sh3, sh4, lnsh;
            double mh, mh2, mh3, mh4, mh6, mh8, mh10, mh12, lnmh;
            double lnmqmu;
            double rho, rho2, rho3, rho4, rho5, rho6, rho7, lnrho, lnrhom1;
            double radixrho, radix4mh2;
            double lnradixrho, lndeltarho4mh2;
            double atanrho, atanh4mh2, atan4mh2rho, atannu, lnsigma;
            complex<double> aminus, aplus, lnam;
            double bminus, bplus, lnbm;
            complex<double> lntau;
            complex<double> dilogx4;
            complex<double> dilogx5;
            complex<double> diloginvx7;
            complex<double> diloginvx9;
            complex<double> dilogx13;
            double redilogx12;
            double redilog2ap;
            complex<double> trilogx4;
            complex<double> trilogx5;
            double retrilogx12;

            DileptonIntegralsCharm(const double & sh, const double & mh, const double & mB, const double & mu) :
                sh(sh), sh2(sh * sh), sh3(sh2 * sh), sh4(sh2 * sh2),
                lnsh(std::log(sh)),
                mh(mh), mh2(mh * mh), mh3(mh2 * mh), mh4(mh2 * mh2), mh6(mh4 * mh2), mh8(mh4 * mh4), mh10(mh8 * mh2), mh12(mh8 * mh4),
                lnmh(std::log(mh)), lnmqmu(2.0 * std::log(mh * mB / mu)),
                rho(4.0 * mh * mh / sh), rho2(rho * rho), rho3(rho2 * rho), rho4(rho2 * rho2), rho5(rho3 * rho2), rho6(rho3 * rho3), rho7(rho4 * rho3),
                lnrho(std::log(rho)), lnrhom1(std::log(rho - 1.0)),
                radixrho(std::sqrt(rho - 1.0)), radix4mh2(std::sqrt(1.0 - 4.0 * mh2)),
                lnradixrho(0.5 * lnrhom1), lndeltarho4mh2(std::log(rho - 4.0 * mh2)),
                atanrho(std::atan(radixrho)), atanh4mh2(std::atanh(radix4mh2)),
                atan4mh2rho(std::atan(radix4mh2 / radixrho)),
                atannu(std::atan((-2 + (2*rho)/(1 - radix4mh2))/(2.0*radixrho))),
                lnsigma(std::log(mh2 * rho / (rho - 4.0 * mh2))),
                aminus(0.5 * complex<double>(1.0, -radixrho)), aplus(1.0 - aminus), lnam(std::log(aminus)),
                bminus(0.5 * (1.0 + radix4mh2)), bplus(1.0 - bminus), lnbm(std::log(bminus)),
                lntau(std::log(bminus / mh)),
                dilogx4(dilog(power_of<2>(aminus / aplus))), dilogx5(dilog(-1.0 * bminus / bplus)),
                diloginvx7(dilog(aminus / bplus)), diloginvx9(dilog(aplus / bplus)),
                dilogx13(dilog((aplus * bminus) / (aminus * bplus))),
                redilogx12(real(dilog((aminus * bminus) / (aplus * bplus)))), redilog2ap(real(dilog(2.0 * aplus))),
                trilogx4(trilog(power_of<2>(aminus / aplus))), trilogx5(trilog(-1.0 * bminus / bplus)),
                retrilogx12(real(trilog((aminus * bminus) / (aplus * bplus))))
            {
            }

            complex<double> j1(const double & a1, const double & a2) const;
            complex<double> j2(const double & a1, const double & a2) const;
            complex<double> j3(const double & a1, const double & a2) const;
            complex<double> j4(const double & a1, const double & a2) const;
            complex<double> j5(const double & a1, const double & a2) const;
            complex<double> j6(const double & a1, const double & a2) const;
            complex<double> j7(const double & a1, const double & a2) const;
        };

        // J1
        complex<double>
        DileptonIntegralsCharm::j1(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi, pi3 = pi2 * pi;
            static const double ln2 = std::log(2.0);
            static const double zeta3 = 1.2020569031595942854;

            // Asymptotic part
            complex<double> asymp = (-80*power_of<3>(atanh4mh2)*mh2*rho)/(4*mh2 - rho) + (complex<double>(0,32)*power_of<3>(atanrho)*mh2*rho)/(4*mh2 - rho) -
   (12*zeta3*mh2*rho)/(4*mh2 - rho) - (complex<double>(0,24)*ln2*lnrhom1*mh2*pi*rho)/(4*mh2 - rho) +
   (complex<double>(0,24)*lnmh*lnsigma*mh2*pi*rho)/(4*mh2 - rho) +
   (complex<double>(0,24)*lnmh*mh2*pi*rho*(8*ln2*mh2 - 4*lnrho*mh2 - rho - 2*ln2*rho + lnrho*rho))/power_of<2>(-4*mh2 + rho) +
   lntau*((complex<double>(0,48)*lnmh*mh2*pi*rho)/(4*mh2 - rho) -
      (complex<double>(0,48)*mh2*pi*rho*(mh2*rho + ln2*(-4*mh2 + rho)))/power_of<2>(-4*mh2 + rho)) +
   atanh4mh2*((complex<double>(0,-96)*lnmh*mh2*pi*rho)/(4*mh2 - rho) - (complex<double>(0,48)*lntau*mh2*pi*rho)/(4*mh2 - rho) -
      (complex<double>(0,24)*mh2*rho*(4*ln2*pi*(4*mh2 - rho) + complex<double>(0,1)*(radix4mh2 + complex<double>(0,2)*mh2*pi)*rho + lnrho*pi*(-4*mh2 + rho)))/
       power_of<2>(-4*mh2 + rho)) + power_of<2>(atanrho)*((-48*atanh4mh2*mh2*rho)/(4*mh2 - rho) + (24*lnsigma*mh2*rho)/(4*mh2 - rho) +
      (48.0*lntau*mh2*rho)/(4*mh2 - rho) + (24.0*mh2*rho*(complex<double>(0,1)*(complex<double>(0,1) + pi)*rho + 2*mh2*(complex<double>(0,-2)*pi + rho)))/
       power_of<2>(-4*mh2 + rho)) + power_of<2>(atanh4mh2)*((24*lnsigma*mh2*rho)/(4*mh2 - rho) + (48.0*lntau*mh2*rho)/(4*mh2 - rho) +
      (24*mh2*rho*((-1.0 - complex<double>(0,3)*pi)*rho + 2*mh2*(complex<double>(0,6)*pi + rho)))/power_of<2>(-4*mh2 + rho)) +
   atanrho*((48*atanh4mh2*mh2*pi*rho)/(4*mh2 - rho) - (24*lnsigma*mh2*pi*rho)/(4*mh2 - rho) - (48.0*lntau*mh2*pi*rho)/(4*mh2 - rho) +
      (4*mh2*rho*(24*radixrho*mh2 + pi*(complex<double>(0,-4)*mh2*(pi - complex<double>(0,3)*rho) + (6.0 + complex<double>(0,1)*pi)*rho)))/power_of<2>(-4*mh2 + rho)) +
   (complex<double>(0,-4)*mh2*rho*(complex<double>(0,-6) + (complex<double>(0,-3) + 3*pi + pi3)*rho) + 3*rho2 +
      16*mh4*(3.0 + (3.0 - 3*radixrho*pi + complex<double>(0,1)*pi3)*rho - complex<double>(0,3)*(-1 + 2*ln2)*pi*rho2))/power_of<2>(-4*mh2 + rho) +
   (complex<double>(0,24)*mh2*pi*rho*diloginvx7)/(4*mh2 - rho) + (complex<double>(0,24)*mh2*pi*rho*diloginvx9)/(4*mh2 - rho) +
   ((complex<double>(0,-24)*atanrho*mh2*rho)/(4*mh2 - rho) + (complex<double>(0,12)*mh2*pi*rho)/(4*mh2 - rho))*dilogx4 +
   (complex<double>(0,24)*mh2*pi*rho*dilogx5)/(4*mh2 - rho) - (complex<double>(0,48)*mh2*pi*rho*redilog2ap)/(4*mh2 - rho) +
   ((-48*atanh4mh2*mh2*rho)/(4*mh2 - rho) - (complex<double>(0,24)*mh2*pi*rho)/(4*mh2 - rho))*redilogx12 +
   (24*mh2*rho*retrilogx12)/(4*mh2 - rho) - (12*mh2*rho*trilogx4)/(4*mh2 - rho);
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (-240*power_of<3>(atanh4mh2)*mh2*rho)/(4*mh2 - rho) + (complex<double>(0,96)*power_of<3>(atanrho)*mh2*rho)/(4*mh2 - rho) -
   (36*zeta3*mh2*rho)/(4*mh2 - rho) - (complex<double>(0,72)*ln2*lnrhom1*mh2*pi*rho)/(4*mh2 - rho) +
   (complex<double>(0,72)*lnmh*lnsigma*mh2*pi*rho)/(4*mh2 - rho) +
   (complex<double>(0,72)*lnmh*mh2*pi*rho*(2*ln2*power_of<2>(-4*mh2 + rho) - lnrho*power_of<2>(-4*mh2 + rho) + 2*rho*(-2*mh2 + rho)))/
    power_of<3>(4*mh2 - rho) + lntau*((complex<double>(0,144)*lnmh*mh2*pi*rho)/(4*mh2 - rho) +
      (complex<double>(0,144)*mh2*pi*rho*(ln2*power_of<2>(-4*mh2 + rho) + mh2*rho*(3*rho - mh2*(4 + 3*rho))))/power_of<3>(4*mh2 - rho)) +
   atanh4mh2*((complex<double>(0,-288)*lnmh*mh2*pi*rho)/(4*mh2 - rho) - (complex<double>(0,144)*lntau*mh2*pi*rho)/(4*mh2 - rho) -
      (complex<double>(0,36)*mh2*rho*(8*ln2*pi*power_of<2>(-4*mh2 + rho) - 2*lnrho*pi*power_of<2>(-4*mh2 + rho) +
           rho*(-4*mh2*pi*(-3*rho + mh2*(4 + 3*rho)) + complex<double>(0,1)*radix4mh2*(-5*rho + mh2*(8 + 6*rho)))))/power_of<3>(4*mh2 - rho)) +
   power_of<2>(atanrho)*((-144*atanh4mh2*mh2*rho)/(4*mh2 - rho) + (72*lnsigma*mh2*rho)/(4*mh2 - rho) + (144.0*lntau*mh2*rho)/(4*mh2 - rho) +
      (72*mh2*rho*(2*mh2*(-2.0 + complex<double>(0,4)*pi - 3*rho)*rho + 2*mh4*(complex<double>(0,-8)*pi + rho*(4 + 3*rho)) + (2.0 - complex<double>(0,1)*pi)*rho2))/
       power_of<3>(4*mh2 - rho)) + power_of<2>(atanh4mh2)*((72*lnsigma*mh2*rho)/(4*mh2 - rho) + (144.0*lntau*mh2*rho)/(4*mh2 - rho) +
      (72*mh2*rho*(-2*mh2*rho*(2.0 + complex<double>(0,12)*pi + 3*rho) + mh4*(complex<double>(0,48)*pi + 2*rho*(4 + 3*rho)) + (2.0 + complex<double>(0,3)*pi)*rho2))/
       power_of<3>(4*mh2 - rho)) + atanrho*((144*atanh4mh2*mh2*pi*rho)/(4*mh2 - rho) - (72*lnsigma*mh2*pi*rho)/(4*mh2 - rho) -
      (144.0*lntau*mh2*pi*rho)/(4*mh2 - rho) - (12*mh2*rho*
         (-72*radixrho*(-(mh2*rho) + mh4*(2 + rho)) + pi*
            (-4*mh2*rho*(6.0 + complex<double>(0,2)*pi + 9*rho) + 4*mh4*(complex<double>(0,4)*pi + 3*rho*(4 + 3*rho)) + (12.0 + complex<double>(0,1)*pi)*rho2)))/
       power_of<3>(4*mh2 - rho)) + (3.0*(mh2*(12.0 + (33.0 + complex<double>(0,30)*pi + complex<double>(0,4)*pi3)*rho)*rho2 +
        4*mh4*rho*(-12.0 + 4.0*(-12.0 - complex<double>(0,3)*pi + 9*radixrho*pi - complex<double>(0,2)*pi3)*rho +
           complex<double>(0,9)*(complex<double>(0,1) - 4*pi + 8*ln2*pi)*rho2) - rho3 +
        16*mh6*(4.0 + (15.0 - 18*radixrho*pi + complex<double>(0,4)*pi3)*rho +
           3.0*(3.0 + (complex<double>(0,4) - complex<double>(0,8)*ln2 - 3*radixrho)*pi)*rho2 - complex<double>(0,6)*(-1 + 4*ln2)*pi*rho3)))/power_of<3>(4*mh2 - rho) +
   (complex<double>(0,72)*mh2*pi*rho*diloginvx7)/(4*mh2 - rho) + (complex<double>(0,72)*mh2*pi*rho*diloginvx9)/(4*mh2 - rho) +
   ((complex<double>(0,-72)*atanrho*mh2*rho)/(4*mh2 - rho) + (complex<double>(0,36)*mh2*pi*rho)/(4*mh2 - rho))*dilogx4 +
   (complex<double>(0,72)*mh2*pi*rho*dilogx5)/(4*mh2 - rho) - (complex<double>(0,144)*mh2*pi*rho*redilog2ap)/(4*mh2 - rho) +
   ((-144*atanh4mh2*mh2*rho)/(4*mh2 - rho) - (complex<double>(0,72)*mh2*pi*rho)/(4*mh2 - rho))*redilogx12 +
   (72*mh2*rho*retrilogx12)/(4*mh2 - rho) - (36*mh2*rho*trilogx4)/(4*mh2 - rho);
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (complex<double>(0,192)*power_of<3>(atanrho)*mh2*rho)/(4*mh2 - rho) - (72*zeta3*mh2*rho)/(4*mh2 - rho) +
   (480*power_of<3>(atanh4mh2)*mh2*rho)/(-4*mh2 + rho) + (complex<double>(0,144)*ln2*lnrhom1*mh2*pi*rho)/(-4*mh2 + rho) -
   (complex<double>(0,144)*lnmh*lnsigma*mh2*pi*rho)/(-4*mh2 + rho) +
   (complex<double>(0,48)*lnmh*mh2*pi*rho*(6*ln2*power_of<3>(4*mh2 - rho) - 3*lnrho*power_of<3>(4*mh2 - rho) - 8*rho*(6*mh4 - 3*mh2*rho + rho2)))/
    power_of<4>(-4*mh2 + rho) + lntau*((complex<double>(0,-288)*lnmh*mh2*pi*rho)/(-4*mh2 + rho) +
      (complex<double>(0,96)*mh2*pi*rho*(3*ln2*power_of<3>(4*mh2 - rho) + mh2*rho*(3*mh2*rho*(8 + 15*rho) - 18*rho2 - 2*mh4*(24 + 25*rho2))))/
       power_of<4>(-4*mh2 + rho)) + atanh4mh2*((complex<double>(0,576)*lnmh*mh2*pi*rho)/(-4*mh2 + rho) +
      (complex<double>(0,288)*lntau*mh2*pi*rho)/(-4*mh2 + rho) +
      (8*mh2*rho*(complex<double>(0,-72)*ln2*pi*power_of<3>(4*mh2 - rho) + complex<double>(0,18)*lnrho*pi*power_of<3>(4*mh2 - rho) +
           rho*(-4*mh2*rho*(36*radix4mh2 + 55*radix4mh2*rho - complex<double>(0,54)*pi*rho) + 73*radix4mh2*rho2 +
              complex<double>(0,24)*mh6*pi*(24 + 25*rho2) + 12*mh4*(complex<double>(0,-3)*pi*rho*(8 + 15*rho) + radix4mh2*(24 + 25*rho2)))))/
       power_of<4>(-4*mh2 + rho)) + atanrho*((288*atanh4mh2*mh2*pi*rho)/(4*mh2 - rho) + (144*lnsigma*mh2*pi*rho)/(-4*mh2 + rho) +
      (288.0*lntau*mh2*pi*rho)/(-4*mh2 + rho) - (8*mh2*rho*
         (36*mh2*(complex<double>(0,1)*radixrho*power_of<2>(pi) - 12*(-1 + rho) + 2*radixrho*pi*(2 + 3*rho))*rho2 -
           36*mh4*rho*(complex<double>(0,4)*radixrho*power_of<2>(pi) - 6*(-6 + rho + 5*rho2) + radixrho*pi*(8 + 8*rho + 15*rho2)) +
           8*mh6*(224.0 + complex<double>(0,24)*radixrho*power_of<2>(pi) - 124*rho + 50*rho2 + 3*radixrho*pi*rho*(24 + 25*rho2) - 150*rho3) -
           complex<double>(0,3)*radixrho*pi*(complex<double>(0,-16) + pi)*rho3))/(radixrho*power_of<4>(-4*mh2 + rho))) +
   power_of<2>(atanrho)*((288*atanh4mh2*mh2*rho)/(-4*mh2 + rho) - (144*lnsigma*mh2*rho)/(-4*mh2 + rho) -
      (288.0*lntau*mh2*rho)/(-4*mh2 + rho) + (48*mh2*rho*(12*mh2*(2.0 - complex<double>(0,3)*pi + 3*rho)*rho2 -
           6*mh4*rho*(8.0 - complex<double>(0,24)*pi + 8*rho + 15*rho2) + (-8.0 + complex<double>(0,3)*pi)*rho3 + 4*mh6*(complex<double>(0,-48)*pi + 24*rho + 25*rho3)))/
       power_of<4>(-4*mh2 + rho)) + power_of<2>(atanh4mh2)*((-144*lnsigma*mh2*rho)/(-4*mh2 + rho) - (288.0*lntau*mh2*rho)/(-4*mh2 + rho) +
      (48*mh2*rho*(12*mh2*(2.0 + complex<double>(0,9)*pi + 3*rho)*rho2 - 6*mh4*rho*(8.0 + complex<double>(0,72)*pi + 8*rho + 15*rho2) +
           (-8.0 - complex<double>(0,9)*pi)*rho3 + 4*mh6*(complex<double>(0,144)*pi + 24*rho + 25*rho3)))/power_of<4>(-4*mh2 + rho)) +
   (9*radixrho*rho4 + 12*mh4*rho2*(complex<double>(0,72)*radixrho*pi3*rho +
         144*pi*rho*(3.0 + complex<double>(0,1)*radixrho + (-3.0 + complex<double>(0,3)*radixrho - complex<double>(0,6)*ln2*radixrho)*rho) +
         radixrho*(72 + 576*rho + 245*rho2)) - complex<double>(0,2)*radixrho*mh2*(complex<double>(0,-72) + (complex<double>(0,-533) + 438*pi + 36*pi3)*rho)*
       rho3 + 144*mh6*rho*(complex<double>(0,-24)*radixrho*pi3*rho +
         6*pi*rho*(-18.0 - complex<double>(0,4)*radixrho + (3.0 - complex<double>(0,8)*radixrho + complex<double>(0,16)*ln2*radixrho)*rho +
            5.0*(3.0 - complex<double>(0,2)*radixrho + complex<double>(0,8)*ln2*radixrho)*rho2) - radixrho*(16 + 102*rho + 90*rho2 + 25*rho3)) +
      64*mh8*(complex<double>(0,72)*radixrho*pi3*rho + radixrho*(36 + 256*rho + 75*rho2 + 225*rho3) +
         3.0*pi*rho*(112.0 + (-62.0 + complex<double>(0,72)*radixrho - complex<double>(0,144)*ln2*radixrho)*rho + 25*rho2 -
            complex<double>(0,5)*(complex<double>(0,-15) - 8*radixrho + 48*ln2*radixrho)*rho3)))/(3.*radixrho*power_of<4>(-4*mh2 + rho)) +
   (complex<double>(0,144)*mh2*pi*rho*diloginvx7)/(4*mh2 - rho) + (complex<double>(0,144)*mh2*pi*rho*diloginvx9)/(4*mh2 - rho) +
   ((complex<double>(0,-144)*atanrho*mh2*rho)/(4*mh2 - rho) + (complex<double>(0,72)*mh2*pi*rho)/(4*mh2 - rho))*dilogx4 -
   (complex<double>(0,144)*mh2*pi*rho*dilogx5)/(-4*mh2 + rho) - (complex<double>(0,288)*mh2*pi*rho*redilog2ap)/(4*mh2 - rho) +
   ((-288*atanh4mh2*mh2*rho)/(4*mh2 - rho) - (complex<double>(0,144)*mh2*pi*rho)/(4*mh2 - rho))*redilogx12 +
   (144*mh2*rho*retrilogx12)/(4*mh2 - rho) - (72*mh2*rho*trilogx4)/(4*mh2 - rho);
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J2
        complex<double>
        DileptonIntegralsCharm::j2(const double & a1, const double & a2) const
        {
            static const double pi = M_PI;
            static const double ln2 = std::log(2.0);

            // Asymptotic part
            complex<double> asymp = (complex<double>(0,-12)*atannu*pi*(-1 + rho))/radixrho - (6*lnrhom1*pi*(-1 + rho))/radixrho - (6*lnsigma*pi*(-1 + rho))/radixrho -
   (12.0*lntau*pi*(-1 + rho))/radixrho - (12*power_of<2>(atanh4mh2)*(2*mh2*(-2 + rho) + rho))/(-4*mh2 + rho) +
   (12*power_of<2>(atanrho)*(2*mh2*(radixrho*(-2 + rho) - complex<double>(0,4)*(-1 + rho)) + (radixrho + complex<double>(0,2)*(-1 + rho))*rho))/
    (radixrho*(4*mh2 - rho)) + (6*radixrho*(-4*mh2 + rho + complex<double>(0,1)*radix4mh2*pi*rho) +
      pi*(-1 + rho)*(48*lnrho*mh2 - complex<double>(0,4)*mh2*pi + 6*rho - 12*lnrho*rho + complex<double>(0,1)*pi*rho + 24*ln2*(-4*mh2 + rho)))/
    (radixrho*(4*mh2 - rho)) + (12*atanh4mh2*(2*mh2*pi*(4.0 + complex<double>(0,1)*radixrho*(-2 + rho) - 4*rho) -
        2*atan4mh2rho*(4*mh2 - rho)*(-1 + rho) + rho*(radix4mh2*radixrho + pi*(-2.0 + complex<double>(0,1)*radixrho + 2*rho))))/
    (radixrho*(-4*mh2 + rho)) + atanrho*((-48*atanh4mh2*(-1 + rho))/radixrho + (12*lnrhom1*(-1 + rho))/radixrho +
      (12*lnsigma*(-1 + rho))/radixrho + (24.0*lntau*(-1 + rho))/radixrho +
      (12.0*(complex<double>(0,8)*mh2*pi - 4*radixrho*mh2*pi - 4*ln2*(4*mh2 - rho)*(-1 + rho) + 2*lnrho*(4*mh2 - rho)*(-1 + rho) - rho -
           complex<double>(0,2)*pi*rho + radixrho*pi*rho - complex<double>(0,8)*mh2*pi*rho + 2*radixrho*mh2*pi*rho + rho2 + complex<double>(0,2)*pi*rho2))/
       (radixrho*(-4*mh2 + rho))) + (complex<double>(0,12)*(-1 + rho)*dilogx13)/radixrho + (complex<double>(0,6)*(-1 + rho)*dilogx4)/radixrho -
   (complex<double>(0,12)*(-1 + rho)*redilogx12)/radixrho;
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (complex<double>(0,-36)*atannu*pi*(-1 + rho))/radixrho - (18*lnrhom1*pi*(-1 + rho))/radixrho - (18*lnsigma*pi*(-1 + rho))/radixrho -
   (36.0*lntau*pi*(-1 + rho))/radixrho + (36*power_of<2>(atanh4mh2)*(2*mh2*(4 - 3*rho)*rho - rho2 + 4*mh4*(-4 + 2*rho + rho2)))/
    power_of<2>(-4*mh2 + rho) - (3.0*(pi*(-1 + rho)*(complex<double>(0,16)*mh4*pi - 24*mh2*rho + 48*mh4*rho - complex<double>(0,8)*mh2*pi*rho +
           24*ln2*power_of<2>(-4*mh2 + rho) - 12*lnrho*power_of<2>(-4*mh2 + rho) + 12*rho2 + complex<double>(0,1)*pi*rho2) -
        3.0*radixrho*(16*mh4*(-3 + rho) + 4*mh2*(8.0 - complex<double>(0,1)*radix4mh2*pi*(-2 + rho) - rho)*rho +
           (-5.0 - complex<double>(0,4)*radix4mh2*pi)*rho2)))/(radixrho*power_of<2>(-4*mh2 + rho)) +
   (36*atanh4mh2*(2*atan4mh2rho*(-1 + rho)*power_of<2>(-4*mh2 + rho) +
        2*mh2*rho*(radix4mh2*radixrho*(-2 + rho) + pi*(8 - 8*rho + complex<double>(0,1)*radixrho*(-4 + 3*rho))) +
        (2*radix4mh2*radixrho + pi*(-2.0 + complex<double>(0,1)*radixrho + 2*rho))*rho2 +
        4*mh4*pi*(8.0*(-1 + rho) - complex<double>(0,1)*radixrho*(-4 + 2*rho + rho2))))/(radixrho*power_of<2>(-4*mh2 + rho)) +
   (36*power_of<2>(atanrho)*(-2*mh2*rho*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-4 + 3*rho)) - (radixrho + complex<double>(0,2)*(-1 + rho))*rho2 +
        4*mh4*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-4 + 2*rho + rho2))))/(radixrho*power_of<2>(-4*mh2 + rho)) +
   atanrho*((-144*atanh4mh2*(-1 + rho))/radixrho + (36*lnrhom1*(-1 + rho))/radixrho + (36*lnsigma*(-1 + rho))/radixrho +
      (72.0*lntau*(-1 + rho))/radixrho + (36.0*(complex<double>(0,-32)*mh4*pi + 16*radixrho*mh4*pi + 4*mh2*rho - 8*mh4*rho +
           complex<double>(0,16)*mh2*pi*rho - 8*radixrho*mh2*pi*rho + complex<double>(0,32)*mh4*pi*rho - 8*radixrho*mh4*pi*rho +
           4*ln2*(-1 + rho)*power_of<2>(-4*mh2 + rho) - 2*lnrho*(-1 + rho)*power_of<2>(-4*mh2 + rho) - 2*rho2 - 4*mh2*rho2 + 8*mh4*rho2 -
           complex<double>(0,2)*pi*rho2 + radixrho*pi*rho2 - complex<double>(0,16)*mh2*pi*rho2 + 6*radixrho*mh2*pi*rho2 - 4*radixrho*mh4*pi*rho2 +
           2*rho3 + complex<double>(0,2)*pi*rho3))/(radixrho*power_of<2>(-4*mh2 + rho))) + (complex<double>(0,36)*(-1 + rho)*dilogx13)/radixrho +
   (complex<double>(0,18)*(-1 + rho)*dilogx4)/radixrho - (complex<double>(0,36)*(-1 + rho)*redilogx12)/radixrho;
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (-72*atannu*pi*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) +
   (complex<double>(0,36)*lnrhom1*pi*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) +
   (complex<double>(0,36)*lnsigma*pi*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) +
   (complex<double>(0,72)*lntau*pi*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) +
   (complex<double>(0,72)*power_of<2>(atanh4mh2)*rho2*(12*mh2*(-1 + rho)*rho2 - 4*mh4*rho*(-12 + 4*rho + 5*rho2) + rho3 + 4*mh6*(-16 + 8*rho + 5*rho3)))/
    (power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))*power_of<3>(-4*mh2 + rho)) +
   (complex<double>(0,2)*rho2*(-3*pi*(-1 + rho)*(complex<double>(0,64)*mh6*pi + 24*ln2*power_of<3>(4*mh2 - rho) - 12*lnrho*power_of<3>(4*mh2 - rho) - 96*mh4*rho +
           160*mh6*rho - complex<double>(0,48)*mh4*pi*rho + 48*mh2*rho2 - 240*mh4*rho2 + 240*mh6*rho2 + complex<double>(0,12)*mh2*pi*rho2 - 16*rho3 -
           complex<double>(0,1)*pi*rho3) + radixrho*(3*mh2*(-192 + 55*rho + complex<double>(0,2)*radix4mh2*pi*(-24 + 25*rho))*rho2 +
           16*mh6*(-112 + 15*rho + 45*rho2) + 36*mh4*rho*(44 - 20*rho - 5*rho2 - complex<double>(0,1)*radix4mh2*pi*(-8 + 5*rho2)) +
           (73.0 + complex<double>(0,48)*radix4mh2*pi)*rho3)))/
    (radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))*power_of<3>(-4*mh2 + rho)) -
   (complex<double>(0,24)*atanh4mh2*rho2*(-6*atan4mh2rho*power_of<3>(4*mh2 - rho)*(-1 + rho) +
        mh2*(complex<double>(0,36)*(complex<double>(0,2) + radixrho)*pi*(-1 + rho) + radix4mh2*radixrho*(-24 + 25*rho))*rho2 -
        6*mh4*rho*(radix4mh2*radixrho*(-8 + 5*rho2) + 2*pi*(-24*(-1 + rho) + complex<double>(0,1)*radixrho*(-12 + 4*rho + 5*rho2))) +
        (8*radix4mh2*radixrho + 3*pi*(-2.0 + complex<double>(0,1)*radixrho + 2*rho))*rho3 +
        12*mh6*pi*(-32*(-1 + rho) + complex<double>(0,1)*radixrho*(-16 + 8*rho + 5*rho3))))/
    (radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))*power_of<3>(-4*mh2 + rho)) +
   (complex<double>(0,72)*power_of<2>(atanrho)*rho2*(12.0*(complex<double>(0,-2) + radixrho)*mh2*(-1 + rho)*rho2 -
        4*mh4*rho*(complex<double>(0,-24)*(-1 + rho) + radixrho*(-12 + 4*rho + 5*rho2)) + (radixrho + complex<double>(0,2)*(-1 + rho))*rho3 +
        4*mh6*(complex<double>(0,-32)*(-1 + rho) + radixrho*(-16 + 8*rho + 5*rho3))))/
    (radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))*power_of<3>(-4*mh2 + rho)) +
   atanrho*((complex<double>(0,288)*atanh4mh2*(-1 + rho)*rho2)/
       (radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) -
      (complex<double>(0,72)*lnrhom1*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) -
      (complex<double>(0,72)*lnsigma*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) -
      (complex<double>(0,144)*lntau*(-1 + rho)*rho2)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) +
      (complex<double>(0,24)*rho2*(complex<double>(0,-384)*mh6*pi + 192*radixrho*mh6*pi + 12*ln2*power_of<3>(4*mh2 - rho)*(-1 + rho) -
           6*lnrho*power_of<3>(4*mh2 - rho)*(-1 + rho) + 48*mh4*rho - 80*mh6*rho + complex<double>(0,288)*mh4*pi*rho - 144*radixrho*mh4*pi*rho +
           complex<double>(0,384)*mh6*pi*rho - 96*radixrho*mh6*pi*rho - 24*mh2*rho2 + 72*mh4*rho2 - 40*mh6*rho2 - complex<double>(0,72)*mh2*pi*rho2 +
           36*radixrho*mh2*pi*rho2 - complex<double>(0,288)*mh4*pi*rho2 + 48*radixrho*mh4*pi*rho2 + 8*rho3 + 24*mh2*rho3 - 120*mh4*rho3 +
           120*mh6*rho3 + complex<double>(0,6)*pi*rho3 - 3*radixrho*pi*rho3 + complex<double>(0,72)*mh2*pi*rho3 - 36*radixrho*mh2*pi*rho3 +
           60*radixrho*mh4*pi*rho3 - 60*radixrho*mh6*pi*rho3 - 8*rho4 - complex<double>(0,6)*pi*rho4))/
       (radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))*power_of<3>(-4*mh2 + rho))) +
   (72*(-1 + rho)*rho2*dilogx13)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) +
   (36*(-1 + rho)*rho2*dilogx4)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho))) -
   (72*(-1 + rho)*rho2*redilogx12)/(radixrho*power_of<2>(complex<double>(0,-1) + radixrho)*(2*radixrho - complex<double>(0,1)*(-2 + rho)));
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J3
        complex<double>
        DileptonIntegralsCharm::j3(const double & a1, const double & a2) const
        {
            static const double pi = M_PI;
            static const double ln2 = std::log(2.0);

            // Asymptotic part
            complex<double> asymp = (complex<double>(0,-48)*atannu*mh2*pi*(-1 + rho))/(radixrho*rho) - (24*lnrhom1*mh2*pi*(-1 + rho))/(radixrho*rho) -
   (24*lnsigma*mh2*pi*(-1 + rho))/(radixrho*rho) - (48.0*lntau*mh2*pi*(-1 + rho))/(radixrho*rho) -
   (24*power_of<2>(atanh4mh2)*mh2*(-((-2 + rho)*rho) + mh2*(-8 + 4*rho + rho2)))/(rho*(-4*mh2 + rho)) +
   (2*pi*(-1 + rho)*(-192*lnrho*mh4 + complex<double>(0,16)*mh4*pi - 24*mh2*rho + 48*lnrho*mh2*rho + 24*mh4*rho - complex<double>(0,4)*mh2*pi*rho +
         96*ln2*(4*mh4 - mh2*rho) + 3*rho2) + 3*radixrho*
       (-16*mh4*(-5 + rho) + 4*mh2*rho*(-8.0 + complex<double>(0,1)*radix4mh2*pi*(-4 + rho) + rho) + (3.0 + complex<double>(0,2)*radix4mh2*pi)*rho2))/
    (2.*radixrho*rho*(-4*mh2 + rho)) - (6*atanh4mh2*(16*atan4mh2rho*mh2*(4*mh2 - rho)*(-1 + rho) +
        2*mh2*(2*pi*(4.0 + complex<double>(0,1)*radixrho*(-2 + rho) - 4*rho) + radix4mh2*radixrho*(-4 + rho))*rho +
        radix4mh2*radixrho*rho2 + 4*mh4*pi*(16*(-1 + rho) - complex<double>(0,1)*radixrho*(-8 + 4*rho + rho2))))/
    (radixrho*rho*(-4*mh2 + rho)) - (24*power_of<2>(atanrho)*mh2*
      ((-(radixrho*(-2 + rho)) + complex<double>(0,4)*(-1 + rho))*rho + mh2*(complex<double>(0,-16)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2))))/
    (radixrho*rho*(-4*mh2 + rho)) + atanrho*((-192*atanh4mh2*mh2*(-1 + rho))/(radixrho*rho) +
      (48*lnrhom1*mh2*(-1 + rho))/(radixrho*rho) + (48*lnsigma*mh2*(-1 + rho))/(radixrho*rho) +
      (96.0*lntau*mh2*(-1 + rho))/(radixrho*rho) - (6.0*(complex<double>(0,-64)*mh4*pi + 32*radixrho*mh4*pi +
           32*ln2*mh2*(4*mh2 - rho)*(-1 + rho) - 16*lnrho*mh2*(4*mh2 - rho)*(-1 + rho) + 8*mh2*rho - 8*mh4*rho +
           complex<double>(0,16)*mh2*pi*rho - 8*radixrho*mh2*pi*rho + complex<double>(0,64)*mh4*pi*rho - 16*radixrho*mh4*pi*rho - rho2 - 8*mh2*rho2 +
           8*mh4*rho2 - complex<double>(0,16)*mh2*pi*rho2 + 4*radixrho*mh2*pi*rho2 - 4*radixrho*mh4*pi*rho2 + rho3))/
       (radixrho*rho*(-4*mh2 + rho))) + (complex<double>(0,48)*mh2*(-1 + rho)*dilogx13)/(radixrho*rho) +
   (complex<double>(0,24)*mh2*(-1 + rho)*dilogx4)/(radixrho*rho) - (complex<double>(0,48)*mh2*(-1 + rho)*redilogx12)/(radixrho*rho);
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (complex<double>(0,-144)*atannu*mh2*pi*(-1 + rho))/(radixrho*rho) - (72*lnrhom1*mh2*pi*(-1 + rho))/(radixrho*rho) -
   (72*lnsigma*mh2*pi*(-1 + rho))/(radixrho*rho) - (144.0*lntau*mh2*pi*(-1 + rho))/(radixrho*rho) +
   (72*power_of<2>(atanh4mh2)*mh2*(mh2*rho*(16 - 8*rho - 3*rho2) + (-2 + rho)*rho2 + 4*mh4*(-8 + 4*rho + rho2 + rho3)))/
    (rho*power_of<2>(-4*mh2 + rho)) + (-6*pi*(-1 + rho)*(complex<double>(0,64)*mh6*pi - 96*mh4*rho + 160*mh6*rho - complex<double>(0,32)*mh4*pi*rho +
         96*ln2*power_of<2>(-4*mh3 + mh*rho) - 48*lnrho*power_of<2>(-4*mh3 + mh*rho) + 36*mh2*rho2 - 72*mh4*rho2 + 96*mh6*rho2 +
         complex<double>(0,4)*mh2*pi*rho2 - rho3) + radixrho*(12*mh2*(-33 + 8*rho + complex<double>(0,1)*radix4mh2*pi*(-18 + 7*rho))*rho2 +
         64*mh6*(-53 + 12*rho + 9*rho2) + 144*mh4*rho*(15 - 4*rho - rho2 - complex<double>(0,1)*radix4mh2*pi*(-4 + rho + rho2)) +
         (17.0 + complex<double>(0,6)*radix4mh2*pi)*rho3))/(2.*radixrho*rho*power_of<2>(-4*mh2 + rho)) +
   (72*power_of<2>(atanrho)*mh2*(mh2*rho*(complex<double>(0,32)*(-1 + rho) + radixrho*(16 - 8*rho - 3*rho2)) +
        (radixrho*(-2 + rho) - complex<double>(0,4)*(-1 + rho))*rho2 + 4*mh4*(complex<double>(0,-16)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2 + rho3))))/
    (radixrho*rho*power_of<2>(-4*mh2 + rho)) + (6*atanh4mh2*
      (48*atan4mh2rho*mh2*(-1 + rho)*power_of<2>(-4*mh2 + rho) +
        2*mh2*(6*pi*(complex<double>(0,-1)*radixrho*(-2 + rho) + 4*(-1 + rho)) + radix4mh2*radixrho*(18 - 7*rho))*rho2 +
        12*mh4*rho*(2*radix4mh2*radixrho*(-4 + rho + rho2) + pi*(-32*(-1 + rho) + complex<double>(0,1)*radixrho*(-16 + 8*rho + 3*rho2))) -
        radix4mh2*radixrho*rho3 - complex<double>(0,48)*mh6*pi*(complex<double>(0,16)*(-1 + rho) + radixrho*(-8 + 4*rho + rho2 + rho3))))/
    (radixrho*rho*power_of<2>(-4*mh2 + rho)) + atanrho*((-576*atanh4mh2*mh2*(-1 + rho))/(radixrho*rho) +
      (144*lnrhom1*mh2*(-1 + rho))/(radixrho*rho) + (144*lnsigma*mh2*(-1 + rho))/(radixrho*rho) +
      (288.0*lntau*mh2*(-1 + rho))/(radixrho*rho) + (6.0*
         (complex<double>(0,-768)*mh6*pi + 384*radixrho*mh6*pi + 96*mh4*rho - 160*mh6*rho + complex<double>(0,384)*mh4*pi*rho -
           192*radixrho*mh4*pi*rho + complex<double>(0,768)*mh6*pi*rho - 192*radixrho*mh6*pi*rho +
           96*ln2*mh2*(-1 + rho)*power_of<2>(-4*mh2 + rho) - 48*lnrho*mh2*(-1 + rho)*power_of<2>(-4*mh2 + rho) - 36*mh2*rho2 - 24*mh4*rho2 +
           64*mh6*rho2 - complex<double>(0,48)*mh2*pi*rho2 + 24*radixrho*mh2*pi*rho2 - complex<double>(0,384)*mh4*pi*rho2 + 96*radixrho*mh4*pi*rho2 -
           48*radixrho*mh6*pi*rho2 + rho3 + 36*mh2*rho3 - 72*mh4*rho3 + 96*mh6*rho3 + complex<double>(0,48)*mh2*pi*rho3 -
           12*radixrho*mh2*pi*rho3 + 36*radixrho*mh4*pi*rho3 - 48*radixrho*mh6*pi*rho3 - rho4))/(radixrho*rho*power_of<2>(-4*mh2 + rho)))
     + (complex<double>(0,144)*mh2*(-1 + rho)*dilogx13)/(radixrho*rho) + (complex<double>(0,72)*mh2*(-1 + rho)*dilogx4)/(radixrho*rho) -
   (complex<double>(0,144)*mh2*(-1 + rho)*redilogx12)/(radixrho*rho);
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (complex<double>(0,-288)*atannu*mh2*pi*(-1 + rho))/(radixrho*rho) - (144*lnrhom1*mh2*pi*(-1 + rho))/(radixrho*rho) -
   (144*lnsigma*mh2*pi*(-1 + rho))/(radixrho*rho) - (288.0*lntau*mh2*pi*(-1 + rho))/(radixrho*rho) +
   (144*power_of<2>(atanh4mh2)*mh2*(6*mh2*rho2*(-4 + 2*rho + rho2) - (-2 + rho)*rho3 - 4*mh4*rho*(-24 + 12*rho + 2*rho2 + 5*rho3) +
        mh6*(-128 + 64*rho + 16*rho2 + 25*rho4)))/(power_of<3>(4*mh2 - rho)*rho) +
   atanrho*((-1152*atanh4mh2*mh2*(-1 + rho))/(radixrho*rho) + (288*lnrhom1*mh2*(-1 + rho))/(radixrho*rho) +
      (288*lnsigma*mh2*(-1 + rho))/(radixrho*rho) + (576.0*lntau*mh2*(-1 + rho))/(radixrho*rho) +
      (6.0*(complex<double>(0,-6144)*mh8*pi + 3072*radixrho*mh8*pi + 192*ln2*mh2*power_of<3>(4*mh2 - rho)*(-1 + rho) -
           96*lnrho*mh2*power_of<3>(4*mh2 - rho)*(-1 + rho) + 768*mh6*rho - 1408*mh8*rho + complex<double>(0,4608)*mh6*pi*rho -
           2304*radixrho*mh6*pi*rho + complex<double>(0,6144)*mh8*pi*rho - 1536*radixrho*mh8*pi*rho + rho5 - 480*mh4*rho2 +
           256*mh6*rho2 + 608*mh8*rho2 - complex<double>(0,1152)*mh4*pi*rho2 + 576*radixrho*mh4*pi*rho2 - complex<double>(0,4608)*mh6*pi*rho2 +
           1152*radixrho*mh6*pi*rho2 - 384*radixrho*mh8*pi*rho2 + 96*mh2*rho3 + 192*mh4*rho3 - 64*mh6*rho3 - 400*mh8*rho3 +
           complex<double>(0,96)*mh2*pi*rho3 - 48*radixrho*mh2*pi*rho3 + complex<double>(0,1152)*mh4*pi*rho3 - 288*radixrho*mh4*pi*rho3 +
           192*radixrho*mh6*pi*rho3 - rho4 - 96*mh2*rho4 + 288*mh4*rho4 - 960*mh6*rho4 + 1200*mh8*rho4 - complex<double>(0,96)*mh2*pi*rho4 +
           24*radixrho*mh2*pi*rho4 - 144*radixrho*mh4*pi*rho4 + 480*radixrho*mh6*pi*rho4 - 600*radixrho*mh8*pi*rho4))/
       (radixrho*power_of<3>(4*mh2 - rho)*rho)) + (-12*pi*(-1 + rho)*
       (complex<double>(0,512)*mh8*pi + 192*ln2*mh2*power_of<3>(4*mh2 - rho) - 96*lnrho*mh2*power_of<3>(4*mh2 - rho) - 768*mh6*rho + 1408*mh8*rho -
         complex<double>(0,384)*mh6*pi*rho + 480*mh4*rho2 - 1024*mh6*rho2 + 800*mh8*rho2 + complex<double>(0,96)*mh4*pi*rho2 - 96*mh2*rho3 + 288*mh4*rho3 -
         960*mh6*rho3 + 1200*mh8*rho3 - complex<double>(0,8)*mh2*pi*rho3 + rho4) +
      radixrho*(12*mh4*rho2*(-1392 + 384*rho + 215*rho2 + complex<double>(0,2)*radix4mh2*pi*(-240 + 48*rho + 95*rho2)) +
         8*mh2*(252 - 83*rho - complex<double>(0,9)*radix4mh2*pi*(-16 + 7*rho))*rho3 + 64*mh8*(-900 + 184*rho + 75*rho2 + 225*rho3) +
         16*mh6*rho*(3232 - 672*rho - 720*rho2 - 225*rho3 - complex<double>(0,9)*radix4mh2*pi*(-64 + 16*rho + 25*rho3)) +
         (-43.0 - complex<double>(0,12)*radix4mh2*pi)*rho4))/(4.*radixrho*power_of<3>(4*mh2 - rho)*rho) +
   (144*power_of<2>(atanrho)*mh2*(6*mh2*rho2*(complex<double>(0,-8)*(-1 + rho) + radixrho*(-4 + 2*rho + rho2)) +
        (-(radixrho*(-2 + rho)) + complex<double>(0,4)*(-1 + rho))*rho3 -
        4*mh4*rho*(complex<double>(0,-48)*(-1 + rho) + radixrho*(-24 + 12*rho + 2*rho2 + 5*rho3)) +
        mh6*(complex<double>(0,-256)*(-1 + rho) + radixrho*(-128 + 64*rho + 16*rho2 + 25*rho4))))/(radixrho*power_of<3>(4*mh2 - rho)*rho) +
   (6*atanh4mh2*(96*atan4mh2rho*mh2*power_of<3>(4*mh2 - rho)*(-1 + rho) -
        2*mh4*rho2*(radix4mh2*radixrho*(-240.0 + 48*rho + 95*rho2) +
           complex<double>(0,72)*pi*(complex<double>(0,8)*(-1 + rho) + radixrho*(-4 + 2*rho + rho2))) +
        6*mh2*(4*pi*(4.0 + complex<double>(0,1)*radixrho*(-2 + rho) - 4*rho) + radix4mh2*radixrho*(-16 + 7*rho))*rho3 +
        12*mh6*rho*(radix4mh2*radixrho*(-64 + 16*rho + 25*rho3) +
           complex<double>(0,8)*pi*(complex<double>(0,48)*(-1 + rho) + radixrho*(-24.0 + 12*rho + 2*rho2 + 5*rho3))) + radix4mh2*radixrho*rho4 -
        complex<double>(0,24)*mh8*pi*(complex<double>(0,256)*(-1 + rho) + radixrho*(-128 + 64*rho + 16*rho2 + 25*rho4))))/
    (radixrho*power_of<3>(4*mh2 - rho)*rho) + (complex<double>(0,288)*mh2*(-1 + rho)*dilogx13)/(radixrho*rho) +
   (complex<double>(0,144)*mh2*(-1 + rho)*dilogx4)/(radixrho*rho) - (complex<double>(0,288)*mh2*(-1 + rho)*redilogx12)/(radixrho*rho);
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J4
        complex<double>
        DileptonIntegralsCharm::j4(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            const double acotrho = pi / 2.0 - atanrho;

            // Begin of asymptotic part
            complex<double> asymp = (4*mh4*(-312 + 768*mh2 + (288*mh2)/rho - (3*lnmqmu*power_of<3>(4*mh2 - rho))/(mh4*rho) + 72*rho - (144*rho)/(-1 + radix4mh2) +
       720*atanh4mh2*radix4mh2*rho + (144*rho)/(1 + radix4mh2) + (24*rho)/mh2 - (72*atanh4mh2*radix4mh2*rho)/mh2 +
       288*mh2*rho - 864*power_of<2>(atanh4mh2)*mh2*rho - (72*mh2*rho)/power_of<2>(-1 + radix4mh2) + (216*mh2*rho)/(-1 + radix4mh2) -
       (72*mh2*rho)/power_of<2>(1 + radix4mh2) - (216*mh2*rho)/(1 + radix4mh2) - complex<double>(0,360)*radix4mh2*pi*rho +
       (complex<double>(0,36)*radix4mh2*pi*rho)/mh2 + complex<double>(0,864)*atanh4mh2*mh2*pi*rho + 216*mh2*pi2*rho +
       72*power_of<2>(acotrho)*(4*mh2*(-3 + rho) - 3*rho)*rho - 216*power_of<2>(atanh4mh2)*rho2 - (18*rho2)/power_of<2>(-1 + radix4mh2) +
       (54*rho2)/(-1 + radix4mh2) + 144*atanh4mh2*radix4mh2*rho2 - (18*rho2)/power_of<2>(1 + radix4mh2) - (54*rho2)/(1 + radix4mh2) -
       (2*rho2)/mh4 + (6*atanh4mh2*radix4mh2*rho2)/mh4 - (36*rho2)/mh2 +
       (12*atanh4mh2*radix4mh2*rho2)/mh2 + 288*power_of<2>(atanh4mh2)*mh2*rho2 - (16*mh2*rho2)/power_of<3>(-1 + radix4mh2) +
       (48*mh2*rho2)/power_of<2>(-1 + radix4mh2) + (96*mh2*rho2)/(-1 + radix4mh2) + (16*mh2*rho2)/power_of<3>(1 + radix4mh2) +
       (48*mh2*rho2)/power_of<2>(1 + radix4mh2) - (96*mh2*rho2)/(1 + radix4mh2) + complex<double>(0,216)*atanh4mh2*pi*rho2 -
       complex<double>(0,72)*radix4mh2*pi*rho2 - (complex<double>(0,3)*radix4mh2*pi*rho2)/mh4 - (complex<double>(0,6)*radix4mh2*pi*rho2)/mh2 -
       complex<double>(0,288)*atanh4mh2*mh2*pi*rho2 + 54*pi2*rho2 - 72*mh2*pi2*rho2 -
       (48*acotrho*radixrho*(3*rho*(-2 + 5*rho) + 4*mh2*(2 + rho + 3*rho2)))/rho))/(27.*power_of<3>(-1 + (4*mh2)/rho)*rho2);
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (-64*acotrho*radixrho*mh4*(2 + 18*mh4*(-2 + rho) - 5*rho - 12*mh2*(2 + rho))*rho2)/power_of<4>(-4*mh2 + rho) +
   (96*power_of<2>(acotrho)*mh4*rho2*(-4*mh2*(-4 + rho)*rho + rho2 + 2*mh4*(8 - 8*rho + 3*rho2)))/power_of<4>(-4*mh2 + rho) +
   (96*power_of<2>(atanh4mh2)*mh4*rho2*(-4*mh2*(-4 + rho)*rho + rho2 + 2*mh4*(8 - 8*rho + 3*rho2)))/power_of<4>(-4*mh2 + rho) +
   atanh4mh2*((complex<double>(0,-96)*mh4*pi*rho2*(-4*mh2*(-4 + rho)*rho + rho2 + 2*mh4*(8 - 8*rho + 3*rho2)))/power_of<4>(-4*mh2 + rho) +
      (16*radix4mh2*mh4*rho2*(8 - 24*rho - 9*rho2 + 2*mh2*(-40 - 24*rho + 9*rho2)))/power_of<4>(-4*mh2 + rho)) +
   (8*mh4*((2*rho2)/mh2 + 3.0*(complex<double>(0,1)*radix4mh2*pi*(-8 + 24*rho + 9*rho2) +
           (12 + 4*rho + (8 - 3*pi2)*rho2 + power_of<2>(1.0 - 4.0 * mh2)*(12 + 20*rho - 3*pi2*rho2) +
              (1 - 4*mh2)*(-24 - 24*rho + (-4 + 6*pi2)*rho2))/(16.*mh4)) +
        2*mh2*(complex<double>(0,-3)*radix4mh2*pi*(-40 - 24*rho + 9*rho2) -
           (-16 + 156*rho + 12*(-5 + 3*pi2)*rho2 - (8 + 9*pi2)*rho3 +
              power_of<3>(1.0 - 4.0 * mh2)*(16 - 108*rho - 36*(1 + pi2)*rho2 + 9*pi2*rho3) -
              3*power_of<2>(1.0 - 4.0 * mh2)*(16 - 124*rho - 4*(5 + 9*pi2)*rho2 + 3*(4 + 3*pi2)*rho3) +
              3*(1 - 4*mh2)*(16 - 140*rho + (12 - 36*pi2)*rho2 + (20 + 9*pi2)*rho3))/(32.*mh6*rho)) +
        mh4*(-18*pi2*(8 - 8*rho + 3*rho2) + (-4 - 24*rho - 15*rho2 + 43*rho3 +
              power_of<4>(1.0 - 4.0 * mh2)*(-4 - 24*rho - 63*rho2 + 27*rho3) - 16*rho4 -
              6*power_of<2>(1.0 - 4.0 * mh2)*(4 + 24*rho + 47*rho2 - 59*rho3 + 15*rho4) +
              power_of<3>(1.0 - 4.0 * mh2)*(16 + 96*rho + 228*rho2 - 180*rho3 + 27*rho4) + (1 - 4*mh2)*(16 + 96*rho + 132*rho2 - 244*rho3 + 91*rho4))
             /(32.*mh8*rho2))))/(3.*power_of<4>(1 - (4*mh2)/rho)*rho2);
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (-128*acotrho*radixrho*mh4*rho2*(rho*(-2 + 5*rho) + mh4*(96 + 228*rho - 90*rho2) + 4*mh2*(-2 + 17*rho + 6*rho2) +
        8*mh6*(26 - 35*rho + 15*rho2)))/power_of<5>(4*mh2 - rho) +
   (192*power_of<2>(acotrho)*mh4*rho2*(4*mh2*(-9 + 2*rho)*rho2 - 6*mh4*rho*(24 - 16*rho + 5*rho2) - rho3 +
        8*mh6*(-8 + 16*rho - 15*rho2 + 5*rho3)))/power_of<5>(4*mh2 - rho) +
   (192*power_of<2>(atanh4mh2)*mh4*rho2*(4*mh2*(-9 + 2*rho)*rho2 - 6*mh4*rho*(24 - 16*rho + 5*rho2) - rho3 +
        8*mh6*(-8.0 + 16*rho - 15*rho2 + 5*rho3)))/power_of<5>(4*mh2 - rho) +
   atanh4mh2*((complex<double>(0,-192)*mh4*pi*rho2*(4*mh2*(-9 + 2*rho)*rho2 - 6*mh4*rho*(24 - 16*rho + 5*rho2) - rho3 +
           8*mh6*(-8.0 + 16*rho - 15*rho2 + 5*rho3)))/power_of<5>(4*mh2 - rho) +
      (32*radix4mh2*mh4*rho2*(rho*(-8 + 24*rho + 13*rho2) + mh2*(-32 + 272*rho + 228*rho2 - 70*rho3) +
           8*mh4*(40.0 + 48*rho - 45*rho2 + 15*rho3)))/power_of<5>(4*mh2 - rho)) +
   (16*mh4*((-5*rho2)/mh2 + 5.0*(complex<double>(0,-3)*radix4mh2*pi*(-8 + 24*rho + 13*rho2) +
           (-36 - 52*rho + 3*(-8 + 3*pi2)*rho2 + power_of<2>(1.0 - 4.0 * mh2)*(-36 - 100*rho + 9*pi2*rho2) +
              2*(1 - 4*mh2)*(36.0 + 76*rho + (6.0 - 9*pi2)*rho2))/(16.*mh4)) +
        (10.0*mh2*(complex<double>(0,3)*radix4mh2*pi*(16 - 136*rho - 114*rho2 + 35*rho3) +
             (-68 + 228*rho + 9*(-16 + 9*pi2)*rho2 - 2*(8 + 9*pi2)*rho3 +
                power_of<3>(1.0 - 4.0 * mh2)*(68 - 84*rho - 9*(8 + 9*pi2)*rho2 + 18*pi2*rho3) -
                3*power_of<2>(1.0 - 4.0 * mh2)*(68 - 132*rho - 9*(4 + 9*pi2)*rho2 + 6*(4 + 3*pi2)*rho3) +
                3*(1 - 4*mh2)*(68 - 180*rho - 9*(-4 + 9*pi2)*rho2 + 2*(20 + 9*pi2)*rho3))/(32.*mh6)))/rho +
        4*mh6*((-90*pi2*(-8 + 16*rho - 15*rho2 + 5*rho3))/rho -
           (-4 - 60*rho + 184*rho5 - 295*rho2 + 760*rho3 +
              5*power_of<2>(1.0 - 4.0 * mh2)*(-8 - 120*rho + 323*rho5 - 854*rho2 + 1856*rho3 - 1355*rho4) -
              5*(1 - 4*mh2)*(-4 - 60*rho + 229*rho5 - 367*rho2 + 904*rho3 - 760*rho4) +
              5*power_of<4>(1.0 - 4.0 * mh2)*(-4 - 60*rho + 45*rho5 - 511*rho2 + 744*rho3 - 360*rho4) - 625*rho4 +
              power_of<5>(1.0 - 4.0 * mh2)*(4 + 60*rho + 535*rho2 - 600*rho3 + 225*rho4) +
              power_of<3>(1.0 - 4.0 * mh2)*(40 + 600*rho - 975*rho5 + 4750*rho2 - 8640*rho3 + 5175*rho4))/(128.*mh10*rho3))\
         + (10.0*mh4*(complex<double>(0,-12)*radix4mh2*pi*rho*(40 + 48*rho - 45*rho2 + 15*rho3) +
             (-48 + 1728*rho + 36*(-5 + 18*pi2)*rho2 - 12*(77 + 36*pi2)*rho3 + 5*(64 + 27*pi2)*rho4 +
                3*power_of<4>(1.0 - 4.0 * mh2)*(-16 + 448*rho + 12*(43 + 18*pi2)*rho2 - 36*(5 + 4*pi2)*rho3 + 45*pi2*rho4) -
                12*power_of<3>(1.0 - 4.0 * mh2)*(-16 + 480*rho + 12*(37 + 18*pi2)*rho2 - 36*(9 + 4*pi2)*rho3 +
                   45*(1 + pi2)*rho4) + 18*power_of<2>(1.0 - 4.0 * mh2)*
                 (-16 + 512*rho + 108*(3 + 2*pi2)*rho2 - 4*(109 + 36*pi2)*rho3 + 5*(20 + 9*pi2)*rho4) -
                4*(1 - 4*mh2)*(-48 + 1632*rho + 36*(13 + 18*pi2)*rho2 - 12*(113 + 36*pi2)*rho3 +
                   5*(91 + 27*pi2)*rho4))/(256.*mh8)))/rho2))/(15.*power_of<5>(-1 + (4*mh2)/rho)*rho2);
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J5
        complex<double>
        DileptonIntegralsCharm::j5(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            static const double ln2 = std::log(2.0);
            const double acotrho = pi / 2.0 - atanrho;
            // There has been an error in the computation of bminus for this integral.
            // Replacing bminus by bplus fixes this error:
            const double lnbm = std::log(0.5 * (1.0 - radix4mh2));

            // Begin of asymptotic part
            complex<double> asymp = (4*mh2*(24 - 24*lnrho - 24*lnmqmu*lnrho + 48*ln2 + 48*lnmqmu*ln2 - 184*mh2 + (48*mh2)/(-1 + radix4mh2) -
        (48*mh2)/(1 + radix4mh2) - complex<double>(0,64)*radix4mh2*pi + complex<double>(0,64)*radix4mh2*mh2*pi + 12*pi2 + (104*mh2)/rho -
        (48*lnmqmu*mh2)/rho + 40*rho + (12*rho)/(-1 + radix4mh2) - (12*rho)/(1 + radix4mh2) - (2*rho)/mh2 +
        (3*lnmqmu*rho)/mh2 + (6*mh2*rho)/power_of<2>(-1 + radix4mh2) - (18*mh2*rho)/(-1 + radix4mh2) +
        (6*mh2*rho)/power_of<2>(1 + radix4mh2) + (18*mh2*rho)/(1 + radix4mh2) + complex<double>(0,30)*radix4mh2*pi*rho -
        (complex<double>(0,3)*radix4mh2*pi*rho)/mh2 - 18*mh2*pi2*rho - 16*lnbm*(-5 + 12*mh2 + 3*rho) +
        16*lnmh*(-2 + 3*lnmqmu + 12*mh2 + 3*rho)))/(9.*power_of<3>(-1 + (4*mh2)/rho)*rho) +
   (64*acotrho*radixrho*mh2*rho*(2*(-4 + rho)*rho + 3*mh2*(-2 + 5*rho)))/(9.*power_of<3>(4*mh2 - rho)) +
   (32*power_of<2>(acotrho)*mh2*(-2 + 3*mh2*rho)*rho2)/(3.*power_of<3>(4*mh2 - rho)) +
   (32*power_of<2>(atanh4mh2)*mh2*(-2 + 3*mh2*rho)*rho2)/(3.*power_of<3>(4*mh2 - rho)) +
   atanh4mh2*((8*radix4mh2*(-64*mh4 + mh2*(64 - 30*rho) + 3*rho)*rho2)/(9.*power_of<3>(4*mh2 - rho)) +
      (32.0*mh2*(10.0 + complex<double>(0,6)*pi - 6*rho + mh2*(-24.0 - complex<double>(0,9)*pi*rho))*rho2)/(9.*power_of<3>(4*mh2 - rho)));
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (4*mh2*(-108*lnmqmu + 72*lnrho + 72*lnmqmu*lnrho - 144*ln2 - 144*lnmqmu*ln2 + 216*mh2 - (576*mh2)/(-1 + radix4mh2) +
        (576*mh2)/(1 + radix4mh2) + 576*mh4 - (216*mh4)/power_of<2>(-1 + radix4mh2) + (648*mh4)/(-1 + radix4mh2) -
        (216*mh4)/power_of<2>(1 + radix4mh2) - (648*mh4)/(1 + radix4mh2) + complex<double>(0,300)*radix4mh2*pi -
        complex<double>(0,1272)*radix4mh2*mh2*pi - 36*pi2 + 648*mh4*pi2 - (672*mh4)/power_of<2>(rho) +
        (192*lnmqmu*mh4)/power_of<2>(rho) - (1224*mh2)/rho + (432*lnmqmu*mh2)/rho + (288*lnrho*mh2)/rho + (288*lnmqmu*lnrho*mh2)/rho -
        (576*ln2*mh2)/rho - (576*lnmqmu*ln2*mh2)/rho + (3744*mh4)/rho + (complex<double>(0,768)*radix4mh2*mh2*pi)/rho -
        (complex<double>(0,768)*radix4mh2*mh4*pi)/rho - (144*mh2*pi2)/rho - 192*rho - (36*rho)/(-1 + radix4mh2) +
        (36*rho)/(1 + radix4mh2) + (2*rho)/mh2 - (3*lnmqmu*rho)/mh2 - (54*mh2*rho)/power_of<2>(-1 + radix4mh2) +
        (162*mh2*rho)/(-1 + radix4mh2) - (54*mh2*rho)/power_of<2>(1 + radix4mh2) - (162*mh2*rho)/(1 + radix4mh2) -
        (32*mh4*rho)/power_of<3>(-1 + radix4mh2) + (96*mh4*rho)/power_of<2>(-1 + radix4mh2) + (192*mh4*rho)/(-1 + radix4mh2) +
        (32*mh4*rho)/power_of<3>(1 + radix4mh2) + (96*mh4*rho)/power_of<2>(1 + radix4mh2) - (192*mh4*rho)/(1 + radix4mh2) -
        complex<double>(0,102)*radix4mh2*pi*rho + (complex<double>(0,3)*radix4mh2*pi*rho)/mh2 - complex<double>(0,144)*radix4mh2*mh2*pi*rho +
        162*mh2*pi2*rho - 144*mh4*pi2*rho + (576*mh4)/(rho - radix4mh2*rho) + (576*mh4)/(rho + radix4mh2*rho) +
        (48*lnbm*(48*mh4 + rho*(-5 + 3*rho) + 4*mh2*(-5 + 12*rho)))/rho -
        (48*lnmh*(48*mh4 + rho*(-2 + 3*lnmqmu + 3*rho) + 4*mh2*(-2 + 3*lnmqmu + 12*rho)))/rho))/(9.*power_of<4>(1 - (4*mh2)/rho)*rho) +
   (32*power_of<2>(acotrho)*mh2*(2*rho + 4*mh4*rho*(-9 + 2*rho) + mh2*(8 - 9*rho2))*rho2)/power_of<4>(-4*mh2 + rho) +
   (32*power_of<2>(atanh4mh2)*mh2*(2*rho + 4*mh4*rho*(-9 + 2*rho) + mh2*(8 - 9*rho2))*rho2)/power_of<4>(-4*mh2 + rho) -
   (64*acotrho*radixrho*mh2*rho*(mh2*rho*(-50 + 53*rho) + 2*(-4 + rho)*rho2 + 4*mh4*(-2 + 17*rho + 6*rho2)))/(3.*power_of<4>(-4*mh2 + rho)) +
   atanh4mh2*((8*radix4mh2*rho2*(256*mh6 + 2*mh2*rho*(-50 + 17*rho) - rho2 + 8*mh4*(-32 + 53*rho + 6*rho2)))/(3.*power_of<4>(-4*mh2 + rho)) +
      (32*mh2*rho2*(2.0*rho*(-5.0 - complex<double>(0,3)*pi + 3*rho) + 12.0*mh4*(8.0 - complex<double>(0,1)*pi*rho*(-9 + 2*rho)) +
           mh2*(-40.0 + 96*rho + complex<double>(0,3)*pi*(-8 + 9*rho2))))/(3.*power_of<4>(-4*mh2 + rho)));
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (4*mh2*(-80 + 336*lnmqmu - 144*lnrho - 144*lnmqmu*lnrho + 288*ln2 + 288*lnmqmu*ln2 + 48*radix4mh2*ln2 -
        48*ln2*radix4mh2 + 3936*mh2 + (2592*mh2)/(-1 + radix4mh2) - (2592*mh2)/(1 + radix4mh2) - 192*radix4mh2*ln2*mh2 +
        192*ln2*radix4mh2*mh2 - 5760*mh4 + (2592*mh4)/power_of<2>(-1 + radix4mh2) - (7776*mh4)/(-1 + radix4mh2) +
        (2592*mh4)/power_of<2>(1 + radix4mh2) + (7776*mh4)/(1 + radix4mh2) + 6480*mh6 + (1280*mh6)/power_of<3>(-1 + radix4mh2) -
        (3840*mh6)/power_of<2>(-1 + radix4mh2) - (7680*mh6)/(-1 + radix4mh2) - (1280*mh6)/power_of<3>(1 + radix4mh2) -
        (3840*mh6)/power_of<2>(1 + radix4mh2) + (7680*mh6)/(1 + radix4mh2) - complex<double>(0,768)*radix4mh2*pi +
        complex<double>(0,48)*radix4mh2*pi + complex<double>(0,6816)*radix4mh2*mh2*pi - complex<double>(0,192)*radix4mh2*mh2*pi +
        complex<double>(0,5760)*radix4mh2*mh4*pi + 72*pi2 - 7776*mh4*pi2 + 5760*mh6*pi2 + (3264*mh6)/power_of<3>(rho) -
        (768*lnmqmu*mh6)/power_of<3>(rho) + (16512*mh4)/power_of<2>(rho) - (5376*lnmqmu*mh4)/power_of<2>(rho) - (2304*lnrho*mh4)/power_of<2>(rho) -
        (2304*lnmqmu*lnrho*mh4)/power_of<2>(rho) + (4608*ln2*mh4)/power_of<2>(rho) + (4608*lnmqmu*ln2*mh4)/power_of<2>(rho) +
        (768*radix4mh2*ln2*mh4)/power_of<2>(rho) - (768*ln2*radix4mh2*mh4)/power_of<2>(rho) - (38784*mh6)/power_of<2>(rho) +
        (4608*mh6)/((-1 + radix4mh2)*power_of<2>(rho)) - (4608*mh6)/((1 + radix4mh2)*power_of<2>(rho)) -
        (3072*radix4mh2*ln2*mh6)/power_of<2>(rho) + (3072*ln2*radix4mh2*mh6)/power_of<2>(rho) -
        (complex<double>(0,6912)*radix4mh2*mh4*pi)/power_of<2>(rho) + (complex<double>(0,768)*radix4mh2*mh4*pi)/power_of<2>(rho) +
        (complex<double>(0,9216)*radix4mh2*mh6*pi)/power_of<2>(rho) - (complex<double>(0,3072)*radix4mh2*mh6*pi)/power_of<2>(rho) +
        (1152*mh4*pi2)/power_of<2>(rho) + (4320*mh2)/rho - (1728*lnrho*mh2)/rho - (1728*lnmqmu*lnrho*mh2)/rho +
        (3456*ln2*mh2)/rho + (3456*lnmqmu*ln2*mh2)/rho + (576*radix4mh2*ln2*mh2)/rho - (576*ln2*radix4mh2*mh2)/rho -
        (27648*mh4)/rho + (10368*mh4)/((-1 + radix4mh2)*rho) - (2304*radix4mh2*ln2*mh4)/rho + (2304*ln2*radix4mh2*mh4)/rho -
        (20880*mh6)/rho + (3456*mh6)/(power_of<2>(-1 + radix4mh2)*rho) + (3456*mh6)/(power_of<2>(1 + radix4mh2)*rho) -
        (complex<double>(0,6912)*radix4mh2*mh2*pi)/rho + (complex<double>(0,576)*radix4mh2*mh2*pi)/rho + (complex<double>(0,24192)*radix4mh2*mh4*pi)/rho -
        (complex<double>(0,2304)*radix4mh2*mh4*pi)/rho + (864*mh2*pi2)/rho - (10368*mh6*pi2)/rho + 480*rho +
        (72*rho)/(-1 + radix4mh2) - (72*rho)/(1 + radix4mh2) - (2*rho)/mh2 + (3*lnmqmu*rho)/mh2 +
        (216*mh2*rho)/power_of<2>(-1 + radix4mh2) - (648*mh2*rho)/(-1 + radix4mh2) + (216*mh2*rho)/power_of<2>(1 + radix4mh2) +
        (648*mh2*rho)/(1 + radix4mh2) + (320*mh4*rho)/power_of<3>(-1 + radix4mh2) - (960*mh4*rho)/power_of<2>(-1 + radix4mh2) -
        (1920*mh4*rho)/(-1 + radix4mh2) - (320*mh4*rho)/power_of<3>(1 + radix4mh2) - (960*mh4*rho)/power_of<2>(1 + radix4mh2) +
        (1920*mh4*rho)/(1 + radix4mh2) + (180*mh6*rho)/power_of<4>(-1 + radix4mh2) - (600*mh6*rho)/power_of<3>(-1 + radix4mh2) +
        (270*mh6*rho)/power_of<2>(-1 + radix4mh2) + (2970*mh6*rho)/(-1 + radix4mh2) + (180*mh6*rho)/power_of<4>(1 + radix4mh2) +
        (600*mh6*rho)/power_of<3>(1 + radix4mh2) + (270*mh6*rho)/power_of<2>(1 + radix4mh2) - (2970*mh6*rho)/(1 + radix4mh2) +
        complex<double>(0,210)*radix4mh2*pi*rho - (complex<double>(0,3)*radix4mh2*pi*rho)/mh2 + complex<double>(0,1170)*radix4mh2*mh2*pi*rho -
        complex<double>(0,1620)*radix4mh2*mh4*pi*rho - 648*mh2*pi2*rho + 1440*mh4*pi2*rho - 1620*mh6*pi2*rho +
        (10368*mh6)/(rho - radix4mh2*rho) - (10368*mh4)/(rho + radix4mh2*rho) + (10368*mh6)/(rho + radix4mh2*rho) +
        (48*lnbm*(64*(-6 - radix4mh2 + radix4mh2)*mh6 +
             4*mh2*rho*(30 - 3*radix4mh2 - radix4mh2*(-3 + rho) + (-54 + radix4mh2)*rho) +
             16*mh4*(10 + radix4mh2 - radix4mh2 - 3*radix4mh2*rho + 3*(-18 + radix4mh2)*rho) +
             (10 + radix4mh2 - radix4mh2 - 6*rho)*rho2))/power_of<2>(rho) +
        (96*lnmh*(192*mh6 + 12*mh2*rho*(-2 + 3*lnmqmu + 9*rho) + 16*mh4*(-2 + 3*lnmqmu + 27*rho) + (-2 + 3*lnmqmu + 3*rho)*rho2))/
         power_of<2>(rho)))/(9.*power_of<5>(-1 + (4*mh2)/rho)*rho) +
   (64*power_of<2>(acotrho)*mh2*rho2*(-2*rho2 + 6*mh2*rho*(-4 + 3*rho2) + mh6*rho*(288 - 160*rho + 45*rho2) - 8*mh4*(4 - 27*rho2 + 5*rho3)))/
    power_of<5>(4*mh2 - rho) + (64*power_of<2>(atanh4mh2)*mh2*rho2*
      (-2*rho2 + 6*mh2*rho*(-4 + 3*rho2) + mh6*rho*(288 - 160*rho + 45*rho2) - 8*mh4*(4 - 27*rho2 + 5*rho3)))/power_of<5>(4*mh2 - rho) -
   (128*acotrho*radixrho*mh2*rho*(-6*mh2*(-22 + 19*rho)*rho2 - 24*mh4*rho*(-10 + 23*rho + 5*rho2) - 2*(-4 + rho)*rho3 +
        mh6*(16 - 280*rho - 390*rho2 + 135*rho3)))/(3.*power_of<5>(4*mh2 - rho)) +
   atanh4mh2*((8*radix4mh2*rho2*(-3072*mh8 + 2*mh2*(128 - 35*rho)*rho2 - 2*mh4*rho*(-1152 + 1136*rho + 195*rho2) + rho3 +
           12*mh6*(192 - 672*rho - 160*rho2 + 45*rho3)))/(3.*power_of<5>(4*mh2 - rho)) +
      (64*mh2*rho2*(-2.0*(-5.0 + radix4mh2 - complex<double>(0,3)*pi + 3*rho)*rho2 +
           2*mh2*rho*(4.0*(15.0 - 3*radix4mh2 - 27*rho + radix4mh2*rho) - complex<double>(0,9)*pi*(-4 + 3*rho2)) +
           mh6*(128*(-3 + radix4mh2) - complex<double>(0,3)*pi*rho*(288 - 160*rho + 45*rho2)) +
           8*mh4*(4*(5 - radix4mh2 + 3.0*(-9.0 + radix4mh2)*rho) + complex<double>(0,3)*pi*(4 - 27*rho2 + 5*rho3))))/(3.*power_of<5>(4*mh2 - rho))
      );
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // J6
        complex<double>
        DileptonIntegralsCharm::j6(const double & a1, const double & a2) const
        {
            static const double pi = M_PI, pi2 = pi * pi;
            const double acotrho = pi / 2.0 - atanrho;

            // Begin of asymptotic part
            complex<double> asymp = (-32*power_of<2>(acotrho)*(-3 + 2*mh2)*mh4)/(3.*power_of<3>(-1 + (4*mh2)/rho)) +
   (4*mh2*((-3*lnmqmu*power_of<3>(4*mh2 - rho))/(power_of<2>(mh)*rho2) + ((-2.0 - complex<double>(0,3)*radix4mh2*pi)*rho)/power_of<2>(mh) +
        6.0*(10.0 + (-(-3.0 + 9*(1 - 4*mh2))/(4.*mh2) + complex<double>(0,8)*radix4mh2*pi)*rho) +
        6.0*mh2*(-60.0 - 28/rho + 3.0*(-(-1.0 - 4*mh2)/(4.*power_of<2>(mh2)) + complex<double>(0,2)*radix4mh2*pi - 3*pi2)*rho) +
        4*mh4*(-36.0 + 48/rho2 + 84/rho - ((-8 - 9*pi2 + 9*power_of<3>(1.0 - 4.0 * mh2)*pi2 -
                9*power_of<2>(1 - 4*mh2)*(4 + 3*pi2) + 3*(1 - 4*mh2)*(20 + 9*pi2))*rho)/(64.*power_of<3>(mh2)))))/
    (27.*power_of<3>(-1 + (4*mh2)/rho)*rho) + (32*acotrho*radixrho*mh2*
      (6*mh2*rho*(2 + rho) + 3*(-1 + rho)*rho2 + 4*mh4*(-4 - 2*rho + 3*rho2)))/(9.*power_of<3>(4*mh2 - rho)) -
   (32*power_of<2>(atanh4mh2)*(-3 + 2*mh2)*mh4*rho3)/(3.*power_of<3>(4*mh2 - rho)) +
   atanh4mh2*((-8*radix4mh2*(-1 + 16*mh2 + 12*mh4)*rho3)/(9.*power_of<3>(4*mh2 - rho)) +
      (complex<double>(0,10.666666666666666)*(-3 + 2*mh2)*mh4*pi*rho3)/power_of<3>(4*mh2 - rho));
            // End of asymptotic part

            // Begin of 1st Gegenbauer moment
            complex<double> gb1 = (32*acotrho*radixrho*mh2*(-18*mh2*rho - 36*mh4*rho - (-1 + rho)*rho + 12*mh6*(-2 + 3*rho))*rho2)/(3.*power_of<4>(-4*mh2 + rho)) -
   (32*power_of<2>(acotrho)*mh4*(-6*mh2*(-2 + rho) + 3*rho + mh4*(-8 + 6*rho))*rho3)/power_of<4>(-4*mh2 + rho) -
   (32*power_of<2>(atanh4mh2)*mh4*(-6*mh2*(-2 + rho) + 3*rho + mh4*(-8 + 6*rho))*rho3)/power_of<4>(-4*mh2 + rho) +
   atanh4mh2*((complex<double>(0,32)*mh4*pi*(-6*mh2*(-2 + rho) + 3*rho + mh4*(-8 + 6*rho))*rho3)/power_of<4>(-4*mh2 + rho) +
      (16*radix4mh2*mh2*(-2 + 3*rho - 6*mh4*(-4 + 3*rho) + mh2*(32 + 15*rho))*rho3)/(3.*power_of<4>(-4*mh2 + rho))) +
   (8*mh2*(-6 - 17*rho - (3*rho)/(-1 + radix4mh2) + (3*rho)/(1 + radix4mh2) - complex<double>(0,3)*radix4mh2*pi*(-2 + 3*rho) +
        3*mh2*(complex<double>(0,-1)*radix4mh2*pi*(32 + 15*rho) +
           3*((-2*(-1 - 4*mh2))/mh2 - 4/rho + ((-8 + 3*pi2 + 3*power_of<2>(1 - 4*mh2)*pi2 + (1 - 4*mh2)*(4 - 6*pi2))*
                 rho)/(16.*power_of<2>(mh2)))) + 6*mh4*(-9*pi2*(-2 + rho) + complex<double>(0,3)*radix4mh2*pi*(-4 + 3*rho) -
           (-2 - 9*rho + 3*rho2 + power_of<3>(1.0 - 4.0 * mh2)*(2 + 9*rho + 9*rho2) + (1 - 4*mh2)*(6 + 27*rho + 9*rho2 - 15*rho3) + 2*rho3 +
              3*power_of<2>(1 - 4*mh2)*(-2 - 9*rho - 7*rho2 + 3*rho3))/(16.*power_of<3>(mh2)*rho2)) +
        mh6*(18*pi2*(-4 + 3*rho) - (4 + 12*rho - 27*rho2 + 35*rho3 + power_of<4>(1.0 - 4.0 * mh2)*(4 + 12*rho - 27*rho2 + 27*rho3) -
              16*rho4 - 6*power_of<2>(1 - 4*mh2)*(-4 - 12*rho + 27*rho2 - 43*rho3 + 15*rho4) +
              power_of<3>(1.0 - 4.0 * mh2)*(-16 - 48*rho + 108*rho2 - 144*rho3 + 27*rho4) +
              (1 - 4*mh2)*(-16 - 48*rho + 108*rho2 - 176*rho3 + 91*rho4))/(32.*power_of<4>(mh2)*rho3))))/
    (9.*power_of<4>(1 - (4*mh2)/rho)*rho);
            // End of 1st Gegenbauer moment

            // Begin of 2nd Gegenbauer moment
            complex<double> gb2 = (32*acotrho*radixrho*mh2*rho2*(24*mh6*(22 - 15*rho)*rho + 144*mh4*rho*(1 + rho) + 4*mh2*rho*(-2 + 11*rho) + (-1 + rho)*rho2 +
        24*mh8*(4 - 20*rho + 15*rho2)))/(3.*power_of<5>(4*mh2 - rho)) -
   (64*power_of<2>(acotrho)*mh4*(12*mh2*(-3 + rho)*rho + mh4*(-48 + 64*rho - 30*rho2) - 3*rho2 + mh6*(32 - 60*rho + 30*rho2))*rho3)/
    power_of<5>(4*mh2 - rho) - (64*power_of<2>(atanh4mh2)*mh4*(12*mh2*(-3 + rho)*rho + mh4*(-48 + 64*rho - 30*rho2) - 3*rho2 +
        mh6*(32 - 60*rho + 30*rho2))*rho3)/power_of<5>(4*mh2 - rho) +
   atanh4mh2*((-16*radix4mh2*mh2*(rho*(-2 + 3*rho) + mh4*(256 + 324*rho - 150*rho2) + 4*mh2*(-4 + 23*rho + 12*rho2) +
           12*mh6*(16 - 30*rho + 15*rho2))*rho3)/(3.*power_of<5>(4*mh2 - rho)) +
      (complex<double>(0,64)*mh4*pi*(12*mh2*(-3 + rho)*rho + mh4*(-48 + 64*rho - 30*rho2) - 3*rho2 + mh6*(32 - 60*rho + 30*rho2))*rho3)/
       power_of<5>(4*mh2 - rho)) + (4.0*mh2*(5.0*(12 - ((-31 + 43*(1 - 4*mh2))*rho)/(4.*mh2) + complex<double>(0,6)*radix4mh2*pi*(-2 + 3*rho)) +
        20*mh2*(-(52 + 44*(1 - 4*mh2))/(4.*mh2) + 60/rho -
           (9*(-8 + 3*pi2 + 3*power_of<2>(1 - 4*mh2)*pi2 + (1 - 4*mh2)*(4 - 6*pi2))*rho)/(16.*power_of<2>(mh2)) +
           (complex<double>(0,6)*radix4mh2*pi*(-4 + 23*rho + 12*rho2))/rho) +
        (60*mh4*(complex<double>(0,-1)*radix4mh2*pi*rho*(-128 - 162*rho + 75*rho2) -
             (-4 + 84*rho + 9*(-4 + 3*pi2)*rho2 - (8 + 9*pi2)*rho3 +
                power_of<3>(1.0 - 4.0 * mh2)*(4 - 60*rho - 9*(4 + 3*pi2)*rho2 + 9*pi2*rho3) -
                3*power_of<2>(1 - 4*mh2)*(4 - 68*rho - 3*(8 + 9*pi2)*rho2 + 3*(4 + 3*pi2)*rho3) +
                3*(1 - 4*mh2)*(4 - 76*rho - 27*pi2*rho2 + (20 + 9*pi2)*rho3))/(16.*power_of<3>(mh2))))/rho2 +
        20*mh6*(pi2*(576 - 432/rho - 270*rho) + (complex<double>(0,18)*radix4mh2*pi*(16 - 30*rho + 15*rho2))/rho +
           (-4 - 48*rho - 99*rho2 + 199*rho3 + power_of<4>(1.0 - 4.0 * mh2)*(-4 - 48*rho - 243*rho2 + 135*rho3) - 80*rho4 -
              6*power_of<2>(1 - 4*mh2)*(4 + 48*rho + 195*rho2 - 263*rho3 + 75*rho4) +
              power_of<3>(1.0 - 4.0 * mh2)*(16 + 192*rho + 900*rho2 - 828*rho3 + 135*rho4) +
              (1 - 4*mh2)*(16 + 192*rho + 612*rho2 - 1084*rho3 + 455*rho4))/(32.*power_of<4>(mh2)*rho3)) +
        8*mh8*(45*pi2*(-30 + 16/rho + 15*rho) + (8 + 60*rho + 552*rho5 - 390*rho2 + 1285*rho3 +
              5*power_of<2>(1 - 4*mh2)*(16 + 120*rho + 969*rho5 - 780*rho2 + 2906*rho3 - 3160*rho4) -
              5*power_of<3>(1.0 - 4.0 * mh2)*(16 + 120*rho + 585*rho5 - 780*rho2 + 2778*rho3 - 2520*rho4) -
              5*(1 - 4*mh2)*(8 + 60*rho + 687*rho5 - 390*rho2 + 1429*rho3 - 1745*rho4) +
              5*power_of<4>(1.0 - 4.0 * mh2)*(8 + 60*rho + 135*rho5 - 390*rho2 + 1269*rho3 - 945*rho4) - 1475*rho4 +
              power_of<10>(radix4mh2)*(-8 - 60*rho + 390*rho2 - 1125*rho3 + 675*rho4))/(256.*power_of<5>(mh2)*rho4))))/
    (45.*power_of<5>(-1 + (4*mh2)/rho)*rho);
            // End of 2nd Gegenbauer moment

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // Massless case

        // cf. [vD:2011A], Eq. (xx), p. ?
        inline complex<double> j2_massless(const double & sh, const double & a1, const double & a2)
        {
            static const double pi2 = M_PI * M_PI;

            double lnsh = std::log(sh), ln1msh = std::log(1.0 - sh);
            double atanhsh = std::atanh(1.0 - 2.0 * sh);
            complex<double> dilogsh = dilog(sh);

            // asymptotic part
            complex<double> asymp = ((6 + pi2) * (1.0 - sh) + 3.0 * lnsh * (2.0 - 2.0 * (1.0 - sh) * ln1msh + (1 - sh) * lnsh) - 6.0 * (1.0 - sh) * dilogsh)
                / (sh - 1.0);

            // first Gegenbauer moment
            complex<double> gb1 = -3.0 * ((-1.0 + sh) * (-15.0 + pi2 * (-1.0 + sh) + 9.0 * sh)
                    - 6.0 * (-2.0 + sh + power_of<2>(1.0 - sh) * (ln1msh - lnsh)) * lnsh
                    - 3.0 * power_of<2>(1.0 - sh) * (lnsh * lnsh + 2.0 * dilogsh))
                / power_of<2>(1.0 - sh);

            // second Gegenbauer moment
            complex<double> gb2 = -2.0 * ((-1.0 + sh) * (73.0 + 3.0 * pi2 * power_of<2>(1.0 - sh) + sh * (-71.0 + 28.0 * sh))
                    - 6.0 * (8.0 + 3.0 * (-2.0 + sh) * sh) * lnsh
                    - 9.0 * power_of<3>(-1.0 + sh) * (lnsh * (4.0 * atanhsh + lnsh) + 2.0 * dilogsh))
                / power_of<3>(-1.0 + sh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // cf. [vD:2011A], Eq. (xx), p. ?
        inline complex<double> j3_massless(const double & sh, const double & a1, const double & a2)
        {
            static const double pi = M_PI, pi2 = pi * pi;

            double sh2 = sh * sh;
            double lnsh = std::log(sh), ln1msh = std::log(1.0 - sh);
            complex<double> dilogsh = dilog(sh);

            // asymptotic part
            complex<double> asymp = ((1.0 - sh) * (-9.0 + (15.0 + 2.0 * pi2) * sh)
                    + 6.0 * (-1.0 + 2.0 * sh + (-1.0 + sh) * sh * (2.0 * ln1msh - lnsh)) * lnsh
                    + 12.0 * (-1.0 + sh) * sh * dilogsh)
                / 2.0 / (-1.0 + sh);

            // first Gegenbauer moment
            complex<double> gb1 = ((1.0 - sh) * (17.0 + sh * (-82.0 + 6.0 * pi2 * (-1.0 + sh) + 53.0 * sh))
                    + 6.0 * (1.0 - 9.0 * sh + 6.0 * sh2 + 3.0 * power_of<2>(-1.0 + sh) * sh * (2.0 * ln1msh - lnsh)) * lnsh
                    + 36.0 * power_of<2>(-1.0 + sh) * sh * dilogsh)
                / 2.0 / power_of<2>(-1.0 + sh);

            // second Gegenbauer moment
            complex<double> gb2 = ((1.0 - sh) * (-43.0 + sh * (461.0 + 24.0 * pi2 * power_of<2>(1.0 - sh) + sh * (-583.0 + 225.0 * sh)))
                    + 12.0 * (-1.0 + 6.0 * sh * (4.0 + sh * (-5.0 + 2.0 * sh)) + 6.0 * power_of<3>(-1.0 + sh) * sh * (2.0 * ln1msh - lnsh)) * lnsh
                    + 144.0 * power_of<3>(-1.0 + sh) * sh * dilogsh)
                / 4.0 / power_of<3>(-1.0 + sh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // cf. [vD:2011A], Eq. (xx), p. ?
        inline complex<double> j4_massless(const double & sh, const double & mB, const double & mu, const double & a1, const double & a2)
        {
            static const double pi = M_PI;

            double sh2 = sh * sh, sh3 = sh * sh2, sh4 = sh2 * sh2;
            double lnsh = std::log(sh), lnmbmu = std::log(mB / mu);
//            complex<double> dilogsh = dilog(sh);

            complex<double> asymp = 2.0 / 9.0 * ((1 - sh) * complex<double>(3.0 - 10.0 * sh + 3.0 * sh2, 2.0 * pi * power_of<2>(1.0 - sh))
                    - 4.0 * power_of<3>(1.0 - sh) * lnmbmu - 2.0 * (3.0 - sh) * sh2 * lnsh)
                / power_of<3>(1.0 - sh);
            complex<double> gb1 = (1.0 - 8.0 * sh + 8.0 * sh3 - sh4 - 12.0 * sh2 * lnsh)
                / 3.0 / power_of<4>(1.0 - sh);
            complex<double> gb2 = 2.0 * ((-1.0 + sh) * (1.0 + sh * (-14.0 + sh * (-94.0 + (-14.0 + sh) * sh)))
                    + 60.0 * sh2 * (1.0 + sh) * lnsh)
                / 15.0 / power_of<5>(-1.0 + sh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j5_massless(const double & sh, const double & mB, const double & mu, const double & a1, const double & a2)
        {
            static const double pi = M_PI;

            double sh2 = sh * sh, sh3 = sh * sh2, sh4 = sh2 * sh2;
            double lnsh = std::log(sh), lnmbmu = std::log(mB / mu);
            //complex<double> dilogsh = dilog(sh);

            complex<double> asymp = 2.0 / 9.0 * (13.0 * (1.0 - sh2) + 2.0 * sh * (10.0 + 3.0 * sh) * lnsh - 6.0 * sh * lnsh * lnsh
                    - 6.0 * (1.0 - sh2 + 2.0 * sh * lnsh) * complex<double>(2.0 * lnmbmu, -pi))
                / power_of<3>(1.0 - sh);
            complex<double> gb1 = 2.0 / 3.0 * (7.0 + 39.0 * sh - 39.0 * sh2 - 7.0 * sh3 + 2.0 * sh * (10.0 + 19.0 * sh + sh2) * lnsh
                    - 6.0 * sh * (1.0 + sh) * lnsh * lnsh
                    - 2.0 * (1.0 + 9.0 * sh - 9.0 * sh2 - sh3 + 6.0 * sh * (1.0 + sh) * lnsh) * complex<double>(2.0 * lnmbmu, -pi))
                / power_of<4>(1.0 - sh);
            complex<double> gb2 = 1.0 / 3.0 * (17.0 + 296.0 * sh - 296.0 * sh3 - 17.0 * sh4 + 4.0 * sh * (20.0 + 96.0 * sh + 48.0 * sh2 + sh3) * lnsh
                    - 24.0 * sh * (1.0 + 3.0 * sh + sh2) * lnsh * lnsh
                    - 4.0 * (1.0 + 28.0 * sh - 28.0 * sh3 - sh4 + 12.0 * sh * (1.0 + 3.0 * sh + sh2) * lnsh) * complex<double>(2.0 * lnmbmu, -pi))
                / power_of<5>(1.0 - sh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        inline complex<double> j6_massless(const double & sh, const double & mB, const double & mu, const double & a1, const double & a2)
        {
            static const double pi = M_PI;

            double sh2 = sh * sh, sh3 = sh * sh2, sh4 = sh2 * sh2, sh5 = sh * sh4;
            double lnsh = std::log(sh), lnmbmu = std::log(mB / mu);
            //complex<double> dilogsh = dilog(sh);

            complex<double> asymp = 2.0 / 9.0 * (5.0 - 10.0 * sh + 7.0 * sh2 - 2.0 * sh3
                    - 2.0 * power_of<3>(1.0 - sh) * complex<double>(2.0 * lnmbmu, -pi)
                    + 2.0 * sh * (3.0 - 3.0 * sh + sh2) * lnsh)
                / power_of<3>(1.0 - sh);
            complex<double> gb1 = 1.0 / 9.0 * (3.0 + 10.0 * sh - 18.0 * sh2 + 6.0 * sh3 - sh4 + 12.0 * sh * lnsh)
                / power_of<4>(1.0 - sh);
            complex<double> gb2 = 1.0 / 45.0 * (6.0 + 125.0 * sh - 80.0 * sh2 - 60.0 * sh3 + 10.0 * sh4 - sh5 + 60.0 * sh * (1.0 + 2.0 * sh) * lnsh)
                / power_of<5>(1.0 - sh);

            return asymp + a1 * gb1 + a2 * gb2;
        }

        // We use the same regularising cut-off x ~= Lambda / m_B as in j7_zero as to ensure
        // a smooth transition B->K^*ll -> B->K^*gamma for s -> 0.
        // The relative error for j7  in the QCDF region 1 <= q^2 <= 6 is less than 25%.
        // Since j7 enters only via subleading terms, it amounts to a relative error of A_FB
        // in the SM of < 0.3%.
        inline double j7_massless(const double & sh, const double & x, const double & a1, const double & a2)
        {
            double lnsh = std::log(sh + x - sh * x), sh2 = sh * sh, sh3 = sh2 * sh, sh4 = sh2 * sh2;
            double x2 = x * x, x3 = x2 * x, x4 = x2 * x2;

            double asymp = 6.0 * (-1.0 + x + sh * (2.0 - x - 1.0 / (sh + x - sh * x)) - (1.0 + sh) * lnsh)
                / power_of<3>(1.0 - sh);
            double gb1 = 18.0 * (-3.0 * sh - lnsh * sh - 4.0 * lnsh * sh2 + 3.0 * sh3
                    - lnsh * sh3 - 2.0 * x - lnsh * x + 6.0 * sh * x - 3.0 * lnsh * sh * x
                    + 3.0 * lnsh * sh2 * x - 4.0 * sh3 * x + lnsh * sh3 * x + 3.0 * x2 - 6.0 * sh * x2
                    + 3.0 * sh2 * x2 - x3 + 3.0 * sh * x3 - 3.0 * sh2 * x3
                    + sh3 * x3)
                / (power_of<4>(1.0 - sh) * (sh + x - sh * x));
            double gb2 = 12.0 * (-11.0 * sh - 3.0 * lnsh * sh - 27.0 * sh2 - 27.0 * lnsh * sh2
                    + 27.0 * sh3 - 27.0 * lnsh * sh3 + 11.0 * sh4 - 3.0 * lnsh * sh4 - 8.0 * x
                    - 3.0 * lnsh * x + 8.0 * sh * x - 24.0 * lnsh * sh * x + 54.0 * sh2 * x
                    - 40.0 * sh3 * x + 24.0 * lnsh * sh3 * x - 14.0 * sh4 * x + 3.0 * lnsh * sh4 * x
                    + 18.0 * x2 - 27.0 * sh * x2 + 3.0 * sh2 * x2 + 3.0 * sh3 * x2 + 3.0 * sh4 * x2
                    - 15.0 * x3 + 50.0 * sh * x3 - 60.0 * sh2 * x3 + 30.0 * sh3 * x3 - 5.0 * sh4 * x3
                    + 5.0 * x4 - 20.0 * sh * x4 + 30.0 * sh2 * x4 - 20.0 * sh3 * x4 + 5.0 * sh4 * x4)
                / (power_of<5>(1.0 - sh) * (sh + x - sh * x));

            return asymp + a1 * gb1 + a2 * gb2;
        }
    }

    /* s = 0, case for B->V gamma */

    // bottom case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_bottom_case(const double & m_b, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        double mh = m_b / m_B;
        double eh = (1.0 + power_of<2>(m_V / m_B)) / 2.0;

        /*
         * J2 itself is divergent for s -> 0, however, it enters via s * J2.
         * J3 is divergent for s -> 0 but does not emerge in B->Vgamma decays.
         * Set both integrals to NaN.
         */

        // perpendicular amplitude
        results.j0_perp = impl::moment_inverse_ubar(a_1_perp, a_2_perp);
        results.j0bar_perp = impl::moment_inverse_ubar(-a_1_perp, a_2_perp);
        results.j1_perp = impl::j1_szero_bottom(mh, a_1_perp, a_2_perp);
        results.j2_perp = complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        results.j4_perp = impl::j4_szero_bottom(m_b, m_B, mu, a_1_perp, a_2_perp);
        results.j5_perp = impl::j5_szero_bottom(m_b, m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        results.j6_perp = impl::j6_szero_bottom(m_b, m_B, mu, a_1_parallel, a_2_parallel);
        results.j7_perp = impl::j7_szero(0.5 / m_B, a_1_perp, a_2_perp);

        // parallel amplitude
        results.j0_parallel = impl::moment_inverse_ubar(a_1_parallel, a_2_parallel);
        results.j1_parallel = impl::j1_szero_bottom(mh, a_1_parallel, a_2_parallel);
        results.j3_parallel = complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        results.j4_parallel = impl::j4_szero_bottom(m_b, m_B, mu, a_1_parallel, a_2_parallel);

        // composite results
        results.jtilde1_perp = 2.0 / eh * results.j1_perp;
        results.jtilde2_parallel = std::numeric_limits<double>::signaling_NaN();

        return results;
    }

    // charm case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_charm_case(const double & m_c, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        double mh = m_c / m_B;
        double eh = (1.0 + power_of<2>(m_V / m_B)) / 2.0;

        /*
         * J2 itself is divergent for s -> 0, however, it enters via s * J2.
         * J3 is divergent for s -> 0 but does not emerge in B->Vgamma decays.
         * Set both integrals to NaN.
         */

        // perpendicular amplitude
        results.j0_perp = impl::moment_inverse_ubar(a_1_perp, a_2_perp);
        results.j0bar_perp = impl::moment_inverse_ubar(-a_1_perp, a_2_perp);
        results.j1_perp = impl::j1_szero_charm(mh, a_1_perp, a_2_perp);
        results.j2_perp = complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        results.j4_perp = impl::j4_szero_charm(m_c, m_B, mu, a_1_perp, a_2_perp);
        results.j5_perp = impl::j5_szero_charm(m_c, m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        results.j6_perp = impl::j6_szero_charm(m_c, m_B, mu, a_1_parallel, a_2_parallel);
        results.j7_perp = impl::j7_szero(0.5 / m_B, a_1_perp, a_2_perp);

        // parallel amplitude
        results.j0_parallel = impl::moment_inverse_ubar(a_1_parallel, a_2_parallel);
        results.j1_parallel = impl::j1_szero_charm(mh, a_1_parallel, a_2_parallel);
        results.j3_parallel = complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        results.j4_parallel = impl::j4_szero_charm(m_c, m_B, mu, a_1_parallel, a_2_parallel);

        // composite results
        results.jtilde1_perp = 2.0 / eh * results.j1_perp;
        results.jtilde2_parallel = std::numeric_limits<double>::signaling_NaN();

        return results;
    }

    // massless case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::photon_massless_case(const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        double eh = (1.0 + power_of<2>(m_V / m_B)) / 2.0;

        /*
         * J2 itself is divergent for s -> 0, however, it enters via s * J2.
         * J3 is divergent for s -> 0 but does not emerge in B->Vgamma decays.
         * Set both integrals to NaN.
         */

        // perpendicular amplitude
        results.j0_perp = impl::moment_inverse_ubar(a_1_perp, a_2_perp);
        results.j0bar_perp = impl::moment_inverse_ubar(-a_1_perp, a_2_perp);
        results.j1_perp = results.j0_perp;
        results.j2_perp = complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        results.j4_perp = impl::j4_szero_massless(m_B, mu, a_1_perp, a_2_perp);
        results.j5_perp = impl::j5_szero_massless(m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        results.j6_perp = impl::j6_szero_massless(m_B, mu, a_1_parallel, a_2_parallel);
        results.j7_perp = impl::j7_szero(0.5 / m_B, a_1_perp, a_2_perp);

        // parallel amplitude
        results.j0_parallel = impl::moment_inverse_ubar(a_1_parallel, a_2_parallel);
        results.j1_parallel = results.j0_parallel;
        results.j3_parallel = complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        results.j4_parallel = impl::j4_szero_massless(m_B, mu, a_1_parallel, a_2_parallel);

        // composite results
        results.jtilde1_perp = 2.0 / eh * results.j1_perp;
        results.jtilde2_parallel = std::numeric_limits<double>::signaling_NaN();

        return results;
    }

    /* s > 0, case for B->V l^+ l^- */

    // bottom case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_bottom_case(const double & s, const double & m_b, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        double sh = s / m_B / m_B, mh = m_b / m_B;
        double eh = (1.0 + power_of<2>(m_V / m_B) - sh) / 2.0;

        impl::DileptonIntegralsBottom integrals(sh, mh, m_B, mu);

        // perpendicular amplitude
        results.j0_perp = impl::j0(sh, a_1_perp, a_2_perp);
        results.j0bar_perp = impl::j0bar(sh, a_1_perp, a_2_perp);
        results.j1_perp = integrals.j1(a_1_perp, a_2_perp);
        results.j2_perp = integrals.j2(a_1_perp, a_2_perp);
        results.j4_perp = integrals.j4(a_1_perp, a_2_perp);
        results.j5_perp = integrals.j5(a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        results.j6_perp = integrals.j6(a_1_parallel, a_2_parallel);
        results.j7_perp = impl::j7_massless(sh, 0.5 / m_B, a_1_perp, a_2_perp);

        // parallel amplitude
        results.j0_parallel = impl::j0(sh, a_1_parallel, a_2_parallel);
        results.j1_parallel = integrals.j1(a_1_parallel, a_2_parallel);
        results.j3_parallel = integrals.j3(a_1_parallel, a_2_parallel);
        results.j4_parallel = integrals.j4(a_1_parallel, a_2_parallel);

        // composite results
        results.jtilde1_perp = 2.0 / eh * results.j1_perp + sh * results.j2_perp / (eh * eh);
        results.jtilde2_parallel = 2.0 / eh * results.j1_parallel + results.j3_parallel / (eh * eh);

        return results;
    }

    // charm case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_charm_case(const double & s, const double & m_c, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        double sh = s / m_B / m_B, rho = 4.0 * m_c * m_c / s, mh = m_c / m_B;
        double eh = (1.0 + power_of<2>(m_V / m_B) - sh) / 2.0;

        if ((rho > 0) && (rho < 1.0))
            throw InternalError("QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_charm_case: charm mass too small, rho = " + stringify(rho) + ", m_c = " + stringify(m_c) + ", s = " + stringify(s));

        impl::DileptonIntegralsCharm integrals(sh, mh, m_B, mu);

        // perpendicular amplitude
        results.j0_perp = impl::j0(sh, a_1_perp, a_2_perp);
        results.j0bar_perp = impl::j0bar(sh, a_1_perp, a_2_perp);
        results.j1_perp = integrals.j1(a_1_perp, a_2_perp);
        results.j2_perp = integrals.j2(a_1_perp, a_2_perp);
        results.j4_perp = integrals.j4(a_1_perp, a_2_perp);
        results.j5_perp = integrals.j5(a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        results.j6_perp = integrals.j6(a_1_parallel, a_2_parallel);
        results.j7_perp = impl::j7_massless(sh, 0.5 / m_B, a_1_perp, a_2_perp);

        // parallel amplitude
        results.j0_parallel = impl::j0(sh, a_1_parallel, a_2_parallel);
        results.j1_parallel = integrals.j1(a_1_parallel, a_2_parallel);
        results.j3_parallel = integrals.j3(a_1_parallel, a_2_parallel);
        results.j4_parallel = integrals.j4(a_1_parallel, a_2_parallel);

        // composite results
        results.jtilde1_perp = 2.0 / eh * results.j1_perp + sh * results.j2_perp / (eh * eh);
        results.jtilde2_parallel = 2.0 / eh * results.j1_parallel + results.j3_parallel / (eh * eh);

        return results;
    }

    // massless case
    template <>
    QCDFIntegrals<BToKstarDilepton>
    QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_massless_case(const double & s, const double & m_B, const double & m_V, const double & mu,
                    const double & a_1_perp, const double & a_2_perp,
                    const double & a_1_parallel, const double & a_2_parallel)
    {
        QCDFIntegrals<BToKstarDilepton> results;
        double sh = s / m_B / m_B;
        double eh = (1.0 + power_of<2>(m_V / m_B) - sh) / 2.0;

        // perpendicular amplitude
        results.j0_perp = impl::j0(sh, a_1_perp, a_2_perp);
        results.j0bar_perp = impl::j0bar(sh, a_1_perp, a_2_perp);
        results.j1_perp = impl::moment_inverse_ubar(a_1_perp, a_2_perp);
        results.j2_perp = impl::j2_massless(sh, a_1_perp, a_2_perp);
        results.j4_perp = impl::j4_massless(sh, m_B, mu, a_1_perp, a_2_perp);
        results.j5_perp = impl::j5_massless(sh, m_B, mu, a_1_perp, a_2_perp);
        // This integral arises in perpendicular amplitudes, but depends on parallel Gegenbauer moments!
        results.j6_perp = impl::j6_massless(sh, m_B, mu, a_1_parallel, a_2_parallel);
        results.j7_perp = impl::j7_massless(sh, 0.5 / m_B, a_1_perp, a_2_perp);

        // parallel amplitude
        results.j0_parallel = impl::j0(sh, a_1_parallel, a_2_parallel);
        results.j1_parallel = impl::moment_inverse_ubar(a_1_parallel, a_2_parallel);
        results.j3_parallel = impl::j3_massless(sh, a_1_parallel, a_2_parallel);
        results.j4_parallel = impl::j4_massless(sh, m_B, mu, a_1_parallel, a_2_parallel);

        // composite results
        results.jtilde1_perp = 2.0 / eh * results.j1_perp + sh * results.j2_perp / (eh * eh);
        results.jtilde2_parallel = 2.0 / eh * results.j1_parallel + results.j3_parallel / (eh * eh);

        return results;
    }
}
