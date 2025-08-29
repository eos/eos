/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2013, 2015 Danny van Dyk
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

#include <eos/maths/integrate.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>
#include <eos/rare-b-decays/bremsstrahlung.hh>

#include <cmath>
#include <limits>
#include <functional>

namespace eos
{
    // cf. [AAGW:2002], Eq. (30), p. 12
    complex<double>
    Bremsstrahlung::G_m1(const double & t)
    {
        complex<double> result;

        if (t < 4)
        {
            double x = atan(sqrt((4.0 - t) / t));

            result = 2.0 * M_PI * x - power_of<2>(M_PI) / 2 - 2.0 * power_of<2>(x);
        }
        else
        {
            double x = log((sqrt(t) + sqrt(t - 4.0)) / 2.0);

            result = complex<double>(-power_of<2>(M_PI) / 2.0 + 2.0 * power_of<2>(x),
                    -2.0 * M_PI * x);
        }

        return result;
    }

    // cf. [AAGW:2002], Eq. (31), p. 12
    complex<double>
    Bremsstrahlung::G_0(const double & t)
    {
        complex<double> result;

        if (t < 4)
        {
            double x = sqrt((4.0 - t) / t);

            result = M_PI * x - 2 - 2 * x * atan(x);
        }
        else
        {
            double x = sqrt((t - 4.0) / t);
            double y = log((sqrt(t) + sqrt(t - 4.0)) / 2.0);

            result = complex<double>(-2.0 + 2.0 * x * y,
                    -M_PI * x);
        }

        return result;
    }

    // cf. [AAGW:2002], Eq. (28), p. 11
    complex<double>
    Bremsstrahlung::Deltai_23(const double & s_hat, const double & w, const double & z)
    {
        return -2.0 + 4.0 / (w - s_hat) * (
                z * (Bremsstrahlung::G_m1(s_hat / z) - Bremsstrahlung::G_m1(w / z))
                - s_hat / 2.0 * (Bremsstrahlung::G_0(s_hat / z) - Bremsstrahlung::G_0(w / z)));
    }

    // cf. [AAGW:2002], Eq. (29), p. 11
    complex<double>
    Bremsstrahlung::Deltai_27(const double & s_hat, const double & w, const double & z)
    {
        return 2.0 * (Bremsstrahlung::G_0(s_hat / z) - Bremsstrahlung::G_0(w / z));
    }

    // cf. [AAGW:2002], Eq. (23), p. 10
    complex<double>
    Bremsstrahlung::tau_22(const double & s_hat, const double & w, const double & z)
    {
        double s_hat2 = s_hat * s_hat, w2 = w * w, w3 = w2 * w;

        return 8.0 / 27.0 * (w - s_hat) * power_of<2>(1 - w) / s_hat / w3 * (
                (3.0 * w2 + 2 * s_hat2 * (2.0 + w) - s_hat * w * (5 - 2.0 * w)) * norm(Deltai_23(s_hat, w, z))
                + (2.0 * s_hat2 * (2.0 + w) + s_hat * w * (1.0 + 2.0 * w)) * norm(Deltai_27(s_hat, w, z))
                + 4.0 * s_hat * (w * (1.0 - w) - s_hat * (2.0 + w)) * real(Deltai_23(s_hat, w, z) * conj(Deltai_27(s_hat, w, z))));
    }

    // cf. [AAGW:2002], Eq. (24), p. 10
    complex<double>
    Bremsstrahlung::tau_27(const double & s_hat, const double & w, const double & z)
    {
        double s_hat2 = s_hat * s_hat, w2 = w * w;

        return 8.0 / 3.0 / (s_hat * w) * (
                ((1.0 - w) * (4.0 * s_hat2 - s_hat * w + w2) + s_hat * w * (4.0 + s_hat - w) * log(w)) * Deltai_23(s_hat, w, z)
                - (4.0 * s_hat2 * (1.0 - w) + s_hat * w * (4.0 + s_hat - w) * log(w)) * Deltai_27(s_hat, w, z));
    }

    // cf. [AAGW:2002], Eq. (25), p. 10
    complex<double>
    Bremsstrahlung::tau_28(const double & s_hat, const double & w, const double & z)
    {
        double w2 = w * w;
        double x = s_hat / (1.0 + s_hat - w) / (w2 + s_hat * (1.0 - w));

        return 8.0 / 9.0 / (s_hat * w * (w - s_hat)) * (
                (power_of<2>(w - s_hat) * (2.0 * s_hat - w) * (1.0 - w)) * Deltai_23(s_hat, w, z)
                - (2.0 * s_hat * power_of<2>(w - s_hat) * (1.0 - w)) * Deltai_27(s_hat, w, z)
                + s_hat * w * ((1.0 + 2.0 * s_hat - 2.0 * w) * Deltai_23(s_hat, w, z)
                    - 2.0 * (1.0 + s_hat - w) * Deltai_27(s_hat, w, z)) * log(x));
    }

    // cf. [AAGW:2002], Eq. (24), p. 10
    complex<double>
    Bremsstrahlung::tau_29(const double & s_hat, const double & w, const double & z)
    {
        return 4.0 / 3.0 / w * (
                (2.0 * s_hat * (1.0 - w) * (s_hat + w) + 4.0 * s_hat * w * log(w)) * Deltai_23(s_hat, w, z)
                - (2.0 * s_hat * (1.0 - w) * (s_hat + w) + w * (3.0 * s_hat + w) * log(w)) * Deltai_27(s_hat, w, z));
    }

    // cf. [AAGW:2002], Eq. (15), p. 8
    double
    Bremsstrahlung::tau_78(const double & s_hat)
    {
        static const double pi = M_PI, pi2 = pi * pi;
        double ln_s_hat = std::log(s_hat), sqrt_s_hat = std::sqrt(s_hat), sqrt_4_m_s_hat = std::sqrt(4.0 - s_hat);
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
        double atan1 = std::atan((2.0 - 4.0 * s_hat + s_hat2) / ((2.0 - s_hat) * sqrt_s_hat * sqrt_4_m_s_hat));
        double atan2 = std::atan(sqrt_s_hat * sqrt_4_m_s_hat / (2.0 - s_hat));
        double atan3 = std::atan(sqrt_4_m_s_hat / sqrt_s_hat);
        double reli2 = real(dilog(complex<double>(s_hat / 2.0, -0.5 * sqrt_s_hat * sqrt_4_m_s_hat)));

        return 8.0 / (9.0 * s_hat) * (
                25.0 - 2.0 * pi2 - 27.0 * s_hat + 3.0 * s_hat2 - s_hat3 + 12.0 * (s_hat + s_hat2) * ln_s_hat
                + 6.0 * power_of<2>(pi / 2.0 - atan1) - 24.0 * reli2
                - 12.0 * ((1.0 - s_hat) * sqrt_s_hat * sqrt_4_m_s_hat - 2.0 * atan2) * (atan3 - atan2));
    }

    // cf. [AAGW:2002], Eq. (16), p. 8
    double
    Bremsstrahlung::tau_88(const double & s_hat)
    {
        static const double pi = M_PI, pi2 = pi * pi;
        double ln_s_hat = std::log(s_hat), sqrt_s_hat = std::sqrt(s_hat), sqrt_4_m_s_hat = std::sqrt(4.0 - s_hat);
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
        double atan1 = std::atan(sqrt_4_m_s_hat / sqrt_s_hat);
        double atan2 = std::atan(sqrt_s_hat * sqrt_4_m_s_hat / (2.0 - s_hat));
        double reli1 = real(dilog(1.0 - s_hat));
        double reli2 = real(dilog(complex<double>((3.0 - s_hat) / 2.0, (1.0 - s_hat) * sqrt_4_m_s_hat / (2.0 * sqrt_s_hat))));

        return 4.0 / (27.0 * s_hat) * (
                -8.0 * pi2 + (1.0 - s_hat) * (77.0 - s_hat - 4.0 * s_hat2) - 24.0 * reli1
                + 3.0 * (10.0 - 4.0 * s_hat - 9.0 * s_hat2 + 8.0 * std::log(sqrt_s_hat / (1.0 - s_hat))) * ln_s_hat
                + 48.00 * reli2
                - 6.0 * ((20.0 * s_hat + 10.0 * s_hat2 - 3.0 * s_hat3) / (sqrt_s_hat * sqrt_4_m_s_hat) - 8.0 * pi + 8.0 * atan1)
                    * (atan1 - atan2));
    }

    // cf. [AAGW:2002], Eq. (17), p. 8
    double
    Bremsstrahlung::tau_89(const double & s_hat)
    {
        double ln_s_hat = std::log(s_hat), sqrt_s_hat = std::sqrt(s_hat), sqrt_4_m_s_hat = std::sqrt(4.0 - s_hat);
        double s_hat2 = s_hat * s_hat;
        double arctan1 = std::atan(sqrt_s_hat * sqrt_4_m_s_hat / (2.0 - s_hat));
        double arctan2 = std::atan(sqrt_4_m_s_hat / sqrt_s_hat);
        double reli1 = real(dilog(complex<double>(s_hat / 2.0, sqrt_s_hat * sqrt_4_m_s_hat / 2.0)));
        double reli2 = real(dilog(complex<double>((-2.0 + s_hat * (4.0 - s_hat)) / 2.0, (2.0 - s_hat) * sqrt_s_hat * sqrt_4_m_s_hat / 2.0)));

        return 2.0 / 3.0 * (
                s_hat * (4.0 - s_hat) - 3.0 - 4.0 * ln_s_hat * (1.0 - s_hat - s_hat2)
                - 8.0 * (reli1 - reli2) + 4.0 * (s_hat2 * sqrt_4_m_s_hat / sqrt_s_hat + 2.0 * arctan1) * (arctan2 - arctan1));
    }

    // Integrals of tau_2x from w = s_hat to w = 1
    complex<double>
    Bremsstrahlung::itau_22(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate1D(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_22, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }

    complex<double>
    Bremsstrahlung::itau_27(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate1D(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_27, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }

    complex<double>
    Bremsstrahlung::itau_28(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate1D(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_28, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }

    complex<double>
    Bremsstrahlung::itau_29(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate1D(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_29, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }
}
