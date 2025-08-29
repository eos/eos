/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2016, 2017 Danny van Dyk
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
#include <eos/nonlocal-form-factors/hard-scattering.hh>

#include <gsl/gsl_sf_dilog.h>

#include <limits>

namespace eos
{
    inline double lcda_tw2(const double & u, const double & a_1, const double & a_2)
    {
        return 6.0 * u * (1.0 - u) * (1.0 + a_1 * 3.0 * (2.0 * u - 1.0) + a_2 * 3.0 / 2.0 * (5.0 * power_of<2>(2.0 * u - 1.0) - 1.0));
    }

    static complex<double> z_ratio(const double & s, const double mq2)
    {
        static const complex<double> i{0.0, 1.0};

        if (s >= 4.0 * mq2)
        {
            const auto root = sqrt(s / (s - 4.0 * mq2));
            return (root + 1.0) / (root - 1.0) + i * sqrt(std::numeric_limits<double>::epsilon());
        }
        else if (s > 0)
        {
            const auto root = sqrt((4.0 * mq2 - s) / s);
            return (i - root) / (i + root) + i * sqrt(std::numeric_limits<double>::epsilon());
        }
        else if (s == 0)
        {
            return -1.0 + i * sqrt(std::numeric_limits<double>::epsilon());
        }
        else
        {
            const auto root = sqrt(-s / (-s + 4.0 * mq2));
            return (root + 1.0) / (root - 1.0) + i * sqrt(std::numeric_limits<double>::epsilon());
        }
    }

    complex<double>
    HardScattering::I1(const double & q2, const double & u, const double & m_q, const double & m_B)
    {
        if (m_q == 0.0)
            return complex<double>(1.0, 0.0);

        const auto ubar = 1.0 - u;
        const auto s    = ubar * m_B * m_B + u * q2;

        const auto x_ratio = z_ratio(s,  m_q * m_q);
        const auto y_ratio = z_ratio(q2, m_q * m_q);

        const auto IauxX = -0.5 * M_PI * M_PI - 0.5 * pow(log(x_ratio), 2) + log(x_ratio) * log(-x_ratio);
        const auto IauxY = -0.5 * M_PI * M_PI - 0.5 * pow(log(y_ratio), 2) + log(y_ratio) * log(-y_ratio);

        return 1.0 + 2.0 * m_q * m_q / (ubar * (m_B * m_B - q2)) * (IauxX - IauxY);
    }
/*
    complex<double>
    HardScattering::I1(const double & s, const double & u, const double & m_q, const double & m_B)
    {
        if (m_q == 0.0)
            return complex<double>(1.0, 0.0);

        gsl_sf_result res_re, res_im;

        double ubar = 1.0 - u;
        double m_q2 = m_q * m_q;
        double m_B2 = m_B * m_B;

        double a, a2, sign;
        complex<double> dilogArg, dilog1, dilog2;
        complex<double> LxpLxm, LypLym;

        if (1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) > 0)
        {
            a = (1 - std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))))
              / (1 + std::sqrt(1.0 - 4.0 * m_q2 / (m_B2 - u * (m_B2 - s))));
            LxpLxm = -M_PI* M_PI / 3.0
                + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                + gsl_sf_dilog(-a) + gsl_sf_dilog(-1.0/a);
        }
        else
        {
            a  = sqrt(4.0 * m_q2 / (m_B2 - u * (m_B2 - s)) - 1);
            a2 = a * a;

            if (a2 - 1 > 0)
                sign = +1.0;
            else
                sign = -1.0;

            dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), -2.0 * a / (a2 + 1.0));
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog1 = complex<double>( res_re.val, res_im.val);

            dilogArg = complex<double>((a2 - 1.0)/(a2 + 1.0), +2.0 * a / (a2 + 1.0));
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog2 = complex<double>(res_re.val, res_im.val);

            LxpLxm = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                + dilog1 + dilog2;
        }

        if (1.0 - 4.0 * m_q2 / s > 0)
        {
            a = (1.0 - std::sqrt(1.0 - 4.0 * m_q2 / s))
              / (1.0 + std::sqrt(1.0 - 4.0 * m_q2 / s));
            LypLym = -1.0 / 3.0 * M_PI * M_PI + std::log(a) * (std::log(a) + complex<double>(0.0, M_PI))
                + gsl_sf_dilog(-a) + gsl_sf_dilog(-1./a);
        }
        else
        {
            a  = std::sqrt(4.0 * m_q2 / s - 1.0);
            a2 = a * a;
            if (a2 - 1.0 > 0)
                sign = +1.0;
            else
                sign = -1.0;

            dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), -2.0 * a / (a2 + 1.0));
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog1 = complex<double>(res_re.val, res_im.val);

            dilogArg = complex<double>((a2 - 1.0) / (a2 + 1.0), +2.0 * a / (a2 + 1.0));
            gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
            dilog2 = complex<double>(res_re.val, res_im.val);

            LypLym = -1.0 / 3.0 * M_PI * M_PI - std::atan(2.0 * a / (a2 - 1.0)) * (std::atan(2.0 * a / (a2 - 1.0)) - M_PI * sign)
                + dilog1 + dilog2;
        }

        return 1.0 + 2.0 * m_q2 / (ubar * (m_B2 - s)) * (LxpLxm - LypLym);
    }
*/

    complex<double>
    HardScattering::t_perp_s0(const double & u, const double & m_q, const double & m_B)
    {
        const double ub = 1 - u;
        const double m2 = power_of<2>(m_q);
        const double m_B2 = power_of<2>(m_B);
        double a, a2, sign, sq;
        complex<double> dilogArg, dilog1, dilog2;
        complex<double> LxpLxm;

        gsl_sf_result res_re, res_im;

        if (m_q > 0) { // m != 0
            if (1 - 4 * m2/(m_B2 - u * m_B2) > 0) {
                sq = std::sqrt(1 - 4 * m2/(m_B2 - u * m_B2));
                a  = (1 - sq)/(1 + sq);
                LxpLxm = -1./3 * power_of<2>(M_PI) + std::log(a) * (std::log(a) + complex<double>(0, M_PI)) +
                        gsl_sf_dilog(-a) + gsl_sf_dilog(-1./a);
            } else {
                a2 = 4. * m2/(m_B2 - u * m_B2) - 1;
                a  = std::sqrt(a2);
                if (a2 - 1 > 0)
                    sign = +1.;
                else
                    sign = -1.;
                dilogArg = complex<double>((a2 - 1)/(a2 + 1), -2 * a/(a2 + 1));
                gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog1 = complex<double>(res_re.val, res_im.val);
                dilogArg = complex<double>((a2 - 1)/(a2 + 1), +2. * a/(a2 + 1));
                gsl_sf_complex_dilog_e(abs(dilogArg), arg(dilogArg), &res_re, &res_im);
                dilog2 = complex<double>(res_re.val, res_im.val);
                sq = std::atan(2 *a/(a2 - 1));
                LxpLxm = -1./3 * power_of<2>(M_PI) - sq * (sq - M_PI * sign) + dilog1 + dilog2;
            }
            return 4.0 / ub * (1.0 + 2.0 * m2 / ub / m_B2 * LxpLxm);
        } else {     // m == 0
            return complex<double>(4.0 / ub, 0);
        }
    }

    complex<double>
    HardScattering::t_perp(const double & s, const double & u, const double & m_q,
                           const double & m_B, const double & m_M)
    {
        if (s == 0)
            return t_perp_s0(u, m_q, m_B);

        const double ub = 1.0 - u;
        const double E = 1.0 / (2.0 * m_B) * (power_of<2>(m_B)+ power_of<2>(m_M) - s);

        if (m_q > 0)
            return 2 * m_B / (ub * E) * I1(s, u, m_q, m_B) +
                    s / (ub * ub * E * E) * (CharmLoops::B0(ub * power_of<2>(m_B)+ u * s, m_q) - CharmLoops::B0(s, m_q));
        return 2 * m_B / (ub * E) * I1(s, u, m_q, m_B);
    }

    complex<double>
    HardScattering::t_par(const double & s, const double & u, const double & m_q,
                          const double & m_B, const double & m_M)
    {
        const double ub = 1.0 - u;
        const double E = 1. / (2 * m_B) * (power_of<2>(m_B)+ power_of<2>(m_M) - s);

        if (m_q > 0)
            return 2 * m_B / (ub * E) * I1(s, u, m_q, m_B) +
                   (ub * power_of<2>(m_B) + u * s) / (ub * ub * E * E)
                   * (CharmLoops::B0(ub * power_of<2>(m_B)+ u * s, m_q) - CharmLoops::B0(s, m_q));
        return 2 * m_B / (ub * E) * I1(s, u, m_q, m_B);
    }

    double
    HardScattering::j0(const double & s, const double & u, const double & m_B, const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u;
        const double s_hat = s / power_of<2>(m_B);

        return lcda_tw2(u, a_1, a_2) / (ubar + u * s_hat);
    }

    complex<double>
    HardScattering::j1(const double & s, const double & u, const double & m_q, const double & m_B,
            const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u;

        return lcda_tw2(u, a_1, a_2) / ubar * I1(s, u, m_q, m_B);
    }

    complex<double>
    HardScattering::j2(const double & s, const double & u, const double & m_q, const double & m_B,
            const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u, ubar2 = ubar * ubar;
        const double m_B2 = m_B * m_B;

        return lcda_tw2(u, a_1, a_2) * (CharmLoops::B0(ubar * m_B2 + u * s, m_q) - CharmLoops::B0(s, m_q)) / ubar2;
    }

    complex<double>
    HardScattering::j2_massless(const double & s, const double & u, const double & m_B,
            const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u, ubar2 = ubar * ubar;
        const double s_hat = s / power_of<2>(m_B);

        return lcda_tw2(u, a_1, a_2) * std::log(s_hat / (ubar + u * s_hat)) / ubar2;
    }

    complex<double>
    HardScattering::j3(const double & s, const double & u, const double & m_q, const double & m_B,
            const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u;
        const double s_hat = s / (m_B * m_B);

        return HardScattering::j2(s, u, m_q, m_B, a_1, a_2) * (ubar + u * s_hat);
    }

    complex<double>
    HardScattering::j3_massless(const double & s, const double & u, const double & m_B,
            const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u;
        const double s_hat = s / power_of<2>(m_B);

        return HardScattering::j2_massless(s, u, m_B, a_1, a_2) * (ubar + u * s_hat);
    }

    complex<double>
    HardScattering::j4(const double & s, const double & u, const double & m_q, const double & m_B,
            const double & mu, const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u;

        return lcda_tw2(u, a_1, a_2) * CharmLoops::h(mu, ubar * power_of<2>(m_B) + u * s, m_q);
    }

    complex<double>
    HardScattering::j5(const double & s, const double & u, const double & m_q, const double & m_B,
            const double & mu, const double & a_1, const double & a_2)
    {
        const double ubar = 1.0 - u;
        const double s_hat = s / power_of<2>(m_B);

        return lcda_tw2(u, a_1, a_2) / (ubar + u * s_hat) * CharmLoops::h(mu, ubar * power_of<2>(m_B) + u * s, m_q);
    }

    complex<double>
    HardScattering::j6(const double & s, const double & u, const double & m_q, const double & m_B,
            const double & mu, const double & a_1, const double & a_2)
    {
        // cf. [BFS:2004], eq. (52): this integral does not involve the lcda
        // itself, but its first inverse partial moment, as the weight function.
        const double weight = power_of<2>(u) * (3.0 + a_1 * (-9.0 + 12.0 * u) +
                a_2 * (18.0 - 60.0 * u + 45.0 * power_of<2>(u)));
        const double ubar = 1.0 - u;

        return weight * CharmLoops::h(mu, ubar * power_of<2>(m_B) + u * s, m_q);
    }

    double
    HardScattering::j7(const double & s, const double & u, const double & m_B, const double & a_1,
            const double & a_2)
    {
        const double ubar = 1.0 - u;

        return lcda_tw2(u, a_1, a_2) / power_of<2>(ubar + u * s / power_of<2>(m_B));
    }
}
