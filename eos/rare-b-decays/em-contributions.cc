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

#include <eos/maths/power-of.hh>
#include <eos/rare-b-decays/em-contributions.hh>

#include <cmath>

#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    // cf. [HLMW:2005A], Eq. (94), p. 23
    double
    EMContributions::omegaem_99(const double & s_hat, const double & log_m_l_hat)
    {
        double li2 = gsl_sf_dilog(s_hat);
        double ln = std::log(s_hat), ln1 = std::log(1.0 - s_hat);
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;

        return -2.0 * log_m_l_hat * (
                - (1.0 + 4.0 * s_hat - 8.0 * s_hat2) / (6.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
                + ln1
                - (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3) / (2.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)) * ln
                )
            - 1.0 / 9.0 * li2
            + 4.0 * power_of<2>(M_PI) / 27.0
            - (37.0 - 3.0 * s_hat - 6.0 * s_hat2) / (72.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
            - (41.0 + 76.0 * s_hat) / (36.0 * (1.0 + 2.0 * s_hat)) * ln1
            + (14.0 * s_hat3 - 17.0 * s_hat2 - 10.0 * s_hat + 6.0) / (18.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)) * ln
            + 17.0 / 18.0 * ln1 * ln
            - (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3) / (2.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)) * ln * ln;
    }

    // cf. [HLMW:2005A], Eq. (100), p. 24
    double
    EMContributions::omegaem_1010(const double & s_hat, const double & log_m_l_hat)
    {
        double ln = std::log(s_hat), ln1 = std::log(1.0 - s_hat);
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;

        return -2.0 * log_m_l_hat * (
                - (1.0 + 4.0 * s_hat - 8.0 * s_hat2) / (6.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
                + ln1
                - (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3) / (2.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)) * ln
                );
    }

    // cf. [HLMW:2005A], Eq. (101), p. 25
    double
    EMContributions::omegaem_77(const double & s_hat, const double & log_m_l_hat)
    {
        double ln = std::log(s_hat), ln1 = std::log(1.0 - s_hat);
        double s_hat2 = s_hat * s_hat;

        return -2.0 * log_m_l_hat * (
                s_hat / (2.0 * (1.0 - s_hat) * (2.0 + s_hat))
                + ln1
                - s_hat * (-3.0 + 2.0 * s_hat2) / (2.0 * power_of<2>(1.0 - s_hat) * (2.0 + s_hat)) * ln
                );
    }

    // cf. [HLMW:2005A], Eq. (102), p. 25
    double
    EMContributions::omegaem_79(const double & s_hat, const double & log_m_l_hat)
    {
        double ln = std::log(s_hat), ln1 = std::log(1.0 - s_hat);
        double s_hat2 = s_hat * s_hat;

        return -2.0 * log_m_l_hat * (
                -1.0 / (2.0 * (1.0 - s_hat))
                + ln1
                + (-1.0 + 2.0 * s_hat - 2.0 * s_hat2) / (2.0 * power_of<2>(1.0 - s_hat)) * ln
                );
    }

    // cf. [HLMW:2005A], Eq. (103), p. 25
    complex<double>
    EMContributions::omegaem_29(const double & s_hat, const double & log_m_l_hat, const double & mu)
    {
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
        double sigma_1 = 23.787 - 120.948 * s_hat + 365.373 * s_hat2 - 584.206 * s_hat3;
        double sigma_1_I = 1.653 + 6.009 * s_hat - 17.080 * s_hat2 + 115.880 * s_hat3;

        return -2.0 * log_m_l_hat * complex<double>(sigma_1, sigma_1_I) / (8.0 * (1.0 - s_hat) * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
            + 16.0 / 9.0 * omegaem_1010(s_hat, log_m_l_hat) * log(mu / 5.0);
    }

    // cf. [HLMW:2005A], Eq. (104), p. 25
    double
    EMContributions::omegaem_22(const double & s_hat, const double & log_m_l_hat, const double & mu)
    {
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat, s_hat4 = s_hat2 * s_hat2;
        double sigma_1 = 23.787 - 120.948 * s_hat + 365.373 * s_hat2 - 584.206 * s_hat3;
        double sigma_2 = 11.488 - 36.987 * s_hat + 255.330 * s_hat2 - 812.388 * s_hat3 + 1011.791 * s_hat4;

        return -2.0 * log_m_l_hat * (
                sigma_2 / (8.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat))
                + sigma_1 / (9.0 * power_of<2>(1.0 - s_hat) * (1.0 + 2.0 * s_hat)) * log(mu / 5.0)
                )
            + 64.0 / 81.0 * omegaem_1010(s_hat, log_m_l_hat) * power_of<2>(log(mu / 5.0));
    }

    // cf. [HLMW:2005A], Eq. (105), p. 25
    complex<double>
    EMContributions::omegaem_27(const double & s_hat, const double & log_m_l_hat, const double & mu)
    {
        double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
        double sigma_3 = 109.311 - 846.039 * s_hat + 2890.115 * s_hat2 - 4179.072 * s_hat3;
        double sigma_3_I = 4.606 + 17.650 * s_hat - 53.244 * s_hat2 + 348.069 * s_hat3;

        return -2.0 * log_m_l_hat * complex<double>(sigma_3, sigma_3_I) / (96.0 * (1.0 - s_hat) * (1.0 - s_hat))
            + 8.0 / 9.0 * omegaem_79(s_hat, log_m_l_hat) * log(mu / 5.0);
    }
}
