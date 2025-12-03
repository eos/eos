/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#include <eos/nonlocal-form-factors/long-distance.hh>

#include <vector>

namespace eos
{
    complex<double>
    LongDistance::g_had_ccbar(const double & s, const double & m_c)
    {
        // cf. [KS:1996A], Eqs. (3.3) and (3.4), p. 5
        static const double alpha = 1.0 / 133.0;
        // cf. [KS:1996A], Appendix: Input Parameters, p. 9
        static const double m_b = 4.8;
        static const double m_D = 1.865;
        static const double s_0_hat = 4.0 * m_D * m_D / m_b / m_b;
        // cf. [KS:1996A], Eq. (A2), p. 9
        static const double s_1_hat = 0.60;
        static const double s_2_hat = 0.69;

        // We us a universal fudge factor, kappa_V = kappa, adjusted
        // so that C_0({C_i}) * kappa ~= 0.72. Using C_0^NLL = 0.61.
        // This yields fudge O(1), so we simply set:
        static const double fudge = 1.0;

        static const std::vector<double> m_cc_hat
        {
            /* Masses in GeV */
            3.0969 / m_b,
            3.6861 / m_b,
            3.771 / m_b,
            4.039 / m_b,
            4.153 / m_b,
            4.421 / m_b,
        };
        static const std::vector<double> gamma_cc_hat
        {
            /* Decay widths in GeV */
            9.34e-5 / m_b,
            3.37e-4 / m_b,
            2.30e-2 / m_b,
            8.00e-2 / m_b,
            1.03e-1 / m_b,
            6.20e-2 / m_b,
        };
        static const std::vector<double> br_cc
        {
            5.935e-2,
            7.325e-3,
            1.050e-5,
            1.070e-5,
            8.100e-6,
            9.400e-6,
        };

        double s_hat = s / m_b / m_b, alpha2 = alpha * alpha;

        double imag_res = 0.0, real_res = 0.0;
        for (auto m = m_cc_hat.cbegin(), g = gamma_cc_hat.cbegin(), b = br_cc.cbegin() ; m_cc_hat.cend() != m ; ++m, ++g, ++b)
        {
            double aa = (9.0 / alpha2) * (*b) * (*g) * (*g);
            double bb = (*m) * (*m);
            double cc = (*m) * (*g);

            imag_res += aa * s_hat / ((s_hat - bb) * (s_hat - bb) + cc * cc);
            real_res += aa / (2.0 * cc) * (
                        (s_hat - bb) * (M_PI + 2.0 * std::atan((bb - s_0_hat) / cc))
                        - cc * std::log((s_0_hat - s_hat) * (s_0_hat - s_hat) / ((bb - s_0_hat) * (bb - s_0_hat) + cc * cc))
                    ) / ((bb - s_hat) * (bb - s_hat) + cc * cc);
        }

        double imag_cont = 0.0;
        if ((s_hat >= 0.60) && (s_hat < 0.69))
        {
            imag_cont = 11.33 * s_hat - 6.80;
        }
        else if (s_hat >= 0.69)
        {
            imag_cont = 1.02;
        }

        double real_cont = (
                    0.571896
                    + 34.0 / 3.0 * ((s_1_hat - s_hat) * std::log(std::abs(s_1_hat - s_hat)) - (s_2_hat - s_hat) * std::log(std::abs(s_2_hat - s_hat)))
                ) / s_hat;

        return complex<double>(s_hat * (real_cont + fudge * real_res) / 3.0, M_PI / 3.0 * (fudge * imag_res + imag_cont))
            -8.0 / 9.0 * log(m_c / m_b) - 4.0 / 9.0;
    }
}
