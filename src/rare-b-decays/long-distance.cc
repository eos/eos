/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/long-distance.hh>

#include <vector>

namespace wf
{
    complex<double>
    LongDistance::g_had_ccbar(const double & s, const double & m_c)
    {
        // cf. [KS1996], Eqs. (3.3) and (3.4), p. 5
        static const double alpha = 1.0 / 133.0;
        static const double m_B = 5.279;
        static const double m_D = 1.865;
        static const double s_0_hat = 4.0 * m_D * m_D / m_B / m_B;
        // cf. [KS1996], Eq. (A2), p. 9
        static const double s_1_hat = 0.60;
        static const double s_2_hat = 0.69;

        // We us a universal fudge factor, kappa_V = kappa, adjusted
        // so that C_0({C_i}) * kappa ~= 0.72. Using C_0^NLL = 0.61.
        static const double fudge = 1.2;

        static const std::vector<double> m_cc_hat
        {
            /* Masses in GeV */
            3.0969 / m_B,
            3.6861 / m_B,
            3.771 / m_B,
            4.039 / m_B,
            4.153 / m_B,
            4.421 / m_B,
        };
        static const std::vector<double> gamma_cc_hat
        {
            /* Decay widths in GeV */
            9.34e-5 / m_B,
            3.37e-4 / m_B,
            2.30e-2 / m_B,
            8.00e-2 / m_B,
            1.03e-1 / m_B,
            6.20e-2 / m_B,
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

        double s_hat = s / m_B / m_B, alpha2 = alpha * alpha;

        double imag_res = 0.0, real_res = 0.0;
        for (auto m = m_cc_hat.cbegin(), g = gamma_cc_hat.cbegin(), b = br_cc.cbegin() ; m_cc_hat.cend() != m ; ++m, ++g, ++b)
        {
            double aa = (9.0 / alpha2) * (*b) * (*g) * (*g);
            double bb = (*m) * (*m);
            double cc = (*m) * (*g);

            imag_res += aa * s_hat / ((s_hat - bb) * (s_hat - bb) + cc * cc);
            real_res += aa * s_hat / (6.0 * cc) * (
                    (bb - s_hat) * (M_PI + 2.0 * atan((bb - s_0_hat) / cc))
                    + cc * log((s_0_hat - s_hat) * (s_0_hat - s_hat) / ((bb - s_0_hat) * (bb - s_0_hat) + cc * cc)));
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

        double real_cont = s_hat / 3.0 * (
                0.0 // Contribution from s_hat_prime < 0.60
                + log(std::abs((s_2_hat - s_hat) / (s_1_hat - s_hat))) // Contribution from 0.60 < s_hat_prime < 0.69
                - log(std::abs((1.0 - s_hat / s_2_hat))) / s_hat); // Contribution from 0.69 < s_hat_prime < 1

        return complex<double>(real_cont + fudge * real_res, M_PI / 3.0 * (fudge * imag_res + imag_cont))
            -8.0 / 9.0 * log(m_c / m_B) - 4.0 / 9.0;
    }
}

