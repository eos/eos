/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014 Danny van Dyk
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
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

#include <algorithm>
#include <cmath>

namespace eos
{
    const double QCD::casimir_f = 4.0 / 3.0;

    /* 6 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 6
    const QCD::BetaFunction QCD::beta_function_nf_6{
        {
         21.0 / 3.0,
         78.0 / 3.0,
         -65.0 / 2.0,
         2472.2837425797160, }
    };

    /* 5 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 5
    const QCD::BetaFunction QCD::beta_function_nf_5{
        {
         23.0 / 3.0,
         116.0 / 3.0,
         9769.0 / 54.0,
         4826.1563287908967, }
    };

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 5
    const QCD::AnomalousMassDimension QCD::gamma_m_nf_5{
        {
         1.0, 506.0 / 9.0,
         474.87124557719461, 2824.7862379694232,
         }
    };

    /* 4 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 4
    const QCD::BetaFunction QCD::beta_function_nf_4{
        {
         25.0 / 3.0,
         154.0 / 3.0,
         21943.0 / 54.0,
         8035.1864197901160, }
    };

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 4
    const QCD::AnomalousMassDimension QCD::gamma_m_nf_4{
        {
         1.0, 526.0 / 9.0,
         636.61057670866927, 6989.5510103599477,
         }
    };

    /* 3 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 3
    const QCD::BetaFunction QCD::beta_function_nf_3{
        {
         9.0, 64.0,
         3863.0 / 6.0,
         12090.378130803711, }
    };

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 3
    const QCD::AnomalousMassDimension QCD::gamma_m_nf_3{
        { 1.0, 182.0 / 3.0, 794.89311771668714, 11331.304567227756 }
    };

    double
    QCD::alpha_s(const double & mu, const double & alpha_s_0, const double & mu_0, const BetaFunction & beta, unsigned int loop_order)
    {
        if (loop_order == 0)
        {
            throw InternalError("QCD::alpha_s: zero loop order not implemented");
        }
        if (loop_order > 4)
        {
            throw InternalError("QCD::alpha_s: loop order " + stringify(loop_order) + " not yet implemented");
        }

        double                a = alpha_s_0 / M_PI;
        std::array<double, 4> switches;
        std::generate(switches.begin(), switches.end(), [i = 0u, l = loop_order]() mutable { return (i++ < l) ? 1.0 : 0.0; });

        // Adjust for a different convention on beta function coefficients
        double beta0 = beta[0] / 4.0;
        double beta1 = switches[1] * beta[1] / 16.0;
        double beta2 = switches[2] * beta[2] / 64.0;
        double beta3 = switches[3] * beta[3] / 256.0;
        double b1    = beta1 / beta0;
        double b2    = beta2 / beta0;
        double b3    = beta3 / beta0;

        // cf. [CKS2000], Eq. (4), p. 3
        double ln_lambda2 = 2.0 * log(mu_0)
                            - (1.0 / a + b1 * log(a) + (b2 - b1 * b1) * a + (b3 / 2.0 - b1 * b2 + b1 * b1 * b1 / 2.0) * a * a) / beta0
                            // Use C for MSbar definition
                            - b1 / beta0 * log(beta0);

        double L = 2.0 * log(mu) - ln_lambda2, lnL = log(L);
        double denom = beta0 * L, denom2 = denom * denom, denom3 = denom2 * denom, denom4 = denom2 * denom2;

        // cf. [CKS2000], Eq. (5), p. 3
        double result = 1.0 / denom - b1 * lnL / denom2 + (b1 * b1 * (lnL * lnL - lnL - 1.0) + b2) / denom3
                        + (b1 * b1 * b1 * (-1.0 * lnL * lnL * lnL + 5.0 / 2.0 * lnL * lnL + 2.0 * lnL - 0.5) - 3.0 * b1 * b2 * lnL + b3 / 2.0) / denom4;

        return M_PI * result;
    }

    double
    QCD::m_q_msbar(const double & m_q_0, const double & alpha_s_0, const double & alpha_s_mu, const BetaFunction & beta, const AnomalousMassDimension & gamma_m,
                   unsigned int loop_order)
    {
        if (loop_order == 0)
        {
            throw InternalError("QCD::alpha_s: zero loop order not implemented");
        }
        if (loop_order > 4)
        {
            throw InternalError("QCD::alpha_s: loop order " + stringify(loop_order) + " not yet implemented");
        }

        double                a_mu0 = alpha_s_0 / M_PI;
        double                a_mu  = alpha_s_mu / M_PI;
        std::array<double, 4> switches;
        std::generate(switches.begin(), switches.end(), [i = 0u, l = loop_order]() mutable { return (i++ < l) ? 1.0 : 0.0; });

        // Adjust for a different convention on beta function coefficients
        double beta0 = beta[0] / 4.0;
        double beta1 = switches[1] * beta[1] / 16.0;
        double beta2 = switches[2] * beta[2] / 64.0;
        double beta3 = switches[3] * beta[3] / 256.0;
        double b1    = beta1 / beta0;
        double b2    = beta2 / beta0;
        double b3    = beta3 / beta0;

        // Adjust for a different convention on gamma function coefficients
        double gamma0_m = gamma_m[0];
        double gamma1_m = switches[1] * gamma_m[1] / 16.0;
        double gamma2_m = switches[2] * gamma_m[2] / 64.0;
        double gamma3_m = switches[3] * gamma_m[3] / 256.0;
        double c0       = gamma0_m / beta0;
        double c1       = gamma1_m / beta0;
        double c2       = gamma2_m / beta0;
        double c3       = gamma3_m / beta0;

        // cf. [CKS2000], Eq. (10), p. 6
        double c_mu0 = pow(a_mu0, c0)
                       * (1.0 + a_mu0 * (c1 - b1 * c0) + a_mu0 * a_mu0 * 0.5 * (pow(c1 - b1 * c0, 2) + c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                          + a_mu0 * a_mu0 * a_mu0
                                    * (pow(c1 - b1 * c0, 3) / 6.0 + (c1 - b1 * c0) / 2.0 * (c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                                       + (c3 - b1 * c2 + b1 * b1 * c1 - b2 * c1 - b1 * b1 * b1 * c0 + 2.0 * b1 * b2 * c0 - b3 * c0) / 3.0));

        // cf. [CKS2000], Eq. (10), p. 6
        double c_mu = pow(a_mu, c0)
                      * (1.0 + a_mu * (c1 - b1 * c0) + a_mu * a_mu * 0.5 * (pow(c1 - b1 * c0, 2) + c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                         + a_mu * a_mu * a_mu
                                   * (pow(c1 - b1 * c0, 3) / 6.0 + (c1 - b1 * c0) / 2.0 * (c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                                      + (c3 - b1 * c2 + b1 * b1 * c1 - b2 * c1 - b1 * b1 * b1 * c0 + 2.0 * b1 * b2 * c0 - b3 * c0) / 3.0));

        // cf. [CKS2000], Eq. (9), p. 6
        return m_q_0 * c_mu / c_mu0;
    }

    double
    QCD::m_q_msbar(const double & m_q_pole, const double & alpha_s, const double & nf)
    {
        double a_s = alpha_s / M_PI;

        // cf. [MvR1999], Eq. (12), pp. 4-5 for alpha_s = alpha_s(m_q_pole)
        // thus we return m_b(mu)
        return m_q_pole * (1.0 + a_s * (-4.0 / 3.0 + a_s * (1.04 * nf - 14.3323 + a_s * (-0.65269 * nf * nf + 26.9239 * nf - 198.8068))));
    }

    double
    QCD::m_q_pole(const double & m_q_MSbar, const double & alpha_s_mb, const double & nf, unsigned int loop_order)
    {
        double a_s = alpha_s_mb / M_PI;

        // cf. [CERN2003-002], Eq. (16), p. 45
        // Collect result from the inside out: m_q_MSbar * (1.0 + a_s ... * (... + a_s ... ) )
        double result = 0.0;
        switch (loop_order)
        {
            case 3:
                result = 190.8 - 26.7 * nf + 0.65 * power_of<2>(nf);
                // fall through
            case 2:
                result = 13.44 - 1.04 * nf + a_s * result;
                // fall through
            case 1:
                result = 4.0 / 3.0 + a_s * result;
                // fall through
            case 0:  result = 1.0 + a_s * result; break;
            default: throw InternalError("QCD::m_q_pole: loop order " + stringify(loop_order) + " not yet implemented");
        }
        return m_q_MSbar * result;
    }

    double
    QCD::m_q_ps(const double & m_q_MSbar, const double & alpha_s_mb, const double & mu_f, const double & nf, const QCD::BetaFunction & beta)
    {
        double a_s = alpha_s_mb / M_PI;
        double K   = 13.44 - 1.04 * nf;
        double a_1 = 10.33 - 1.11 * nf;
        double b_0 = beta[0];
        double L   = log(mu_f / m_q_MSbar);

        // cf. [B1998], Eq. (25), p. 12
        return m_q_MSbar * (1.0 + a_s * (4.0 / 3.0 * (1.0 - mu_f / m_q_MSbar) + a_s * (K - (mu_f / 3.0 / m_q_MSbar) * (a_1 - 2.0 * b_0 * (L - 1.0)))));
    }

    double
    QCD::m_q_kin(const double & m_q_MSbar, const double & alpha_s_mq, const double & mu, const BetaFunction & beta)
    {
        static const double zeta3 = 1.20206;
        static const double ln2   = log(2.0);
        static const double pi = M_PI, pi2 = pi * pi;

        double a_s = alpha_s_mq / M_PI, r = mu / m_q_MSbar;
        double b_0 = beta[0]; // We do not need to adjust for a factor of 4 when using [BBMU2003], Eq. (A.8).
        double L   = log(m_q_MSbar / (2.0 * mu));

        // cf. [BBMU2003], Eq. (A.8) and the underlying work [MvR2000]. Note that
        // the latter does use 4 beta_0 = beta_0|here.
        return m_q_MSbar
               * (1.0
                  + a_s
                            * (4.0 / 3.0 * (1.0 - 4.0 / 3.0 * r - r * r / 2.0)
                               + a_s
                                         * (b_0 / 2.0 * (pi2 / 6.0 + 71.0 / 48.0) + 665.0 / 144.0 + pi2 / 18.0 * (2.0 * ln2 - 19.0 / 2.0) - zeta3 / 6.0 - 8.0 / 3.0
                                            - r * (8.0 * b_0 / 9.0 * (L + 8.0 / 3.0) - 8.0 * pi2 / 9.0 + 52.0 / 9.0)
                                            - r * r * (b_0 / 3.0 * (L + 13.0 / 6.0) - pi2 / 3.0 + 23.0 / 18.0)
                                            + a_s * b_0 * b_0 / 4.0
                                                      * (2353.0 / 2592.0 + 13.0 / 36.0 * pi2 + 7.0 / 6.0 * zeta3
                                                         - r * 16.0 / 9.0 * (power_of<2>(L + 8.0 / 3.0) + 67.0 / 36.0 - pi2 / 6.0)
                                                         - r * r * 2.0 / 3.0 * (power_of<2>(L + 13.0 / 6.0) + 10.0 / 9.0 - pi2 / 6.0)))));
    }
} // namespace eos
