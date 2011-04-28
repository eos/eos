/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <src/utils/qcd.hh>

#include <cmath>

namespace eos
{
    const double QCD::casimir_f = 4.0 / 3.0;

    /* 6 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 6
    const QCD::BetaFunction QCD::beta_function_nf_6
    {{
        21.0 / 3.0,
        78.0 / 3.0,
        -65.0 / 2.0,
        2472.2837425797160,
    }};

    /* 5 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 5
    const QCD::BetaFunction QCD::beta_function_nf_5
    {{
        23.0 / 3.0,
        116.0 / 3.0,
        9769.0 / 54.0,
        4826.1563287908967,
    }};

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 5
    const QCD::AnomalousMassDimension QCD::gamma_m_nf_5
    {{
        1.0,
        506.0 / 9.0,
        474.87124557719461,
        2824.7862379694232,
    }};

    /* 4 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 4
    const QCD::BetaFunction QCD::beta_function_nf_4
    {{
        25.0 / 3.0,
        154.0 / 3.0,
        21943.0 / 54.0,
        8035.1864197901160,
    }};

    // cf. [CKS2000], Eq. (7), p. 5 with n_f = 4
    const QCD::AnomalousMassDimension QCD::gamma_m_nf_4
    {{
        1.0,
        526.0 / 9.0,
        636.61057670866927,
        6989.5510103599477,
    }};

    /* 3 flavor QCD constants */

    // cf. [CKS2000], Eq. (2), p. 2 with n_f = 3
    const QCD::BetaFunction QCD::beta_function_nf_3
    {{
        9.0,
        64.0,
        3863.0 / 6.0,
        12090.378130803711,
    }};

    double
    QCD::alpha_s(const double & mu, const double & alpha_s_0, const double & mu_0, const BetaFunction & beta)
    {
        double a = alpha_s_0 / M_PI;
        // Adjust for a different convention on beta function coefficients
        double beta0 = beta[0] / 4.0;
        double beta1 = beta[1] / 16.0;
        double beta2 = beta[2] / 64.0;
        double beta3 = beta[3] / 256.0;
        double b1 = beta1 / beta0;
        double b2 = beta2 / beta0;
        double b3 = beta3 / beta0;

        // cf. [CKS2000], Eq. (4), p. 3
        double ln_lambda2 = 2.0 * log(mu_0)
            - (1.0 / a + b1 * log(a) + (b2 - b1 * b1) * a + (b3 / 2.0 - b1 * b2 + b1 * b1 * b1 / 2.0) * a * a) / beta0
            // Use C for MSbar definition
            - b1 / beta0 * log(beta0);

        double L = 2.0 * log(mu) - ln_lambda2, lnL = log(L);
        double denom = beta0 * L, denom2 = denom * denom, denom3 = denom2 * denom, denom4 = denom2 * denom2;

        // cf. [CKS2000], Eq. (5), p. 3
        double result = 1.0 / denom
            - b1 * lnL / denom2
            + (b1 * b1 * (lnL * lnL - lnL - 1.0) + b2) / denom3
            + (b1 * b1 * b1 * (-1.0 * lnL * lnL * lnL + 5.0 / 2.0 * lnL * lnL + 2.0 * lnL - 0.5) - 3.0 * b1 * b2 * lnL + b3 / 2.0) / denom4;

        return M_PI * result;
    }

    double
    QCD::m_q_msbar(const double & m_q_0, const double & alpha_s_0, const double & alpha_s_mu, const BetaFunction & beta, const AnomalousMassDimension & gamma_m)
    {
        double a_mu0 = alpha_s_0 / M_PI;
        double a_mu = alpha_s_mu / M_PI;

        // Adjust for a different convention on beta function coefficients
        double beta0 = beta[0] / 4.0;
        double beta1 = beta[1] / 16.0;
        double beta2 = beta[2] / 64.0;
        double beta3 = beta[3] / 256.0;
        double b1 = beta1 / beta0;
        double b2 = beta2 / beta0;
        double b3 = beta3 / beta0;

        // Adjust for a different convention on gamma function coefficients
        double gamma0_m = gamma_m[0];
        double gamma1_m = gamma_m[1] / 16.0;
        double gamma2_m = gamma_m[2] / 64.0;
        double gamma3_m = gamma_m[3] / 256.0;
        double c0 = gamma0_m / beta0;
        double c1 = gamma1_m / beta0;
        double c2 = gamma2_m / beta0;
        double c3 = gamma3_m / beta0;

        // cf. [CKS2000], Eq. (10), p. 6
        double c_mu0 = pow(a_mu0, c0) * (
                1.0
                + a_mu0 * (c1 - b1 * c0)
                + a_mu0 * a_mu0 * 0.5 * (pow(c1 - b1 * c0, 2) + c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                + a_mu0 * a_mu0 * a_mu0 * (pow(c1 - b1 * c0, 3) / 6.0 + (c1 - b1 * c0) / 2.0 * (c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                    + (c3 - b1 * c2 + b1 * b1 * c1 - b2 * c1 - b1 * b1 * b1 * c0 + 2.0 * b1 * b2 * c0 - b3 * c0)));

        // cf. [CKS2000], Eq. (10), p. 6
        double c_mu = pow(a_mu, c0) * (
                1.0
                + a_mu * (c1 - b1 * c0)
                + a_mu * a_mu * 0.5 * (pow(c1 - b1 * c0, 2) + c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                + a_mu * a_mu * a_mu * (pow(c1 - b1 * c0, 3) / 6.0 + (c1 - b1 * c0) / 2.0 * (c2 - b1 * c1 + b1 * b1 * c0 - b2 * c0)
                    + (c3 - b1 * c2 + b1 * b1 * c1 - b2 * c1 - b1 * b1 * b1 * c0 + 2.0 * b1 * b2 * c0 - b3 * c0)));

        // cf. [CKS2000], Eq. (9), p. 6
        return m_q_0 * c_mu / c_mu0;
    }

    double
    QCD::m_q_pole(const double & m_q_MSbar, const double & alpha_s_mb, const double & nf)
    {
        double a_s = alpha_s_mb / M_PI;

        // cf. [CERN2003-002], Eq. (16), p. 45
        return m_q_MSbar * (1.0 + a_s * (4.0/3.0 + a_s * (13.44 - 1.04 * nf + a_s * (190.8 - 26.7 * nf + 0.65 * nf * nf))));
    }

    double
    QCD::m_q_ps(const double & m_q_MSbar, const double & alpha_s_mb, const double & mu_f, const double & nf, const QCD::BetaFunction & beta)
    {
        double a_s = alpha_s_mb / M_PI;
        double K = 13.44 - 1.04 * nf;
        double a_1 = 10.33 - 1.11 * nf;
        double b_0 = beta[0];
        double L = log(mu_f / m_q_MSbar);

        // cf. [B1998], Eq. (25), p. 12
        return m_q_MSbar * (1.0 + a_s * (4.0/3.0 * (1.0 - mu_f / m_q_MSbar) + a_s * (K - (mu_f / 3.0 / m_q_MSbar) * (a_1 - 2.0 * b_0 * (L - 1.0)))));
    }
}
