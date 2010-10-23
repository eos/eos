/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/qcd.hh>

#include <cmath>

namespace eos
{
    const double QCD::casimir_f = 4.0 / 3.0;

    double
    QCD::alpha_s(const double & mu, const double & alpha_s_0, const double & mu_0, const Parameters & params)
    {
        double a = alpha_s_0 / M_PI;
        // Adjust for a different convention on beta function coefficients
        double beta0 = params._beta0 / 4.0;
        double beta1 = params._beta1 / 16.0;
        double beta2 = params._beta2 / 64.0;
        double beta3 = params._beta3 / 256.0;
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
    QCD::m_q_msbar(const double & m_q_0, const double & alpha_s_0, const double & alpha_s_mu, const QCD::Parameters & params)
    {
        double a_mu0 = alpha_s_0 / M_PI;
        double a_mu = alpha_s_mu / M_PI;

        // Adjust for a different convention on beta function coefficients
        double beta0 = params._beta0 / 4.0;
        double beta1 = params._beta1 / 16.0;
        double beta2 = params._beta2 / 64.0;
        double beta3 = params._beta3 / 256.0;
        double b1 = beta1 / beta0;
        double b2 = beta2 / beta0;
        double b3 = beta3 / beta0;

        double gamma0_m = params._gamma0_m;
        // Adjust for a different convention on gamma function coefficients
        double gamma1_m = params._gamma1_m / 16.0;
        double gamma2_m = params._gamma2_m / 64.0;
        double gamma3_m = params._gamma3_m / 256.0;
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
    QCD::m_q_pole(const double & m_q_MSbar, const double & alpha_s_mb, const QCD::Parameters & params)
    {
        double a_s = alpha_s_mb / M_PI;
        double nf = params._nf;

        // cf. [CERN2003-002], Eq. (16), p. 45
        return m_q_MSbar * (1.0 + a_s * (4.0/3.0 + a_s * (13.44 - 1.04 * nf + a_s * (190.8 - 26.7 * nf + 0.65 * nf * nf))));
    }

    double
    QCD::m_q_ps(const double & m_q_MSbar, const double & alpha_s_mb, const double & mu_f, const QCD::Parameters & params)
    {
        double a_s = alpha_s_mb / M_PI;
        double K = 13.44 - 1.04 * params._nf;
        double a_1 = 10.33 - 1.11 * params._nf;
        double b_0 = params._beta0;
        double L = log(mu_f / m_q_MSbar);

        // cf. [CERN2003-002], Eq. (16), p. 45
        return m_q_MSbar * (1.0 + a_s * (4.0/3.0 * (1.0 - mu_f / m_q_MSbar) + a_s * (K - (mu_f / 3.0 / m_q_MSbar) * (a_1 - 2.0 * b_0 * (L - 1.0)))));
    }
}
