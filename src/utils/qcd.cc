/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/qcd.hh>

#include <cmath>

#include <iostream>

namespace wf
{
    double
    QCD::alpha_s(const double & mu)
    {
        /* QCD Running of alpha_s by 5 quark active quark flavors */
        static const double alpha_s_MZ = 0.1176; // cf.[PDG2008] p. 5, Table 1.1
        static const double MZ = 91.1876; // cf. [PDG2008] p. 9

        // cf. [vRVL1997] Eq. (11), n_f = 5
        static const double b_0 = 23.0/3.0;
        static const double b_1 = 116.0/3.0;
        static const double b_2 = 6519.0/54.0;

        // cf. [AEMSS2002], Eq. (5.7)
        double a = alpha_s_MZ / (4.0 * M_PI);
        double v0 = 1.0 + b_0 * 2.0 * a * std::log(mu / MZ);
        double v1 = -b_1 / b_0 * std::log(v0) / v0;
        double v2 = ((b_1 * b_1 / b_0 / b_0) * (1.0 - std::log(v0) + std::log(v0) * std::log(v0)) - b_2 / b_0 * (v0 - 1.0)) / (v0 * v0);
        double result = alpha_s_MZ / v0 * (1.0 + a * v1 + a * a * v2);

        return result;
    }

    double
    QCD::mb_pole(const double & mb_MSbar)
    {
        double a_s = alpha_s(mb_MSbar) / M_PI;

        // cf. [CERN2003-002], Eq. (16), p. 45
        return mb_MSbar * (1.0 + a_s * (4.0/3.0 + a_s * (8.24 + a_s * (73.55))));
    }

    double
    QCD::mb_PS(const double & mb_MSbar, const double & mu, const double & mu_PS)
    {
        double a_s = alpha_s(mu) / M_PI;

        // cf. [CERN2003-002], below Eq. (20), p. 46
        static const double a_1 = 43.0/3.0;
        static const double a_2 = 155.85;
        static const double b_0 = 23.0/3.0;
        static const double b_1 = 116.0/3.0;

        double casimir_f = 4.0 / 3.0;
        double L = std::log(mu_PS / mu);

        // cf. [CERN2003-002], Eq. (20), p. 46
        return mb_pole(mb_MSbar) -  casimir_f * a_s * mu_PS * (
                1.0 + a_s / 4.0 * (a_1 - b_0 * 2.0 * (L - 1.0))
                + 1.0 / 4.0 / 4.0 * a_s * a_s * (a_2 - (2.0 * a_1 * b_0 + b_1) * 2.0 * (L - 1.0)
                    + b_0 * b_0 * (L * L - 4.0 * L + 8.0)));
    }

    double
    QCD::mb_MSbar(const double & mb_MSbar0, const double & mu)
    {
        // 5 active flavors!

        static const double b_0 = 23.0/3.0;
        static const double b_1 = 116.0/3.0;
        static const double g_0 = 8.0;
        static const double g_1 = 1012/9.0;
        static const double rho = g_0 / (2.0 * b_0);

        double alpha_m = alpha_s(mb_MSbar0), alpha_mu = alpha_s(mu);
        double eta = alpha_mu / alpha_m;

        // cf. [BBL1995], Eq. (III.20), p. 17
        return mb_MSbar0 * pow(eta, rho) * (1.0 + alpha_m / (4.0 * M_PI) * rho * (g_1 / g_0 - b_1 / b_0) * (eta - 1));
    }

    double
    QCD::mc_MSbar(const double & mc_MSbar0, const double & mu)
    {
        // 4 active flavors!

        static const double b_0 = 25.0/3.0;
        static const double b_1 = 154.0/3.0;
        static const double g_0 = 8.0;
        static const double g_1 = 1052.0/9.0;
        static const double rho = g_0 / (2.0 * b_0);

        double alpha_m = alpha_s(mc_MSbar0), alpha_mu = alpha_s(mu);
        double eta = alpha_mu / alpha_m;

        // cf. [BBL1995], Eq. (III.20), p. 17
        return mc_MSbar0 * pow(eta, rho) * (1.0 + alpha_m / (4.0 * M_PI) * rho * (g_1 / g_0 - b_1 / b_0) * (eta - 1));
    }
}
