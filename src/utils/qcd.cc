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
        static const double b_3 = 4826.1;

        // cf. [AEMSS2002], Eq. (5.7)
        double a = alpha_s_MZ / (4.0 * M_PI);
        double v0 = 1.0 + b_0 * 2.0 * a * std::log(mu / MZ);
        double v1 = -b_1 / b_0 * std::log(v0) / v0;
        double v2 = ((b_1 * b_1 / b_0 / b_0) * (1.0 - std::log(v0) + std::log(v0) * std::log(v0)) - b_2 / b_0 * (v0 - 1.0)) / (v0 * v0);
        double result = alpha_s_MZ / v0 * (1.0 + a * v1 + a * a * v2);

        return result;
    }

    double
    QCD::mb_PS(const double & mb_MSbar, const double & mu)
    {
        double alpha = alpha_s(mu);

        return mb_pole(mb_MSbar, mu) - alpha * mu * 4.0 / 3.0 / M_PI;
    }

    double
    QCD::mb_pole(const double & mb_MSbar, const double & mu)
    {
        double alpha = alpha_s(mb_MSbar);

        return mb_MSbar * (1.0 + 4.0 / 3.0 * alpha);
    }
}
