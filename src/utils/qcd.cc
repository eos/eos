/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/qcd.hh>

#include <cmath>

#include <iostream>

namespace wf
{
    double
    QCD::alpha_s(const double & mu)
    {
        static const double alpha_s_MZ = 0.1176; // cf.[PDG2008] p. 5, Table 1.1
        static const double b_0 = 11.0 - 10.0 / 3; // cf. [PS] Eq. (17.13), n_f = 5
        static const double MZ = 91.1876; // cf. [PDG2008] p. 9


        double result = alpha_s_MZ / (1 + b_0 * alpha_s_MZ / (2.0 * M_PI) * std::log(mu / MZ)); // cf. [PS] Eq. (17.14)

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
