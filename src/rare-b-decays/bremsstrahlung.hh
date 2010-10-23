/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_BREMSSTRAHLUNG_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_BREMSSTRAHLUNG_HH 1

#include <src/utils/complex.hh>

namespace eos
{
    struct Bremsstrahlung
    {
        // cf. [AAGW2002], Eqs. (30), (31), p. 12
        static complex<double> G_m1(const double & t);
        static complex<double> G_0(const double & t);

        // cf. [AAGW2002], Eqs. (28), (29), p. 11
        static complex<double> Deltai_23(const double & s_hat, const double & w, const double & z);
        static complex<double> Deltai_27(const double & s_hat, const double & w, const double & z);

        // cf. [AAGW2002], Eqs. (23)-(26), p. 10
        static complex<double> tau_22(const double & s_hat, const double & w, const double & z);
        static complex<double> tau_27(const double & s_hat, const double & w, const double & z);
        static complex<double> tau_28(const double & s_hat, const double & w, const double & z);
        static complex<double> tau_29(const double & s_hat, const double & w, const double & z);

        // Integrals of tau_2x from w = s_hat to w = 1, cf [AAGW2002], Eq. (22)
        static complex<double> itau_22(const double & s_hat, const double & z);
        static complex<double> itau_27(const double & s_hat, const double & z);
        static complex<double> itau_28(const double & s_hat, const double & z);
        static complex<double> itau_29(const double & s_hat, const double & z);
    };
}

#endif
