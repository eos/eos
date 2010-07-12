/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/bremsstrahlung.hh>
#include <src/utils/integrate.hh>

#include <cmath>
#include <limits>
#include <tr1/functional>

namespace wf
{
    // cf. [AAGW2002], Eq. (30), p. 12
    complex<double>
    Bremsstrahlung::G_m1(const double & t)
    {
        complex<double> result;

        if (t < 4)
        {
            double x = atan(sqrt((4.0 - t) / t));

            result = 2.0 * M_PI * x - pow(M_PI, 2) / 2 - 2.0 * pow(x, 2);
        }
        else
        {
            double x = log((sqrt(t) + sqrt(t - 4.0)) / 2.0);

            result = complex<double>(-pow(M_PI, 2) / 2.0 + 2.0 * pow(x, 2),
                    -2.0 * M_PI * x);
        }

        return result;
    }

    // cf. [AAGW2002], Eq. (31), p. 12
    complex<double>
    Bremsstrahlung::G_0(const double & t)
    {
        complex<double> result;

        if (t < 4)
        {
            double x = sqrt((4.0 - t) / t);

            result = M_PI * x - 2 - 2 * x * atan(x);
        }
        else
        {
            double x = sqrt((t - 4.0) / t);
            double y = log((sqrt(t) + sqrt(t - 4.0)) / 2.0);

            result = complex<double>(-2.0 + 2.0 * x * y,
                    -M_PI * x);
        }

        return result;
    }

    // cf. [AAGW2002], Eq. (28), p. 11
    complex<double>
    Bremsstrahlung::Deltai_23(const double & s_hat, const double & w, const double & z)
    {
        return -2.0 + 4.0 / (w - s_hat) * (
                z * (Bremsstrahlung::G_m1(s_hat / z) - Bremsstrahlung::G_m1(w / z))
                - s_hat / 2.0 * (Bremsstrahlung::G_0(s_hat / z) - Bremsstrahlung::G_0(w / z)));
    }

    // cf. [AAGW2002], Eq. (29), p. 11
    complex<double>
    Bremsstrahlung::Deltai_27(const double & s_hat, const double & w, const double & z)
    {
        return 2.0 * (Bremsstrahlung::G_0(s_hat / z) - Bremsstrahlung::G_0(w / z));
    }

    // cf. [AAGW2002], Eq. (23), p. 10
    complex<double>
    Bremsstrahlung::tau_22(const double & s_hat, const double & w, const double & z)
    {
        double s_hat2 = s_hat * s_hat, w2 = w * w, w3 = w2 * w;

        return 8.0 / 27.0 * (w - s_hat) * pow(1 - w, 2) / s_hat / w3 * (
                (3.0 * w2 + 2 * s_hat2 * (2.0 + w) - s_hat * w * (5 - 2.0 * w)) * norm(Deltai_23(s_hat, w, z))
                + (2.0 * s_hat2 * (2.0 + w) + s_hat * w * (1.0 + 2.0 * w)) * norm(Deltai_27(s_hat, w, z))
                + 4.0 * s_hat * (w * (1.0 - w) - s_hat * (2.0 + w)) * real(Deltai_23(s_hat, w, z) * conj(Deltai_27(s_hat, w, z))));
    }

    // cf. [AAGW2002], Eq. (24), p. 10
    complex<double>
    Bremsstrahlung::tau_27(const double & s_hat, const double & w, const double & z)
    {
        double s_hat2 = s_hat * s_hat, w2 = w * w;

        return 8.0 / 3.0 / (s_hat * w) * (
                ((1.0 - w) * (4.0 * s_hat2 - s_hat * w + w2) + s_hat * w * (4.0 + s_hat - w) * log(w)) * Deltai_23(s_hat, w, z)
                - (4.0 * s_hat2 * (1.0 - w) + s_hat * w * (4.0 + s_hat - w) * log(w)) * Deltai_27(s_hat, w, z));
    }

    // cf. [AAGW2002], Eq. (25), p. 10
    complex<double>
    Bremsstrahlung::tau_28(const double & s_hat, const double & w, const double & z)
    {
        double w2 = w * w;
        double x = s_hat / (1.0 + s_hat - w) / (w2 + s_hat * (1.0 - w));

        return 8.0 / 9.0 / (s_hat * w * (w - s_hat)) * (
                (pow(w - s_hat, 2) * (2.0 * s_hat - w) * (1.0 - w)) * Deltai_23(s_hat, w, z)
                - (2.0 * s_hat * pow(w - s_hat, 2) * (1.0 - w)) * Deltai_27(s_hat, w, z)
                + s_hat * w * ((1.0 + 2.0 * s_hat - 2.0 * w) * Deltai_23(s_hat, w, z)
                    - 2.0 * (1.0 + s_hat - w) * Deltai_27(s_hat, w, z)) * log(x));
    }

    // cf. [AAGW2002], Eq. (24), p. 10
    complex<double>
    Bremsstrahlung::tau_29(const double & s_hat, const double & w, const double & z)
    {
        return 4.0 / 3.0 / w * (
                (2.0 * s_hat * (1.0 - w) * (s_hat + w) + 4.0 * s_hat * w * log(w)) * Deltai_23(s_hat, w, z)
                - (2.0 * s_hat * (1.0 - w) * (s_hat + w) + w * (3.0 * s_hat + w) * log(w)) * Deltai_27(s_hat, w, z));
    }

    // Integrals of tau_2x from w = s_hat to w = 1
    complex<double>
    Bremsstrahlung::itau_22(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_22, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }

    complex<double>
    Bremsstrahlung::itau_27(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_27, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }

    complex<double>
    Bremsstrahlung::itau_28(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_28, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }

    complex<double>
    Bremsstrahlung::itau_29(const double & s_hat, const double & z)
    {
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        if (1.0 - s_hat < eps)
            return 0.0;

        return integrate(std::function<complex<double> (const double &)>(
                    std::bind(&Bremsstrahlung::tau_29, s_hat, std::placeholders::_1, z)),
                128, s_hat + eps, 1.0);
    }
}
