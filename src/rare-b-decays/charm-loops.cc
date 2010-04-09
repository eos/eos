/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/charm-loops.hh>

#include <cmath>
#include <complex>

#include <gsl/gsl_sf_dilog.h>

namespace wf
{
    using std::atan;
    using std::complex;
    using std::log;
    using std::pow;
    using std::sqrt;

    complex<double>
    CharmLoops::h(const double & mu, const double & s)
    {
        // cf. [BFS2001], Eq. (11), p. 4 in the limit m_q -> 0
        return 4.0 / 9.0 * complex<double>(2.0 / 3.0 + 2.0 * std::log(mu) - std::log(s), M_PI);
    }

    complex<double>
    CharmLoops::h(const double & mu, const double & s, const double & m_q)
    {
        if (m_q < 1e-4)
            return h(mu, s);

        const double z = 4.0 * m_q * m_q / s;
        if (z < 1e-10)
            return complex<double>(-4.0/9.0 * (2.0 * std::log(m_q / mu) + 1.0), 0.0);

        const double sqrt1z = std::sqrt(std::abs(z - 1.0));

        double a = 2.0 * std::log(m_q / mu) - 2.0 / 3.0 - z;
        double b = (2.0 + z) * sqrt1z;
        double rc, ic;
        if (z > 1.0)
        {
            ic = 0.0;
            rc = std::atan(1.0 / sqrt1z);
        }
        else
        {
            ic = -M_PI / 2.0;
            rc = std::log((1.0 + sqrt1z) / std::sqrt(z));
        }

        // cf. [BFS2001], Eq. (11), p. 4
        return -4.0 / 9.0 * (a + b * complex<double>(rc, ic));
    }

    complex<double>
    CharmLoops::A(const double & mu, const double & s, const double & m_b)
    {
            double s_hat = s / m_b / m_b, s_hat2 = s_hat * s_hat, denom = 1 - s_hat;
            double ln = log(s_hat), ln2 = ln * ln, li_2 = gsl_sf_dilog(s_hat);
            double z = 4.0 / s_hat, sqrt1z = sqrt(z - 1.0);

            double a = -104.0 / 243.0 * 2.0 * log(m_b / mu);
            double b = +4.0 * s_hat / 27.0 / denom * (li_2 + ln * log(1.0 - s_hat));
            complex<double> c = +1.0 / 729.0 / pow(denom, 2) * complex<double>(
                    6.0 * s_hat * (29.0 - 47.0 * s_hat) * ln + 785.0 - 1600.0 * s_hat + 833.0 * s_hat2,
                    6.0 * M_PI * (20.0 - 49.0 * s_hat + 47.0 * s_hat2));
            complex<double> d = -2.0 / 243.0 / pow(denom, 3) * complex<double>(
                    2.0 * sqrt1z * (-4.0 + 9.0 * s_hat - 15.0 * s_hat2 + 4.0 * s_hat * s_hat2) * (M_PI / 2.0 - atan(sqrt1z))
                    + 9.0 * s_hat * s_hat2 * ln2,
                    18.0 * M_PI * s_hat * (1.0 - 2.0 * s_hat) * ln);
            double e = +2.0 * s_hat / 243.0 / pow(denom, 4) * (36.0 * pow(M_PI / 2.0 - atan(sqrt1z), 2) +
                    M_PI * M_PI * (-4.0 + 9.0 * s_hat - 9.0 * s_hat2 + 3.0 * s_hat * s_hat2));

            // cf. [S2004], Eq. (29), p. 8
            complex<double> result = a + b + c + d + e;

            return result;
    }

    complex<double>
    CharmLoops::B(const double & mu, const double & s, const double & m_b)
    {
            double s_hat = s / m_b / m_b, s_hat2 = s_hat * s_hat, denom = 1 - s_hat;
            double ln = log(s_hat), ln2 = ln * ln, li_2 = gsl_sf_dilog(s_hat);
            double z = 4.0 / s_hat, sqrt1z = sqrt(z - 1.0);
            double lnmu = 2.0 * log(m_b / mu);

            complex<double> x1(0.5, 0.5 * sqrt1z), x2(0.5, -0.5 * sqrt1z), x3(0.5, 0.5 / sqrt1z), x4(0.5, -0.5 / sqrt1z);
            gsl_sf_result rp, ip;
            gsl_sf_complex_dilog_e(1.0, arg(-x2 / x1), &rp, &ip);
            complex<double> dilog(rp.val, ip.val);
            complex<double> lx1 = log(x1), lx2 = log(x2), lx3 = log(x3), lx4 = log(x4);


            complex<double> a = 8.0 / 243.0 / s_hat * (complex<double>(4.0 - 34.0 * s_hat, -17.0 * M_PI * s_hat) * lnmu
                    + 8.0 * s_hat * lnmu * lnmu + 17.0 * s_hat * ln * lnmu);
            complex<double> b = (2.0 + s_hat) * sqrt1z / 729.0 / s_hat * (
                        -48.0 * lnmu * (M_PI / 2.0 - atan(sqrt1z)) - 18.0 * M_PI * 2.0 * log(sqrt1z) -12.0 * M_PI * (2.0 * lx1 + lx3 + lx4)
                        + complex<double>(0.0, 1.0) * (
                            3.0 * 4.0 * pow(log(sqrt1z), 2.0) - 5.0 * M_PI * M_PI - 24.0 * dilog
                            +6.0 * (-9.0 * pow(lx1, 2.0) + pow(lx2, 2.0) - 2.0 * pow(lx4, 2.0) + 6.0 * lx1 * lx2
                            -4.0 * lx1 * lx3 + 8.0 * lx1 * lx4))
                        );
            complex<double> c = -2.0 / 243.0 / s_hat / denom * (
                    4.0 * s_hat * (-8.0 + 17.0 * s_hat) * (li_2 + ln * log(1.0 - s_hat))
                    + 3.0 * (2.0 + s_hat) * (3.0 - s_hat) * pow(lx2 - lx1, 2.0)
                    + 12.0 * M_PI * (-6.0 - s_hat + s_hat2) * (M_PI / 2.0 - atan(sqrt1z)));
            complex<double> d = 2.0 / (2187.0 * s_hat * pow(denom, 2.0)) * complex<double>(
                    -18.0 * s_hat * (120.0 - 211.0 * s_hat + 73.0 * s_hat2) * ln
                    -288.0 - 8.0 * s_hat + 934.0 * s_hat2 - 692.0 * s_hat * s_hat2,
                    18.0 * M_PI * s_hat * (82.0 - 173.0 * s_hat + 73.0 * s_hat2));
            complex<double> e = -4.0 / (243.0 * s_hat * pow(denom, 3.0)) * std::complex<double>(
                    -2.0 * sqrt1z * (4.0 - 3.0 * s_hat - 18.0 * s_hat2 + 16.0 * s_hat * s_hat2 - 5.0 * s_hat2 * s_hat2) * (M_PI / 2 - atan(sqrt1z))
                    -9.0 * s_hat * s_hat2 * ln2,
                    2.0 * M_PI * s_hat * (8.0 - 33.0 * s_hat + 51.0 * s_hat2 - 17.0 * s_hat * s_hat2) * ln);
            complex<double> f = 2.0 / (729.0 * s_hat * pow(denom, 4.0)) * (72.0 * (3.0 - 8.0 * s_hat + 2.0 * s_hat2)
                    * pow(M_PI / 2.0 - atan(sqrt1z), 2.0)
                    - M_PI * M_PI * (54.0 - 53.0 * s_hat - 286.0 * s_hat2 + 612.0 * s_hat * s_hat2 - 446 * s_hat2 * s_hat2 + 113.0 * s_hat2 * s_hat2 * s_hat));

            // cf. [S2004], Eq. (29), p. 8
            return a + b + c + d + e + f;
    }

    complex<double>
    CharmLoops::C(const double & mu, const double & s, const double & m_b)
    {
        static const double zeta3 = 1.20206;

        return complex<double>(16.0 / 81.0 * log(mu * mu / s) + 428.0 / 243.0 - 64.0 / 27.0 * zeta3,
                16.0 / 81.0 * M_PI);
    }

    complex<double>
    CharmLoops::F17(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (22), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return -A(mu, s, m_b);
    }

    complex<double>
    CharmLoops::F19(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (24), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return -B(mu, s, m_b) - 4.0 * C(mu, s, m_b);
    }

    complex<double>
    CharmLoops::F27(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (23), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return 6.0 * A(mu, s, m_b);
    }

    complex<double>
    CharmLoops::F29(const double & mu, const double & s, const double & m_b)
    {
        // cf. [S2004], Eq. (25), p. 7 and consider a global sign (compare [ABGW2003], Eq. (7), p. 8 with [S2004], Eq. (16), p. 6)
        return 6.0 * B(mu, s, m_b) - 3.0 * C(mu, s, m_b);
    }
}
