/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/bremsstrahlung.hh>
#include <src/rare-b-decays/charm-loops.hh>
#include <src/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <tr1/functional>
#include <utility>
#include <map>
#include <vector>

#include <gsl/gsl_sf_dilog.h>

namespace wf
{
    /* HLMW2005 */

    template <>
    struct Implementation<BToXsDilepton<HLMW2005>>
    {
        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c7;

        Parameter c8;

        Parameter c9;

        Parameter c10;

        Parameter m_b_MSbar;

        Parameter m_c_MSbar;

        Parameter m_s;

        Parameter m_tau;

        Parameter m_Z;

        Parameter br_clnu;

        Parameter lambda_2;

        Parameter mu;

        Parameter ckm;

        Parameter C;

        double m_l;

        Implementation(const Parameters & p) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["Re{c7}"]),
            c8(p["c8"]),
            c9(p["Re{c9}"]),
            c10(p["Re{c10}"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c_MSbar(p["mass::c"]),
            m_s(p["mass::s"]),
            m_tau(p["mass::tau"]),
            m_Z(p["mass::Z"]),
            br_clnu(p["exp::BR(B->X_clnu)"]),
            lambda_2(p["B->X_s::lambda_2"]),
            mu(p["mu"]),
            ckm(p["exp::CKM(B->X_sll, B->X_clnu)"]),
            C(p["exp::C(B->X_clnu, B->X_ulnu)"])
        {
            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.10565836; // (GeV), cf. [PDG2008], p. 13
        }

        double m_b_pole() const
        {
            return 4.8;
        }

        double m_c_pole() const
        {
            return 1.4;
        }

        double s_hat(const double & s) const
        {
            double m_b = m_b_pole();

            return s / m_b / m_b;
        }

        /* NLO functions */
        // cf. [BFS2001], Eq. (29), p. 8
        complex<double> B0(const double & s, const double & m_q) const
        {
            double z = 4.0 * m_q * m_q / s;
            double rp, ip;

            if (m_q < 0.0)
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q < 0!");
            };

            if ((0.0 == m_q) && (0.0 == s))
            {
                throw InternalError("Implementation<BToKstarDilepton<LargeRecoil>>::B0: m_q == 0 & s == 0");
            }

            if (0 == s)
                return complex<double>(-2.0, 0.0);

            if (z > 1.0)
            {
                rp = -2.0 * std::sqrt(z - 1.0) * std::atan(1.0 / std::sqrt(z - 1.0));
                ip = 0.0;
            }
            else
            {
                rp = std::sqrt(1.0 - z) * std::log((1.0 - std::sqrt(1 - z)) / (1.0 + std::sqrt(1.0 - z)));
                ip = std::sqrt(1.0 - z) * M_PI;
            }

            return complex<double>(rp, ip);
        }

        // cf. [BFS2001], Eq. (84), p. 30
        double C0(const double & s_hat) const
        {
            static const int points = 40;
            // Integration boundaries of u = (0, u_max]
            static const double dx = (1.0 - 0.0) / points;
            static const double g[2] =
            {
                (1.0 + std::sqrt(3.0 / 5.0)) / 2.0,
                (1.0 - std::sqrt(3.0 / 5.0)) / 2.0
            };

            long double x1, x2, x;
            long double y[3];
            long double result = 0;

            for (int i = 0; i < points; i++)
            {
                // 0 is lower integration boundary

                x1 = 0.0 + i * dx;
                x2 = 0.0 + (i + 1) * dx;
                for (int j = 0; j < 3; j++)
                {
                    switch (j)
                    {
                        case 0: x = x1 * g[0] + x2 * g[1];
                                break;

                        case 1: x = (x1 + x2) / 2.0;
                                break;

                        case 2: x = x1 * g[1] + x2 * g[0];
                                break;
                    }

                    y[j] = std::log(x * x / (1.0 - s_hat * x * (1.0 - x))) / (1.0 + x * (1.0 - s_hat));
                }

                result += (5.0 * y[0] + 8.0 * y[1] + 5.0 * y[2]) * dx / 18.0;
            };

            return result;
        }

        // cf. [BFS2001], Eq. (82), p. 30
        complex<double> F87(const double & s_hat) const
        {
            double m_b = m_b_pole();
            double s_hat2 = s_hat * s_hat;
            double s = s_hat * m_b * m_b;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            complex<double> a = complex<double>(-32.0 * std::log(mu / m_b) - 8.0 * s_hat / denom * std::log(s_hat)
                    - 4.0 * (11.0 - 16.0 * s_hat + 8.0 * s_hat2) / denom2,
                    -8.0 * M_PI);
            complex<double> b = (4.0 / denom / denom2)
                * ((9.0 * s_hat - 5.0 * s_hat2 + 2.0 * s_hat * s_hat2) * B0(s, m_b) - (4.0 + 2.0 * s_hat) * C0(s_hat));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BFS2001], Eq. (83), p. 30
        complex<double> F89(const double & s_hat) const
        {
            double m_b = m_b_pole();
            double s = s_hat * m_b * m_b;
            double denom = (1.0 - s_hat);
            double denom2 = denom * denom;

            double a = 16.0 * std::log(s_hat) / denom + 8.0 * (5.0 - 2.0 * s_hat) / denom2;
            complex<double> b = (-8.0 * (4.0 - s_hat) / denom / denom2) * ((1.0 + s_hat) * B0(s, m_b) - 2.0 * C0(s_hat));

            return (1.0 / 9.0) * (a + b);
        }

        // cf. [BMU1999], Eq. (34), p. 9 and [HLMW2005], Eq. (127), p. 29
        double omega1_99(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat;

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (5.0 + 4.0 * s_hat) / (3.0 * (1.0 + 2.0 * s_hat)) * ln1
                - (2.0 * s_hat * (1.0 + s_hat) * (1.0 - 2.0 * s_hat)) / (3.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * ln
                + (5.0 + 9.0 * s_hat - 6.0 * s_hat2) / (6.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat));
        }

        // cf. [HLMW2005], Eq. (128), p. 29
        // only valid for 0 < s_hat < 0.4
        double omega2_99(const double & s_hat) const
        {
            double ln = log(s_hat);;
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;

            return -19.2 + 6.1 * s_hat + (37.9 + 17.2 * ln) * s_hat2 - 18.7 * s_hat3;
        }

        // cf. [HLMW2005], Eq. (130), p. 29
        double omega1_77(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat;

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (8.0 + s_hat) / (3.0 * (2.0 + s_hat)) * ln1
                - (2.0 * s_hat * (2.0 - 2.0 * s_hat - s_hat2)) / (3.0 * pow(1.0 - s_hat, 2) * (2.0 + s_hat)) * ln
                - (16.0 - 11.0 * s_hat - 17.0 * s_hat2) / (8.0 * (1.0 - s_hat) * (2.0 + s_hat))
                // We use mu_b in MSbar scheme globally, so use m_b_MSbar here instead of m_b_pole
                - 8.0/3.0 * log(mu / m_b_MSbar());
        }

        // cf.[HLMW2005], Eq. (131), p. 29
        double omega1_79(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (2.0 + 7.0 * s_hat) / (9.0 * s_hat) * ln1
                - (2.0 * s_hat * (3.0 - 2.0 * s_hat)) / (9.0 * pow(1.0 - s_hat, 2)) * ln
                + (5.0 - 9.0 * s_hat) / (18.0 * (1.0 - s_hat))
                // We use mu_b in MSbar scheme globally, so use m_b_MSbar here instead of m_b_pole
                - 4.0/3.0 * log(mu / m_b_MSbar());
        }

        // cf. [HLMW2005], Eq. (94), p. 23
        double omegaem_99(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;

            return 2.0 * log(m_b_pole() / m_l) * (
                    - (1.0 + 4.0 * s_hat - 8.0 * s_hat2) / (6.0 * (1.0 - s_hat) + (1.0 + 2.0 * s_hat))
                    + ln1
                    - (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3) / (2.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * ln
                    - 1.0 / 9.0 * li2
                    + 4.0 * pow(M_PI, 2) / 27.0
                    - (121.0 - 27.0 * s_hat - 30.0 * s_hat2) / (72.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
                    - (41.0 + 76.0 * s_hat) / (36.0 * (1.0 + 2.0 * s_hat)) * ln1
                    + (14.0 * s_hat3 - 17.0 * s_hat2 - 10.0 * s_hat - 3.0) / (18.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * ln
                    + 17.0 / 18.0 * ln1 * ln
                    - (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3) / (2.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * ln * ln
                    );
        }

        // cf. [HLMW2005], Eq. (100), p. 24
        double omegaem_1010(const double & s_hat) const
        {
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;

            return 2.0 * log(m_b_pole() / m_l) * (
                    - (1.0 + 4.0 * s_hat - 8.0 * s_hat2) / (6.0 * (1.0 - s_hat) + (1.0 + 2.0 * s_hat))
                    + ln1
                    - (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3) / (2.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * ln
                    );
        }

        // cf. [HLMW2005], Eq. (101), p. 25
        double omegaem_77(const double & s_hat) const
        {
            return 0.0;
        }

        // cf. [HLMW2005], Eq. (102), p. 25
        double omegaem_79(const double & s_hat) const
        {
            return 0.0;
        }

        // cf. [HLMW2005], Eq. (103), p. 25
        double omegaem_29(const double & s_hat) const
        {
            return 0.0;
        }

        // cf. [HLMW2005], Eq. (104), p. 25
        double omegaem_22(const double & s_hat) const
        {
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat, s_hat4 = s_hat2 * s_hat2;
            double sigma_1 = 23.787 - 120.948 * s_hat + 365.373 * s_hat2 - 584.206 * s_hat3;
            double sigma_2 = 11.488 - 36.987 * s_hat + 255.330 * s_hat2 - 812.388 * s_hat3 + 1011.791 * s_hat4;

            return 2.0 * log(m_b_pole() / m_l) * (
                        sigma_2 / (8.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat))
                        + sigma_1 / (9.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * log(mu / 5.0)
                        + 64.0 / 81.0 * omegaem_1010(s_hat) * pow(log(mu / 5.0), 2)
                    );
        }

        // cf. [HLMW2005], Eq. (105), p. 25
        double omegaem_27(const double & s_hat) const
        {
            return 0.0;
        }

        // cf. [HLMW2005], Eq. (126), p. 28
        complex<double> g(const double & y) const
        {
            if (y < 1)
                throw InternalError("[HLMW2005] g(y) only valid for y >= 1");

            return 20.0 / 27.0 + 4.0 / 9.0 * y - 2.0 / 9.0 * (2.0 + y) * sqrt(y - 1) * 2.0 * std::atan(1.0 / sqrt(y - 1));
        }

        // cf. [HLMW2005], Eq. (72), p. 17
        // i = { 1 .. 10, Q3 .. Q6, b }
        complex<double> f(unsigned i, const double & s_hat) const
        {
            if ((7 == i) || (8 == i))
                throw InternalError("[HLMW2005] f_i not defined for i = 7,8!");

            static const std::vector<double> rho_b = {
                /* 1 .. 6 */
                0.0, 0.0, -7.0 / 2.0, -2.0 / 3.0, -38.0, -32.0 / 3.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                7.0 / 6.0, 2.0 / 9.0, 38.0 / 3.0, 32.0 / 9.0, -2.0
            };

            static const std::vector<double> rho_c = {
                /* 1 .. 6 */
                4.0 / 3.0, 1.0, 6.0, 0.0, 60.0, 0.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                4.0, 0.0, 40.0, 0.0, 0.0
            };

            static const std::vector<double> rho_0 = {
                /* 1 .. 6 */
                0.0, 0.0, 2.0/9.0, 8.0 / 27.0, 32.0 / 9.0, 128.0 / 27.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                -74.0 / 27.0, -8.0 / 81.0, -752.0 / 27.0, -128.0 / 81.0, 0.0
            };

            static const std::vector<double> rho_sharp = {
                /* 1 .. 6 */
                -16.0 / 27.0, -4.0 / 9.0, 2.0 / 27.0, 8.0 / 81.0, -136.0 / 27.0, 320.0 / 81.0,
                /* 7 .. 10 */
                0.0, 0.0, 0.0, 0.0,
                /* Q3 .. Q6 */
                358.0 / 81.0, -8.0 / 243.0, 1144.0 / 81.0, -320.0 / 243.0, 26.0 / 27.0
            };

            static const std::vector<double> gamma9 = {
                /* 1 .. 6 */
                -32.0 / 27.0, -8.0 / 9.0, -16.0 / 9.0, 32.0 / 27.0, -112.0 / 9.0, 512.0 / 27.0,
                /* 7 .. 10 */
                0.0, 0.0, 8.0, -4.0,
                /* Q3 .. Q6 */
                -272.0 / 27.0, -32.0 / 81.0, 2768.0 / 27.0, -512.0 / 81.0, 16.0 / 9.0
            };

            double m_b = m_b_pole(), m_c = m_c_pole();
            double s = s_hat * pow(m_b, 2);

            complex<double> g_b = g(4.0 / s_hat);
            complex<double> g_c = g(4.0 * pow(m_c, 2) / s);

            /* mu == mu in MSbar scheme, so use m_b_MSbar here */
            return gamma9[i-1] * log(m_b_MSbar / mu)
                + rho_c[i-1] * (g_c + 8.0 / 9.0 * log(m_b / m_c))
                + rho_b[i-1] * g_b
                + rho_0[i-1] * complex<double>(log(s_hat), -M_PI)
                + rho_sharp[i-1];
        }

        complex<double> f9pen(const double & s_hat) const
        {
            complex<double> g_tau = g(4.0 * pow(m_tau / m_b_pole(), 2) / s_hat);

            return 8.0 * log(m_b_MSbar / mu)
                - 3.0 * (g_tau + 8.0 / 9.0 * log(m_b_MSbar / m_tau))
                + 8.0 / 3.0 * std::complex<double>(log(s_hat), -M_PI)
                - 40.0 / 9.0;
        }

        // cf. [HLMW2005], Eq. (132), p. 29
        complex<double> F(const double & s_hat) const
        {
            double r = s_hat * pow(m_b_pole() / m_c_pole(), 2) / 4.0;
            complex<double> result;

            if ((r < 0) || (r > 1))
                throw InternalError("[HLMW] F only implemented for 0 < r < 1");

            result = 3.0 / (2.0 * r) * (std::atan(std::sqrt(r / (1.0 - r))) / std::sqrt(r * (1.0 - r)) - 1.0);

            return result;
        }

        // cf. [HLMW2005], Eq. (6), p. 4
        double branching_ratio(const double & s) const
        {
            static const double alpha_e = 1.0/133.0;
            double m_c = m_c_pole(), m_b = m_b_pole();
            double s_hat = s / pow(m_b, 2), s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
            double lambda_2_hat = lambda_2 / pow(m_b, 2);
            double alpha_s = QCD::alpha_s(mu);
            double kappa = alpha_e / alpha_s, alpha_s_tilde = alpha_s / (4.0 * M_PI);

            static const double u1 = (4.0 * M_PI * M_PI - 25.0) / 12.0;
            double u2 = 27.1 + 23.0 / 3.0 * u1 * log(mu / m_b);
            double uem = 12.0 / 23.0 * (QCD::alpha_s(m_Z) / alpha_s - 1.0);

            // cf. [HLMW2005], Eq. (69), p. 16
            double c7eff = c7() - c3() / 3.0 - 4.0 * c4() / 9.0 - 20.0 * c5() / 3.0 - 80.0 * c6() / 9.0;

            /* S_{AB} */
            // cf. [HLMW2005], Eqs. (112)-(115), p. 26
            double s77 = pow(1.0 - s_hat, 2) * (4.0 + 8.0 / s_hat) * (
                    1.0
                    + 8.0 * alpha_s_tilde * (omega1_77(s_hat) + u1)
                    + kappa * uem
                    + 8.0 * alpha_s_tilde * kappa * omegaem_77(s_hat)
                ) + 24.0 * lambda_2_hat * (2.0 * s_hat2 - 3.0);

            double s79 = 12.0 * pow(1.0 - s_hat, 2) * (
                    1.0
                    + 8.0 * alpha_s_tilde * (omega1_79(s_hat) + u1)
                    + kappa * uem
                    + 8.0 * alpha_s_tilde * omegaem_79(s_hat)
                ) + 24.0 * lambda_2_hat * (1.0 - 6.0 * s_hat + 4.0 * s_hat2);

            double s99 = pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat) * (
                    1.0
                    + 8.0 * alpha_s_tilde * (omega1_99(s_hat) + u1)
                    + kappa * uem
                    + 8.0 * alpha_s_tilde * kappa * omegaem_99(s_hat)
                    + 16.0 * pow(alpha_s_tilde, 2) * (omega2_99(s_hat) + u2 + 4.0 * u1 * omega1_99(s_hat))
                ) + 6.0 * lambda_2_hat * (1.0 - 6.0 * s_hat2 + 4.0 * s_hat3);

            double s1010 = s99
                + 8.0 * alpha_s_tilde * kappa * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat) * (omegaem_1010(s_hat) - omegaem_99(s_hat));

            /* Wilson coefficients */
            std::vector<double> wc = {
                c1, c2, c3, c4, c5, c6, c7eff, c8,
                // We use a different basis of operators: O_{9,10} = alpha_e_tilde * P_{9,10} */
                alpha_s_tilde * kappa * c9,
                alpha_s_tilde * kappa * c10,
                // cf. [HLMW2005], Table 3, p. 17. Using values at mu = 5.0 GeV
                alpha_s_tilde * kappa * -3.72e-2,
                alpha_s_tilde * kappa * -1.04e-2,
                alpha_s_tilde * kappa * -1.71e-6,
                alpha_s_tilde * kappa * -1.03e-3,
                0.0
            };

            /* Corrections, cf. [HLMW2005], Table 6, p. 18 */
            std::vector<complex<double>> m7 = {
                -pow(alpha_s_tilde, 2) * kappa * CharmLoops::F17_massive(mu, s, m_b, m_c),
                -pow(alpha_s_tilde, 2) * kappa * CharmLoops::F27_massive(mu, s, m_b, m_c),
                0.0,
                0.0,
                0.0,
                0.0,
                alpha_s_tilde * kappa,
                -pow(alpha_s_tilde, 2) * kappa * F87(s_hat),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            };

            std::vector<complex<double>> m9 = {
                alpha_s_tilde * kappa * f(1, s_hat) - pow(alpha_s_tilde, 2) * kappa * CharmLoops::F19_massive(mu, s, m_b, m_c),
                alpha_s_tilde * kappa * f(2, s_hat) - pow(alpha_s_tilde, 2) * kappa * CharmLoops::F29_massive(mu, s, m_b, m_c),
                alpha_s_tilde * kappa * f(3, s_hat),
                alpha_s_tilde * kappa * f(4, s_hat),
                alpha_s_tilde * kappa * f(5, s_hat),
                alpha_s_tilde * kappa * f(6, s_hat),
                0.0,
                -pow(alpha_s_tilde, 2) * kappa * F89(s_hat),
                1.0 + alpha_s_tilde * kappa * f9pen(s_hat),
                0.0,
                alpha_s_tilde * kappa * f(3, s_hat),
                alpha_s_tilde * kappa * f(4, s_hat),
                alpha_s_tilde * kappa * f(5, s_hat),
                alpha_s_tilde * kappa * f(6, s_hat),
            };

            std::vector<complex<double>> m10(14, 0.0);
            m10[9] = 1.0; // M^10_i = delta_{10,i}

            // cf. [HLMW2005], Eq. (111)
            double ratio_phi = 0.0;
            for (unsigned i(0) ; i < 14 ; ++i)
            {
                /* diagonal */
                ratio_phi += pow(wc[i], 2) * real(
                        s77 * norm(m7[i])
                        + s99 * norm(m9[i])
                        + s1010 * norm(m10[i])
                        + s79 * m7[i] * conj(m9[i])
                    );

                /* upper */
                for (unsigned j(i + 1) ; j < 14 ; ++j)
                {
                    ratio_phi += wc[i] * wc[j] * real(
                            2.0 * s77 * m7[i] * conj(m7[j])
                            + 2.0 * s99 * m9[i] * conj(m9[j])
                            + 2.0 * s1010 * m10[i] * conj(m10[j])
                            + s79 * (m7[i] * conj(m9[j]) + m9[i] * conj(m7[j]))
                        );
                }
            }

            /* bremsstrahlung */
            static const double c_tau1 = 1.0 / 27.0;
            static const double c_tau2 = - 2.0 / 9.0;
            double z = pow(m_c / m_b, 2);
            double itau_22 = real(Bremsstrahlung::itau_22(s_hat, z));
            double itau_27 = real(Bremsstrahlung::itau_27(s_hat, z));
            double itau_28 = real(Bremsstrahlung::itau_28(s_hat, z));
            double itau_29 = real(Bremsstrahlung::itau_29(s_hat, z));
            double b11 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_22 * c_tau1;
            double b12 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_22 * c_tau2 * 2.0;
            double b22 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_22 * QCD::casimir_f;
            double b17 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_27 * c_tau2 * 2.0;
            double b27 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_27 * QCD::casimir_f * 2.0;
            double b18 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_28 * c_tau2 * 2.0;
            double b28 = pow(alpha_s_tilde, 3) * pow(kappa, 2) * itau_28 * QCD::casimir_f * 2.0;
            double b19 = pow(alpha_s_tilde, 2) * kappa * itau_29 * c_tau2 * 2.0;
            double b29 = pow(alpha_s_tilde, 2) * kappa * itau_29 * QCD::casimir_f * 2.0;
            double add = 0.0;
            add += wc[0] * wc[0] * b11 + wc[0] * wc[1] * b12 + wc[1] * wc[1] * b22;
            add += wc[6] * (wc[0] * b17 + wc[1] * b27);
            add += wc[7] * (wc[0] * b18 + wc[1] * b28);
            add += wc[8] * (wc[0] * b19 + wc[1] * b29);
            ratio_phi += add;

            /* non-perturbative 1/m_c^2 */
            complex<double> cF = F(s_hat);
            double c27 = - pow(alpha_s_tilde * kappa, 2) * 8.0 * lambda_2 / (9.0 * pow(m_c, 2)) * pow(1.0 - s_hat, 2)
                * (1.0 + 6.0 * s_hat - s_hat2) / s_hat * real(cF);
            double c29 = - alpha_s_tilde * kappa * 8.0 * lambda_2 / (9.0 * pow(m_c, 2)) * pow(1.0 - s_hat, 2) * (2.0 + s_hat) * real(cF);
            double c22 = - alpha_s_tilde * kappa * 8.0 * lambda_2 / (9.0 * pow(m_c, 2)) * pow(1.0 - s_hat, 2) * (2.0 + s_hat) * real(cF * conj(m9[1]));
            ratio_phi += c22 * (-2.0 / 9.0 * wc[0] * wc[0] + 7.0 / 6.0 * wc[0] * wc[1] + wc[1] * wc[1]);
            ratio_phi += c27 * (-1.0 / 6.0 * wc[0] + wc[1]) * wc[6];
            ratio_phi += c29 * (-1.0 / 6.0 * wc[0] + wc[1]) * wc[8];

            /* log enhanced em */
            double e22 = 8.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat) * pow(alpha_s_tilde * kappa, 3) * omegaem_22(s_hat);
            double e27 = 96.0 * pow(1.0 - s_hat, 2) * pow(alpha_s_tilde * kappa, 3) * omegaem_27(s_hat);
            double e29 = 8.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat) * pow(alpha_s_tilde * kappa, 2) * omegaem_29(s_hat);
            ratio_phi += e22 * (16.0 / 9.0 * wc[0] * wc[0] + 8.0 / 3.0 * wc[0] * wc[1] + wc[1] * wc[1]);
            ratio_phi += e27 * (4.0 / 3.0 * wc[0] + wc[1]) * wc[6];
            ratio_phi += e29 * (4.0 / 3.0 * wc[0] + wc[1]) * wc[8];

            return br_clnu * ckm * 4.0 / C * ratio_phi;
        }
    };

    BToXsDilepton<HLMW2005>::BToXsDilepton(const Parameters & parameters, const ObservableOptions &) :
        PrivateImplementationPattern<BToXsDilepton<HLMW2005>>(new Implementation<BToXsDilepton<HLMW2005>>(parameters))
    {
    }

    BToXsDilepton<HLMW2005>::~BToXsDilepton()
    {
    }

    double
    BToXsDilepton<HLMW2005>::differential_branching_ratio(const double & s) const
    {
        return _imp->branching_ratio(s) / pow(_imp->m_b_pole(), 2);
    }

    double
    BToXsDilepton<HLMW2005>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrate(std::function<double (const double &)>(
                    std::bind(&BToXsDilepton<HLMW2005>::differential_branching_ratio, this, std::placeholders::_1)),
                128, s_min, s_max);
    }
}
