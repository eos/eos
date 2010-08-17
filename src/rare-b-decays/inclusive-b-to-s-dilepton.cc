/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/bremsstrahlung.hh>
#include <src/rare-b-decays/charm-loops.hh>
#include <src/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/memoise.hh>
#include <src/utils/model.hh>
#include <src/utils/observable.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <functional>
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
        std::shared_ptr<Model> model;

        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter abs_c7;

        Parameter arg_c7;

        Parameter c8;

        Parameter abs_c9;

        Parameter arg_c9;

        Parameter abs_c10;

        Parameter arg_c10;

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
            model(new StandardModel(p)),
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            abs_c7(p["Abs{c7}"]),
            arg_c7(p["Arg{c7}"]),
            c8(p["c8"]),
            abs_c9(p["Abs{c9}"]),
            arg_c9(p["Arg{c9}"]),
            abs_c10(p["Abs{c10}"]),
            arg_c10(p["Arg{c10}"]),
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
            return model->m_b_pole();
        }

        double m_c_pole() const
        {
            return 1.4;
        }

        complex<double> c7() const
        {
            return abs_c7() * complex<double>(cos(arg_c7), sin(arg_c7));
        }

        complex<double> c9() const
        {
            return abs_c9() * complex<double>(cos(arg_c9), sin(arg_c9));
        }

        complex<double> c10() const
        {
            return abs_c10() * complex<double>(cos(arg_c10), sin(arg_c10));
        }

        double s_hat(const double & s) const
        {
            double m_b = m_b_pole();

            return s / m_b / m_b;
        }

        /* NLO functions */
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
                    - (1.0 + 4.0 * s_hat - 8.0 * s_hat2) / (6.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
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
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat;

            return 2.0 * log(m_b_pole() / m_l) * (
                    s_hat / (2.0 * (1.0 - s_hat) * (2.0 + s_hat))
                    + ln1
                    - s_hat * (-3.0 + 2.0 * s_hat2) / (2.0 * pow(1.0 - s_hat, 2) * (2.0 + s_hat)) * ln
                    );
        }

        // cf. [HLMW2005], Eq. (102), p. 25
        double omegaem_79(const double & s_hat) const
        {
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat * s_hat;

            return 2.0 * log(m_b_pole() / m_l) * (
                    -1.0 / (2.0 * (1.0 - s_hat))
                    + ln1
                    + (-1.0 + 2.0 * s_hat - 2.0 * s_hat2) / (2.0 * pow(1.0 - s_hat, 2)) * ln
                    );
        }

        // cf. [HLMW2005], Eq. (103), p. 25
        complex<double> omegaem_29(const double & s_hat) const
        {
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
            double sigma_1 = 23.787 - 120.948 * s_hat + 365.373 * s_hat2 - 584.206 * s_hat3;
            double sigma_1_I = 1.653 + 6.009 * s_hat - 17.080 * s_hat2 + 115.880 * s_hat3;

            return 2.0 * log(m_b_pole() / m_l) * (
                    complex<double>(sigma_1, sigma_1_I) / (8.0 * (1.0 - s_hat) * (1.0 - s_hat) * (1.0 + 2.0 * s_hat))
                    + 16.0 / 9.0 * omegaem_1010(s_hat) * log(mu / 5.0)
                    );
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
        complex<double> omegaem_27(const double & s_hat) const
        {
            double s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
            double sigma_3 = 109.311 - 846.039 * s_hat + 2890.115 * s_hat2 - 4179.072 * s_hat3;
            double sigma_3_I = 4.606 + 17.650 * s_hat - 53.244 * s_hat2 + 348.069 * s_hat3;

            return 2.0 * log(m_b_pole() / m_l) * (
                    complex<double>(sigma_3, sigma_3_I) / (96.0 * (1.0 - s_hat) * (1.0 - s_hat))
                    + 8.0 / 9.0 * omegaem_79(s_hat) * log(mu / 5.0)
                    );
        }

        // cf. [HLMW2005], Eq. (126), p. 28
        complex<double> g(const double & y) const
        {
            complex<double> x;

            if (y < 1)
                x = complex<double>(log(abs(1.0 + sqrt(1.0 - y)) / (1.0 - sqrt(1.0 - y))), -M_PI);
            else
                x = 2.0 * atan(1.0 / sqrt(y - 1.0));

            return 20.0 / 27.0 + 4.0 / 9.0 * y - 2.0 / 9.0 * (2.0 + y) * sqrt(abs(y - 1)) * x;
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
            double x = sqrt(1.0 - 1.0 / r);
            complex<double> result;

            if ((0 < r) && (r < 1))
                result = 3.0 / (2.0 * r) * (std::atan(std::sqrt(r / (1.0 - r))) / std::sqrt(r * (1.0 - r)) - 1.0);
            else
                result = 3.0 / (2.0 * r) * (
                    complex<double>(log((1.0 - x) / (1.0 + x)), M_PI) / (2.0 * sqrt(r * (r - 1.0)))
                    - 1.0);

            return result;
        }

        // cf. [HLMW2005], Eq. (6), p. 4
        double branching_ratio(const double & s) const
        {
            static const double alpha_e = 1.0/133.0;
            double m_c = m_c_pole(), m_b = m_b_pole();
            double s_hat = s / pow(m_b, 2), s_hat2 = s_hat * s_hat, s_hat3 = s_hat2 * s_hat;
            double lambda_2_hat = lambda_2 / pow(m_b, 2);
            double alpha_s = model->alpha_s(mu);
            double kappa = alpha_e / alpha_s, alpha_s_tilde = alpha_s / (4.0 * M_PI);

            static const double u1 = (4.0 * M_PI * M_PI - 25.0) / 12.0;
            double u2 = 27.1 + 23.0 / 3.0 * u1 * log(mu / m_b);
            double uem = 12.0 / 23.0 * (model->alpha_s(m_Z) / alpha_s - 1.0);

            // cf. [HLMW2005], Eq. (69), p. 16
            complex<double> c7eff = c7() - c3() / 3.0 - 4.0 * c4() / 9.0 - 20.0 * c5() / 3.0 - 80.0 * c6() / 9.0;

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
                    + 8.0 * alpha_s_tilde * kappa * omegaem_79(s_hat)
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
            std::vector<complex<double>> wc = {
                c1(), c2(), c3(), c4(), c5(), c6(), c7eff, c8(),
                // We use a different basis of operators: O_{9,10} = alpha_e_tilde * P_{9,10} */
                alpha_s_tilde * kappa * c9(),
                alpha_s_tilde * kappa * c10(),
                // cf. [HLMW2005], Table 3, p. 17. Using values at mu = 5.0 GeV
                alpha_s_tilde * kappa * -3.72e-2,
                alpha_s_tilde * kappa * -1.04e-2,
                alpha_s_tilde * kappa * -1.71e-6,
                alpha_s_tilde * kappa * -1.03e-3,
                0.0
            };

            /* Corrections, cf. [HLMW2005], Table 6, p. 18 */
            std::vector<complex<double>> m7 = {
                -pow(alpha_s_tilde, 2) * kappa * memoise(CharmLoops::F17_massive, mu(), s, m_b, m_c),
                -pow(alpha_s_tilde, 2) * kappa * memoise(CharmLoops::F27_massive, mu(), s, m_b, m_c),
                0.0,
                0.0,
                0.0,
                0.0,
                alpha_s_tilde * kappa,
                -pow(alpha_s_tilde, 2) * kappa * CharmLoops::F87_massless(mu, s, m_b),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            };

            std::vector<complex<double>> m9 = {
                alpha_s_tilde * kappa * f(1, s_hat) - pow(alpha_s_tilde, 2) * kappa * memoise(CharmLoops::F19_massive, mu(), s, m_b, m_c),
                alpha_s_tilde * kappa * f(2, s_hat) - pow(alpha_s_tilde, 2) * kappa * memoise(CharmLoops::F29_massive, mu(), s, m_b, m_c),
                alpha_s_tilde * kappa * f(3, s_hat),
                alpha_s_tilde * kappa * f(4, s_hat),
                alpha_s_tilde * kappa * f(5, s_hat),
                alpha_s_tilde * kappa * f(6, s_hat),
                0.0,
                -pow(alpha_s_tilde, 2) * kappa * CharmLoops::F89_massless(mu, s, m_b),
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
                ratio_phi += norm(wc[i]) * real(
                        s77 * norm(m7[i])
                        + s99 * norm(m9[i])
                        + s1010 * norm(m10[i])
                        + s79 * m7[i] * conj(m9[i])
                    );

                /* upper */
                for (unsigned j(i + 1) ; j < 14 ; ++j)
                {
                    ratio_phi += real(wc[i] * conj(wc[j]) * (
                            2.0 * s77 * m7[i] * conj(m7[j])
                            + 2.0 * s99 * m9[i] * conj(m9[j])
                            + 2.0 * s1010 * m10[i] * conj(m10[j])
                            + s79 * (m7[i] * conj(m9[j]) + m9[i] * conj(m7[j])))
                        );
                }
            }

            /* bremsstrahlung */
            static const double c_tau1 = 1.0 / 27.0;
            static const double c_tau2 = - 2.0 / 9.0;
            double z = pow(m_c / m_b, 2);
            double itau_22 = real(memoise(Bremsstrahlung::itau_22, s_hat, z));
            double itau_27 = real(memoise(Bremsstrahlung::itau_27, s_hat, z));
            double itau_28 = real(memoise(Bremsstrahlung::itau_28, s_hat, z));
            double itau_29 = real(memoise(Bremsstrahlung::itau_29, s_hat, z));
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
            add += norm(wc[0]) * b11 + real(wc[0] * conj(wc[1])) * b12 + norm(wc[1]) * b22;
            add += real(conj(wc[6]) * (wc[0] * b17 + wc[1] * b27));
            add += real(conj(wc[7]) * (wc[0] * b18 + wc[1] * b28));
            add += real(conj(wc[8]) * (wc[0] * b19 + wc[1] * b29));
            ratio_phi += add;

            /* non-perturbative 1/m_c^2 */
            complex<double> cF = F(s_hat);
            double c27 = - pow(alpha_s_tilde * kappa, 2) * 8.0 * lambda_2 / (9.0 * pow(m_c, 2)) * pow(1.0 - s_hat, 2)
                * (1.0 + 6.0 * s_hat - s_hat2) / s_hat * real(cF);
            double c29 = - alpha_s_tilde * kappa * 8.0 * lambda_2 / (9.0 * pow(m_c, 2)) * pow(1.0 - s_hat, 2) * (2.0 + s_hat) * real(cF);
            double c22 = - alpha_s_tilde * kappa * 8.0 * lambda_2 / (9.0 * pow(m_c, 2)) * pow(1.0 - s_hat, 2) * (2.0 + s_hat) * real(cF * conj(m9[1]));
            ratio_phi += c22 * (-2.0 / 9.0 * norm(wc[0]) + 7.0 / 6.0 * real(wc[0] * conj(wc[1])) + norm(wc[1]));
            ratio_phi += c27 * real((-1.0 / 6.0 * wc[0] + wc[1]) * conj(wc[6]));
            ratio_phi += c29 * real((-1.0 / 6.0 * wc[0] + wc[1]) * conj(wc[8]));

            /* log enhanced em */
            double e22 = 8.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat) * pow(alpha_s_tilde * kappa, 3) * omegaem_22(s_hat);
            complex<double> e27 = 96.0 * pow(1.0 - s_hat, 2) * pow(alpha_s_tilde * kappa, 3) * omegaem_27(s_hat);
            complex<double> e29 = 8.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat) * pow(alpha_s_tilde * kappa, 2) * omegaem_29(s_hat);
            ratio_phi += e22 * (16.0 / 9.0 * norm(wc[0]) + 8.0 / 3.0 * real(wc[0] * conj(wc[1])) + norm(wc[1]));
            ratio_phi += real(e27 * (4.0 / 3.0 * wc[0] + wc[1]) * conj(wc[6]));
            ratio_phi += real(e29 * (4.0 / 3.0 * wc[0] + wc[1]) * conj(wc[8]));

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
