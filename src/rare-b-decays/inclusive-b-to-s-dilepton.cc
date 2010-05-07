/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <iostream>

namespace wf
{
    /* ALGH2001 */

    template <>
    struct Implementation<BToXsDilepton<ALGH2001>>
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

        Parameter br_bcsl;

        Parameter lambda_1;

        Parameter lambda_2;

        Parameter mu;

        double m_l;

        Implementation(const Parameters & p) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["c7"]),
            c8(p["c8"]),
            c9(p["c9"]),
            c10(p["c10"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c_MSbar(p["mass::c"]),
            m_s(p["mass::s"]),
            br_bcsl(p["exp::BR(B->X_clnu)"]),
            lambda_1(p["B->X_s::lambda_1"]),
            lambda_2(p["B->X_s::lambda_2"]),
            mu(p["mu"])
        {
            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.0;//0.10565836; // (GeV), cf. [PDG2008], p. 13
#if 0
            double sh = 1.0 / pow(m_b_pole(), 2);
            std::cout << "s_hat  = " << sh << std::endl;
            std::cout << "c7eff  = " << c7eff(sh) << std::endl;
            std::cout << "c9eff  = " << c9eff(sh) << std::endl;
            std::cout << "c10eff = " << c10eff(sh) << std::endl;
            std::cout << "dBR    = " << branching_ratio(sh) << std::endl;
            std::cout << "omega_7= " << omega_7(sh) << std::endl;
            std::cout << "omega_9= " << omega_9(sh) << std::endl;
            std::cout << "gc     = " << gc(sh) << std::endl;
            std::cout << "g1     = " << g1(sh) << std::endl;
            std::cout << "g2     = " << g2(sh) << std::endl;
            std::cout << "==================" << std::endl << std::endl;

            sh = 6.0 / pow(m_b_pole(), 2);
            std::cout << "s_hat  = " << sh << std::endl;
            std::cout << "c7eff  = " << c7eff(sh) << std::endl;
            std::cout << "c9eff  = " << c9eff(sh) << std::endl;
            std::cout << "c10eff = " << c10eff(sh) << std::endl;
            std::cout << "dBR    = " << branching_ratio(sh) << std::endl;
            std::cout << "omega_7= " << omega_7(sh) << std::endl;
            std::cout << "omega_9= " << omega_9(sh) << std::endl;
            std::cout << "gc     = " << gc(sh) << std::endl;
            std::cout << "g1     = " << g1(sh) << std::endl;
            std::cout << "g2     = " << g2(sh) << std::endl;
            std::cout << "==================" << std::endl << std::endl;
#endif
        }

        double m_b_pole() const
        {
            double result = QCD::mb_pole(m_b_MSbar);

            return result;
        }

        double m_c_pole() const
        {
            double result = 1.4;

            return result;
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
        // cf. [AAGW2002], Eq. (81), p. 26
        double omega_7(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat2;

            return - 8.0/3.0 * log(mu / m_b_pole())
                - 4.0/3.0 * li2 - 2.0/9.0 * M_PI * M_PI - 2.0/3.0 * ln * ln1
                - (8.0 + s_hat) / (3.0 * (2.0 + s_hat)) * ln1
                - (2.0 * s_hat * (2.0 - 2.0 * s_hat + s_hat2)) / (3.0 * pow(1.0 - s_hat, 2) * (2.0 + s_hat)) * ln
                - (16.0 - 11.0 * s_hat - 17.0 * s_hat2) / (18.0 * (2.0 + s_hat) * (1.0 - s_hat));
        }

        /* Effective wilson coefficients */
        // cf. [ALGH2001], Eq. (13), p. 6
        complex<double> c7eff(const double & s_hat) const
        {
            double m_b = m_b_pole();
            double s = s_hat * m_b;
            double alpha_s = QCD::alpha_s(mu());

            complex<double> a8 = c8() + c3() - 1.0/6.0 * c4() + 20.0 * c5() - 10.0/3.0 * c6();
            complex<double> lo = c7() - 1.0/3.0 * c3() - 4.0/9.0 * c4() - 20.0/3.0 * c5() - 80.0/9.0 * c6();
            complex<double> prefactor = 1.0 + alpha_s / M_PI * omega_7(s_hat);
            complex<double> nlo = -1.0 * alpha_s / (4.0 * M_PI)
                * (c1() * CharmLoops::F17(mu, s, m_b) + c2() * CharmLoops::F27(mu, s, m_b) + a8 * F87(s_hat));

            return prefactor * lo + nlo;
        }

        // cf. [BMU1999], Eq. (34), p. 9
        double omega_9(const double & s_hat) const
        {
            double li2 = gsl_sf_dilog(s_hat);
            double ln = log(s_hat), ln1 = log(1.0 - s_hat);
            double s_hat2 = s_hat2;

            return -4.0/3.0 * li2 - 2.0/3.0 * ln1 * ln - 2.0/9.0 * M_PI * M_PI
                - (5.0 + 4.0 * s_hat) / (3.0 * (1.0 + 2.0 * s_hat)) * ln1
                - (2.0 * s_hat * (1.0 + s_hat) * (1.0 - 2.0 * s_hat)) / (3.0 * pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * ln
                + (5.0 + 9.0 * s_hat - 6.0 * s_hat2) / (6.0 * (1.0 - s_hat) * (1.0 + 2.0 * s_hat));
        }

        // cf. [ALGH2001], Eq. (14), p. 6
        complex<double> c9eff(const double & s_hat) const
        {
            double m_b = m_b_pole();
            double s = s_hat * m_b;
            double alpha_s = QCD::alpha_s(mu());

            complex<double> a8 = c8() + c3() - 1.0/6.0 * c4() + 20.0 * c5() - 10.0/3.0 * c6();
            complex<double> lo = c9 + 4.0/3.0 * c3() + 64.0 / 9.0 * c5() + 64./27.0 * c6(); // + gamma^0_{i9} terms
            complex<double> lo_loops = (4.0/3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5()) * CharmLoops::h(mu, s, m_c_pole())
                - 0.5 * (7.0 * c3() + 4.0/3.0 * c4() + 76.0 * c5() + 64.0/3.0 * c6()) * CharmLoops::h(mu, s, m_b_pole())
                - 0.5 * (c3() + 4.0/3.0 * c4() + 16.0 * c5() + 64.0/3.0 * c6()) * CharmLoops::h(mu, s);
            complex<double> prefactor = 1.0 + alpha_s / M_PI * omega_9(s_hat);
            complex<double> nlo = -1.0 * alpha_s / (4.0 * M_PI)
                * (c1() + CharmLoops::F19(mu, s, m_b) + c2() * CharmLoops::F29(mu, s, m_b) + a8 * F89(s_hat));

            return prefactor * (lo + lo_loops) + nlo;
        }

        // cf. [ALGH2001], Eq. (15), p. 6
        complex<double> c10eff(const double & s_hat) const
        {
            double alpha_s = QCD::alpha_s(mu());

            return (1.0 + alpha_s / M_PI * omega_9(s_hat)) * c10();
        }

        // cf. [ALGH2001], Eq. (29), p. 8
        double g1(const double & s_hat) const
        {
            double m_b2 = pow(m_b_pole(), 2);
            double s_hat2 = pow(s_hat, 2);
            double s_hat3 = pow(s_hat, 3);

            return 1.0
                + lambda_1 / (2.0 * m_b2)
                + 3.0 * (1.0 - 15.0 * s_hat2 + 10.0 * s_hat3) / (pow(1.0 - s_hat, 2) * (1.0 + 2.0 * s_hat)) * lambda_2 / (2.0 * m_b2);
        }

        // cf. [ALGH2001], Eq. (30), p. 8
        double g2(const double & s_hat) const
        {
            double m_b2 = pow(m_b_pole(), 2);
            double s_hat2 = pow(s_hat, 2);
            double s_hat3 = pow(s_hat, 3);

            return 1.0
                + lambda_1 / (2.0 * m_b2)
                - 3.0 * (6.0 + 3.0 * s_hat - 5.0 * s_hat3) / (pow(1.0 - s_hat, 2) * (2.0 + s_hat)) * lambda_2 / (2.0 * m_b2);
        }

        // cf. [ALGH2001], Eq. (30), p. 8
        double g3(const double & s_hat) const
        {
            double m_b2 = pow(m_b_pole(), 2);
            double s_hat2 = pow(s_hat, 2);
            double s_hat3 = pow(s_hat, 3);

            return 1.0
                + lambda_1 / (2.0 * m_b2)
                - (5.0 + 6.0 * s_hat - 7.0 * s_hat2) / pow(1.0 - s_hat, 2) * lambda_2 / (2.0 * m_b2);
        }

        // cf. [ALGH2001], Eq. (32), p. 8
        double gc(const double & s_hat) const
        {
            double m_c2 = pow(m_c_pole(), 2);
            double r = s_hat / (4.0 * m_c2 / pow(m_b_pole(), 2));

            // cf. [ALGH2001], Eq. (33), p. 9
            complex<double> f;
            if (r > 1.0)
            {
                double a = 1.0 / (2.0 * sqrt(r * (r - 1.0)));
                double b = sqrt(1.0 - 1.0 / r);
                f = a * std::complex<double>(log((1.0 - b) / (1.0 + b)), M_PI) - 1.0;
            }
            else
            {
                double a = sqrt(r / (1.0 - r));
                f = a / r * atan(a) - 1.0;
            }
            f *= 1.5 / r;

            return -8.0/9.0 * (c2() - c1() / 6.0) * lambda_2 / m_c2
                * real(conj(f) * (c9eff(s_hat) * (2.0 + s_hat) + c7eff(s_hat) * (1.0 + 6.0 * s_hat - s_hat * s_hat) / s_hat));
        }

        /* We compute the branching ratio as given in [ALGH2001], Eqs. (28),(35) and
         * [BMU1999], Eq. (45):
         *
         *  dBr(B->X_s l^+ l^-)/dsHat = Br(B->X_c l nu) / Gamma(b->X_c l nu) * dGamma(b->X_s l^+ l^-)/dsHat,
         *  dBr/ds = mb^{-2} dBr/dsHat
         *
         */
        double branching_ratio(const double & s_hat) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]

            double z = pow(m_c_pole() / m_b_pole(), 2);
            double z2 = z * z, z3 = z * z2, z4 = z2 * z2, lnz = log(z), ln2z = lnz * lnz;
            double sqrtz = m_c_pole() / m_b_pole();
            double li2z = gsl_sf_dilog(z), li2psqrtz = gsl_sf_dilog(sqrtz), li2msqrtz = gsl_sf_dilog(-sqrtz);

            // cf. [BMU1999], Eq. (46), p. 16
            double g = 1.0 - 8.0 * z + 8.0 * z3 - z4 - 12.0 * z2 * log(z);
            // cf. [BMU1999], Eq. (48), p. 16
            double h = -(1.0 - z2) * (25.0/4.0 - 239.0/3.0 * z + 25.0/4.0 * z2)
                + z * lnz * (20.0 + 90.0 * z - 4.0/3.0 * z2 + 17.0/3.0 * z3)
                + z2 * ln2z * (36.0 + z2)
                + (1 - z2) * (17.0/3.0 - 64.0/3.0 * z + 17.0/3.0 * z2) * log(1 - z)
                - 4.0 * (1.0 + 30.0 * z2 + z4) * lnz * log(1 - z)
                - (1.0 + 16.0 * z2 + z4) * (6.0 * li2z - M_PI * M_PI)
                - 32.0 * pow(z, 1.5) * (1.0 + z) * (M_PI * M_PI - 4.0 * (li2psqrtz - li2msqrtz) - 2.0 * lnz * log((1.0 - sqrtz) / (1.0 + sqrtz)));
            // cf. [BMU1999], Eq. (47), pp. 16-17
            double kappa = 1.0 - 2.0/3.0/M_PI * QCD::alpha_s(m_b_pole()) * h / g;

            double ckm = 1.0; // lambda_t^2/V_cb^2 = 1.0 + O(lambda^4)
            double a = br_bcsl * pow(alpha_e, 2) / (4.0 * M_PI * M_PI) / pow(m_b_pole(), 2) * ckm / g / kappa;
            double b = pow(1.0 - s_hat, 2) * ((1.0 + 2.0 * s_hat) * (norm(c9eff(s_hat)) + norm(c10eff(s_hat))) * g1(s_hat)
                    + 4.0 * (1.0 + 2.0 / s_hat) * norm(c7eff(s_hat)) * g2(s_hat)
                    + 12.0 * real(c7eff(s_hat) * conj(c9eff(s_hat))) * g3(s_hat)
                    + gc(s_hat));

            return a * b;
        }
    };

    BToXsDilepton<ALGH2001>::BToXsDilepton(const Parameters & parameters, const ObservableOptions &) :
        PrivateImplementationPattern<BToXsDilepton<ALGH2001>>(new Implementation<BToXsDilepton<ALGH2001>>(parameters))
    {
    }

    BToXsDilepton<ALGH2001>::~BToXsDilepton()
    {
    }

    double
    BToXsDilepton<ALGH2001>::differential_branching_ratio(const double & s) const
    {
        return _imp->branching_ratio(_imp->s_hat(s));
    }

    double
    BToXsDilepton<ALGH2001>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrate(std::tr1::function<double (const double &)>(
                    std::tr1::bind(&BToXsDilepton<ALGH2001>::differential_branching_ratio, this, std::tr1::placeholders::_1)),
                128, s_min, s_max);
    }
}
