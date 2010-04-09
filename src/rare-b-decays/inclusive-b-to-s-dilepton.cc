/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <iostream>

namespace wf
{
    /* GN1997 */

    template <>
    struct Implementation<BToXsDilepton<GN1997>>
    {
        Parameter c1;

        Parameter c2;

        Parameter c3;

        Parameter c4;

        Parameter c5;

        Parameter c6;

        Parameter c7;

        Parameter c7prime;

        Parameter c9;

        Parameter c9prime;

        Parameter c10;

        Parameter c10prime;

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_s;

        Parameter br_bcsl;

        double m_l;

        Implementation(const Parameters & p) :
            c1(p["c1"]),
            c2(p["c2"]),
            c3(p["c3"]),
            c4(p["c4"]),
            c5(p["c5"]),
            c6(p["c6"]),
            c7(p["c7"]),
            c9(p["c9"]),
            c10(p["c10"]),
            c7prime(p["c7prime"]),
            c9prime(p["c9prime"]),
            c10prime(p["c10prime"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            m_c(p["mass::c"]),
            m_s(p["mass::s"]),
            br_bcsl(p["exp::BR(B->X_clnu)"])
        {
            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.0;//0.10565836; // (GeV), cf. [PDG2008], p. 13
        }

        double s_hat(const double & s) const
        {
            double m_b_pole = QCD::mb_pole(m_b_MSbar);

            return s / m_b_pole / m_b_pole;
        }

        /* Effective wilson coefficients */
        // cf. [GN19997], Eq. (3.7), p. 8
        complex<double> c9eff(const double & s_hat) const
        {
            return c9() + (3.0 * c1() + c2()) * g(s_hat);
        }

        complex<double> c9primeeff(const double & s_hat) const
        {
            return complex<double>(c9prime(), 0.0);
        }

        // cf. [BM1995], Eqs. (3.9), (2.29)
        complex<double> g(const double & s_hat) const
        {
            const double m_b_pole = QCD::mb_pole(m_b_MSbar);
            const double m_c_hat = m_c / m_b_pole;
            const double x = 4.0 * m_c_hat * m_c_hat / s_hat;
            double y = std::sqrt(std::abs(1.0 - x));
            complex<double> z;

            if (x <= 1.0)
                z = complex<double>(std::log(std::abs((y + 1.0) / (y - 1.0))), -M_PI);
            else
                z = complex<double>(2.0 * std::atan(1.0 / y), 0.0);

            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]

            // k (no unit) and mass, partial decay width, total decay width (all in GeV)
            std::vector<std::tuple<double, double, double, double>> resonances =
            {
                // J/psi(1S)
                std::make_tuple(2.3, 3.096, 5.55e-6, 93.2e-6),
                // psi(2S)
                std::make_tuple(1.0, 3.686, 2.38e-6, 317e-6),
                // psi(3770)
                std::make_tuple(1.0, 3.772, 0.27e-6, 27.3e-6),
                // psi(4040)
                std::make_tuple(1.0, 4.040, 0.86e-6, 80e-6),
            };

            complex<double> resonance_parts = complex<double>(0.0, 0.0);
            for (auto r(resonances.begin()) ; r != resonances.end() ; ++r)
            {
                double k = std::get<0>(*r);
                double m_hat = std::get<1>(*r) / m_b_pole;
                double m_hat2 = m_hat * m_hat;
                double gamma_ee_hat = std::get<2>(*r) / m_b_pole;
                double gamma_hat = std::get<3>(*r) / m_b_pole;
                double denom = (s_hat - m_hat2) * (s_hat - m_hat2) + m_hat2 * gamma_hat * gamma_hat;

                resonance_parts = resonance_parts
                    + (k / denom) * complex<double>((s_hat - m_hat2) * m_hat * gamma_ee_hat, -m_hat * m_hat * gamma_ee_hat * gamma_hat);
            }

            return -1.0 / 3.0 * (4.0 * std::log(m_c_hat) - 16.0 / 3.0 - 2.0 * x
                + (2.0 + x) * y * z) - (3.0 * M_PI / alpha_e / alpha_e) * resonance_parts;
        }

        double m1(const double & u_hat, const double & s_hat, const double & m_s_hat2, const double & m_l_hat2) const
        {
            return 32.0/3.0 * u_hat / s_hat * (2.0 * (1.0 + m_s_hat2) * std::pow(1.0 - - m_s_hat2, 2.0)
                    - s_hat * (1.0 + 14.0 * m_s_hat2 + std::pow(m_s_hat2, 2.0))
                    - std::pow(s_hat, 2.0) * (1.0 + m_s_hat2));
        }

        double m2(const double & u_hat, const double & s_hat, const double & m_s_hat2) const
        {
            return 128.0/3.0 * u_hat / s_hat * std::sqrt(m_s_hat2) * (2.0 * std::pow(1.0 - m_s_hat2, 2.0)
                    - s_hat * (4.0 + 4.0 * m_s_hat2 + s_hat));
        }

        double m3(const double & u_hat, const double & s_hat, const double & m_s_hat2) const
        {
            return 32.0 * u_hat * (std::pow(1.0 - m_s_hat2, 2.0) - s_hat * (1.0 + m_s_hat2));
        }

        double m4(const double & u_hat, const double & s_hat, const double & m_s_hat2, const double & m_l_hat2) const
        {
            return 8.0/3.0 * u_hat * (std::pow(1.0 - m_s_hat2, 2.0) + s_hat * (1.0 + m_s_hat2)
                    - 2.0 * std::pow(s_hat, 2.0) + std::pow(u_hat, 2.0) * 2.0 * m_l_hat2 / s_hat);
        }

        double m5(const double & u_hat, const double & s_hat, const double & m_s_hat2) const
        {
            return -32.0 * u_hat * std::sqrt(m_s_hat2) * s_hat;
        }

        /* We compute the branching ratio similar to [GN1997], Eq. (3.9) and
         * [FKMY1998], Eq. (21):
         *
         *  dBR/dsHat = B_0 / (2 m_b^6) * Sigma(sHat),
         *
         *  Sigma(sHat) = W_i({C_j}) M_i(sHat)
         *
         */
        double sigma(const double & s_hat) const
        {
            complex<double> c9_eff(c9eff(s_hat));
            complex<double> c9prime_eff(c9primeeff(s_hat));
            double m_b_pole = QCD::mb_pole(m_b_MSbar);
            double m_l_hat = m_l / m_b_pole, m_l_hat2 = std::pow(m_l_hat, 2.0);
            double m_s_hat = m_s / m_b_pole, m_s_hat2 = std::pow(m_s_hat, 2.0);
            double u_hat = lambda(1.0, m_s_hat, s_hat);

            // cf. [GN1997], Eq,. (3.10), p. 9-10. The s dependent terms and the
            // F_i have been absorbed into M_i for i = 1, 2, 3, 4. M_5 ~ -3
            // m_s * s.
            double w1 = std::pow(c7, 2.0) + std::pow(c7prime, 2.0);
            double w2 = c7 * c7prime;
            double w3 = real(c7() * c9_eff + c7prime() * c9prime_eff);
            double w4 = norm(c9_eff) + norm(c9prime_eff) + std::pow(c10, 2.0) + std::pow(c10prime, 2.0);
            double w5 = real(2.0 * c7() * c9prime_eff + 2.0 * conj(c9_eff) * c7prime() + conj(c9_eff) * c9prime_eff);

            // At the time being we do not consider scalar or pseudo-scalar wilson coefficients, so no additional term ~f5.
            double result = w1 * m1(u_hat, s_hat, m_s_hat2, m_l_hat2)
                + w2 * m2(u_hat, s_hat, m_s_hat2)
                + w3 * m3(u_hat, s_hat, m_s_hat2)
                + w4 * m4(u_hat, s_hat, m_s_hat2, m_l_hat2)
                + w5 * m5(u_hat, s_hat, m_s_hat2);
            result *= (1.0 + 2.0 * m_l_hat2 / s_hat);

            return result;
        }

        double norm2() const
        {
            static const double alpha_e = 1.0/133.0;
            static const double ckm = 1.0; // |V_tb V_ts^*|/|V_cb| = 1.0 + O(lambda^4) according to Wolfenstein

            double m_b_pole = QCD::mb_pole(m_b_MSbar);
            double m_c_hat = m_c / m_b_pole;
            double f = 1.0 - 8.0 * std::pow(m_c_hat, 2.0) + 8.0 * std::pow(m_c_hat, 6.0)
                - std::pow(m_c_hat, 8.0) - 24.0 * std::pow(m_c_hat, 4.0) * std::log(m_c_hat);
            double kappa = 1.0 - 2.0 * QCD::alpha_s(m_b_MSbar) / (3.0 * M_PI)
                * ((M_PI * M_PI - 32.0/4.0) * std::pow(1.0 - m_c_hat, 2.0) + 1.5);

            return br_bcsl * alpha_e * alpha_e * 3.0 / (32.0 * M_PI * M_PI * kappa * f) / std::pow(m_b_pole, 2.0);
        }
    };

    BToXsDilepton<GN1997>::BToXsDilepton(const Parameters & parameters, const ObservableOptions &) :
        PrivateImplementationPattern<BToXsDilepton<GN1997>>(new Implementation<BToXsDilepton<GN1997>>(parameters))
    {
    }

    BToXsDilepton<GN1997>::~BToXsDilepton()
    {
    }

    double
    BToXsDilepton<GN1997>::differential_branching_ratio(const double & s) const
    {
        return _imp->norm2() * _imp->sigma(_imp->s_hat(s));
    }

    double
    BToXsDilepton<GN1997>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrate(std::tr1::function<double (const double &)>(
                    std::tr1::bind(&BToXsDilepton<GN1997>::differential_branching_ratio, this, std::tr1::placeholders::_1)),
                400, s_min, s_max);
    }
}
