/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/low-recoil/low-recoil.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>

#include <iostream>

namespace wf
{
    // TODO: Move to low_recoil/form_factors ?
    // cf. [ABHH1999], p. 8, Table 3
    struct FormFactors
    {
        static double v(const double & s_hat)
        {
            return 0.457 * std::exp(1.482 * s_hat + 1.014 * s_hat * s_hat);
        }

        static double a_1(const double & s_hat)
        {
            return 0.337 * std::exp(0.602 * s_hat + 0.258 * s_hat * s_hat);
        }

        static double a_2(const double & s_hat)
        {
            return 0.282 * std::exp(1.172 * s_hat + 0.567 * s_hat * s_hat);
        }
    };

    template <>
    struct Implementation<Decay<BToKstarDilepton>>
    {
        double c1;

        double c2;

        double c3;

        double c4;

        double c5;

        double c6;

        double c7;

        double c9;

        double c10;

        double m_b_MSbar;

        double m_c;

        double m_B;

        double m_Kstar;


        double mu;

        Implementation(double mu) :
            c1(-0.248), // TODO
            c2(+1.107), // TODO
            c3(+0.011), // TODO
            c4(-0.026), // TODO
            c5(+0.007), // TODO
            c6(-0.031), // TODO
            c7(-0.313), // TODO: Compute correct value from [GP]
            c9(+4.344), // TODO: Compute correct value from [GP]
            c10(-4.669), // TODO: Compute correct value from [GP]
            m_b_MSbar(4.2), // (GeV), cf. [PDG2006], p. 24
            m_c(1.27), // TODO
            m_B(5.279), // (GeV), cf. [PDG2006], p. 87
            m_Kstar(0.896), // (GeV), cf. [PDG2006], p. 51
            mu(mu)
        {
        }

        // Helpers
        Complex<double> c7eff(double s) const
        {
            // TODO: Neglecting contributions ~alpha_s / 4.0 / M_PI. These do involve spectator scattering,
            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            double c7eff0 = c7 - 4.0/9.0 * c3 - 4.0/3.0 * c4 + 1.0/9.0 * c5 + 1.0/3.0 * c6;
            return Complex<double>::Cartesian(c7eff0, 0); // cf. [GP2004] Eq. (56), p. 10
        }

        Complex<double> c9eff(double s) const
        {
            // For r_i, g_i cf. [GP2004], Eqs. (27)-(29), p. 6
            // TODO: I think all sqrt(r_c) needs to be replaced by r_c. Asked Christoph for checking it.
            double r_b = std::sqrt(4.0 * m_b_MSbar * m_b_MSbar / s - 1.0);
            double r_c = std::sqrt(1.0 - 4.0 * s / m_c / m_c);
            Complex<double> g_0 = Complex<double>::Cartesian(std::log(s / mu / mu), -M_PI) * (1.0 / 6.0) - 5.0 / 18.0;
            double g_m_b = std::log(m_b_MSbar * m_b_MSbar / mu / mu) / 6.0 - 5.0 / 18.0 - 2.0 * m_b_MSbar * m_b_MSbar / 3.0 / s
                + r_b / 3.0 * (1.0 + 2.0 * m_b_MSbar * m_b_MSbar / s) * std::atan(1.0 / r_b);
            Complex<double> g_m_c = std::log(m_c * m_c / mu / mu) / 6.0 - 5.0 / 18.0 - 2.0 * m_c * m_c / 3.0 / s
                + std::sqrt(r_c) / 6.0 * (1.0 + 2.0 * m_c * m_c / s) * Complex<double>::Cartesian(std::log((1.0 + std::sqrt(r_c)) / (1.0 - std::sqrt(r_c))), -M_PI);

            Complex<double> c9eff0 = c9 - (c1 + c2 / 3.0) * (8.0 * g_0 - 4.0 / 3.0) - c3 * (20.0 / 3.0 * g_0 - 16.0 / 3.0 * g_m_b + 2.0 / 27.0)
                + c4 * (4.0 / 3.0 * g_0 + 16.0 / 3.0 * g_m_b + 14.0 / 9.0) - c5 * (8.0 * g_0 - 4.0 * g_m_b - 14.0 / 27.0)
                - c6 * (8.0 / 3.0 * g_0 - 4.0 / 3.0 * g_m_b + 2.0 / 9.0);

            // TODO: Neglecting contributions ~alpha_s / 4.0 / M_PI. These do involve spectator scattering,
            // cf. [BFS2001] Eq. (29), p. 8, and Eqs. (82)-(84), p. 30
            return c9eff0; // cf. [GP2004] Eq. (55), p. 10
        }

        double kappa_1() const
        {
            static const double c0v // cf. [GP] Eq. (48)
                = 1.0 - QCD::alpha_s(mu) * QCD::casimir_f / (4.0 * M_PI) * (3.0 * std::log(mu / m_b_MSbar) + 4.0);
            static const double d0v // cf. [GP] Eq. (A30)
                = QCD::alpha_s(mu) * QCD::casimir_f / (2.0 * M_PI) * (std::log(mu / m_b_MSbar) + 1.0);

            // TODO: [GP] uses m_b(\mu). Which m_b? Using m_b_MSbar for the time
            // being.
            return (1.0 + 2.0 * d0v / c0v) * m_b_MSbar / m_B; // cf. [GP] Eq. (A24)
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2006], p.5
            static const double lambda_t = 0.2272 * 0.2272 * 0.818; // cf. [PDG2006], Eqs. (11.2, 11.25), pp. 183,189

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(const double & s) const
        {
            return s / m_B / m_B;
        }

        // Amplitudes
        Complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            Complex<double> wilson = c9eff(s) + h * c10 + kappa_1() * c7eff(s) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -0.5 * norm(s) * m_B * m_B / m_Kstar / std::sqrt(s));
            double formfactor = lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s)) * FormFactors::a_1(s_hat(s)) - (1 - s_hat(s)) * FormFactors::a_2(s_hat(s));

            return prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            Complex<double> wilson = c9eff(s) + h * c10 + kappa_1() * c7eff(s) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, std::sqrt(2 * lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s))) * norm(s) * m_B);

            return prefactor * wilson * FormFactors::v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = c9eff(s) + h * c10 + kappa_1() * c7eff(s) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -std::sqrt(2) * norm(s) * m_B);

            return prefactor * wilson * FormFactors::a_1(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }
    };

    Decay<BToKstarDilepton>::Decay(const double & mu) :
        PrivateImplementationPattern<Decay<BToKstarDilepton>>(new Implementation<Decay<BToKstarDilepton>>(mu))
    {
    }

    Decay<BToKstarDilepton>::~Decay()
    {
    }

    Complex<double>
    Decay<BToKstarDilepton>::a_long(const Helicity & h, const double & s)
    {
        return _imp->a_long(h, s);
    }

    Complex<double>
    Decay<BToKstarDilepton>::a_perp(const Helicity & h, const double & s)
    {
        return _imp->a_perp(h, s);
    }

    Complex<double>
    Decay<BToKstarDilepton>::a_par(const Helicity & h, const double & s)
    {
        return _imp->a_par(h, s);
    }
}
