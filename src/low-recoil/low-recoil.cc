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
        Complex<double> c7eff;

        Complex<double> c9eff;

        Complex<double> c10;

        double m_b_MSbar;

        double m_B;

        double m_Kstar;


        double mu;

        Implementation(double mu) :
            c7eff(Complex<double>::Cartesian(-0.313, 0)), // TODO: Compute correct value from [GP]
            c9eff(Complex<double>::Cartesian(+4.344, 0)), // TODO: Compute correct value from [GP]
            c10(Complex<double>::Cartesian(-4.669, 0)), // TODO: Compute correct value from [GP]
            m_b_MSbar(4.2), // (GeV), cf. [PDG2006], p. 24
            m_B(5.279), // (GeV), cf. [PDG2006], p. 87
            m_Kstar(0.896), // (GeV), cf. [PDG2006], p. 51
            mu(mu)
        {
        }

        // Helpers
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
                    * lambda(m_B * m_B, m_Kstar * m_Kstar, s)); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(const double & s) const
        {
            return s / m_B / m_B;
        }

        // Amplitudes
        Complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = c9eff + h * c10 + kappa_1() * c7eff * (2 * m_b_MSbar * m_b_MSbar / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -0.5 * norm(s) * m_B * m_B / m_Kstar / std::sqrt(s));
            double formfactor = lambda(1, 0, s_hat(s)) * FormFactors::a_1(s_hat(s)) - (1 - s_hat(s)) * FormFactors::a_2(s_hat(s));

            return prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = c9eff + h * c10 + kappa_1() * c7eff * c7eff * (2 * m_b_MSbar * m_b_MSbar / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, std::sqrt(2 * lambda(1, 0, s_hat(s))) * norm(s) * m_B);

            return prefactor * wilson * FormFactors::v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = c9eff + h * c10 + kappa_1() * c7eff * c7eff * (2 * m_b_MSbar * m_b_MSbar / s);
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
        return _imp->a_long(h, s);
    }

    Complex<double>
    Decay<BToKstarDilepton>::a_par(const Helicity & h, const double & s)
    {
        return _imp->a_long(h, s);
    }
}
