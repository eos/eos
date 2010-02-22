/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <tr1/functional>
#include <utility>
#include <map>

#include <iostream>

namespace wf
{
    // TODO: Move to low_recoil/form_factors ?
    // cf. [ABHH1999], p. 8, Table 3
    struct FormFactors
    {
        static double v(const double & s_hat)
        {
            return 0.457 * std::exp(1.482 * s_hat + 1.015 * s_hat * s_hat);
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
    struct Implementation<BToKstarDilepton<LowRecoil>>
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

        double m_b_MSbar;

        double m_c;

        double m_B;

        double m_Kstar;


        Parameter mu;

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
            m_b_MSbar(4.2), // (GeV), cf. [PDG2008], p. 21
            m_c(1.27), // (GeV), cf. [PDG2008], p. 21
            m_B(5.279), // (GeV), cf. [PDG2008], p. 81
            m_Kstar(0.896), // (GeV), cf. [PDG2008], p. 44
            mu(p["mu"])
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
            // All occurences of sqrt(r_c) in Eq. need to be replaced by r_c. Cf. also the footnote
            // on p. 6.
            double r_b = std::sqrt(4.0 * m_b_MSbar * m_b_MSbar / s - 1.0);
            double r_c = std::sqrt(1.0 - 4.0 * s / m_c / m_c);
            Complex<double> g_0 = Complex<double>::Cartesian(std::log(s / mu / mu), -M_PI) * (1.0 / 6.0) - 5.0 / 18.0;
            double g_m_b = std::log(m_b_MSbar * m_b_MSbar / mu / mu) / 6.0 - 5.0 / 18.0 - 2.0 * m_b_MSbar * m_b_MSbar / 3.0 / s
                + r_b / 3.0 * (1.0 + 2.0 * m_b_MSbar * m_b_MSbar / s) * std::atan(1.0 / r_b);
            Complex<double> g_m_c = std::log(m_c * m_c / mu / mu) / 6.0 - 5.0 / 18.0 - 2.0 * m_c * m_c / 3.0 / s
                + r_c / 6.0 * (1.0 + 2.0 * m_c * m_c / s) * Complex<double>::Cartesian(std::log((1.0 + r_c) / (1.0 - r_c)), -M_PI);

            Complex<double> c9eff0 = c9() - (c1 + c2 / 3.0) * (8.0 * g_0 - 4.0 / 3.0) - c3() * (20.0 / 3.0 * g_0 - 16.0 / 3.0 * g_m_b + 2.0 / 27.0)
                + c4() * (4.0 / 3.0 * g_0 + 16.0 / 3.0 * g_m_b + 14.0 / 9.0) - c5() * (8.0 * g_0 - 4.0 * g_m_b - 14.0 / 27.0)
                - c6() * (8.0 / 3.0 * g_0 - 4.0 / 3.0 * g_m_b + 2.0 / 9.0);

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
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5
            static const double lambda_t = -0.2257 * 0.2257 * 0.814; // cf. [PDG2008], Eqs. (11.5, 11.26), pp. 169,174

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
            Complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) + c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -0.5 * norm(s) * m_B * m_B * m_B / m_Kstar / std::sqrt(s));
            double formfactor = lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s)) * FormFactors::a_1(s_hat(s)) - (1 - s_hat(s)) * FormFactors::a_2(s_hat(s));

            return prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            Complex<double> wilson = (c9eff(s) + c9prime()) + h * (c10() + c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, std::sqrt(2 * lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s))) * norm(s) * m_B);

            return prefactor * wilson * FormFactors::v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -std::sqrt(2) * norm(s) * m_B);

            return prefactor * wilson * FormFactors::a_1(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }
    };

    BToKstarDilepton<LowRecoil>::BToKstarDilepton(const Parameters & parameters) :
        PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>(new Implementation<BToKstarDilepton<LowRecoil>>(parameters))
    {
    }

    BToKstarDilepton<LowRecoil>::~BToKstarDilepton()
    {
    }

    Complex<double>
    BToKstarDilepton<LowRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    Complex<double>
    BToKstarDilepton<LowRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    Complex<double>
    BToKstarDilepton<LowRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) / Gamma;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_decay_width(const double & s) const
    {
        return a_long(left_handed, s).absolute_squared()
            + a_long(right_handed, s).absolute_squared()
            + a_perp(left_handed, s).absolute_squared()
            + a_perp(right_handed, s).absolute_squared()
            + a_par(left_handed, s).absolute_squared()
            + a_par(right_handed, s).absolute_squared();
    }

    double
    BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 / differential_decay_width(s) * (
                (a_par(left_handed, s) * a_perp(left_handed, s).conjugate()).real()
                -(a_par(right_handed, s) * a_perp(right_handed, s).conjugate()).real()
            );
    }

    double
    BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
            / differential_decay_width(s);
    }

    double integrate(const std::tr1::function<double (const double &)> & f, unsigned n, const double & a, const double & b)
    {
        if (n & 0x1)
            n += 1;

        double h = (b - a) / n;

        double result = f(a) + f(b);
        for (unsigned i = 1 ; i <= n / 2 - 1 ; ++i)
        {
            result += 2.0 * f(a + h * (2 * i));
        }
        for (unsigned i = 1 ; i <= n / 2 ; ++i)
        {
            result += 4.0 * f(a + h * (2 * i - 1));
        }

        result *= h / 3;

        return result;
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrate(std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_branching_ratio),
                    this, std::tr1::placeholders::_1), 100, s_min, s_max);
    }

    double
    BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        return integrate(std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry),
                    this, std::tr1::placeholders::_1), 100, s_min, s_max);
    }

    std::tr1::shared_ptr<Observable>
    BToKstarDileptonFactory::make(const std::string & name, const Parameters & parameters)
    {
        static const std::map<std::string, ObservableFactory *> simple_observables
        {
            std::make_pair("dBR/ds", new ConcreteObservableFactory<BToKstarDilepton<LowRecoil>, 1>(
                        ConcreteObservableData<BToKstarDilepton<LowRecoil>, 1>("dBR/ds",
                            std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_branching_ratio),
                            "s"))),
            std::make_pair("A_FB(s)", new ConcreteObservableFactory<BToKstarDilepton<LowRecoil>, 1>(
                        ConcreteObservableData<BToKstarDilepton<LowRecoil>, 1>("dA_FB/ds",
                            std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry),
                            "s"))),
            std::make_pair("F_L(s)", new ConcreteObservableFactory<BToKstarDilepton<LowRecoil>, 1>(
                        ConcreteObservableData<BToKstarDilepton<LowRecoil>, 1>("dF_L/ds",
                            std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation),
                            "s"))),
            std::make_pair("A_FB", new ConcreteObservableFactory<BToKstarDilepton<LowRecoil>, 2>(
                        ConcreteObservableData<BToKstarDilepton<LowRecoil>, 2>("A_FB",
                            std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry),
                            "s_min", "s_max"))),
            std::make_pair("BR", new ConcreteObservableFactory<BToKstarDilepton<LowRecoil>, 2>(
                        ConcreteObservableData<BToKstarDilepton<LowRecoil>, 2>("BR",
                            std::tr1::mem_fn(&BToKstarDilepton<LowRecoil>::integrated_branching_ratio),
                            "s_min", "s_max")))
        };

        auto i(simple_observables.find(name));

        if (simple_observables.end() == i)
            return std::tr1::shared_ptr<Observable>();

        return i->second->make(parameters);
    }
}
