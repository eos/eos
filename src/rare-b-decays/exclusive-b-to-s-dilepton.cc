/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <src/rare-b-decays/form_factors.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

#include <cmath>
#include <tr1/functional>
#include <utility>
#include <map>

namespace wf
{
    // Large Recoil

    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
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

        double m_b_pole;

        double m_c;

        double m_B;

        double m_Kstar;

        double m_l;

        Parameter mu;

        double lambda_t2;

        Complex<double> lambda_u_hat;

        std::tr1::shared_ptr<FormFactors<BToKstar>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
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
            m_b_pole(4.8), // (GeV)
            m_c(1.27), // (GeV), cf. [PDG2008], p. 21
            m_B(5.279), // (GeV), cf. [PDG2008], p. 81
            m_Kstar(0.896), // (GeV), cf. [PDG2008], p. 44
            mu(p["mu"])
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"]);
            if (! form_factors.get())
                throw std::string("InternalError");

            // TODO: Lepton masses, m_l = m_mu
            m_l = 0.10565836; // (GeV), cf. [PDG2008], p. 13

            double A = p["CKM::A"];
            double A2 = A * A;
            double A4 = A2 * A2;
            double lambda = p["CKM::lambda"];
            double lambda2 = lambda * lambda;
            double lambda4 = lambda2 * lambda2;
            double etabar = p["CKM::etabar"];
            double rhobar = p["CKM::rhobar"];

            // cf. [PDG2008], Eq. (11.4)
            lambda_u_hat = -1.0 * lambda2 * std::sqrt(1.0 - A2 * lambda4) / std::sqrt(1.0 - lambda2)
                / (1.0 - 2 * A2 * lambda4 * rhobar + A4 * lambda4 * lambda4 * (rhobar * rhobar + etabar * etabar))
                * Complex<double>::Cartesian(rhobar - A2 * lambda4 * (rhobar * rhobar + etabar * etabar), -etabar);
            lambda_t2 = A2 * lambda4;
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double norm(const double & s) const
        {
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]
            static const double g_fermi = 1.16637e-5; // (Gev^-2 (hbar c)^3), cf. [PDG2008], p.5
            static const double lambda_t = -0.2257 * 0.2257 * 0.814; // cf. [PDG2008], Eqs. (11.5, 11.26), pp. 169,174

            return std::sqrt(g_fermi * g_fermi * alpha_e * alpha_e / 3.0 / 1024 / std::pow(M_PI, 5.0) / m_B
                    * lambda_t2 * s_hat(s)
                    * std::sqrt(lambda(m_B * m_B, m_Kstar * m_Kstar, s))); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(double s) const
        {
            return s / m_B / m_B;
        }

        /* Effective wilson coefficients */
        double c7eff() const
        {
            return c7() - 1.0/3.0 * c3() - 4.0/9.0 * c4() - 20.0/3.0 * c5() - 80.0/9.0 * c6();
        }

        /* Charm loop parts */
        Complex<double> h(const double & s, const double & m) const
        {
            const double z = 4.0 * m * m / s;
            const double sqrt1z = std::sqrt(std::abs(z - 1.0));

            double a = 2.0 * std::log(m / mu()) - 2.0 / 3.0 - z;
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

            return -4.0 / 9.0 * (a + b * Complex<double>::Cartesian(rc, ic));
        }

        Complex<double> h0(const double & s) const
        {
            return 4.0 / 9.0 * Complex<double>::Cartesian(2.0 / 3.0 + 2.0 * std::log(2.0 * mu()) - std::log(s), M_PI);
        }

        Complex<double> Y0(const double & s) const
        {
            double Y_c = 4.0 / 3.0 * c1() + c2() + 6.0 * c3() + 60.0 * c5();
            double Y_b = -0.5 * (7.0 * c3() + 4.0 / 3.0 * (c4() + 76.0 * c5()) + 64.0 / 3.0 * c6());
            double Y_0 = -0.5 * (c3() + 4.0 / 3.0 * c4() + 16.0 * c5() - 64 / 3.0 * c6());
            double Y = 2.0 / 9.0 * (6.0 * c3() + 32.0 * c5() + 32.0 / 3.0 * c6());

            // Uses b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses
            return Y_c * h(s, m_c) + Y_b * h(s, m_b_pole) + Y_0 * h0(s) + Y;
        }

        Complex<double> Y0u(const double & s) const
        {
            double a = 4.0 / 3.0 * c1() + c2();

            return a * (h(s, m_c) - h0(s));
        }

        // Form factors
        double xi_perp(const double & s) const
        {
            const double factor = m_B / (m_B + m_Kstar);

            return factor * form_factors->v(s_hat(s));
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B * (m_B + m_Kstar)) / (m_B * m_B + m_Kstar * m_Kstar - s);
            const double factor2 = (1.0 - m_Kstar / m_B);

            return factor1 * form_factors->a_1(s_hat(s)) - factor2 * form_factors->a_2(s_hat(s));
        }

        Complex<double> tensor_perp(const double & h, const double & s) const
        {
            return c7eff() + h * c7prime() + s / (2.0 * QCD::mb_PS(m_b_pole, mu()) * m_B) * (Y0(s) + lambda_u_hat * Y0u(s));
        }

        Complex<double> tensor_par(const double & s) const
        {
            return -1.0 * (c7eff() - c7prime() + m_B / (2.0 * QCD::mb_PS(m_b_pole, mu())) * (Y0(s) + lambda_u_hat * Y0u(s)));
        }

        // Amplitudes
        Complex<double> a_long(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = QCD::mb_PS(m_b_pole, mu()) / m_B;

            double prefactor = -1.0 * m_B * m_B * (1.0 - shat) * (1.0 - shat) / (2.0 * m_Kstar * std::sqrt(shat));
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return prefactor * (wilson - 2.0 * mbhat * tensor_par(s)) * xi_par(s);
        }

        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = QCD::mb_PS(m_b_pole, mu()) / m_B;

            double prefactor = +std::sqrt(2.0) * m_B * (1.0 - shat);
            double wilson = (c9() + c9prime()) + h * (c10() + c10prime());

            return prefactor * (wilson + 2.0 * mbhat / shat * tensor_perp(+1.0, s)) * xi_perp(s);
        }

        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double shat = s_hat(s);
            double mbhat = QCD::mb_PS(m_b_pole, mu()) / m_B;

            double prefactor = -std::sqrt(2.0) * m_B * (1.0 - shat);
            double wilson = (c9() - c9prime()) + h * (c10() - c10prime());

            return prefactor * (wilson + 2.0 * mbhat / shat * tensor_perp(-1.0, s)) * xi_perp(s);
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    Complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        return _imp->a_long(h, s);
    }

    Complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        return _imp->a_perp(h, s);
    }

    Complex<double>
    BToKstarDilepton<LargeRecoil>::a_par(const Helicity & h, const double & s) const
    {
        return _imp->a_par(h, s);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        return differential_decay_width(s) / Gamma;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_decay_width(const double & s) const
    {
        return _imp->norm(s) * _imp->norm(s) * (a_long(left_handed, s).absolute_squared()
            + a_long(right_handed, s).absolute_squared()
            + a_perp(left_handed, s).absolute_squared()
            + a_perp(right_handed, s).absolute_squared()
            + a_par(left_handed, s).absolute_squared()
            + a_par(right_handed, s).absolute_squared());
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return 1.5 * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s) * (
                (a_par(left_handed, s) * a_perp(left_handed, s).conjugate()).real()
                -(a_par(right_handed, s) * a_perp(right_handed, s).conjugate()).real()
            );
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared();
        double b = a_par(left_handed, s).absolute_squared() + a_par(right_handed, s).absolute_squared();

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_par(left_handed, s).conjugate() + a_long(right_handed, s).conjugate() * a_par(right_handed, s)).absolute();
        double b = std::sqrt(a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
                * std::sqrt(a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared());

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_perp(left_handed, s).conjugate() - a_long(right_handed, s).conjugate() * a_perp(right_handed, s)).absolute();
        double b = (a_long(left_handed, s).conjugate() * a_par(left_handed, s) + a_long(right_handed, s) * a_par(right_handed, s).conjugate()).absolute();

        return a / b;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
            * _imp->norm(s) * _imp->norm(s) / differential_decay_width(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrate(std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_branching_ratio),
                    this, std::tr1::placeholders::_1), 100, s_min, s_max);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        return integrate(std::tr1::bind(std::tr1::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry),
                    this, std::tr1::placeholders::_1), 100, s_min, s_max);
    }


    // Low Recoil

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

        Parameter m_b_MSbar;

        Parameter m_c;

        Parameter m_B;

        Parameter m_Kstar;

        Parameter mu;

        std::tr1::shared_ptr<FormFactors<BToKstar>> form_factors;

        Implementation(const Parameters & p, const ObservableOptions & o) :
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
            m_B(p["mass::B0"]),
            m_Kstar(p["mass::K^*0"]),
            mu(p["mu"])
        {
            form_factors = FormFactorFactory<BToKstar>::create(o["form-factors"]);
            if (! form_factors.get())
                throw std::string("InternalError");
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
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -0.5 * m_B * m_B * m_B / m_Kstar / std::sqrt(s));
            double formfactor = lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s)) * form_factors->a_1(s_hat(s)) - (1 - s_hat(s)) * form_factors->a_2(s_hat(s));

            return prefactor * wilson * formfactor; // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_perp(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            double m_Kstar_hat = m_Kstar / m_B;
            Complex<double> wilson = (c9eff(s) + c9prime()) + h * (c10() + c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, std::sqrt(2 * lambda(1.0, m_Kstar_hat * m_Kstar_hat, s_hat(s))) * m_B);

            return prefactor * wilson * form_factors->v(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }

        Complex<double> a_par(const Helicity & helicity, const double & s) const
        {
            double h = helicity;
            Complex<double> wilson = (c9eff(s) - c9prime()) + h * (c10() - c10prime()) + kappa_1() * (c7eff(s) - c7prime()) * (2 * m_b_MSbar * m_B / s);
            Complex<double> prefactor = Complex<double>::Cartesian(0.0, -std::sqrt(2) * m_B);

            return prefactor * wilson * form_factors->a_1(s_hat(s)); // cf. [BHvD2010], Eq. (??)
        }
    };

    BToKstarDilepton<LowRecoil>::BToKstarDilepton(const Parameters & parameters, const ObservableOptions & options) :
        PrivateImplementationPattern<BToKstarDilepton<LowRecoil>>(new Implementation<BToKstarDilepton<LowRecoil>>(parameters, options))
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

        return differential_decay_width(s) * (_imp->norm(s) * _imp->norm(s) / Gamma);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_decay_width(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared()
            + a_long(right_handed, s).absolute_squared()
            + a_perp(left_handed, s).absolute_squared()
            + a_perp(right_handed, s).absolute_squared()
            + a_par(left_handed, s).absolute_squared()
            + a_par(right_handed, s).absolute_squared());
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
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        double a = a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared();
        double b = a_par(left_handed, s).absolute_squared() + a_par(right_handed, s).absolute_squared();

        return (a - b) / (a + b);
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_par(left_handed, s).conjugate() + a_long(right_handed, s).conjugate() * a_par(right_handed, s)).absolute();
        double b = std::sqrt((a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
                * (a_perp(left_handed, s).absolute_squared() + a_perp(right_handed, s).absolute_squared()));

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        double a = (a_long(left_handed, s) * a_perp(left_handed, s).conjugate() - a_long(right_handed, s).conjugate() * a_perp(right_handed, s)).absolute();
        double b = (a_long(left_handed, s).conjugate() * a_par(left_handed, s) + a_long(right_handed, s) * a_par(right_handed, s).conjugate()).absolute();

        return a / b;
    }

    double
    BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        return (a_long(left_handed, s).absolute_squared() + a_long(right_handed, s).absolute_squared())
            / differential_decay_width(s);
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
}
