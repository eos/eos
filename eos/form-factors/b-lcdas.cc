/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
 * Copyright (c) 2018 Nico Gubernari
 * Copyright (c) 2018 Ahmet Kokulu
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/form-factors/b-lcdas.hh>
#include <eos/utils/options-impl.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/qualified-name.hh>

#include <gsl/gsl_sf_expint.h>

namespace eos
{
    template <>
    struct Implementation<BMesonLCDAs>
    {
        SwitchOption opt_q;

        UsedParameter lambda_B_inv;
        UsedParameter lambda_E2;
        UsedParameter lambda_H2;

        SwitchOption opt_gminus;

        double switch_gminus;

        inline
        QualifiedName parameter(const char * _name) const
        {
            qnp::Name name(_name);

            if (opt_q.value() == "s")
                return QualifiedName(qnp::Prefix("B_s"), name);

            return QualifiedName(qnp::Prefix("B"), name);
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_q(o, "q", { "u", "d", "s" }, "u"),
            lambda_B_inv(p[parameter("1/lambda_B_p").str()], u),
            lambda_E2(p[parameter("lambda_E^2").str()], u),
            lambda_H2(p[parameter("lambda_H^2").str()], u),
            opt_gminus(o, "gminus", { "zero", "WW-limit" }, "WW-limit"),
            switch_gminus(1.0)
        {
            if (opt_gminus.value() == "zero")
            {
                switch_gminus = 0.0;
            }
        }

        /* the inverse moment of phi_+ */
        inline double lambda_B() const
        {
            return 1.0 / lambda_B_inv();
        }

        /* Leading twist two-particle LCDAs */

        inline double phi_plus(const double & omega) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B();

            return omega / (omega_0 * omega_0) * std::exp(-omega / omega_0);
        }

        inline double phi_minus(const double & omega) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B();

            const double limitWW = 1.0 / omega_0 * std::exp(-omega / omega_0);
            const double nonWW   = -(lambda_E2 - lambda_H2) / (18.0 * pow(omega_0, 5)) *
                (
                    2.0 * omega_0 * omega_0 - 4.0 * omega_0 * omega + omega * omega
                ) * std::exp(-omega / omega_0);

            return limitWW + nonWW;
        }

        inline double phi_bar(const double & omega) const
        {
            const double omega_0 = lambda_B();

            const double limitWW = -omega / omega_0 * std::exp(-omega / omega_0);
            const double nonWW   = (lambda_E2 - lambda_H2) / (18.0 * pow(omega_0, 4))
                * (2.0 * omega_0 - omega) * omega * std::exp(-omega / omega_0);

            return limitWW + nonWW;
        }

        inline double phi_bar_d1(const double & omega) const
        {
            return phi_plus(omega) - phi_minus(omega);
        }

        /* Next-to-leading twist two-particle LCDAs */

        inline double g_minusWW(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // Wandzura-Wilcek limit of g_minus

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            return switch_gminus * (3.0 / 4.0 ) * omega * exp;
        }

       inline double g_minusWW_d1(const double & omega) const
        {
            // Wandzura-Wilcek limit of g_minus

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            return switch_gminus * -1.0 * (3.0 / (4.0 * omega_0)) * (omega - omega_0) * exp;
        }

       inline double g_minusWW_d2(const double & omega) const
        {
            // Wandzura-Wilcek limit of g_minus

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            return switch_gminus * (3.0 / (4.0 * omega_0 * omega_0)) * (omega - 2.0 * omega_0) * exp;
        }

        inline double g_plus(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double Ei = gsl_sf_expint_Ei(-omega / omega_0);
            const double exp = std::exp(-omega / omega_0);

            const double termA = -lambda_E2 / (6.0 * pow(omega_0, 2)) *
                (
                    (omega - 2.0 * omega_0) * Ei
                    + (omega + 2.0 * omega_0) * exp * (std::log(omega / omega_0) + gamma_E)
                    - 2.0 * omega * exp
                );
            const double termB = exp / (2.0 * omega_0) * pow(omega, 2) * (
                    1.0 - (lambda_E2 - lambda_H2) / (36.0 * pow(omega_0, 2))
                );

            return termA + termB;
        }

        inline double g_plus_d1(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double Ei = gsl_sf_expint_Ei(-omega / omega_0);
            const double exp = std::exp(-omega / omega_0);

            const double termA = lambda_E2 / (6.0 * pow(omega_0, 3)) *
                (
                    - omega_0 * Ei
                    + (omega + omega_0) * exp * (std::log(omega / omega_0) + gamma_E)
                    - 2.0 * omega * exp
                );
            const double termB = exp / (2.0 * pow(omega_0, 2)) * (2.0 * omega_0 - omega) * omega * (
                    1.0 - (lambda_E2 - lambda_H2) / (36.0 * pow(omega_0, 2))
                );

            return termA + termB;
        }

        inline double g_plus_d2(const double & omega) const
        {
            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            const double termA = lambda_E2 / (6.0 * pow(omega_0, 4)) * exp *
                (
                    - omega_0
                    - omega * (std::log(omega / omega_0) + gamma_E - 2.0)
                );
            const double termB = exp / (2.0 * pow(omega_0, 3)) * (2.0 * pow(omega_0, 2) - 4.0 * omega_0 * omega + pow(omega, 2)) * (
                    1.0 - (lambda_E2 - lambda_H2) / (36.0 * pow(omega_0, 2))
                );

            return termA + termB;
        }

        inline double g_bar(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))

            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double Ei = gsl_sf_expint_Ei(-omega / omega_0);
            const double exp = std::exp(-omega / omega_0);
            const double exp_plus = std::exp(omega / omega_0);

            // integral of g_plus
            const double  termA = -lambda_E2 / (12.0 * pow(omega_0, 2)) *
            ((pow(omega, 2) - 4.0 * omega_0 * omega + 6.0 * pow(omega_0, 2)) * Ei - omega_0 * exp * (std::log(omega / omega_0) + gamma_E) * 2.0 * (3.0 * omega_0 + omega)
             - omega_0 * exp * (omega_0 - 5.0 * omega));
            const double  termB = -exp / 2.0 * (2.0 * pow(omega_0, 2) + 2.0 * omega_0 * omega + pow(omega, 2)) * (1.0 - (lambda_E2 - lambda_H2) / (36.0 * pow(omega_0, 2)));
            const double  int_glus = (termA - lambda_E2 / 12.0) + (termB + pow(omega_0, 2) - (lambda_E2 - lambda_H2) / 36.0);
            // integral of g_minusWW
            const double int_gminusWW = (3.0 / 4.0) * exp * omega_0 * (exp_plus * omega_0 - omega - omega_0);
            return       int_glus - switch_gminus * int_gminusWW;
        }

        inline double g_bar_d1(const double & omega) const
        {
            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))
            return g_plus(omega) - g_minusWW(omega);
        }

        inline double g_bar_d2(const double & omega) const
        {
            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))
            return g_plus_d1(omega) - g_minusWW_d1(omega);
        }

        inline double g_bar_d3(const double & omega) const
        {
            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))
            return g_plus_d2(omega) - g_minusWW_d2(omega);
        }

        /* Leading twist three-particle LCDAs */

        inline double phi_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            // cf. [1703.02446], eq. (5.8), p. 17
            return (lambda_E2 - lambda_H2) / (6.0 * pow(omega_0, 5)) * omega_1 * pow(omega_2, 2) * std::exp(-(omega_1 + omega_2) / omega_0);
        }

        inline double phi_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            // cf. [1703.02446], eq. (5.8), p. 17
            return (lambda_E2 + lambda_H2) / (6.0 * pow(omega_0, 4)) * pow(omega_2, 2) * std::exp(-(omega_1 + omega_2) / omega_0);
        }

        inline double phi_bar_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 - lambda_H2) / (6.0 * pow(omega_0, 4)) * (omega_0 + omega_1) * omega_2 * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 - lambda_H2) / (6.0 * pow(omega_0, 3)) * omega_2 * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        inline double phi_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 + lambda_H2) / (6.0 * pow(omega_0, 3)) * omega_2 * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 + lambda_H2) / (6.0 * pow(omega_0, 3)) * omega_2 * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        inline double phi_bar2_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 - lambda_H2) / (6.0 * pow(omega_0, 4)) * omega_1 * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 - lambda_H2) / (3.0 * pow(omega_0, 2)) * omega_1 * std::exp(- omega_1 / omega_0);

            return termA + termB;
        }

        inline double phi_bar2_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 + lambda_H2) / (6.0 * pow(omega_0, 3)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 + lambda_H2) / (3.0 * omega_0) * std::exp(- omega_1 / omega_0);

            return termA + termB;
        }

        inline double phi_bar_bar_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = (lambda_E2 - lambda_H2) / (6.0 * pow(omega_0, 3)) * (omega_0 + omega_1) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = - (lambda_E2 - lambda_H2) / (3.0 * omega_0) * (omega_0 + omega_1) * std::exp(- omega_1 / omega_0);
            const double termC = - (lambda_E2 - lambda_H2) / (6.0 * pow(omega_0, 2)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(- omega_2 / omega_0);
            const double termD = - 1.0 / 3.0 * (- lambda_E2 + lambda_H2);

            return termA + termB + termC + termD;
        }

        inline double phi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = (lambda_E2 + lambda_H2) / (6.0 * pow(omega_0, 2)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = - 1.0 / 3.0 * (lambda_E2 + lambda_H2) * std::exp(- omega_1 / omega_0);
            const double termC = - (lambda_E2 + lambda_H2) / (6.0 * pow(omega_0, 2)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(- omega_2 / omega_0);
            const double termD = - 1.0 / 3.0 * (- lambda_E2 - lambda_H2);

            return termA + termB + termC + termD;
        }

        inline double psi_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_E2 / (3.0 * pow(omega_0, 3)) * (omega_0 + omega_1) * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_E2 / (3.0 * pow(omega_0, 2)) * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        inline double psi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_E2 / (3.0 * pow(omega_0, 2)) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * (omega_0 + omega_2)  * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_E2 / (3.0 * omega_0) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * std::exp(- omega_1/ omega_0);

            return termA + termB;
        }


        inline double chi_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_H2 / (3.0 * pow(omega_0, 3)) * (omega_0 + omega_1) * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_H2 / (3.0 * pow(omega_0, 2)) * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        inline double chi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_H2 / (3.0 * pow(omega_0, 2)) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * (omega_0 + omega_2)  * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_H2 / (3.0 * omega_0) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * std::exp(- omega_1/ omega_0);

            return termA + termB;
        }

        inline double psi_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return lambda_E_2 / (6.0 * omega_0_4) * xi * xi * std::exp(-(omega + xi) / omega_0);
        }

        inline double psi_V(const double & omega, const double & xi) const
        {
            return psi_A(omega, xi);
        }

        inline double X_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return lambda_E_2 / (6.0 * omega_0_4) * xi * (2.0 * omega - xi) * std::exp(-(omega + xi) / omega_0);
        }

        inline double Y_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return -lambda_E_2 / (24.0 * omega_0_4) * xi * (7.0 * omega_0 - 13.0 * omega + 3.0 * xi) * std::exp(-(omega + xi) / omega_0);
        }

        inline double Xbar_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_3 = omega_0_2 * omega_0;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            // obtained by analytica integrating Y_A(tau, xi) over 0 <= tau <= omega.
            return lambda_E_2 / (6.0 * omega_0_3) * xi * std::exp(-(xi + omega) / omega_0)
                * (xi - 2.0 * (omega + omega_0) + std::exp(omega / omega_0) * (2.0 * omega_0 - xi));
        }

        inline double Ybar_A(const double & omega, const double & xi) const
        {
            // cf. [KMO2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_3 = omega_0_2 * omega_0;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            // obtained by analytica integrating Y_A(tau, xi) over 0 <= tau <= omega.
            return -lambda_E_2 / (24.0 * omega_0_3) * xi * std::exp(-(xi + omega) / omega_0)
                * (-3.0 * xi + 13.0 * omega + 6.0 * omega_0 + 3.0 * std::exp(omega / omega_0) * (xi - 2.0 * omega_0));
        }
    };

    BMesonLCDAs::BMesonLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BMesonLCDAs>(new Implementation<BMesonLCDAs>(p, o, *this))
    {
    }

    BMesonLCDAs::~BMesonLCDAs() = default;

    double
    BMesonLCDAs::phi_plus(const double & omega) const
    {
        return _imp->phi_plus(omega);
    }

    double
    BMesonLCDAs::phi_minus(const double & omega) const
    {
        return _imp->phi_minus(omega);
    }

    double
    BMesonLCDAs::phi_bar(const double & omega) const
    {
        return _imp->phi_bar(omega);
    }

    double
    BMesonLCDAs::phi_bar_d1(const double & omega) const
    {
        return _imp->phi_bar_d1(omega);
    }

    double
    BMesonLCDAs::g_plus(const double & omega) const
    {
        return _imp->g_plus(omega);
    }

    double
    BMesonLCDAs::g_plus_d1(const double & omega) const
    {
        return _imp->g_plus_d1(omega);
    }

    double
    BMesonLCDAs::g_plus_d2(const double & omega) const
    {
        return _imp->g_plus_d2(omega);
    }

    double
    BMesonLCDAs::g_minusWW(const double & omega) const
    {
        return _imp->g_minusWW(omega);
    }

    double
    BMesonLCDAs::g_minusWW_d1(const double & omega) const
    {
        return _imp->g_minusWW_d1(omega);
    }

    double
    BMesonLCDAs::g_minusWW_d2(const double & omega) const
    {
        return _imp->g_minusWW_d2(omega);
    }

    double
    BMesonLCDAs::g_bar(const double & omega) const
    {
        return _imp->g_bar(omega);
    }

    double
    BMesonLCDAs::g_bar_d1(const double & omega) const
    {
        return _imp->g_bar_d1(omega);
    }

    double
    BMesonLCDAs::g_bar_d2(const double & omega) const
    {
        return _imp->g_bar_d2(omega);
    }

    double
    BMesonLCDAs::g_bar_d3(const double & omega) const
    {
        return _imp->g_bar_d3(omega);
    }

    double
    BMesonLCDAs::phi_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_3(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_bar_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_3(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_bar2_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar2_3(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_bar2_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar2_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_bar_bar_3(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_bar_3(omega_1, omega_2);
    }

    double
    BMesonLCDAs::phi_bar_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->phi_bar_bar_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::psi_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->psi_bar_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::psi_bar_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->psi_bar_bar_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::chi_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->chi_bar_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::chi_bar_bar_4(const double & omega_1, const double & omega_2) const
    {
        return _imp->chi_bar_bar_4(omega_1, omega_2);
    }

    double
    BMesonLCDAs::inverse_lambda_plus() const
    {
        return 1.0 / _imp->lambda_B();
    }

    double
    BMesonLCDAs::psi_A(const double & omega, const double & xi) const
    {
        return _imp->psi_A(omega, xi);
    }

    double
    BMesonLCDAs::psi_V(const double & omega, const double & xi) const
    {
        return _imp->psi_V(omega, xi);
    }

    double
    BMesonLCDAs::X_A(const double & omega, const double & xi) const
    {
        return _imp->X_A(omega, xi);
    }

    double
    BMesonLCDAs::Y_A(const double & omega, const double & xi) const
    {
        return _imp->Y_A(omega, xi);
    }

    double
    BMesonLCDAs::Xbar_A(const double & omega, const double & xi) const
    {
        return _imp->Xbar_A(omega, xi);
    }

    double
    BMesonLCDAs::Ybar_A(const double & omega, const double & xi) const
    {
        return _imp->Ybar_A(omega, xi);
    }
}
