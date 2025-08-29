/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2025 Danny van Dyk
 * Copyright (c) 2018      Nico Gubernari
 * Copyright (c) 2018      Ahmet Kokulu
 * Copyright (c) 2022      Philip Lüghausen
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

#include <eos/form-factors/heavy-meson-lcdas-exponential.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <gsl/gsl_sf_expint.h>

namespace eos
{
    namespace heavy_meson_lcdas
    {
        using namespace std::literals::string_literals;

        const std::vector<OptionSpecification> Exponential::options
        {
            { "Q"_ok,      { "b"s },                 "b"s        },
            { "q"_ok,      { "u"s, "d"s, "s"s },     "u"s        },
            { "gminus"_ok, { "zero"s, "WW-limit"s }, "WW-limit"s }
        };

        Exponential::Exponential(const Parameters & p, const Options & o) :
            opt_Q(o, options, "Q"_ok),
            opt_q(o, options, "q"_ok),
            opt_gminus(o, options, "gminus"_ok),
            lambda_B_inv(p[parameter("1/lambda_B_p")], *this),
            lambda_E2(p[parameter("lambda_E^2")], *this),
            lambda_H2(p[parameter("lambda_H^2")], *this),
            switch_gminus(1.0)
        {
            if (opt_gminus.value() == "zero")
            {
                switch_gminus = 0.0;
            }
        }

        std::string
        Exponential::parameter(const char * _name) const
        {
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, qnp::Prefix> prefixes
            {
                { { QuarkFlavor::bottom, QuarkFlavor::up },      qnp::Prefix("B")   },
                { { QuarkFlavor::bottom, QuarkFlavor::down },    qnp::Prefix("B")   },
                { { QuarkFlavor::bottom, QuarkFlavor::strange }, qnp::Prefix("B_s") }
            };

            auto it = prefixes.find(std::make_tuple(opt_Q.value(), opt_q.value()));
            if (it == prefixes.end())
                throw InternalError("Combination of options Q=" + opt_Q.str() + ", q=" + opt_q.str() + " is not supported");

            return QualifiedName(it->second, qnp::Name(_name)).str();
        }

        HeavyMesonLCDAs *
        Exponential::make(const Parameters & p, const Options & o)
        {
            return new Exponential(p, o);
        }

        /* the inverse moment of phi_+ */
        double
        Exponential::lambda_B() const
        {
            return 1.0 / lambda_B_inv();
        }

        double
        Exponential::inverse_lambda_plus() const
        {
            return lambda_B_inv();
        }

        /* Leading twist two-particle LCDAs */

        double
        Exponential::phi_plus(const double & omega) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B();

            return omega / (omega_0 * omega_0) * std::exp(-omega / omega_0);
        }

        double
        Exponential::phi_minus(const double & omega) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B();

            const double limitWW = 1.0 / omega_0 * std::exp(-omega / omega_0);
            const double nonWW   = -(lambda_E2 - lambda_H2) / (18.0 * power_of<5>(omega_0)) *
                (
                    2.0 * omega_0 * omega_0 - 4.0 * omega_0 * omega + omega * omega
                ) * std::exp(-omega / omega_0);

            return limitWW + nonWW;
        }

        double
        Exponential::phi_bar(const double & omega) const
        {
            const double omega_0 = lambda_B();
            const double limitWW = -omega / omega_0 * std::exp(-omega / omega_0);
            const double nonWW   = (lambda_E2 - lambda_H2) / (18.0 * power_of<4>(omega_0))
                * (2.0 * omega_0 - omega) * omega * std::exp(-omega / omega_0);

            return limitWW + nonWW;
        }

        double
        Exponential::phi_bar_d1(const double & omega) const
        {
            return phi_plus(omega) - phi_minus(omega);
        }

        /* Next-to-leading twist two-particle LCDAs */

        double
        Exponential::g_minusWW(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // Wandzura-Wilcek limit of g_minus

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            return switch_gminus * (3.0 / 4.0 ) * omega * exp;
        }

        double
        Exponential::g_minusWW_d1(const double & omega) const
        {
            // Wandzura-Wilcek limit of g_minus

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            return switch_gminus * -1.0 * (3.0 / (4.0 * omega_0)) * (omega - omega_0) * exp;
        }

        double
        Exponential::g_minusWW_d2(const double & omega) const
        {
            // Wandzura-Wilcek limit of g_minus

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            return switch_gminus * (3.0 / (4.0 * omega_0 * omega_0)) * (omega - 2.0 * omega_0) * exp;
        }

        double
        Exponential::g_plus(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double Ei = gsl_sf_expint_Ei(-omega / omega_0);
            const double exp = std::exp(-omega / omega_0);

            const double termA = -lambda_E2 / (6.0 * power_of<2>(omega_0)) *
                (
                    (omega - 2.0 * omega_0) * Ei
                    + (omega + 2.0 * omega_0) * exp * (std::log(omega / omega_0) + gamma_E)
                    - 2.0 * omega * exp
                );
            const double termB = exp / (2.0 * omega_0) * power_of<2>(omega) * (
                    1.0 - (lambda_E2 - lambda_H2) / (36.0 * power_of<2>(omega_0))
                );

            return termA + termB;
        }

        double
        Exponential::g_plus_d1(const double & omega) const
        {
            if (omega < 1.0e-5)
                return 0.0;

            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double Ei = gsl_sf_expint_Ei(-omega / omega_0);
            const double exp = std::exp(-omega / omega_0);

            const double termA = lambda_E2 / (6.0 * power_of<3>(omega_0)) *
                (
                    - omega_0 * Ei
                    + (omega + omega_0) * exp * (std::log(omega / omega_0) + gamma_E)
                    - 2.0 * omega * exp
                );
            const double termB = exp / (2.0 * power_of<2>(omega_0)) * (2.0 * omega_0 - omega) * omega * (
                    1.0 - (lambda_E2 - lambda_H2) / (36.0 * power_of<2>(omega_0))
                );

            return termA + termB;
        }

        double
        Exponential::g_plus_d2(const double & omega) const
        {
            // Euler-Mascheroni gamma constant
            constexpr static double gamma_E = 0.57721566490153286;

            const double omega_0 = lambda_B();
            const double exp = std::exp(-omega / omega_0);

            const double termA = lambda_E2 / (6.0 * power_of<4>(omega_0)) * exp *
                (
                    - omega_0
                    - omega * (std::log(omega / omega_0) + gamma_E - 2.0)
                );
            const double termB = exp / (2.0 * power_of<3>(omega_0)) * (2.0 * power_of<2>(omega_0) - 4.0 * omega_0 * omega + power_of<2>(omega)) * (
                    1.0 - (lambda_E2 - lambda_H2) / (36.0 * power_of<2>(omega_0))
                );

            return termA + termB;
        }

        double
        Exponential::g_bar(const double & omega) const
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
            const double  termA = -lambda_E2 / (12.0 * power_of<2>(omega_0)) *
            ((power_of<2>(omega) - 4.0 * omega_0 * omega + 6.0 * power_of<2>(omega_0)) * Ei - omega_0 * exp * (std::log(omega / omega_0) + gamma_E) * 2.0 * (3.0 * omega_0 + omega)
             - omega_0 * exp * (omega_0 - 5.0 * omega));
            const double  termB = -exp / 2.0 * (2.0 * power_of<2>(omega_0) + 2.0 * omega_0 * omega + power_of<2>(omega)) * (1.0 - (lambda_E2 - lambda_H2) / (36.0 * power_of<2>(omega_0)));
            const double  int_glus = (termA - lambda_E2 / 12.0) + (termB + power_of<2>(omega_0) - (lambda_E2 - lambda_H2) / 36.0);
            // integral of g_minusWW
            const double int_gminusWW = (3.0 / 4.0) * exp * omega_0 * (exp_plus * omega_0 - omega - omega_0);
            return       int_glus - switch_gminus * int_gminusWW;
        }

        double
        Exponential::g_bar_d1(const double & omega) const
        {
            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))
            return g_plus(omega) - g_minusWW(omega);
        }

        double
        Exponential::g_bar_d2(const double & omega) const
        {
            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))
            return g_plus_d1(omega) - g_minusWW_d1(omega);
        }

        double
        Exponential::g_bar_d3(const double & omega) const
        {
            // including the WW-limit of g_minus
            // in this case: g_bar = \int_0^omega d(eta) (g_plus(eta) - g_minusWW(eta))
            return g_plus_d2(omega) - g_minusWW_d2(omega);
        }

        /* Leading twist three-particle LCDAs */

        double
        Exponential::phi_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            // cf. [1703.02446], eq. (5.8), p. 17
            return (lambda_E2 - lambda_H2) / (6.0 * power_of<5>(omega_0)) * omega_1 * power_of<2>(omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
        }

        double
        Exponential::phi_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            // cf. [1703.02446], eq. (5.8), p. 17
            return (lambda_E2 + lambda_H2) / (6.0 * power_of<4>(omega_0)) * power_of<2>(omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
        }

        double
        Exponential::phi_bar_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 - lambda_H2) / (6.0 * power_of<4>(omega_0)) * (omega_0 + omega_1) * omega_2 * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 - lambda_H2) / (6.0 * power_of<3>(omega_0)) * omega_2 * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        double
        Exponential::phi_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 + lambda_H2) / (6.0 * power_of<3>(omega_0)) * omega_2 * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 + lambda_H2) / (6.0 * power_of<3>(omega_0)) * omega_2 * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        double
        Exponential::phi_bar2_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 - lambda_H2) / (6.0 * power_of<4>(omega_0)) * omega_1 * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 - lambda_H2) / (3.0 * power_of<2>(omega_0)) * omega_1 * std::exp(- omega_1 / omega_0);

            return termA + termB;
        }

        double
        Exponential::phi_bar2_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - (lambda_E2 + lambda_H2) / (6.0 * power_of<3>(omega_0)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = (lambda_E2 + lambda_H2) / (3.0 * omega_0) * std::exp(- omega_1 / omega_0);

            return termA + termB;
        }

        double
        Exponential::phi_bar_bar_3(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = (lambda_E2 - lambda_H2) / (6.0 * power_of<3>(omega_0)) * (omega_0 + omega_1) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = - (lambda_E2 - lambda_H2) / (3.0 * omega_0) * (omega_0 + omega_1) * std::exp(- omega_1 / omega_0);
            const double termC = - (lambda_E2 - lambda_H2) / (6.0 * power_of<2>(omega_0)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(- omega_2 / omega_0);
            const double termD = - 1.0 / 3.0 * (- lambda_E2 + lambda_H2);

            return termA + termB + termC + termD;
        }

        double
        Exponential::phi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = (lambda_E2 + lambda_H2) / (6.0 * power_of<2>(omega_0)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = - 1.0 / 3.0 * (lambda_E2 + lambda_H2) * std::exp(- omega_1 / omega_0);
            const double termC = - (lambda_E2 + lambda_H2) / (6.0 * power_of<2>(omega_0)) * (2.0 * omega_0 * omega_0 + 2.0 * omega_0 * omega_2 + omega_2 * omega_2) * std::exp(- omega_2 / omega_0);
            const double termD = - 1.0 / 3.0 * (- lambda_E2 - lambda_H2);

            return termA + termB + termC + termD;
        }

        double
        Exponential::psi_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_E2 / (3.0 * power_of<3>(omega_0)) * (omega_0 + omega_1) * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_E2 / (3.0 * power_of<2>(omega_0)) * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        double
        Exponential::psi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_E2 / (3.0 * power_of<2>(omega_0)) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * (omega_0 + omega_2)  * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_E2 / (3.0 * omega_0) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * std::exp(- omega_1/ omega_0);

            return termA + termB;
        }


        double
        Exponential::chi_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_H2 / (3.0 * power_of<3>(omega_0)) * (omega_0 + omega_1) * omega_2 * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_H2 / (3.0 * power_of<2>(omega_0)) * omega_2 * std::exp(- omega_2 / omega_0);

            return termA + termB;
        }

        double
        Exponential::chi_bar_bar_4(const double & omega_1, const double & omega_2) const
        {
            const double omega_0 = lambda_B();

            const double termA = - lambda_H2 / (3.0 * power_of<2>(omega_0)) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * (omega_0 + omega_2)  * std::exp(-(omega_1 + omega_2) / omega_0);
            const double termB = lambda_H2 / (3.0 * omega_0) * ((-1.0 +  std::exp( omega_1 / omega_0))
                               * omega_0 - omega_1) * std::exp(- omega_1/ omega_0);

            return termA + termB;
        }

        double
        Exponential::psi_A(const double & omega, const double & xi) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return lambda_E_2 / (6.0 * omega_0_4) * xi * xi * std::exp(-(omega + xi) / omega_0);
        }

        double
        Exponential::psi_V(const double & omega, const double & xi) const
        {
            return psi_A(omega, xi);
        }

        double
        Exponential::X_A(const double & omega, const double & xi) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return lambda_E_2 / (6.0 * omega_0_4) * xi * (2.0 * omega - xi) * std::exp(-(omega + xi) / omega_0);
        }

        double
        Exponential::Y_A(const double & omega, const double & xi) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_4 = omega_0_2 * omega_0_2;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            return -lambda_E_2 / (24.0 * omega_0_4) * xi * (7.0 * omega_0 - 13.0 * omega + 3.0 * xi) * std::exp(-(omega + xi) / omega_0);
        }

        double
        Exponential::Xbar_A(const double & omega, const double & xi) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_3 = omega_0_2 * omega_0;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            // obtained by analytica integrating Y_A(tau, xi) over 0 <= tau <= omega.
            return lambda_E_2 / (6.0 * omega_0_3) * xi * std::exp(-(xi + omega) / omega_0)
                * (xi - 2.0 * (omega + omega_0) + std::exp(omega / omega_0) * (2.0 * omega_0 - xi));
        }

        double
        Exponential::Ybar_A(const double & omega, const double & xi) const
        {
            // cf. [KMO:2006], eq. (53), p. 16
            const double omega_0 = lambda_B(), omega_0_2 = omega_0 * omega_0, omega_0_3 = omega_0_2 * omega_0;
            const double lambda_E_2 = 3.0 / 2.0 * omega_0_2;

            // obtained by analytica integrating Y_A(tau, xi) over 0 <= tau <= omega.
            return -lambda_E_2 / (24.0 * omega_0_3) * xi * std::exp(-(xi + omega) / omega_0)
                * (-3.0 * xi + 13.0 * omega + 6.0 * omega_0 + 3.0 * std::exp(omega / omega_0) * (xi - 2.0 * omega_0));
        }

        std::tuple<HeavyMesonLCDAs::CoefficientIterator, HeavyMesonLCDAs::CoefficientIterator>
        Exponential::coefficient_range(const double & /* mu */) const
        {
            static const std::array<double, 9> cs = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            return {cs.begin(), cs.end()};
        }

        Diagnostics Exponential::diagnostics() const
        {
            Diagnostics results;
            // add diagnostic results
            return results;
        }
    }

    template <>
    struct WrappedForwardIteratorTraits<HeavyMesonLCDAs::CoefficientIteratorTag>
    {
        using UnderlyingIterator = std::array<double, 9>::const_iterator;
    };
    template class WrappedForwardIterator<HeavyMesonLCDAs::CoefficientIteratorTag, const double &>;
}
