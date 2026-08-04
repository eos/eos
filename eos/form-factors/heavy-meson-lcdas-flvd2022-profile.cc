/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/form-factors/heavy-meson-lcdas-flvd2022-profile.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <map>
#include <numeric>

namespace eos
{
    namespace heavy_meson_lcdas
    {
        using namespace std::literals::string_literals;

        FLvD2022Profile::FLvD2022Profile(const Parameters & p, const Options & o) :
            leading_twist(p, o),
            opt_Q(o, options, "Q"_ok),
            opt_q(o, options, "q"_ok),
            mu_0(p[flvd2022_parameter("mu_0")], *this),
            omega_0(p[flvd2022_parameter("omega_0")], *this),
            lambda_bar(p[parameter("LambdaBar")], *this),
            lambda_E2(p[parameter("lambda_E^2")], *this),
            lambda_H2(p[parameter("lambda_H^2")], *this)
        {
            // Verify the size of Weights used internally
            Weights weights;
            auto [a_begin, a_end] = leading_twist.coefficient_range(mu_0);
            if (weights.size() < static_cast<std::size_t>(std::distance(a_begin, a_end)))
            {
                throw InternalError("The number of weights implemented is smaller than the number of coefficients of phi_+");
            }

            this->uses(leading_twist);
        }

        std::string
        FLvD2022Profile::parameter(const char * _name) const
        {
            // The hadronic higher-twist parameters use the same names as in the exponential
            // model, cf. Exponential::parameter().
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, qnp::Prefix> prefixes
            {
                { { QuarkFlavor::bottom, QuarkFlavor::up },      qnp::Prefix("B")   },
                { { QuarkFlavor::bottom, QuarkFlavor::strange }, qnp::Prefix("B_s") }
            };

            auto it = prefixes.find(std::make_tuple(opt_Q.value(), opt_q.value()));
            if (it == prefixes.end())
                throw InternalError("Combination of options Q=" + opt_Q.str() + ", q=" + opt_q.str() + " is not supported");

            return QualifiedName(it->second, qnp::Name(_name)).str();
        }

        std::string
        FLvD2022Profile::flvd2022_parameter(const char * _name) const
        {
            // The leading-twist parameters are shared with the FLvD2022 model,
            // cf. FLvD2022::parameter().
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, qnp::Prefix> prefixes
            {
                { { QuarkFlavor::bottom, QuarkFlavor::up },      qnp::Prefix("B_u") },
                { { QuarkFlavor::bottom, QuarkFlavor::strange }, qnp::Prefix("B_s") }
            };

            auto it = prefixes.find(std::make_tuple(opt_Q.value(), opt_q.value()));
            if (it == prefixes.end())
                throw InternalError("Combination of options Q=" + opt_Q.str() + ", q=" + opt_q.str() + " is not supported");

            return QualifiedName(it->second, qnp::Name(_name), qnp::Suffix("FLvD2022")).str();
        }

        HeavyMesonLCDAs *
        FLvD2022Profile::make(const Parameters & parameters, const Options & options)
        {
            return new FLvD2022Profile(parameters, options);
        }

        std::tuple<HeavyMesonLCDAs::CoefficientIterator, HeavyMesonLCDAs::CoefficientIterator>
        FLvD2022Profile::coefficient_range(const double & mu) const
        {
            // The leading-twist convolutions are unchanged with respect to FLvD2022, including
            // the RG evolution of the coefficients from mu_0 to mu.
            return leading_twist.coefficient_range(mu);
        }

        /* Leading twist two-particle LCDAs */

        double
        FLvD2022Profile::phi_plus(const double & omega) const
        {
            // cf. [FLvD:2022A], the momentum-space expansion labelled (paramom),
            //
            //   phi_+(omega) = omega exp(-x) / omega_0^2 sum_k a_k / (1 + k) L_k^(1)(2 x) ,
            //
            // with x = omega / omega_0 and L_k^(1) the associated Laguerre polynomials.
            // Evaluated at mu_0, cf. the comment on xi_1() for the choice of scale.
            return omega * this->profile_function(omega).f;
        }

        FLvD2022Profile::ProfileFunction
        FLvD2022Profile::profile_function(const double & omega) const
        {
            // f(omega) = phi_+(omega) / omega, the profile function of [BBJW:2018A],
            // Eq. (T4ansatz), and
            //
            //   f'(omega) = exp(-x) / omega_0^3 sum_k a_k / (1 + k)
            //               [-L_k^(1)(2 x) - 2 L_{k-1}^(2)(2 x)] ,
            //
            // using d/dy L_n^(a)(y) = -L_{n-1}^(a+1)(y) and L_{-1}^(2) = 0. At
            // a_k = (1, 0, ..., 0) they reduce to exp(-omega / omega_0) / omega_0^2 and
            // -f(omega) / omega_0, i.e. to the exponential model.
            const double x = omega / omega_0;

            // L_k^(1)(2 x) and L_k^(2)(2 x) by the recurrence
            // (k + 1) L_{k+1}^(a)(y) = (2 k + 1 + a - y) L_k^(a)(y) - (k + a) L_{k-1}^(a)(y).
            Weights L1, L2;
            L1[0] = 1.0;
            L1[1] = 2.0 - 2.0 * x;
            L2[0] = 1.0;
            L2[1] = 3.0 - 2.0 * x;
            for (unsigned k = 1 ; k + 1 < L1.size() ; ++k)
            {
                L1[k + 1] = ((2.0 * k + 2.0 - 2.0 * x) * L1[k] - (k + 1.0) * L1[k - 1]) / (k + 1.0);
                L2[k + 1] = ((2.0 * k + 3.0 - 2.0 * x) * L2[k] - (k + 2.0) * L2[k - 1]) / (k + 1.0);
            }

            Weights c_f, c_f_d1;
            for (unsigned k = 0 ; k < c_f.size() ; ++k)
            {
                const double L2_km1 = (k == 0) ? 0.0 : L2[k - 1];
                c_f[k]    = L1[k] / (1.0 + k);
                c_f_d1[k] = (-L1[k] - 2.0 * L2_km1) / (1.0 + k);
            }

            auto [a_begin, a_end] = leading_twist.coefficient_range(mu_0);
            const double exp_x = std::exp(-x);

            return {
                exp_x / power_of<2>(omega_0()) * std::inner_product(a_begin, a_end, c_f.begin(),    0.0),
                exp_x / power_of<3>(omega_0()) * std::inner_product(a_begin, a_end, c_f_d1.begin(), 0.0)
            };
        }

        double
        FLvD2022Profile::eom_moment() const
        {
            // (1/6) int_0^infinity domega omega^2 phi_+(omega), the left-hand side of the EOM
            // relation that fixes varkappa in [BBJW:2018A], Eq. (T4ansatz). With
            // I_k = int_0^infinity dx x^3 exp(-x) L_k^(1)(2 x) the generating function is
            //
            //   sum_k I_k t^k = 6 (1 - t)^2 (1 + t)^-4 ,
            //
            // so I_k / 6 = b_k - 2 b_{k-1} + b_{k-2} with b_m = (-1)^m binom(m + 3, 3), and
            //
            //   (1/6) int domega omega^2 phi_+ = omega_0^2 sum_k a_k (I_k / 6) / (1 + k) .
            //
            // At a_k = (1, 0, ..., 0) this is omega_0^2, the value of the exponential model.
            const Weights c = {
                1.0, -6.0 / 2.0, 19.0 / 3.0, -44.0 / 4.0, 85.0 / 5.0, -146.0 / 6.0,
                231.0 / 7.0, -344.0 / 8.0, 489.0 / 9.0
            };

            auto [a_begin, a_end] = leading_twist.coefficient_range(mu_0);
            return power_of<2>(omega_0()) * std::inner_product(a_begin, a_end, c.begin(), 0.0);
        }

        double
        FLvD2022Profile::eom_parameter_combination() const
        {
            // LambdaBar^2 + (2 lambda_E^2 + lambda_H^2) / 6, the right-hand side of the same
            // relation.
            return power_of<2>(lambda_bar()) + (2.0 * lambda_E2 + lambda_H2) / 6.0;
        }

        double
        FLvD2022Profile::varkappa() const
        {
            // Taken from the parameter combination rather than from the moment: in the
            // [FLvD:2022A] basis the moment is a linear functional of the a_k and need not be
            // positive, so a moment-defined varkappa does not always exist. Where the moment is
            // not positive the ansatz is inconsistent with its own construction, and returning a
            // number would hide that; see the class documentation.
            const double moment = this->eom_moment();

            if (moment <= 0.0)
            {
                auto [a_begin, a_end] = leading_twist.coefficient_range(mu_0);
                std::string coefficients;
                for (auto a = a_begin ; a != a_end ; ++a)
                {
                    coefficients += (coefficients.empty() ? "" : ", ") + stringify(*a, 8);
                }

                throw InternalError("The profile-function ansatz for the higher-twist B-meson LCDAs"
                        " requires a positive moment (1/6) int domega omega^2 phi_+(omega), but it"
                        " evaluates to " + stringify(moment, 8) + " GeV^2 for the coefficients a_k = ("
                        + coefficients + ") at mu_0 = " + stringify(mu_0(), 4) + " GeV and omega_0 = "
                        + stringify(omega_0(), 4) + " GeV");
            }

            return 1.0 / this->eom_parameter_combination();
        }

        double
        FLvD2022Profile::phi_minusWW(const double & omega) const
        {
            // The Wandzura-Wilczek relation is a linear functional of phi_+ and is therefore
            // defined by the [FLvD:2022A] parametrisation itself; nothing is modelled here.
            // In the ansatz it is phi_-^WW(omega) = int_omega^infinity drho f(rho),
            // cf. [BBJW:2018A], Eq. (eq:wwpartsinf).
            return leading_twist.phi_minusWW(omega);
        }

        double
        FLvD2022Profile::inverse_lambda_plus() const
        {
            return leading_twist.inverse_lambda_plus();
        }

        double
        FLvD2022Profile::xi_1(const double & omega) const
        {
            // cf. [BBJW:2018A], Eq. (4.15), in the equivalent profile-function form following
            // Eq. (A.39). Evaluated at mu_0: [BBJW:2018A], Sec. 5, evaluate the higher-twist
            // contributions and the related soft corrections at the reference scale.
            const auto [f, f_prime]   = this->profile_function(omega);
            const double phi_minus_ww = this->phi_minusWW(omega);
            const double kappa        = this->varkappa();
            const double bracket      = omega * omega * f - 2.0 * omega * phi_minus_ww;

            return 2.0 / 3.0 * kappa * (lambda_E2 + 2.0 * lambda_H2) * bracket
                - 2.0 * omega * phi_minus_ww + 3.0 * omega * omega * f + power_of<3>(omega) * f_prime;
        }

        double
        FLvD2022Profile::xi_2(const double & omega) const
        {
            // cf. [BBJW:2018A], Eq. (4.16), in the equivalent profile-function form following
            // Eq. (A.39). See xi_1() for the choice of scale.
            const double f            = this->profile_function(omega).f;
            const double phi_minus_ww = this->phi_minusWW(omega);
            const double kappa        = this->varkappa();
            const double bracket      = omega * omega * f - 2.0 * omega * phi_minus_ww;

            return -2.0 / 3.0 * kappa * (lambda_E2 - lambda_H2) * bracket
                + (lambda_bar - omega) * omega * f - omega * phi_minus_ww;
        }

        /*
         * The remaining accessors are not supplied by the ansatz for a generic profile
         * function; cf. the class documentation. They throw rather than fall back to the
         * exponential model, which would mix two different phi_+.
         */

        double
        FLvD2022Profile::phi_minus(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_minusWW(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_minusWW_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_minusWW_d2(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_plus(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_plus_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_plus_d2(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_bar(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_bar_d1(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_bar_d2(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::g_bar_d3(const double & /* omega */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            // The ansatz does give phi_3 in closed form, but in the normalisation convention of
            // [BBJW:2018A], Eq. (T4ansatz), which does not coincide with the [KMO:2006A]
            // convention used by this interface elsewhere. Only the integrated combinations
            // entering xi_1() and xi_2() are established here, so the individual three-particle
            // LCDAs are left unimplemented rather than converted on an unverified assumption.
            throw InternalError("The profile-function ansatz supplies the three-particle LCDAs only"
                    " through the kernels Xi_1 and Xi_2; phi_3 is not implemented");
        }

        double
        FLvD2022Profile::phi_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("The profile-function ansatz supplies the three-particle LCDAs only"
                    " through the kernels Xi_1 and Xi_2; phi_4 is not implemented");
        }

        double
        FLvD2022Profile::phi_bar_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar2_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar2_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar_bar_3(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::phi_bar_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::psi_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::psi_bar_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::chi_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::chi_bar_bar_4(const double & /* omega_1 */, const double & /* omega_2 */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::psi_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::psi_V(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::X_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::Y_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::Xbar_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        double
        FLvD2022Profile::Ybar_A(const double & /* omega */, const double & /* xi */) const
        {
            throw InternalError("Function not yet implemented");
        }

        Diagnostics
        FLvD2022Profile::diagnostics() const
        {
            Diagnostics results;

            // The two sides of the EOM relation that fixes varkappa, and their ratio. The ratio
            // is 1 exactly on the codimension-one surface where the ansatz is consistent, and
            // deviates from 1 by the amount by which this model violates it. Reported rather
            // than enforced, so that the violation is testable instead of silent.
            const double moment      = this->eom_moment();
            const double combination = this->eom_parameter_combination();

            results.add({ moment,                 "eom_moment()"                    });
            results.add({ combination,            "eom_parameter_combination()"     });
            results.add({ moment / combination,   "eom_moment() / eom_parameter_combination()" });

            return results;
        }

        const std::set<ReferenceName>
        FLvD2022Profile::references
        {
            "BBJW:2018A"_rn,
            "FLvD:2022A"_rn
        };

        const std::vector<OptionSpecification>
        FLvD2022Profile::options
        {
            { "Q"_ok,       { "b"_ov },                   "b"_ov        },
            { "q"_ok,       { "u"_ov, "s"_ov },           "u"_ov        },
            { "gminus"_ok,  { "zero"_ov, "WW-limit"_ov }, "WW-limit"_ov },
            { "alpha-s"_ok, { "naive"_ov, "full"_ov  },   "full"_ov     },
        };

        std::vector<OptionSpecification>::const_iterator
        FLvD2022Profile::begin_options()
        {
            return options.cbegin();
        }

        std::vector<OptionSpecification>::const_iterator
        FLvD2022Profile::end_options()
        {
            return options.cend();
        }
    }
}
