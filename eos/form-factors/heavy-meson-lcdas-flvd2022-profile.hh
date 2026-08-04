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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_HEAVY_MESON_LCDAS_FLVD2022_PROFILE_HH
#define EOS_GUARD_EOS_FORM_FACTORS_HEAVY_MESON_LCDAS_FLVD2022_PROFILE_HH 1

#include <eos/form-factors/heavy-meson-lcdas-flvd2022.hh>
#include <eos/form-factors/heavy-meson-lcdas.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <array>
#include <string>
#include <tuple>

namespace eos
{
    namespace heavy_meson_lcdas
    {
        /*!
         * A [FLvD:2022A] leading-twist LCDA phi_+ combined with three-particle LCDAs
         * *modelled* from it through the profile-function ansatz of [BBJW:2018A],
         * Eq. (T4ansatz) and the closed forms following Eq. (A.39).
         *
         * This is a model of its own, not a completion of [FLvD:2022A]. That
         * parametrisation defines phi_+ (and its Wandzura-Wilczek descendant phi_-^WW)
         * only; it defines neither phi_3 nor [psi_4 + psitilde_4]. Here those are
         * *generated* from phi_+ by identifying the [BBJW:2018A] profile function with
         * f(omega) = phi_+(omega) / omega. The resulting xi^soft_(tw-3,4) therefore
         * carries an unquantifiable *modelling* error, not the controlled truncation
         * error that is the point of the [FLvD:2022A] expansion.
         *
         * The ansatz determines the three-particle input for a generic profile function
         * only in the integrated combinations that enter the kernels Xi_1 and Xi_2, cf.
         * [BBJW:2018A] below Eq. (A.39): "for [psi_4 - psitilde_4] only the integral
         * [...] can be determined for a generic profile function f(omega) without
         * additional assumptions. Luckily, only this integral is relevant for
         * B -> gamma l nu_l." Accordingly xi_1() and xi_2() are implemented while the
         * accessors for the individual three-particle LCDAs are not.
         *
         * *The EOM caveat.* The ansatz fixes its normalisation constant varkappa
         * through the equation of motion to
         *
         *     varkappa^-1 = (1/6) int domega omega^2 phi_+(omega)
         *                 = LambdaBar^2 + (2 lambda_E^2 + lambda_H^2) / 6 .
         *
         * In the [FLvD:2022A] basis the left-hand side is a linear functional of the
         * coefficients a_k while the right-hand side is a fixed combination of
         * parameters, so the two agree only on a codimension-one surface in coefficient
         * space; away from it this model violates the very EOM relation the ansatz was
         * built to satisfy. varkappa is taken from the parameter combination -- the only
         * choice that is always defined, and the choice the exponential model makes --
         * and both sides are reported by diagnostics() so that their ratio is visible
         * and testable, cf. eom_moment() and eom_parameter_combination(). Where the
         * moment is not positive the ansatz is not merely inaccurate but inconsistent
         * with its own construction, and xi_1() and xi_2() throw rather than return a
         * number.
         *
         * All higher-twist accessors are evaluated at mu_0, cf. [BBJW:2018A], Sec. 5.
         */
        class FLvD2022Profile : public HeavyMesonLCDAs
        {
            private:
                static const unsigned int number_of_parameters = 9u;
                using Weights                                  = std::array<double, number_of_parameters>; // We implement the weights as fixed-size arrays

                // The leading-twist half of the model, held by composition. It supplies the
                // coefficients a_k including their RG evolution, phi_-^WW, and 1/lambda_B.
                // Composition rather than inheritance: FLvD2022 declares its overrides final,
                // and this class must not present itself as a variant of it.
                FLvD2022 leading_twist;

                QuarkFlavorOption opt_Q;
                QuarkFlavorOption opt_q;

                UsedParameter mu_0;
                UsedParameter omega_0;

                // The higher-twist parameters use the same names as the exponential model, so
                // that at the default coefficients the two models coincide term by term.
                UsedParameter lambda_bar;
                UsedParameter lambda_E2;
                UsedParameter lambda_H2;

                std::string parameter(const char * _name) const;
                std::string flvd2022_parameter(const char * _name) const;

                /*!
                 * The profile function f(omega) = phi_+(omega) / omega of [BBJW:2018A],
                 * Eq. (T4ansatz), and its first derivative, at mu_0.
                 *
                 * Both are linear functionals of the coefficients a_k(mu_0) and both are needed at
                 * every point by xi_1(), so they are evaluated together from a single read of the
                 * coefficients and a single pass of the Laguerre recurrences.
                 */
                struct ProfileFunction
                {
                        double f;
                        double f_d1;
                };

                ProfileFunction profile_function(const double & omega) const;

                /*!
                 * The normalisation constant varkappa of the ansatz, taken from the parameter
                 * combination. Throws if eom_moment() is not positive.
                 */
                double varkappa() const;

                static const std::vector<OptionSpecification> options;

            public:
                FLvD2022Profile(const Parameters & parameters, const Options & options);
                virtual ~FLvD2022Profile() = default;

                static HeavyMesonLCDAs * make(const Parameters & parameters, const Options & options);

                /*!
                 * Parmeters of the B-Meson LCDA phi+ as defined in Ref. [FLvD:2022A]
                 *
                 * mu: the renormalization scale
                 */
                virtual std::tuple<HeavyMesonLCDAs::CoefficientIterator, HeavyMesonLCDAs::CoefficientIterator> coefficient_range(const double & mu) const final override;

                /*!
                 * Leading twist two-particle LCDAs
                 *
                 * omega: plus-component of the spectator momentum
                 */
                virtual double phi_plus(const double & omega) const final override;

                /*!
                 * Twist-three two-particle LCDAs
                 *
                 * omega: plus-component of the spectator momentum
                 */
                virtual double phi_minus(const double & omega) const final override;
                virtual double phi_minusWW(const double & omega) const final override;
                virtual double phi_bar(const double & omega) const final override;
                virtual double phi_bar_d1(const double & omega) const final override;

                /*!
                 * Next-to-leading twist two-particle LCDAs
                 *
                 * omega: plus-component of the spectator momentum
                 */
                virtual double g_plus(const double & omega) const final override;
                virtual double g_plus_d1(const double & omega) const final override;
                virtual double g_plus_d2(const double & omega) const final override;

                virtual double g_minusWW(const double & omega) const final override;
                virtual double g_minusWW_d1(const double & omega) const final override;
                virtual double g_minusWW_d2(const double & omega) const final override;

                virtual double g_bar(const double & omega) const final override;
                virtual double g_bar_d1(const double & omega) const final override;
                virtual double g_bar_d2(const double & omega) const final override;
                virtual double g_bar_d3(const double & omega) const final override;

                /*!
                 * Leading power three-particle LCDAs
                 *
                 * omega_1: plus-component of the spectator momentum
                 * omega_2: plus-component of the gluon momentum
                 */
                virtual double phi_3(const double & omega_1, const double & omega_2) const final override;
                virtual double phi_4(const double & omega_1, const double & omega_2) const final override;

                virtual double phi_bar_3(const double & omega_1, const double & omega_2) const final override;
                virtual double phi_bar_4(const double & omega_1, const double & omega_2) const final override;

                virtual double phi_bar2_3(const double & omega_1, const double & omega_2) const final override;
                virtual double phi_bar2_4(const double & omega_1, const double & omega_2) const final override;

                virtual double phi_bar_bar_3(const double & omega_1, const double & omega_2) const final override;
                virtual double phi_bar_bar_4(const double & omega_1, const double & omega_2) const final override;

                virtual double psi_bar_4(const double & omega_1, const double & omega_2) const final override;
                virtual double chi_bar_4(const double & omega_1, const double & omega_2) const final override;

                virtual double psi_bar_bar_4(const double & omega_1, const double & omega_2) const final override;
                virtual double chi_bar_bar_4(const double & omega_1, const double & omega_2) const final override;

                /*!
                 * Higher-twist kernels entering B -> gamma form factors.
                 *
                 * cf. [BBJW:2018A], Eqs. (4.15) and (4.16), in the profile-function form
                 * following Eq. (A.39). Modelled, see the class documentation.
                 */
                virtual double xi_1(const double & omega) const final override;
                virtual double xi_2(const double & omega) const final override;

                /*!
                 * Pseudo observables for the two-particle LCDAs
                 */
                virtual double inverse_lambda_plus() const final override;

                /*!
                 * The two sides of the EOM relation that fixes varkappa.
                 *
                 * eom_moment() is (1/6) int domega omega^2 phi_+(omega), a linear functional
                 * of the coefficients a_k at mu_0; eom_parameter_combination() is
                 * LambdaBar^2 + (2 lambda_E^2 + lambda_H^2) / 6. The ansatz is consistent only
                 * where the two agree.
                 */
                double eom_moment() const;
                double eom_parameter_combination() const;

                /*!
                 * Leading power three-particle LCDAs
                 *
                 * omega: plus-component of the spectator momentum
                 * xi:    plus-component of the gluon momentum
                 */
                virtual double psi_A(const double & omega, const double & xi) const final override;
                virtual double psi_V(const double & omega, const double & xi) const final override;
                virtual double X_A(const double & omega, const double & xi) const final override;
                virtual double Y_A(const double & omega, const double & xi) const final override;

                /*!
                 * Auxiliary functions for the three-particle LCDAs
                 *
                 * See [KMO:2006A], below eq. (72), p. 28 for their definition.
                 */
                virtual double Xbar_A(const double & omega, const double & xi) const final override;
                virtual double Ybar_A(const double & omega, const double & xi) const final override;

                /* Internal diagnostics */
                virtual Diagnostics diagnostics() const final override;

                /*!
                 * References used in the computation of our (pseudo)observables.
                 */
                static const std::set<ReferenceName> references;

                /*!
                 * Options used in the computation of our (pseudo)observables.
                 */
                static std::vector<OptionSpecification>::const_iterator begin_options();
                static std::vector<OptionSpecification>::const_iterator end_options();
        };
    } // namespace heavy_meson_lcdas
} // namespace eos

#endif
