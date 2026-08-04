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

#include <eos/form-factors/heavy-meson-lcdas-exponential.hh>
#include <eos/form-factors/heavy-meson-lcdas-flvd2022-profile.hh>
#include <eos/form-factors/heavy-meson-lcdas-flvd2022.hh>
#include <eos/utils/diagnostics.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <tuple>
#include <vector>

using namespace test;
using namespace eos;

class HeavyMesonLCDAsFLvD2022ProfileTest : public TestCase
{
    public:
        HeavyMesonLCDAsFLvD2022ProfileTest() :
            TestCase("heavy_meson_lcdas_flvd2022_profile_test")
        {
        }

        static void
        set_generic_coefficients(Parameters & p)
        {
            // The "generic B" vector of issues/1206/reference-profile-ansatz.py. Chosen because
            // its EOM moment is positive, so the ansatz can be evaluated at all; cf. the
            // negative-moment block below for the other case.
            p["B_u::a^phi+_0@FLvD2022"] = 0.80;
            p["B_u::a^phi+_1@FLvD2022"] = 0.25;
            p["B_u::a^phi+_2@FLvD2022"] = -0.15;
            p["B_u::a^phi+_3@FLvD2022"] = 0.05;
            p["B_u::a^phi+_4@FLvD2022"] = 0.10;
            p["B_u::a^phi+_5@FLvD2022"] = -0.08;
            p["B_u::a^phi+_6@FLvD2022"] = 0.04;
            p["B_u::a^phi+_7@FLvD2022"] = -0.02;
            p["B_u::a^phi+_8@FLvD2022"] = 0.01;
        }

        virtual void
        run() const
        {
            static const std::array<double, 7> omegas{ 0.0, 1.0e-6, 0.05, 0.15, 0.30, 0.60, 1.20 };

            /*
             * At the default coefficients a_k = (1, 0, ..., 0) the [FLvD:2022A] phi_+ collapses to
             * the exponential shape, so every accessor of this model must reproduce the
             * exponential model's. **This is a wiring test, not a physics test**: per
             * issues/1206/STEP-07-FINDING.md Sec. 2.3 it passes by construction, because the
             * profile-function kernels reduce identically to the exponential ones at these
             * coefficients. It pins the implementation, and says nothing about whether the
             * modelled three-particle LCDAs are right.
             */
            {
                Parameters p               = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"] = 0.35;
                p["B_u::mu_0@FLvD2022"]    = 1.0;

                heavy_meson_lcdas::FLvD2022Profile profile(p, Options{});
                // lambda-b-source=FLvD2022 makes the exponential model read the same
                // inverse-moment scale, B_u::omega_0@FLvD2022.
                heavy_meson_lcdas::Exponential     exponential(p,
                                                           Options{
                                                                   { "lambda-b-source"_ok, "FLvD2022"_ov }
                });

                for (const double & omega : omegas)
                {
                    TEST_CHECK_NEARLY_EQUAL(profile.phi_plus(omega), exponential.phi_plus(omega), 1e-14);
                    TEST_CHECK_NEARLY_EQUAL(profile.phi_minusWW(omega), exponential.phi_minusWW(omega), 1e-14);
                    TEST_CHECK_NEARLY_EQUAL(profile.xi_1(omega), exponential.xi_1(omega), 1e-14);
                    TEST_CHECK_NEARLY_EQUAL(profile.xi_2(omega), exponential.xi_2(omega), 1e-14);
                }

                TEST_CHECK_NEARLY_EQUAL(profile.inverse_lambda_plus(), exponential.inverse_lambda_plus(), 1e-14);

                // At the default coefficients the EOM moment is omega_0^2, the exponential-model
                // value, independently of the parameter combination it is compared against.
                TEST_CHECK_NEARLY_EQUAL(profile.eom_moment(), 0.35 * 0.35, 1e-15);
            }

            /*
             * The closed forms at a generic coefficient vector, against
             * issues/1206/reference-profile-ansatz.py, which implements phi_+ from
             * [FLvD:2022A], Eq. (paramom), and Xi_1, Xi_2 from [BBJW:2018A], Eqs. (4.15), (4.16)
             * in the profile form following Eq. (A.39), without importing EOS. Its f' is checked
             * there against finite differences of phi_+(omega)/omega and its EOM moment against
             * quadrature.
             */
            {
                Parameters p               = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"] = 0.30;
                p["B_u::mu_0@FLvD2022"]    = 1.0;
                set_generic_coefficients(p);

                heavy_meson_lcdas::FLvD2022Profile profile(p, Options{});

                for (const auto & [omega, phi_plus, phi_minusWW, Xi_1, Xi_2] : std::array<std::tuple<double, double, double, double, double>, 3>{
                         { { 0.05, +0.450374811056351, +2.088567238043244, -0.214969597373140, +0.084437489499010 },
                          { 0.30, +1.004027557485672, +0.666456537763827, +0.115393320317877, -0.006440016899613 },
                          { 1.20, -0.163491357955376, -0.050092895495690, -1.086989832160827, +0.168928263655681 } }
                })
                {
                    TEST_CHECK_NEARLY_EQUAL(profile.phi_plus(omega), phi_plus, 1e-13);
                    TEST_CHECK_NEARLY_EQUAL(profile.phi_minusWW(omega), phi_minusWW, 1e-13);
                    TEST_CHECK_NEARLY_EQUAL(profile.xi_1(omega), Xi_1, 1e-13);
                    TEST_CHECK_NEARLY_EQUAL(profile.xi_2(omega), Xi_2, 1e-13);
                }

                // The EOM caveat, made numerical. The moment is a linear functional of the a_k,
                // the parameter combination is not, and here they disagree by a factor 1.64: away
                // from the codimension-one surface where they agree, this model violates the very
                // EOM relation the ansatz was constructed to satisfy. varkappa is nevertheless
                // taken from the parameter combination, the only choice that always exists.
                TEST_CHECK_NEARLY_EQUAL(profile.eom_moment(), 0.4428, 1e-13);
                TEST_CHECK_NEARLY_EQUAL(profile.eom_parameter_combination(), 0.2700, 1e-15);

                const Diagnostics d     = profile.diagnostics();
                auto              entry = d.begin();
                TEST_CHECK_NEARLY_EQUAL((entry++)->value, 0.4428, 1e-13);          // eom_moment()
                TEST_CHECK_NEARLY_EQUAL((entry++)->value, 0.2700, 1e-15);          // eom_parameter_combination()
                TEST_CHECK_NEARLY_EQUAL((entry++)->value, 0.4428 / 0.2700, 1e-12); // their ratio
                TEST_CHECK(entry == d.end());
            }

            /*
             * The consistent case: the ratio is 1 exactly on the surface where the moment equals
             * the parameter combination. At the default coefficients the moment is omega_0^2, so
             * choosing LambdaBar^2 = omega_0^2 - (2 lambda_E^2 + lambda_H^2)/6 puts the model
             * there. This is the only configuration in which the ansatz is EOM-consistent.
             */
            {
                Parameters p               = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"] = 0.50;
                p["B_u::mu_0@FLvD2022"]    = 1.0;
                p["B::lambda_E^2"]         = 0.03;
                p["B::lambda_H^2"]         = 0.06;
                p["B::LambdaBar"]          = std::sqrt(0.50 * 0.50 - (2.0 * 0.03 + 0.06) / 6.0);

                heavy_meson_lcdas::FLvD2022Profile profile(p, Options{});

                TEST_CHECK_NEARLY_EQUAL(profile.eom_moment() / profile.eom_parameter_combination(), 1.0, 1e-14);
            }

            /*
             * Where the moment is not positive, a moment-defined varkappa does not exist and the
             * ansatz is inconsistent with its own construction rather than merely inaccurate. The
             * kernels must throw instead of returning a number. The coefficients are the
             * "generic A" vector of issues/1206/reference-profile-ansatz.py, for which the moment
             * is -0.06525 GeV^2 at omega_0 = 0.30 GeV.
             */
            {
                Parameters p                = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"]  = 0.30;
                p["B_u::mu_0@FLvD2022"]     = 1.0;
                p["B_u::a^phi+_0@FLvD2022"] = 1.000;
                p["B_u::a^phi+_1@FLvD2022"] = -0.350;
                p["B_u::a^phi+_2@FLvD2022"] = 0.200;
                p["B_u::a^phi+_3@FLvD2022"] = 0.100;
                p["B_u::a^phi+_4@FLvD2022"] = -0.050;
                p["B_u::a^phi+_5@FLvD2022"] = 0.030;
                p["B_u::a^phi+_6@FLvD2022"] = -0.020;
                p["B_u::a^phi+_7@FLvD2022"] = 0.010;
                p["B_u::a^phi+_8@FLvD2022"] = -0.005;

                heavy_meson_lcdas::FLvD2022Profile profile(p, Options{});

                TEST_CHECK_NEARLY_EQUAL(profile.eom_moment(), -0.06525, 1e-13);
                // phi_+ and phi_-^WW are unaffected: they are defined by the parametrisation, not
                // by the ansatz.
                TEST_CHECK_NO_THROW(profile.phi_plus(0.30));
                TEST_CHECK_NO_THROW(profile.phi_minusWW(0.30));
                // The modelled kernels are not.
                TEST_CHECK_THROWS(InternalError, profile.xi_1(0.30));
                TEST_CHECK_THROWS(InternalError, profile.xi_2(0.30));
            }

            /*
             * The leading-twist half is shared with FLvD2022 verbatim, evolution included, while
             * FLvD2022 itself must keep throwing for the modelled kernels: this step adds a model,
             * it does not complete that parametrisation.
             */
            {
                Parameters p               = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"] = 0.30;
                p["B_u::mu_0@FLvD2022"]    = 1.0;
                set_generic_coefficients(p);

                heavy_meson_lcdas::FLvD2022Profile profile(p, Options{});
                heavy_meson_lcdas::FLvD2022        flvd2022(p, Options{});

                for (const double & mu : std::array<double, 3>{ 1.0, 1.5, 2.5 })
                {
                    // The two ranges must be copied out before the second call: FLvD2022 returns
                    // iterators into a function-local thread_local buffer, so holding both ranges
                    // at once would compare that buffer with itself and pass vacuously.
                    std::vector<double> from_profile;
                    {
                        auto [begin, end] = profile.coefficient_range(mu);
                        from_profile.assign(begin, end);
                    }
                    std::vector<double> from_flvd2022;
                    {
                        auto [begin, end] = flvd2022.coefficient_range(mu);
                        from_flvd2022.assign(begin, end);
                    }

                    TEST_CHECK_EQUAL(from_profile.size(), from_flvd2022.size());
                    TEST_CHECK(from_profile.size() > 0u);
                    for (unsigned k = 0; k < from_profile.size(); ++k)
                    {
                        TEST_CHECK_NEARLY_EQUAL(from_profile[k], from_flvd2022[k], 1e-15);
                    }
                }

                for (const double & omega : omegas)
                {
                    TEST_CHECK_NEARLY_EQUAL(profile.phi_minusWW(omega), flvd2022.phi_minusWW(omega), 1e-15);
                }
                TEST_CHECK_NEARLY_EQUAL(profile.inverse_lambda_plus(), flvd2022.inverse_lambda_plus(), 1e-15);

                TEST_CHECK_THROWS(InternalError, flvd2022.xi_1(0.30));
                TEST_CHECK_THROWS(InternalError, flvd2022.xi_2(0.30));
            }

            /* The model is reachable through the factory under its own name. */
            {
                Parameters p               = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"] = 0.35;
                p["B_u::mu_0@FLvD2022"]    = 1.0;

                auto lcdas = HeavyMesonLCDAs::make("FLvD2022+profile", p, Options{});
                TEST_CHECK(nullptr != lcdas.get());

                heavy_meson_lcdas::Exponential exponential(p,
                                                           Options{
                                                               { "lambda-b-source"_ok, "FLvD2022"_ov }
                });
                TEST_CHECK_NEARLY_EQUAL(lcdas->xi_1(0.30), exponential.xi_1(0.30), 1e-14);
                TEST_CHECK_NEARLY_EQUAL(lcdas->xi_2(0.30), exponential.xi_2(0.30), 1e-14);
            }
        }
} heavy_meson_lcdas_flvd2022_profile_test;
