/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022-2026 Danny van Dyk
 * Copyright (c) 2022-2024 Philip Lüghausen
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

#include <eos/form-factors/analytic-p-to-gamma-qcdf.hh>
#include <eos/models/model.hh>
#include <eos/observable.hh>
#include <eos/utils/diagnostics.hh>

#include <test/test.hh>

#include <array>
#include <cmath>

using namespace test;
using namespace eos;

class AnalyticFormFactorBToGammaQCDFTest : public TestCase
{
    public:
        AnalyticFormFactorBToGammaQCDFTest() :
            TestCase("analytic_b_to_gamma_qcdf_test")
        {
        }

        virtual void
        run() const
        {
            // Regression tests at the inputs of [BBJW:2018A], Table 1, with sigmahat_1 = 0.
            {
                Parameters p                      = Parameters::Defaults();
                p["QCD::alpha_s(MZ)"]             = 0.117400;
                p["mass::b(MSbar)"]               = 4.458251328586803;
                p["mass::B_u"]                    = 5.27929;
                p["mass::rho^+"]                  = 0.77526;
                p["decay-constant::B_u"]          = 0.192;
                p["B::LambdaBar"]                 = 0.47929;
                p["B::lambda_E^2"]                = 0.0625;
                p["B::lambda_H^2"]                = 0.125;
                p["B->gamma::mu@FLvD2022QCDF"]    = 1.5;
                p["B->gamma::mu_h1@FLvD2022QCDF"] = 4.8;
                p["B->gamma::mu_h2@FLvD2022QCDF"] = 4.8;
                p["B->gamma::s_0@FLvD2022QCDF"]   = 1.5;
                p["B->gamma::M^2@FLvD2022QCDF"]   = 1.25;

                auto model = Model::make("SM"_ov, p, Options());
                TEST_CHECK_NEARLY_EQUAL(model->m_b_pole(1), 4.8, 1e-6);

                // Positional access into the diagnostics of AnalyticFormFactorPToGammaQCDF. The
                // indices used below are, in the order in which diagnostics() emits them:
                //
                //    5 -> L0_effective(3.0)      6 -> L0_effective(2.16)
                //    9 -> C at Egamma=2.16      10 -> K_inv at Egamma=2.16    11 -> U at Egamma=2.16
                //
                // Inserting an entry into diagnostics() silently shifts these. That is caught by the
                // TEST_CHECK_DIAGNOSTICS blocks below, which pin the full ordered list and would fail
                // first -- keep them in sync when adding a diagnostic.
                const auto diagnostic = [](const Diagnostics & d, const unsigned i)
                {
                    auto entry = d.begin();
                    for (unsigned j = 0; j < i; ++j)
                    {
                        ++entry;
                    }
                    return entry->value;
                };

                const auto average = [](const AnalyticFormFactorPToGammaQCDF<BToGamma> & ff, const double & E_gamma) { return 0.5 * (ff.F_V(E_gamma) + ff.F_A(E_gamma)); };

                // [BR:2011A], Eq. (2.14): J = L_hat^2 - 1 for the exponential LCDA.
                for (const double & lambda_B : std::array<double, 2>{ 0.20, 0.35 })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p,
                                                                Options{
                                                                    { "lcda-model"_ok, "exponential"_ov }
                    });
                    const Diagnostics                        d = ff.diagnostics();

                    for (const auto & [E_gamma, index] : std::array<std::pair<double, unsigned>, 2>{
                             { { 3.0, 5 }, { 2.16, 6 } }
                    })
                    {
                        const double L_hat = 0.577215664901533 + std::log(1.5 * 1.5 / (2.0 * E_gamma * lambda_B));
                        TEST_CHECK_NEARLY_EQUAL(diagnostic(d, index) * lambda_B, L_hat * L_hat - 1.0, 1e-12);
                    }
                }

                // R / J = C K^-1 U: only J may depend on the LCDA.
                double perturbative_factor = 0.0;
                for (const double & lambda_B : std::array<double, 3>{ 0.20, 0.35, 0.50 })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p,
                                                                Options{
                                                                    { "lcda-model"_ok, "exponential"_ov }
                    });
                    const Diagnostics                        d     = ff.diagnostics();
                    const double                             value = diagnostic(d, 9) * diagnostic(d, 10) * diagnostic(d, 11);
                    if (lambda_B == 0.20)
                    {
                        perturbative_factor = value;
                    }
                    else
                    {
                        TEST_CHECK_NEARLY_EQUAL(value, perturbative_factor, 1e-12);
                    }
                }

                // [BBJW:2018A], Eq. (5.11) and Fig. 4: xi-hat^soft_(LO) at E_gamma = 2 GeV.
                for (const auto & [lambda_B, published] : std::array<std::pair<double, double>, 3>{
                         { { 0.20, -0.7875 }, { 0.35, -0.4085 }, { 0.50, -0.2625 } }
                })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    const Options common{
                        {      "lcda-model"_ok, "exponential"_ov },
                        { "evolution-order"_ok,          "LL"_ov }
                    };
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_none(p,
                                                                     common
                                                                             + Options{
                                                                                 { "contributions"_ok, "none"_ov }
                    });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_soft(p,
                                                                     common
                                                                             + Options{
                                                                                 { "contributions"_ok, "soft-nlo"_ov }
                    });
                    const double                             F_LP    = average(ff_none, 2.0);
                    const double                             xi_soft = average(ff_soft, 2.0) - F_LP;
                    TEST_CHECK_NEARLY_EQUAL(2.0 * 2.0 * xi_soft / F_LP, published, 1e-3);
                }

                // [BBJW:2018A], plots/sigma0.pdf, dash-dotted curve: xi^ht at E_gamma = 2 GeV.
                for (const auto & [lambda_B, published] : std::array<std::pair<double, double>, 3>{
                         { { 0.20, -0.07039 }, { 0.35, -0.08846 }, { 0.50, -0.09570 } }
                })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    const Options common{
                        { "lcda-model"_ok, "exponential"_ov }
                    };
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_none(p,
                                                                     common
                                                                             + Options{
                                                                                 { "contributions"_ok, "none"_ov }
                    });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_ht(p,
                                                                   common
                                                                           + Options{
                                                                               { "contributions"_ok, "ht"_ov }
                    });
                    TEST_CHECK_NEARLY_EQUAL(average(ff_ht, 2.0) - average(ff_none, 2.0), published, 1e-4);
                }

                // [BBJW:2018A], Eqs. (4.13), (4.15), (4.16): xi^soft_(tw-3,4) at E_gamma = 2 GeV.
                //
                // Checked for both LCDA models that supply the kernels Xi_1, Xi_2. For
                // "exponential" this pins the closed forms against the direct numerical
                // evaluation of Eqs. (4.15), (4.16) reported in issues/1206/DESIGN.md, STEP-04.
                // For "FLvD2022+profile" it is a **wiring test, not a physics test**: at the
                // default coefficients a_k = (1, 0, ..., 0) the profile-function ansatz reduces
                // identically to the exponential kernels, so the agreement holds by construction
                // (issues/1206/STEP-07-FINDING.md, Sec. 2.3). It pins the implementation and says
                // nothing about the modelled three-particle LCDAs.
                for (const auto & [lambda_B, expected] : std::array<std::pair<double, double>, 3>{
                         { { 0.20, +0.02511170 }, { 0.35, +0.01717038 }, { 0.50, +0.01198488 } }
                })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    p["B::1/lambda_B_p"]       = 1.0 / 0.46; // deliberately different: QCDF uses omega_0 above
                    for (const char * lcda_model : { "exponential", "FLvD2022+profile" })
                    {
                        const Options common{
                            { "lcda-model"_ok, qnp::OptionValue(lcda_model) }
                        };
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff_none(p,
                                                                         common
                                                                                 + Options{
                                                                                     { "contributions"_ok, "none"_ov }
                        });
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff_tw_3_4(p,
                                                                           common
                                                                                   + Options{
                                                                                       { "contributions"_ok, "soft-tw-3+4"_ov }
                        });
                        TEST_CHECK_NEARLY_EQUAL(average(ff_tw_3_4, 2.0) - average(ff_none, 2.0), expected, 1e-6);
                    }
                }

                // xi^soft_(tw-3,4) needs the kernels Xi_1, Xi_2 and hence the three-particle LCDAs,
                // which [FLvD:2022A] does not parametrise. Naming it under FLvD2022 must be
                // rejected rather than ignored. xi^soft_(tw-5,6) needs only phi_-^WW and
                // 1/lambda_B, which that parametrisation does define, so it is available for both
                // models and must not throw.
                TEST_CHECK_THROWS(InternalError,
                                  AnalyticFormFactorPToGammaQCDF<BToGamma>(p,
                                                                           Options{
                                                                               {    "lcda-model"_ok,    "FLvD2022"_ov },
                                                                               { "contributions"_ok, "soft-tw-3+4"_ov }
                }));
                // FLvD2022+profile models those LCDAs by the [BBJW:2018A] profile-function ansatz,
                // so the same request must be accepted there. That is a modelling choice, not a
                // completion of [FLvD:2022A]; see heavy-meson-lcdas-flvd2022-profile.hh.
                TEST_CHECK_NO_THROW(AnalyticFormFactorPToGammaQCDF<BToGamma>(p,
                                                                             Options{
                                                                                 {    "lcda-model"_ok, "FLvD2022+profile"_ov },
                                                                                 { "contributions"_ok,      "soft-tw-3+4"_ov }
                }));

                // Equal-model anchor: with matched lambda_B, default coefficients and mu_0 = mu, so
                // that no evolution occurs, both LCDA models must give the same xi^soft_(tw-5,6).
                // Unlike a kernel generated from phi_+ by a model ansatz, phi_-^WW is independently
                // defined by each model, so this compares two implementations of one object.
                {
                    Parameters q                   = p;
                    q["B_u::mu_0@FLvD2022"]        = 1.5;
                    q["B->gamma::mu@FLvD2022QCDF"] = 1.5;
                    q["B_u::omega_0@FLvD2022"]     = 0.35;
                    q["B::1/lambda_B_p"]           = 1.0 / 0.35;

                    const auto tw_5_6 = [&q](const char * lcda_model) -> double
                    {
                        const Options common{
                            { "lcda-model"_ok, qnp::OptionValue(lcda_model) }
                        };
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff_none(q,
                                                                         common
                                                                                 + Options{
                                                                                     { "contributions"_ok, "none"_ov }
                        });
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff_tw(q,
                                                                       common
                                                                               + Options{
                                                                                   { "contributions"_ok, "soft-tw-5+6"_ov }
                        });
                        return ff_tw.F_V(2.0) - ff_none.F_V(2.0);
                    };

                    TEST_CHECK_NEARLY_EQUAL(tw_5_6("FLvD2022"), tw_5_6("exponential"), 1e-12);
                }

                // FLvD2022 is the default LCDA model and must reproduce an explicit selection.
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_flvd_default(p, Options{});
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_flvd_all(p,
                                                                     Options{
                                                                         {    "lcda-model"_ok, "FLvD2022"_ov },
                                                                         { "contributions"_ok,      "all"_ov }
                });
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_flvd_ht(p,
                                                                    Options{
                                                                        {    "lcda-model"_ok, "FLvD2022"_ov },
                                                                        { "contributions"_ok,       "ht"_ov }
                });
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_flvd_soft(p,
                                                                      Options{
                                                                          {    "lcda-model"_ok, "FLvD2022"_ov },
                                                                          { "contributions"_ok,     "soft"_ov }
                });
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_flvd_none(p,
                                                                      Options{
                                                                          {    "lcda-model"_ok, "FLvD2022"_ov },
                                                                          { "contributions"_ok,     "none"_ov }
                });
                TEST_CHECK_NEARLY_EQUAL(ff_flvd_default.F_V(2.0), ff_flvd_all.F_V(2.0), 1e-12);
                TEST_CHECK_NEARLY_EQUAL(ff_flvd_default.F_A(2.0), ff_flvd_all.F_A(2.0), 1e-12);
                TEST_CHECK_NEARLY_EQUAL(average(ff_flvd_all, 2.0) - average(ff_flvd_ht, 2.0) - average(ff_flvd_soft, 2.0) + average(ff_flvd_none, 2.0), 0.0, 1e-12);

                // FLvD2022 must not have changed at all by the addition of FLvD2022+profile: the
                // two models share phi_+, its RG evolution, phi_-^WW and 1/lambda_B, and differ
                // *only* in the modelled twist-3,4 soft term. Asserted rather than assumed, and in
                // both form factors, since delta_xi carries the term too.
                {
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_profile_all(p,
                                                                            Options{
                                                                                {    "lcda-model"_ok, "FLvD2022+profile"_ov },
                                                                                { "contributions"_ok,              "all"_ov }
                    });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_profile_none(p,
                                                                             Options{
                                                                                 {    "lcda-model"_ok, "FLvD2022+profile"_ov },
                                                                                 { "contributions"_ok,             "none"_ov }
                    });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_profile_tw_3_4(p,
                                                                               Options{
                                                                                   {    "lcda-model"_ok, "FLvD2022+profile"_ov },
                                                                                   { "contributions"_ok,      "soft-tw-3+4"_ov }
                    });

                    // The leading power is untouched by the ansatz.
                    TEST_CHECK_NEARLY_EQUAL(ff_profile_none.F_V(2.0), ff_flvd_none.F_V(2.0), 1e-12);
                    TEST_CHECK_NEARLY_EQUAL(ff_profile_none.F_A(2.0), ff_flvd_none.F_A(2.0), 1e-12);
                    // ... and the whole difference in "all" is the twist-3,4 term.
                    TEST_CHECK_NEARLY_EQUAL(ff_profile_all.F_V(2.0) - ff_flvd_all.F_V(2.0), ff_profile_tw_3_4.F_V(2.0) - ff_profile_none.F_V(2.0), 1e-12);
                    TEST_CHECK_NEARLY_EQUAL(ff_profile_all.F_A(2.0) - ff_flvd_all.F_A(2.0), ff_profile_tw_3_4.F_A(2.0) - ff_profile_none.F_A(2.0), 1e-12);
                    // The term is non-zero, so the checks above are not vacuous.
                    TEST_CHECK(std::abs(ff_profile_all.F_V(2.0) - ff_flvd_all.F_V(2.0)) > 1e-4);
                }

                // Every contribution included in "all" must also be independently selectable, and
                // the composites must be exactly the sum of their parts. This is the guard that
                // catches a term being added to "all" without being added to the decomposition.
                // Checked for both models that support all four terms.
                p["B_u::omega_0@FLvD2022"] = 0.35;
                for (const char * lcda_model : { "exponential", "FLvD2022+profile" })
                {
                    const Options common{
                        { "lcda-model"_ok, qnp::OptionValue(lcda_model) }
                    };
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_sum_none(p,
                                                                         common
                                                                                 + Options{
                                                                                     { "contributions"_ok, "none"_ov }
                    });
                    const double                             F_V_none = ff_sum_none.F_V(2.0);

                    const auto term = [&](const char * contribution) -> double
                    {
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p,
                                                                    common
                                                                            + Options{
                                                                                { "contributions"_ok, qnp::OptionValue(contribution) }
                        });
                        return ff.F_V(2.0) - F_V_none;
                    };

                    // "soft" = xi^soft_(NLO) + xi^soft_(tw-3,4) + xi^soft_(tw-5,6)
                    TEST_CHECK_NEARLY_EQUAL(term("soft"), term("soft-nlo") + term("soft-tw-3+4") + term("soft-tw-5+6"), 1e-12);
                    // "all" = xi^ht + "soft"
                    TEST_CHECK_NEARLY_EQUAL(term("all"), term("ht") + term("soft"), 1e-12);
                    // ... and therefore "all" is the sum of the four individual terms.
                    TEST_CHECK_NEARLY_EQUAL(term("all"), term("ht") + term("soft-nlo") + term("soft-tw-3+4") + term("soft-tw-5+6"), 1e-12);
                }

                // Which terms each (lcda-model, contributions) pair actually enables. The four
                // switches are the last entries of diagnostics(). Without this, "soft" and
                // "soft-nlo" would be silently indistinguishable under FLvD2022. Three models
                // times seven values of "contributions" is twenty-one combinations; the one
                // remaining cell is FLvD2022 with "soft-tw-3+4", which throws and is covered
                // separately above, so twenty rows are pinned here.
                for (const auto & [lcda_model, contributions, expected] : std::array<std::tuple<const char *, const char *, std::array<double, 4>>, 20>{
                         { // lcda-model          contributions    ht   soft  tw3+4 tw5+6
                           { "exponential", "none", { 0.0, 0.0, 0.0, 0.0 } },
                          { "exponential", "ht", { 1.0, 0.0, 0.0, 0.0 } },
                          { "exponential", "soft-nlo", { 0.0, 1.0, 0.0, 0.0 } },
                          { "exponential", "soft-tw-3+4", { 0.0, 0.0, 1.0, 0.0 } },
                          { "exponential", "soft-tw-5+6", { 0.0, 0.0, 0.0, 1.0 } },
                          { "exponential", "soft", { 0.0, 1.0, 1.0, 1.0 } },
                          { "exponential", "all", { 1.0, 1.0, 1.0, 1.0 } },
                          // [FLvD:2022A] does not parametrise the three-particle LCDAs that Xi_1 and
                           // Xi_2 require, so the twist-3,4 switch stays off under that model and its
                           // composites deliver one term less than under the exponential model. It
                           // does define phi_-^WW, so the twist-5,6 switch behaves identically. Both
                           // facts are pinned here so neither can change unnoticed.
                           { "FLvD2022", "none", { 0.0, 0.0, 0.0, 0.0 } },
                          { "FLvD2022", "ht", { 1.0, 0.0, 0.0, 0.0 } },
                          { "FLvD2022", "soft-nlo", { 0.0, 1.0, 0.0, 0.0 } },
                          { "FLvD2022", "soft-tw-5+6", { 0.0, 0.0, 0.0, 1.0 } },
                          { "FLvD2022", "soft", { 0.0, 1.0, 0.0, 1.0 } },
                          { "FLvD2022", "all", { 1.0, 1.0, 0.0, 1.0 } },
                          // FLvD2022+profile *models* those LCDAs by the [BBJW:2018A] ansatz, so all
                           // four switches are available and its rows match the exponential model's.
                           // The identical switch pattern is exactly why the pattern alone cannot tell
                           // the two models apart; the numerical anchors above do.
                           { "FLvD2022+profile", "none", { 0.0, 0.0, 0.0, 0.0 } },
                          { "FLvD2022+profile", "ht", { 1.0, 0.0, 0.0, 0.0 } },
                          { "FLvD2022+profile", "soft-nlo", { 0.0, 1.0, 0.0, 0.0 } },
                          { "FLvD2022+profile", "soft-tw-3+4", { 0.0, 0.0, 1.0, 0.0 } },
                          { "FLvD2022+profile", "soft-tw-5+6", { 0.0, 0.0, 0.0, 1.0 } },
                          { "FLvD2022+profile", "soft", { 0.0, 1.0, 1.0, 1.0 } },
                          { "FLvD2022+profile", "all", { 1.0, 1.0, 1.0, 1.0 } } }
                })
                {
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p,
                                                                Options{
                                                                    {    "lcda-model"_ok,    qnp::OptionValue(lcda_model) },
                                                                    { "contributions"_ok, qnp::OptionValue(contributions) }
                    });
                    const Diagnostics                        d     = ff.diagnostics();
                    const unsigned                           first = d.size() - 4u;
                    for (unsigned i = 0; i < 4u; ++i)
                    {
                        TEST_CHECK_NEARLY_EQUAL(diagnostic(d, first + i), expected[i], 1e-15);
                    }
                }

                // [BBJW:2018A], plots/Deltaxi.pdf, black dash-dotted curve: Delta_xi^ht.
                for (const auto & [E_gamma, published] : std::array<std::pair<double, double>, 3>{
                         { { 1.5, +0.05162 }, { 2.0, +0.02464 }, { 2.6, +0.01145 } }
                })
                {
                    double first = 0.0;
                    for (const double & lambda_B : std::array<double, 3>{ 0.20, 0.35, 0.50 })
                    {
                        p["B_u::omega_0@FLvD2022"] = lambda_B;
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p,
                                                                    Options{
                                                                        {    "lcda-model"_ok, "exponential"_ov },
                                                                        { "contributions"_ok,          "ht"_ov }
                        });
                        const double                             delta_xi = 0.5 * (ff.F_V(E_gamma) - ff.F_A(E_gamma));
                        TEST_CHECK_NEARLY_EQUAL(delta_xi, published, 1e-4);
                        if (lambda_B == 0.20)
                        {
                            first = delta_xi;
                        }
                        else
                        {
                            TEST_CHECK_NEARLY_EQUAL(delta_xi, first, 1e-12);
                        }
                    }
                }

                // Complete twist-3--6 soft correction at E_gamma = 2 GeV. With the BBJW Table 1
                // setup used throughout this block but EOS's threshold-crossing alpha_s and GMOR
                // condensate, these values exceed the dashed curves in plots/Deltaxi.pdf by 31--36%.
                // BBJW instead use fixed-n_f=4 alpha_s(1 GeV)=0.348929 and a slightly smaller
                // condensate. The deterministic result is pinned tightly and the test fails if the
                // twist-5,6 term is removed.
                for (const auto & [lambda_B, expected] : std::array<std::pair<double, double>, 3>{
                         { { 0.20, +0.006630321674412021 }, { 0.35, +0.004272383522116230 }, { 0.50, +0.003003337577852400 } }
                })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    const Options common{
                        { "lcda-model"_ok, "exponential"_ov }
                    };
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_ht(p,
                                                                   common
                                                                           + Options{
                                                                               { "contributions"_ok, "ht"_ov }
                    });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_all(p,
                                                                    common
                                                                            + Options{
                                                                                { "contributions"_ok, "all"_ov }
                    });
                    const double                             delta_xi_soft = 0.5 * ((ff_all.F_V(2.0) - ff_all.F_A(2.0)) - (ff_ht.F_V(2.0) - ff_ht.F_A(2.0)));
                    TEST_CHECK_NEARLY_EQUAL(delta_xi_soft, expected, 1e-8);
                }
            }

            // Independent numerical reference
            {
                Parameters p                      = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"]        = 0.7;
                p["B_u::mu_0@FLvD2022"]           = 1.3;
                p["B->gamma::mu@FLvD2022QCDF"]    = 1.3;
                p["B->gamma::mu_h1@FLvD2022QCDF"] = 4.7;
                p["B->gamma::mu_h2@FLvD2022QCDF"] = 4.5;
                p["B->gamma::s_0@FLvD2022QCDF"]   = 1.59;
                p["B->gamma::M^2@FLvD2022QCDF"]   = 1.35;
                p["decay-constant::B_u"]          = 0.129;
                p["mass::B_u"]                    = 5.27929;
                p["mass::b(MSbar)"]               = 4.45432371854873; // fix m_b_pole@1-loop to 4.8
                p["mass::rho^+"]                  = 0.7;
                p["B::lambda_E^2"]                = 0.0625;
                p["B::lambda_H^2"]                = 0.125;
                p["B::LambdaBar"]                 = 0.5;
                p["B_u::a^phi+_0@FLvD2022"]       = 1.868119356054707;
                p["B_u::a^phi+_1@FLvD2022"]       = 0.151143197362311;
                p["B_u::a^phi+_2@FLvD2022"]       = 1.203196552637887;
                p["B_u::a^phi+_3@FLvD2022"]       = 0.429631987348729;
                p["B_u::a^phi+_4@FLvD2022"]       = 0.304198191688109;
                p["B_u::a^phi+_5@FLvD2022"]       = -0.324469147908141;
                p["B_u::a^phi+_6@FLvD2022"]       = 0.381019563820993;
                p["B_u::a^phi+_7@FLvD2022"]       = -0.246884872397705;
                p["B_u::a^phi+_8@FLvD2022"]       = -0.058121797086248;


                // Diagnostics: check pieces against Mathematica implementation

                AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p, Options());

                static const std::vector<std::pair<double, double>> reference = {
                    std::make_pair(3.39713985820215, 1e-9),   // L0()
                    std::make_pair(3.32067923218836, 1e-9),   // L0_incomplete(8.0)
                    std::make_pair(2.9528207810186, 1e-9),    // norm_incomplete(8.0)
                    std::make_pair(0.156908479594529, 1e-9),  // lapltr_incomplete(8.0, 4.0)
                    std::make_pair(-0.186641425933295, 1e-9), // lapltr_incomplete_dsigma(8.0, 4.0)
                    std::make_pair(0.272486748944117, 1e-3),  // L0_effective(3.0); numerical reference is imprecise
                    std::make_pair(0.0534704590581136, 1e-3), // L0_effective(2.16)
                    std::make_pair(10.4492075178413 + -10.4781709714967 + 6.58190087562423 + -8.92720937287174,
                                   1e-3), // L0_incomplete_effective(3.0, 8.0); numerical reference is imprecise
                    std::make_pair(0.101195623867872 + -0.188854545332271 + 0.334217768087141 + -0.315849218854024, 1e-6), // lapltr_effective_incomplete(3.0, 8.0, 4.0)
                    std::make_pair(0.856077632552228, 1e-8),                                                               // C at Egamma=2.16
                    std::make_pair(1.04437026215851, 1e-8),                                                                // K_inv at Egamma=2.16
                    std::make_pair(0.948778129043308, 1e-8),                                                               // U at Egamma=2.16
                    std::make_pair(0.303055287938432, 1e-8),                                                               // F_leading_power(2.16)
                    std::make_pair(-0.0261633822, 1e-5),                                                                   // xi(2.16)
                    std::make_pair(0.013380027478431 + 0.002539709, 1e-6), // delta_xi(2.16); the second term is Delta_xi^soft_(tw-5,6)
                    // Active contributions for the default options, i.e. lcda-model=FLvD2022 and
                    // contributions=all. That model defines phi_-^WW, so the twist-5,6 term is on;
                    // it does not parametrise the three-particle LCDAs that Xi_1 and Xi_2 need, so
                    // the twist-3,4 term is off. Cf. the option matrix above.
                    std::make_pair(1.0, 1e-15), // switch_ht
                    std::make_pair(1.0, 1e-15), // switch_soft
                    std::make_pair(0.0, 1e-15), // switch_soft_tw_3_4
                    std::make_pair(1.0, 1e-15), // switch_soft_tw_5_6
                };

                Diagnostics diagnostics = ff.diagnostics();
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                // Integration test: observable evaluation

                Kinematics k = Kinematics({
                    { "E_gamma", 2.16 }
                });
                Options    o{
                       { "form-factors"_ok, "FLvD2022QCDF"_ov }
                };
                auto obs_F_V = Observable::make("B->gamma::F_V(E_gamma)", p, k, o);
                auto obs_F_A = Observable::make("B->gamma::F_A(E_gamma)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL(ff.F_V(k["E_gamma"]), 0.292811641821400, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(obs_F_V->evaluate(), 0.292811641821400, 1e-8);

                TEST_CHECK_NEARLY_EQUAL(ff.F_A(k["E_gamma"]), 0.260972168867438, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(obs_F_A->evaluate(), 0.260972168867438, 1e-8);


                // Math integrity test: cross-check complete form factors against the corrected implementation

                TEST_CHECK_NEARLY_EQUAL(ff.F_V(4.0), 0.154216592881089, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.F_V(12.0), 0.046511232564673, 1e-8);

                TEST_CHECK_NEARLY_EQUAL(ff.F_A(4.0), 0.150371157254392, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(ff.F_A(12.0), 0.048711476280875, 1e-8);


                // The same generic coefficients under lcda-model=FLvD2022+profile, i.e. with the
                // modelled twist-3,4 soft term switched on. Regression values.
                //
                // These coefficients are far from the surface on which the profile-function ansatz
                // is EOM-consistent: (1/6) int domega omega^2 phi_+ = 18.3297 GeV^2 against
                // LambdaBar^2 + (2 lambda_E^2 + lambda_H^2)/6 = 0.2917 GeV^2, a ratio of 62.8. The
                // ansatz is evaluated anyway, with varkappa taken from the parameter combination,
                // which is what makes the resulting term a modelling assumption rather than a
                // controlled correction. The ratio itself is pinned in
                // heavy-meson-lcdas-flvd2022-profile_TEST.cc.
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_profile(p,
                                                                    Options{
                                                                        { "lcda-model"_ok, "FLvD2022+profile"_ov }
                });

                TEST_CHECK_NEARLY_EQUAL(ff_profile.F_V(k["E_gamma"]), 0.295761908079351, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(ff_profile.F_A(k["E_gamma"]), 0.260663359999444, 1e-8);

                // FLvD2022 itself is unchanged, and the entire difference is the twist-3,4 term.
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_profile_tw_3_4(p,
                                                                           Options{
                                                                               {    "lcda-model"_ok, "FLvD2022+profile"_ov },
                                                                               { "contributions"_ok,      "soft-tw-3+4"_ov }
                });
                AnalyticFormFactorPToGammaQCDF<BToGamma> ff_profile_none(p,
                                                                         Options{
                                                                             {    "lcda-model"_ok, "FLvD2022+profile"_ov },
                                                                             { "contributions"_ok,             "none"_ov }
                });
                TEST_CHECK_NEARLY_EQUAL(ff_profile.F_V(k["E_gamma"]) - ff.F_V(k["E_gamma"]), ff_profile_tw_3_4.F_V(k["E_gamma"]) - ff_profile_none.F_V(k["E_gamma"]), 1e-12);
                TEST_CHECK_NEARLY_EQUAL(ff_profile.F_A(k["E_gamma"]) - ff.F_A(k["E_gamma"]), ff_profile_tw_3_4.F_A(k["E_gamma"]) - ff_profile_none.F_A(k["E_gamma"]), 1e-12);
            }

            // "evolution-order" = "LL": check the truncation implemented via switch_nll
            {
                Parameters p                      = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"]        = 0.7;
                p["B_u::mu_0@FLvD2022"]           = 1.3;
                p["B->gamma::mu@FLvD2022QCDF"]    = 1.3;
                p["B->gamma::mu_h1@FLvD2022QCDF"] = 4.7;
                p["B->gamma::mu_h2@FLvD2022QCDF"] = 4.5;
                p["B->gamma::s_0@FLvD2022QCDF"]   = 1.59;
                p["B->gamma::M^2@FLvD2022QCDF"]   = 1.35;
                p["decay-constant::B_u"]          = 0.129;
                p["mass::B_u"]                    = 5.27929;
                p["mass::b(MSbar)"]               = 4.45432371854873; // fix m_b_pole@1-loop to 4.8
                p["mass::rho^+"]                  = 0.7;
                p["B::lambda_E^2"]                = 0.0625;
                p["B::lambda_H^2"]                = 0.125;
                p["B::LambdaBar"]                 = 0.5;
                p["B_u::a^phi+_0@FLvD2022"]       = 1.868119356054707;
                p["B_u::a^phi+_1@FLvD2022"]       = 0.151143197362311;
                p["B_u::a^phi+_2@FLvD2022"]       = 1.203196552637887;
                p["B_u::a^phi+_3@FLvD2022"]       = 0.429631987348729;
                p["B_u::a^phi+_4@FLvD2022"]       = 0.304198191688109;
                p["B_u::a^phi+_5@FLvD2022"]       = -0.324469147908141;
                p["B_u::a^phi+_6@FLvD2022"]       = 0.381019563820993;
                p["B_u::a^phi+_7@FLvD2022"]       = -0.246884872397705;
                p["B_u::a^phi+_8@FLvD2022"]       = -0.058121797086248;

                AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p,
                                                            Options{
                                                                { "evolution-order"_ok, "LL"_ov }
                });

                // Only C, K_inv, U, F_leading_power, and xi depend on "evolution-order"; the remaining
                // functionals reproduce the "NLL" (default) reference values from above unchanged.
                static const std::vector<std::pair<double, double>> reference = {
                    std::make_pair(3.39713985820215, 1e-9),   // L0()
                    std::make_pair(3.32067923218836, 1e-9),   // L0_incomplete(8.0)
                    std::make_pair(2.9528207810186, 1e-9),    // norm_incomplete(8.0)
                    std::make_pair(0.156908479594529, 1e-9),  // lapltr_incomplete(8.0, 4.0)
                    std::make_pair(-0.186641425933295, 1e-9), // lapltr_incomplete_dsigma(8.0, 4.0)
                    std::make_pair(0.272486748944117, 1e-3),  // L0_effective(3.0); numerical reference is imprecise
                    std::make_pair(0.0534704590581136, 1e-3), // L0_effective(2.16)
                    std::make_pair(10.4492075178413 + -10.4781709714967 + 6.58190087562423 + -8.92720937287174,
                                   1e-3), // L0_incomplete_effective(3.0, 8.0); numerical reference is imprecise
                    std::make_pair(0.101195623867872 + -0.188854545332271 + 0.334217768087141 + -0.315849218854024, 1e-6), // lapltr_effective_incomplete(3.0, 8.0, 4.0)
                    std::make_pair(1.0, 1e-14),                            // C at Egamma=2.16: no O(alpha_s(mu_h1)) matching correction at LL
                    std::make_pair(1.0, 1e-14),                            // K_inv at Egamma=2.16: K = 1 identically at LL
                    std::make_pair(0.968013353629127, 1e-8),               // U at Egamma=2.16
                    std::make_pair(0.345608964545694, 1e-8),               // F_leading_power(2.16)
                    std::make_pair(-0.0335807458, 1e-5),                   // xi(2.16)
                    std::make_pair(0.013380027478431 + 0.002539709, 1e-6), // delta_xi(2.16); the second term is Delta_xi^soft_(tw-5,6); does not depend on "evolution-order"
                    // Active contributions; "contributions" is at its default "all" here, and the
                    // twist switches are off because this block uses the default lcda-model=FLvD2022.
                    std::make_pair(1.0, 1e-15), // switch_ht
                    std::make_pair(1.0, 1e-15), // switch_soft
                    std::make_pair(0.0, 1e-15), // switch_soft_tw_3_4
                    std::make_pair(1.0, 1e-15), // switch_soft_tw_5_6
                };

                Diagnostics diagnostics = ff.diagnostics();
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                Kinematics k = Kinematics({
                    { "E_gamma", 2.16 }
                });
                Options    o{
                       {    "form-factors"_ok, "FLvD2022QCDF"_ov },
                       { "evolution-order"_ok,           "LL"_ov }
                };
                auto obs_F_V = Observable::make("B->gamma::F_V(E_gamma)", p, k, o);
                auto obs_F_A = Observable::make("B->gamma::F_A(E_gamma)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL(ff.F_V(k["E_gamma"]), 0.327947957642603, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(obs_F_V->evaluate(), 0.327947957642603, 1e-8);

                TEST_CHECK_NEARLY_EQUAL(ff.F_A(k["E_gamma"]), 0.296108484688641, 1e-8);
                TEST_CHECK_NEARLY_EQUAL(obs_F_A->evaluate(), 0.296108484688641, 1e-8);
            }
        }
} analytic_b_to_gamma_qcdf_test;
