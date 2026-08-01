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

#include <test/test.hh>
#include <eos/form-factors/analytic-p-to-gamma-qcdf.hh>
#include <eos/models/model.hh>
#include <eos/observable.hh>
#include <eos/utils/diagnostics.hh>

#include <array>
#include <cmath>

using namespace test;
using namespace eos;

class AnalyticFormFactorBToGammaQCDFTest :
    public TestCase
{
    public:
        AnalyticFormFactorBToGammaQCDFTest():
            TestCase("analytic_b_to_gamma_qcdf_test")
        {
        }

        virtual void run() const
        {
            // Regression tests at the inputs of [BBJW:2018A], Table 1, with sigmahat_1 = 0.
            {
                Parameters p = Parameters::Defaults();
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
                const auto diagnostic = [] (const Diagnostics & d, const unsigned i)
                {
                    auto entry = d.begin();
                    for (unsigned j = 0 ; j < i ; ++j)
                        ++entry;
                    return entry->value;
                };

                const auto average = [] (const AnalyticFormFactorPToGammaQCDF<BToGamma> & ff, const double & E_gamma)
                {
                    return 0.5 * (ff.F_V(E_gamma) + ff.F_A(E_gamma));
                };

                // [BR:2011A], Eq. (2.14): J = L_hat^2 - 1 for the exponential LCDA.
                for (const double & lambda_B : std::array<double, 2>{ 0.20, 0.35 })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p, Options{ { "lcda-model"_ok, "exponential"_ov } });
                    const Diagnostics d = ff.diagnostics();

                    for (const auto & [E_gamma, index] : std::array<std::pair<double, unsigned>, 2>{ { { 3.0, 5 }, { 2.16, 6 } } })
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
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p, Options{ { "lcda-model"_ok, "exponential"_ov } });
                    const Diagnostics d = ff.diagnostics();
                    const double value = diagnostic(d, 9) * diagnostic(d, 10) * diagnostic(d, 11);
                    if (lambda_B == 0.20)
                        perturbative_factor = value;
                    else
                        TEST_CHECK_NEARLY_EQUAL(value, perturbative_factor, 1e-12);
                }

                // [BBJW:2018A], Eq. (5.11) and Fig. 4: xi-hat^soft_(LO) at E_gamma = 2 GeV.
                for (const auto & [lambda_B, published] : std::array<std::pair<double, double>, 3>{ {
                        { 0.20, -0.7875 }, { 0.35, -0.4085 }, { 0.50, -0.2625 }
                    } })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    const Options common{ { "lcda-model"_ok, "exponential"_ov }, { "evolution-order"_ok, "LL"_ov } };
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_none(p, common + Options{ { "contributions"_ok, "none"_ov } });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_soft(p, common + Options{ { "contributions"_ok, "soft"_ov } });
                    const double F_LP = average(ff_none, 2.0);
                    const double xi_soft = average(ff_soft, 2.0) - F_LP;
                    TEST_CHECK_NEARLY_EQUAL(2.0 * 2.0 * xi_soft / F_LP, published, 1e-3);
                }

                // [BBJW:2018A], plots/sigma0.pdf, dash-dotted curve: xi^ht at E_gamma = 2 GeV.
                for (const auto & [lambda_B, published] : std::array<std::pair<double, double>, 3>{ {
                        { 0.20, -0.07039 }, { 0.35, -0.08846 }, { 0.50, -0.09570 }
                    } })
                {
                    p["B_u::omega_0@FLvD2022"] = lambda_B;
                    const Options common{ { "lcda-model"_ok, "exponential"_ov } };
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_none(p, common + Options{ { "contributions"_ok, "none"_ov } });
                    AnalyticFormFactorPToGammaQCDF<BToGamma> ff_ht(p, common + Options{ { "contributions"_ok, "ht"_ov } });
                    TEST_CHECK_NEARLY_EQUAL(average(ff_ht, 2.0) - average(ff_none, 2.0), published, 1e-4);
                }

                // [BBJW:2018A], plots/Deltaxi.pdf, black dash-dotted curve: Delta_xi^ht.
                for (const auto & [E_gamma, published] : std::array<std::pair<double, double>, 3>{ {
                        { 1.5, +0.05162 }, { 2.0, +0.02464 }, { 2.6, +0.01145 }
                    } })
                {
                    double first = 0.0;
                    for (const double & lambda_B : std::array<double, 3>{ 0.20, 0.35, 0.50 })
                    {
                        p["B_u::omega_0@FLvD2022"] = lambda_B;
                        AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p, Options{
                                { "lcda-model"_ok, "exponential"_ov }, { "contributions"_ok, "ht"_ov }
                            });
                        const double delta_xi = 0.5 * (ff.F_V(E_gamma) - ff.F_A(E_gamma));
                        TEST_CHECK_NEARLY_EQUAL(delta_xi, published, 1e-4);
                        if (lambda_B == 0.20)
                            first = delta_xi;
                        else
                            TEST_CHECK_NEARLY_EQUAL(delta_xi, first, 1e-12);
                    }
                }
            }

            // Independent numerical reference
            {
                Parameters p = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"]  = 0.7;
                p["B_u::mu_0@FLvD2022"]  = 1.3;
                p["B->gamma::mu@FLvD2022QCDF"]  = 1.3;
                p["B->gamma::mu_h1@FLvD2022QCDF"]  = 4.7;
                p["B->gamma::mu_h2@FLvD2022QCDF"]  = 4.5;
                p["B->gamma::s_0@FLvD2022QCDF"] = 1.59;
                p["B->gamma::M^2@FLvD2022QCDF"] = 1.35;
                p["decay-constant::B_u"] = 0.129;
                p["mass::B_u"] = 5.27929;
                p["mass::b(MSbar)"] = 4.45432371854873; // fix m_b_pole@1-loop to 4.8
                p["mass::rho^+"] = 0.7;
                p["B::lambda_E^2"] = 0.0625;
                p["B::lambda_H^2"] = 0.125;
                p["B::LambdaBar"] = 0.5;
                p["B_u::a^phi+_0@FLvD2022"] =  1.868119356054707;
                p["B_u::a^phi+_1@FLvD2022"] =  0.151143197362311;
                p["B_u::a^phi+_2@FLvD2022"] =  1.203196552637887;
                p["B_u::a^phi+_3@FLvD2022"] =  0.429631987348729;
                p["B_u::a^phi+_4@FLvD2022"] =  0.304198191688109;
                p["B_u::a^phi+_5@FLvD2022"] = -0.324469147908141;
                p["B_u::a^phi+_6@FLvD2022"] =  0.381019563820993;
                p["B_u::a^phi+_7@FLvD2022"] = -0.246884872397705;
                p["B_u::a^phi+_8@FLvD2022"] = -0.058121797086248;


                // Diagnostics: check pieces against Mathematica implementation

                AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p, Options());

                static const std::vector<std::pair<double, double>> reference = {
                    std::make_pair(3.39713985820215,    1e-9), // L0()
                    std::make_pair(3.32067923218836,    1e-9), // L0_incomplete(8.0)
                    std::make_pair(2.9528207810186,     1e-9), // norm_incomplete(8.0)
                    std::make_pair(0.156908479594529,   1e-9), // lapltr_incomplete(8.0, 4.0)
                    std::make_pair(-0.186641425933295,  1e-9), // lapltr_incomplete_dsigma(8.0, 4.0)
                    std::make_pair(0.272486748944117,   1e-3), // L0_effective(3.0); numerical reference is imprecise
                    std::make_pair(0.0534704590581136,  1e-3), // L0_effective(2.16)
                    std::make_pair(10.4492075178413 + -10.4781709714967 + 6.58190087562423 + -8.92720937287174,     1e-3), // L0_incomplete_effective(3.0, 8.0); numerical reference is imprecise
                    std::make_pair(0.101195623867872 + -0.188854545332271 + 0.334217768087141 + -0.315849218854024, 1e-6), // lapltr_effective_incomplete(3.0, 8.0, 4.0)
                    std::make_pair(0.856077632552228,  1e-8), // C at Egamma=2.16
                    std::make_pair(1.04437026215851,   1e-8), // K_inv at Egamma=2.16
                    std::make_pair(0.948778129043308,  1e-8), // U at Egamma=2.16
                    std::make_pair(0.303055287938432,  1e-8), // F_leading_power(2.16)
                    std::make_pair(-0.0419868760541205, 1e-5), // xi(2.16)
                    std::make_pair(0.013380027478431 + 0.0, 1e-6), // delta_xi(2.16)
                };

                Diagnostics diagnostics = ff.diagnostics();
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                // Integration test: observable evaluation

                Kinematics k = Kinematics({ { "E_gamma", 2.16 } });
                Options o { {"form-factors"_ok, "FLvD2022QCDF"_ov} };
                auto obs_F_V = Observable::make("B->gamma::F_V(E_gamma)", p, k, o);
                auto obs_F_A = Observable::make("B->gamma::F_A(E_gamma)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL( ff.F_V(k["E_gamma"]),      0.274448439362743,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( obs_F_V->evaluate(),       0.274448439362743,     1e-8 );

                TEST_CHECK_NEARLY_EQUAL( ff.F_A(k["E_gamma"]),      0.247688384405881,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( obs_F_A->evaluate(),       0.247688384405881,     1e-8 );


                // Math integrity test: cross-check complete form factors against the corrected implementation

                TEST_CHECK_NEARLY_EQUAL( ff.F_V(4.0),               0.14854710850733,      1e-8 );
                TEST_CHECK_NEARLY_EQUAL( ff.F_V(12.0),              0.0458074023590614,    1e-8 );

                TEST_CHECK_NEARLY_EQUAL( ff.F_A(4.0),               0.146182831168587,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( ff.F_A(12.0),              0.0481722192183694,    1e-8 );
            }

            // "evolution-order" = "LL": check the truncation implemented via switch_nll
            {
                Parameters p = Parameters::Defaults();
                p["B_u::omega_0@FLvD2022"]  = 0.7;
                p["B_u::mu_0@FLvD2022"]  = 1.3;
                p["B->gamma::mu@FLvD2022QCDF"]  = 1.3;
                p["B->gamma::mu_h1@FLvD2022QCDF"]  = 4.7;
                p["B->gamma::mu_h2@FLvD2022QCDF"]  = 4.5;
                p["B->gamma::s_0@FLvD2022QCDF"] = 1.59;
                p["B->gamma::M^2@FLvD2022QCDF"] = 1.35;
                p["decay-constant::B_u"] = 0.129;
                p["mass::B_u"] = 5.27929;
                p["mass::b(MSbar)"] = 4.45432371854873; // fix m_b_pole@1-loop to 4.8
                p["mass::rho^+"] = 0.7;
                p["B::lambda_E^2"] = 0.0625;
                p["B::lambda_H^2"] = 0.125;
                p["B::LambdaBar"] = 0.5;
                p["B_u::a^phi+_0@FLvD2022"] =  1.868119356054707;
                p["B_u::a^phi+_1@FLvD2022"] =  0.151143197362311;
                p["B_u::a^phi+_2@FLvD2022"] =  1.203196552637887;
                p["B_u::a^phi+_3@FLvD2022"] =  0.429631987348729;
                p["B_u::a^phi+_4@FLvD2022"] =  0.304198191688109;
                p["B_u::a^phi+_5@FLvD2022"] = -0.324469147908141;
                p["B_u::a^phi+_6@FLvD2022"] =  0.381019563820993;
                p["B_u::a^phi+_7@FLvD2022"] = -0.246884872397705;
                p["B_u::a^phi+_8@FLvD2022"] = -0.058121797086248;

                AnalyticFormFactorPToGammaQCDF<BToGamma> ff(p, Options{ { "evolution-order"_ok, "LL"_ov } });

                // Only C, K_inv, U, F_leading_power, and xi depend on "evolution-order"; the remaining
                // functionals reproduce the "NLL" (default) reference values from above unchanged.
                static const std::vector<std::pair<double, double>> reference = {
                    std::make_pair(3.39713985820215,    1e-9), // L0()
                    std::make_pair(3.32067923218836,    1e-9), // L0_incomplete(8.0)
                    std::make_pair(2.9528207810186,     1e-9), // norm_incomplete(8.0)
                    std::make_pair(0.156908479594529,   1e-9), // lapltr_incomplete(8.0, 4.0)
                    std::make_pair(-0.186641425933295,  1e-9), // lapltr_incomplete_dsigma(8.0, 4.0)
                    std::make_pair(0.272486748944117,   1e-3), // L0_effective(3.0); numerical reference is imprecise
                    std::make_pair(0.0534704590581136,  1e-3), // L0_effective(2.16)
                    std::make_pair(10.4492075178413 + -10.4781709714967 + 6.58190087562423 + -8.92720937287174,     1e-3), // L0_incomplete_effective(3.0, 8.0); numerical reference is imprecise
                    std::make_pair(0.101195623867872 + -0.188854545332271 + 0.334217768087141 + -0.315849218854024, 1e-6), // lapltr_effective_incomplete(3.0, 8.0, 4.0)
                    std::make_pair(1.0,                1e-14), // C at Egamma=2.16: no O(alpha_s(mu_h1)) matching correction at LL
                    std::make_pair(1.0,                1e-14), // K_inv at Egamma=2.16: K = 1 identically at LL
                    std::make_pair(0.968013353629127,  1e-8), // U at Egamma=2.16
                    std::make_pair(0.345608964545694,  1e-8), // F_leading_power(2.16)
                    std::make_pair(-0.0494042368401784, 1e-5), // xi(2.16)
                    std::make_pair(0.013380027478431 + 0.0, 1e-6), // delta_xi(2.16); does not depend on "evolution-order"
                };

                Diagnostics diagnostics = ff.diagnostics();
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                Kinematics k = Kinematics({ { "E_gamma", 2.16 } });
                Options o { {"form-factors"_ok, "FLvD2022QCDF"_ov}, {"evolution-order"_ok, "LL"_ov} };
                auto obs_F_V = Observable::make("B->gamma::F_V(E_gamma)", p, k, o);
                auto obs_F_A = Observable::make("B->gamma::F_A(E_gamma)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL( ff.F_V(k["E_gamma"]),      0.309584755183946,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( obs_F_V->evaluate(),       0.309584755183946,     1e-8 );

                TEST_CHECK_NEARLY_EQUAL( ff.F_A(k["E_gamma"]),      0.282824700227084,     1e-8 );
                TEST_CHECK_NEARLY_EQUAL( obs_F_A->evaluate(),       0.282824700227084,     1e-8 );
            }
        }
} analytic_b_to_gamma_qcdf_test;
