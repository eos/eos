/*
 * Copyright (c) 2021      Méril Reboud
 * Copyright (c) 2021-2026 Danny van Dyk
 * Copyright (c) 2024      Lorenz Gärtner
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

#include <eos/observable.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/options.hh>
#include <eos/utils/units.hh>

#include <test/test.hh>

#include <format>
#include <iostream>
#include <ranges>
#include <set>

using namespace test;
using namespace eos;

class ObservableTest : public TestCase
{
    public:
        ObservableTest() :
            TestCase("observable_test")
        {
        }

        virtual void
        run() const
        {
            /* Test insertion of a new observable */
            {
                auto observables = Observables();

                TEST_CHECK_THROWS(ParsingError, observables.insert("mass::ratio", "m_r", Unit::None(), Options(), "<<mass::mu>> /* <<mass::tau>>"));

                TEST_CHECK_THROWS(UnknownObservableError, observables.insert("mass::ratio", "m_r", Unit::None(), Options(), "<<mass::qwerty>>"));

                observables.insert("mass::ratio", "m_r", Unit::None(), Options(), "<<mass::mu>> / <<mass::tau>>");

                Parameters p   = Parameters::Defaults();
                p["mass::mu"]  = 0.105658;
                p["mass::tau"] = 1.77682;
                Kinematics k;
                Options    o;

                TEST_CHECK(observables["mass::ratio"]->make(p, k, o));
                TEST_CHECK_NEARLY_EQUAL(observables["mass::ratio"]->make(p, k, o)->evaluate(), 0.059464662, 10e-5);

                auto observable = Observable::make("mass::ratio", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(observable->evaluate(), 0.059464662, 10e-5);
            }

            /* Test insertion and evaluation of a cacheable observable */
            {
                auto observables = Observables();

                observables.insert("B->D^*lnu::2*S_1c",
                                   R"()",
                                   Unit::Undefined(),
                                   Options(),
                                   R"(
                                2 * <<B->D^*lnu::S_1c;l=mu>>
                                )");

                Parameters p = Parameters::Defaults();
                Kinematics k{
                    { "q2_min",  4.00 },
                    { "q2_max", 10.68 }
                };
                Options o;

                TEST_CHECK(observables["B->D^*lnu::2*S_1c"]->make(p, k, o));
                TEST_CHECK_NO_THROW(observables["B->D^*lnu::2*S_1c"]->make(p, k, o)->evaluate());

                auto observable = Observable::make("B->D^*lnu::2*S_1c", p, k, o);
                TEST_CHECK_NO_THROW(observable->evaluate());

                o.declare("l"_ok, "mu"_ov);
                auto obs_S1c = Observable::make("B->D^*lnu::S_1c", p, k, o);
                TEST_CHECK_EQUAL(observable->evaluate(), 2 * obs_S1c->evaluate());
            }

            /* Test insertion of nested cacheable observables */
            {
                auto observables = Observables();

                observables.insert("B->D^*lnu::4*S_1c",
                                   R"()",
                                   Unit::Undefined(),
                                   Options(),
                                   R"(
                                2 * <<B->D^*lnu::2*S_1c;l=mu>>
                                )");

                Parameters p = Parameters::Defaults();
                Kinematics k{
                    { "q2_min",  4.00 },
                    { "q2_max", 10.68 }
                };
                Options o;

                TEST_CHECK(observables["B->D^*lnu::4*S_1c"]->make(p, k, o));
                TEST_CHECK_NO_THROW(observables["B->D^*lnu::4*S_1c"]->make(p, k, o)->evaluate());

                auto observable = Observable::make("B->D^*lnu::4*S_1c", p, k, o);
                TEST_CHECK_NO_THROW(observable->evaluate());

                auto obs_2S1c = Observable::make("B->D^*lnu::2*S_1c", p, k, o);
                TEST_CHECK_EQUAL(observable->evaluate(), 2 * obs_2S1c->evaluate());
            }

            /* Test insertion of nested multiple cacheable observables with different kinematics */
            {
                auto observables = Observables();
                std::printf("\n\n Starting insertion of nested \n");
                observables.insert("B->D^*lnu::S_1c_q2_min_num/S_1c_q2_min_denom",
                                   R"()",
                                   Unit::Undefined(),
                                   Options(),
                                   R"(
                                <<B->D^*lnu::S_1c;l=mu>>[q2_min=>q2_min_num] / <<B->D^*lnu::S_1c;l=mu>>[q2_min=>q2_min_denom]
                                )");

                Parameters p = Parameters::Defaults();
                Kinematics k{
                    {   "q2_min_num",  4.00 },
                    { "q2_min_denom",  8.00 },
                    {       "q2_max", 10.68 }
                };
                Kinematics k_num{
                    { "q2_min",  4.00 },
                    { "q2_max", 10.68 }
                };
                Kinematics k_denom{
                    { "q2_min",  8.00 },
                    { "q2_max", 10.68 }
                };
                Options o;

                auto obs_num   = Observable::make("B->D^*lnu::S_1c", p, k_num, o);
                auto obs_denom = Observable::make("B->D^*lnu::S_1c", p, k_denom, o);
                auto obs       = Observable::make("B->D^*lnu::S_1c_q2_min_num/S_1c_q2_min_denom", p, k, o);

                ObservableCache cache(p);
                TEST_CHECK_NO_THROW(cache.add(obs_num));
                TEST_CHECK_NO_THROW(cache.add(obs_denom));

                TEST_CHECK_RELATIVE_ERROR(obs->evaluate(), obs_num->evaluate() / obs_denom->evaluate(), 1e-3);

                ObservableCache::ObservableId obs_id;
                TEST_CHECK_NO_THROW(obs_id = cache.add(obs));

                TEST_CHECK_NO_THROW(cache.update());
                TEST_CHECK_RELATIVE_ERROR(obs->evaluate(), cache[obs_id], 1e-3);
            }

            /* Test a fixed kinematic value of a variable that the kinematics also provide */
            {
                auto observables = Observables();

                // P(z) = z (4 - z), so P(1) = 3 and P(2) = 4
                observables.insert("test::P(z)/P(1)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>> / <<TestLegendre1D::UnnormalizedPDF(z)>>[z=1.0]");

                Parameters p = Parameters::Defaults();
                Kinematics k{
                    { "z", 2.0 }
                };
                Options o;

                auto obs = Observable::make("test::P(z)/P(1)", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(obs->evaluate(), 4.0 / 3.0, 1e-12);

                k.set("z", 1.0);
                TEST_CHECK_NEARLY_EQUAL(obs->evaluate(), 1.0, 1e-12);

                // an outer reference must not rename a variable that the inner expression fixes already
                observables.insert("test::P(z=3)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>>[z=3.0]");
                observables.insert("test::P(z=3)[z=w]", R"()", Unit::None(), Options(), "<<test::P(z=3)>>[z=>w]");

                Kinematics k_w{
                    { "w", 1.0 }
                };
                TEST_CHECK_THROWS(DuplicateKinematicAliasError, Observable::make("test::P(z=3)[z=w]", p, k_w, o));
                TEST_CHECK(! k_w.has("z"));
            }

            /* Test fixed kinematic values of several, also nested, expressions that share one set of kinematics */
            {
                auto observables = Observables();

                // P(z) = z (4 - z), so P(1) = 3 and P(2) = 4
                observables.insert("test::P(z=1)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>>[z=1.0]");
                observables.insert("test::P(z=2)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>>[z=2.0]");
                observables.insert("test::P(z=1)*P(z=2)", R"()", Unit::None(), Options(), "<<test::P(z=1)>> * <<TestLegendre1D::UnnormalizedPDF(z)>>[z=2.0]");

                Parameters p = Parameters::Defaults();
                Kinematics k;
                Options    o;

                auto obs_1       = Observable::make("test::P(z=1)", p, k, o);
                auto obs_2       = Observable::make("test::P(z=2)", p, k, o);
                auto obs_product = Observable::make("test::P(z=1)*P(z=2)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL(obs_1->evaluate(), 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(obs_2->evaluate(), 4.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(obs_product->evaluate(), 12.0, 1e-12);

                TEST_CHECK_NEARLY_EQUAL(obs_1->clone()->evaluate(), 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(obs_product->clone()->evaluate(), 12.0, 1e-12);

                ObservableCache cache(p);
                auto            id_1       = cache.add(obs_1);
                auto            id_2       = cache.add(obs_2);
                auto            id_product = cache.add(obs_product);

                TEST_CHECK_NO_THROW(cache.update());
                TEST_CHECK_NEARLY_EQUAL(cache[id_1], 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache[id_2], 4.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache[id_product], 12.0, 1e-12);

                ObservableCache cache_clone = cache.clone(p.clone());
                TEST_CHECK_NO_THROW(cache_clone.update());
                TEST_CHECK_NEARLY_EQUAL(cache_clone[id_1], 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache_clone[id_2], 4.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache_clone[id_product], 12.0, 1e-12);
            }

            /* Test kinematic substitutions */
            {
                auto observables = Observables();

                // P(z) = z (4 - z), so P(1) = 3, P(2) = 4, and P(3) = 3
                observables.insert("test::P(z=2w)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>>[z=2 * {w}]");
                observables.insert("test::P(z=2E)", R"()", Unit::None(), Options(), "<<test::P(z=2w)>>[w=>E]");
                observables.insert("test::P(z=2z)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>>[z=2 * {z}]");

                const auto kinematic_variables = [&observables](const char * name)
                {
                    const auto entry = observables[name];
                    return std::set<std::string>(entry->begin_kinematic_variables(), entry->end_kinematic_variables());
                };
                TEST_CHECK_EQUAL(kinematic_variables("test::P(z=2w)"), std::set<std::string>{ "w" });
                TEST_CHECK_EQUAL(kinematic_variables("test::P(z=2E)"), std::set<std::string>{ "E" });
                TEST_CHECK_EQUAL(kinematic_variables("test::P(z=2z)"), std::set<std::string>{ "z" });

                TEST_CHECK_THROWS(exp::ExpressionError,
                                  observables.insert("test::P(z=c)", R"()", Unit::None(), Options(), "<<TestLegendre1D::UnnormalizedPDF(z)>>[z=[[mass::c]] * {w}]"));

                Parameters p = Parameters::Defaults();
                Kinematics k{
                    { "E", 1.5 },
                    { "w", 1.0 },
                    { "z", 0.5 },
                };
                Options o;

                auto obs_w = Observable::make("test::P(z=2w)", p, k, o);
                auto obs_E = Observable::make("test::P(z=2E)", p, k, o);
                auto obs_z = Observable::make("test::P(z=2z)", p, k, o);

                TEST_CHECK_NEARLY_EQUAL(obs_w->evaluate(), 4.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(obs_E->evaluate(), 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(obs_z->evaluate(), 3.0, 1e-12);

                // the observables report the variables that their substitutions read
                const auto used_kinematics = [](const ObservablePtr & observable)
                {
                    const KinematicUser & user = static_cast<const KinematicUser &>(*observable);
                    return std::set<KinematicVariable::Id>(user.begin_kinematics(), user.end_kinematics());
                };
                TEST_CHECK(used_kinematics(obs_w).contains(k["w"].id()));
                TEST_CHECK(used_kinematics(obs_E).contains(k["E"].id()));
                TEST_CHECK(used_kinematics(obs_z).contains(k["z"].id()));

                k.set("w", 0.5);
                TEST_CHECK_NEARLY_EQUAL(obs_w->evaluate(), 3.0, 1e-12);

                // a clone follows its own kinematics
                auto clone = obs_w->clone();
                clone->kinematics().set("w", 1.5);
                TEST_CHECK_NEARLY_EQUAL(clone->evaluate(), 3.0, 1e-12);
                k.set("w", 1.0);
                TEST_CHECK_NEARLY_EQUAL(obs_w->evaluate(), 4.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(clone->evaluate(), 3.0, 1e-12);

                // the cache reads the substituted variable when it is updated
                ObservableCache cache(p);
                auto            id_w = cache.add(obs_w);
                auto            id_E = cache.add(obs_E);

                TEST_CHECK_NO_THROW(cache.update());
                TEST_CHECK_NEARLY_EQUAL(cache[id_w], 4.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache[id_E], 3.0, 1e-12);

                k.set("w", 0.5);
                k.set("E", 1.0);
                TEST_CHECK_NO_THROW(cache.update());
                TEST_CHECK_NEARLY_EQUAL(cache[id_w], 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache[id_E], 4.0, 1e-12);

                ObservableCache cache_clone = cache.clone(p.clone());
                TEST_CHECK_NO_THROW(cache_clone.update());
                TEST_CHECK_NEARLY_EQUAL(cache_clone[id_w], 3.0, 1e-12);
                TEST_CHECK_NEARLY_EQUAL(cache_clone[id_E], 4.0, 1e-12);
            }

            /* Test insertion of a problematic observable */
            {
                auto observables = Observables();

                TEST_CHECK_THROWS(InternalError,
                                  observables.insert("test::problematic_observable",
                                                     R"()",
                                                     Unit::Undefined(),
                                                     Options(),
                                                     R"(
                        <<test::obs1>>[q2_min=>q2_min_num] * {q2_min}
                    )"));
            }

            /* Test the names of the kinematic variables against a whitelist */
            {
                auto return_allowed_kinematic_variables = []() -> std::set<std::string>
                {
                    std::set<std::string> result{ "q2",
                                                  "q2_min",
                                                  "q2_max",
                                                  "q2_num",
                                                  "q2_min_num",
                                                  "q2_max_num",
                                                  "q2_denom",
                                                  "q2_min_denom",
                                                  "q2_max_denom",
                                                  "Re{q2}",
                                                  "Im{q2}",
                                                  "k2",
                                                  "k2_min",
                                                  "k2_max",
                                                  "sqrt(k2)",
                                                  "sqrt(k2)_min",
                                                  "sqrt(k2)_max",
                                                  "k2_min_num",
                                                  "k2_max_num",
                                                  "k2_min_denom",
                                                  "k2_max_denom",
                                                  "kperp",
                                                  "kperp_min",
                                                  "kperp_max",
                                                  "phi",
                                                  "phi_min",
                                                  "phi_max",
                                                  "Re{psi}",
                                                  "Im{psi}",
                                                  "E",
                                                  "E_min",
                                                  "E_max",
                                                  "Re{E}",
                                                  "Im{E}",
                                                  "E_gamma",
                                                  "E_gamma_min",
                                                  "E_pi",
                                                  "w",
                                                  "w_min",
                                                  "w_max",
                                                  "w_min_num",
                                                  "w_max_num",
                                                  "w_min_denom",
                                                  "w_max_denom",
                                                  "mu",
                                                  "tau",
                                                  "z",
                                                  "z_min",
                                                  "z_max",
                                                  "z_min_num",
                                                  "z_max_num",
                                                  "z_min_denom",
                                                  "z_max_denom",
                                                  "z_gamma",
                                                  "z_w",
                                                  // Temporary hack for discrete itemization
                                                  "k" };

                    for (std::string l : { "e", "mu", "tau" })
                    {
                        result.insert(std::format("q2_{}", l));
                        result.insert(std::format("q2_{}_min", l));
                        result.insert(std::format("q2_{}_max", l));
                    }

                    for (std::string M : { "l", "pi", "K", "D", "D_s", "L" })
                    {
                        result.insert(std::format("cos(theta_{})", M));
                        result.insert(std::format("cos(theta_{})_min", M));
                        result.insert(std::format("cos(theta_{})_max", M));
                    }

                    return result;
                };

                static const auto allowed_kinematic_variables = return_allowed_kinematic_variables();

                auto observables = Observables();

                bool found_problematic_observable = false;

                for (const auto & [name, entry] : observables)
                {
                    // The whitelist governs the vocabulary that EOS's physics observables expose to
                    // users. The ``Test...`` fixtures that the test runner registers in every test
                    // executable are not part of it and sample abstract variables such as 'x' and 'y'.
                    if (name.prefix_part().str().starts_with("Test"))
                    {
                        continue;
                    }

                    for (const auto & kinematic_variable : std::ranges::subrange(entry->begin_kinematic_variables(), entry->end_kinematic_variables()))
                    {
                        if (allowed_kinematic_variables.find(kinematic_variable) != allowed_kinematic_variables.end())
                        {
                            continue;
                        }

                        found_problematic_observable = true;
                        std::cerr << "Found observable: " << name << ", with invalid kinematic variable: " << kinematic_variable << std::endl;
                    }
                }

                TEST_CHECK(found_problematic_observable == false);
            }

            /* Test the LaTeX description of observables for unneeded dollar signs */
            {
                auto observables = Observables();

                bool found_problematic_observable = false;

                for (const auto & [name, entry] : observables)
                {
                    if (entry->latex().find("$") == std::string::npos)
                    {
                        continue;
                    }

                    found_problematic_observable = true;
                    std::cerr << "Found problematic observable: " << name << ", with latex description: " << entry->latex() << "that has (likely) unneeded dollar signs"
                              << std::endl;
                }

                TEST_CHECK(found_problematic_observable == false);
            }

            // Observables::has
            {
                Observables o = Observables();

                TEST_CHECK_EQUAL(o.has("B_c->lnu::BR"), true);
                TEST_CHECK_EQUAL(o.has("B_c->lnu::TEST"), false);

                // an observable inserted at run time must be found, just as operator[] finds it
                o.insert("mass::ratio@has", "m_r", Unit::None(), Options(), "<<mass::mu>> / <<mass::tau>>");

                TEST_CHECK_EQUAL(o.has("mass::ratio@has"), true);
                TEST_CHECK_NO_THROW(o["mass::ratio@has"]);
            }
        }

} observable_test;
