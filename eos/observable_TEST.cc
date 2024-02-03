/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2021-2023 Danny van Dyk
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
#include <eos/observable.hh>
#include <eos/utils/options.hh>
#include <eos/utils/units.hh>

#include <ranges>
#include <iostream>

using namespace test;
using namespace eos;

class ObservableTest :
    public TestCase
{
    public:
        ObservableTest() :
            TestCase("observable_test")
        {
        }

        virtual void run() const
        {
            /* Test insertion of a new observable */
            {
                auto observables = Observables();

                TEST_CHECK_THROWS(
                    ParsingError,
                    observables.insert("mass::ratio", "m_r", Unit::None(), Options(), "<<mass::mu>> /* <<mass::tau>>")
                );

                TEST_CHECK_THROWS(
                    UnknownObservableError,
                    observables.insert("mass::ratio", "m_r", Unit::None(), Options(), "<<mass::qwerty>>")
                );

                observables.insert("mass::ratio", "m_r", Unit::None(), Options(), "<<mass::mu>> / <<mass::tau>>");

                Parameters p = Parameters::Defaults();
                p["mass::mu"]     =  0.105658;
                p["mass::tau"]    =  1.77682;
                Kinematics k;
                Options o;

                TEST_CHECK(observables["mass::ratio"]->make(p, k, o));
                TEST_CHECK_NEARLY_EQUAL(observables["mass::ratio"]->make(p, k, o)->evaluate(), 0.059464662, 10e-5);

                auto observable = Observable::make("mass::ratio", p, k, o);
                TEST_CHECK_NEARLY_EQUAL(observable->evaluate(), 0.059464662, 10e-5);
            }

            /* Test insertion and evaluation of a cacheable observable */
            {
                auto observables = Observables();

                observables.insert("B->D^*lnu::2*S_1c", R"()", Unit::Undefined(), Options(),
                                R"(
                                2 * <<B->D^*lnu::S_1c;l=mu>>
                                )");

                Parameters p = Parameters::Defaults();
                Kinematics k
                {
                    { "q2_min",   4.00 }, { "q2_max",  10.68 }
                };
                Options o;

                TEST_CHECK(observables["B->D^*lnu::2*S_1c"]->make(p, k, o));
                TEST_CHECK_NO_THROW(observables["B->D^*lnu::2*S_1c"]->make(p, k, o)->evaluate());

                auto observable = Observable::make("B->D^*lnu::2*S_1c", p, k, o);
                TEST_CHECK_NO_THROW(observable->evaluate());

                o.declare("l", "mu");
                auto obs_S1c = Observable::make("B->D^*lnu::S_1c", p, k, o);
                TEST_CHECK_EQUAL(observable->evaluate(), 2 * obs_S1c->evaluate());
            }

            /* Test insertion of nested cacheable observables */
            {
                auto observables = Observables();

                observables.insert("B->D^*lnu::4*S_1c", R"()", Unit::Undefined(), Options(),
                                R"(
                                2 * <<B->D^*lnu::2*S_1c;l=mu>>
                                )");

                Parameters p = Parameters::Defaults();
                Kinematics k
                {
                    { "q2_min",   4.00 }, { "q2_max",  10.68 }
                };
                Options o;

                TEST_CHECK(observables["B->D^*lnu::4*S_1c"]->make(p, k, o));
                TEST_CHECK_NO_THROW(observables["B->D^*lnu::4*S_1c"]->make(p, k, o)->evaluate());

                auto observable = Observable::make("B->D^*lnu::4*S_1c", p, k, o);
                TEST_CHECK_NO_THROW(observable->evaluate());

                auto obs_2S1c = Observable::make("B->D^*lnu::2*S_1c", p, k, o);
                TEST_CHECK_EQUAL(observable->evaluate(), 2 * obs_2S1c->evaluate());
            }

            /* Test inserion of a problematic observable */
            {
                auto observables = Observables();

                TEST_CHECK_THROWS(InternalError,
                    observables.insert("test::problematic_observable", R"()", Unit::Undefined(), Options(),
                    R"(
                        <<test::obs1>>[q2_min=>q2_min_num] * {q2_min}
                    )")
                );
            }

            /* Test the names of the kinematic variables against a whitelist */
            {
                static const std::set<std::string> allowed_kinematic_variables
                {
                    "q2", "q2_min", "q2_max",
                    "q2_num", "q2_min_num", "q2_max_num",
                    "q2_denom", "q2_min_denom", "q2_max_denom",
                    "Re{q2}", "Im{q2}",
                    "q2_e", "q2_e_min", "q2_e_max",
                    "q2_mu", "q2_mu_min", "q2_mu_max",
                    "q2_tau", "q2_tau_min", "q2_tau_max",
                    "k2", "k2_min", "k2_max", "k2_min_num", "k2_max_num", "k2_min_denom", "k2_max_denom",
                    "cos(theta_l)",
                    "cos(theta_pi)",
                    "cos(theta_k)",
                    "phi",
                    "E", "E_min", "E_max",
                    "Re{E}", "Im{E}",
                    "E_gamma", "E_gamma_min",
                    "E_pi",
                    "w", "w_min", "w_max", "w_min_num", "w_max_num", "w_min_denom", "w_max_denom",
                    "mu", "tau",
                    "z", "z_min", "z_max", "z_min_num", "z_max_num", "z_min_denom", "z_max_denom",
                    // needs to be unified with the notation of Re{E} and Im{E}
                    "q2_real", "q2_imag",
                    // deprecated
                    "s", "s_min", "s_max",
                    "i",
                    "c_D", "c_D_min", "c_D_max",
                    "c_l", "c_l_min", "c_l_max",
                    "chi", "chi_min", "chi_max",
                };

                auto observables = Observables();

                bool found_problematic_observable = false;

                for (auto const & [name, entry] : observables)
                {
                    for (auto const & kinematic_variable : std::ranges::subrange(entry->begin_kinematic_variables(), entry->end_kinematic_variables()))
                    {
                        if (allowed_kinematic_variables.find(kinematic_variable) != allowed_kinematic_variables.end())
                            continue;

                        found_problematic_observable = true;
                        std::cerr << "Found problematic observable: " << name << ", with kinematic variable: " << kinematic_variable << std::endl;
                    }
                }

                TEST_CHECK(found_problematic_observable == false);
            }

            /* Test the LaTeX description of observables for unneeded dollar signs */
            {
                auto observables = Observables();

                bool found_problematic_observable = false;

                for (auto const & [name, entry] : observables)
                {
		            if (entry->latex().find("$") == std::string::npos)
                        continue;

                    found_problematic_observable = true;
                    std::cerr << "Found problematic observable: " << name << ", with latex description: " << entry->latex() << "that has (likely) unneeded dollar signs" << std::endl;
                }

                TEST_CHECK(found_problematic_observable == false);

            }
        }

} observable_test;
