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

#include <eos/nonlocal-form-factors/observables.hh>
#include <eos/observable.hh>

#include <test/test.hh>

#include <iostream>
#include <map>
#include <ranges>
#include <set>
#include <string>
#include <variant>

using namespace test;
using namespace eos;

static const std::set<std::string> expected_failures{};

class NonlocalFormFactorsObservablesTest : public TestCase
{
    public:
        NonlocalFormFactorsObservablesTest() :
            TestCase("nonlocal_form_factors_observables_test")
        {
        }

        virtual void
        run() const
        {
            /* Test that every observable in this section can be constructed */
            {
                Parameters parameters = Parameters::Defaults();

                // the nonlocal form factors take 'q' directly, so the entries do not declare it
                const Options base_options{
                    { "q"_ok, "d"_ov }
                };

                std::map<std::string, std::string> failures;

                for (const auto & group : make_nonlocal_form_factors_section())
                {
                    for (const auto & [name, entry] : group)
                    {
                        Kinematics kinematics;

                        for (const auto & kinematic_variable : std::ranges::subrange(entry->begin_kinematic_variables(), entry->end_kinematic_variables()))
                        {
                            kinematics.declare(kinematic_variable, 1.0);
                        }

                        try
                        {
                            try
                            {
                                entry->make(parameters, kinematics, base_options);
                            }
                            catch (UnspecifiedOptionError &)
                            {
                                // retry with the first allowed value for every option that has no default
                                Options options;

                                for (const auto & specification : std::ranges::subrange(entry->begin_options(), entry->end_options()))
                                {
                                    if (specification.default_value.has_value())
                                    {
                                        continue;
                                    }

                                    if (auto * value = std::get_if<qnp::OptionValue>(&specification.allowed_values))
                                    {
                                        options.declare(specification.key, *value);
                                    }
                                    else
                                    {
                                        const auto & values = std::get<std::vector<qnp::OptionValue>>(specification.allowed_values);

                                        if (! values.empty())
                                        {
                                            options.declare(specification.key, values.front());
                                        }
                                    }
                                }

                                entry->make(parameters, kinematics, base_options + options);
                            }
                        }
                        catch (std::exception & e)
                        {
                            failures.insert({ name.str(), e.what() });
                        }
                    }
                }

                bool found_problematic_observable = false;

                for (const auto & [name, message] : failures)
                {
                    if (expected_failures.find(name) != expected_failures.end())
                    {
                        continue;
                    }

                    found_problematic_observable = true;
                    std::cerr << "Found problematic observable: " << name << ", which cannot be constructed: " << message << std::endl;
                }

                for (const auto & name : expected_failures)
                {
                    if (failures.find(name) != failures.end())
                    {
                        continue;
                    }

                    found_problematic_observable = true;
                    std::cerr << "Found problematic observable: " << name << ", which is listed as an expected failure but can now be constructed" << std::endl;
                }

                TEST_CHECK(found_problematic_observable == false);
            }
        }
} nonlocal_form_factors_observables_test;
