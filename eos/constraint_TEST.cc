/* vim: set sw=4 sts=4 et foldmethod=marker foldmarker={{{,}}} : */

/*
 * Copyright (c) 2011, 2013, 2014, 2015, 2017 Danny van Dyk
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
#include <eos/constraint.hh>
#include <eos/statistics/log-likelihood.hh>

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <vector>

using namespace test;
using namespace eos;

class ConstraintDeserializationTest :
    public TestCase
{
    public:
        ConstraintDeserializationTest() :
            TestCase("constraint_deserialization_test")
        {
        }

        virtual void run() const
        {
            // {{{ Gaussian (correct order)
            {
                static const std::string input(
                    "type: Gaussian\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_max: 20, s_min: 15}\n"
                    "options: {form-factors: DM2016, l: mu}\n"
                    "mean: 6e-07\n"
                    "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}\n"
                    "dof: 1"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input == output);
            }
            /// }}}

            // {{{ Gaussian (incorrect order)
            {
                static const std::string input(
                    "type: Gaussian\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_max: 20, s_min: 15}\n"
                    "options: {l: mu, form-factors: DM2016}\n"
                    "mean: 6e-07\n"
                    "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}\n"
                    "dof: 1"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input != output);
            }
            // }}}

            // {{{ Gaussian (double entries 1: kinematics, 2: options)
            {
                static const std::string input1(
                    "type: Gaussian\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_max: 20, s_min: 15, s_max: 71}\n"
                    "options: {form-factors: DM2016, l: mu}\n"
                    "mean: 6e-07\n"
                    "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}\n"
                    "dof: 1"
                );

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node1)));

                static const std::string input2(
                    "type: Gaussian\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_max: 20, s_min: 15}\n"
                    "options: {form-factors: DM2016, l: mu, l: tau}\n"
                    "mean: 6e-07\n"
                    "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}\n"
                    "dof: 1"
                );

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node2)));
            }
            /// }}}

            // {{{ Gaussian (values)
            {
                // {{{ only stat errors
                {
                    static const std::string input(
                        "type: Gaussian\n"
                        "observable: mass::b(MSbar)\n"
                        "kinematics: {}\n"
                        "options: {}\n"
                        "mean: 4.3\n"
                        "sigma-stat: {hi: 0.1, lo: -0.1}\n"
                        "sigma-sys: {hi: 0.0, lo: -0.0}\n"
                        "dof: 1"
                    );

                    YAML::Node node = YAML::Load(input);

                    std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                    TEST_CHECK(nullptr != entry.get());

                    Constraint c = entry->make("Test::Gaussian", Options{ });
                    std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                    TEST_CHECK_EQUAL(1, blocks.size());

                    Parameters p = Parameters::Defaults();
                    LogLikelihood llh(p);
                    llh.add(c);

                    static const double eps = 1e-14;
                    p["mass::b(MSbar)"] = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656, eps);

                    p["mass::b(MSbar)"] = 4.4;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656, eps);
                }
                // }}}

                // {{{ combined sys+stat errors
                {
                    static const std::string input(
                        "type: Gaussian\n"
                        "observable: mass::b(MSbar)\n"
                        "kinematics: {}\n"
                        "options: {}\n"
                        "mean: 4.3\n"
                        "sigma-stat: {hi: 0.4, lo: -0.4}\n"
                        "sigma-sys: {hi: 0.3, lo: -0.3}\n"
                        "dof: 1"
                    );

                    YAML::Node node = YAML::Load(input);

                    std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                    TEST_CHECK(nullptr != entry.get());

                    Constraint c = entry->make("Test::Gaussian", Options{ });
                    std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                    TEST_CHECK_EQUAL(1, blocks.size());

                    Parameters p = Parameters::Defaults();
                    LogLikelihood llh(p);
                    llh.add(c);

                    static const double eps = 1e-14;
                    p["mass::b(MSbar)"] = 3.8;
                    TEST_CHECK_NEARLY_EQUAL(llh(), -0.72579135264472743, eps);

                    p["mass::b(MSbar)"] = 4.8;
                    TEST_CHECK_NEARLY_EQUAL(llh(), -0.72579135264472743, eps);
                }
                // }}}
            }
            // }}}

            // {{{ LogGamma (correct order)
            {
                static const std::string input(
                    "type: LogGamma\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_max: 20, s_min: 15}\n"
                    "options: {form-factors: DM2016, l: mu}\n"
                    "mode: 6e-07\n"
                    "sigma: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "dof: 1"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ LogGamma (incorrect order)
            {
                static const std::string input(
                    "type: LogGamma\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_min: 15, s_max: 20}\n"
                    "options: {form-factors: DM2016, l: mu}\n"
                    "mode: 6e-07\n"
                    "sigma: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "dof: 1"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input != output);
            }
            // }}}

            // {{{ LogGamma (double entries 1: kinematics, 2: options)
            {
                static const std::string input1(
                    "type: LogGamma\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_min: 15, s_max: 20, s_max: 13.5}\n"
                    "options: {form-factors: DM2016, l: mu}\n"
                    "mode: 6e-07\n"
                    "sigma: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "dof: 1"
                );

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node1)));

                static const std::string input2(
                    "type: LogGamma\n"
                    "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                    "kinematics: {s_min: 15, s_max: 20}\n"
                    "options: {form-factors: DM2016, l: mu, l: tau}\n"
                    "mode: 6e-07\n"
                    "sigma: {hi: 1.2e-08, lo: -3.4e-08}\n"
                    "dof: 1"
                );

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node2)));
            }
            // }}}

            // {{{ LogGamma (values)
            {
                static const std::string input(
                    "type: LogGamma\n"
                    "observable: mass::b(MSbar)\n"
                    "kinematics: {}\n"
                    "options: {}\n"
                    "mode: 0.53\n"
                    "sigma: {hi: 0.1, lo: 0.19}\n"
                    "dof: 1"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint c = entry->make("Test::LogGamma", Options{ });
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                static const double eps = 5e-4;
                p["mass::b(MSbar)"] = 0.57;
                TEST_CHECK_RELATIVE_ERROR(llh(), +1.005543554, eps);
            }
            // }}}

            // {{{ Amoroso (correct order)
            {
                static const std::string input(
                    "type: Amoroso\n"
                    "observable: B_q->ll::BR@Untagged\n"
                    "kinematics: {}\n"
                    "options: {l: mu, q: s}\n"
                    "physical-limit: 0\n"
                    "theta: 1.54243e-10\n"
                    "alpha: 19.4025\n"
                    "beta: 1.0048"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ Amoroso (incorrect order)
            {
                static const std::string input(
                    "type: Amoroso\n"
                    "observable: B_q->ll::BR@Untagged\n"
                    "kinematics: {}\n"
                    "options: {q: s, l: mu}\n"
                    "physical-limit: 0\n"
                    "theta: 1.54243e-10\n"
                    "alpha: 19.4025\n"
                    "beta: 1.0048"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input != output);
            }
            // }}}

            // {{{ Amoroso (double entries 1: kinematics, 2: options)
            {
                static const std::string input1(

                    "type: Amoroso\n"
                    "observable: B_q->ll::BR@Untagged\n"
                    "kinematics: {none: 1, none: 2}\n"
                    "options: {q: s, l: mu}\n"
                    "physical-limit: 0\n"
                    "theta: 1.54243e-10\n"
                    "alpha: 19.4025\n"
                    "beta: 1.0048"
                );

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node1)));

                static const std::string input2(
                    "type: Amoroso\n"
                    "observable: B_q->ll::BR@Untagged\n"
                    "kinematics: {}\n"
                    "options: {l: mu, q: s, q: s}\n"
                    "physical-limit: 0\n"
                    "theta: 1.54243e-10\n"
                    "alpha: 19.4025\n"
                    "beta: 1.0048"
                );

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node2)));
            }
            // }}}

            // {{{ Amoroso (value)
            {
                static const std::string input(
                    "type: Amoroso\n"
                    "observable: mass::b(MSbar)\n"
                    "kinematics: {}\n"
                    "options: {}\n"
                    "physical-limit: 0.0\n"
                    "theta: 2.9708273062\n"
                    "alpha: 8.2392613044e-01\n"
                    "beta: 1.6993290032"
                );

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint c = entry->make("Test::Amoroso", Options{ });
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation at mode
                p["mass::b(MSbar)"] = 1.268439;
                TEST_CHECK_RELATIVE_ERROR(llh(), std::log(2.824624787700217e-01), 1e-8);

                p["mass::b(MSbar)"] = 2.536877;
                TEST_CHECK_RELATIVE_ERROR(llh(), -1.516059, 1e-6);

                p["mass::b(MSbar)"] = 3.805315;
                TEST_CHECK_RELATIVE_ERROR(llh(), -2.112174, 1e-6);
            }
            // }}}
        }
} constraint_deserialization_test;

class ConstraintTest :
    public TestCase
{
    public:
        ConstraintTest() :
            TestCase("constraint_test")
        {
        }

        virtual void run() const
        {
            /* Test making constraints */
            {
                std::cout << "# Constraints :" << std::endl;

                Options o;
                auto constraints = Constraints();
                unsigned n = 0;

                for (auto cf = constraints.begin(); cf != constraints.end(); ++cf, ++n)
                {
                    std::cout << "#  " << cf->first.full() << ": ";

                    Constraint c = Constraint::make(cf->first, o);
                    TEST_CHECK_EQUAL(c.name(), cf->first);
                    TEST_CHECK(std::distance(c.begin_observables(), c.end_observables()) > 0);
                    TEST_CHECK(std::distance(c.begin_blocks(), c.end_blocks()) > 0);

                    for (auto o = c.begin_observables(), o_end = c.end_observables(); o != o_end ; ++o)
                    {
                        std::cout << (**o).name() << '['
                                << (**o).kinematics().as_string() << ']'
                                << " with options: " << (**o).options().as_string();
                    }
                    for (auto b = c.begin_blocks(), b_end = c.end_blocks(); b != b_end ; ++b)
                    {
                        std::cout << ", " << (**b).as_string();
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
                std::cout << "# Found " << n << " constraints" << std::endl;
            }
        }
} constraint_test;
