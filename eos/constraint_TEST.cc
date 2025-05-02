/* vim: set sw=4 sts=4 et foldmethod=marker foldmarker={{{,}}} : */

/*
 * Copyright (c) 2011-2025 Danny van Dyk
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

#include <eos/constraint.hh>
#include <eos/maths/power-of.hh>
#include <eos/observable.hh>
#include <eos/statistics/log-likelihood.hh>

#include <test/test.hh>

#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

using namespace test;
using namespace eos;

class ConstraintDeserializationTest : public TestCase
{
    public:
        ConstraintDeserializationTest() :
            TestCase("constraint_deserialization_test")
        {
        }

        virtual void
        run() const
        {
            // {{{ Gaussian (correct order, non-ASCI character present)
            {
                static const std::string input("type: Gaussian\n"
                                               "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                               "kinematics: {s_max: 20, s_min: 15}\n"
                                               "options: {form-factors: DM2016, l: mu}\n"
                                               "mean: 6e-07\n"
                                               "sigma-stat: {hi: 1.2e-08, lo: −3.4e-08}\n"
                                               "sigma-sys: {hi: 4.6e-08, lo: −7.8e-08}"); // sigma-stat.lo and sigma-sys.lo contain a U+2212 (minus sign, '−'), rather than a '-'.

                TEST_CHECK_THROWS(ConstraintEntryEncodingError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", input)));
            }
            /// }}}

            // {{{ Gaussian (correct order)
            {
                static const std::string input("type: Gaussian\n"
                                               "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                               "kinematics: {s_max: 20, s_min: 15}\n"
                                               "options: {form-factors: DM2016, l: mu}\n"
                                               "mean: 6e-07\n"
                                               "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                                               "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}\n"
                                               "references:\n"
                                               "  []");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                std::cout << output << std::endl;

                TEST_CHECK(input == output);
            }
            /// }}}

            // {{{ Gaussian (incorrect order)
            {
                static const std::string input("type: Gaussian\n"
                                               "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                               "kinematics: {s_max: 20, s_min: 15}\n"
                                               "options: {l: mu, form-factors: DM2016}\n"
                                               "mean: 6e-07\n"
                                               "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                                               "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}");

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
                static const std::string input1("type: Gaussian\n"
                                                "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                                "kinematics: {s_max: 20, s_min: 15, s_max: 71}\n"
                                                "options: {form-factors: DM2016, l: mu}\n"
                                                "mean: 6e-07\n"
                                                "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                                                "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}");

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node1)));

                static const std::string input2("type: Gaussian\n"
                                                "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                                "kinematics: {s_max: 20, s_min: 15}\n"
                                                "options: {form-factors: DM2016, l: mu, l: tau}\n"
                                                "mean: 6e-07\n"
                                                "sigma-stat: {hi: 1.2e-08, lo: -3.4e-08}\n"
                                                "sigma-sys: {hi: 4.6e-08, lo: -7.8e-08}");

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node2)));
            }
            /// }}}

            // {{{ Gaussian (values)
            {
                // {{{ only stat errors
                {
                    static const std::string input("type: Gaussian\n"
                                                   "observable: mass::b(MSbar)\n"
                                                   "kinematics: {}\n"
                                                   "options: {}\n"
                                                   "mean: 4.3\n"
                                                   "sigma-stat: {hi: 0.1, lo: -0.1}\n"
                                                   "sigma-sys: {hi: 0.0, lo: -0.0}");

                    YAML::Node node = YAML::Load(input);

                    std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                    TEST_CHECK(nullptr != entry.get());

                    Constraint                         c = entry->make("Test::Gaussian", Options{});
                    std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                    TEST_CHECK_EQUAL(1, blocks.size());

                    Parameters    p = Parameters::Defaults();
                    LogLikelihood llh(p);
                    llh.add(c);

                    static const double eps = 1e-14;
                    p["mass::b(MSbar)"]     = 4.2;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656, eps);

                    p["mass::b(MSbar)"] = 4.4;
                    TEST_CHECK_NEARLY_EQUAL(llh(), +0.88364655978937656, eps);
                }
                // }}}

                // {{{ combined sys+stat errors
                {
                    static const std::string input("type: Gaussian\n"
                                                   "observable: mass::b(MSbar)\n"
                                                   "kinematics: {}\n"
                                                   "options: {}\n"
                                                   "mean: 4.3\n"
                                                   "sigma-stat: {hi: 0.4, lo: -0.4}\n"
                                                   "sigma-sys: {hi: 0.3, lo: -0.3}");

                    YAML::Node node = YAML::Load(input);

                    std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Gaussian", node));
                    TEST_CHECK(nullptr != entry.get());

                    Constraint                         c = entry->make("Test::Gaussian", Options{});
                    std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                    TEST_CHECK_EQUAL(1, blocks.size());

                    Parameters    p = Parameters::Defaults();
                    LogLikelihood llh(p);
                    llh.add(c);

                    static const double eps = 1e-14;
                    p["mass::b(MSbar)"]     = 3.8;
                    TEST_CHECK_NEARLY_EQUAL(llh(), -0.72579135264472743, eps);

                    p["mass::b(MSbar)"] = 4.8;
                    TEST_CHECK_NEARLY_EQUAL(llh(), -0.72579135264472743, eps);
                }
                // }}}
            }
            // }}}

            // {{{ LogGamma (correct order)
            {
                static const std::string input("type: LogGamma\n"
                                               "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                               "kinematics: {q2_max: 20, q2_min: 15}\n"
                                               "options: {form-factors: DM2016, l: mu}\n"
                                               "mode: 6e-07\n"
                                               "sigma: {hi: 1.2e-08, lo: 3.4e-08}\n"
                                               "alpha: 0.17132\n"
                                               "lambda: 5.78825e-09\n"
                                               "references:\n"
                                               "  []");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                std::cerr << output << std::endl;
                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ LogGamma (incorrect order)
            {
                static const std::string input("type: LogGamma\n"
                                               "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                               "kinematics: {q2_min: 15, q2_max: 20}\n"
                                               "options: {form-factors: DM2016, l: mu}\n"
                                               "mode: 6e-07\n"
                                               "sigma: {hi: 1.2e-08, lo: 3.4e-08}\n"
                                               "alpha: 0.17132\n"
                                               "lambda: 5.78825e-09");

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
                static const std::string input1("type: LogGamma\n"
                                                "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                                "kinematics: {q2_min: 15, q2_max: 20, q2_max: 13.5}\n"
                                                "options: {form-factors: DM2016, l: mu}\n"
                                                "mode: 6e-07\n"
                                                "sigma: {hi: 1.2e-08, lo: 3.4e-08}\n"
                                                "alpha: 0.17132\n"
                                                "lambda: 5.78825e-09");

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node1)));

                static const std::string input2("type: LogGamma\n"
                                                "observable: Lambda_b->Lambdall::BR@LowRecoil\n"
                                                "kinematics: {q2_min: 15, q2_max: 20}\n"
                                                "options: {form-factors: DM2016, l: mu, l: tau}\n"
                                                "mode: 6e-07\n"
                                                "sigma: {hi: 1.2e-08, lo: 3.4e-08}\n"
                                                "alpha: 0.17132\n"
                                                "lambda: 5.78825e-09");

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node2)));
            }
            // }}}

            // {{{ LogGamma (values)
            {
                static const std::string input("type: LogGamma\n"
                                               "observable: mass::b(MSbar)\n"
                                               "kinematics: {}\n"
                                               "options: {}\n"
                                               "mode: 0.53\n"
                                               "sigma: {hi: 0.1, lo: 0.19}\n"
                                               "alpha: 0.383056\n"
                                               "lambda: 0.0687907");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::LogGamma", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::LogGamma", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                static const double eps = 5e-4;
                p["mass::b(MSbar)"]     = 0.57;
                TEST_CHECK_RELATIVE_ERROR(llh(), +1.005543554, eps);
            }
            // }}}

            // {{{ Amoroso (correct order)
            {
                static const std::string input("type: Amoroso\n"
                                               "observable: B_q->ll::BR@Untagged\n"
                                               "kinematics: {}\n"
                                               "options: {l: mu, q: s}\n"
                                               "physical-limit: 0\n"
                                               "theta: 1.54243e-10\n"
                                               "alpha: 19.4025\n"
                                               "beta: 1.0048\n"
                                               "references:\n"
                                               "  []");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                std::cerr << output << std::endl;
                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ Amoroso (incorrect order)
            {
                static const std::string input("type: Amoroso\n"
                                               "observable: B_q->ll::BR@Untagged\n"
                                               "kinematics: {}\n"
                                               "options: {q: s, l: mu}\n"
                                               "physical-limit: 0\n"
                                               "theta: 1.54243e-10\n"
                                               "alpha: 19.4025\n"
                                               "beta: 1.0048");

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
                        "beta: 1.0048");

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node1)));

                static const std::string input2("type: Amoroso\n"
                                                "observable: B_q->ll::BR@Untagged\n"
                                                "kinematics: {}\n"
                                                "options: {l: mu, q: s, q: s}\n"
                                                "physical-limit: 0\n"
                                                "theta: 1.54243e-10\n"
                                                "alpha: 19.4025\n"
                                                "beta: 1.0048");

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node2)));
            }
            // }}}

            // {{{ Amoroso (value)
            {
                static const std::string input("type: Amoroso\n"
                                               "observable: mass::b(MSbar)\n"
                                               "kinematics: {}\n"
                                               "options: {}\n"
                                               "physical-limit: 0.0\n"
                                               "theta: 2.9708273062\n"
                                               "alpha: 8.2392613044e-01\n"
                                               "beta: 1.6993290032");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Amoroso", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::Amoroso", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
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

            // {{{ MultivariateGaussian (correct order)
            {
                static const std::string input("type: MultivariateGaussian\n"
                                               "observables:\n"
                                               "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                               "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                               "kinematics:\n"
                                               "  - {s: 13.5}\n"
                                               "  - {s: 20.5}\n"
                                               "options:\n"
                                               "  - {form-factors: BFvD2014}\n"
                                               "  - {form-factors: BFvD2014}\n"
                                               "means: [0.73, 1.4]\n"
                                               "sigma-stat-hi: [0.2, 0.2]\n"
                                               "sigma-stat-lo: [0.2, 0.2]\n"
                                               "sigma-sys: [0, 0]\n"
                                               "correlations:\n"
                                               "  - [1, 0]\n"
                                               "  - [0, 1]\n"
                                               "references:\n"
                                               "  []\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                std::cerr << output << std::endl;
                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ MultivariateGaussian (incorrect order)
            {
                static const std::string input("type: MultivariateGaussian\n"
                                               "observables:\n"
                                               "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                               "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                               "kinematics:\n"
                                               "  - {s: 13.5}\n"
                                               "  - {s: 20.5}\n"
                                               "options:\n"
                                               "  - {l: mu, form-factors: BFvD2014}\n"
                                               "  - {q: s, form-factors: BFvD2014}\n"
                                               "means: [0.73, 1.4]\n"
                                               "sigma-stat-hi: [0.2, 0.2]\n"
                                               "sigma-stat-lo: [0.2, 0.2]\n"
                                               "sigma-sys: [0, 0]\n"
                                               "correlations:\n"
                                               "  - [1, 0]\n"
                                               "  - [0, 1]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input != output);
            }
            // }}}

            // {{{ MultivariateGaussian (double entries 1: kinematics, 2: options)
            {
                static const std::string input1("type: MultivariateGaussian\n"
                                                "observables:\n"
                                                "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                                "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                                "kinematics:\n"
                                                "  - {s: 13.5, s: 20.5}\n"
                                                "  - {s: 20.5, s: 13.5}\n"
                                                "options:\n"
                                                "  - {form-factors: BFvD2014}\n"
                                                "  - {form-factors: BFvD2014}\n"
                                                "means: [0.73, 1.4]\n"
                                                "sigma-stat-hi: [0.2, 0.2]\n"
                                                "sigma-stat-lo: [0.2, 0.2]\n"
                                                "sigma-sys: [0, 0]\n"
                                                "correlations:\n"
                                                "  - [1, 0]\n"
                                                "  - [0, 1]\n"
                                                "dof: 2");

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian", node1)));

                static const std::string input2("type: MultivariateGaussian\n"
                                                "observables:\n"
                                                "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                                "  - Lambda_b->Lambda::f_perp^V(s)\n"
                                                "kinematics:\n"
                                                "  - {s: 13.5}\n"
                                                "  - {s: 20.5}\n"
                                                "options:\n"
                                                "  - {form-factors: BFvD2014, form-factors: DM2016}\n"
                                                "  - {form-factors: BFvD2014, form-factors: DM2016}\n"
                                                "means: [0.73, 1.4]\n"
                                                "sigma-stat-hi: [0.2, 0.2]\n"
                                                "sigma-stat-lo: [0.2, 0.2]\n"
                                                "sigma-sys: [0, 0]\n"
                                                "correlations:\n"
                                                "  - [1, 0]\n"
                                                "  - [0, 1]\n"
                                                "dof: 2");

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError, std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian", node2)));
            }
            // }}}

            // {{{ MultivariateGaussian (value)
            {
                static const std::string input("type: MultivariateGaussian\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "means: [4.3, 1.1]\n"
                                               "sigma-stat-hi: [0.1, 0.05]\n"
                                               "sigma-stat-lo: [0.1, 0.05]\n"
                                               "sigma-sys: [0, 0]\n"
                                               "correlations:\n"
                                               "  - [1, 0.6]\n"
                                               "  - [0.6, 1]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::MultivariateGaussian", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation at mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ MultivariateGaussian (begin and end options)
            {
                static const std::string input("type: MultivariateGaussian\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "  - mass::u\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "means: [4.3, 1.1, 0.002]\n"
                                               "sigma-stat-hi: [0.1, 0.05, 0.001]\n"
                                               "sigma-stat-lo: [0.1, 0.05, 0.001]\n"
                                               "sigma-sys: [0, 0, 0]\n"
                                               "correlations:\n"
                                               "  - [ 1.0,  0.6,  -0.1 ]\n"
                                               "  - [ 0.6,  1.0,   0.01]\n"
                                               "  - [-0.1,  0.01,  1.0 ]\n"
                                               "dof: 3");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::MultivariateGaussian",
                                           Options{
                                                                       { "begin"_ok, "0" },
                                                                       {   "end"_ok, "2" }
                });
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation at mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (correct order)
            {
                static const std::string input("type: MultivariateGaussian(Covariance)\n"
                                               "observables:\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_0(s)\n"
                                               "  - B->D::f_0(s)\n"
                                               "  - B->D::f_0(s)\n"
                                               "kinematics:\n"
                                               "  - {s: 0}\n"
                                               "  - {s: 4}\n"
                                               "  - {s: 8}\n"
                                               "  - {s: 11.62}\n"
                                               "  - {s: 4}\n"
                                               "  - {s: 8}\n"
                                               "  - {s: 11.62}\n"
                                               "options:\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "means: [0.665, 0.798, 0.972, 1.177, 0.729, 0.81, 0.901]\n"
                                               "covariance:\n"
                                               "  - [0.001128, 0.001042, 0.000923, 0.0007727, 0.001093, 0.001063, 0.001045]\n"
                                               "  - [0.001042, 0.001079, 0.001108, 0.001123, 0.001026, 0.001017, 0.001021]\n"
                                               "  - [0.000923, 0.001108, 0.001331, 0.001576, 0.000931, 0.0009511, 0.0009865]\n"
                                               "  - [0.000773, 0.001123, 0.001576, 0.002112, 0.000811, 0.0008681, 0.0009425]\n"
                                               "  - [0.001093, 0.001026, 0.0009307, 0.0008108, 0.001126, 0.001165, 0.00121]\n"
                                               "  - [0.001063, 0.001017, 0.0009511, 0.0008681, 0.001165, 0.001283, 0.00141]\n"
                                               "  - [0.001045, 0.001021, 0.0009865, 0.0009425, 0.00121, 0.00141, 0.001635]\n"
                                               "references:\n"
                                               "  []\n"
                                               "dof: 0");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());
                std::cerr << output << std::endl;
                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (incorrect order)
            {
                static const std::string input("type: MultivariateGaussian(Covariance)\n"
                                               "observables:\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_+(s)\n"
                                               "  - B->D::f_0(s)\n"
                                               "  - B->D::f_0(s)\n"
                                               "  - B->D::f_0(s)\n"
                                               "kinematics:\n"
                                               "  - {s: 0}\n"
                                               "  - {s: 4}\n"
                                               "  - {s: 8}\n"
                                               "  - {s: 11.62}\n"
                                               "  - {s: 4}\n"
                                               "  - {s: 8}\n"
                                               "  - {s: 11.62}\n"
                                               "options:\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {form-factors: BCL2008}\n"
                                               "  - {q: s, form-factors: BCL2008}\n"
                                               "means: [0.665, 0.798, 0.972, 1.177, 0.729, 0.8100000000000001, 0.901]\n"
                                               "covariance:\n"
                                               "  - [0.001128, 0.001042, 0.000923, 0.0007727, 0.001093, 0.001063, 0.001045]\n"
                                               "  - [0.001042, 0.001079, 0.001108, 0.001123, 0.001026, 0.001017, 0.001021]\n"
                                               "  - [0.000923, 0.001108, 0.001331, 0.001576, 0.0009307, 0.0009511, 0.0009865]\n"
                                               "  - [0.0007727, 0.001123, 0.001576, 0.002112, 0.0008108, 0.0008681, 0.0009425]\n"
                                               "  - [0.001093, 0.001026, 0.0009307, 0.0008108, 0.001126, 0.001165, 0.00121]\n"
                                               "  - [0.001063, 0.001017, 0.0009511, 0.0008681, 0.001165, 0.001283, 0.00141]\n"
                                               "  - [0.001045, 0.001021, 0.0009865, 0.0009425, 0.00121, 0.00141, 0.001635]\n"
                                               "dof: 0");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input != output);
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (double entries 1: kinematics, 2: options)
            {
                static const std::string input1("type: MultivariateGaussian(Covariance)\n"
                                                "observables:\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_0(s)\n"
                                                "  - B->D::f_0(s)\n"
                                                "  - B->D::f_0(s)\n"
                                                "kinematics:\n"
                                                "  - {s: 0}\n"
                                                "  - {s: 4}\n"
                                                "  - {s: 8}\n"
                                                "  - {s: 11.62}\n"
                                                "  - {s: 4, s: 8}\n"
                                                "  - {s: 8}\n"
                                                "  - {s: 11.62}\n"
                                                "options:\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "means: [0.665, 0.798, 0.972, 1.177, 0.729, 0.8100000000000001, 0.901]\n"
                                                "covariance:\n"
                                                "  - [0.001128, 0.001042, 0.000923, 0.0007727, 0.001093, 0.001063, 0.001045]\n"
                                                "  - [0.001042, 0.001079, 0.001108, 0.001123, 0.001026, 0.001017, 0.001021]\n"
                                                "  - [0.000923, 0.001108, 0.001331, 0.001576, 0.0009307, 0.0009511, 0.0009865]\n"
                                                "  - [0.0007727, 0.001123, 0.001576, 0.002112, 0.0008108, 0.0008681, 0.0009425]\n"
                                                "  - [0.001093, 0.001026, 0.0009307, 0.0008108, 0.001126, 0.001165, 0.00121]\n"
                                                "  - [0.001063, 0.001017, 0.0009511, 0.0008681, 0.001165, 0.001283, 0.00141]\n"
                                                "  - [0.001045, 0.001021, 0.0009865, 0.0009425, 0.00121, 0.00141, 0.001635]\n"
                                                "dof: 0");

                YAML::Node node1 = YAML::Load(input1);

                TEST_CHECK_THROWS(ConstraintDeserializationError,
                                  std::shared_ptr<ConstraintEntry>
                                          entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node1)));

                static const std::string input2("type: MultivariateGaussian(Covariance)\n"
                                                "observables:\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_+(s)\n"
                                                "  - B->D::f_0(s)\n"
                                                "  - B->D::f_0(s)\n"
                                                "  - B->D::f_0(s)\n"
                                                "kinematics:\n"
                                                "  - {s: 0}\n"
                                                "  - {s: 4}\n"
                                                "  - {s: 8}\n"
                                                "  - {s: 11.62}\n"
                                                "  - {s: 4}\n"
                                                "  - {s: 8}\n"
                                                "  - {s: 11.62}\n"
                                                "options:\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008, form-factors: invalid}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "  - {form-factors: BCL2008}\n"
                                                "means: [0.665, 0.798, 0.972, 1.177, 0.729, 0.8100000000000001, 0.901]\n"
                                                "covariance:\n"
                                                "  - [0.001128, 0.001042, 0.000923, 0.0007727, 0.001093, 0.001063, 0.001045]\n"
                                                "  - [0.001042, 0.001079, 0.001108, 0.001123, 0.001026, 0.001017, 0.001021]\n"
                                                "  - [0.000923, 0.001108, 0.001331, 0.001576, 0.0009307, 0.0009511, 0.0009865]\n"
                                                "  - [0.0007727, 0.001123, 0.001576, 0.002112, 0.0008108, 0.0008681, 0.0009425]\n"
                                                "  - [0.001093, 0.001026, 0.0009307, 0.0008108, 0.001126, 0.001165, 0.00121]\n"
                                                "  - [0.001063, 0.001017, 0.0009511, 0.0008681, 0.001165, 0.001283, 0.00141]\n"
                                                "  - [0.001045, 0.001021, 0.0009865, 0.0009425, 0.00121, 0.00141, 0.001635]\n"
                                                "dof: 0");

                YAML::Node node2 = YAML::Load(input2);

                TEST_CHECK_THROWS(ConstraintDeserializationError,
                                  std::shared_ptr<ConstraintEntry>
                                          entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node2)));
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (value)
            {
                static const std::string input("type: MultivariateGaussian(Covariance)\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "means: [4.3, 1.1]\n"
                                               "covariance:\n"
                                               "  - [0.0100, 0.0030]\n"
                                               "  - [0.0030, 0.0025]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::MultivariateGaussian(Covariance)", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation off mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (value; trivial response matrix)
            {
                static const std::string input("type: MultivariateGaussian(Covariance)\n"
                                               "observables:\n"
                                               "  - mass::s(2GeV)\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "means: [4.3, 1.1]\n"
                                               "covariance:\n"
                                               "  - [0.0100, 0.0030]\n"
                                               "  - [0.0030, 0.0025]\n"
                                               "response:\n"
                                               "  - [0, 1, 0]\n"
                                               "  - [0, 0, 1]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::MultivariateGaussian(Covariance)", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation off mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (value; non-trivial response matrix)
            {
                static const std::string input("type: MultivariateGaussian(Covariance)\n"
                                               "observables:\n"
                                               "  - mass::s(2GeV)\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "means: [4.3, 1.1]\n"
                                               "covariance:\n"
                                               "  - [0.0100, 0.0030]\n"
                                               "  - [0.0030, 0.0025]\n"
                                               "response:\n"
                                               "  - [0, +0.5, +0.5]\n"
                                               "  - [0, -0.5, +0.5]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::MultivariateGaussian(Covariance)", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation off mode
                p["mass::b(MSbar)"] = 3.3;
                p["mass::c"]        = 5.9;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ MultivariateGaussian(Covariance) (value; sub-sample)
            {
                static const std::string input("type: MultivariateGaussian(Covariance)\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "  - mass::K_d\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "means: [4.3, 1.1, 6.2]\n"
                                               "covariance:\n"
                                               "  - [0.0100, 0.0030, 0.0060]\n"
                                               "  - [0.0030, 0.0025, 0.0015]\n"
                                               "  - [0.0060, 0.0015, 0.0105]\n"
                                               "dof: 3");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::MultivariateGaussian(Covariance)", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::MultivariateGaussian(Covariance)",
                                           Options{
                                                                       { "begin"_ok, "0" },
                                                                       {   "end"_ok, "2" }
                });
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation off mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ UniformBound (correct order)
            {
                static const std::string input("type: UniformBound\n"
                                               "observables:\n"
                                               "  - b->c::Bound[0^+]@CLN\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "bound: 0.7\n"
                                               "uncertainty: 0.1\n"
                                               "references:\n"
                                               "  []");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::UniformBound", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());
                std::cerr << output << std::endl;
                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ UniformBound (incorrect order)
            {
                static const std::string input("type: UniformBound\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "observables:\n"
                                               "  - b->c::Bound[0^+]@CLN\n"
                                               "bound: 0.7\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "uncertainty: 0.1");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::UniformBound", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                TEST_CHECK(input != output);
            }
            // }}}

            // {{{ UniformBound (value)
            {
                static const std::string input("type: UniformBound\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "bound: 1.0\n"
                                               "uncertainty: 0.1");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::UniformBound", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::UniformBound", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                p["mass::b(MSbar)"] = 0.0;
                TEST_CHECK_NEARLY_EQUAL(llh(), 0.0, 1e-7);

                p["mass::b(MSbar)"] = 1.1;
                TEST_CHECK_NEARLY_EQUAL(llh(), -0.5, 1e-7);

                p["mass::b(MSbar)"] = 11.;
                TEST_CHECK_NEARLY_EQUAL(llh(), -5000.0, 1e-7);
            }
            // }}}

            // {{{ UniformBound (value with non-trivial kinematics & options)
            {
                // abusing the B->pi form factor to have an observable with non-trivial kinematics and options.
                static const std::string input("type: UniformBound\n"
                                               "observables:\n"
                                               "  - B->pi::f_+(q2)\n"
                                               "kinematics:\n"
                                               "  - {q2: 1.0}\n"
                                               "options:\n"
                                               "  - {form-factors: BSZ2015}\n"
                                               "bound: 0.1\n"
                                               "uncertainty: 0.01");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::UniformBound", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::UniformBound", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                ObservablePtr o = Observable::make("B->pi::f_+(q2);form-factors=BSZ2015",
                                                   p,
                                                   Kinematics{
                                                       { "q2", 1.0 }
                },
                                                   Options{});
                TEST_CHECK_NEARLY_EQUAL(llh(), -0.5 * power_of<2>((o->evaluate() - 0.1) / 0.01), 1e-17);
            }
            // }}}

            // {{{ Mixture (correct order)
            {
                static const std::string input("type: Mixture\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "components:\n"
                                               "  - means: [4.3, 1.1]\n"
                                               "    covariance:\n"
                                               "      - [0.01, 0.003]\n"
                                               "      - [0.003, 0.0025]\n"
                                               "  - means: [4.2, 1.1]\n"
                                               "    covariance:\n"
                                               "      - [0.04, 0.006]\n"
                                               "      - [0.006, 0.0025]\n"
                                               "weights: [0.6, 0.4]\n"
                                               "test statistics:\n"
                                               "  sigma: []\n"
                                               "  densities: []\n"
                                               "references:\n"
                                               "  []\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Mixture", node));
                TEST_CHECK(nullptr != entry.get());

                YAML::Emitter out;
                entry->serialize(out);

                std::string output(out.c_str());

                std::cerr << output << std::endl;
                TEST_CHECK(input == output);
            }
            // }}}

            // {{{ Mixture (value; trivial cross check against MultivariateGaussian(Covariance))
            {
                static const std::string input("type: Mixture\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "components:\n"
                                               "  - means: [4.3, 1.1]\n"
                                               "    covariance:\n"
                                               "      - [0.0100, 0.0030]\n"
                                               "      - [0.0030, 0.0025]\n"
                                               "test statistics:\n"
                                               "  sigma: []\n"
                                               "  densities: []\n"
                                               "weights: [1.0]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Mixture", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::Mixture", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation off mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.597666149, 1e-8);
            }
            // }}}

            // {{{ Mixture (value; non-trivial cross check using two components and naive test statistics values)
            {
                static const std::string input("type: Mixture\n"
                                               "observables:\n"
                                               "  - mass::b(MSbar)\n"
                                               "  - mass::c\n"
                                               "kinematics:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "options:\n"
                                               "  - {}\n"
                                               "  - {}\n"
                                               "components:\n"
                                               "  - means: [4.3, 1.1]\n"
                                               "    covariance:\n"
                                               "      - [0.01, 0.003]\n"
                                               "      - [0.003, 0.0025]\n"
                                               "  - means: [4.2, 1.1]\n"
                                               "    covariance:\n"
                                               "      - [0.04, 0.006]\n"
                                               "      - [0.006, 0.0025]\n"
                                               "test statistics:\n"
                                               "  sigma: [0.5, 1., 1.5]\n"
                                               "  densities: [-13., -8., -1.]\n"
                                               "weights: [0.6, 0.4]\n"
                                               "dof: 2");

                YAML::Node node = YAML::Load(input);

                std::shared_ptr<ConstraintEntry> entry(ConstraintEntry::FromYAML("Test::Mixture", node));
                TEST_CHECK(nullptr != entry.get());

                Constraint                         c = entry->make("Test::Mixture", Options{});
                std::vector<LogLikelihoodBlockPtr> blocks(c.begin_blocks(), c.end_blocks());
                TEST_CHECK_EQUAL(1, blocks.size());

                Parameters    p = Parameters::Defaults();
                LogLikelihood llh(p);
                llh.add(c);

                // evaluation off mode
                p["mass::b(MSbar)"] = 4.6;
                p["mass::c"]        = 1.3;
                TEST_CHECK_NEARLY_EQUAL(llh(), -4.779399451, 1e-8);
            }
            // }}}
        }
} constraint_deserialization_test;

class ConstraintTest : public TestCase
{
    public:
        ConstraintTest() :
            TestCase("constraint_test")
        {
        }

        virtual void
        run() const
        {
            /* Test making constraints */
            /*{
                std::cout << "# Constraints :" << std::endl;

                Options o;
                auto constraints = Constraints();
                unsigned n = 0;

                for (auto cf = constraints.begin(); cf != constraints.end(); ++cf, ++n)
                {
                    std::cout << "#  " << cf->first.full() << ": ";

                    auto suffix = cf->first.suffix_part();
                    if (0 < suffix.str().size())
                    {
                        TEST_CHECK_NO_THROW(auto rn = ReferenceName(suffix.str()));
                    }

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
            }*/

            /* Test retrieving ConstraintEntry by name */
            {
                try
                {
                    auto constraints = Constraints();

                    static const std::vector<QualifiedName> names{ "B->pi::f_+@IKMvD:2014A", "B->K::f_0+f_++f_T@HPQCD:2013A" };

                    for (auto & n : names)
                    {
                        std::shared_ptr<const ConstraintEntry> c;
                        TEST_CHECK_NO_THROW(c = constraints[n]);
                        TEST_CHECK(c.get() != nullptr);
                    }
                }
                catch (std::exception & e)
                {
                    std::cerr << "Caught unexpected exception: " << e.what() << std::endl;
                    throw e;
                }
                catch (...)
                {
                    std::cerr << "Caught unexpected unknown exception" << std::endl;
                }
            }
        }
} constraint_test;
