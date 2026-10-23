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

#include <eos/utils/expression-kinematic-hook.hh>
#include <eos/utils/expression-parser.hh>
#include <eos/utils/expression.hh>

#include <test/test.hh>

#include <cmath>
#include <memory>
#include <set>

using namespace test;
using namespace eos;
using namespace eos::exp;

class ExpressionKinematicHookTest : public TestCase
{
    public:
        ExpressionKinematicHookTest() :
            TestCase("expression_kinematic_hook_test")
        {
        }

        virtual void
        run() const
        {
            // a derived variable follows the variables it depends on
            {
                Kinematics k{
                    { "m", 2.0 },
                    { "E", 0.0 },
                };
                auto q2 = k.declare_derived("q2", std::make_shared<ExpressionKinematicHook>(parse_expression("{m}^2 + exp({E})"), k));

                TEST_CHECK_NEARLY_EQUAL(q2.evaluate(), 5.0, 1.0e-14);

                // the hook reports the variables that it reads
                const auto hook = k.hook("q2");
                TEST_CHECK_EQUAL(std::set<KinematicVariable::Id>(hook->begin_kinematics(), hook->end_kinematics()), (std::set<KinematicVariable::Id>{ k["m"].id(), k["E"].id() }));

                k.set("m", 3.0);
                k.set("E", 1.0);
                TEST_CHECK_NEARLY_EQUAL(q2.evaluate(), 9.0 + std::exp(1.0), 1.0e-14);

                // a clone follows its own variables
                Kinematics clone = k.clone();
                clone.set("m", 4.0);
                TEST_CHECK_NEARLY_EQUAL(clone["q2"].evaluate(), 16.0 + std::exp(1.0), 1.0e-14);
                TEST_CHECK_NEARLY_EQUAL(k["q2"].evaluate(), 9.0 + std::exp(1.0), 1.0e-14);
            }

            // a variable reached through an alias is bound to the aliased variable
            {
                Kinematics k{
                    { "m_mu", 2.0 }
                };
                k.alias("m", "m_mu");
                auto hook = std::make_shared<ExpressionKinematicHook>(parse_expression("{m}^2"), k);
                k.remove_alias("m");

                auto q2 = k.declare_derived("q2", hook);
                TEST_CHECK_NEARLY_EQUAL(q2.evaluate(), 4.0, 1.0e-14);

                Kinematics clone = k.clone();
                clone.set("m_mu", 3.0);
                TEST_CHECK_NEARLY_EQUAL(clone["q2"].evaluate(), 9.0, 1.0e-14);
            }

            // a clone binds to the variable, even if an alias of the same name is installed
            {
                Kinematics k{
                    {   "m", 2.0 },
                    { "m_e", 5.0 },
                };
                k.declare_derived("q2", std::make_shared<ExpressionKinematicHook>(parse_expression("{m}^2"), k));
                k.alias("m", "m_e");

                Kinematics clone = k.clone();
                TEST_CHECK_NEARLY_EQUAL(clone["q2"].evaluate(), 4.0, 1.0e-14);
            }

            // only kinematic variables, constants, and functions are permitted
            {
                Kinematics k{
                    { "m", 2.0 }
                };

                TEST_CHECK_THROWS(ExpressionError, ExpressionKinematicHook(parse_expression("[[mass::mu]] * {m}"), k));
                TEST_CHECK_THROWS(ExpressionError, ExpressionKinematicHook(parse_expression("<<mass::mu>> * {m}"), k));
                TEST_CHECK_THROWS(UnknownKinematicVariableError, ExpressionKinematicHook(parse_expression("{x}^2"), k));
            }
        }
} expression_kinematic_hook_test;
