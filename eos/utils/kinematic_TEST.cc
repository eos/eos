/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2026 Danny van Dyk
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

#include <eos/utils/kinematic.hh>

#include <test/test.hh>

#include <cmath>
#include <iterator>
#include <memory>

using namespace test;
using namespace eos;

namespace
{
    // Squares the value of another kinematic variable.
    class SquareHook : public KinematicHook
    {
        private:
            KinematicVariable _variable;

        public:
            SquareHook(const KinematicVariable & variable) :
                _variable(variable)
            {
                uses_kinematic(variable.id());
            }

            virtual double
            evaluate() const
            {
                return _variable * _variable;
            }

            virtual std::shared_ptr<const KinematicHook>
            clone(const Kinematics & kinematics) const
            {
                return std::make_shared<SquareHook>(kinematics[_variable.name()]);
            }
    };
} // namespace

class KinematicsTest : public TestCase
{
    public:
        KinematicsTest() :
            TestCase("kinematics_test")
        {
        }

        virtual void
        run() const
        {
            // Creation from initializer list
            {
                Kinematics kinematics{
                    { "s_min", 1.0 },
                    { "s_max", 6.0 },
                };

                TEST_CHECK_EQUAL(1.0, kinematics["s_min"]);
                TEST_CHECK_EQUAL(6.0, kinematics["s_max"]);
            }

            // Creation from initializer list and aliasing
            {
                Kinematics kinematics{
                    { "s_min", 1.0 },
                    { "s_max", 6.0 },
                };
                kinematics.alias("q2_min", "s_min");
                kinematics.alias("q2_max", "s_max");

                TEST_CHECK_EQUAL(1.0, kinematics["s_min"]);
                TEST_CHECK_EQUAL(6.0, kinematics["s_max"]);
                TEST_CHECK_EQUAL(1.0, kinematics["q2_min"]);
                TEST_CHECK_EQUAL(6.0, kinematics["q2_max"]);

                TEST_CHECK_EQUAL("s_min", kinematics["s_min"].name());
                TEST_CHECK_EQUAL("s_max", kinematics["s_max"].name());

                // Test clearing the alias map and access values
                auto aliased = kinematics["s_min"];
                auto alias   = kinematics["q2_min"];
                kinematics.clear_aliases();
                TEST_CHECK_EQUAL(1.0, aliased);
                TEST_CHECK_EQUAL(1.0, alias);
            }

            // Aliases take precedence over variables of the same name
            {
                Kinematics kinematics{
                    {   "q2", 1.0 },
                    { "q2_e", 2.0 },
                };
                kinematics.alias("q2", "q2_e");

                TEST_CHECK_EQUAL(2.0, kinematics["q2"].evaluate());
                TEST_CHECK_EQUAL("q2_e", kinematics["q2"].name());

                kinematics.set("q2", 3.0);
                TEST_CHECK_EQUAL(3.0, kinematics["q2_e"].evaluate());

                kinematics.declare("q2", 4.0);
                TEST_CHECK_EQUAL(4.0, kinematics["q2_e"].evaluate());

                TEST_CHECK_THROWS(UnknownKinematicVariableError, kinematics.set("q2_mu", 5.0));

                kinematics.remove_alias("q2");
                TEST_CHECK_EQUAL(1.0, kinematics["q2"].evaluate());
                TEST_CHECK_EQUAL(4.0, kinematics["q2_e"].evaluate());
            }

            // Testing has(), for variables as well as for aliases
            {
                Kinematics kinematics{
                    { "s_min", 1.0 },
                    { "s_max", 6.0 },
                };
                kinematics.alias("q2_min", "s_min");

                TEST_CHECK_EQUAL(kinematics.has("s_min"), true);
                TEST_CHECK_EQUAL(kinematics.has("q2_min"), true);
                TEST_CHECK_EQUAL(kinematics.has("q2_max"), false);

                // an alias is no longer found once the alias map has been cleared
                kinematics.clear_aliases();

                TEST_CHECK_EQUAL(kinematics.has("s_min"), true);
                TEST_CHECK_EQUAL(kinematics.has("q2_min"), false);
                TEST_CHECK_THROWS(UnknownKinematicVariableError, kinematics["q2_min"]);
            }

            // Access
            {
                Kinematics kinematics;
                kinematics.declare("foo");
                kinematics.set("foo", 17.0);

                TEST_CHECK_EQUAL(17.0, kinematics["foo"]);
            }

            // Equality/Inequality
            {
                Kinematics a, b;

                // check self equality
                TEST_CHECK(a == a);
                TEST_CHECK(b == b);

                // check empty kinematics
                TEST_CHECK(a == b);

                // populate a
                a.declare("foo");
                a.set("foo", 19.0);
                TEST_CHECK(! (a == b));
                TEST_CHECK(a != b);

                // populate b false
                b.declare("foo");
                b.set("foo", 21.3);
                TEST_CHECK(! (a == b));

                // populate b correctly
                b.set("foo", 19.0);
                TEST_CHECK(a == b);

                // add alias to a, but not b
                a.alias("baz", "foo");
                TEST_CHECK(! (a == b));

                // add alias to b, too
                b.alias("baz", "foo");
                TEST_CHECK(a == b);

                // copy
                Kinematics c = a;
                TEST_CHECK(a == c);
                TEST_CHECK(b == c);
            }

            // Iteration (check for names, values, and order)
            {
                Kinematics k{
                    {      "s_min",  1.0 },
                    {      "s_max",  6.0 },
                    { "cos(theta)", -0.5 },
                };

                auto i = k.begin();
                TEST_CHECK("s_min" == i->name());
                TEST_CHECK(1.0 == i->evaluate());

                ++i;
                TEST_CHECK("s_max" == i->name());
                TEST_CHECK(6.0 == i->evaluate());

                ++i;
                TEST_CHECK("cos(theta)" == i->name());
                TEST_CHECK(-0.5 == i->evaluate());

                ++i;
                TEST_CHECK(k.end() == i);
            }

            // Iteration (check for names, values, and order, in presence of an alias)
            {
                Kinematics k{
                    {      "s_min",  1.0 },
                    {      "s_max",  6.0 },
                    { "cos(theta)", -0.5 },
                };
                k.alias("z", "cos(theta)");

                auto i = k.begin();
                TEST_CHECK("s_min" == i->name());
                TEST_CHECK(1.0 == i->evaluate());

                ++i;
                TEST_CHECK("s_max" == i->name());
                TEST_CHECK(6.0 == i->evaluate());

                ++i;
                TEST_CHECK("cos(theta)" == i->name());
                TEST_CHECK(-0.5 == i->evaluate());

                ++i;
                TEST_CHECK(k.end() == i);

                TEST_CHECK_NO_THROW(-0.5 == k["z"].evaluate());
            }

            // Derived variables
            {
                Kinematics k{
                    { "m", 2.0 }
                };
                auto q2 = k.declare_derived("q2", std::make_shared<SquareHook>(k["m"]));

                TEST_CHECK_EQUAL(4.0, q2.evaluate());
                TEST_CHECK_EQUAL(4.0, static_cast<double>(k["q2"]));
                TEST_CHECK_EQUAL(4.0, q2());

                k.set("m", 3.0);
                TEST_CHECK_EQUAL(9.0, q2.evaluate());

                // read through an alias
                k.alias("s", "q2");
                TEST_CHECK_EQUAL(9.0, k["s"].evaluate());
                k.remove_alias("s");

                // derived variables cannot be set or declared anew
                TEST_CHECK_THROWS(InternalError, k.set("q2", 1.0));
                TEST_CHECK_THROWS(InternalError, q2.set(1.0));
                TEST_CHECK_THROWS(InternalError, k.declare("q2", 1.0));
                TEST_CHECK_THROWS(InternalError, k.declare_derived("q2", std::make_shared<SquareHook>(k["m"])));
                TEST_CHECK_THROWS(InternalError, k.declare_derived("m", std::make_shared<SquareHook>(k["m"])));

                // a clone reads its own variables
                Kinematics clone = k.clone();
                clone.set("m", 4.0);
                TEST_CHECK_EQUAL(16.0, clone["q2"].evaluate());
                TEST_CHECK_EQUAL(9.0, k["q2"].evaluate());

                // comparison and printing use the derived values
                TEST_CHECK(k != clone);
                clone.set("m", 3.0);
                TEST_CHECK(k == clone);
                TEST_CHECK_EQUAL_STR("m=3, q2=9", k.as_string());

                // the sum reads its own variables, too
                Kinematics sum = Kinematics{
                    { "E", 1.0 }
                } + k;
                sum.set("m", 5.0);
                TEST_CHECK_EQUAL(25.0, sum["q2"].evaluate());
                TEST_CHECK_EQUAL(9.0, k["q2"].evaluate());

                // a derived variable must not share its name with a variable of the other summand
                const Kinematics other{
                    { "q2", 1.0 }
                };
                TEST_CHECK_THROWS(InternalError, k + other);
                TEST_CHECK_THROWS(InternalError, other + k);

                // the hook reports the variables that it reads
                TEST_CHECK(nullptr == k.hook("m"));
                TEST_CHECK(nullptr != k.hook("q2"));
                TEST_CHECK_EQUAL(1, std::distance(k.hook("q2")->begin_kinematics(), k.hook("q2")->end_kinematics()));
                TEST_CHECK_EQUAL(k["m"].id(), *k.hook("q2")->begin_kinematics());
            }
        }
} kinematics_test;
