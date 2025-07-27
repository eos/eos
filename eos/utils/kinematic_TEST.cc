/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2018 Danny van Dyk
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

using namespace test;
using namespace eos;

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
        }
} kinematics_test;
