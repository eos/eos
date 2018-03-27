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

#include <test/test.hh>
#include <eos/utils/kinematic.hh>

#include <cmath>

using namespace test;
using namespace eos;

class KinematicsTest :
    public TestCase
{
    public:
        KinematicsTest() :
            TestCase("kinematics_test")
        {
        }

        virtual void run() const
        {
            // Creation from initializer list
            {
                Kinematics kinematics
                {
                    { "s_min", 1.0 },
                    { "s_max", 6.0 },
                };

                TEST_CHECK_EQUAL(1.0, kinematics["s_min"]);
                TEST_CHECK_EQUAL(6.0, kinematics["s_max"]);
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

                // copy
                Kinematics c = a;
                TEST_CHECK(a == c);
                TEST_CHECK(b == c);
            }

            // Iteration (check for names, values, and order)
            {
                Kinematics k
                {
                    { "s_min",       1.0 },
                    { "s_max",       6.0 },
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
        }
} kinematics_test;
