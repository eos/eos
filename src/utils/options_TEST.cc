/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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
#include <src/utils/options.hh>

#include <cmath>

using namespace test;
using namespace eos;

class OptionsTest :
    public TestCase
{
    public:
        OptionsTest() :
            TestCase("options_test")
        {
        }

        virtual void run() const
        {
            // Access
            {
                Options options;
                options.set("foo", "bar");

                TEST_CHECK(options.has("foo"));
                TEST_CHECK(! options.has("baz"));
                TEST_CHECK_EQUAL("bar", options.get("foo", ""));
            }

            // Equality/Inequality
            {
                Options a, b;

                // check self equality
                TEST_CHECK(a == a);
                TEST_CHECK(b == b);

                // check empty options
                TEST_CHECK(a == b);

                // populate a
                a.set("foo", "bar");
                TEST_CHECK(! (a == b));
                TEST_CHECK(a != b);

                // populate b false
                b.set("foo", "baz");
                TEST_CHECK(! (a == b));

                // populate b correctly
                b.set("foo", "bar");
                TEST_CHECK(a == b);

                // copy
                Options c = a;
                TEST_CHECK(a == c);
                TEST_CHECK(b == c);
            }
        }
} options_test;
