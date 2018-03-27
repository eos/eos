/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2016, 2018 Danny van Dyk
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
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>

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
            // Creation from initializer list
            {
                Options options
                {
                    { "q", "d"  },
                    { "l", "mu" },
                };

                TEST_CHECK_EQUAL("d",  options.get("q", "s"));
                TEST_CHECK_EQUAL("mu", options.get("l", "tau"));
            }

            // Access
            {
                Options options;
                options.set("foo", "bar");

                TEST_CHECK(options.has("foo"));
                TEST_CHECK(! options.has("baz"));
                TEST_CHECK_EQUAL("bar", options.get("foo", ""));
            }

            // Merging w/o duplicates
            {
                Options o1
                {
                    { "q", "d" }
                };

                Options o2
                {
                    { "l", "mu" }
                };

                Options o3 = o1 + o2;

                TEST_CHECK_EQUAL("d",  o3.get("q", "s"));
                TEST_CHECK_EQUAL("mu", o3.get("l", "tau"));
            }

            // Merging w/ duplicates
            {
                Options o1
                {
                    { "q", "d" }
                };

                Options o2
                {
                    { "q", "s" },
                    { "l", "mu" }
                };

                Options o3 = o1 + o2;

                TEST_CHECK_EQUAL("s",  o3.get("q", "foo"));
                TEST_CHECK_EQUAL("mu", o3.get("l", "tau"));
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

            // Iteration (check for names, values, and lexicographical order)
            {
                Options o
                {
                    { "l",     "tau"     },
                    { "q",     "u"       },
                    { "model", "CKMScan" },
                };

                auto i = o.begin();
                TEST_CHECK("l" == i->first);
                TEST_CHECK("tau" == i->second);

                ++i;
                TEST_CHECK("model" == i->first);
                TEST_CHECK("CKMScan" == i->second);

                ++i;
                TEST_CHECK("q" == i->first);
                TEST_CHECK("u" == i->second);

                ++i;
                TEST_CHECK(o.end() == i);
            }
        }
} options_test;

class SwitchOptionTest :
    public TestCase
{
    public:
        SwitchOptionTest() :
            TestCase("switch_option_test")
        {
        }

        virtual void run() const
        {
            // Creation with valid default value, value = default value, non-empty list
            {
                SwitchOption so
                {
                    Options{ { "key", "value1" }, { "unused", "foo" } },
                    "key",
                    { "value1", "value2", "value4" },
                    "value1"
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation with unspecified value, non-empty list
            {
                SwitchOption so
                {
                    Options{ { "unused", "foo" } },
                    "key",
                    { "value1", "value2", "value4" },
                    "value1"
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation with non-empty list, no default value
            {
                SwitchOption so
                {
                    Options{ { "key", "value4" }, { "unused", "foo" } },
                    "key",
                    { "value1", "value2", "value4" }
                };
                TEST_CHECK_EQUAL(so.value(), "value4");
            }

            // Creation with unspecified value, non-empty list, no default value
            {
                auto test = [] ()
                {
                    SwitchOption so
                    {
                        Options{ { "unused", "foo" } },
                        "key",
                        { "value1", "value2", "value4" }
                    };
                };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }

            // Creation with value not in list of allowed values
            {
                auto test = [] ()
                {
                    SwitchOption so
                    {
                        Options{ { "key", "value3"}, { "unused", "foo" } },
                        "key",
                        { "value1", "value2", "value4" },
                        "value1"
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creating with empty list
            {
                auto test = [] ()
                {
                    SwitchOption so
                    {
                        Options{ { "key", "value1" }, { "unused", "foo" } },
                        "key",
                        { },
                        "value1"
                    };
                };
                TEST_CHECK_THROWS(InternalError, test());
            }

            // Creating with empty list
            {
                auto test = [] ()
                {
                    SwitchOption so
                    {
                        Options{ { "key", "value1" }, { "unused", "foo" } },
                        "key",
                        { },
                        "value1"
                    };
                };
                TEST_CHECK_THROWS(InternalError, test());
            }
        }
} switch_option_test;
