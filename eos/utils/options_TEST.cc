/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
                    { "q"_ok, "d"  },
                    { "l"_ok, "mu" },
                };

                TEST_CHECK_EQUAL("d",  options.get("q"_ok, "s"));
                TEST_CHECK_EQUAL("mu", options.get("l"_ok, "tau"));
            }

            // Access
            {
                Options options;
                options.declare("foo"_ok, "bar");

                TEST_CHECK(options.has("foo"_ok));
                TEST_CHECK(! options.has("baz"_ok));
                TEST_CHECK_EQUAL("bar", options.get("foo"_ok, ""));
            }

            // Merging w/o duplicates
            {
                Options o1
                {
                    { "q"_ok, "d" }
                };

                Options o2
                {
                    { "l"_ok, "mu" }
                };

                Options o3 = o1 + o2;

                TEST_CHECK_EQUAL("d"_ok,  o3.get("q"_ok, "s"));
                TEST_CHECK_EQUAL("mu"_ok, o3.get("l"_ok, "tau"));
            }

            // Merging w/ duplicates
            {
                Options o1
                {
                    { "q"_ok, "d" }
                };

                Options o2
                {
                    { "q"_ok, "s" },
                    { "l"_ok, "mu" }
                };

                Options o3 = o1 + o2;

                TEST_CHECK_EQUAL("s"_ok,  o3.get("q"_ok, "foo"));
                TEST_CHECK_EQUAL("mu"_ok, o3.get("l"_ok, "tau"));
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
                a.declare("foo"_ok, "bar");
                TEST_CHECK(! (a == b));
                TEST_CHECK(a != b);

                // populate b false
                b.declare("foo"_ok, "baz");
                TEST_CHECK(! (a == b));

                // populate b correctly
                b.declare("foo"_ok, "bar");
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
                    { "l"_ok,     "tau"     },
                    { "q"_ok,     "u"       },
                    { "model"_ok, "CKM" },
                };

                auto i = o.begin();
                TEST_CHECK("l"_ok == i->first);
                TEST_CHECK("tau" == i->second);

                ++i;
                TEST_CHECK("model"_ok == i->first);
                TEST_CHECK("CKM" == i->second);

                ++i;
                TEST_CHECK("q"_ok == i->first);
                TEST_CHECK("u" == i->second);

                ++i;
                TEST_CHECK(o.end() == i);
            }
        }
} options_test;

class NameOptionTest :
    public TestCase
{
    public:
        NameOptionTest() :
            TestCase("name_option_test")
        {
        }

        virtual void run() const
        {
            // Creation with valid default value, value = default value
            {
                NameOption no
                {
                    Options{ { "key"_ok, "value1" }, { "unused"_ok, "foo" } },
                    "key"_ok,
                    qnp::Name("value1")
                };
                TEST_CHECK_EQUAL(no.value(), qnp::Name("value1"));
            }

            // Creation with unspecified value
            {
                NameOption no
                {
                    Options{ { "unused"_ok, "foo" } },
                    "key"_ok,
                    qnp::Name("value1")
                };
                TEST_CHECK_EQUAL(no.value(), qnp::Name("value1"));
            }

            // Creation without default value
            {
                NameOption no
                {
                    Options{ { "key"_ok, "value4" }, { "unused"_ok, "foo" } },
                    "key"_ok,
                };
                TEST_CHECK_EQUAL(no.value(), qnp::Name("value4"));
            }

            // Creation with unspecified value, non-empty list, no default value
            {
                auto test = [] ()
                {
                    NameOption no
                    {
                        Options{ { "unused"_ok, "foo" } },
                        "key"_ok
                    };
                };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }

            // Creation with invalid value
            {
                auto test = [] ()
                {
                    NameOption so
                    {
                        Options{ { "key"_ok, "invalid@value"}, { "unused"_ok, "foo" } },
                        "key"_ok,
                        qnp::Name("value1")
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} name_option_test;

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
                    Options{ { "key"_ok, "value1" }, { "unused"_ok, "foo" } },
                    "key"_ok,
                    { "value1", "value2", "value4" },
                    "value1"
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation with unspecified value, non-empty list
            {
                SwitchOption so
                {
                    Options{ { "unused"_ok, "foo" } },
                    "key"_ok,
                    { "value1", "value2", "value4" },
                    "value1"
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation with non-empty list, no default value
            {
                SwitchOption so
                {
                    Options{ { "key"_ok, "value4" }, { "unused"_ok, "foo" } },
                    "key"_ok,
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
                        Options{ { "unused"_ok, "foo" } },
                        "key"_ok,
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
                        Options{ { "key"_ok, "value3"}, { "unused"_ok, "foo" } },
                        "key"_ok,
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
                        Options{ { "key"_ok, "value1" }, { "unused"_ok, "foo" } },
                        "key"_ok,
                        { },
                        "value1"
                    };
                };
                TEST_CHECK_THROWS(InternalError, test());
            }
        }
} switch_option_test;

class SpecifiedOptionTest :
    public TestCase
{
    public:
        SpecifiedOptionTest() :
            TestCase("specified_option_test")
        {
        }

        virtual void run() const
        {
            // specify permitted options
            std::vector<OptionSpecification> specifications
            {
                { "key-with-full-specification"_ok,   { "value1", "value2", "value4" }, "value1" },
                { "key-without-default"_ok,           { "value1", "value2", "value4" }           }
            };

            // Creation of fully specified option, with value
            {
                SpecifiedOption so
                {
                    Options{
                        { "key-with-full-specification"_ok, "value1" }, { "unused"_ok, "foo" }
                    },
                    specifications,
                    "key-with-full-specification"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation of fully specified option, without value
            {
                SpecifiedOption so
                {
                    Options{
                        { "unused"_ok, "foo" }
                    },
                    specifications,
                    "key-with-full-specification"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation of specified option without default value, with value
            {
                SpecifiedOption so
                {
                    Options{
                        { "key-without-default"_ok, "value4" }, { "unused"_ok, "foo" }
                    },
                    specifications,
                    "key-without-default"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value4");
            }
            // Creation of specified option without default value, without value
            {
                auto test = [specifications] ()
                {
                    SpecifiedOption so
                    {
                        Options{
                            { "unused"_ok, "foo" }
                        },
                        specifications,
                        "key-without-default"_ok
                    };
                };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }
        }
} specified_option_test;

class RestrictedOptionTest :
    public TestCase
{
    public:
        RestrictedOptionTest() :
            TestCase("restricted_option_test")
        {
        }

        virtual void run() const
        {
            // specify permitted options
            std::vector<OptionSpecification> specifications
            {
                { "key-with-default"_ok,              { "value1", "value2", "value4" }, "value1" },
                { "key-without-value"_ok,             { "foo", "bar" },                 "foo"    },
                { "key-without-default"_ok,           { "value1", "value2", "value4" }           },
                { "key-without-value-and-default"_ok, { "foo", "bar" }                           },
                { "key-with-invalid-default"_ok,      { "value1", "value2", "value4" }, "value3" },
                { "key-with-invalid-value"_ok,        { "value1", "value2", "value4" }, "value1" },
            };

            // Creation of option with valid default value, with value = default value
            {
                RestrictedOption so
                {
                    Options{
                        { "key-with-default"_ok, "value1" }, { "unused"_ok, "foo" }
                    },
                    specifications,
                    "key-with-default"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation of option with valid default value, with unspecified value
            {
                RestrictedOption so
                {
                    Options{
                        { "unused"_ok, "foo" }
                    },
                    specifications,
                    "key-without-value"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "foo");
            }

            // Creation of option without default valuee, with valid value
            {
                RestrictedOption so
                {
                    Options{
                        { "key-without-default"_ok, "value4" }, { "unused"_ok, "foo" }
                    },
                    specifications,
                    "key-without-default"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value4");
            }

            // Creation of option without default value, with unspecified valu
            {
                auto test = [specifications] ()
                {
                    RestrictedOption so
                    {
                        Options{
                            { "unused"_ok, "foo" }
                        },
                        specifications,
                        "key-without-value-and-default"_ok
                    };
                };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }

            // Creation with invalid value
            {
                auto test = [specifications] ()
                {
                    RestrictedOption so
                    {
                        Options{
                            { "unused"_ok, "foo" }
                        },
                        specifications,
                        "key-with-invalid-default"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creation with invalid value
            {
                auto test = [specifications] ()
                {
                    RestrictedOption so
                    {
                        Options{
                            { "key-with-invalid-value"_ok, "value3"}, { "unused"_ok, "foo" }
                        },
                        specifications,
                        "key-with-invalid-value"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} restricted_option_test;

class IsospinOptionTest :
    public TestCase
{
    public:
        IsospinOptionTest() :
            TestCase("isospin_option_test")
        {
        }

        virtual void run() const
        {
            // specify permitted options
            std::vector<OptionSpecification> specifications
            {
                { "I"_ok,              { "0|2" }, "0" },
            };

            // Creation of option with valid default value, with value = default value
            {
                IsospinOption so
                {
                    Options{
                        { "I"_ok, "0" }
                    },
                    specifications,
                    "I"_ok
                };
                TEST_CHECK_EQUAL(so.value(), Isospin::zero);
            }

            // Creation of option with valid default value, with value != default value
            {
                IsospinOption so
                {
                    Options{
                        { "I"_ok, "2" }
                    },
                    specifications,
                    "I"_ok
                };
                TEST_CHECK_EQUAL(so.value(), Isospin::two);
            }

            // Creation of option with valid default value, with value != default value but containing default value
            {
                IsospinOption so
                {
                    Options{
                        { "I"_ok, "0|2" }
                    },
                    specifications,
                    "I"_ok
                };
                TEST_CHECK_EQUAL(so.value(), Isospin::zero | Isospin::two);
            }

            // Creation with invalid value
            {
                auto test = [specifications] ()
                {
                    IsospinOption so
                    {
                        Options{
                            { "I"_ok, "1"}
                        },
                        specifications,
                        "I"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creation with invalid value but containing a valid value
            {
                auto test = [specifications] ()
                {
                    IsospinOption so
                    {
                        Options{
                            { "I"_ok, "0|1"}
                        },
                        specifications,
                        "I"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} isospin_option_test;

class PartialWaveOptionTest :
    public TestCase
{
    public:
        PartialWaveOptionTest() :
            TestCase("partial_wave_option_test")
        {
        }

        virtual void run() const
        {
            // specify permitted options
            std::vector<OptionSpecification> specifications
            {
                { "L"_ok,              { "S|P|D" }, "S|P" },
            };

            // Creation of option with valid default value, with value = default value
            {
                PartialWaveOption so
                {
                    Options{
                        { "L"_ok, "S|P" }
                    },
                    specifications,
                    "L"_ok
                };
                TEST_CHECK_EQUAL(so.value(), PartialWave::S | PartialWave::P);
            }

            // Creation of option with valid default value, with value != default value but part of default value
            {
                PartialWaveOption so
                {
                    Options{
                        { "L"_ok, "P" }
                    },
                    specifications,
                    "L"_ok
                };
                TEST_CHECK_EQUAL(so.value(), PartialWave::P);
            }

            // Creation of option with valid default value, with value != default value but containing default value
            {
                PartialWaveOption so
                {
                    Options{
                        { "L"_ok, "S|P|D" }
                    },
                    specifications,
                    "L"_ok
                };
                TEST_CHECK_EQUAL(so.value(), PartialWave::S | PartialWave::P | PartialWave::D);
            }

            // Creation with invalid value
            {
                auto test = [specifications] ()
                {
                    PartialWaveOption so
                    {
                        Options{
                            { "L"_ok, "F"}
                        },
                        specifications,
                        "L"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creation with invalid value but containing a valid value
            {
                auto test = [specifications] ()
                {
                    PartialWaveOption so
                    {
                        Options{
                            { "L"_ok, "P|F"}
                        },
                        specifications,
                        "L"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} partial_wave_option_test;
