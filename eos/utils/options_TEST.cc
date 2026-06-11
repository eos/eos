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

#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>

#include <test/test.hh>

#include <cmath>

using namespace test;
using namespace eos;
using namespace std::literals::string_literals;

class OptionsTest : public TestCase
{
    public:
        OptionsTest() :
            TestCase("options_test")
        {
        }

        virtual void
        run() const
        {
            // Creation from initializer list
            {
                Options options{
                    { "q"_ok,  "d"_ov },
                    { "l"_ok, "mu"_ov },
                };

                TEST_CHECK_EQUAL("d"_ov, options.get("q"_ok, "s"_ov));
                TEST_CHECK_EQUAL("mu"_ov, options.get("l"_ok, "tau"_ov));
            }

            // Access
            {
                Options options;
                options.declare("foo"_ok, "bar"_ov);

                TEST_CHECK(options.has("foo"_ok));
                TEST_CHECK(! options.has("baz"_ok));
                TEST_CHECK_EQUAL("bar"_ov, options.get("foo"_ok, "baz"_ov));
            }

            // Merging w/o duplicates
            {
                Options o1{
                    { "q"_ok, "d"_ov }
                };

                Options o2{
                    { "l"_ok, "mu"_ov }
                };

                Options o3 = o1 + o2;

                TEST_CHECK_EQUAL("d"_ov, o3.get("q"_ok, "s"_ov));
                TEST_CHECK_EQUAL("mu"_ov, o3.get("l"_ok, "tau"_ov));
            }

            // Merging w/ duplicates
            {
                Options o1{
                    { "q"_ok, "d"_ov }
                };

                Options o2{
                    { "q"_ok,  "s"_ov },
                    { "l"_ok, "mu"_ov }
                };

                Options o3 = o1 + o2;

                TEST_CHECK_EQUAL("s"_ov, o3.get("q"_ok, "foo"_ov));
                TEST_CHECK_EQUAL("mu"_ov, o3.get("l"_ok, "tau"_ov));
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
                a.declare("foo"_ok, "bar"_ov);
                TEST_CHECK(! (a == b));
                TEST_CHECK(a != b);

                // populate b false
                b.declare("foo"_ok, "baz"_ov);
                TEST_CHECK(! (a == b));

                // populate b correctly
                b.declare("foo"_ok, "bar"_ov);
                TEST_CHECK(a == b);

                // copy
                Options c = a;
                TEST_CHECK(a == c);
                TEST_CHECK(b == c);
            }

            // Iteration (check for names, values, and lexicographical order)
            {
                Options o{
                    {     "l"_ok, "tau"_ov },
                    {     "q"_ok,   "u"_ov },
                    { "model"_ok, "CKM"_ov },
                };

                auto i = o.begin();
                TEST_CHECK("l"_ok == i->first);
                TEST_CHECK("tau"_ov == i->second);

                ++i;
                TEST_CHECK("model"_ok == i->first);
                TEST_CHECK("CKM"_ov == i->second);

                ++i;
                TEST_CHECK("q"_ok == i->first);
                TEST_CHECK("u"_ov == i->second);

                ++i;
                TEST_CHECK(o.end() == i);
            }
        }
} options_test;

class NameOptionTest : public TestCase
{
    public:
        NameOptionTest() :
            TestCase("name_option_test")
        {
        }

        virtual void
        run() const
        {
            // Creation with valid default value, value = default value
            {
                NameOption no{
                    Options{ { "key"_ok, "value1"_ov }, { "unused"_ok, "foo"_ov } },
                    "key"_ok,
                    qnp::Name("value1")
                };
                TEST_CHECK_EQUAL(no.value(), qnp::Name("value1"));
            }

            // Creation with unspecified value
            {
                NameOption no{ Options{ { "unused"_ok, "foo"_ov } }, "key"_ok, qnp::Name("value1") };
                TEST_CHECK_EQUAL(no.value(), qnp::Name("value1"));
            }

            // Creation without default value
            {
                NameOption no{
                    Options{ { "key"_ok, "value4"_ov }, { "unused"_ok, "foo"_ov } },
                    "key"_ok,
                };
                TEST_CHECK_EQUAL(no.value(), qnp::Name("value4"));
            }

            // Creation with unspecified value, non-empty list, no default value
            {
                auto test = []() { NameOption no{ Options{ { "unused"_ok, "foo"_ov } }, "key"_ok }; };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }
        }
} name_option_test;

class SwitchOptionTest : public TestCase
{
    public:
        SwitchOptionTest() :
            TestCase("switch_option_test")
        {
        }

        virtual void
        run() const
        {
            // Creation with valid default value, value = default value, non-empty list
            {
                SwitchOption so{
                    Options{ { "key"_ok, "value1"_ov }, { "unused"_ok, "foo"_ov } },
                    "key"_ok,
                    { "value1"_ov, "value2"_ov, "value4"_ov },
                    "value1"_ov
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation with unspecified value, non-empty list
            {
                SwitchOption so{
                    Options{ { "unused"_ok, "foo"_ov } },
                    "key"_ok,
                    { "value1"_ov, "value2"_ov, "value4"_ov },
                    "value1"_ov
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation with non-empty list, no default value
            {
                SwitchOption so{
                    Options{ { "key"_ok, "value4"_ov }, { "unused"_ok, "foo"_ov } },
                    "key"_ok,
                    { "value1"_ov, "value2"_ov, "value4"_ov }
                };
                TEST_CHECK_EQUAL(so.value(), "value4");
            }

            // Creation with unspecified value, non-empty list, no default value
            {
                auto test = []()
                {
                    SwitchOption so{
                        Options{ { "unused"_ok, "foo"_ov } },
                        "key"_ok,
                        { "value1"_ov, "value2"_ov, "value4"_ov }
                    };
                };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }

            // Creation with value not in list of allowed values
            {
                auto test = []()
                {
                    SwitchOption so{
                        Options{ { "key"_ok, "value3"_ov }, { "unused"_ok, "foo"_ov } },
                        "key"_ok,
                        { "value1"_ov, "value2"_ov, "value4"_ov },
                        "value1"_ov
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creating with empty list
            {
                auto test = []()
                {
                    SwitchOption so{
                        Options{ { "key"_ok, "value1"_ov }, { "unused"_ok, "foo"_ov } },
                        "key"_ok,
                        {},
                        "value1"_ov
                    };
                };
                TEST_CHECK_THROWS(InternalError, test());
            }
        }
} switch_option_test;

class SpecifiedOptionTest : public TestCase
{
    public:
        SpecifiedOptionTest() :
            TestCase("specified_option_test")
        {
        }

        virtual void
        run() const
        {
            // specify permitted options
            std::vector<OptionSpecification> specifications{
                { "key-with-full-specification"_ok, { "value1"_ov, "value2"_ov, "value4"_ov }, "value1"_ov },
                { "key-without-default"_ok, { "value1"_ov, "value2"_ov, "value4"_ov } }
            };

            // Creation of fully specified option, with value
            {
                SpecifiedOption so{
                    Options{ { "key-with-full-specification"_ok, "value1"_ov }, { "unused"_ok, "foo"_ov } },
                    specifications,
                    "key-with-full-specification"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation of fully specified option, without value
            {
                SpecifiedOption so{ Options{ { "unused"_ok, "foo"_ov } }, specifications, "key-with-full-specification"_ok };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation of specified option without default value, with value
            {
                SpecifiedOption so{
                    Options{ { "key-without-default"_ok, "value4"_ov }, { "unused"_ok, "foo"_ov } },
                    specifications,
                    "key-without-default"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value4");
            }
            // Creation of specified option without default value, without value
            {
                auto test = [specifications]() { SpecifiedOption so{ Options{ { "unused"_ok, "foo"_ov } }, specifications, "key-without-default"_ok }; };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }
        }
} specified_option_test;

class RestrictedOptionTest : public TestCase
{
    public:
        RestrictedOptionTest() :
            TestCase("restricted_option_test")
        {
        }

        virtual void
        run() const
        {
            // specify permitted options
            std::vector<OptionSpecification> specifications{
                { "key-with-default"_ok, { "value1"_ov, "value2"_ov, "value4"_ov }, "value1"_ov },
                { "key-without-value"_ok, { "foo"_ov, "bar"_ov }, "foo"_ov },
                { "key-without-default"_ok, { "value1"_ov, "value2"_ov, "value4"_ov } },
                { "key-without-value-and-default"_ok, { "foo"_ov, "bar"_ov } },
                { "key-with-invalid-default"_ok, { "value1"_ov, "value2"_ov, "value4"_ov }, "value3"_ov },
                { "key-with-invalid-value"_ok, { "value1"_ov, "value2"_ov, "value4"_ov }, "value1"_ov },
            };

            // Creation of option with valid default value, with value = default value
            {
                RestrictedOption so{
                    Options{ { "key-with-default"_ok, "value1"_ov }, { "unused"_ok, "foo"_ov } },
                    specifications,
                    "key-with-default"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value1");
            }

            // Creation of option with valid default value, with unspecified value
            {
                RestrictedOption so{ Options{ { "unused"_ok, "foo"_ov } }, specifications, "key-without-value"_ok };
                TEST_CHECK_EQUAL(so.value(), "foo");
            }

            // Creation of option without default valuee, with valid value
            {
                RestrictedOption so{
                    Options{ { "key-without-default"_ok, "value4"_ov }, { "unused"_ok, "foo"_ov } },
                    specifications,
                    "key-without-default"_ok
                };
                TEST_CHECK_EQUAL(so.value(), "value4");
            }

            // Creation of option without default value, with unspecified valu
            {
                auto test = [specifications]() { RestrictedOption so{ Options{ { "unused"_ok, "foo"_ov } }, specifications, "key-without-value-and-default"_ok }; };
                TEST_CHECK_THROWS(UnspecifiedOptionError, test());
            }

            // Creation with invalid value
            {
                auto test = [specifications]() { RestrictedOption so{ Options{ { "unused"_ok, "foo"_ov } }, specifications, "key-with-invalid-default"_ok }; };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creation with invalid value
            {
                auto test = [specifications]()
                {
                    RestrictedOption so{
                        Options{ { "key-with-invalid-value"_ok, "value3"_ov }, { "unused"_ok, "foo"_ov } },
                        specifications,
                        "key-with-invalid-value"_ok
                    };
                };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} restricted_option_test;

class IsospinOptionTest : public TestCase
{
    public:
        IsospinOptionTest() :
            TestCase("isospin_option_test")
        {
        }

        virtual void
        run() const
        {
            // correctly specify permitted options
            std::vector<OptionSpecification> specifications{
                { "I"_ok, { "0|2"_ov }, "0"_ov },
            };

            // incorrectly specify permitted options
            std::vector<OptionSpecification> invalid_specifications{
                { "I"_ok, std::vector<qnp::OptionValue>{ "0|2"_ov }, "0"_ov },
            };

            // Creation of option with invalid specification
            {
                auto test = [invalid_specifications]() { IsospinOption so{ Options{ { "I"_ok, "0"_ov } }, invalid_specifications, "I"_ok }; };
                TEST_CHECK_THROWS(InternalError, test());
            }

            // Creation of option with valid default value, with value = default value
            {
                IsospinOption so{ Options{ { "I"_ok, "0"_ov } }, specifications, "I"_ok };
                TEST_CHECK_EQUAL(so.value(), Isospin::zero);
            }

            // Creation of option with valid default value, with value != default value
            {
                IsospinOption so{ Options{ { "I"_ok, "2"_ov } }, specifications, "I"_ok };
                TEST_CHECK_EQUAL(so.value(), Isospin::two);
            }

            // Creation of option with valid default value, with value != default value but containing default value
            {
                IsospinOption so{ Options{ { "I"_ok, "0|2"_ov } }, specifications, "I"_ok };
                TEST_CHECK_EQUAL(so.value(), Isospin::zero | Isospin::two);
            }

            // Creation with invalid value
            {
                auto test = [specifications]() { IsospinOption so{ Options{ { "I"_ok, "1"_ov } }, specifications, "I"_ok }; };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creation with invalid value but containing a valid value
            {
                auto test = [specifications]() { IsospinOption so{ Options{ { "I"_ok, "0|1"_ov } }, specifications, "I"_ok }; };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} isospin_option_test;

class PartialWaveOptionTest : public TestCase
{
    public:
        PartialWaveOptionTest() :
            TestCase("partial_wave_option_test")
        {
        }

        virtual void
        run() const
        {
            // correctly specify permitted options
            std::vector<OptionSpecification> specifications{
                { "L"_ok, "S|P|D"_ov, "S|P"_ov },
            };

            // incorrectly specify permitted options
            std::vector<OptionSpecification> invalid_specifications{
                { "L"_ok, std::vector<qnp::OptionValue>{ "S|P|D"_ov }, "S|P"_ov },
            };

            // Creation of option with invalid specification
            {
                auto test = [invalid_specifications]() { PartialWaveOption so{ Options{ { "L"_ok, "S|P"_ov } }, invalid_specifications, "L"_ok }; };
                TEST_CHECK_THROWS(InternalError, test());
            }

            // Creation of option with valid default value, with value = default value
            {
                PartialWaveOption so{ Options{ { "L"_ok, "S|P"_ov } }, specifications, "L"_ok };
                TEST_CHECK_EQUAL(so.value(), PartialWave::S | PartialWave::P);
            }

            // Creation of option with valid default value, with value != default value but part of default value
            {
                PartialWaveOption so{ Options{ { "L"_ok, "P"_ov } }, specifications, "L"_ok };
                TEST_CHECK_EQUAL(so.value(), PartialWave::P);
            }

            // Creation of option with valid default value, with value != default value but containing default value
            {
                PartialWaveOption so{ Options{ { "L"_ok, "S|P|D"_ov } }, specifications, "L"_ok };
                TEST_CHECK_EQUAL(so.value(), PartialWave::S | PartialWave::P | PartialWave::D);
            }

            // Creation with invalid value
            {
                auto test = [specifications]() { PartialWaveOption so{ Options{ { "L"_ok, "F"_ov } }, specifications, "L"_ok }; };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }

            // Creation with invalid value but containing a valid value
            {
                auto test = [specifications]() { PartialWaveOption so{ Options{ { "L"_ok, "P|F"_ov } }, specifications, "L"_ok }; };
                TEST_CHECK_THROWS(InvalidOptionValueError, test());
            }
        }
} partial_wave_option_test;
