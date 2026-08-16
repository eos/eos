/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <eos/utils/yaml-schema.hh>

#include <test/test.hh>

#include <optional>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

using namespace test;
using namespace eos;

class YamlSchemaTest : public TestCase
{
    public:
        YamlSchemaTest() :
            TestCase("yaml_schema_test")
        {
        }

        static std::optional<std::string>
        first_violation(const YAML::Node & n, const std::vector<yaml_schema::Field> & fields)
        {
            std::optional<std::string> violation;
            yaml_schema::validate(n, fields, [&](const std::string & msg) { violation = msg; });

            return violation;
        }

        virtual void
        run() const
        {
            static const std::vector<yaml_schema::Field> fields{
                { "a",   yaml_schema::Kind::Scalar,  true, "a required scalar" },
                { "b",      yaml_schema::Kind::Map,  true,    "a required map" },
                { "c", yaml_schema::Kind::Sequence, false,  "an optional list" },
            };

            // {{{ well-formed input, all fields present
            {
                YAML::Node n = YAML::Load("a: 1\nb: {x: 1}\nc: [1, 2]");

                TEST_CHECK(! first_violation(n, fields).has_value());
            }
            /// }}}

            // {{{ well-formed input, optional field absent
            {
                YAML::Node n = YAML::Load("a: 1\nb: {x: 1}");

                TEST_CHECK(! first_violation(n, fields).has_value());
            }
            /// }}}

            // {{{ a stray key not listed in the schema is ignored
            {
                YAML::Node n = YAML::Load("a: 1\nb: {x: 1}\nstray: 1");

                TEST_CHECK(! first_violation(n, fields).has_value());
            }
            /// }}}

            // {{{ required scalar missing
            {
                YAML::Node n = YAML::Load("b: {x: 1}");

                auto violation = first_violation(n, fields);
                TEST_CHECK(violation.has_value());
                TEST_CHECK_EQUAL_STR(*violation, "required key 'a' not specified");
            }
            /// }}}

            // {{{ required map missing
            {
                YAML::Node n = YAML::Load("a: 1");

                auto violation = first_violation(n, fields);
                TEST_CHECK(violation.has_value());
                TEST_CHECK_EQUAL_STR(*violation, "required key 'b' not specified");
            }
            /// }}}

            // {{{ required scalar present but wrong kind
            {
                YAML::Node n = YAML::Load("a: [1, 2]\nb: {x: 1}");

                auto violation = first_violation(n, fields);
                TEST_CHECK(violation.has_value());
                TEST_CHECK_EQUAL_STR(*violation, "key 'a' not mapped to a scalar value");
            }
            /// }}}

            // {{{ required map present but wrong kind
            {
                YAML::Node n = YAML::Load("a: 1\nb: [1, 2]");

                auto violation = first_violation(n, fields);
                TEST_CHECK(violation.has_value());
                TEST_CHECK_EQUAL_STR(*violation, "key 'b' not mapped to a map");
            }
            /// }}}

            // {{{ optional sequence present but wrong kind
            {
                YAML::Node n = YAML::Load("a: 1\nb: {x: 1}\nc: 1");

                auto violation = first_violation(n, fields);
                TEST_CHECK(violation.has_value());
                TEST_CHECK_EQUAL_STR(*violation, "key 'c' not mapped to a sequence");
            }
            /// }}}

            // {{{ raise() may throw; validate() propagates the exception
            {
                YAML::Node n = YAML::Load("b: {x: 1}");

                TEST_CHECK_THROWS(std::runtime_error, yaml_schema::validate(n, fields, [](const std::string & msg) { throw std::runtime_error(msg); }));
            }
            /// }}}
        }
} yaml_schema_test;
