/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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
#include <eos/utils/one-of.hh>

#include <array>

using namespace test;
using namespace eos;

class OneOfTest :
    public TestCase
{
    public:
        OneOfTest() :
            TestCase("one_of_test")
        {
        }

        virtual void run() const
        {
            typedef OneOf<int, std::string> Type;

            Type x = 0;
            Type y = std::string("foo");
            Type z = 10;
        }
} one_of_test;

struct Foo {};
struct Bar {};
struct Baz {};

class TestVisitorReturningVoid
{
    public:
        std::string result;

        void visit(const Foo &)
        {
            result += "Foo";
        }

        void visit(Bar &)
        {
            result += "Bar";
        }

        void visit(Baz)
        {
            result += "Baz";
        }
};

class OneOfVisitorReturningVoidTest :
    public TestCase
{
    public:
        OneOfVisitorReturningVoidTest() :
            TestCase("one_of_visitor_returning_void_test")
        {
        }

        virtual void run() const
        {
            typedef OneOf<Foo, Bar, Baz> Type;
            std::array<Type, 5> items
            {{
                Type{Foo()},
                Type{Bar()},
                Type{Bar()},
                Type{Foo()},
                Type{Baz()}
            }};

            TestVisitorReturningVoid visitor;
            for (auto i = items.cbegin(), i_end = items.cend() ; i != i_end ; ++i)
            {
                (*i).accept(visitor);
            }

            TEST_CHECK_EQUAL("FooBarBarFooBaz", visitor.result);
        }
} one_of_visitor_returning_void_test;

class TestVisitorReturningString
{
    public:
        std::string visit(const Foo &)
        {
            return "Foo";
        }

        std::string visit(Bar &)
        {
            return "Bar";
        }

        std::string visit(Baz)
        {
            return "Baz";
        }
};

class OneOfVisitorReturningStringTest :
    public TestCase
{
    public:
        OneOfVisitorReturningStringTest() :
            TestCase("one_of_visitor_returning_string_test")
        {
        }

        virtual void run() const
        {
            typedef OneOf<Foo, Bar, Baz> Type;
            std::array<Type, 5> items
            {{
                Type{Foo()},
                Type{Bar()},
                Type{Bar()},
                Type{Foo()},
                Type{Baz()}
            }};

            TestVisitorReturningString visitor;
            std::string result;
            for (auto i = items.cbegin(), i_end = items.cend() ; i != i_end ; ++i)
            {
                result += (*i).accept_returning<std::string>(visitor);
            }

            TEST_CHECK_EQUAL("FooBarBarFooBaz", result);
        }
} one_of_visitor_returning_string_test;
