/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#ifndef WILSON_FITTER_GUARD_TEST_TEST_HH
#define WILSON_FITTER_GUARD_TEST_TEST_HH 1

#include <limits>
#include <string>
#include <src/utils/stringify.hh>

namespace test
{
    class TestCase
    {
        private:
            std::string _name;

        public:
            TestCase(const std::string & name);

            virtual ~TestCase();

            std::string name() const;

            virtual void run() const = 0;
    };

    class TestCaseFailedException
    {
        private:
            int _line;

            std::string _file;

            std::string _reason;

        public:
            TestCaseFailedException(int line, const std::string & file, const std::string & reason);

            const std::string & reason() const;

            std::string where() const;
    };

#define TEST_CHECK(a) \
    do \
    { \
        if (! (a)) \
            throw TestCaseFailedException(__LINE__, __FILE__, "'" #a "' is false"); \
    } \
    while (false)

#define TEST_CHECK_EQUAL(a, b) \
    do \
    { \
        if (! ((a) == (b))) \
            throw TestCaseFailedException(__LINE__, __FILE__, "'" #a "' is not equal to '" #b "'"); \
    } \
    while (false)

#define TEST_CHECK_FAILED(s) \
    do \
    { \
        throw TestCaseFailedException(__LINE__, __FILE__, s); \
    } \
    while (false)

#define TEST_CHECK_NEARLY_EQUAL(a, b, eps) \
    do \
    { \
        if (std::abs((a - b)) <= eps) \
            break; \
        else \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #a "' = " + stringify(a) + " is not nearly-equal to '" #b "' = " + stringify(b) + " within '" + stringify(eps) + "'" \
                    + ", difference is '" + stringify(a - b) + "'"); \
    } \
    while (false)

#define TEST_CHECK_RELATIVE_ERROR(a, b, eps) \
    do \
    { \
        if (std::sqrt(std::abs(a)) < std::numeric_limits<decltype(a)>::epsilon()) \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #a "' has been evaluated to the zero within computational accuracy, result = " + stringify(a)); \
         \
        if (std::sqrt(std::abs(b)) < std::numeric_limits<decltype(b)>::epsilon()) \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #b "' has been evaluated to the zero within computational accuracy, result = " + stringify(b)); \
         \
        if (((std::abs((a - b) / a)) <= eps) && ((std::abs((a - b) / b)) <= eps)) \
            break; \
        else \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "One relative error of '" #a "' = '" + stringify(a) + "' and '" #b "' = '" + stringify(b) + "' is greater than " + stringify(eps) + ". The results are " + \
                    stringify(std::abs((a - b) / a)) + " and " + stringify(std::abs((a - b) / b))); \
    } \
    while (false)

#define TEST_CHECK_NO_THROW(expression) \
    do \
    { \
        try \
        { \
            expression; \
        } \
        catch (...) \
        { \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "Caught unexpected exception in '" #expression "'"); \
        } \
    } \
    while (false)

#define TEST_CHECK_THROWS(exception, expression) \
    do \
    { \
        try \
        { \
            expression; \
        } \
        catch (exception & e) \
        { \
            break; \
        } \
        catch (...) \
        { \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "Caught unexpected exception when expecting " #exception " in '" #expression "'"); \
        } \
        \
        throw TestCaseFailedException(__LINE__, __FILE__, \
                "Caught no exception in " #expression " when expecting '" #exception "'"); \
    } \
    while (false)
}

#endif
