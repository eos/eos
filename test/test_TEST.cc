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

#include <eos/maths/complex.hh>

#include <test/test.hh>

#include <cmath>

using namespace test;
using namespace eos;

class NoThrowTest : public TestCase
{
    public:
        NoThrowTest() :
            TestCase("no_throw_test")
        {
        }

        virtual void
        run() const
        {
            try
            {
                TEST_CHECK_NO_THROW(throw std::string("failed"));
            }
            catch (TestCaseFailedException & e)
            {
                // as should be
            }
            catch (std::string & e)
            {
                TEST_CHECK_FAILED("failed");
            }
        }
} no_throw_test;

class EqualTest : public TestCase
{
    public:
        EqualTest() :
            TestCase("equal_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(TEST_CHECK_EQUAL(0, 0));
            TEST_CHECK_NO_THROW(TEST_CHECK_EQUAL(std::string("foo"), std::string("foo")));
            TEST_CHECK_NO_THROW(TEST_CHECK_EQUAL(0.0, 0.0));

            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_EQUAL(0, 1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_EQUAL(std::string("foo"), std::string("bar")));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_EQUAL(17.0, 23.0));
        }
} equal_test;

class RelativeErrorTest : public TestCase
{
    public:
        RelativeErrorTest() :
            TestCase("relative_error_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(TEST_CHECK_RELATIVE_ERROR(1.0, 1.09, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(+1.0, +2.0, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(+1.0, -2.0, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(-1.0, +2.0, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(-1.0, -2.0, 0.1));

            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(-0.1, 0.0, 0.2));
        }
} relative_error_test;

class RelativeErrorComplexTest : public TestCase
{
    public:
        RelativeErrorComplexTest() :
            TestCase("relative_error_complex_test")
        {
        }

        virtual void
        run() const
        {
            complex<double> x(1.0, 2.0);
            TEST_CHECK_NO_THROW(TEST_CHECK_RELATIVE_ERROR_C(x, complex<double>(1.03, 2.1), 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR_C(x, complex<double>(1.15, 2.1), 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR_C(x, complex<double>(1.0, 2.5), 0.1));
        }
} relative_error_complex_test;
