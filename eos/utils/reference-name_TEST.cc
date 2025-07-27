/* vim: set sw=4 sts=4 et tw=140 foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/utils/reference-name.hh>

#include <test/test.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class NameTest : public TestCase
{
    public:
        NameTest() :
            TestCase("name_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto p = rnp::Name("IKMvD"));
            TEST_CHECK_NO_THROW(auto p = rnp::Name("KMPW"));
            TEST_CHECK_NO_THROW(auto p = rnp::Name("CMS"));
            TEST_CHECK_NO_THROW(auto p = rnp::Name("BES2"));
            TEST_CHECK_NO_THROW(auto p = rnp::Name("LHCb"));

            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Name(""));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Name("["));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Name("[KMPW"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Name("IKMvD:"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Name("BES2:"));
        }
} name_test;

class YearTest : public TestCase
{
    public:
        YearTest() :
            TestCase("year_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto n = rnp::Year("0000"));
            TEST_CHECK_NO_THROW(auto n = rnp::Year("1905"));
            TEST_CHECK_NO_THROW(auto n = rnp::Year("1999"));
            TEST_CHECK_NO_THROW(auto n = rnp::Year("2010"));
            TEST_CHECK_NO_THROW(auto n = rnp::Year("2019"));

            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto n = rnp::Year(":2014]"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto n = rnp::Year(":2014"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto n = rnp::Year("2014A"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto n = rnp::Year(""));
        }
} year_test;

class IndexTest : public TestCase
{
    public:
        IndexTest() :
            TestCase("index_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto p = rnp::Index("A"));
            TEST_CHECK_NO_THROW(auto p = rnp::Index("Z"));
            TEST_CHECK_NO_THROW(auto p = rnp::Index("AA"));

            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Index(""));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Index("0A"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto p = rnp::Index("A]"));
        }
} index_test;

class ReferenceNameTest : public TestCase
{
    public:
        ReferenceNameTest() :
            TestCase("reference_name_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto qn = ReferenceName("A:2010A"));
            TEST_CHECK_NO_THROW(auto qn = ReferenceName("IKMvD:2014A"));
            TEST_CHECK_NO_THROW(auto qn = ReferenceName("BES2:2006A"));
            TEST_CHECK_NO_THROW(auto qn = ReferenceName("LHCb:2010A"));

            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName(""));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName("A"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName("A199"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName("A:199"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName("A:1999"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName("A:1999-"));
            TEST_CHECK_THROWS(ReferenceNameSyntaxError, auto qn = ReferenceName("[A:1999-B]"));
        }
} reference_name_test;
