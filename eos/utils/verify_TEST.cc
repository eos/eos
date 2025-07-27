/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
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

#include <eos/utils/verify.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class VerifiedRangeTest : public TestCase
{
    public:
        VerifiedRangeTest() :
            TestCase("verified_range_test")
        {
        }

        virtual void
        run() const
        {
            // proper construction
            TEST_CHECK_NO_THROW(VerifiedRange<double> test_range(0.0, 1.0, 0.5));
            VerifiedRange<double> test_range(0.0, 1.0, 0.5);

            // false assignments
            TEST_CHECK_THROWS(VerifiedRangeUnderflow, test_range = -800);
            TEST_CHECK_THROWS(VerifiedRangeOverflow, test_range = +100);

            // retrieve value
            TEST_CHECK_EQUAL(test_range, 0.5);
            TEST_CHECK_NO_THROW(test_range = 0.9);
            TEST_CHECK_EQUAL(test_range, 0.9);

            // constructor exceptions
            TEST_CHECK_THROWS(VerifiedRangeUnderflow, VerifiedRange<double> test_range(0.0, 1.0, -0.1));
            TEST_CHECK_THROWS(VerifiedRangeOverflow, VerifiedRange<double> test_range(0.0, 1.0, 1.1));
        }
} verified_range_test;
