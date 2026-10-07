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

#include <eos/statistics/unbinned-observations.hh>
#include <eos/utils/exception.hh>

#include <test/test.hh>

#include <vector>

using namespace test;
using namespace eos;

class UnbinnedObservationsTest : public TestCase
{
    public:
        UnbinnedObservationsTest() :
            TestCase("unbinned_observations_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_THROWS(InternalError, UnbinnedObservations(std::vector<double>{}, 1));
            TEST_CHECK_THROWS(InternalError, UnbinnedObservations(std::vector<double>{ 1.0, 2.0, 3.0 }, 2));
            TEST_CHECK_THROWS(InternalError, UnbinnedObservations(std::vector<double>{ 1.0 }, 0));

            UnbinnedObservations observations(std::vector<double>{ 1.0, 2.0, 3.0, 4.0 }, 2);
            TEST_CHECK_EQUAL(observations.rank(), 2u);
            TEST_CHECK_EQUAL(observations.size(), 2u);
            TEST_CHECK_EQUAL(observations.values()[3], 4.0);

            // copies share the events, and see their replacement
            UnbinnedObservations                   copy       = observations;
            const UnbinnedObservations::Generation generation = observations.generation();
            TEST_CHECK(copy.generation() == generation);

            observations.set(std::vector<double>{ 5.0, 6.0 });
            TEST_CHECK_EQUAL(copy.size(), 1u);
            TEST_CHECK_EQUAL(copy.values()[1], 6.0);
            TEST_CHECK(! (copy.generation() == generation));
            TEST_CHECK(copy.generation() == observations.generation());

            // a rejected replacement leaves the events and the generation unchanged
            const UnbinnedObservations::Generation current = observations.generation();
            TEST_CHECK_THROWS(InternalError, copy.set(std::vector<double>{ 7.0 }));
            TEST_CHECK_EQUAL(observations.values()[0], 5.0);
            TEST_CHECK(observations.generation() == current);

            // independent sets of observations never share a generation
            UnbinnedObservations other(std::vector<double>{ 5.0, 6.0 }, 2);
            TEST_CHECK(! (other.generation() == UnbinnedObservations(std::vector<double>{ 5.0, 6.0 }, 2).generation()));
            TEST_CHECK(! (other.generation() == UnbinnedObservations::Generation()));
        }
} unbinned_observations_test;
