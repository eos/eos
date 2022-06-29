/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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
#include <interpolation.hh>
#include <eos/utils/parameters.hh>

using namespace test;
using namespace eos;

class InterpolationTest :
    public TestCase
{
    public:
        InterpolationTest() :
            TestCase("interpolation_test")
        {
        }

        virtual void run() const
        {
            // standard initialisation
            {
                CSplineInterpolation interp = CSplineInterpolation({0.0, 1.0, 2.0}, {0.0, 1.0, 2.0});
                TEST_CHECK_NEARLY_EQUAL(interp(0.5), 0.5, 1e-14);
            }

            // uniform initialization syntax
            {
                CSplineInterpolation interp = {
                    {0.0, 1.0, 2.0},
                    {0.0, 1.0, 2.0}
                };
                TEST_CHECK_NEARLY_EQUAL(interp(0.5), 0.5, 1e-14);
            }

            // evaluate outside of data range: must throw
            {
                CSplineInterpolation interp = {
                    {0.0, 1.0, 2.0},
                    {0.0, 1.0, 2.0}
                };
                TEST_CHECK_THROWS(GSLError, interp(3.0));
            }

            // dimension mismatch: must throw
            {
                TEST_CHECK_THROWS(InternalError, CSplineInterpolation({
                    {0.0, 1.0, 2.0},
                    {0.0, 1.0, 2.0, 3.0}
                }));
            }
        }
} interpolation_test;

