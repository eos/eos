/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Frederik Beaujean
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
#include <eos/statistics/simple-parameters.hh>
#include <test/test.hh>

using namespace eos;
using namespace test;

class SimpleParametersTest :
    public TestCase
{
    public:
        SimpleParametersTest() :
            TestCase("simple_parameters_test")
        {
        }

        virtual void run() const
        {
            // create, access, and modify
            {
                SimpleParameters p;
                auto mH = p.declare("mH", 120, 130);
                TEST_CHECK_EQUAL(mH.name(), "mH");

                mH.set(125);
                TEST_CHECK_EQUAL(double(mH), 125);

                p["mH"] = 129;
                TEST_CHECK_EQUAL(mH, 129);

                p[0] = 128;
                TEST_CHECK_EQUAL(mH, 128);

                TEST_CHECK_EQUAL(p.values().size(), 1);
                TEST_CHECK_EQUAL(p.values().front(), 128);

                TEST_CHECK( (p != p) == false);

                TEST_CHECK_EQUAL(p.begin()->min, 120);
                TEST_CHECK_EQUAL(p.begin()->max, 130);
                TEST_CHECK_EQUAL(p.begin()->nuisance, false);
            }

            // cloning
            {
                SimpleParameters p1;
                p1.declare("mH", 120, 130);
                p1.declare("mt", 170, 180);

                p1[0] = 125;
                p1[1] = 174;

                SimpleParameters p2 = p1.clone();

                TEST_CHECK(p1 != p2);
                TEST_CHECK_EQUAL(p1[0], p2[0]);

                // now modify p1, does p2 change?
                p1[0] = 126;
                TEST_CHECK_EQUAL(p2[0], 125);

                p2[1] = 173;
                TEST_CHECK_EQUAL(p1[1], 174);
            }
        }
} simple_parameters_test;
