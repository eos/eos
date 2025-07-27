/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#include <eos/utils/observable_stub.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class ObservableStubTest : public TestCase
{
    public:
        ObservableStubTest() :
            TestCase("observable_stub_test")
        {
        }

        virtual void
        run() const
        {
            // Setting and retrieval
            {
                Parameters    parameters = Parameters::Defaults();
                Parameter     p          = parameters["mass::c"];
                ObservablePtr o          = ObservablePtr(new ObservableStub(parameters, "mass::c"));

                TEST_CHECK_EQUAL(p(), o->evaluate());

                p = 0.0;
                TEST_CHECK_EQUAL(p(), o->evaluate());

                p = p.central();
                TEST_CHECK_EQUAL(p(), o->evaluate());
            }

            // Cloning w/ anonymous parameters
            {
                Parameters    parameters = Parameters::Defaults();
                Parameter     p          = parameters["mass::c"];
                ObservablePtr o1         = ObservablePtr(new ObservableStub(parameters, "mass::c"));
                ObservablePtr o2         = o1->clone();

                TEST_CHECK_EQUAL(o1->evaluate(), o2->evaluate());

                p = 0.0;
                TEST_CHECK(o1->evaluate() != o2->evaluate());

                p = p.central();
                TEST_CHECK_EQUAL(o1->evaluate(), o2->evaluate());
            }

            // Check name, kinematics, options after cloning
            {
                Parameters    parameters = Parameters::Defaults();
                ObservablePtr o1         = ObservablePtr(new ObservableStub(parameters, "mass::c"));
                ObservablePtr o2         = o1->clone();

                TEST_CHECK_EQUAL(o1->name().full(), o2->name().full());
                TEST_CHECK_EQUAL(o1->kinematics(), o2->kinematics());
                TEST_CHECK_EQUAL(o1->options(), o2->options());
            }
        }
} observable_stub_test;
