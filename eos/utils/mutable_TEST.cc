/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Danny van Dyk
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

#include <eos/utils/mutable.hh>
#include <eos/utils/parameters.hh>

#include <test/test.hh>

using namespace eos;
using namespace test;

class MutableAccessTest : public TestCase
{
    public:
        MutableAccessTest() :
            TestCase("mutable_access_test")
        {
        }

        virtual void
        run() const
        {
            // Test Parameter
            {
                Parameters parameters = Parameters::Defaults();
                Parameter  p          = parameters["mass::b(MSbar)"];

                Mutable & m = p;

                TEST_CHECK_EQUAL(static_cast<double>(m), static_cast<double>(p));
                TEST_CHECK_EQUAL(m(), p());
                TEST_CHECK_EQUAL(m.name(), p.name());
            }

            // Test Parameter Clone
            {
                Parameters parameters = Parameters::Defaults();
                Parameter  p          = parameters["mass::b(MSbar)"];

                MutablePtr m1(new Parameter(p));
                MutablePtr m2(m1->clone());

                TEST_CHECK_EQUAL(static_cast<double>(*m1), static_cast<double>(p));
                TEST_CHECK_EQUAL((*m1)(), p());
                TEST_CHECK_EQUAL(m1->name(), p.name());

                TEST_CHECK_EQUAL(static_cast<double>(*m2), static_cast<double>(p));
                TEST_CHECK_EQUAL((*m2)(), p());
                TEST_CHECK_EQUAL(m2->name(), p.name());
            }
        }
} mutable_access_test;
