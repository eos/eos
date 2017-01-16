/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2014, 2015, 2017 Danny van Dyk
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
#include <eos/constraint.hh>
#include <eos/statistics/log-likelihood.hh>

#include <iostream>
#include <vector>

using namespace test;
using namespace eos;

class ConstraintTest :
    public TestCase
{
    public:
        ConstraintTest() :
            TestCase("constraint_test")
        {
        }

        virtual void run() const
        {
            /* Test making constraints */
            {
                std::cout << "# Constraints :" << std::endl;

                Options o;
                auto constraints = Constraints();
                unsigned n = 0;

                for (auto cf = constraints.begin(); cf != constraints.end(); ++cf, ++n)
                {
                    std::cout << "#  " << cf->first.full() << ": ";

                    Constraint c = Constraint::make(cf->first, o);
                    TEST_CHECK_EQUAL(c.name(), cf->first);
                    TEST_CHECK(std::distance(c.begin_observables(), c.end_observables()) > 0);
                    TEST_CHECK(std::distance(c.begin_blocks(), c.end_blocks()) > 0);

                    for (auto o = c.begin_observables(), o_end = c.end_observables(); o != o_end ; ++o)
                    {
                        std::cout << (**o).name() << '['
                                << (**o).kinematics().as_string() << ']'
                                << " with options: " << (**o).options().as_string();
                    }
                    for (auto b = c.begin_blocks(), b_end = c.end_blocks(); b != b_end ; ++b)
                    {
                        std::cout << ", " << (**b).as_string();
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
                std::cout << "# Found " << n << " constraints" << std::endl;
            }
        }
} constraint_test;
