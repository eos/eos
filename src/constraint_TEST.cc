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

#include <test/test.hh>
#include <src/constraint.hh>

#include <iostream>

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
                Constraint c = Constraint::make("B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011", Options());
            }
        }
} constraint_test;
