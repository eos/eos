/*
 * Copyright (c) 2023 Stefan Meiser
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
#include <eos/form-factors/vec-lcdas.hh>

#include <eos/models/model.hh>
#include <eos/maths/power-of.hh>

using namespace test;
using namespace eos;

class VectorLCDAsTest :
    public TestCase
{
    public:
        VectorLCDAsTest() :
            TestCase("vector_lcdas_test")
        {
        }

        virtual void run() const
        {
            // creation
            {
                auto parameter = Parameters::Defaults();
                auto options   = Options();

                TEST_CHECK_THROWS(InternalError, VectorLCDAs::make("FooBar", parameter, options));
            }
        }
} vector_lcdas_test;
