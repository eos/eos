/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#include <eos/utils/exception.hh>

#include <test/test.hh>

#include <gsl/gsl_sf.h>

using namespace test;
using namespace eos;

class GSLErrorHandlerTest : public TestCase
{
    public:
        GSLErrorHandlerTest() :
            TestCase("gsl_error_handler_test")
        {
        }

        virtual void
        run() const
        {
            gsl_sf_result result;
            int           status;

            TEST_CHECK_THROWS(GSLError, { status = gsl_sf_expint_3_e(-1.0, &result); });
        }
} gsl_error_handler_test;
