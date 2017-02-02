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

#include <test/test.hh>
#include <eos/utils/gsl-interface.hh>

using namespace test;
using namespace eos;

class GSLMatrixPtrTestCase :
    public TestCase
{
    public:
        GSLMatrixPtrTestCase() :
            TestCase("gsl_matrix_ptr_test_case")
        {
        }

        virtual void run() const
        {
            // default construction
            {
                GSLMatrixPtr();
            }

            // allocation with alloc
            {
                GSLMatrixPtr p(gsl_matrix_alloc(10, 10));
            }

            // allocation with calloc
            {
                GSLMatrixPtr p(gsl_matrix_calloc(10, 10));
            }

            // allocation with make_gsl_matrix
            {
                GSLMatrixPtr p(make_gsl_matrix(10, 10));
            }
        }
} gsl_matrix_ptr_test_case;

class GSLVectorPtrTestCase :
    public TestCase
{
    public:
        GSLVectorPtrTestCase() :
            TestCase("gsl_vector_ptr_test_case")
        {
        }

        virtual void run() const
        {
            // default construction
            {
                GSLVectorPtr();
            }

            // allocation with alloc
            {
                GSLVectorPtr p(gsl_vector_alloc(10));
            }

            // allocation with calloc
            {
                GSLVectorPtr p(gsl_vector_calloc(10));
            }

            // allocation with make_gsl_vector
            {
                GSLVectorPtr p(make_gsl_vector(10));
            }
        }
} gsl_vector_ptr_test_case;
