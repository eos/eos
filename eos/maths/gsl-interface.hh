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

#ifndef EOS_GUARD_EOS_MATHS_GSL_INTERFACE_HH
#define EOS_GUARD_EOS_MATHS_GSL_INTERFACE_HH 1

#include <memory>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

namespace eos
{
    namespace gsl
    {
        /**
         * Deletes a gsl_matrix upon release from std::unique_ptr<>.
         */
        struct MatrixDeleter
        {
            MatrixDeleter() { }

            void operator() (gsl_matrix * m)
            {
                if (m)
                {
                    gsl_matrix_free(m);
                }
            }
        };
    }

    using GSLMatrixPtr = std::unique_ptr<gsl_matrix, gsl::MatrixDeleter>;

    GSLMatrixPtr make_gsl_matrix(size_t n1, size_t n2)
    {
        gsl_matrix * result = gsl_matrix_calloc(n1, n2);

        return GSLMatrixPtr(result);
    }

    namespace gsl
    {
        /**
         * Deletes a gsl_vector upon release from std::unique_ptr<>.
         */
        struct VectorDeleter
        {
            VectorDeleter() { }

            void operator() (gsl_vector * v)
            {
                if (v)
                {
                    gsl_vector_free(v);
                }
            }
        };
    }

    using GSLVectorPtr = std::unique_ptr<gsl_vector, gsl::VectorDeleter>;

    GSLVectorPtr make_gsl_vector(size_t n)
    {
        gsl_vector * result = gsl_vector_calloc(n);

        return GSLVectorPtr(result);
    }
}

#endif
