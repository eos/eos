/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

/*
 * Under Debian/Ubuntu, libgsl.so has no DT_NEEDED flag pointing to
 * libgslcblas.so. Force a need on the latter within libeosutils.so.
 *
 * Contents of eos::gsl_cblas_hack::gsl_cblas_hack() is an example from the
 * GSL CBLAS Manual.
 */

#include <gsl/gsl_cblas.h>

namespace eos
{
    namespace gsl_cblas_hack
    {
        void gsl_cblas_hack()
        {
            int lda = 3;

            float A[6] = { 0.11, 0.12, 0.13, 0.21, 0.22, 0.23 };

            int ldb = 2;

            float B[6] = { 1011, 1012, 1021, 1022, 1031, 1032 };

            int ldc = 2;

            float C[4] = { 0.00, 0.00, 0.00, 0.00 };

            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 3, 1.0, A, lda, B, ldb, 0.0, C, ldc);
        }
    }
}
