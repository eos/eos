/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2022 Danny van Dyk
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
#include <eos/utils/stringify.hh>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>

namespace eos
{
    /*
     * Under Debian/Ubuntu, libgsl.so has no DT_NEEDED flag pointing to
     * libgslcblas.so. Force a need on the latter within libeosutils.so.
     *
     * Contents of eos::gsl_cblas_hack::gsl_cblas_hack() is an example from the
     * GSL CBLAS Manual.
     */

    namespace gsl_cblas_hack
    {
        void
        gsl_cblas_hack()
        {
            int lda = 3;

            float A[6] = { 0.11, 0.12, 0.13, 0.21, 0.22, 0.23 };

            int ldb = 2;

            float B[6] = { 1011, 1012, 1021, 1022, 1031, 1032 };

            int ldc = 2;

            float C[4] = { 0.00, 0.00, 0.00, 0.00 };

            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 3, 1.0, A, lda, B, ldb, 0.0, C, ldc);
        }
    } // namespace gsl_cblas_hack

    /*
     * Disable GSL default error handler and register our own handler, which
     * throws an eos::Exception of type eos::GSLError.
     */
    namespace gsl_error_handler
    {
        void
        gsl_error_handler(const char * reason, const char * /*file*/, int /*line*/, int gsl_errno)
        {
            throw GSLError(stringify(reason) + " (error code: " + stringify(gsl_errno) + ")");
        }

        void gsl_error_handler_hack() __attribute__((constructor));

        void
        gsl_error_handler_hack()
        {
            gsl_error_handler_t * previous_gsl_error_handler __attribute__((unused)) = gsl_set_error_handler(&gsl_error_handler);
        }
    } // namespace gsl_error_handler
} // namespace eos
