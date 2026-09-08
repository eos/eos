/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_TEST_TEST_PDFS_HH
#define EOS_GUARD_TEST_TEST_PDFS_HH 1

namespace eos::test
{
    /*
     * Truth and resolution PDFs whose convolutions have exact closed forms, for use in
     * unit tests only:
     *   Gaussian (x) Gaussian = N(mu, sqrt(sigma_1^2 + sigma_2^2))
     *   box(w1) (x) box(w2) = a trapezoid, a triangle when w1 == w2
     *   the 2D cases are separable products of the 1D ones
     *
     * The test runner's main() calls this before any test case runs, which matters because
     * Parameters::Defaults() snapshots the parameter defaults: a set constructed before the
     * registration would not know the PDFs' parameters. It is idempotent, so a test case may
     * call it again without harm.
     */
    void register_test_pdfs();
} // namespace eos::test

#endif
