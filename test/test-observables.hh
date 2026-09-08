/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2026 Danny van Dyk
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

#ifndef EOS_GUARD_TEST_TEST_OBSERVABLES_HH
#define EOS_GUARD_TEST_TEST_OBSERVABLES_HH 1

namespace eos::test
{
    /*
     * Observables for use in unit tests only:
     *   test::obs1 = mass::c * multiplier * (q2_max - q2_min)
     *
     * The test runner's main() calls this before any test case runs; it is idempotent, so a
     * test case may call it again without harm.
     */
    void register_test_observables();
} // namespace eos::test

#endif
