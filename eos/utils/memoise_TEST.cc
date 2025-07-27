/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <eos/utils/memoise.hh>

#include <test/test.hh>

#include <complex>

using namespace test;
using namespace eos;

class MemoiseTest : public TestCase
{
    public:
        MemoiseTest() :
            TestCase("memoise_test")
        {
        }

        static double
        f1(const double & x, const double & y)
        {
            return x / y;
        }

        static std::complex<double>
        f2(const double & x, const double & y)
        {
            return std::complex<double>(x, y);
        }

        virtual void
        run() const
        {
            /* f1 */
            {
                TEST_CHECK_EQUAL(0.5, f1(1.0, 2.0));
                TEST_CHECK_EQUAL(2.0, f1(2.0, 1.0));

                // First round of memoisation
                TEST_CHECK_EQUAL(0, number_of_memoisations(f1, 0.0, 0.0));
                TEST_CHECK_EQUAL(0.5, memoise(f1, 1.0, 2.0));
                TEST_CHECK_EQUAL(1, number_of_memoisations(f1, 0.0, 0.0));
                TEST_CHECK_EQUAL(2.0, memoise(f1, 2.0, 1.0));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f1, 0.0, 0.0));

                // Second round of memoisation
                TEST_CHECK_EQUAL(0.5, memoise(f1, 1.0, 2.0));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f1, 0.0, 0.0));
                TEST_CHECK_EQUAL(2.0, memoise(f1, 2.0, 1.0));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f1, 0.0, 0.0));
            }

            /* f2 */
            {
                TEST_CHECK_EQUAL(1.0, real(f2(1.0, 2.0)));
                TEST_CHECK_EQUAL(2.0, imag(f2(1.0, 2.0)));

                // First round of memoisation
                TEST_CHECK_EQUAL(0, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(1.0, real(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(1, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(2.0, imag(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(1, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(2.0, real(memoise(f2, 2.0, 1.0)));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(1.0, imag(memoise(f2, 2.0, 1.0)));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));

                // Second round of memoisation
                TEST_CHECK_EQUAL(1.0, real(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(2.0, imag(memoise(f2, 1.0, 2.0)));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(2.0, real(memoise(f2, 2.0, 1.0)));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));
                TEST_CHECK_EQUAL(1.0, imag(memoise(f2, 2.0, 1.0)));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));
            }

            /* Test clearing all memoisations */
            {
                // There should be 2 memoisations per function
                TEST_CHECK_EQUAL(2, number_of_memoisations(f1, 0.0, 0.0));
                TEST_CHECK_EQUAL(2, number_of_memoisations(f2, 0.0, 0.0));

                MemoisationControl::instance()->clear();

                // There should be no memoisations left
                TEST_CHECK_EQUAL(0, number_of_memoisations(f1, 0.0, 0.0));
                TEST_CHECK_EQUAL(0, number_of_memoisations(f2, 0.0, 0.0));
            }
        }
} memoise_test;
