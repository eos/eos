/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Christian Wacker
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
#include <eos/utils/polylog.hh>

#include <fstream>

#include <iomanip>

using namespace test;
using namespace eos;

class PolylogarithmTest :
    public TestCase
{
    public:
        PolylogarithmTest() :
            TestCase("polylogarithm_test")
        {
        }

        static complex<double> read_complex(std::istream& input)
        {
            complex<double> result;

            input.read(reinterpret_cast<char*>(&result), sizeof(result));

            return result;
        }

        virtual void run() const
        {
            std::ifstream z_file(EOS_SRCDIR "/eos/utils/polylog_TEST_z.bin");
            std::ifstream dilog_file(EOS_SRCDIR "/eos/utils/polylog_TEST_dilog.bin");
            std::ifstream trilog_file(EOS_SRCDIR "/eos/utils/polylog_TEST_trilog.bin");

            static const double eps = 1e-12;

            // check that the C implementation uses the principal branch of the (complex) logarithm
            TEST_CHECK_RELATIVE_ERROR(+M_PI, imag(std::log(complex<double>(-1.0, 0.0))), eps);

            while (! z_file.eof())
            {
                complex<double> z = read_complex(z_file);
                complex<double> dilog_reference = read_complex(dilog_file);
                complex<double> dilog_value = dilog(z);
                complex<double> trilog_reference = read_complex(trilog_file);
                complex<double> trilog_value = trilog(z);

                TEST_CHECK_NEARLY_EQUAL(real(dilog_reference),  real(dilog_value),  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(dilog_reference),  imag(dilog_value),  eps);
                TEST_CHECK_NEARLY_EQUAL(real(trilog_reference), real(trilog_value), eps);
                TEST_CHECK_NEARLY_EQUAL(imag(trilog_reference), imag(trilog_value), eps);
            }

            // check that the C implementation handles z = -2.0 - i 0.0
            // correctly
            static const complex<double> c1(1.0, 0.0), c05(0.5, 0.0), c2(2.0, 0.0);
            static const complex<double> z    = (c2 - c1) / (c05 - c1); // which turns out to be (-2.0,-0.0)
            static const complex<double> zbar = (c05 - c1) / (c2 - c1); // which turns out to be (-0.5,-0.0)
            TEST_CHECK_RELATIVE_ERROR(real(dilog(-c2)),   +real(dilog(z)),     eps); // has no imaginary part
            TEST_CHECK_RELATIVE_ERROR(real(trilog(-c2)),  +real(trilog(z)),    eps); // has no imaginary part
            TEST_CHECK_RELATIVE_ERROR(real(dilog(-c05)),  +real(dilog(zbar)),  eps); // has no imaginary part
            TEST_CHECK_RELATIVE_ERROR(real(trilog(-c05)), +real(trilog(zbar)), eps); // has no imaginary part
        }
} polylogarithm_test;
