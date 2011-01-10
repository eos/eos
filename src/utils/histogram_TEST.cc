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

#include <test/test.hh>
#include <src/utils/histogram.hh>

#include <iostream>

using namespace test;
using namespace eos;

class Histogram1Test :
    public TestCase
{
    public:
        Histogram1Test() :
            TestCase("histogram<1>_test")
        {
        }

        static std::string _stringify(const Histogram<1> & histogram)
        {
            std::stringstream stream;

            for (auto b = histogram.begin(), b_end = histogram.end() ; b != b_end ; ++b)
            {
                stream << '[' << b->lower << ',' << b->upper << ',' << b->value << ']';
            }

            return stream.str();
        }

        virtual void run() const
        {
            // Test insertion of values
            {
                Histogram<1> histogram = Histogram<1>::WithEqualBinning(0.0, 7.0, 7);

                TEST_CHECK_EQUAL(_stringify(histogram), "[0,1,0][1,2,0][2,3,0][3,4,0][4,5,0][5,6,0][6,7,0]");

                histogram.insert(3.1415);
                TEST_CHECK_EQUAL(_stringify(histogram), "[0,1,0][1,2,0][2,3,0][3,4,1][4,5,0][5,6,0][6,7,0]");
                TEST_CHECK_EQUAL(histogram.entries(), 1);

                histogram.insert(5);
                TEST_CHECK_EQUAL(_stringify(histogram), "[0,1,0][1,2,0][2,3,0][3,4,1][4,5,0][5,6,1][6,7,0]");
                TEST_CHECK_EQUAL(histogram.entries(), 2);

                histogram.insert(3);
                TEST_CHECK_EQUAL(_stringify(histogram), "[0,1,0][1,2,0][2,3,0][3,4,2][4,5,0][5,6,1][6,7,0]");
                TEST_CHECK_EQUAL(histogram.entries(), 3);
            }

            // ToDo: Test ECDF
            {
            }
        }
} histogram_1_test;
