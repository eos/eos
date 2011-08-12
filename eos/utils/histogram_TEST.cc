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

#include <test/test.hh>
#include <eos/utils/histogram.hh>

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

class Histogram2Test :
    public TestCase
{
    public:
        Histogram2Test() :
            TestCase("histogram<2>_test")
        {
        }

        static std::string _stringify(const Histogram<2> & histogram)
        {
            std::stringstream stream;

            for (auto b = histogram.cbegin(), b_end = histogram.cend() ; b != b_end ; ++b)
            {
                stream << "[<" << b->lower[0] << ',' << b->upper[0] << ">,<" << b->lower[1] << ',' << b->upper[1] << ">," << b->value << ']';
            }

            return stream.str();
        }

        virtual void run() const
        {
            // Test insertion of values
            {
                Histogram<2> histogram = Histogram<2>::WithEqualBinning(std::array<double, 2>{{0.0, 0.0}}, std::array<double, 2>{{6.0, 6.0}},
                        std::array<unsigned, 2>{{2, 3}});

                TEST_CHECK_EQUAL(_stringify(histogram), "[<0,3>,<0,2>,0][<0,3>,<2,4>,0][<0,3>,<4,6>,0][<3,6>,<0,2>,0][<3,6>,<2,4>,0][<3,6>,<4,6>,0]");

                auto b = histogram.find(std::array<double, 2>{{2.0, 3.0}});
                TEST_CHECK(histogram.end() != b);
                b->value = 7.0;

                TEST_CHECK_EQUAL(_stringify(histogram), "[<0,3>,<0,2>,0][<0,3>,<2,4>,7][<0,3>,<4,6>,0][<3,6>,<0,2>,0][<3,6>,<2,4>,0][<3,6>,<4,6>,0]");

                b = histogram.find(std::array<double, 2>{{6.0, 6.0}});
                TEST_CHECK(histogram.end() == b);

                b = histogram.find(std::array<double, 2>{{5.9, 5.9}});
                TEST_CHECK(histogram.end() != b);
                b->value = -10.0;

                TEST_CHECK_EQUAL(_stringify(histogram), "[<0,3>,<0,2>,0][<0,3>,<2,4>,7][<0,3>,<4,6>,0][<3,6>,<0,2>,0][<3,6>,<2,4>,0][<3,6>,<4,6>,-10]");
            }

            // Test finding of bins
            {
                Histogram<2> histogram = Histogram<2>::WithEqualBinning(std::array<double, 2>{{0.0, 0.0}}, std::array<double, 2>{{6.0, 6.0}},
                        std::array<unsigned, 2>{{2, 3}});

                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.0, 0.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.0, 0.1}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.1, 0.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{2.9, 1.9}}));

                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.0, 2.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.0, 2.1}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.1, 2.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{2.9, 3.9}}));

                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.0, 4.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.0, 4.1}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{0.1, 4.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{2.9, 5.9}}));

                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.0, 0.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.0, 0.1}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.1, 0.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{5.9, 1.9}}));

                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.0, 2.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.0, 2.1}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.1, 2.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{5.9, 3.9}}));

                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.0, 4.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.0, 4.1}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{3.1, 4.0}}));
                TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{5.9, 5.9}}));
            }

            // Test not finding of bins
            {
                Histogram<2> histogram = Histogram<2>::WithEqualBinning(std::array<double, 2>{{0.0, 0.0}}, std::array<double, 2>{{6.0, 6.0}},
                        std::array<unsigned, 2>{{2, 3}});

                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, -0.1}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, +0.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, 1.9}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, 2.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, 3.9}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, 4.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, 5.9}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{-0.1, 6.0}}));

                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, -0.1}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, +0.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 1.9}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 2.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 3.9}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 4.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 5.9}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 6.0}}));

                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{0.0, -0.1}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{3.9, -0.1}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{3.0, -0.1}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{5.9, -0.1}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, -0.1}}));

                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{0.0, 6.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{3.9, 6.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{3.0, 6.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{5.9, 6.0}}));
                TEST_CHECK(histogram.end() == histogram.find(std::array<double, 2>{{6.0, 6.0}}));
            }

            // Test a common use case
            {
                static const std::array<double, 2> start{{0.0, 0.0}}, end{{15.25, 15.25}};
                static const std::array<unsigned, 2> count{{61, 61}};
                Histogram<2> histogram = Histogram<2>::WithEqualBinning(start, end, count);

                // check lattice points
                for (unsigned i = 0 ; i < 61 ; ++i)
                {
                    for (unsigned j = 0 ; j < 61 ; ++j)
                    {
                        TEST_CHECK(histogram.end() != histogram.find(std::array<double, 2>{{i * 0.25, j * 0.25}}));
                    }
                }
            }
        }
} histogram_2_test;
