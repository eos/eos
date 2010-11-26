/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
        }
} histogram_1_test;
