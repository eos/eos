/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/charm-loops.hh>
#include <src/rare-b-decays/long-distance.hh>

#include <cmath>
#include <vector>

#include <iostream>

using namespace test;
using namespace wf;

class LongDistanceTest :
    public TestCase
{
    public:
        LongDistanceTest() :
            TestCase("long_distance_test")
        {
        }

        virtual void run() const
        {
            static const double m_c = 1.2;
            static const double eps = 1e-6;
            static const std::vector<double> inputs{ 14.00, 15.00, 16.00, 19.21 };
            static const std::vector<complex<double>> results{
                complex<double>(6.716354, 0.381747),
                complex<double>(4.064683, 0.766175),
                complex<double>(3.230018, 2.113815),
                complex<double>(2.195244, 1.780107),
            };
            auto r = results.cbegin();
            for (auto i = inputs.cbegin() ; inputs.cend() != i ; ++i, ++r)
            {
                complex<double> g = LongDistance::g_had_ccbar(*i, m_c);

                TEST_CHECK_NEARLY_EQUAL(real(*r), real(g), eps);
                TEST_CHECK_NEARLY_EQUAL(imag(*r), imag(g), eps);
            }
        }
} long_distance_test;
