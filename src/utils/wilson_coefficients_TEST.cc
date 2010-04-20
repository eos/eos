/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/wilson_coefficients.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace wf;

class WilsonCoefficientsTest :
    public TestCase
{
    public:
        WilsonCoefficientsTest() :
            TestCase("wilson_coefficients_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            Parameter c1(p["c1"]), c2(p["c2"]), c3(p["c3"]), c4(p["c4"]), c5(p["c5"]), c6(p["c6"]);

            for (auto i(0) ; i < 20 ; ++i)
            {
                double mu = 2.4 + i * (9.6 - 2.4) / 20.0;

                calculate_wilson_coefficients(mu, p);

                std::cout
                    << mu << '\t'
                    << c1 << '\t'
                    << c2 << '\t'
                    << c3 << '\t'
                    << c4 << '\t'
                    << c5 << '\t'
                    << c6 << '\t'
                    << std::endl;
            }
        }
} wilson_coefficients_test;
