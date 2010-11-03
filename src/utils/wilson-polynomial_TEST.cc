/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/complex.hh>
#include <src/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <vector>

#include <iostream>

using namespace test;
using namespace eos;

struct WilsonPolynomialTestObservable :
    public Observable
{
    std::string n;
    Parameters p;
    Parameter c1;
    Parameter c2;
    Parameter abs_c7;
    Parameter arg_c7;
    Parameter abs_c9;
    Parameter arg_c9;
    Parameter abs_c10;
    Parameter arg_c10;

    WilsonPolynomialTestObservable(const Parameters & p, const ObservableOptions &) :
        n("WilsonPolynomialTestObservable"),
        p(p),
        c1(p["c1"]),
        c2(p["c2"]),
        abs_c7(p["Abs{c7}"]),
        arg_c7(p["Arg{c7}"]),
        abs_c9(p["Abs{c9}"]),
        arg_c9(p["Arg{c9}"]),
        abs_c10(p["Abs{c10}"]),
        arg_c10(p["Arg{c10}"])
    {
    }

    virtual const std::string & name() const { return n; }
    virtual Parameters parameters() { return p; }
    virtual ObservablePtr clone() const { return ObservablePtr(new WilsonPolynomialTestObservable(p.clone(), ObservableOptions())); }

    virtual double evaluate(const Kinematics &) const
    {
        complex<double> c7 = abs_c7() * complex<double>(std::cos(arg_c7()), std::sin(arg_c7));
        complex<double> c9 = abs_c9() * complex<double>(std::cos(arg_c9()), std::sin(arg_c9));
        complex<double> c10 = abs_c10() * complex<double>(std::cos(arg_c10()), std::sin(arg_c10));

        return real(
                    0.01234
                    + c7 * complex<double>(0.321, 1.000)
                    + c9 * complex<double>(0.731, 1.000)
                    + conj(c7) * c7 * 0.6
                    + conj(c7) * c9 * complex<double>(1.300, 0.123)
                    + conj(c9) * c9 * 2.1
                    + conj(c10) * c10 * 1.23
                );
    }
};

class WilsonPolynomialTest :
    public TestCase
{
    public:
        WilsonPolynomialTest() :
            TestCase("wilson_polynomial_test")
        {
        }

        void run_one(const ObservablePtr & o, const WilsonPolynomial & p, const std::array<double, 6> & values) const
        {
            Parameters parameters = o->parameters();
            Parameter abs_c7(parameters["Abs{c7}"]);
            Parameter arg_c7(parameters["Arg{c7}"]);
            Parameter abs_c9(parameters["Abs{c9}"]);
            Parameter arg_c9(parameters["Arg{c9}"]);
            Parameter abs_c10(parameters["Abs{c10}"]);
            Parameter arg_c10(parameters["Arg{c10}"]);

            abs_c7 = values[0];
            arg_c7 = values[1];
            abs_c9 = values[2];
            arg_c9 = values[3];
            abs_c10 = values[4];
            arg_c10 = values[5];

            static const double eps = 1e-10;
            WilsonPolynomialEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(Kinematics()), p.accept_returning<double>(evaluator), eps);
        }

        virtual void run() const
        {
            Parameters parameters = Parameters::Defaults();
            Kinematics kinematics;

            ObservablePtr o = ObservablePtr(new WilsonPolynomialTestObservable(parameters, ObservableOptions()));
            WilsonPolynomial p = make_polynomial(o, kinematics, std::list<std::string>{ "c7", "c9", "c10" });

            WilsonPolynomialPrinter printer;
            std::cout << p.accept_returning<std::string>(printer) << std::endl;

            static const std::vector<std::array<double, 6>> inputs
            {
                std::array<double, 6>{{0.0,       0.0,       0.0,       0.0,       0.0,       0.0      }},
                std::array<double, 6>{{1.0,       0.0,       1.0,       0.0,       1.0,       0.0      }},
                std::array<double, 6>{{0.7808414, 0.8487257, 0.7735165, 0.5383695, 0.6649164, 0.7235497}},
                std::array<double, 6>{{0.5860642, 0.9830907, 0.7644369, 0.8330194, 0.4935018, 0.4492084}},
                std::array<double, 6>{{0.2177456, 0.5062894, 0.6463376, 0.3624364, 0.6770480, 0.0718421}},
                std::array<double, 6>{{0.0088306, 0.9441413, 0.8721501, 0.2984633, 0.2961408, 0.9145809}},
                std::array<double, 6>{{0.7967655, 0.2427081, 0.8403112, 0.3351082, 0.6477823, 0.5569495}},
                std::array<double, 6>{{0.7607454, 0.5025871, 0.5877762, 0.5516025, 0.2930899, 0.4882813}},
            };

            for (auto i = inputs.cbegin(), i_end = inputs.cend() ; i != i_end ; ++i)
            {
                run_one(o, p, *i);
            }
        }
} wilson_polynomial_test;
