/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_WILSON_POLYNOMIAL_HH
#define EOS_GUARD_SRC_UTILS_WILSON_POLYNOMIAL_HH 1

#include <src/utils/one-of.hh>
#include <src/utils/observable.hh>

#include <list>
#include <string>

namespace eos
{
    class Constant;

    class Sum;

    class Product;

    class Sine;

    class Cosine;

    typedef OneOf<Constant, Sum, Product, Sine, Cosine, Parameter> WilsonPolynomial;

    WilsonPolynomial make_polynomial(const ObservablePtr &, const Kinematics &, const std::list<std::string> &);

    class WilsonPolynomialPrinter
    {
        private:
            bool _pretty;

        public:
            WilsonPolynomialPrinter(const bool & pretty = true);

            std::string visit(const Constant & c);
            std::string visit(const Sum & s);
            std::string visit(const Product & p);
            std::string visit(const Sine & s);
            std::string visit(const Cosine & s);
            std::string visit(const Parameter & p);
    };

    class WilsonPolynomialEvaluator
    {
        public:
            double visit(const Constant & c);
            double visit(const Parameter & p);
            double visit(const Sum & s);
            double visit(const Product & p);
            double visit(const Sine & s);
            double visit(const Cosine & c);
    };
}


#endif
