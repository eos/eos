/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/observable.hh>
#include <src/utils/stringify.hh>
#include <src/utils/wilson-polynomial.hh>

#include <cmath>

namespace eos
{
    /* WilsonPolynomial elements */
    struct Constant
    {
        double value;

        Constant(const double & value) :
            value(value)
        {
        }
    };

    struct Sum
    {
        std::list<WilsonPolynomial> summands;

        Sum()
        {
        }

        Sum(const WilsonPolynomial & x, const WilsonPolynomial & y)
        {
            add(x);
            add(y);
        }

        void add(const WilsonPolynomial & summand)
        {
            summands.push_back(summand);
        }
    };

    struct Product
    {
        WilsonPolynomial x, y;

        Product() :
            x(Constant(0)),
            y(Constant(0))
        {
        }

        Product(const WilsonPolynomial & x, const WilsonPolynomial & y) :
            x(x),
            y(y)
        {
        };
    };

    struct Sine
    {
        WilsonPolynomial phi;

        Sine(const WilsonPolynomial & phi) :
            phi(phi)
        {
        }
    };

    struct Cosine
    {
        WilsonPolynomial phi;

        Cosine(const WilsonPolynomial & phi) :
            phi(phi)
        {
        }
    };

    /* Build a WilsonPolynomial from an observable */
    WilsonPolynomial make_polynomial(const ObservablePtr & o, const Kinematics & k, const std::list<std::string> & _coefficients)
    {
        Sum result;

        std::list<std::tuple<Parameter, Parameter, double, double, double>> coefficients; // <r, phi, q_i, c_i, s_i>
        for (auto i = _coefficients.cbegin(), i_end = _coefficients.cend() ; i != i_end ; ++i)
        {
            Parameter modulus = o->parameters()["Abs{" + *i + "}"];
            Parameter phase = o->parameters()["Arg{" + *i + "}"];

            coefficients.push_back(std::make_tuple(modulus, phase, 0.0, 0.0, 0.0));
        }

        /*
         * Wilson-Polynomials have the form
         *
         *   p = n
         *     + \sum_i q_i Abs{C_i}^2 + c_i Abs{C_i} Cos[Arg{C_i}] + s_i Abs{C_i} Sin[Arg{C_i}]
         *     + \sum_{i, j > i} Abs{C_i} Abs{C_j} (bl_{ij} + c_{ij} Cos[Arg{C_i} - Arg{C_j}] + s_{ij} Sin[Arg{C_i} - Arg{C_j}])
         */

        // Set all parameters to zero
        for (auto i = coefficients.begin(), i_end = coefficients.end() ; i != i_end ; ++i)
        {
            std::get<0>(*i) = 0.0; std::get<1>(*i) = 0.0;
        }

        // Determine the constant part 'n'
        double n = o->evaluate(k);
        result.add(Constant(n));

        // Determine the true quadratic terms 'q_i' and linear terms 'l_i'
        for (auto i = coefficients.begin(), i_end = coefficients.end() ; i != i_end ; ++i)
        {
            Parameter r_i = std::get<0>(*i), phi_i = std::get<1>(*i);

            r_i = r_i.central();
            // calculate observables
            phi_i = 0.0;
            double pzero = o->evaluate(k);
            phi_i = M_PI;
            double ppluspi = o->evaluate(k);
            phi_i = +M_PI / 2.0;
            double ppluspihalf = o->evaluate(k);
            phi_i = -M_PI / 2.0;
            double pminuspihalf = o->evaluate(k);

            double q_i = 0.5 * ((ppluspi + pzero) - 2.0 * n) / r_i() / r_i();
            result.add(Product(Constant(q_i), Product(r_i, r_i)));
            std::get<2>(*i) = q_i;

            double c_i = 0.5 * ((pzero - ppluspi)) / r_i();
            result.add(Product(Constant(c_i), Product(r_i, Cosine(phi_i))));
            std::get<3>(*i) = c_i;

            double s_i = 0.5 * (ppluspihalf - pminuspihalf) / r_i();
            result.add(Product(Constant(s_i), Product(r_i, Sine(phi_i))));
            std::get<4>(*i) = s_i;

            // reset parameters
            r_i = 0.0; phi_i = 0.0;
        }

        // Determine the bilinear terms 'bl_{ij}', 'c_{ij}' and 's_{ij}'
        for (auto i = coefficients.begin(), i_end = coefficients.end() ; i != i_end ; ++i)
        {
            Parameter r_i = std::get<0>(*i), phi_i = std::get<1>(*i);
            r_i = r_i.central(); phi_i = 0.0;
            double q_i = std::get<2>(*i), c_i = std::get<3>(*i);

            auto j = i;
            if (j != i_end)
                ++j;

            for ( ; j != i_end; ++j)
            {
                Parameter r_j = std::get<0>(*j), phi_j = std::get<1>(*j);
                r_j = r_j.central();
                double q_j = std::get<2>(*j), c_j = std::get<3>(*j), s_j = std::get<4>(*j);

                // calculate observables
                phi_j = +M_PI;
                double ppluspi = o->evaluate(k);
                phi_j = 0.0;
                double pzero = o->evaluate(k);
                phi_j = +M_PI / 2.0;
                double ppluspihalf = o->evaluate(k);
                phi_j = -M_PI / 2.0;
                double pminuspihalf = o->evaluate(k);

                // extract bilinear terms
                double c_ij = 0.5 * ((pzero - ppluspi) - 2.0 * r_j() * c_j) / r_i() / r_j();
                double b_ij = 0.5 * ((pzero + ppluspi) - 2.0 * (n + r_i() * r_i() * q_i + r_j() * r_j() * q_j + r_i() * c_i)) / r_i() / r_j();
                double s_ij = 0.5 * ((pminuspihalf - ppluspihalf) + 2.0 * r_j() * s_j) / r_i() / r_j();

                Product bilinear(r_i, r_j);
                result.add(Product(bilinear, Constant(b_ij)));
                result.add(Product(bilinear, Product(Constant(c_ij), Cosine(Sum(phi_i, Product(Constant(-1.0), phi_j))))));
                result.add(Product(bilinear, Product(Constant(s_ij), Sine(Sum(phi_i, Product(Constant(-1.0), phi_j))))));

                r_j = 0.0; phi_j = 0.0;
            }

            r_i = 0.0; phi_i = 0.0;
        }

        // Reset parameters to defaults
        for (auto i = coefficients.cbegin(), i_end = coefficients.cend() ; i != i_end ; ++i)
        {
            Parameter r_i = std::get<0>(*i), phi_i = std::get<1>(*i);

            r_i = r_i.central();
            phi_i = phi_i.central();
        }

        return result;
    }

    /* WilsonPolynomialPrinter */
    WilsonPolynomialPrinter::WilsonPolynomialPrinter(const bool & pretty) :
        _pretty(pretty)
    {
    }

    std::string
    WilsonPolynomialPrinter::visit(const Constant & c)
    {
        return stringify(c.value);
    }

    std::string
    WilsonPolynomialPrinter::visit(const Sum & s)
    {
        std::string result = "(";

        if (_pretty)
            result += "\n   ";

        auto i = s.summands.cbegin(), i_end = s.summands.cend();

        if (i != i_end)
        {
            result += (*i).accept_returning<std::string>(*this);
            ++i;
        }

        for ( ; i != i_end ; ++i)
        {
            result += std::string(_pretty ? "\n" : "") + " + " + (*i).accept_returning<std::string>(*this);
        }

        if (_pretty)
            result += "\n";

        result += ")";

        return result;
    }

    std::string
    WilsonPolynomialPrinter::visit(const Product & p)
    {
        return p.x.accept_returning<std::string>(*this) + " * " + p.y.accept_returning<std::string>(*this);
    }

    std::string
    WilsonPolynomialPrinter::visit(const Sine & s)
    {
        WilsonPolynomialPrinter printer(false);

        return "Sin[" + s.phi.accept_returning<std::string>(printer) + "]";
    }

    std::string
    WilsonPolynomialPrinter::visit(const Cosine & s)
    {
        WilsonPolynomialPrinter printer(false);

        return "Cos[" + s.phi.accept_returning<std::string>(printer) + "]";
    }

    std::string
    WilsonPolynomialPrinter::visit(const Parameter & p)
    {
        return p.name();
    }

    /* WilsonPolynomialEvaluator */
    double
    WilsonPolynomialEvaluator::visit(const Constant & c) { return c.value; }
    double
    WilsonPolynomialEvaluator::visit(const Parameter & p) { return p(); }

    double
    WilsonPolynomialEvaluator::visit(const Sum & s)
    {
        double result = 0.0;

        for (auto i = s.summands.cbegin(), i_end = s.summands.cend() ; i != i_end ; ++i)
        {
            result += (*i).accept_returning<double>(*this);
        }

        return result;
    }

    double
    WilsonPolynomialEvaluator::visit(const Product & p)
    {
        return p.x.accept_returning<double>(*this) * p.y.accept_returning<double>(*this);
    }

    double
    WilsonPolynomialEvaluator::visit(const Sine & s)
    {
        return std::sin(s.phi.accept_returning<double>(*this));
    }

    double
    WilsonPolynomialEvaluator::visit(const Cosine & c)
    {
        return std::cos(c.phi.accept_returning<double>(*this));
    }
}
