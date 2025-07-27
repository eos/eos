/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2015, 2016 Danny van Dyk
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

#include <eos/observable.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wilson-polynomial.hh>

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

            Sum() {}

            Sum(const WilsonPolynomial & x, const WilsonPolynomial & y)
            {
                add(x);
                add(y);
            }

            void
            add(const WilsonPolynomial & summand)
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
            }
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
    WilsonPolynomial
    make_polynomial(const ObservablePtr & o, const std::list<std::string> & _coefficients)
    {
        Sum result;

        std::list<std::tuple<Parameter, double, double>> coefficients; // <p_i, q_i, l_i>
        for (const auto & _coefficient : _coefficients)
        {
            coefficients.push_back(std::make_tuple(o->parameters()[_coefficient], 0.0, 0.0));
        }

        /*
         * Wilson-Polynomials have the form
         *
         *   p = n
         *     + \sum_i q_i P_i^2 + l_i P_i
         *     + \sum_{i, j > i} c_ij P_i P_j
         */

        // Set all parameters to zero
        for (auto & coefficient : coefficients)
        {
            std::get<0>(coefficient) = 0.0;
        }

        // Determine the constant part 'n'
        double n = o->evaluate();
        result.add(Constant(n));

        // Determine the true quadratic terms 'q_i' and linear terms 'l_i'
        for (auto & coefficient : coefficients)
        {
            Parameter p_i = std::get<0>(coefficient);


            // calculate observables
            p_i               = +1.0;
            double o_plus_one = o->evaluate();

            p_i                = -1.0;
            double o_minus_one = o->evaluate();

            double q_i = 0.5 * ((o_plus_one + o_minus_one) - 2.0 * n);
            result.add(Product(Constant(q_i), Product(p_i, p_i)));
            std::get<1>(coefficient) = q_i;

            double l_i = 0.5 * (o_plus_one - o_minus_one);
            result.add(Product(Constant(l_i), p_i));
            std::get<2>(coefficient) = l_i;

            // reset parameter to zero
            p_i = 0.0;
        }

        // Determine the bilinear terms 'b_{ij}'
        for (auto i = coefficients.begin(), i_end = coefficients.end(); i != i_end; ++i)
        {
            Parameter p_i = std::get<0>(*i);
            double    q_i = std::get<1>(*i), l_i = std::get<2>(*i);
            p_i = 1.0;

            auto j = i;
            if (j != i_end)
            {
                ++j;
            }

            for (; j != i_end; ++j)
            {
                Parameter p_j = std::get<0>(*j);
                double    q_j = std::get<1>(*j), l_j = std::get<2>(*j);
                p_j = 1.0;

                // extract bilinear term
                double b_ij = o->evaluate() - n - q_i - l_i - q_j - l_j;

                result.add(Product(Constant(b_ij), Product(p_i, p_j)));

                p_j = 0.0;
            }

            p_i = 0.0;
        }

        // Reset parameters to defaults
        for (const auto & coefficient : coefficients)
        {
            Parameter p_i = std::get<0>(coefficient);

            p_i = p_i.central();
        }

        return result;
    }

    class WilsonPolynomialRatio : public Observable
    {
        private:
            WilsonPolynomial _numerator, _denominator;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            QualifiedName _name;

        public:
            WilsonPolynomialRatio(const WilsonPolynomial & numerator, const WilsonPolynomial & denominator, const Parameters & parameters) :
                _numerator(numerator),
                _denominator(denominator),
                _parameters(parameters),
                _name("WilsonPolynomial::Ratio")
            {
            }

            ~WilsonPolynomialRatio() {}

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual Kinematics
            kinematics()
            {
                return _kinematics;
            }

            virtual Parameters
            parameters()
            {
                return _parameters;
            }

            virtual Options
            options()
            {
                return _options;
            }

            virtual ObservablePtr
            clone() const
            {
                throw InternalError("Cloning WilsonPolynomialRatio without external parameters");
            }

            virtual double
            evaluate() const
            {
                WilsonPolynomialEvaluator evaluator;

                return _numerator.accept_returning<double>(evaluator) / _denominator.accept_returning<double>(evaluator);
            }

            virtual ObservablePtr
            clone(const Parameters & parameters) const
            {
                WilsonPolynomialCloner cloner(parameters);

                return ObservablePtr(
                        new WilsonPolynomialRatio(_numerator.accept_returning<WilsonPolynomial>(cloner), _denominator.accept_returning<WilsonPolynomial>(cloner), parameters));
            }
    };

    class WilsonPolynomialHTLikeRatio : public Observable
    {
        private:
            WilsonPolynomial _numerator, _denominator1, _denominator2;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            QualifiedName _name;

        public:
            WilsonPolynomialHTLikeRatio(const WilsonPolynomial & numerator, const WilsonPolynomial & denominator1, const WilsonPolynomial & denominator2,
                                        const Parameters & parameters) :
                _numerator(numerator),
                _denominator1(denominator1),
                _denominator2(denominator2),
                _parameters(parameters),
                _name("WilsonPolynomial::HTLikeRatio")
            {
            }

            ~WilsonPolynomialHTLikeRatio() {}

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual Kinematics
            kinematics()
            {
                return _kinematics;
            }

            virtual Parameters
            parameters()
            {
                return _parameters;
            }

            virtual Options
            options()
            {
                return _options;
            }

            virtual ObservablePtr
            clone() const
            {
                throw InternalError("Cloning WilsonPolynomialHTLikeRatio without external parameters");
            }

            virtual double
            evaluate() const
            {
                WilsonPolynomialEvaluator evaluator;

                return _numerator.accept_returning<double>(evaluator)
                       / std::sqrt(_denominator1.accept_returning<double>(evaluator) * _denominator2.accept_returning<double>(evaluator));
            }

            virtual ObservablePtr
            clone(const Parameters & parameters) const
            {
                WilsonPolynomialCloner cloner(parameters);

                return ObservablePtr(new WilsonPolynomialHTLikeRatio(_numerator.accept_returning<WilsonPolynomial>(cloner),
                                                                     _denominator1.accept_returning<WilsonPolynomial>(cloner),
                                                                     _denominator2.accept_returning<WilsonPolynomial>(cloner),
                                                                     parameters));
            }
    };

    ObservablePtr
    make_polynomial_observable(const WilsonPolynomial & polynomial, const Parameters & parameters)
    {
        return ObservablePtr(new WilsonPolynomialRatio(polynomial, Constant(1.0), parameters));
    }

    ObservablePtr
    make_polynomial_ratio(const WilsonPolynomial & numerator, const WilsonPolynomial & denominator, const Parameters & parameters)
    {
        return ObservablePtr(new WilsonPolynomialRatio(numerator, denominator, parameters));
    }

    ObservablePtr
    make_polynomial_ht_like_ratio(const WilsonPolynomial & numerator, const WilsonPolynomial & denominator1, const WilsonPolynomial & denominator2, const Parameters & parameters)
    {
        return ObservablePtr(new WilsonPolynomialHTLikeRatio(numerator, denominator1, denominator2, parameters));
    }

    /* WilsonPolynomialCloner */
    WilsonPolynomialCloner::WilsonPolynomialCloner(const Parameters & parameters) :
        _parameters(parameters)
    {
    }

    WilsonPolynomial
    WilsonPolynomialCloner::visit(const Constant & c)
    {
        return WilsonPolynomial(Constant(c));
    }

    WilsonPolynomial
    WilsonPolynomialCloner::visit(const Sum & s)
    {
        Sum result;

        for (const auto & summand : s.summands)
        {
            result.add(summand.accept_returning<WilsonPolynomial>(*this));
        }

        return WilsonPolynomial(result);
    }

    WilsonPolynomial
    WilsonPolynomialCloner::visit(const Product & p)
    {
        return WilsonPolynomial(Product(p.x.accept_returning<WilsonPolynomial>(*this), p.y.accept_returning<WilsonPolynomial>(*this)));
    }

    WilsonPolynomial
    WilsonPolynomialCloner::visit(const Sine & s)
    {
        return WilsonPolynomial(Sine(s.phi.accept_returning<WilsonPolynomial>(*this)));
    }

    WilsonPolynomial
    WilsonPolynomialCloner::visit(const Cosine & c)
    {
        return WilsonPolynomial(Cosine(c.phi.accept_returning<WilsonPolynomial>(*this)));
    }

    WilsonPolynomial
    WilsonPolynomialCloner::visit(const Parameter & p)
    {
        return WilsonPolynomial(_parameters[p.name()]);
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
        {
            result += "\n   ";
        }

        auto i = s.summands.cbegin(), i_end = s.summands.cend();

        if (i != i_end)
        {
            result += (*i).accept_returning<std::string>(*this);
            ++i;
        }

        for (; i != i_end; ++i)
        {
            result += std::string(_pretty ? "\n" : "") + " + " + (*i).accept_returning<std::string>(*this);
        }

        if (_pretty)
        {
            result += "\n";
        }

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
    WilsonPolynomialEvaluator::visit(const Constant & c)
    {
        return c.value;
    }

    double
    WilsonPolynomialEvaluator::visit(const Parameter & p)
    {
        return p();
    }

    double
    WilsonPolynomialEvaluator::visit(const Sum & s)
    {
        double result = 0.0;

        for (const auto & summand : s.summands)
        {
            result += summand.accept_returning<double>(*this);
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
} // namespace eos
