/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
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

#include <eos/maths/complex.hh>
#include <eos/utils/expression-cloner.hh>
#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-printer.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using namespace test;
using namespace eos;

struct WilsonPolynomialTestObservable : public Observable
{
        QualifiedName n;
        Parameters    p;
        Kinematics    k;
        Parameter     c1;
        Parameter     c2;
        Parameter     re_c7;
        Parameter     im_c7;
        Parameter     re_c9;
        Parameter     im_c9;
        Parameter     re_c10;
        Parameter     im_c10;

        WilsonPolynomialTestObservable(const Parameters & p, const Kinematics & k, const Options &) :
            n("WilsonPolynomial::TestObservable"),
            p(p),
            k(k),
            c1(p["b->s::c1"]),
            c2(p["b->s::c2"]),
            re_c7(p["b->s::Re{c7}"]),
            im_c7(p["b->s::Im{c7}"]),
            re_c9(p["b->smumu::Re{c9}"]),
            im_c9(p["b->smumu::Im{c9}"]),
            re_c10(p["b->smumu::Re{c10}"]),
            im_c10(p["b->smumu::Im{c10}"])
        {
        }

        virtual const QualifiedName &
        name() const
        {
            return n;
        }

        virtual Parameters
        parameters()
        {
            return p;
        }

        virtual Kinematics
        kinematics()
        {
            return k;
        }

        virtual Options
        options()
        {
            return Options();
        }

        virtual ObservablePtr
        clone() const
        {
            return ObservablePtr(new WilsonPolynomialTestObservable(p.clone(), k.clone(), Options()));
        }

        virtual ObservablePtr
        clone(const Parameters & p) const
        {
            return ObservablePtr(new WilsonPolynomialTestObservable(p, k.clone(), Options()));
        }

        virtual double
        evaluate() const
        {
            complex<double> c7(re_c7(), im_c7());
            complex<double> c9(re_c9(), im_c9());
            complex<double> c10(re_c10(), im_c10());

            return real(0.01234 + c7 * complex<double>(0.321, 1.000) + c9 * complex<double>(0.731, 1.000) + conj(c7) * c7 * 0.6 + conj(c7) * c9 * complex<double>(1.300, 0.123)
                        + conj(c9) * c9 * 2.1 + conj(c10) * c10 * 1.23);
        }
};

class WilsonPolynomialTest : public TestCase
{
    public:
        WilsonPolynomialTest() :
            TestCase("wilson_polynomial_test")
        {
        }

        void
        run_one(const ObservablePtr & o, const exp::Expression & p, const std::array<double, 6> & values) const
        {
            Parameters parameters = o->parameters();
            Parameter  re_c7(parameters["b->s::Re{c7}"]);
            Parameter  im_c7(parameters["b->s::Im{c7}"]);
            Parameter  re_c9(parameters["b->smumu::Re{c9}"]);
            Parameter  im_c9(parameters["b->smumu::Im{c9}"]);
            Parameter  re_c10(parameters["b->smumu::Re{c10}"]);
            Parameter  im_c10(parameters["b->smumu::Im{c10}"]);

            re_c7  = values[0];
            im_c7  = values[1];
            re_c9  = values[2];
            im_c9  = values[3];
            re_c10 = values[4];
            im_c10 = values[5];

            static const double      eps = 1e-10;
            exp::ExpressionEvaluator evaluator;
            TEST_CHECK_NEARLY_EQUAL(o->evaluate(), std::visit(evaluator, p), eps);
        }

        virtual void
        run() const
        {
            Parameters parameters = Parameters::Defaults();
            Kinematics kinematics;

            auto o = ObservablePtr(new WilsonPolynomialTestObservable(parameters, kinematics, Options()));
            auto p = make_polynomial(o, std::list<std::string>{ "b->s::Re{c7}", "b->s::Im{c7}", "b->smumu::Re{c9}", "b->smumu::Im{c9}", "b->smumu::Re{c10}", "b->smumu::Im{c10}" });

            exp::ExpressionPrinter printer(std::cout);
            std::visit(printer, p);

            static const std::vector<std::array<double, 6>> inputs{
                std::array<double, 6>{ { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                std::array<double, 6>{ { 1.0, 0.0, 1.0, 0.0, 1.0, 0.0 } },
                std::array<double, 6>{ { 0.7808414, 0.8487257, 0.7735165, 0.5383695, 0.6649164, 0.7235497 } },
                std::array<double, 6>{ { 0.5860642, 0.9830907, 0.7644369, 0.8330194, 0.4935018, 0.4492084 } },
                std::array<double, 6>{ { 0.2177456, 0.5062894, 0.6463376, 0.3624364, 0.6770480, 0.0718421 } },
                std::array<double, 6>{ { 0.0088306, 0.9441413, 0.8721501, 0.2984633, 0.2961408, 0.9145809 } },
                std::array<double, 6>{ { 0.7967655, 0.2427081, 0.8403112, 0.3351082, 0.6477823, 0.5569495 } },
                std::array<double, 6>{ { 0.7607454, 0.5025871, 0.5877762, 0.5516025, 0.2930899, 0.4882813 } },
            };

            for (const auto & input : inputs)
            {
                run_one(o, p, input);
            }
        }
} wilson_polynomial_test;

class WilsonPolynomialClonerTest : public TestCase
{
    public:
        WilsonPolynomialClonerTest() :
            TestCase("wilson_polynomial_cloner_test")
        {
        }

        virtual void
        run() const
        {
            Parameters parameters = Parameters::Defaults();
            Kinematics kinematics;

            auto o = ObservablePtr(new WilsonPolynomialTestObservable(parameters, kinematics, Options()));
            auto p = make_polynomial(o, std::list<std::string>{ "b->s::Re{c7}", "b->smumu::Re{c9}", "b->smumu::Re{c10}" });

            Parameters            clone_parameters = Parameters::Defaults();
            exp::ExpressionCloner cloner(clone_parameters, kinematics, Options());
            exp::Expression       c = std::visit(cloner, p);

            std::string rep_original, rep_clone;
            TEST_CHECK_NO_THROW({
                std::stringstream      output;
                exp::ExpressionPrinter printer(output);
                std::visit(printer, p);
                rep_original = output.str();
            });
            TEST_CHECK_NO_THROW({
                std::stringstream      output;
                exp::ExpressionPrinter printer(output);
                std::visit(printer, c);
                rep_clone = output.str();
            });
            TEST_CHECK_EQUAL(rep_original, rep_clone);

            exp::ExpressionEvaluator evaluator;
            TEST_CHECK_EQUAL(std::visit(evaluator, p), std::visit(evaluator, c));

            parameters["b->smumu::Re{c10}"] = 10;
            TEST_CHECK(std::visit(evaluator, p) != std::visit(evaluator, c));

            clone_parameters["b->smumu::Re{c10}"] = 10;
            TEST_CHECK_EQUAL(std::visit(evaluator, p), std::visit(evaluator, c));
        }
} wilson_polynomial_cloner_test;
