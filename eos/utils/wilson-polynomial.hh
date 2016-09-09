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

#ifndef EOS_GUARD_SRC_UTILS_WILSON_POLYNOMIAL_HH
#define EOS_GUARD_SRC_UTILS_WILSON_POLYNOMIAL_HH 1

#include <eos/observable.hh>
#include <eos/utils/one-of.hh>

#include <list>
#include <string>

namespace eos
{
    struct Constant;

    struct Sum;

    struct Product;

    struct Sine;

    struct Cosine;

    typedef OneOf<Constant, Sum, Product, Sine, Cosine, Parameter> WilsonPolynomial;

    WilsonPolynomial make_polynomial(const ObservablePtr &, const std::list<std::string> &);

    /*!
     * Return an Observable that wraps a WilsonPolynomial object.
     *
     * @param polynomial  The polynomial that shall be wrapped.
     * @param parameters  The Parameters object of the polynomial.
     */
    ObservablePtr make_polynomial_observable(const WilsonPolynomial & polynomial,
            const Parameters & parameters);

    /*!
     * Return an Observable that is a ratio of two WilsonPolynomial objects.
     *
     * @param numerator   The numerator of the ratio.
     * @param denominator The denominator of the ratio.
     * @param parameters  The common Parameters object of both numerator and denominator.
     */
    ObservablePtr make_polynomial_ratio(const WilsonPolynomial & numerator, const WilsonPolynomial & denominator,
            const Parameters & parameters);

    /*!
     * Return an Observable that is a ratio similar to H_T^(i) of three WilsonPolynomial objects N,D1,D2:
     * @f[N / \sqrt{D1 \cdot D2}@f]
     *
     * @param numerator    The numerator of the ratio.
     * @param denominator1 The first denominator component of the ratio.
     * @param denominator2 The first denominator component of the ratio.
     * @param parameters   The common Parameters object of both numerator and denominator.
     */
    ObservablePtr make_polynomial_ht_like_ratio(const WilsonPolynomial & numerator, const WilsonPolynomial & denominator,
            const WilsonPolynomial & denominator2, const Parameters & parameters);

    class WilsonPolynomialCloner
    {
        private:
            Parameters _parameters;

        public:
            WilsonPolynomialCloner(const Parameters & parameters);

            WilsonPolynomial visit(const Constant & c);
            WilsonPolynomial visit(const Sum & s);
            WilsonPolynomial visit(const Product & p);
            WilsonPolynomial visit(const Sine & s);
            WilsonPolynomial visit(const Cosine & s);
            WilsonPolynomial visit(const Parameter & p);
    };

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
