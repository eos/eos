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

#include <eos/observable.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <cmath>

namespace eos
{
    /* Build a Wilson polynomial as an exp::Expression from an ObservablePtr */
    exp::Expression
    make_polynomial(const ObservablePtr & o, const std::list<std::string> & _coefficients)
    {
        exp::Expression result;

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
        result   = exp::ConstantExpression(n);

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

            exp::Expression p_i_expr         = exp::ParameterExpression(p_i);
            exp::Expression q_i_expr         = exp::ConstantExpression(q_i);
            exp::Expression p_i_squared_expr = exp::BinaryExpression('*', p_i_expr, p_i_expr);
            result                           = exp::BinaryExpression('+', result, exp::BinaryExpression('*', q_i_expr, p_i_squared_expr));
            std::get<1>(coefficient)         = q_i;

            double          l_i      = 0.5 * (o_plus_one - o_minus_one);
            exp::Expression l_i_expr = exp::ConstantExpression(l_i);
            result                   = exp::BinaryExpression('+', result, exp::BinaryExpression('*', l_i_expr, p_i_expr));
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

                exp::Expression p_i_expr  = exp::ParameterExpression(p_i);
                exp::Expression p_j_expr  = exp::ParameterExpression(p_j);
                exp::Expression b_ij_expr = exp::ConstantExpression(b_ij);
                exp::Expression p_ij_expr = exp::BinaryExpression('*', p_i_expr, p_j_expr);
                result                    = exp::BinaryExpression('+', result, exp::BinaryExpression('*', b_ij_expr, p_ij_expr));

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

    exp::Expression
    make_polynomial_ratio(const ObservablePtr & numerator, const ObservablePtr & denominator, const std::list<std::string> & coefficients)
    {
        const auto numerator_exp   = make_polynomial(numerator, coefficients);
        const auto denominator_exp = make_polynomial(denominator, coefficients);

        return exp::BinaryExpression('/', numerator_exp, denominator_exp);
    }

} // namespace eos
