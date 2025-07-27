/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_PARSER_IMPL_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_PARSER_IMPL_HH 1

#include <eos/utils/expression-parser.hh>
#include <eos/utils/expression.hh>

#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>

#include <utility>

namespace eos
{
    template <typename Iterator>
    ExpressionParser<Iterator>::ExpressionParser() :
        ExpressionParser::base_type(expression)
    {
        using namespace boost::spirit::qi;

        expression = additive_expr[_val = _1];

        auto make_binary = [](char op, const auto & lhs, const auto & rhs) { return eos::exp::BinaryExpression(op, lhs, rhs); };

        additive_expr = multiplicative_expr[_val = _1] >> *((char_("+") >> multiplicative_expr)[_val = phx::bind(make_binary, _1, _val, _2)]
                                                            | (char_("-") >> multiplicative_expr)[_val = phx::bind(make_binary, _1, _val, _2)]);

        multiplicative_expr = exponential_expr[_val = _1] >> *((char_("*") >> exponential_expr)[_val = phx::bind(make_binary, _1, _val, _2)]
                                                               | (char_("/") >> exponential_expr)[_val = phx::bind(make_binary, _1, _val, _2)]);

        exponential_expr = primary_expr[_val = _1] >> -(char_("^") >> primary_expr)[_val = phx::bind(make_binary, _1, _val, _2)];

        auto make_observable = [](const std::string & name, const KinematicsSpecification & spec) { return eos::exp::Expression(eos::exp::ObservableNameExpression(name, spec)); };

        auto make_parameter = [](const std::string & name) { return eos::exp::Expression(eos::exp::ParameterNameExpression(name)); };

        auto make_kinematic_variable = [](const std::string & name) { return eos::exp::Expression(eos::exp::KinematicVariableNameExpression(name)); };

        primary_expr = ('(' >> expression >> ')')[_val = _1] | constant[_val = _1] | (observable_name >> kinematics)[_val = phx::bind(make_observable, _1, _2)]
                       | observable_name[_val = phx::bind(make_observable, _1, KinematicsSpecification())] | parameter_name[_val = phx::bind(make_parameter, _1)]
                       | kinematic_variable_name[_val = phx::bind(make_kinematic_variable, _1)] | function_expr[_val = _1];

        constant = double_ | int_;

        kinematic_variable_name = "{" >> as_string[lexeme[*(char_ - "}")]][_val = _1] >> "}";

        observable_name = "<<" >> as_string[lexeme[*(char_ - ">>")]][_val = _1] >> ">>";

        parameter_name = "[[" >> as_string[lexeme[*(char_ - "]]")]][_val = _1] >> "]]";

        kinematics = '[' >> (kinematics_alias[phx::bind(_val, _1)] | kinematics_value[phx::bind(_val, _1)]) % ',' >> ']';

        kinematics_alias = (as_string[lexeme[*~char_(",=>]")]] >> "=>" >> as_string[lexeme[*~char_(",=]")]]);

        kinematics_value = (as_string[lexeme[*~char_(",=>]")]] >> '=' >> double_);

        auto make_function = [](const std::string & f, const auto & arg) { return eos::exp::FunctionExpression(f, arg); };

        function_expr = (function_name >> '(' >> primary_expr >> ')')[_val = phx::bind(make_function, _1, _2)];

        function_name = *(string("exp") | string("cos") | string("sin"));
    }

    template <typename Iterator> ExpressionParser<Iterator>::~ExpressionParser() {}
} // namespace eos

#endif
