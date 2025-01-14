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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_PARSER_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_PARSER_HH 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/adapted.hpp>

#include <eos/utils/expression.hh>

namespace eos
{
    namespace qi    = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx   = boost::phoenix;

    template <typename Iterator>
    struct ExpressionParser :
        qi::grammar<Iterator, eos::exp::Expression(), ascii::space_type>
    {
        // Constructor
        ExpressionParser();

        qi::rule<Iterator, eos::exp::Expression()               , ascii::space_type> expression;
        qi::rule<Iterator, eos::exp::Expression()               , ascii::space_type> additive_expr;
        qi::rule<Iterator, eos::exp::Expression()               , ascii::space_type> multiplicative_expr;
        qi::rule<Iterator, eos::exp::Expression()               , ascii::space_type> exponential_expr;
        qi::rule<Iterator, eos::exp::Expression()               , ascii::space_type> function_expr;

        qi::rule<Iterator, eos::exp::Expression()               , ascii::space_type> primary_expr;
        qi::rule<Iterator, eos::exp::ConstantExpression()       , ascii::space_type> constant;
        qi::rule<Iterator, std::string()                        , ascii::space_type> observable_name;
        qi::rule<Iterator, std::string()                        , ascii::space_type> parameter_name;
        qi::rule<Iterator, std::string()                        , ascii::space_type> kinematic_variable_name;
        qi::rule<Iterator, std::string()                        , ascii::space_type> function_name;

        using KinematicsSpecification = eos::exp::KinematicsSpecification;

        qi::rule<Iterator, KinematicsSpecification()            , ascii::space_type> kinematics;
        qi::rule<Iterator, std::pair<std::string, std::string>(), ascii::space_type> kinematics_alias;
        qi::rule<Iterator, std::pair<std::string, double>       , ascii::space_type> kinematics_value;

        // Destuctor
        ~ExpressionParser();
    };

    extern template struct ExpressionParser<std::string::const_iterator>;
}

#endif
