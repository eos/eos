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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_EVALUATOR_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_EVALUATOR_HH 1

#include <eos/utils/expression-fwd.hh>

namespace eos::exp
{
    // Evaluate the expression tree
    class ExpressionEvaluator
    {
        public:
            ExpressionEvaluator() = default;
            ~ExpressionEvaluator() = default;

            double visit(BinaryExpression & e);

            double visit(FunctionExpression & e);

            double visit(ConstantExpression & e);

            double visit(ObservableNameExpression &);

            double visit(ObservableExpression & e);

            double visit(ParameterNameExpression &);

            double visit(ParameterExpression & e);

            double visit(KinematicVariableNameExpression & e);

            double visit(KinematicVariableExpression & e);

            double visit(CachedObservableExpression & e);

    };
}

#endif
