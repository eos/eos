/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-kinematic-hook.hh>
#include <eos/utils/expression-used-kinematics-reader.hh>
#include <eos/utils/expression.hh>

#include <array>

namespace eos::exp
{
    namespace
    {
        // Copies an expression tree, binding its kinematic variables to the given kinematics.
        class KinematicVariableBinder
        {
            private:
                Kinematics _kinematics;

            public:
                KinematicVariableBinder(const Kinematics & kinematics) :
                    _kinematics(kinematics)
                {
                }

                Expression
                operator() (const BinaryExpression & e)
                {
                    return BinaryExpression(e.op, ExpressionPtr(new Expression(std::visit(*this, *e.lhs))), ExpressionPtr(new Expression(std::visit(*this, *e.rhs))));
                }

                Expression
                operator() (const FunctionExpression & e)
                {
                    std::array<ExpressionPtr, FunctionExpression::max_arguments> args;
                    for (std::size_t i = 0; i < e.number_of_arguments; ++i)
                    {
                        args[i] = ExpressionPtr(new Expression(std::visit(*this, *e.args[i])));
                    }

                    return FunctionExpression(e.fname, std::span<const ExpressionPtr>(args.data(), e.number_of_arguments));
                }

                Expression
                operator() (const ConstantExpression & e)
                {
                    return ConstantExpression(e);
                }

                Expression
                operator() (const KinematicVariableNameExpression & e)
                {
                    return KinematicVariableExpression(_kinematics[e.variable_name]);
                }

                Expression
                operator() (const KinematicVariableExpression & e)
                {
                    // the variable is bound already, so that an alias of the same name must not capture it
                    return KinematicVariableExpression(_kinematics.variable(e.kinematic_variable.name()));
                }

                // Parameters would not follow a clone of the parameters, and observables would not follow the cache.
                template <typename T_>
                Expression
                operator() (const T_ &)
                {
                    throw ExpressionError("a derived kinematic variable may only depend on kinematic variables, constants, and functions");
                }
        };
    } // namespace

    ExpressionKinematicHook::ExpressionKinematicHook(const ExpressionPtr & expression, const Kinematics & kinematics)
    {
        KinematicVariableBinder binder(kinematics);
        _expression = ExpressionPtr(new Expression(std::visit(binder, *expression)));

        ExpressionUsedKinematicsReader reader;
        std::visit(reader, *_expression);
        for (const auto & id : reader.kinematic_variable_ids)
        {
            uses_kinematic(id);
        }
    }

    double
    ExpressionKinematicHook::evaluate() const
    {
        ExpressionEvaluator evaluator;

        return std::visit(evaluator, *_expression);
    }

    std::shared_ptr<const KinematicHook>
    ExpressionKinematicHook::clone(const Kinematics & kinematics) const
    {
        return std::make_shared<ExpressionKinematicHook>(_expression, kinematics);
    }
} // namespace eos::exp
