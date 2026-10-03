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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_KERNEL_DECOMPOSER_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_KERNEL_DECOMPOSER_HH 1

#include <eos/utils/expression.hh>

#include <string>
#include <vector>

namespace eos::exp
{
    // One term coefficient * kernel(u, parameters...) of a linear combination of kernels in the offset u.
    struct KernelTerm
    {
            ExpressionPtr              coefficient;
            std::string                kernel;
            std::vector<ExpressionPtr> parameters;
    };

    /*!
     * Visit the expression tree and decompose it into a linear combination of kernels in the
     * offset variable. The coefficients and the kernel parameters are sub-expressions that do not
     * depend on the offset variable, and may include kernels in other variables. A kernel in the
     * offset variable must take that variable itself as its first argument. Throws ExpressionError
     * if the expression is not of this form.
     */
    class ExpressionKernelDecomposer
    {
        public:
            // The terms found below the visited node; empty if it does not depend on the offset variable.
            using Result = std::vector<KernelTerm>;

            explicit ExpressionKernelDecomposer(const std::string & offset_variable);
            ~ExpressionKernelDecomposer() = default;

            // Decompose the expression; throws ExpressionError if it contains no kernel.
            std::vector<KernelTerm> decompose(const ExpressionPtr & e);

            Result operator() (const BinaryExpression & e);

            Result operator() (const FunctionExpression & e);

            Result operator() (const ConstantExpression &);

            Result operator() (const ObservableNameExpression & e);

            Result operator() (const ObservableExpression & e);

            Result operator() (const ParameterNameExpression &);

            Result operator() (const ParameterExpression &);

            Result operator() (const KinematicVariableNameExpression & e);

            Result operator() (const KinematicVariableExpression & e);

            Result operator() (const CachedObservableExpression & e);

        private:
            std::string _offset_variable;

            // Visit a sub-expression that must not depend on the offset variable.
            void _independent(const ExpressionPtr & e, const std::string & context);

            void _check_aliases(const KinematicsSpecification & spec);
    };
} // namespace eos::exp

#endif
