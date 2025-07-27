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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_PRINTER_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_PRINTER_HH 1

#include <eos/utils/expression-fwd.hh>

#include <ostream>

namespace eos::exp
{
    // Visit the expression tree and print its components
    class ExpressionPrinter
    {
        private:
            std::ostream & _os;

        public:
            ExpressionPrinter(std::ostream & ostream);

            void visit(BinaryExpression & e);

            void visit(FunctionExpression & e);

            void visit(ConstantExpression & e);

            void visit(ObservableNameExpression & e);

            void visit(ObservableExpression & e);

            void visit(ParameterNameExpression & e);

            void visit(ParameterExpression & e);

            void visit(KinematicVariableNameExpression & e);

            void visit(KinematicVariableExpression & e);

            void visit(CachedObservableExpression & e);
    };
} // namespace eos::exp

#endif
