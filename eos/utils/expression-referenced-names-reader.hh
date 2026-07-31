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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_REFERENCED_NAMES_READER_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_REFERENCED_NAMES_READER_HH 1

#include <eos/utils/expression-fwd.hh>
#include <eos/utils/qualified-name.hh>

#include <set>

namespace eos::exp
{
    class ExpressionReferencedNamesReader
    {
        public:
            std::set<QualifiedName> observable_names;
            std::set<QualifiedName> parameter_names;

            ExpressionReferencedNamesReader()  = default;
            ~ExpressionReferencedNamesReader() = default;

            void operator() (const BinaryExpression & e);
            void operator() (const FunctionExpression & e);
            void operator() (const ConstantExpression &);
            void operator() (const ObservableNameExpression & e);
            void operator() (const ObservableExpression & e);
            void operator() (const ParameterNameExpression & e);
            void operator() (const ParameterExpression & e);
            void operator() (const KinematicVariableNameExpression &);
            void operator() (const KinematicVariableExpression &);
            void operator() (const CachedObservableExpression &);
    };
} // namespace eos::exp

#endif
