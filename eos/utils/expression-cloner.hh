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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_CLONER_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_CLONER_HH 1

#include <eos/utils/expression-fwd.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>

namespace eos::exp
{
    // Visit the expression tree and clone its components
    class ExpressionCloner
    {
        private:
            Parameters _parameters;
            Kinematics _kinematics;
            Options    _options;

        public:
            ExpressionCloner(const Parameters & parameters, const Kinematics & kinematics, const Options & options);

            Expression operator() (const BinaryExpression & e);

            Expression operator() (const FunctionExpression & e);

            Expression operator() (const ConstantExpression & e);

            Expression operator() (const ObservableNameExpression & e);

            Expression operator() (const ObservableExpression & e);

            Expression operator() (const ParameterNameExpression & e);

            Expression operator() (const ParameterExpression & e);

            Expression operator() (const KinematicVariableNameExpression & e);

            Expression operator() (const KinematicVariableExpression & e);

            Expression operator() (const CachedObservableExpression & e);
    };
} // namespace eos::exp

#endif
