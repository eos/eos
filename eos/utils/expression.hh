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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_HH 1

#include <eos/observable-fwd.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/observable_cache.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <cassert>
#include <memory>
#include <iostream>
#include <map>

namespace eos::exp
{
    class ExpressionError :
        public eos::Exception
    {
        public:
            ExpressionError(const std::string & msg);
    };

    class BinaryExpression
    {
        public:
            char op;
            Expression lhs, rhs;

            using func = double(*)(const double &, const double &);

            static double sum(const double &, const double &);
            static double difference(const double &, const double &);
            static double product(const double &, const double &);
            static double ratio(const double &, const double &);
            static double power(const double &, const double &);

            static BinaryExpression::func Method(char op);

            BinaryExpression() {}

            BinaryExpression(char op, const Expression & l, const Expression & r) :
                op(op), lhs(l), rhs(r)
            {
            }
    };

    class FunctionExpression
    {
        public:
            using FunctionType = double(*)(const double &);

            FunctionType f;
            std::string fname;
            Expression arg;

            FunctionExpression() {};

            FunctionExpression(const std::string & f, const Expression & arg);
    };

    class ConstantExpression
    {
        public:
            double value;

            ConstantExpression(const double& v = 0) :
                value(v)
            {
            }
    };

    class KinematicsSpecification
    {
        public:
            std::map<std::string, double> values;
            std::map<std::string, std::string> aliases;

            void operator() (const std::pair<std::string, double> & value) { values.insert(value); }
            void operator() (const std::pair<std::string, std::string> & alias) { aliases.insert(alias); }
    };

    class KinematicVariableNameExpression
    {
        public:
            std::string variable_name;

            KinematicVariableNameExpression(const std::string & variable_name) :
                 variable_name(variable_name)
            {
            }
    };

    class KinematicVariableExpression
    {
        public:
            KinematicVariable kinematic_variable;

            KinematicVariableExpression(const KinematicVariable & kinematic_variable) :
                kinematic_variable(kinematic_variable)
            {
            }
    };

    class ObservableNameExpression
    {
        public:
            QualifiedName observable_name;
            KinematicsSpecification kinematics_specification;

            ObservableNameExpression(const QualifiedName & observable_name, const KinematicsSpecification & kinematics_specification) :
                 observable_name(observable_name),
                 kinematics_specification(kinematics_specification)
            {
            }
    };

    class ObservableExpression
    {
        public:
            ObservablePtr observable;
            KinematicsSpecification kinematics_specification;

            ObservableExpression(ObservablePtr observable, const KinematicsSpecification & kinematics_specification) :
                 observable(observable),
                 kinematics_specification(kinematics_specification)
            {
            }
    };

    class CachedObservableExpression
    {
        public:
            ObservableCache cache;
            ObservableCache::Id id;
            KinematicsSpecification kinematics_specification;

            CachedObservableExpression(const ObservableCache & cache, const ObservableCache::Id & id, const KinematicsSpecification & kinematics_specification) :
                cache(cache),
                id(id),
                kinematics_specification(kinematics_specification)
            {
            }
    };

    class ParameterNameExpression
    {
        public:
            QualifiedName parameter_name;

            ParameterNameExpression(const QualifiedName & parameter_name) :
                parameter_name(parameter_name)
            {
            }
    };

    class ParameterExpression
    {
        public:
            Parameter parameter;

            ParameterExpression(const Parameter & parameter) :
                parameter(parameter)
            {
            }
    };
}

#endif
