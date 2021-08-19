/*
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/utils/exception.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/join.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <iostream>
#include <set>

namespace eos::exp
{
    // Visit the expression tree and print its components
    class ExpressionPrinter
    {
        private:
            std::ostream & _os;

        public:
            ExpressionPrinter(std::ostream& ostream) : _os(ostream)
            {
            }

            void visit(BinaryExpression & e)
            {
                _os << "BinaryExpression(";
                e.lhs.accept(*this);
                _os << " " << e.op << " ";
                e.rhs.accept(*this);
                _os << ")";
            }

            void visit(ConstantExpression & e)
            {
                _os << "ConstantExpression(" << e.value << ")";
            }

            void visit(ObservableNameExpression & e)
            {
                _os << "ObservableNameExpression(" << e.observable_name.full();
                if (! e.kinematics_specification.aliases.empty())
                {
                    auto a = e.kinematics_specification.aliases.cbegin();
                    _os << ", aliases=[" << a->first << "=>" << a->second;
                    ++a;
                    for (auto a_end = e.kinematics_specification.aliases.cend() ; a != a_end ; ++a)
                    {
                        _os << "," << a->first << "=>" << a->second;
                    }
                    _os << "]";
                }
                if (! e.kinematics_specification.values.empty())
                {
                    auto v = e.kinematics_specification.values.cbegin();
                    _os << ", values=[" << v->first << "=" << v->second;
                    ++v;
                    for (auto v_end = e.kinematics_specification.values.cend() ; v != v_end ; ++v)
                    {
                        _os << "," << v->first << "=" << v->second;
                    }
                    _os << "]";
                }
                _os << ")";
            }

            void visit(ObservableExpression & e)
            {
                _os << "ObservableExpression(" << e.observable->name() << ")";
                if (! e.kinematics_specification.aliases.empty())
                {
                    auto a = e.kinematics_specification.aliases.cbegin();
                    _os << ", aliases=[" << a->first << "=>" << a->second;
                    ++a;
                    for (auto a_end = e.kinematics_specification.aliases.cend() ; a != a_end ; ++a)
                    {
                        _os << "," << a->first << "=>" << a->second;
                    }
                    _os << "]";
                }
                if (! e.kinematics_specification.values.empty())
                {
                    auto v = e.kinematics_specification.values.cbegin();
                    _os << ", values=[" << v->first << "=" << v->second;
                    ++v;
                    for (auto v_end = e.kinematics_specification.values.cend() ; v != v_end ; ++v)
                    {
                        _os << "," << v->first << "=" << v->second;
                    }
                    _os << "]";
                }
                _os << ")";
            }
    };

    // Evaluate the expression tree
    class ExpressionEvaluator
    {
        public:
            double visit(BinaryExpression & e)
            {
                BinaryExpression::func f = BinaryExpression::Method(e.op);

                return f(e.lhs.accept_returning<double>(*this), e.rhs.accept_returning<double>(*this));
            }

            double visit(ConstantExpression & e)
            {
                return e.value;
            }

            double visit(ObservableExpression & e)
            {
                return e.observable->evaluate();
            }

            double visit(ObservableNameExpression &)
            {
                throw InternalError("Encountered ObserableNameExpression in ExpressionEvaluator::visit");

                return 0.0;
            }
    };

    // Visit the expression tree and clone its components
    class ExpressionCloner
    {
        private:
            Parameters _parameters;
            Kinematics _kinematics;
            Options    _options;

        public:
            ExpressionCloner(const Parameters & parameters, const Kinematics & kinematics, const Options & options) :
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options)
            {
            }

            Expression visit(const BinaryExpression & e)
            {
                return BinaryExpression(e.op, e.lhs.accept_returning<Expression>(*this), e.rhs.accept_returning<Expression>(*this));
            }

            Expression visit(const ConstantExpression & e)
            {
                return ConstantExpression(e);
            }

            Expression visit(const ObservableNameExpression & e)
            {
                return ObservableNameExpression(e);
            }

            Expression visit(const ObservableExpression & e)
            {
                const auto & kinematics_values = e.kinematics_specification.values;
                const auto & kinematics_aliases = e.kinematics_specification.aliases;

                // Set or alias kinematic specifications
                for (const auto & value : kinematics_values)
                {
                    _kinematics.declare(value.first, value.second);
                }
                for (const auto & alias : kinematics_aliases)
                {
                    _kinematics.alias(alias.first, alias.second);
                }

                // Make observable
                auto observable = Observable::make(
                    e.observable->name(),
                    this->_parameters,
                    this->_kinematics,
                    this->_options
                    );
                // Clear alias map
                _kinematics.clear_aliases();

                return ObservableExpression(observable, e.kinematics_specification);
            }
    };

    // Copies one expression tree to a new tree, replacing ObservableNameExpression objects with ObservableExpression objects
    class ExpressionMaker
    {
        private:
            Parameters      _parameters;
            Kinematics      _kinematics;
            Options         _options;
            ParameterUser * _parameter_user;

        public:
            ExpressionMaker(const Parameters & parameters, const Kinematics & kinematics, const Options & options, ParameterUser * parameter_user = nullptr) :
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _parameter_user(parameter_user)
            {
            }

            Expression visit(const BinaryExpression & e)
            {
                return BinaryExpression(e.op, e.lhs.accept_returning<Expression>(*this), e.rhs.accept_returning<Expression>(*this));
            }

            Expression visit(const ConstantExpression & e)
            {
                return ConstantExpression(e);
            }

            Expression visit(const ObservableNameExpression & e)
            {
                const auto & kinematics_values = e.kinematics_specification.values;
                const auto & kinematics_aliases = e.kinematics_specification.aliases;

                // Set or alias kinematic specifications
                for (const auto & value : kinematics_values)
                {
                    _kinematics.declare(value.first, value.second);
                }
                for (const auto & alias : kinematics_aliases)
                {
                    _kinematics.alias(alias.first, alias.second);
                }

                // Make observable
                auto observable = Observable::make(
                    e.observable_name,
                    this->_parameters,
                    this->_kinematics,
                    this->_options
                    );
                // Clear alias map
                _kinematics.clear_aliases();

                // Record the used parameters
                if (_parameter_user)
                {
                    _parameter_user->uses(*observable);
                }

                return ObservableExpression(observable, e.kinematics_specification);
            }

            Expression visit(const ObservableExpression & e)
            {
                throw InternalError("Encountered ObserableExpression in ExpressionMaker::visit");

                return e;
            }
    };

    // Visit the expression tree and return the set of kinematics names used by its components
    class ExpressionKinematicReader
    {
        private:

        public:
            ExpressionKinematicReader()
            {
            }

            std::set<std::string> visit(const BinaryExpression & e)
            {
                std::set<std::string> lhs_set = std::move(e.lhs.accept_returning<std::set<std::string>>(*this));
                std::set<std::string> rhs_set = std::move(e.rhs.accept_returning<std::set<std::string>>(*this));
                lhs_set.insert(rhs_set.begin(), rhs_set.end());

                return lhs_set;
            }

            std::set<std::string> visit(const ConstantExpression &)
            {
                std::set<std::string> empty_set;

                return empty_set;
            }

            std::set<std::string> visit(const ObservableNameExpression & e)
            {
                std::set<std::string> kinematic_set;

                const auto & kinematics_values = e.kinematics_specification.values;
                const auto & kinematics_aliases = e.kinematics_specification.aliases;

                // Set or alias kinematic specifications
                for (const auto & value : kinematics_values)
                {
                    kinematic_set.insert(value.first);
                }
                for (const auto & alias : kinematics_aliases)
                {
                    kinematic_set.insert(alias.second);
                }

                return kinematic_set;
            }

            std::set<std::string> visit(const ObservableExpression & e)
            {
                std::set<std::string> kinematic_set;

                const auto & kinematics_values = e.kinematics_specification.values;
                const auto & kinematics_aliases = e.kinematics_specification.aliases;

                // Set or alias kinematic specifications
                for (const auto & value : kinematics_values)
                {
                    kinematic_set.insert(value.first);
                }
                for (const auto & alias : kinematics_aliases)
                {
                    kinematic_set.insert(alias.second);
                }

                return kinematic_set;
            }
    };
}
