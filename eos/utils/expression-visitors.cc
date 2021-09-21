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
#include <eos/utils/expression-cacher.hh>
#include <eos/utils/expression-cloner.hh>
#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-kinematic-reader.hh>
#include <eos/utils/expression-maker.hh>
#include <eos/utils/expression-printer.hh>
#include <eos/utils/join.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>

#include <iostream>
#include <set>

namespace eos::exp
{
    /*
     * ExpressionPrinter
     */

    ExpressionPrinter::ExpressionPrinter(std::ostream & ostream) :
        _os(ostream)
    {
    }

    void ExpressionPrinter::visit(BinaryExpression & e)
    {
        _os << "BinaryExpression(";
        e.lhs.accept(*this);
        _os << " " << e.op << " ";
        e.rhs.accept(*this);
        _os << ")";
    }

    void ExpressionPrinter::visit(ConstantExpression & e)
    {
        _os << "ConstantExpression(" << e.value << ")";
    }

    void ExpressionPrinter::visit(ObservableNameExpression & e)
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

    void ExpressionPrinter::visit(ObservableExpression & e)
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

    void ExpressionPrinter::visit(CachedObservableExpression & e)
    {
        _os << "CachedObservableExpression(id=" << e.id;
        _os << ", name='" << e.cache.observable(e.id)->name().full() << "'";
        _os << ")";
    }

    /*
     * ExpressionEvaluator
     */

    double
    ExpressionEvaluator::visit(BinaryExpression & e)
    {
        BinaryExpression::func f = BinaryExpression::Method(e.op);

        return f(e.lhs.accept_returning<double>(*this), e.rhs.accept_returning<double>(*this));
    }

    double
    ExpressionEvaluator::visit(ConstantExpression & e)
    {
        return e.value;
    }

    double
    ExpressionEvaluator::visit(ObservableNameExpression &)
    {
        throw InternalError("Encountered ObserableNameExpression in ExpressionEvaluator::visit");

        return 0.0;
    }

    double
    ExpressionEvaluator::visit(ObservableExpression & e)
    {
        return e.observable->evaluate();
    }

    double
    ExpressionEvaluator::visit(CachedObservableExpression & e)
    {
        return e.cache[e.id];
    }

    /*
     * ExpressionCloner
     */

    ExpressionCloner::ExpressionCloner(const Parameters & parameters, const Kinematics & kinematics, const Options & options) :
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options)
    {
    }

    Expression
    ExpressionCloner::visit(const BinaryExpression & e)
    {
        return BinaryExpression(e.op, e.lhs.accept_returning<Expression>(*this), e.rhs.accept_returning<Expression>(*this));
    }

    Expression
    ExpressionCloner::visit(const ConstantExpression & e)
    {
        return ConstantExpression(e);
    }

    Expression
    ExpressionCloner::visit(const ObservableNameExpression & e)
    {
        return ObservableNameExpression(e);
    }

    Expression
    ExpressionCloner::visit(const ObservableExpression & e)
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
            this->_options + e.observable->options()
            );
        // Clear alias map
        _kinematics.clear_aliases();

        return ObservableExpression(observable, e.kinematics_specification);
    }

    Expression
    ExpressionCloner::visit(const CachedObservableExpression & e)
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
            e.cache.observable(e.id)->name(),
            this->_parameters,
            this->_kinematics,
            this->_options + e.cache.observable(e.id)->options()
            );
        // Clear alias map
        _kinematics.clear_aliases();

        return ObservableExpression(observable, e.kinematics_specification);
    }

    /*
     * ExpressionMaker
     */

    ExpressionMaker::ExpressionMaker(const Parameters & parameters, const Kinematics & kinematics, const Options & options, ParameterUser * parameter_user) :
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options),
        _parameter_user(parameter_user)
    {
    }

    Expression
    ExpressionMaker::visit(const BinaryExpression & e)
    {
        return BinaryExpression(e.op, e.lhs.accept_returning<Expression>(*this), e.rhs.accept_returning<Expression>(*this));
    }

    Expression
    ExpressionMaker::visit(const ConstantExpression & e)
    {
        return ConstantExpression(e);
    }

    Expression
    ExpressionMaker::visit(const ObservableNameExpression & e)
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

    Expression
    ExpressionMaker::visit(const ObservableExpression & e)
    {
        // Rebuild the observable expression with the local set of parameters

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
            this->_options + e.observable->options()
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

    Expression
    ExpressionMaker::visit(const CachedObservableExpression & e)
    {
        throw InternalError("Encountered CachedObservableExpression in ExpressionMaker::visit");

        return e;
    }

    /*
     * ExpressionKinematicReader
     */

    std::set<std::string>
    ExpressionKinematicReader::visit(const BinaryExpression & e)
    {
        std::set<std::string> lhs_set = std::move(e.lhs.accept_returning<std::set<std::string>>(*this));
        std::set<std::string> rhs_set = std::move(e.rhs.accept_returning<std::set<std::string>>(*this));
        lhs_set.insert(rhs_set.begin(), rhs_set.end());

        return lhs_set;
    }

    std::set<std::string>
    ExpressionKinematicReader::visit(const ConstantExpression &)
    {
        std::set<std::string> empty_set;

        return empty_set;
    }

    std::set<std::string>
    ExpressionKinematicReader::visit(const ObservableNameExpression & e)
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

    std::set<std::string>
    ExpressionKinematicReader::visit(const ObservableExpression & e)
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

    std::set<std::string>
    ExpressionKinematicReader::visit(const CachedObservableExpression & e)
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

    /*
     * ExpressionCacher
     */

    ExpressionCacher::ExpressionCacher(const ObservableCache & cache) :
        _cache(cache)
    {
    }

    Expression
    ExpressionCacher::visit(const BinaryExpression & e)
    {
        return BinaryExpression(e.op, e.lhs.accept_returning<Expression>(*this), e.rhs.accept_returning<Expression>(*this));
    }

    Expression
    ExpressionCacher::visit(const ConstantExpression & e)
    {
        return ConstantExpression(e);
    }

    Expression
    ExpressionCacher::visit(const ObservableNameExpression & e)
    {
        throw InternalError("Encountered ObservableNameExpression in ExpressionCacher::visit()");

        return e;
    }

    Expression
    ExpressionCacher::visit(const ObservableExpression & e)
    {
        auto id = _cache.add(e.observable);

        return CachedObservableExpression(_cache, id, e.kinematics_specification);
    }

    Expression
    ExpressionCacher::visit(const CachedObservableExpression & e)
    {
        throw InternalError("Encountered CachedObservableExpression in ExpressionCacher::visit");

        return e;
    }
}
