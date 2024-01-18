/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2023-2024 Danny van Dyk
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
#include <eos/observable-impl.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/expression-cacher.hh>
#include <eos/utils/expression-cloner.hh>
#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-kinematic-reader.hh>
#include <eos/utils/expression-maker.hh>
#include <eos/utils/expression-printer.hh>
#include <eos/utils/expression-used-parameter-reader.hh>
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

    void ExpressionPrinter::visit(ParameterNameExpression & e)
    {
        _os << "ParameterNameExpression(" << e.parameter_name.full() << ")";
    }

    void ExpressionPrinter::visit(ParameterExpression & e)
    {
        _os << "ParameterExpression(" << e.parameter.name() << ")";
    }

    void ExpressionPrinter::visit(KinematicVariableNameExpression & e)
    {
        _os << "KinematicVariableNameExpression(" << e.variable_name << ")";
    }

    void ExpressionPrinter::visit(KinematicVariableExpression & e)
    {
        _os << "KinematicVariableExpression(" << e.kinematic_variable.name() << ")";
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
    ExpressionEvaluator::visit(ParameterNameExpression &)
    {
        throw InternalError("Encountered ParameterNameExpression in ExpressionEvaluator::visit");

        return 0.0;
    }

    double
    ExpressionEvaluator::visit(ParameterExpression & e)
    {
        return e.parameter.evaluate();
    }

    double
    ExpressionEvaluator::visit(KinematicVariableNameExpression &)
    {
        throw InternalError("Encountered KinematicVariableNameExpression in ExpressionEvaluator::visit");

        return 0.0;
    }

    double
    ExpressionEvaluator::visit(KinematicVariableExpression & e)
    {
        return e.kinematic_variable.evaluate();
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
    ExpressionCloner::visit(const ParameterNameExpression & e)
    {
        return ParameterNameExpression(e);
    }

    Expression
    ExpressionCloner::visit(const ParameterExpression & e)
    {
        return ParameterExpression(this->_parameters[e.parameter.name()]);
    }

    Expression
    ExpressionCloner::visit(const KinematicVariableNameExpression & e)
    {
        return KinematicVariableNameExpression(e);
    }

    Expression
    ExpressionCloner::visit(const KinematicVariableExpression & e)
    {
        auto kinematic_variable = this->_kinematics[e.kinematic_variable.name()];

        return KinematicVariableExpression(kinematic_variable);
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

        // Clear aliases installed as part of this expression
        for (const auto & alias: kinematics_aliases)
        {
            _kinematics.remove_alias(alias.first);
        }

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

        // Clear aliases installed as part of this expression
        for (const auto & alias: kinematics_aliases)
        {
            _kinematics.remove_alias(alias.first);
        }

        // Record the used parameters
        if (_parameter_user)
        {
            _parameter_user->uses(*observable);
        }

        return ObservableExpression(observable, e.kinematics_specification);
    }

    Expression
    ExpressionMaker::visit(const ParameterNameExpression & e)
    {
        Parameter parameter = this->_parameters[e.parameter_name];

        // Record the used parameters
        if (_parameter_user)
        {
            _parameter_user->uses(parameter.id());
        }

        return ParameterExpression(parameter);
    }

    Expression
    ExpressionMaker::visit(const ParameterExpression & e)
    {
        // Rebuild the parameter expression with the local set of parameters
        Parameter parameter = this->_parameters[e.parameter.name()];

        // Record the used parameters
        if (_parameter_user)
        {
            _parameter_user->uses(parameter.id());
        }

        return ParameterExpression(parameter);
    }

    Expression
    ExpressionMaker::visit(const KinematicVariableNameExpression & e)
    {
        KinematicVariable kinematic_variable = this->_kinematics[e.variable_name];

        return KinematicVariableExpression(kinematic_variable);
    }

    Expression
    ExpressionMaker::visit(const KinematicVariableExpression & e)
    {
        return KinematicVariableExpression(e.kinematic_variable);
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
    void
    ExpressionKinematicReader::clear()
    {
        this->kinematics.clear();
        this->aliases.clear();
    }

    void
    ExpressionKinematicReader::visit(const BinaryExpression & e)
    {
        e.lhs.accept(*this);
        e.rhs.accept(*this);
    }

    void
    ExpressionKinematicReader::visit(const ConstantExpression &)
    {
    }

    void
    ExpressionKinematicReader::visit(const ObservableNameExpression & e)
    {
        std::set<std::string> kinematic_set;
        std::set<std::string> alias_set;

        const auto & entries = impl::observable_entries;
        const auto & it = entries.find(e.observable_name);

        // check if 'e' matches the name of a known observable
        if (it != entries.end())
        {
            const auto entry = it->second;

            kinematic_set.insert(entry->begin_kinematic_variables(), entry->end_kinematic_variables());

            const auto & kinematics_values = e.kinematics_specification.values;
            const auto & kinematics_aliases = e.kinematics_specification.aliases;

            // Set or alias kinematic specifications
            for (const auto & value : kinematics_values)
            {
                kinematic_set.erase(value.first);
                alias_set.insert(value.first);
            }
            for (const auto & alias : kinematics_aliases)
            {
                kinematic_set.erase(alias.first);
                alias_set.insert(alias.first);
                kinematic_set.insert(alias.second);
            }

            this->kinematics.insert(kinematic_set.begin(), kinematic_set.end());
            this->aliases.insert(alias_set.begin(), alias_set.end());

            return;
        }

        // otherwise check if 'e' matches the name of a parameter
        if (kinematic_set.empty())
        {
            auto parameters = eos::Parameters::Defaults();

            // observable_name.full() is used below because if the observable has options, it cannot be a parameter
            auto i = std::find_if(parameters.begin(), parameters.end(), [&] (const Parameter & p) { return p.name() == e.observable_name.full(); });
            if (parameters.end() != i)
            {
                return;
            }
        }

        throw UnknownObservableError("Expression '" + e.observable_name.full() + "' is neither a known Observable nor a Parameter");
    }

    void
    ExpressionKinematicReader::visit(const ObservableExpression & e)
    {
        std::set<std::string> kinematic_set;
        std::set<std::string> alias_set;

        const auto & kinematics_values = e.kinematics_specification.values;
        const auto & kinematics_aliases = e.kinematics_specification.aliases;

        // Set or alias kinematic specifications
        for (const auto & value : kinematics_values)
        {
            kinematic_set.erase(value.first);
            alias_set.insert(value.first);
        }
        for (const auto & alias : kinematics_aliases)
        {
            kinematic_set.erase(alias.first);
            alias_set.insert(alias.first);
            kinematic_set.insert(alias.second);
        }

        this->kinematics.insert(kinematic_set.begin(), kinematic_set.end());
        this->aliases.insert(alias_set.begin(), alias_set.end());
    }

    void
    ExpressionKinematicReader::visit(const ParameterNameExpression &)
    {
    }

    void
    ExpressionKinematicReader::visit(const ParameterExpression &)
    {
    }

    void
    ExpressionKinematicReader::visit(const KinematicVariableNameExpression & e)
    {
        this->kinematics.insert(e.variable_name);
    }

    void
    ExpressionKinematicReader::visit(const KinematicVariableExpression & e)
    {
        this->kinematics.insert(e.kinematic_variable.name());
    }

    void
    ExpressionKinematicReader::visit(const CachedObservableExpression & e)
    {
        std::set<std::string> kinematic_set;
        std::set<std::string> alias_set;

        const auto & kinematics_values = e.kinematics_specification.values;
        const auto & kinematics_aliases = e.kinematics_specification.aliases;

        // Set or alias kinematic specifications
        for (const auto & value : kinematics_values)
        {
            kinematic_set.erase(value.first);
            alias_set.insert(value.first);
        }
        for (const auto & alias : kinematics_aliases)
        {
            kinematic_set.erase(alias.first);
            alias_set.insert(alias.first);
            kinematic_set.insert(alias.second);
        }

        this->kinematics.insert(kinematic_set.begin(), kinematic_set.end());
        this->aliases.insert(alias_set.begin(), alias_set.end());
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
    ExpressionCacher::visit(const ParameterNameExpression & e)
    {
        throw InternalError("Encountered ParameterNameExpression in ExpressionCacher::visit()");

        return e;
    }

    Expression
    ExpressionCacher::visit(const ParameterExpression & e)
    {
        return ParameterExpression(e);
    }

    Expression
    ExpressionCacher::visit(const KinematicVariableNameExpression & e)
    {
        throw InternalError("Encountered KinematicVariableNameExpression in ExpressionCacher::visit()");

        return e;
    }

    Expression
    ExpressionCacher::visit(const KinematicVariableExpression & e)
    {
        return KinematicVariableExpression(e);
    }

    Expression
    ExpressionCacher::visit(const CachedObservableExpression & e)
    {
        throw InternalError("Encountered CachedObservableExpression in ExpressionCacher::visit");

        return e;
    }

    /*
     * ExpressionUsedParameterReader
     */
    void
    ExpressionUsedParameterReader::visit(const BinaryExpression & e)
    {
        e.lhs.accept(*this);
        e.rhs.accept(*this);
    }

    void
    ExpressionUsedParameterReader::visit(const ConstantExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::visit(const ObservableNameExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::visit(const ObservableExpression & e)
    {
        const ParameterUser & parameter_user = static_cast<const ParameterUser &>(*e.observable);
        for (const auto & parameter_id : parameter_user)
        {
            this->parameter_ids.insert(parameter_id);
        }
    }

    void
    ExpressionUsedParameterReader::visit(const ParameterNameExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::visit(const ParameterExpression & e)
    {
        this->parameter_ids.insert(e.parameter.id());
    }

    void
    ExpressionUsedParameterReader::visit(const KinematicVariableNameExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::visit(const KinematicVariableExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::visit(const CachedObservableExpression & e)
    {
        const ParameterUser & parameter_user = static_cast<const ParameterUser &>(*e.cache.observable(e.id));
        for (const auto & parameter_id : parameter_user)
        {
            this->parameter_ids.insert(parameter_id);
        }
    }
}
