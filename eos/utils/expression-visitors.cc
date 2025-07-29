/*
 * Copyright (c) 2021      Méril Reboud
 * Copyright (c) 2023-2025 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/observable.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/expression-cacher.hh>
#include <eos/utils/expression-cloner.hh>
#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-kinematic-reader.hh>
#include <eos/utils/expression-maker.hh>
#include <eos/utils/expression-printer.hh>
#include <eos/utils/expression-used-parameter-reader.hh>
#include <eos/utils/expression.hh>
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

    void
    ExpressionPrinter::operator() (BinaryExpression & e)
    {
        _os << "BinaryExpression(";
        std::visit(*this, *e.lhs);
        _os << " " << e.op << " ";
        std::visit(*this, *e.rhs);
        _os << ")";
    }

    void
    ExpressionPrinter::operator() (FunctionExpression & e)
    {
        _os << "FunctionExpression(";
        _os << e.fname << ", ";
        std::visit(*this, *e.arg);
        _os << ")";
    }

    void
    ExpressionPrinter::operator() (ConstantExpression & e)
    {
        _os << "ConstantExpression(" << e.value << ")";
    }

    void
    ExpressionPrinter::operator() (ObservableNameExpression & e)
    {
        _os << "ObservableNameExpression(" << e.observable_name.full();
        if (! e.kinematics_specification.aliases.empty())
        {
            auto a = e.kinematics_specification.aliases.cbegin();
            _os << ", aliases=[" << a->first << "=>" << a->second;
            ++a;
            for (auto a_end = e.kinematics_specification.aliases.cend(); a != a_end; ++a)
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
            for (auto v_end = e.kinematics_specification.values.cend(); v != v_end; ++v)
            {
                _os << "," << v->first << "=" << v->second;
            }
            _os << "]";
        }
        _os << ")";
    }

    void
    ExpressionPrinter::operator() (ObservableExpression & e)
    {
        _os << "ObservableExpression(" << e.observable->name() << ")";
        if (! e.kinematics_specification.aliases.empty())
        {
            auto a = e.kinematics_specification.aliases.cbegin();
            _os << ", aliases=[" << a->first << "=>" << a->second;
            ++a;
            for (auto a_end = e.kinematics_specification.aliases.cend(); a != a_end; ++a)
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
            for (auto v_end = e.kinematics_specification.values.cend(); v != v_end; ++v)
            {
                _os << "," << v->first << "=" << v->second;
            }
            _os << "]";
        }
        _os << ")";
    }

    void
    ExpressionPrinter::operator() (ParameterNameExpression & e)
    {
        _os << "ParameterNameExpression(" << e.parameter_name.full() << ")";
    }

    void
    ExpressionPrinter::operator() (ParameterExpression & e)
    {
        _os << "ParameterExpression(" << e.parameter.name() << ")";
    }

    void
    ExpressionPrinter::operator() (KinematicVariableNameExpression & e)
    {
        _os << "KinematicVariableNameExpression(" << e.variable_name << ")";
    }

    void
    ExpressionPrinter::operator() (KinematicVariableExpression & e)
    {
        _os << "KinematicVariableExpression(" << e.kinematic_variable.name() << ")";
    }

    void
    ExpressionPrinter::operator() (CachedObservableExpression & e)
    {
        _os << "CachedObservableExpression(id=" << e.id;
        _os << ", name='" << e.cache.observable(e.id)->name().full() << "'";
        _os << ")";
    }

    /*
     * ExpressionEvaluator
     */

    double
    ExpressionEvaluator::operator() (const BinaryExpression & e)
    {
        BinaryExpression::func f = BinaryExpression::Method(e.op);

        return f(std::visit(*this, *e.lhs), std::visit(*this, *e.rhs));
    }

    double
    ExpressionEvaluator::operator() (const FunctionExpression & e)
    {
        return e.f(std::visit(*this, *e.arg));
    }

    double
    ExpressionEvaluator::operator() (const ConstantExpression & e)
    {
        return e.value;
    }

    double
    ExpressionEvaluator::operator() (const ObservableNameExpression &)
    {
        throw InternalError("Encountered ObserableNameExpression in ExpressionEvaluator::operator() ");

        return 0.0;
    }

    double
    ExpressionEvaluator::operator() (const ObservableExpression & e)
    {
        return e.observable->evaluate();
    }

    double
    ExpressionEvaluator::operator() (const ParameterNameExpression &)
    {
        throw InternalError("Encountered ParameterNameExpression in ExpressionEvaluator::operator() ");

        return 0.0;
    }

    double
    ExpressionEvaluator::operator() (const ParameterExpression & e)
    {
        return e.parameter.evaluate();
    }

    double
    ExpressionEvaluator::operator() (const KinematicVariableNameExpression &)
    {
        throw InternalError("Encountered KinematicVariableNameExpression in ExpressionEvaluator::operator() ");

        return 0.0;
    }

    double
    ExpressionEvaluator::operator() (const KinematicVariableExpression & e)
    {
        return e.kinematic_variable.evaluate();
    }

    double
    ExpressionEvaluator::operator() (const CachedObservableExpression & e)
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
    ExpressionCloner::operator() (const BinaryExpression & e)
    {
        return BinaryExpression(e.op, ExpressionPtr(new Expression(std::move(std::visit(*this, *e.lhs)))), ExpressionPtr(new Expression(std::move(std::visit(*this, *e.rhs)))));
    }

    Expression
    ExpressionCloner::operator() (const FunctionExpression & e)
    {
        return FunctionExpression(e.fname, ExpressionPtr(new Expression(std::move(std::visit(*this, *e.arg)))));
    }

    Expression
    ExpressionCloner::operator() (const ConstantExpression & e)
    {
        return ConstantExpression(e);
    }

    Expression
    ExpressionCloner::operator() (const ObservableNameExpression & e)
    {
        return ObservableNameExpression(e);
    }

    Expression
    ExpressionCloner::operator() (const ObservableExpression & e)
    {
        const auto & kinematics_values  = e.kinematics_specification.values;
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
        auto observable = Observable::make(e.observable->name(), this->_parameters, this->_kinematics, this->_options + e.observable->options());
        // Clear alias map
        _kinematics.clear_aliases();

        return ObservableExpression(observable, e.kinematics_specification);
    }

    Expression
    ExpressionCloner::operator() (const ParameterNameExpression & e)
    {
        return ParameterNameExpression(e);
    }

    Expression
    ExpressionCloner::operator() (const ParameterExpression & e)
    {
        return ParameterExpression(this->_parameters[e.parameter.name()]);
    }

    Expression
    ExpressionCloner::operator() (const KinematicVariableNameExpression & e)
    {
        return KinematicVariableNameExpression(e);
    }

    Expression
    ExpressionCloner::operator() (const KinematicVariableExpression & e)
    {
        auto kinematic_variable = this->_kinematics[e.kinematic_variable.name()];

        return KinematicVariableExpression(kinematic_variable);
    }

    Expression
    ExpressionCloner::operator() (const CachedObservableExpression & e)
    {
        const auto & kinematics_values  = e.kinematics_specification.values;
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
        auto observable = Observable::make(e.cache.observable(e.id)->name(), this->_parameters, this->_kinematics, this->_options + e.cache.observable(e.id)->options());
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
    ExpressionMaker::operator() (const BinaryExpression & e)
    {
        return BinaryExpression(e.op, ExpressionPtr(new Expression(std::move(std::visit(*this, *e.lhs)))), ExpressionPtr(new Expression(std::move(std::visit(*this, *e.rhs)))));
    }

    Expression
    ExpressionMaker::operator() (const FunctionExpression & e)
    {
        return FunctionExpression(e.fname, ExpressionPtr(new Expression(std::move(std::visit(*this, *e.arg)))));
    }

    Expression
    ExpressionMaker::operator() (const ConstantExpression & e)
    {
        return ConstantExpression(e);
    }

    Expression
    ExpressionMaker::operator() (const ObservableNameExpression & e)
    {
        const auto & kinematics_values  = e.kinematics_specification.values;
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
        auto observable = Observable::make(e.observable_name, this->_parameters, this->_kinematics, this->_options);

        // Clear aliases installed as part of this expression
        for (const auto & alias : kinematics_aliases)
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
    ExpressionMaker::operator() (const ObservableExpression & e)
    {
        // Rebuild the observable expression with the local set of parameters

        const auto & kinematics_values  = e.kinematics_specification.values;
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
        auto observable = Observable::make(e.observable->name(), this->_parameters, this->_kinematics, this->_options + e.observable->options());

        // Clear aliases installed as part of this expression
        for (const auto & alias : kinematics_aliases)
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
    ExpressionMaker::operator() (const ParameterNameExpression & e)
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
    ExpressionMaker::operator() (const ParameterExpression & e)
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
    ExpressionMaker::operator() (const KinematicVariableNameExpression & e)
    {
        KinematicVariable kinematic_variable = this->_kinematics[e.variable_name];

        return KinematicVariableExpression(kinematic_variable);
    }

    Expression
    ExpressionMaker::operator() (const KinematicVariableExpression & e)
    {
        return KinematicVariableExpression(e.kinematic_variable);
    }

    Expression
    ExpressionMaker::operator() (const CachedObservableExpression & e)
    {
        throw InternalError("Encountered CachedObservableExpression in ExpressionMaker::operator() ");

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
    ExpressionKinematicReader::operator() (const BinaryExpression & e)
    {
        std::visit(*this, *e.lhs);
        std::visit(*this, *e.rhs);
    }

    void
    ExpressionKinematicReader::operator() (const FunctionExpression & e)
    {
        std::visit(*this, *e.arg);
    }

    void
    ExpressionKinematicReader::operator() (const ConstantExpression &)
    {
    }

    void
    ExpressionKinematicReader::operator() (const ObservableNameExpression & e)
    {
        std::set<std::string> kinematic_set;
        std::set<std::string> alias_set;

        const auto & entries = impl::observable_entries;
        const auto & it      = entries.find(e.observable_name);

        // check if 'e' matches the name of a known observable
        if (it != entries.end())
        {
            const auto entry = it->second;

            kinematic_set.insert(entry->begin_kinematic_variables(), entry->end_kinematic_variables());

            const auto & kinematics_values  = e.kinematics_specification.values;
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
            auto i = std::find_if(parameters.begin(), parameters.end(), [&](const Parameter & p) { return p.name() == e.observable_name.full(); });
            if (parameters.end() != i)
            {
                return;
            }
        }

        throw UnknownObservableError("Expression '" + e.observable_name.full() + "' is neither a known Observable nor a Parameter");
    }

    void
    ExpressionKinematicReader::operator() (const ObservableExpression & e)
    {
        std::set<std::string> kinematic_set;
        std::set<std::string> alias_set;

        const auto & kinematics_values  = e.kinematics_specification.values;
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
    ExpressionKinematicReader::operator() (const ParameterNameExpression &)
    {
    }

    void
    ExpressionKinematicReader::operator() (const ParameterExpression &)
    {
    }

    void
    ExpressionKinematicReader::operator() (const KinematicVariableNameExpression & e)
    {
        this->kinematics.insert(e.variable_name);
    }

    void
    ExpressionKinematicReader::operator() (const KinematicVariableExpression & e)
    {
        this->kinematics.insert(e.kinematic_variable.name());
    }

    void
    ExpressionKinematicReader::operator() (const CachedObservableExpression & e)
    {
        std::set<std::string> kinematic_set;
        std::set<std::string> alias_set;

        const auto & kinematics_values  = e.kinematics_specification.values;
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
    ExpressionCacher::operator() (const BinaryExpression & e)
    {
        return BinaryExpression(e.op, ExpressionPtr(new Expression(std::move(std::visit(*this, *e.lhs)))), ExpressionPtr(new Expression(std::move(std::visit(*this, *e.rhs)))));
    }

    Expression
    ExpressionCacher::operator() (const FunctionExpression & e)
    {
        return FunctionExpression(e.fname, ExpressionPtr(new Expression(std::move(std::visit(*this, *e.arg)))));
    }

    Expression
    ExpressionCacher::operator() (const ConstantExpression & e)
    {
        return ConstantExpression(e);
    }

    Expression
    ExpressionCacher::operator() (const ObservableNameExpression & e)
    {
        throw InternalError("Encountered ObservableNameExpression in ExpressionCacher::operator() ()");

        return e;
    }

    Expression
    ExpressionCacher::operator() (const ObservableExpression & e)
    {
        auto id = _cache.add(e.observable);

        return CachedObservableExpression(_cache, id, e.kinematics_specification);
    }

    Expression
    ExpressionCacher::operator() (const ParameterNameExpression & e)
    {
        throw InternalError("Encountered ParameterNameExpression in ExpressionCacher::operator() ()");

        return e;
    }

    Expression
    ExpressionCacher::operator() (const ParameterExpression & e)
    {
        return ParameterExpression(e);
    }

    Expression
    ExpressionCacher::operator() (const KinematicVariableNameExpression & e)
    {
        throw InternalError("Encountered KinematicVariableNameExpression in ExpressionCacher::operator() ()");

        return e;
    }

    Expression
    ExpressionCacher::operator() (const KinematicVariableExpression & e)
    {
        return KinematicVariableExpression(e);
    }

    Expression
    ExpressionCacher::operator() (const CachedObservableExpression & e)
    {
        throw InternalError("Encountered CachedObservableExpression in ExpressionCacher::operator() ");

        return e;
    }

    /*
     * ExpressionUsedParameterReader
     */
    void
    ExpressionUsedParameterReader::operator() (const BinaryExpression & e)
    {
        std::visit(*this, *e.lhs);
        std::visit(*this, *e.rhs);
    }

    void
    ExpressionUsedParameterReader::operator() (const FunctionExpression & e)
    {
        std::visit(*this, *e.arg);
    }

    void
    ExpressionUsedParameterReader::operator() (const ConstantExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::operator() (const ObservableNameExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::operator() (const ObservableExpression & e)
    {
        const ParameterUser & parameter_user = static_cast<const ParameterUser &>(*e.observable);
        for (const auto & parameter_id : parameter_user)
        {
            this->parameter_ids.insert(parameter_id);
        }
    }

    void
    ExpressionUsedParameterReader::operator() (const ParameterNameExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::operator() (const ParameterExpression & e)
    {
        this->parameter_ids.insert(e.parameter.id());
    }

    void
    ExpressionUsedParameterReader::operator() (const KinematicVariableNameExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::operator() (const KinematicVariableExpression &)
    {
    }

    void
    ExpressionUsedParameterReader::operator() (const CachedObservableExpression & e)
    {
        const ParameterUser & parameter_user = static_cast<const ParameterUser &>(*e.cache.observable(e.id));
        for (const auto & parameter_id : parameter_user)
        {
            this->parameter_ids.insert(parameter_id);
        }
    }
} // namespace eos::exp
