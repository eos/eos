/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2024 Danny van Dyk
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
#include <eos/utils/expression-cacher.hh>
#include <eos/utils/expression-cloner.hh>
#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-kinematic-reader.hh>
#include <eos/utils/expression-maker.hh>
#include <eos/utils/expression-observable.hh>
#include <eos/utils/expression-parser.hh>
#include <eos/utils/expression-used-kinematics-reader.hh>
#include <eos/utils/expression-used-parameter-reader.hh>
#include <eos/utils/log.hh>

#include <set>

namespace eos
{
    using eos::exp::Expression;
    using eos::exp::ExpressionPtr;

    ExpressionObservable::ExpressionObservable(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                               const ExpressionPtr & expression) :
        _name(name),
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options)
    {
        if (expression == nullptr)
        {
            throw InternalError("Empty expression encountered in ExpressionObservable!");
        }

        exp::ExpressionMaker maker(parameters, kinematics, options, this);
        _expression = eos::exp::ExpressionPtr(new eos::exp::Expression(std::move(std::visit(maker, *expression))));

        exp::ExpressionUsedParameterReader reader;
        std::visit(reader, *_expression);

        for (Parameter::Id id : reader.parameter_ids)
        {
            this->uses(id);
        }

        exp::ExpressionUsedKinematicsReader used_kinematic_reader;
        std::visit(used_kinematic_reader, *_expression);

        for (KinematicVariable::Id id : used_kinematic_reader.kinematic_variable_ids)
        {
            this->uses_kinematic(id);
        }
    }

    ExpressionObservable::ExpressionObservable(const QualifiedName & name, const ObservableCache & cache, const Kinematics & kinematics, const Options & options,
                                               const ExpressionPtr & expression) :
        _name(name),
        _parameters(cache.parameters()),
        _kinematics(kinematics),
        _options(options)
    {
        if (expression == nullptr)
        {
            throw InternalError("Empty expression encountered in ExpressionObservable!");
        }

        exp::ExpressionCacher cacher(cache);
        _expression = eos::exp::ExpressionPtr(new eos::exp::Expression(std::move(std::visit(cacher, *expression))));

        exp::ExpressionUsedParameterReader reader;
        std::visit(reader, *_expression);

        for (Parameter::Id id : reader.parameter_ids)
        {
            this->uses(id);
        }

        exp::ExpressionUsedKinematicsReader used_kinematic_reader;
        std::visit(used_kinematic_reader, *_expression);

        for (KinematicVariable::Id id : used_kinematic_reader.kinematic_variable_ids)
        {
            this->uses_kinematic(id);
        }
    }

    double
    ExpressionObservable::evaluate() const
    {
        exp::ExpressionEvaluator evaluator;

        return std::visit(evaluator, *_expression);
    }

    ObservablePtr
    ExpressionObservable::clone() const
    {
        auto parameters = _parameters.clone();
        auto kinematics = _kinematics.clone();

        exp::ExpressionCloner cloner(parameters, kinematics, _options);

        return ObservablePtr(new ExpressionObservable(_name, parameters, kinematics, _options, ExpressionPtr(new Expression(std::move(std::visit(cloner, *_expression))))));
    }

    ObservablePtr
    ExpressionObservable::clone(const Parameters & parameters) const
    {
        auto kinematics = _kinematics.clone();

        exp::ExpressionCloner cloner(parameters, kinematics, _options);

        return ObservablePtr(new ExpressionObservable(_name, parameters, kinematics, _options, ExpressionPtr(new Expression(std::move(std::visit(cloner, *_expression))))));
    }

    ExpressionObservableEntry::ExpressionObservableEntry(const QualifiedName & name, const std::string & latex, const Unit & unit, const ExpressionPtr & expression,
                                                         const Options & forced_options) :
        _name(name),
        _latex(latex),
        _unit(unit),
        _expression(expression),
        _forced_options(forced_options)
    {
        if (_expression == nullptr)
        {
            throw InternalError("Empty expression encountered in ExpressionObservableEntry!");
        }

        // Read kinematic variables
        exp::ExpressionKinematicReader kinematic_reader;
        std::visit(kinematic_reader, *_expression);

        // Check the absence of overlap between the used kinematic variables and the aliased variables
        std::set<std::string> intersection;
        std::set_intersection(kinematic_reader.kinematics.begin(),
                              kinematic_reader.kinematics.end(),
                              kinematic_reader.aliases.begin(),
                              kinematic_reader.aliases.end(),
                              std::inserter(intersection, intersection.begin()));

        if (! intersection.empty())
        {
            throw InternalError("An aliased kinematic variable is still present in the expression.");
        }

        _kinematics_names.assign(kinematic_reader.kinematics.begin(), kinematic_reader.kinematics.end());
    }

    ObservableEntry::KinematicVariableIterator
    ExpressionObservableEntry::begin_kinematic_variables() const
    {
        return _kinematics_names.data();
    }

    ObservableEntry::KinematicVariableIterator
    ExpressionObservableEntry::end_kinematic_variables() const
    {
        return _kinematics_names.data() + _kinematics_names.size();
    }

    ObservableEntry::OptionIterator
    ExpressionObservableEntry::begin_options() const
    {
        return _option_specifications.cbegin();
    }

    ObservableEntry::OptionIterator
    ExpressionObservableEntry::end_options() const
    {
        return _option_specifications.cend();
    }

    ObservablePtr
    ExpressionObservableEntry::make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
    {
        if (_expression == nullptr)
        {
            throw InternalError("Empty expression encountered in ExpressionObservableEntry::make!");
        }

        for (const auto & fo : _forced_options)
        {
            const auto & key = std::get<0>(fo);
            if (options.has(key))
            {
                Log::instance()->message("[ExpressionObservableEntry.make]", ll_warning) << "Observable '" << _name << "' forces option key '" << key << "' to value '"
                                                                                         << _forced_options[key] << "', overriding user-provided value '" << options[key] << "'";
            }
        }

        return ObservablePtr(new ExpressionObservable(_name, parameters, kinematics, options + _forced_options, _expression));
    }
} // namespace eos
