/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2025 Danny van Dyk
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

#include "python/_eos/external-observable.hh"

#include <eos/observable-impl.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

using boost::python::extract;
using boost::python::object;
using boost::python::stl_input_iterator;

namespace eos
{
    ExternalObservable::ExternalObservable(const QualifiedName & name, object provider, const Parameters & parameters, const Kinematics & kinematics, const Options & options) :
        _name(name),
        _provider(provider),
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options)
    {
        if (! PyCallable_Check(_provider.ptr()))
        {
            throw InternalError("ExternalObservable encountered an observable provider that is not callable/constructible");
        }

        auto o = boost::python::object(_provider(parameters, kinematics, options));

        _evaluate = o.attr("evaluate");

        if (_evaluate.is_none())
        {
            throw InternalError("ExternalObservable encountered an observable provider that lacks the 'evaluate' attribute");
        }

        if (! PyCallable_Check(_evaluate.ptr()))
        {
            throw InternalError("ExternalObservable encountered an 'evaluate' attribute that is not callable");
        }
    }

    ExternalObservable::~ExternalObservable() = default;

    const QualifiedName &
    ExternalObservable::name() const
    {
        return _name;
    }

    double
    ExternalObservable::evaluate() const
    {
        return extract<double>(_evaluate());
    }

    Parameters
    ExternalObservable::parameters()
    {
        return _parameters;
    }

    Kinematics
    ExternalObservable::kinematics()
    {
        return _kinematics;
    }

    Options
    ExternalObservable::options()
    {
        return _options;
    }

    ObservablePtr
    ExternalObservable::clone() const
    {
        return ObservablePtr(new ExternalObservable(_name, _provider, _parameters, _kinematics, _options));
    }

    ObservablePtr
    ExternalObservable::clone(const Parameters & parameters) const
    {
        return ObservablePtr(new ExternalObservable(_name, _provider, parameters, _kinematics, _options));
    }

    ExternalObservableEntry::ExternalObservableEntry(const QualifiedName & name, object provider, const std::string & latex, const Unit & unit) :
        _name(name),
        _provider(provider),
        _latex(latex),
        _unit(unit)
    {
        object kinematic_variables = _provider.attr("kinematic_variables");
        if (kinematic_variables.is_none())
        {
            throw InternalError("ExternalObservableEntry encountered a factory that posesses no 'kinematic_variables' attribute");
        }

        _kinematic_variables = std::vector<std::string>(stl_input_iterator<std::string>(kinematic_variables), stl_input_iterator<std::string>());
    }

    ExternalObservableEntry::~ExternalObservableEntry() = default;

    ObservablePtr
    ExternalObservableEntry::make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
    {
        return ObservablePtr(new ExternalObservable(_name, _provider, parameters, kinematics, options));
    }

    const QualifiedName &
    ExternalObservableEntry::name() const
    {
        return _name;
    }

    const std::string &
    ExternalObservableEntry::latex() const
    {
        return _latex;
    }

    const Unit &
    ExternalObservableEntry::unit() const
    {
        return _unit;
    }

    ObservableEntry::KinematicVariableIterator
    ExternalObservableEntry::begin_kinematic_variables() const
    {
        return ObservableEntry::KinematicVariableIterator(&(*_kinematic_variables.cbegin()));
    }

    ObservableEntry::KinematicVariableIterator
    ExternalObservableEntry::end_kinematic_variables() const
    {
        return ObservableEntry::KinematicVariableIterator(&(*_kinematic_variables.cend()));
    }

    ObservableEntry::OptionIterator
    ExternalObservableEntry::begin_options() const
    {
        return ObservableEntry::OptionIterator(_used_options.cbegin());
    }

    ObservableEntry::OptionIterator
    ExternalObservableEntry::end_options() const
    {
        return ObservableEntry::OptionIterator(_used_options.cend());
    }

    std::shared_ptr<const ObservableEntry>
    register_python_observable(const QualifiedName & name, boost::python::object provider, const std::string & latex, const Unit & unit)
    {
        std::shared_ptr<const ObservableEntry> entry(new ExternalObservableEntry(name, provider, latex, unit));

        ObservableEntries::instance()->insert_or_assign(name, entry);

        return entry;
    }
} // namespace eos
