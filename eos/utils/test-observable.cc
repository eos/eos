
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2014 Frederik Beaujean
 * Copyright (c) 2022-2024 Danny van Dyk
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
#include <eos/utils/test-observable.hh>

namespace eos
{
    TestObservable::TestObservable(const Parameters & p, const Kinematics & k, const Options & o, const QualifiedName & observable_name,
                                   const std::vector<std::string> &                                                                           kinematic_variable_names,
                                   const std::function<double(const Parameters &, const std::vector<KinematicVariable> &, const Options &)> & function) :
        p(p),
        o(o),
        observable_name(observable_name),
        function(function)
    {
        for (auto kvn : kinematic_variable_names)
        {
            kv.emplace_back(std::move(k[kvn]));
            this->uses_kinematic(kv.back().id());
        }
    }

    TestObservable::~TestObservable() {}

    double
    TestObservable::evaluate() const
    {
        return function(p, kv, o);
    }

    ObservablePtr
    TestObservable::clone() const
    {
        throw InternalError("Not yet implemented");
    }

    ObservablePtr
    TestObservable::clone(const Parameters & /*parameters*/) const
    {
        throw InternalError("Not yet implemented");
    }

    Parameters
    TestObservable::parameters()
    {
        return p;
    }

    Kinematics
    TestObservable::kinematics()
    {
        return k;
    }

    Options
    TestObservable::options()
    {
        return o;
    }

    const QualifiedName &
    TestObservable::name() const
    {
        return observable_name;
    }

    TestObservableEntry::TestObservableEntry(const QualifiedName & name, const std::string & latex, const Unit & unit,
                                             const std::function<double(const Parameters &, const std::vector<KinematicVariable> &, const Options &)> & function,
                                             const std::vector<std::string> &                                                                           kinematics_names) :
        _name(name),
        _latex(latex),
        _unit(unit),
        _function(function),
        _kinematics_names(kinematics_names)
    {
    }

    TestObservableEntry::~TestObservableEntry() {}

    const QualifiedName &
    TestObservableEntry::name() const
    {
        return _name;
    }

    const std::string &
    TestObservableEntry::latex() const
    {
        return _latex;
    }

    const Unit &
    TestObservableEntry::unit() const
    {
        return _unit;
    }

    ObservableEntry::KinematicVariableIterator
    TestObservableEntry::begin_kinematic_variables() const
    {
        return &(*_kinematics_names.begin());
    }

    ObservableEntry::KinematicVariableIterator
    TestObservableEntry::end_kinematic_variables() const
    {
        auto temp = (&(*_kinematics_names.rbegin())) + 1;

        return ObservableEntry::KinematicVariableIterator(temp);
    }

    ObservableEntry::OptionIterator
    TestObservableEntry::begin_options() const
    {
        return _options.begin();
    }

    ObservableEntry::OptionIterator
    TestObservableEntry::end_options() const
    {
        return _options.end();
    }

    ObservablePtr
    TestObservableEntry::make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
    {
        return ObservablePtr(new TestObservable(parameters, kinematics, options, _name, _kinematics_names, _function));
    }

    std::ostream &
    TestObservableEntry::insert(std::ostream & os) const
    {
        os << "    type: test observable (name=" << _name << ")" << std::endl;

        return os;
    }
} // namespace eos
