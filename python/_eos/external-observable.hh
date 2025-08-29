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

#include "eos/observable.hh"

#include <boost/python.hpp>

namespace eos
{
    class ExternalObservable : public Observable
    {
        private:
            eos::QualifiedName    _name;
            boost::python::object _provider;
            Parameters            _parameters;
            Kinematics            _kinematics;
            Options               _options;
            boost::python::object _evaluate;

        public:
            ExternalObservable(const QualifiedName & name, boost::python::object provider, const Parameters & parameters, const Kinematics & kinematics, const Options & options);

            ~ExternalObservable() override;

            const QualifiedName & name() const override;

            double evaluate() const override;

            Kinematics kinematics() override;

            Parameters parameters() override;

            Options options() override;

            ObservablePtr clone() const override;

            ObservablePtr clone(const Parameters & parameters) const override;
    };

    class ExternalObservableEntry : public ObservableEntry
    {
        private:
            eos::QualifiedName               _name;
            boost::python::object            _provider;
            std::string                      _latex;
            Unit                             _unit;
            std::vector<std::string>         _kinematic_variables;
            std::vector<OptionSpecification> _used_options;

        public:
            ExternalObservableEntry(const QualifiedName & name, boost::python::object provider, const std::string & latex, const Unit & unit = Unit::Undefined());

            ~ExternalObservableEntry() override;

            ObservablePtr make(const Parameters &, const Kinematics &, const Options &) const override;

            const QualifiedName & name() const override;

            const std::string & latex() const override;

            const Unit & unit() const override;

            KinematicVariableIterator begin_kinematic_variables() const override;
            KinematicVariableIterator end_kinematic_variables() const override;

            OptionIterator begin_options() const override;
            OptionIterator end_options() const override;
    };

    std::shared_ptr<const ObservableEntry> register_python_observable(const QualifiedName & name, boost::python::object provider, const std::string & latex = "",
                                                                      const Unit & unit = Unit::Undefined());
} // namespace eos
