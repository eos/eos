/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2014 Frederik Beaujean
 * Copyright (c) 2022 Méril Reboud
 * Copyright (c) 2022 Danny van Dyk
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
#include <eos/utils/options-impl.hh>

namespace eos
{
    struct TestObservable : public Observable
    {
            Parameters p;

            Kinematics k;

            Options o;

            std::vector<KinematicVariable> kv;

            const QualifiedName & observable_name;

            const std::function<double(const Parameters &, const std::vector<KinematicVariable> &, const Options &)> & function;

            TestObservable(const Parameters & p, const Kinematics & k, const Options & o, const QualifiedName & observable_name,
                           const std::vector<std::string> &                                                                           kinematic_variable_names,
                           const std::function<double(const Parameters &, const std::vector<KinematicVariable> &, const Options &)> & function);

            virtual ~TestObservable();

            virtual double evaluate() const;

            virtual ObservablePtr clone() const;

            virtual ObservablePtr clone(const Parameters & parameters) const;

            virtual Parameters parameters();

            virtual Kinematics kinematics();

            virtual Options options();

            const QualifiedName & name() const;
    };

    class TestObservableEntry : public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            Unit _unit;

            std::function<double(const Parameters &, const std::vector<KinematicVariable> &, const Options &)> _function;

            std::vector<std::string> _kinematics_names;

            std::vector<OptionSpecification> _options;

        public:
            TestObservableEntry(const QualifiedName & name, const std::string & latex, const Unit & unit,
                                const std::function<double(const Parameters &, const std::vector<KinematicVariable> &, const Options &)> & function,
                                const std::vector<std::string> &                                                                           kinematics_names);

            virtual ~TestObservableEntry();

            virtual const QualifiedName & name() const;

            virtual const std::string & latex() const;

            virtual const Unit & unit() const;

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const;

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const;

            virtual ObservableEntry::OptionIterator begin_options() const;

            virtual ObservableEntry::OptionIterator end_options() const;

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const;

            virtual std::ostream & insert(std::ostream & os) const;
    };
} // namespace eos
