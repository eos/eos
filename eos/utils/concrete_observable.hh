/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2015, 2016, 2017 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_CONCRETE_OBSERVABLE_HH
#define EOS_GUARD_EOS_UTILS_CONCRETE_OBSERVABLE_HH 1

#include <eos/observable-impl.hh>
#include <eos/utils/join.hh>
#include <eos/utils/log.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/units.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <array>
#include <functional>
#include <string>
#include <tuple>

namespace eos
{
    template <typename Decay_, typename ... Args_>
    class ConcreteObservable :
        public Observable
    {
        public:

        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            Decay_ _decay;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple;

        public:
            ConcreteObservable(const QualifiedName & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _decay(parameters, options),
                _function(function),
                _kinematics_names(kinematics_names),
                _argument_tuple(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names, &_decay))
            {
                uses(_decay);
                uses(Decay_::references);
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return std::apply(_function, values);
            };

            virtual Parameters parameters()
            {
                return _parameters;
            };

            virtual Kinematics kinematics()
            {
                return _kinematics;
            };

            virtual Options options()
            {
                return _options;
            }

            virtual ObservablePtr clone() const
            {
                return ObservablePtr(new ConcreteObservable(_name, _parameters.clone(), _kinematics.clone(), _options, _function, _kinematics_names));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteObservable(_name, parameters, _kinematics.clone(), _options, _function, _kinematics_names));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            Unit _unit;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array;

            Options _forced_options;

        public:
            ConcreteObservableEntry(const QualifiedName & name, const std::string & latex,
                    const Unit & unit,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names,
                    const Options & forced_options) :
                _name(name),
                _latex(latex),
                _unit(unit),
                _function(function),
                _kinematics_names(kinematics_names),
                _kinematics_names_array(impl::make_array<const std::string>(kinematics_names)),
                _forced_options(forced_options)
            {
            }

            ~ConcreteObservableEntry()
            {
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual const std::string & latex() const
            {
                return _latex;
            }

            virtual const Unit & unit() const
            {
                return _unit;
            }

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const
            {
                return _kinematics_names_array.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const
            {
                return _kinematics_names_array.end();
            }

            virtual ObservableEntry::OptionIterator begin_options() const
            {
                return Decay_::begin_options();
            }

            virtual ObservableEntry::OptionIterator end_options() const
            {
                return Decay_::end_options();
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                for (const auto & [key, forced_value] : _forced_options)
                {
                    if (options.has(key) && (options[key] != forced_value))
                    {
                        Log::instance()->message("[ConcreteObservableEntry.make]", ll_error)
                            << "Observable '" << _name << "' forces option key '" << key << "' to value '" << forced_value << "', overriding user-provided value '" << options[key] << "'";
                    }
                }
                return ObservablePtr(new ConcreteObservable<Decay_, Args_ ...>(_name, parameters, kinematics, options + _forced_options, _function, _kinematics_names));
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableEntryPtr make_concrete_observable_entry(const QualifiedName & name, const std::string & latex,
            const Unit & unit,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names,
            const Options & forced_options)
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return std::make_shared<ConcreteObservableEntry<Decay_, Args_ ...>>(name, latex,
                unit,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(function)),
                kinematics_names, forced_options);
    }
}


#endif
