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

#ifndef EOS_GUARD_SRC_UTILS_CONCRETE_OBSERVABLE_HH
#define EOS_GUARD_SRC_UTILS_CONCRETE_OBSERVABLE_HH 1

#include <eos/observable.hh>
#include <eos/utils/apply.hh>
#include <eos/utils/join.hh>
#include <eos/utils/tuple-maker.hh>

#include <array>
#include <functional>
#include <string>

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
            }

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return apply(_function, values);
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

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::array<std::string, sizeof...(Args_)> _kinematics_names_array;

        public:
            ConcreteObservableEntry(const QualifiedName & name,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _function(function),
                _kinematics_names(kinematics_names),
                _kinematics_names_array(impl::make_array<std::string>(kinematics_names))
            {
            }

            ~ConcreteObservableEntry()
            {
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new ConcreteObservable<Decay_, Args_ ...>(_name, parameters, kinematics, options, _function, _kinematics_names));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: regular observable" << std::endl;

                if (sizeof...(Args_) > 0)
                {
                    os << "    kinematic variables: " << join(std::begin(_kinematics_names_array), std::end(_kinematics_names_array)) << std::endl;
                }

                return os;
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableEntry * make_concrete_observable_entry(const QualifiedName & name, double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names = std::make_tuple())
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return new ConcreteObservableEntry<Decay_, Args_ ...>(name,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(function)),
                kinematics_names);
    }

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableRatio :
        public Observable
    {
        public:

        private:
            std::string _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            Decay_ _decay;

            std::function<double (const Decay_ *, const Args_ & ...)> _numerator, _denominator;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple;

        public:
            ConcreteObservableRatio(const std::string & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & numerator,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & denominator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _decay(parameters, options),
                _numerator(numerator),
                _denominator(denominator),
                _kinematics_names(kinematics_names),
                _argument_tuple(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names, &_decay))
            {
                uses(_decay);
            }

            ~ConcreteObservableRatio() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values = _argument_tuple;

                return apply(_numerator, values) / apply(_denominator, values);
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
                return ObservablePtr(new ConcreteObservableRatio(_name, _parameters.clone(), _kinematics.clone(), _options, _numerator, _denominator, _kinematics_names));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteObservableRatio(_name, parameters, _kinematics.clone(), _options, _numerator, _denominator, _kinematics_names));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableRatioEntry :
        public ObservableEntry
    {
        private:
            std::string _name;

            std::function<double (const Decay_ *, const Args_ & ...)> _numerator, _denominator;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

        public:
            ConcreteObservableRatioEntry(const std::string & name,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & numerator,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & denominator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _numerator(numerator),
                _denominator(denominator),
                _kinematics_names(kinematics_names)
            {
            }

            ~ConcreteObservableRatioEntry() = default;

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new ConcreteObservableRatio<Decay_, Args_ ...>(_name, parameters, kinematics, options, _numerator, _denominator, _kinematics_names));
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableEntry * make_concrete_observable_ratio_entry(const std::string & name,
            double (Decay_::* numerator)(const Args_ & ...) const,
            double (Decay_::* denominator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names = std::make_tuple())
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return new ConcreteObservableRatioEntry<Decay_, Args_ ...>(name,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(numerator)),
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(denominator)),
                kinematics_names);
    }

}

#endif
