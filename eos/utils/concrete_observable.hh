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

#include <eos/observable-impl.hh>
#include <eos/utils/apply.hh>
#include <eos/utils/join.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

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

            std::string _latex;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array;

            Options _forced_options;

        public:
            ConcreteObservableEntry(const QualifiedName & name, const std::string & latex,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names,
                    const Options & forced_options) :
                _name(name),
                _latex(latex),
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

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const
            {
                return _kinematics_names_array.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const
            {
                return _kinematics_names_array.end();
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new ConcreteObservable<Decay_, Args_ ...>(_name, parameters, kinematics, options + _forced_options, _function, _kinematics_names));
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
    ObservableEntryPtr make_concrete_observable_entry(const QualifiedName & name, const std::string & latex,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names,
            const Options & forced_options)
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return std::make_shared<ConcreteObservableEntry<Decay_, Args_ ...>>(name, latex,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(function)),
                kinematics_names, forced_options);
    }

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableRatio :
        public Observable
    {
        public:

        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options, _forced_options_numerator, _forced_options_denominator;

            Decay_ _decay_numerator, _decay_denominator;

            std::function<double (const Decay_ *, const Args_ & ...)> _numerator, _denominator;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names_numerator, _kinematics_names_denominator;

            std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple_numerator, _argument_tuple_denominator;

        public:
            ConcreteObservableRatio(const QualifiedName & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & numerator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_numerator,
                    const Options & forced_options_numerator,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & denominator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_denominator,
                    const Options & forced_options_denominator) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _forced_options_numerator(forced_options_numerator),
                _forced_options_denominator(forced_options_denominator),
                _decay_numerator(parameters, options + _forced_options_numerator),
                _decay_denominator(parameters, options + _forced_options_denominator),
                _numerator(numerator),
                _denominator(denominator),
                _kinematics_names_numerator(kinematics_names_numerator),
                _kinematics_names_denominator(kinematics_names_denominator),
                _argument_tuple_numerator(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names_numerator, &_decay_numerator)),
                _argument_tuple_denominator(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names_denominator, &_decay_denominator))
            {
                uses(_decay_numerator);
                uses(_decay_denominator);
            }

            ~ConcreteObservableRatio() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values_numerator   = _argument_tuple_numerator;
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values_denominator = _argument_tuple_denominator;

                return apply(_numerator, values_numerator) / apply(_denominator, values_denominator);
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
                return ObservablePtr(new ConcreteObservableRatio(_name, _parameters.clone(), _kinematics.clone(), _options,
                        _numerator,   _kinematics_names_numerator,   _forced_options_numerator,
                        _denominator, _kinematics_names_denominator, _forced_options_denominator));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteObservableRatio(_name, parameters, _kinematics.clone(), _options,
                        _numerator,   _kinematics_names_numerator,   _forced_options_numerator,
                        _denominator, _kinematics_names_denominator, _forced_options_denominator));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableRatioEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            std::function<double (const Decay_ *, const Args_ & ...)> _numerator, _denominator;

            Options _forced_options_numerator, _forced_options_denominator;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names_numerator, _kinematics_names_denominator;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array_numerator, _kinematics_names_array_denominator;

        public:
            ConcreteObservableRatioEntry(const QualifiedName & name, const std::string & latex,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & numerator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_numerator,
                    const Options & forced_options_numerator,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & denominator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_denominator,
                    const Options & forced_options_denominator) :
                _name(name),
                _latex(latex),
                _numerator(numerator),
                _denominator(denominator),
                _forced_options_numerator(forced_options_numerator),
                _forced_options_denominator(forced_options_denominator),
                _kinematics_names_numerator(kinematics_names_numerator),
                _kinematics_names_denominator(kinematics_names_denominator),
                _kinematics_names_array_numerator(impl::make_array<const std::string>(kinematics_names_numerator)),
                _kinematics_names_array_denominator(impl::make_array<const std::string>(kinematics_names_denominator))
            {
            }

            ~ConcreteObservableRatioEntry() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual const std::string & latex() const
            {
                return _latex;
            }

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const
            {
                return _kinematics_names_array_numerator.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const
            {
                return _kinematics_names_array_numerator.end();
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new ConcreteObservableRatio<Decay_, Args_ ...>(_name, parameters, kinematics, options,
                        _numerator,   _kinematics_names_numerator,   _forced_options_numerator,
                        _denominator, _kinematics_names_denominator, _forced_options_denominator));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: observable ratio" << std::endl;

                if (sizeof...(Args_) > 0)
                {
                    os << "    kinematic variables numerator:   " << join(std::begin(_kinematics_names_array_numerator),   std::end(_kinematics_names_array_numerator))   << std::endl;
                    os << "    kinematic variables denominator: " << join(std::begin(_kinematics_names_array_denominator), std::end(_kinematics_names_array_denominator)) << std::endl;
                }

                return os;
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableEntryPtr make_concrete_observable_ratio_entry(const QualifiedName & name, const std::string & latex,
            double (Decay_::* numerator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_numerator,
            const Options & forced_options_numerator,
            double (Decay_::* denominator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_denominator,
            const Options & forced_options_denominator)
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return std::make_shared<ConcreteObservableRatioEntry<Decay_, Args_ ...>>(name, latex,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(numerator)),
                kinematics_names_numerator,
                forced_options_numerator,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(denominator)),
                kinematics_names_denominator,
                forced_options_denominator
                );
    }



    template <typename Decay_, typename ... Args_>
    class ConcreteObservableSum :
        public Observable
    {
        public:

        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options, _forced_options_numerator, _forced_options_denominator;

            Decay_ _decay_numerator, _decay_denominator;

            std::function<double (const Decay_ *, const Args_ & ...)> _numerator, _denominator;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names_numerator, _kinematics_names_denominator;

            std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, KinematicVariable>::Type ...> _argument_tuple_numerator, _argument_tuple_denominator;

            const double _weight_numerator, _weight_denominator;

        public:
            ConcreteObservableSum(const QualifiedName & name,
                    const Parameters & parameters,
                    const Kinematics & kinematics,
                    const Options & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & numerator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_numerator,
                    const Options & forced_options_numerator,
                    const double & weight_numerator,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & denominator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_denominator,
                    const Options & forced_options_denominator,
                    const double & weight_denominator) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _forced_options_numerator(forced_options_numerator),
                _forced_options_denominator(forced_options_denominator),
                _decay_numerator(parameters, options + _forced_options_numerator),
                _decay_denominator(parameters, options + _forced_options_denominator),
                _numerator(numerator),
                _denominator(denominator),
                _kinematics_names_numerator(kinematics_names_numerator),
                _kinematics_names_denominator(kinematics_names_denominator),
                _argument_tuple_numerator(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names_numerator, &_decay_numerator)),
                _argument_tuple_denominator(impl::TupleMaker<sizeof...(Args_)>::make(_kinematics, _kinematics_names_denominator, &_decay_denominator)),
                _weight_numerator(weight_numerator),
                _weight_denominator(weight_denominator)
            {
                uses(_decay_numerator);
                uses(_decay_denominator);
            }

            ~ConcreteObservableSum() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual double evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values_numerator   = _argument_tuple_numerator;
                std::tuple<const Decay_ *, typename impl::ConvertTo<Args_, double>::Type ...> values_denominator = _argument_tuple_denominator;

                return _weight_numerator * apply(_numerator, values_numerator) + _weight_denominator * apply(_denominator, values_denominator);
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
                return ObservablePtr(new ConcreteObservableSum(_name, _parameters.clone(), _kinematics.clone(), _options,
                        _numerator,   _kinematics_names_numerator,   _forced_options_numerator, _weight_numerator,
                        _denominator, _kinematics_names_denominator, _forced_options_denominator, _weight_denominator));
            }

            virtual ObservablePtr clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteObservableSum(_name, parameters, _kinematics.clone(), _options,
                        _numerator,   _kinematics_names_numerator,   _forced_options_numerator, _weight_numerator,
                        _denominator, _kinematics_names_denominator, _forced_options_denominator, _weight_denominator));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableSumEntry :
        public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            std::function<double (const Decay_ *, const Args_ & ...)> _numerator, _denominator;

            Options _forced_options_numerator, _forced_options_denominator;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names_numerator, _kinematics_names_denominator;

            std::array<const std::string, sizeof...(Args_)> _kinematics_names_array_numerator, _kinematics_names_array_denominator;

            const double _weight_numerator, _weight_denominator;

        public:
            ConcreteObservableSumEntry(const QualifiedName & name, const std::string & latex,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & numerator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_numerator,
                    const Options & forced_options_numerator,
                    const double & weight_numerator,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & denominator,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names_denominator,
                    const Options & forced_options_denominator,
                    const double & weight_denominator) :
                _name(name),
                _latex(latex),
                _numerator(numerator),
                _denominator(denominator),
                _forced_options_numerator(forced_options_numerator),
                _forced_options_denominator(forced_options_denominator),
                _kinematics_names_numerator(kinematics_names_numerator),
                _kinematics_names_denominator(kinematics_names_denominator),
                _kinematics_names_array_numerator(impl::make_array<const std::string>(kinematics_names_numerator)),
                _kinematics_names_array_denominator(impl::make_array<const std::string>(kinematics_names_denominator)),
                _weight_numerator(weight_numerator),
                _weight_denominator(weight_denominator)
            {
            }

            ~ConcreteObservableSumEntry() = default;

            virtual const QualifiedName & name() const
            {
                return _name;
            }

            virtual const std::string & latex() const
            {
                return _latex;
            }

            virtual ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const
            {
                return _kinematics_names_array_numerator.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator end_kinematic_variables() const
            {
                return _kinematics_names_array_numerator.end();
            }

            virtual ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                return ObservablePtr(new ConcreteObservableSum<Decay_, Args_ ...>(_name, parameters, kinematics, options,
                        _numerator,   _kinematics_names_numerator,   _forced_options_numerator, _weight_numerator,
                        _denominator, _kinematics_names_denominator, _forced_options_denominator, _weight_denominator));
            }

            virtual std::ostream & insert(std::ostream & os) const
            {
                os << "    type: observable sum" << std::endl;

                if (sizeof...(Args_) > 0)
                {
                    os << "    kinematic variables numerator:   " << join(std::begin(_kinematics_names_array_numerator),   std::end(_kinematics_names_array_numerator))   << std::endl;
                    os << "    kinematic variables denominator: " << join(std::begin(_kinematics_names_array_denominator), std::end(_kinematics_names_array_denominator)) << std::endl;
                }

                return os;
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableEntryPtr make_concrete_observable_sum_entry(const QualifiedName & name, const std::string & latex,
            double (Decay_::* numerator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_numerator,
            const Options & forced_options_numerator,
            const double & weight_numerator,
            double (Decay_::* denominator)(const Args_ & ...) const,
            const Tuple_ & kinematics_names_denominator,
            const Options & forced_options_denominator,
            const double & weight_denominator)
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return std::make_shared<ConcreteObservableSumEntry<Decay_, Args_ ...>>(name, latex,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(numerator)),
                kinematics_names_numerator,
                forced_options_numerator,
                weight_numerator,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(denominator)),
                kinematics_names_denominator,
                forced_options_denominator,
                weight_denominator
                );
    }
}


#endif
