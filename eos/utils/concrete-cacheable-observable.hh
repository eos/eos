/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_CONCRETE_CACHEABLE_OBSERVABLE_HH
#define EOS_GUARD_EOS_UTILS_CONCRETE_CACHEABLE_OBSERVABLE_HH 1

#include <eos/observable-impl.hh>
#include <eos/utils/join.hh>
#include <eos/utils/log.hh>
#include <eos/utils/tuple-maker.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <array>
#include <functional>
#include <string>
#include <tuple>
#include <typeindex>
#include <utility>

namespace eos
{
    template <typename Decay_, typename IntermediateResult_, typename PrepareArgs_, typename EvaluateArgs_> class ConcreteCacheableObservable;

    template <typename Decay_, typename IntermediateResult_, typename PrepareArgs_, typename EvaluateArgs_> class ConcreteCachedObservable;

    /*!
     * CacheablePreparer holds the half of a cacheable observable that produces the intermediate
     * result. Observables that agree on this half share one intermediate result, whichever
     * kinematic variables their evaluation functions consume on top of it.
     */
    template <typename Decay_, typename IntermediateResult_, typename PrepareArgs_> class CacheablePreparer;

    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_>
    class CacheablePreparer<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>> : public CacheableObservable
    {
        public:
            using PrepareFunction = const IntermediateResult_ * (Decay_::*) (const PrepareArgs_ &...) const;

            using PrepareKinematicsNames = std::tuple<typename impl::ConvertTo<PrepareArgs_, const char *>::Type...>;

        protected:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            std::shared_ptr<Decay_> _decay;

            PrepareFunction _prepare_fn;

            PrepareKinematicsNames _prepare_kinematics_names;

            std::tuple<const Decay_ *, typename impl::ConvertTo<PrepareArgs_, KinematicVariable>::Type...> _prepare_argument_tuple;

            CacheablePreparer(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options, PrepareFunction prepare_fn,
                              const PrepareKinematicsNames & prepare_kinematics_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _decay(new Decay_(parameters, options)),
                _prepare_fn(prepare_fn),
                _prepare_kinematics_names(prepare_kinematics_names),
                _prepare_argument_tuple(impl::TupleMaker<sizeof...(PrepareArgs_)>::make(_kinematics, _prepare_kinematics_names, _decay.get()))
            {
                uses(*_decay);
                auto _register_kinematics = [this](const Decay_ *, typename impl::ConvertTo<PrepareArgs_, KinematicVariable>::Type... args)
                {
                    std::array<const KinematicVariable, sizeof...(PrepareArgs_)> kinematics_array = { args... };
                    for (const auto & kinematic_variable : kinematics_array)
                    {
                        this->uses_kinematic(kinematic_variable.id());
                    }
                };
                std::apply(_register_kinematics, _prepare_argument_tuple);
                uses(Decay_::references);
            }

        public:
            ~CacheablePreparer() = default;

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual Parameters
            parameters()
            {
                return _parameters;
            }

            virtual Kinematics
            kinematics()
            {
                return _kinematics;
            }

            virtual Options
            options()
            {
                return _options;
            }

            virtual std::type_index
            prepare_type_index() const
            {
                return std::type_index(typeid(CacheablePreparer));
            }

            virtual const CacheableObservable::IntermediateResult *
            prepare() const
            {
                return this->prepared_result();
            }

            const std::shared_ptr<Decay_> &
            decay() const
            {
                return _decay;
            }

            const IntermediateResult_ *
            prepared_result() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<PrepareArgs_, double>::Type...> values = _prepare_argument_tuple;

                return std::apply(std::mem_fn(_prepare_fn), values);
            }

            /// Can this observable adopt the intermediate result prepared by another one?
            bool
            shares_prepared_result_with(const CacheablePreparer & other) const
            {
                if (_prepare_fn != other._prepare_fn)
                {
                    return false;
                }

                if (_parameters != other._parameters)
                {
                    return false;
                }

                if (_options != other._options)
                {
                    return false;
                }

                return _same_prepare_arguments(other, std::make_index_sequence<sizeof...(PrepareArgs_)>{});
            }

        private:
            template <std::size_t... Indices_>
            bool
            _same_prepare_arguments(const CacheablePreparer & other, std::index_sequence<Indices_...>) const
            {
                return ((std::get<Indices_ + 1u>(_prepare_argument_tuple).evaluate() == std::get<Indices_ + 1u>(other._prepare_argument_tuple).evaluate()) && ...);
            }
    };

    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_, typename... EvaluateArgs_>
    class ConcreteCachedObservable<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>> : public Observable
    {
        public:
            using PrepareFunction  = const IntermediateResult_ * (Decay_::*) (const PrepareArgs_ &...) const;
            using EvaluateFunction = double (Decay_::*)(const IntermediateResult_ *, const EvaluateArgs_ &...) const;

            using PrepareKinematicsNames  = std::tuple<typename impl::ConvertTo<PrepareArgs_, const char *>::Type...>;
            using EvaluateKinematicsNames = std::tuple<typename impl::ConvertTo<EvaluateArgs_, const char *>::Type...>;

        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            std::shared_ptr<Decay_> _decay;

            const IntermediateResult_ * _intermediate_result;

            PrepareFunction _prepare_fn;

            EvaluateFunction _evaluate_fn;

            PrepareKinematicsNames _prepare_kinematics_names;

            EvaluateKinematicsNames _evaluate_kinematics_names;

            std::tuple<const Decay_ *, typename impl::ConvertTo<EvaluateArgs_, KinematicVariable>::Type...> _evaluate_argument_tuple;

        public:
            ConcreteCachedObservable(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                     const std::shared_ptr<Decay_> & decay, const IntermediateResult_ * intermediate_result, PrepareFunction prepare_fn,
                                     EvaluateFunction evaluate_fn, const PrepareKinematicsNames & prepare_kinematics_names,
                                     const EvaluateKinematicsNames & evaluate_kinematics_names) :
                _name(name),
                _parameters(parameters),
                _kinematics(kinematics),
                _options(options),
                _decay(decay),
                _intermediate_result(intermediate_result),
                _prepare_fn(prepare_fn),
                _evaluate_fn(evaluate_fn),
                _prepare_kinematics_names(prepare_kinematics_names),
                _evaluate_kinematics_names(evaluate_kinematics_names),
                _evaluate_argument_tuple(impl::TupleMaker<sizeof...(EvaluateArgs_)>::make(_kinematics, _evaluate_kinematics_names, _decay.get()))
            {
                uses(*_decay);
                auto _register_kinematics = [this](const Decay_ *, const auto &... args) { (this->uses_kinematic(args.id()), ...); };
                // the preparation half's variables are used by the observable that prepared our result, not by us;
                // register them all the same, so that we compare equal to an identical uncached observable
                std::apply(_register_kinematics, impl::TupleMaker<sizeof...(PrepareArgs_)>::make(_kinematics, _prepare_kinematics_names, _decay.get()));
                std::apply(_register_kinematics, _evaluate_argument_tuple);
                uses(Decay_::references);
            }

            ~ConcreteCachedObservable() = default;

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual double
            evaluate() const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<EvaluateArgs_, double>::Type...> values = _evaluate_argument_tuple;

                return std::apply([evaluate_fn = _evaluate_fn, intermediate_result = _intermediate_result](const Decay_ * decay,
                                                                                                           const typename impl::ConvertTo<EvaluateArgs_, double>::Type &... args)
                { return (decay->*evaluate_fn)(intermediate_result, args...); },
                                  values);
            }

            virtual Parameters
            parameters()
            {
                return _parameters;
            }

            virtual Kinematics
            kinematics()
            {
                return _kinematics;
            }

            virtual Options
            options()
            {
                return _options;
            }

            virtual ObservablePtr
            clone() const
            {
                return ObservablePtr(
                        new ConcreteCacheableObservable<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>>(_name,
                                                                                                                                                _parameters.clone(),
                                                                                                                                                _kinematics.clone(),
                                                                                                                                                _options,
                                                                                                                                                _prepare_fn,
                                                                                                                                                _evaluate_fn,
                                                                                                                                                _prepare_kinematics_names,
                                                                                                                                                _evaluate_kinematics_names));
            }

            virtual ObservablePtr
            clone(const Parameters & parameters) const
            {
                return ObservablePtr(
                        new ConcreteCacheableObservable<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>>(_name,
                                                                                                                                                parameters,
                                                                                                                                                _kinematics.clone(),
                                                                                                                                                _options,
                                                                                                                                                _prepare_fn,
                                                                                                                                                _evaluate_fn,
                                                                                                                                                _prepare_kinematics_names,
                                                                                                                                                _evaluate_kinematics_names));
            }
    };

    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_, typename... EvaluateArgs_>
    class ConcreteCacheableObservable<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>> :
        public CacheablePreparer<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>>
    {
        public:
            using Preparer         = CacheablePreparer<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>>;
            using PrepareFunction  = typename Preparer::PrepareFunction;
            using EvaluateFunction = double (Decay_::*)(const IntermediateResult_ *, const EvaluateArgs_ &...) const;

            using PrepareKinematicsNames  = typename Preparer::PrepareKinematicsNames;
            using EvaluateKinematicsNames = std::tuple<typename impl::ConvertTo<EvaluateArgs_, const char *>::Type...>;

        private:
            EvaluateFunction _evaluate_fn;

            EvaluateKinematicsNames _evaluate_kinematics_names;

            std::tuple<const Decay_ *, typename impl::ConvertTo<EvaluateArgs_, KinematicVariable>::Type...> _evaluate_argument_tuple;

        public:
            ConcreteCacheableObservable(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                        PrepareFunction prepare_fn, EvaluateFunction evaluate_fn, const PrepareKinematicsNames & prepare_kinematics_names,
                                        const EvaluateKinematicsNames & evaluate_kinematics_names) :
                Preparer(name, parameters, kinematics, options, prepare_fn, prepare_kinematics_names),
                _evaluate_fn(evaluate_fn),
                _evaluate_kinematics_names(evaluate_kinematics_names),
                _evaluate_argument_tuple(impl::TupleMaker<sizeof...(EvaluateArgs_)>::make(this->_kinematics, _evaluate_kinematics_names, this->_decay.get()))
            {
                auto _register_kinematics = [this](const Decay_ *, typename impl::ConvertTo<EvaluateArgs_, KinematicVariable>::Type... args)
                {
                    std::array<const KinematicVariable, sizeof...(EvaluateArgs_)> kinematics_array = { args... };
                    for (const auto & kinematic_variable : kinematics_array)
                    {
                        this->uses_kinematic(kinematic_variable.id());
                    }
                };
                std::apply(_register_kinematics, _evaluate_argument_tuple);
            }

            ~ConcreteCacheableObservable() = default;

            virtual double
            evaluate() const
            {
                return this->_evaluate(this->prepared_result());
            }

            virtual double
            evaluate(const CacheableObservable::IntermediateResult * intermediate_result) const
            {
                return this->_evaluate(static_cast<const IntermediateResult_ *>(intermediate_result));
            }

            virtual ObservablePtr
            make_cached_observable(const CacheableObservable * _other) const
            {
                // the other observable need not agree on the evaluation half
                auto other = dynamic_cast<const Preparer *>(_other);
                if (nullptr == other)
                {
                    return { nullptr };
                }

                if (! this->shares_prepared_result_with(*other))
                {
                    return { nullptr };
                }

                /*
                 * The intermediate result is owned by the other observable's provider, which we keep
                 * alive by sharing its pointer.
                 */
                return ObservablePtr(
                        new ConcreteCachedObservable<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>>(this->_name,
                                                                                                                                             this->_parameters,
                                                                                                                                             this->_kinematics,
                                                                                                                                             this->_options,
                                                                                                                                             other->decay(),
                                                                                                                                             other->prepared_result(),
                                                                                                                                             this->_prepare_fn,
                                                                                                                                             _evaluate_fn,
                                                                                                                                             this->_prepare_kinematics_names,
                                                                                                                                             _evaluate_kinematics_names));
            }

            virtual ObservablePtr
            clone() const
            {
                return ObservablePtr(new ConcreteCacheableObservable(this->_name,
                                                                     this->_parameters.clone(),
                                                                     this->_kinematics.clone(),
                                                                     this->_options,
                                                                     this->_prepare_fn,
                                                                     _evaluate_fn,
                                                                     this->_prepare_kinematics_names,
                                                                     _evaluate_kinematics_names));
            }

            virtual ObservablePtr
            clone(const Parameters & parameters) const
            {
                return ObservablePtr(new ConcreteCacheableObservable(this->_name,
                                                                     parameters,
                                                                     this->_kinematics.clone(),
                                                                     this->_options,
                                                                     this->_prepare_fn,
                                                                     _evaluate_fn,
                                                                     this->_prepare_kinematics_names,
                                                                     _evaluate_kinematics_names));
            }

        private:
            double
            _evaluate(const IntermediateResult_ * intermediate_result) const
            {
                std::tuple<const Decay_ *, typename impl::ConvertTo<EvaluateArgs_, double>::Type...> values = _evaluate_argument_tuple;

                return std::apply([evaluate_fn = _evaluate_fn, intermediate_result](const Decay_ * decay, const typename impl::ConvertTo<EvaluateArgs_, double>::Type &... args)
                { return (decay->*evaluate_fn)(intermediate_result, args...); },
                                  values);
            }
    };

    template <typename Decay_, typename IntermediateResult_, typename PrepareArgs_, typename EvaluateArgs_> class ConcreteCacheableObservableEntry;

    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_, typename... EvaluateArgs_>
    class ConcreteCacheableObservableEntry<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>> : public ObservableEntry
    {
        public:
            using Observable_      = ConcreteCacheableObservable<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>>;
            using PrepareFunction  = typename Observable_::PrepareFunction;
            using EvaluateFunction = typename Observable_::EvaluateFunction;

            using PrepareKinematicsNames  = typename Observable_::PrepareKinematicsNames;
            using EvaluateKinematicsNames = typename Observable_::EvaluateKinematicsNames;

        private:
            QualifiedName _name;

            std::string _latex;

            Unit _unit;

            PrepareFunction _prepare_fn;

            EvaluateFunction _evaluate_fn;

            PrepareKinematicsNames _prepare_kinematics_names;

            EvaluateKinematicsNames _evaluate_kinematics_names;

            std::array<const std::string, sizeof...(PrepareArgs_)> _prepare_kinematics_names_array;

            std::array<const std::string, sizeof...(EvaluateArgs_)> _evaluate_kinematics_names_array;

            std::array<const std::string, sizeof...(PrepareArgs_) + sizeof...(EvaluateArgs_)> _kinematics_names_array;

            Options _forced_options;

        public:
            ConcreteCacheableObservableEntry(const QualifiedName & name, const std::string & latex, const Unit & unit, PrepareFunction prepare_fn, EvaluateFunction evaluate_fn,
                                             const PrepareKinematicsNames & prepare_kinematics_names, const EvaluateKinematicsNames & evaluate_kinematics_names,
                                             const Options & forced_options) :
                _name(name),
                _latex(latex),
                _unit(unit),
                _prepare_fn(prepare_fn),
                _evaluate_fn(evaluate_fn),
                _prepare_kinematics_names(prepare_kinematics_names),
                _evaluate_kinematics_names(evaluate_kinematics_names),
                _prepare_kinematics_names_array(impl::make_array<const std::string>(prepare_kinematics_names)),
                _evaluate_kinematics_names_array(impl::make_array<const std::string>(evaluate_kinematics_names)),
                _kinematics_names_array(impl::make_array<const std::string>(std::tuple_cat(prepare_kinematics_names, evaluate_kinematics_names))),
                _forced_options(forced_options)
            {
                for (const auto & prepare_name : _prepare_kinematics_names_array)
                {
                    for (const auto & evaluate_name : _evaluate_kinematics_names_array)
                    {
                        if (prepare_name == evaluate_name)
                        {
                            throw InternalError("Observable '" + _name.str() + "' declares the kinematic variable '" + prepare_name
                                                + "' both for its intermediate result and for its evaluation");
                        }
                    }
                }
            }

            ~ConcreteCacheableObservableEntry() {}

            virtual const QualifiedName &
            name() const
            {
                return _name;
            }

            virtual const std::string &
            latex() const
            {
                return _latex;
            }

            virtual const Unit &
            unit() const
            {
                return _unit;
            }

            virtual ObservableEntry::KinematicVariableIterator
            begin_kinematic_variables() const
            {
                return _kinematics_names_array.begin();
            }

            virtual ObservableEntry::KinematicVariableIterator
            end_kinematic_variables() const
            {
                return _kinematics_names_array.end();
            }

            virtual ObservableEntry::OptionIterator
            begin_options() const
            {
                return Decay_::begin_options();
            }

            virtual ObservableEntry::OptionIterator
            end_options() const
            {
                return Decay_::end_options();
            }

            virtual ObservablePtr
            make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
            {
                for (const auto & fo : _forced_options)
                {
                    const auto & key = std::get<0>(fo);
                    if (options.has(key))
                    {
                        Log::instance()->message("[ConcreteCacheableObservableEntry.make]", ll_warning)
                                << "Observable '" << _name << "' forces option key '" << key << "' to value '" << _forced_options[key] << "', overriding user-provided value '"
                                << options[key] << "'";
                    }
                }
                return ObservablePtr(new Observable_(_name,
                                                     parameters,
                                                     kinematics,
                                                     options + _forced_options,
                                                     _prepare_fn,
                                                     _evaluate_fn,
                                                     _prepare_kinematics_names,
                                                     _evaluate_kinematics_names));
            }

            virtual std::ostream &
            insert(std::ostream & os) const
            {
                os << "    type: cacheable observable" << std::endl;

                if (sizeof...(PrepareArgs_) > 0)
                {
                    os << "    cached kinematic variables: " << join(std::begin(_prepare_kinematics_names_array), std::end(_prepare_kinematics_names_array)) << std::endl;
                }

                if (sizeof...(EvaluateArgs_) > 0)
                {
                    os << "    kinematic variables: " << join(std::begin(_evaluate_kinematics_names_array), std::end(_evaluate_kinematics_names_array)) << std::endl;
                }

                return os;
            }
    };

    template <typename Decay_, typename IntermediateResult_, typename... PrepareArgs_, typename... EvaluateArgs_>
    ObservableEntryPtr
    make_concrete_cacheable_observable_entry(const QualifiedName & name, const std::string & latex, const Unit & unit,
                                             const impl::Preparer<Decay_, IntermediateResult_, PrepareArgs_...> &   preparer,
                                             const impl::Evaluator<Decay_, IntermediateResult_, EvaluateArgs_...> & evaluator, const Options & forced_options)
    {
        return std::make_shared<ConcreteCacheableObservableEntry<Decay_, IntermediateResult_, std::tuple<PrepareArgs_...>, std::tuple<EvaluateArgs_...>>>(
                name,
                latex,
                unit,
                preparer.function,
                evaluator.function,
                preparer.kinematics_names,
                evaluator.kinematics_names,
                forced_options);
    }
} // namespace eos


#endif
