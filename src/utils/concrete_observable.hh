/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_CONCRETE_OBSERVABLE_HH
#define WFITTER_GUARD_SRC_UTILS_CONCRETE_OBSERVABLE_HH 1

#include <src/utils/apply.hh>
#include <src/utils/observable.hh>

#include <functional>
#include <string>

namespace wf
{
    namespace impl
    {
        template <typename T_, typename U_> struct ConvertTo { typedef U_ Type; };

        template <unsigned n_> struct TupleMaker
        {
            template <typename Decay_, typename ... TupleElements_, typename ... ResultElements_>
            static auto make(const Kinematics & k, const std::tuple<TupleElements_ ...> & t, const Decay_ * d, ResultElements_ ... r)
                -> std::tuple<const Decay_ *, typename ConvertTo<TupleElements_, double>::Type ...>
            {
                return TupleMaker<n_ - 1>::make(k, t, d, k[std::get<n_ - 1>(t)], r ...);
            }
        };

        template <> struct TupleMaker<0>
        {
            template <typename Decay_, typename ... TupleElements_, typename ... ResultElements_>
            static auto make(const Kinematics & k, const std::tuple<TupleElements_ ...> & t, const Decay_ * d, ResultElements_ ... r)
                -> std::tuple<const Decay_ *, typename ConvertTo<TupleElements_, double>::Type ...>
            {
                return std::make_tuple(d, r ...);
            }
        };

        template <typename T_> struct TupleSize;

        template <typename ... TupleElements_> struct TupleSize<std::tuple<TupleElements_ ...>>
        {
            static const unsigned long size = sizeof...(TupleElements_);
        };
    }

    template <typename Decay_, typename ... Args_>
    class ConcreteObservable :
        public Observable
    {
        public:

        private:
            std::string _name;

            Parameters _parameters;

            ObservableOptions _options;

            Decay_ _decay;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

            std::tuple<const Decay_ *, Args_ ...> _argument_tuple(const Kinematics & k) const
            {
                return impl::TupleMaker<sizeof...(Args_)>::make(k, _kinematics_names, &_decay);
            }

        public:
            ConcreteObservable(const std::string & name,
                    const Parameters & parameters,
                    const ObservableOptions & options,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _parameters(parameters),
                _options(options),
                _decay(parameters, options),
                _function(function),
                _kinematics_names(kinematics_names)
            {
            }

            virtual const std::string & name() const
            {
                return _name;
            }

            virtual double evaluate(const Kinematics & k) const
            {
                return apply(_function, _argument_tuple(k));
            };

            virtual Parameters parameters()
            {
                return _parameters;
            };

            virtual ObservablePtr clone() const
            {
                return ObservablePtr(new ConcreteObservable(_name, _parameters.clone(), _options, _function, _kinematics_names));
            }
    };

    template <typename Decay_, typename ... Args_>
    class ConcreteObservableFactory :
        public ObservableFactory
    {
        private:
            std::string _name;

            std::function<double (const Decay_ *, const Args_ & ...)> _function;

            std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> _kinematics_names;

        public:
            ConcreteObservableFactory(const std::string & name,
                    const std::function<double (const Decay_ *, const Args_ & ...)> & function,
                    const std::tuple<typename impl::ConvertTo<Args_, const char *>::Type ...> & kinematics_names) :
                _name(name),
                _function(function),
                _kinematics_names(kinematics_names)
            {
            }

            ~ConcreteObservableFactory()
            {
            }

            virtual std::tr1::shared_ptr<Observable> make(const Parameters & parameters, const ObservableOptions & options) const
            {
                return std::tr1::shared_ptr<Observable>(new ConcreteObservable<Decay_, Args_ ...>(_name, parameters, options, _function,
                            _kinematics_names));
            }
    };

    template <typename Decay_, typename Tuple_, typename ... Args_>
    ObservableFactory * make_concrete_observable_factory(const std::string & name, double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names = std::make_tuple())
    {
        static_assert(sizeof...(Args_) == impl::TupleSize<Tuple_>::size, "Need as many function arguments as kinematics names!");

        return new ConcreteObservableFactory<Decay_, Args_ ...>(name,
                std::function<double (const Decay_ *, const Args_ & ...)>(std::mem_fn(function)),
                kinematics_names);
    }
}

#endif
