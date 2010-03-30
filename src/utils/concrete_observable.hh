/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_CONCRETE_OBSERVABLE_HH
#define WFITTER_GUARD_SRC_UTILS_CONCRETE_OBSERVABLE_HH 1

#include <src/utils/observable.hh>

#include <string>
#include <tr1/functional>

namespace wf
{
    template <typename Decay_, unsigned n_>
    struct ConcreteObservableData;

    template <typename Decay_>
    struct ConcreteObservableData<Decay_, 0>
    {
        typedef std::tr1::function<double (const Decay_ &)> Type;

        Type function;

        std::string name;

        ConcreteObservableData(const std::string & name, const Type & function) :
            function(function),
            name(name)
        {
        }

        double evaluate(const Decay_ & decay, const Kinematics &) const
        {
            return function(decay);
        }
    };

    template <typename Decay_>
    struct ConcreteObservableData<Decay_, 1>
    {
        typedef std::tr1::function<double (const Decay_ &, const double &)> Type;

        Type function;

        std::string name;

        std::string variable1;

        ConcreteObservableData(const std::string & name, const Type & function, const std::string & variable1) :
            function(function),
            name(name),
            variable1(variable1)
        {
        }

        double evaluate(const Decay_ & decay, const Kinematics & k) const
        {
            return function(decay, k[variable1]);
        }
    };

    template <typename Decay_>
    struct ConcreteObservableData<Decay_, 2>
    {
        typedef std::tr1::function<double (const Decay_ &, const double &, const double &)> Type;

        Type function;

        std::string name;

        std::string variable1;

        std::string variable2;

        ConcreteObservableData(const std::string & name, const Type & function, const std::string & variable1,
                const std::string & variable2) :
            function(function),
            name(name),
            variable1(variable1),
            variable2(variable2)
        {
        }

        double evaluate(const Decay_ & decay, const Kinematics & k) const
        {
            return function(decay, k[variable1], k[variable2]);
        }
    };

    template <typename Decay_, unsigned n_>
    class ConcreteObservable :
        public Observable
    {
        private:
            Parameters _parameters;

            ObservableOptions _options;

            Decay_ _decay;

            ConcreteObservableData<Decay_, n_> _data;

        public:
            ConcreteObservable(const Parameters & parameters, const ObservableOptions & options, const ConcreteObservableData<Decay_, n_> & data) :
                _parameters(parameters),
                _options(options),
                _decay(parameters, options),
                _data(data)
            {
            }

            virtual const std::string & name() const
            {
                return _data.name;
            }

            virtual double evaluate(const Kinematics & k) const
            {
                return _data.evaluate(_decay, k);
            };

            virtual Parameters parameters()
            {
                return _parameters;
            };

            virtual ObservablePtr clone() const
            {
                return ObservablePtr(new ConcreteObservable(_parameters.clone(), _options, _data));
            }
    };

    template <typename Decay_, unsigned n_>
    class ConcreteObservableFactory :
        public ObservableFactory
    {
        private:
            ConcreteObservableData<Decay_, n_> _data;

        public:
            ConcreteObservableFactory(const ConcreteObservableData<Decay_, n_> & data) :
                _data(data)
            {
            }

            ~ConcreteObservableFactory()
            {
            }

            virtual std::tr1::shared_ptr<Observable> make(const Parameters & parameters, const ObservableOptions & options) const
            {
                return std::tr1::shared_ptr<Observable>(new ConcreteObservable<Decay_, n_>(parameters, options, _data));
            }
    };
}

#endif
