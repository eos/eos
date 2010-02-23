/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_OBSERVABLE_HH
#define WFITTER_GUARD_SRC_UTILS_OBSERVABLE_HH 1

#include <src/utils/exception.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/parameters.hh>

#include <string>
#include <tr1/memory>

namespace wf
{
    class Observable
    {
        public:
            Observable();

            virtual ~Observable();

            virtual const std::string & name() const = 0;

            virtual double evaluate(const Kinematics &) const = 0;
    };

    typedef std::tr1::shared_ptr<Observable> ObservablePtr;

    struct UnknownOptionError :
        public Exception
    {
        UnknownOptionError(const std::string & key) throw ();
    };

    class ObservableOptions :
        public PrivateImplementationPattern<ObservableOptions>
    {
        public:
            ObservableOptions();

            ~ObservableOptions();

            const std::string & operator[] (const std::string & key) const;

            bool has(const std::string & key) const;

            void set(const std::string & key, const std::string & value = "");
    };

    class ObservableFactory
    {
        public:
            ObservableFactory();

            virtual ~ObservableFactory();

            virtual ObservablePtr make(const Parameters &) const = 0;
    };
}

#endif
