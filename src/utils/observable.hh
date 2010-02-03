/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_OBSERVABLE_HH
#define WFITTER_GUARD_SRC_UTILS_OBSERVABLE_HH 1

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

    class ObservableFactory
    {
        public:
            ObservableFactory();

            virtual ~ObservableFactory();

            virtual ObservablePtr make(const Parameters &) const = 0;
    };
}

#endif
