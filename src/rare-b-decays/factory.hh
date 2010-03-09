/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_FACTORY_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_FACTORY_HH 1

#include <src/utils/observable.hh>

namespace wf
{
    class RareBFactory
    {
        public:
            static ObservablePtr make(const std::string & name, const Parameters & parameters, const ObservableOptions & options);
    };
}

#endif
