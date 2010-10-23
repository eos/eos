/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_FACTORY_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_FACTORY_HH 1

#include <src/utils/observable.hh>

namespace eos
{
    class RareBFactory
    {
        public:
            static ObservablePtr make(const std::string & name, const Parameters & parameters, const ObservableOptions & options);
    };
}

#endif
