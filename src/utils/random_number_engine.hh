/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WILSON_FITTER_GUARD_SRC_UTILS_RANDOM_NUMBER_ENGINE_HH
#define WILSON_FITTER_GUARD_SRC_UTILS_RANDOM_NUMBER_ENGINE_HH 1

#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    class RandomNumberEngine :
        public PrivateImplementationPattern<RandomNumberEngine>
    {
        public:
            RandomNumberEngine();

            ~RandomNumberEngine();

            unsigned int operator() (void);
    };
}

#endif
