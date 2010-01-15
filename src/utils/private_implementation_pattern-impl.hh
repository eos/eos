/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WILSON_FITTER_GUARD_UTILS_PRIVATE_IMPLEMENTATION_PATTERN_IMPL_HH
#define WILSON_FITTER_GUARD_UTILS_PRIVATE_IMPLEMENTATION_PATTERN_IMPL_HH 1

#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    template <typename T_>
    PrivateImplementationPattern<T_>::PrivateImplementationPattern(Implementation<T_> * imp) :
        _imp(imp)
    {
    }

    template <typename T_>
    PrivateImplementationPattern<T_>::~PrivateImplementationPattern()
    {
    }
}

#endif
