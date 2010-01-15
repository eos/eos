/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WILSON_FITTER_GUARD_UTILS_PRIVATE_IMPLEMENTATION_PATTERN_HH
#define WILSON_FITTER_GUARD_UTILS_PRIVATE_IMPLEMENTATION_PATTERN_HH 1

#include <tr1/memory>

namespace wf
{
    template <typename T_> struct Implementation;

    template <typename T_> class PrivateImplementationPattern
    {
        protected:
            std::tr1::shared_ptr<Implementation<T_> > _imp;

        public:
            explicit PrivateImplementationPattern(Implementation<T_> * imp);

            ~PrivateImplementationPattern();
    };
}

#endif
