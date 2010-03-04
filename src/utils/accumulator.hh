/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_ACCUMULATOR_HH
#define WFITTER_GUARD_SRC_UTILS_ACCUMULATOR_HH 1

#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    class Accumulator :
        public PrivateImplementationPattern<Accumulator>
    {
        public:
            Accumulator(const double & value);

            ~Accumulator();

            double operator+= (const double & addend);

            double value() const;
    };
}

#endif
