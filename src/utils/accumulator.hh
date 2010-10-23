/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_ACCUMULATOR_HH
#define EOS_GUARD_SRC_UTILS_ACCUMULATOR_HH 1

#include <src/utils/private_implementation_pattern.hh>

namespace eos
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
