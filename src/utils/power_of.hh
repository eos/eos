/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_POWER_OF_HH
#define EOS_GUARD_SRC_UTILS_POWER_OF_HH 1

#include <complex>

namespace eos
{
    namespace impl
    {
        template <unsigned n_, typename T_>
        struct PowerOf
        {
            inline static T_ calculate(const T_ & x)
            {
                return x * PowerOf<n_ - 1, T_>::calculate(x);
            }
        };

        template <typename T_>
        struct PowerOf<0, T_>
        {
            inline static T_ calculate(const T_ &)
            {
                return T_() + 1.0;
            }
        };
    }

    template <unsigned n_, typename T_>
    T_ power_of(const T_ & x)
    {
        return impl::PowerOf<n_, T_>::calculate(x);
    }
}

#endif
