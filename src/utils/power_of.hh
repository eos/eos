/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_POWER_OF_HH
#define EOS_GUARD_SRC_UTILS_POWER_OF_HH 1

namespace eos
{
    namespace impl
    {
        template <unsigned n_>
        struct PowerOf
        {
            inline static double calculate(const double & x)
            {
                return x * PowerOf<n_ - 1>::calculate(x);
            }
        };

        template <>
        struct PowerOf<0>
        {
            inline static double calculate(const double &)
            {
                return 1.0;
            }
        };
    }

    template <unsigned n_>
    double power_of(const double & x)
    {
        return impl::PowerOf<n_>::calculate(x);
    }
}

#endif
