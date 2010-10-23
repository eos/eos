/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_DESTRINGIFY_HH
#define EOS_GUARD_SRC_UTILS_DESTRINGIFY_HH 1

#include <sstream>

namespace eos
{
    namespace impl
    {
        template <typename T_>
        struct SimpleDestringify
        {
            static T_ destringify(const std::string & input)
            {
                std::stringstream ss(input);
                T_ value;

                ss >> value;

                return value;
            }
        };

        template <typename T_> struct DoDestringify;

        template <> struct DoDestringify<double> : public SimpleDestringify<double> {};

        template <> struct DoDestringify<unsigned> : public SimpleDestringify<double> {};

        template <> struct DoDestringify<bool>
        {
            static bool destringify(const std::string & input)
            {
                if ("true" == input)
                    return true;

                return false;
            }
        };
    }

    template <typename T_>
    T_
    destringify(const std::string & input)
    {
        return impl::DoDestringify<T_>::destringify(input);
    }
}

#endif
