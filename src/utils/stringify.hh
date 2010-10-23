/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_STRINGIFY_HH
#define EOS_GUARD_SRC_UTILS_STRINGIFY_HH 1

#include <string>
#include <sstream>

namespace eos
{
    namespace implementation
    {
        template <typename T_>
        struct DoStringify
        {
            static std::string stringify(const T_ & x)
            {
                std::stringstream ss;
                ss << x;

                return ss.str();
            }
        };

        template <>
        struct DoStringify<std::string>
        {
            static std::string stringify(const std::string & x)
            {
                return x;
            }
        };
    }

    template <typename T_>
    std::string stringify(const T_ & x)
    {
        return implementation::DoStringify<T_>::stringify(x);
    }
}

#endif
