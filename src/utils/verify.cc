/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/verify.hh>

namespace eos
{
    VerifiedRangeError::VerifiedRangeError(const std::string & message) :
        Exception("VerifiedRange encountered bogus value: " + message)
    {
    }

    VerifiedRangeOverflow::VerifiedRangeOverflow(const std::string & value, const std::string & maximum) :
        VerifiedRangeError("'" + value + "' larger than maximum of '" + maximum + "'")
    {
    }

    VerifiedRangeUnderflow::VerifiedRangeUnderflow(const std::string & value, const std::string & maximum) :
        VerifiedRangeError("'" + value + "' smaller than minimum of '" + maximum + "'")
    {
    }
}

