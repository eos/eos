/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/exception.hh>

namespace wf
{
    Exception::Exception(const std::string & message) throw () :
        _message(message)
    {
    }

    Exception::~Exception() throw ()
    {
    }

    const char *
    Exception::what() const throw ()
    {
        return _message.c_str();
    }
}
