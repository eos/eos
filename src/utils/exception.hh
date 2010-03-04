/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_EXCEPTION_HH
#define WFITTER_GUARD_SRC_UTILS_EXCEPTION_HH 1

#include <exception>
#include <string>

namespace wf
{
    class Exception :
        public std::exception
    {
        private:
            std::string _message;

        protected:
            Exception(const std::string & message) throw ();

        public:
            ~Exception() throw ();

            virtual const char * what() const throw ();
    };

    struct InternalError :
        public Exception
    {
        InternalError(const std::string & message) throw ();
    };
}

#endif
