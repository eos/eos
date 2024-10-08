/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2022 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_UTILS_EXCEPTION_HH
#define EOS_GUARD_EOS_UTILS_EXCEPTION_HH 1

#include <exception>
#include <memory>
#include <string>
#include <source_location>

namespace eos
{
    using source_location = std::source_location;

    class Context
    {
        private:
            Context(const Context &);
            Context & operator= (const Context &);

        public:
            Context(const std::string & entry, const source_location = source_location::current());
            ~Context() noexcept(false);

            std::string backtrace(const std::string & delimiter) const;
    };

    class Exception :
        public std::exception
    {
        private:
            class ContextData;

            std::string _message;
            const std::unique_ptr<ContextData> _context_data;

        protected:
            Exception(const std::string & message) noexcept;
            Exception(const Exception &);

        public:
            ~Exception() noexcept;

            std::string backtrace(const std::string & delimiter) const;

            virtual const char * what() const noexcept override;
    };

    class InternalError :
        public Exception
    {
        public:
            InternalError(const std::string & message) throw ();
    };

    class UnknownObservableError :
        public Exception
    {
        public:
            UnknownObservableError(const std::string & message) throw ();
    };

    class GSLError :
        public Exception
    {
        public:
            GSLError(const std::string & message) throw ();
    };

    class ParsingError :
        public Exception
    {
        public:
            ParsingError(const std::string & message) throw ();
    };
}

#endif
