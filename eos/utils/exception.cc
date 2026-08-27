/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2026 Danny van Dyk
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

#include <eos/utils/exception.hh>

#include <list>
#include <utility>

namespace eos
{
    namespace impl
    {
        static thread_local std::list<std::tuple<source_location, std::string>> context;
    }

    Context::Context(const std::string & entry, const source_location location)
    {
        impl::context.emplace_back(location, entry);
    }

    Context::Context(const Context & other) = default;

    Context & Context::operator= (const Context & other) = default;

    Context::~Context() noexcept(false)
    {
        if (impl::context.empty())
        {
            throw InternalError("empty context");
        }

        impl::context.pop_back();
    }

    std::string
    Context::backtrace(const std::string & delimiter) const
    {
        if (impl::context.empty())
        {
            return "";
        }

        std::string result;

        auto append = [&](const std::tuple<source_location, std::string> & entry)
        {
            result +=
                    std::get<1>(entry) + "[" + std::get<0>(entry).file_name() + ":" + std::to_string(std::get<0>(entry).line()) + " -> " + std::get<0>(entry).function_name() + "]";
            result += delimiter;
        };

        for (const auto & c : impl::context)
        {
            append(c);
        }

        return result;
    }

    class Exception::ContextData
    {
        public:
            std::list<std::tuple<source_location, std::string>> local_context;

            ContextData()
            {
                // use assign to ensure copying each element
                local_context.assign(impl::context.begin(), impl::context.end());
            }

            ContextData(const ContextData &) = default;
            ~ContextData()                   = default;
    };

    Exception::Exception(std::string message) noexcept :
        _message(std::move(message)),
        _context_data(std::make_unique<ContextData>())
    {
    }

    Exception::Exception(const Exception & other) :
        std::exception(other),
        _message(other._message),
        _context_data(std::make_unique<ContextData>(*other._context_data))
    {
    }

    Exception::~Exception() noexcept = default;

    std::string
    Exception::backtrace(const std::string & delimiter) const
    {
        if (_context_data->local_context.empty())
        {
            return "";
        }

        std::string result;

        auto append = [&](const std::tuple<source_location, std::string> & entry)
        {
            result += std::get<1>(entry) + " [" + std::get<0>(entry).file_name() + ":" + std::to_string(std::get<0>(entry).line()) + " -> " + std::get<0>(entry).function_name()
                      + "]";
            result += delimiter;
        };

        for (const auto & c : _context_data->local_context)
        {
            append(c);
        }

        return result;
    }

    const char *
    Exception::what() const noexcept
    {
        return _message.c_str();
    }

    InternalError::InternalError(const std::string & message) noexcept :
        Exception("Internal Error: " + message)
    {
    }

    UnknownObservableError::UnknownObservableError(const std::string & message) noexcept :
        Exception("Unknown Observable Error: " + message)
    {
    }

    GSLError::GSLError(const std::string & message) noexcept :
        Exception("GSL Error: " + message)
    {
    }

    ParsingError::ParsingError(const std::string & message) noexcept :
        Exception("Parsing Error: " + message)
    {
    }
} // namespace eos
