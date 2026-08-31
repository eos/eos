/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2026 Danny van Dyk
 *
 * Based upon 'paludis/util/log.cc', which is
 *
 *   Copyright (c) 2006-2011 Ciaran McCreesh
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
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/log.hh>
#include <eos/utils/mutex.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <iostream>
#include <memory>
#include <set>
#include <time.h>
#include <vector>

namespace eos
{
    /* LogLevel */

    std::ostream &
    operator<< (std::ostream & lhs, const LogLevel & rhs)
    {
        do
        {
            switch (rhs)
            {
                case ll_silent: lhs << "silent"; continue;

                case ll_error: lhs << "error"; continue;

                case ll_warning: lhs << "warning"; continue;

                case ll_success: lhs << "success"; continue;

                case ll_completed: lhs << "completed"; continue;

                case ll_inprogress: lhs << "inprogress"; continue;

                case ll_informational: lhs << "informational"; continue;

                case ll_debug: lhs << "debug"; continue;

                case ll_last: break;
            }

            throw InternalError("LogLevel::operator<<: Bad value for log_level");
        }
        while (false);

        return lhs;
    }

    std::istream &
    operator>> (std::istream & lhs, LogLevel & rhs)
    {
        std::string word;
        lhs >> word;

        do
        {
            if ("silent" == word)
            {
                rhs = ll_silent;
                break;
            }
            else if ("error" == word)
            {
                rhs = ll_error;
                break;
            }
            else if ("warning" == word)
            {
                rhs = ll_warning;
                break;
            }
            else if ("success" == word)
            {
                rhs = ll_success;
                break;
            }
            else if ("completed" == word)
            {
                rhs = ll_completed;
                break;
            }
            else if ("inprogress" == word)
            {
                rhs = ll_inprogress;
                break;
            }
            else if ("informational" == word)
            {
                rhs = ll_informational;
                break;
            }
            else if ("debug" == word)
            {
                rhs = ll_debug;
                break;
            }

            throw InternalError("LogLevel::operator>>: Bad input in stream");
        }
        while (false);

        return lhs;
    }

    /* Log */

    template <> struct Implementation<Log>
    {
            using Callback = std::function<void(const std::string &, const LogLevel &, const std::string &)>;

            // held by value so that an emitting thread can keep reading the list after it has
            // released our mutex; register_callback replaces the list rather than extending it
            using Callbacks = std::shared_ptr<const std::vector<Callback>>;

            Mutex mutex;

            LogLevel log_level;

            std::ostream * stream;

            std::string program_name;

            Callbacks callbacks;

            std::set<std::string> one_time_messages;

            Implementation() :
                log_level(ll_error),
                stream(nullptr),
                callbacks(std::make_shared<const std::vector<Callback>>())
            {
            }

            /*
             * Emits the message to the stream, and returns the callbacks that are yet to be
             * notified, or an empty pointer if the message is filtered out.
             *
             * The caller notifies them once it has released our mutex. A callback runs code that we
             * know nothing about, and may well wait for a resource that another thread holds while
             * emitting a message of its own; holding any lock of ours across the callback would
             * deadlock that pair of threads.
             */
            Callbacks
            message(const std::string & id, const LogLevel & l, const std::string & m)
            {
                Lock ll(mutex);

                if (l > log_level)
                {
                    return {};
                }

                Callbacks result = callbacks;

                if (! stream)
                {
                    return result;
                }

                *stream << program_name << '@' << ::time(0) << ": ";

                do
                {
                    switch (l)
                    {
                        case ll_error: *stream << "[ERROR " << id << "] "; continue;

                        case ll_warning: *stream << "[WARNING " << id << "] "; continue;

                        case ll_informational: *stream << "[INFO " << id << "] "; continue;

                        case ll_success: *stream << "[SUCCESS " << id << "] "; continue;

                        case ll_completed: *stream << "[COMPLETED " << id << "] "; continue;

                        case ll_inprogress: *stream << "[INPROGRESS " << id << "] "; continue;

                        case ll_debug: *stream << "[DEBUG " << id << "] "; continue;

                        case ll_silent: throw InternalError("Implementation<Log>::message: ll_silent used for a message");

                        case ll_last: break;
                    }

                    throw InternalError("Implementation<Log>::message: Bad value for log_level");
                }
                while (false);

                *stream << m << std::endl;

                return result;
            }
    };

    template class InstantiationPolicy<Log, Singleton>;

    Log::Log() :
        PrivateImplementationPattern<Log>(new Implementation<Log>)
    {
    }

    Log::~Log() {}

    LogLevel
    Log::get_log_level() const
    {
        Lock l(_imp->mutex);

        return _imp->log_level;
    }

    void
    Log::set_log_level(const LogLevel & log_level)
    {
        Lock l(_imp->mutex);

        _imp->log_level = log_level;
    }

    void
    Log::set_log_stream(std::ostream * stream)
    {
        Lock l(_imp->mutex);

        _imp->stream = stream;
    }

    void
    Log::set_program_name(const std::string & program_name)
    {
        Lock l(_imp->mutex);

        _imp->program_name = program_name;
    }

    void
    Log::register_callback(const std::function<void(const std::string &, const LogLevel &, const std::string &)> & callback)
    {
        Lock l(_imp->mutex);

        auto callbacks = std::make_shared<std::vector<Implementation<Log>::Callback>>(*_imp->callbacks);
        callbacks->push_back(callback);

        _imp->callbacks = std::move(callbacks);
    }

    void
    Log::_message(const std::string & id, const LogLevel & l, const std::string & m)
    {
        const auto callbacks = _imp->message(id, l, m);

        if (! callbacks)
        {
            return;
        }

        for (const auto & c : *callbacks)
        {
            c(id, l, m);
        }
    }

    LogMessageHandler
    Log::message(const std::string & id, const LogLevel & log_level)
    {
        return LogMessageHandler(this, log_level, id);
    }

    Log::OneTimeMessage::OneTimeMessage(const std::string & id, const LogLevel & log_level, const std::string & message)
    {
        {
            Lock ll(Log::instance()->_imp->mutex);

            if (! Log::instance()->_imp->one_time_messages.insert(id).second)
            {
                return;
            }
        }

        Log::instance()->_message(id, log_level, message + " (Further messages of this type will be suppressed.)");
    }

    /* LogMessageHandler */
    LogMessageHandler::LogMessageHandler(Log * const log, const LogLevel & log_level, const std::string & id) :
        _log(log),
        _log_level(log_level),
        _id(id)
    {
    }

    LogMessageHandler::~LogMessageHandler()
    {
        if ((0 == std::uncaught_exceptions()) && (! _message.empty()))
        {
            _log->_message(_id, _log_level, _message);
        }
    }

    void
    LogMessageHandler::_append(const std::string & s)
    {
        _message.append(s);
    }
} // namespace eos
