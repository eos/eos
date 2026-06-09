/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2023-2025 Danny van Dyk
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

#include "python/_eos/log.hh"

#include <iostream>

using namespace boost::python;

namespace impl
{
    struct LoggingCallbackPayload
    {
            PyObject *          callback;
            const std::string   id;
            const eos::LogLevel log_level;
            const std::string   message;
    };

    int
    thread_safe_logging_callback(void * p)
    {
        // this function is only supposed to be called from the main thread
        // via Py_AddPendingCall, so we can safely call the Python API here
        const LoggingCallbackPayload * const payload = static_cast<const LoggingCallbackPayload *>(p);

        call<void>(payload->callback, payload->id, payload->log_level, payload->message);

        delete payload;

        return 0;
    }

    void
    logging_callback(PyObject * c, const std::string & id, const eos::LogLevel & l, const std::string & m)
    {
        LoggingCallbackPayload * const payload = new LoggingCallbackPayload{ c, id, l, m };

        int result = Py_AddPendingCall(&thread_safe_logging_callback, payload);

        if (result != 0)
        {
            std::cerr << "eos::Log: Warning: Could not schedule log callback via Py_AddPendingCall, error code " << result << "." << std::endl;

            delete payload;
        }
    }

    void
    register_log_callback(PyObject * c)
    {
        eos::Log::instance()->register_callback(std::bind(&logging_callback, c, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
    }

    void
    set_native_log_level(const eos::LogLevel & log_level)
    {
        eos::Log::instance()->set_log_level(log_level);
    }

    // for testing purpose only
    void
    emit_native_log(const std::string & id, const eos::LogLevel & log_level, const std::string & m)
    {
        eos::Log::instance()->message(id, log_level) << m;
    }
} // namespace impl
