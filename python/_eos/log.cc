/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2023-2026 Danny van Dyk
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

#include "python/_eos/gil.hh"

#include <memory>
#include <vector>

using namespace boost::python;

namespace impl
{
    namespace
    {
        // Holds a reference to a registered Python callable for as long as the interpreter is alive.
        struct LogCallback
        {
                object callback;
        };

        std::vector<std::shared_ptr<LogCallback>> &
        registered_callbacks()
        {
            static std::vector<std::shared_ptr<LogCallback>> result;

            return result;
        }

        void
        logging_callback(const std::shared_ptr<LogCallback> & c, const std::string & id, const eos::LogLevel & l, const std::string & m)
        {
            // a message emitted from a static destructor arrives once the interpreter has gone
            if (! Py_IsInitialized())
            {
                return;
            }

            // the Log notifies us from whichever thread emitted the message, holding no lock of its own
            ScopedGILAcquire gil;

            // release_python_log_callbacks() has dropped the callable
            if (c->callback.is_none())
            {
                return;
            }

            try
            {
                call<void>(c->callback.ptr(), id, l, m);
            }
            catch (const error_already_set &)
            {
                // the Log notifies us from a destructor, so nothing may be thrown from here
                PyErr_WriteUnraisable(c->callback.ptr());
            }
        }
    } // namespace

    void
    register_log_callback(PyObject * c)
    {
        auto callback = std::make_shared<LogCallback>(LogCallback{ object(handle<>(borrowed(c))) });

        registered_callbacks().push_back(callback);

        eos::Log::instance()->register_callback(std::bind(&logging_callback, callback, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
    }

    void
    release_python_log_callbacks()
    {
        for (auto & c : registered_callbacks())
        {
            c->callback = object();
        }

        registered_callbacks().clear();
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
