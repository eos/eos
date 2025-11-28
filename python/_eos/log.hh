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

#include "eos/utils/log.hh"

#include <boost/python.hpp>

#ifndef EOS_PYTHON__EOS_LOG_HH
#  define EOS_PYTHON__EOS_LOG_HH 1

using namespace boost::python;

namespace impl
{
    void logging_callback(PyObject * c, const std::string & id, const eos::LogLevel & l, const std::string & m);

    void register_log_callback(PyObject * c);

    void set_native_log_level(const eos::LogLevel & log_level);

    // for testing purposes only
    void emit_native_log(const std::string & id, const eos::LogLevel & log_level, const std::string & m);
} // namespace impl

#endif // EOS_PYTHON__EOS_LOG_HH
