/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2016-2023 Danny van Dyk
 * Copyright (c) 2021-2023 Philip Lüghausen
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

#include "python/_eos/wrappers.hh"

using boost::python::_;
using boost::python::dict;
using boost::python::len;
using boost::python::list;
using boost::python::object;
using boost::python::tuple;

namespace impl
{
    // raw constructor for class Kinematics
    object
    Kinematics_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args        = tuple(args.slice(1, _));

        self.attr("__init__")();

        if (1 < len(args))
        {
            PyErr_SetString(PyExc_TypeError, "eos.Kinematics expects exactly one argument, or keyword arguments, but not both");
            return object();
        }

        dict kinematics;

        if (1 == len(args))
        {
            kinematics = dict(args[0]);
            args       = tuple(args.slice(1, _));
        }
        else
        {
            kinematics = kwargs;
        }

        list items = kinematics.items();
        for (unsigned i = 0; i < len(items); ++i)
        {
            object name  = items[i][0];
            object value = items[i][1];
            self.attr("declare")(name, value);
        }

        return object();
    }

    // raw constructor for class Options
    object
    Options_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args        = tuple(args.slice(1, _));

        self.attr("__init__")();

        if (1 < len(args))
        {
            PyErr_SetString(PyExc_TypeError, "eos.Options expects exactly one argument, or keyword arguments, but not both");
            return object();
        }

        dict options;

        if (1 == len(args))
        {
            options = dict(args[0]);
            args    = tuple(args.slice(1, _));
        }
        else
        {
            options = kwargs;
        }

        list items = options.items();
        for (unsigned i = 0; i < len(items); ++i)
        {
            object name  = items[i][0];
            object value = items[i][1];
            self.attr("declare")(name, value);
        }

        return object();
    }

    // converter for eos::Exception
    void
    translate_exception(const eos::Exception & e)
    {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    }
} // namespace impl
