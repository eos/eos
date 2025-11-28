/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
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

#include "eos/models/model.hh"
#include "eos/utils/exception.hh"

#include <boost/python.hpp>

#ifndef EOS_PYTHON__EOS_WRAPPERS_HH
#  define EOS_PYTHON__EOS_WRAPPERS_HH 1

namespace impl
{
    // raw constructor for class Kinematics
    boost::python::object Kinematics_ctor(boost::python::tuple args, boost::python::dict kwargs);

    // raw constructor for class Options
    boost::python::object Options_ctor(boost::python::tuple args, boost::python::dict kwargs);

    // converter for eos::Exception
    void translate_exception(const eos::Exception & e);

    // wrappers to avoid issues with virtual inheritance and overloading
    inline double
    m_b_pole_wrapper_noargs(const eos::Model & m)
    {
        return m.m_b_pole();
    }
} // namespace impl

#endif // EOS_PYTHON__EOS_WRAPPERS_HH
