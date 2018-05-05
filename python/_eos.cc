/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2016 Danny van Dyk
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

#include "config.h"

#include "eos/constraint.hh"
#include "eos/observable.hh"
#include "eos/signal-pdf.hh"
#include "eos/utils/kinematic.hh"
#include "eos/utils/model.hh"
#include "eos/utils/parameters.hh"
#include "eos/utils/options.hh"
#include "eos/utils/qualified-name.hh"
#include "eos/statistics/log-likelihood.hh"

#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

using namespace boost::python;
using namespace eos;

namespace impl
{
    // raw constructor for class Kinematics
    object
    Kinematics_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args = tuple(args.slice(1,_));

        self.attr("__init__")();

        list items = kwargs.items();
        for (unsigned i = 0 ; i < len(items) ; ++i)
        {
            object name = items[i][0];
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
        args = tuple(args.slice(1,_));

        self.attr("__init__")();

        list items = kwargs.items();
        for (unsigned i = 0 ; i < len(items) ; ++i)
        {
            object name = items[i][0];
            object value = items[i][1];
            self.attr("set")(name, value);
        }

        return object();
    }

    const char *
    version(void)
    {
        static const char version[] = PACKAGE_VERSION;

        return version;
    }
}

BOOST_PYTHON_MODULE(_eos)
{
    using namespace boost::python;
    using namespace eos;

    // {{{ eos/utils
    // qnp::Prefix
    class_<qnp::Prefix>("qnpPrefix", init<std::string>())
        .def("__repr__", &qnp::Prefix::str, return_value_policy<copy_const_reference>())
        .def("__str__", &qnp::Prefix::str, return_value_policy<copy_const_reference>())
        .def("__lt__", &qnp::Prefix::operator<)
        ;

    // qnp::Name
    class_<qnp::Name>("qnpName", init<std::string>())
        .def("__repr__", &qnp::Name::str, return_value_policy<copy_const_reference>())
        .def("__str__", &qnp::Name::str, return_value_policy<copy_const_reference>())
        .def("__lt__", &qnp::Name::operator<)
        ;

    // qnp::Suffix
    class_<qnp::Suffix>("qnpSuffix", init<std::string>())
        .def("__repr__", &qnp::Suffix::str, return_value_policy<copy_const_reference>())
        .def("__str__", &qnp::Suffix::str, return_value_policy<copy_const_reference>())
        .def("__lt__", &qnp::Suffix::operator<)
        ;

    // QualifiedName
    class_<QualifiedName>("QualifiedName", init<std::string>())
        .def("__repr__", &QualifiedName::full, return_value_policy<copy_const_reference>())
        .def("__str__", &QualifiedName::str, return_value_policy<copy_const_reference>())
        .def("__eq__", &QualifiedName::operator==)
        .def("__ne__", &QualifiedName::operator!=)
        .def("__lt__", &QualifiedName::operator<)
        .def("prefix_part", &QualifiedName::prefix_part, return_value_policy<copy_const_reference>())
        .def("name_part", &QualifiedName::name_part, return_value_policy<copy_const_reference>())
        .def("suffix_part", &QualifiedName::suffix_part, return_value_policy<copy_const_reference>())
        ;
    implicitly_convertible<std::string, QualifiedName>();

    // Parameters
    class_<Parameters>("_Parameters", no_init)
        .def("Defaults", &Parameters::Defaults)
        .staticmethod("Defaults")
        .def("__getitem__", (Parameter (Parameters::*)(const std::string &) const) &Parameters::operator[])
        .def("__iter__", range(&Parameters::begin, &Parameters::end))
        .def("set", &Parameters::set)
        .def("override_from_file", &Parameters::override_from_file)
        ;

    // Parameter
    class_<Parameter>("Parameter", no_init)
        .def(float_(self))
        .def("central", &Parameter::central, return_value_policy<copy_const_reference>())
        .def("max", &Parameter::max, return_value_policy<copy_const_reference>())
        .def("min", &Parameter::min, return_value_policy<copy_const_reference>())
        .def("name", &Parameter::name, return_value_policy<copy_const_reference>())
        .def("latex", &Parameter::latex, return_value_policy<copy_const_reference>())
        .def("set", &Parameter::set)
        .def("evaluate", &Parameter::evaluate)
        ;

    // Kinematics
    class_<Kinematics>("Kinematics", no_init)
        .def("__init__", raw_function(&impl::Kinematics_ctor))
        .def(init<>())
        .def("__getitem__", (KinematicVariable (Kinematics::*)(const std::string &) const) &Kinematics::operator[])
        .def("declare", &Kinematics::declare, return_value_policy<return_by_value>())
        .def("as_string", &Kinematics::as_string)
        ;

    // KinematicVariable
    class_<KinematicVariable>("KinematicVariable", no_init)
        .def(float_(self))
        .def("name", &KinematicVariable::name, return_value_policy<copy_const_reference>())
        .def("set", &KinematicVariable::set)
        .def("evaluate", &KinematicVariable::evaluate)
        ;

    // Options
    class_<Options>("Options", no_init)
        .def("__init__", raw_function(&impl::Options_ctor))
        .def(init<>())
        .def("set", &Options::set)
        .def("as_string", &Options::as_string)
        ;

    // Model
    register_ptr_to_python<std::shared_ptr<Model>>();
    class_<Model, boost::noncopyable>("Model", no_init)
        .def("make", &Model::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        // CKM component
        .def("ckm_cd", &Model::ckm_cd)
        .def("ckm_cs", &Model::ckm_cs)
        .def("ckm_cb", &Model::ckm_cb)
        .def("ckm_ud", &Model::ckm_ud)
        .def("ckm_us", &Model::ckm_us)
        .def("ckm_ub", &Model::ckm_ub)
        .def("ckm_td", &Model::ckm_td)
        .def("ckm_ts", &Model::ckm_ts)
        .def("ckm_tb", &Model::ckm_tb)
        // QCD component
        .def("m_t_msbar",  &Model::m_t_msbar)
        .def("m_t_pole",   &Model::m_t_pole)
        .def("m_b_kin",    &Model::m_b_kin)
        .def("m_b_msbar",  &Model::m_b_msbar)
        .def("m_b_pole",   &Model::m_b_pole)
        .def("m_c_kin",    &Model::m_c_kin)
        .def("m_c_msbar",  &Model::m_c_msbar)
        .def("m_c_pole",   &Model::m_c_pole)
        .def("m_s_msbar",  &Model::m_s_msbar)
        .def("m_ud_msbar", &Model::m_ud_msbar)
        ;
    // }}}

    // {{{ eos/statistics
    // LogLikelihoodBlock
    register_ptr_to_python<std::shared_ptr<LogLikelihoodBlock>>();
    class_<LogLikelihoodBlock, boost::noncopyable>("LogLikelihoodBlock", no_init)
        .def("as_string", &LogLikelihoodBlock::as_string)
        ;

    // Constraint
    class_<Constraint>("Constraint", no_init)
        .def("make", &Constraint::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        .def("name", &Constraint::name, return_value_policy<copy_const_reference>())
        .def("blocks", range(&Constraint::begin_blocks, &Constraint::end_blocks))
        .def("observables", range(&Constraint::begin_observables, &Constraint::end_observables))
        ;

    // }}}

    // {{{ eos/
    // Observable
    register_ptr_to_python<std::shared_ptr<Observable>>();
    class_<Observable, boost::noncopyable>("Observable", no_init)
        .def("make", &Observable::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        .def("evaluate", &Observable::evaluate)
        .def("name", &Observable::name, return_value_policy<copy_const_reference>())
        ;

    // SignalPDF
    register_ptr_to_python<std::shared_ptr<SignalPDF>>();
    class_<SignalPDF, boost::noncopyable>("SignalPDF", no_init)
        .def("make", &SignalPDF::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        .def("evaluate", &SignalPDF::evaluate)
        .def("name", &SignalPDF::name, return_value_policy<copy_const_reference>())
        ;
    // }}}

    // EOS version
    def("version", impl::version);
}
