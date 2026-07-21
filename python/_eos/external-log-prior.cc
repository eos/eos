/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#include "python/_eos/external-log-prior.hh"

using boost::python::extract;
using boost::python::object;
using boost::python::stl_input_iterator;

namespace eos
{
    ExternalLogPrior::ExternalLogPrior(const Parameters & parameters, object factory) :
        LogPrior(parameters),
        _factory(factory),
        _python_prior(_factory(parameters)),
        _evaluate(_python_prior.attr("evaluate")),
        _sample(_python_prior.attr("sample")),
        _compute_cdf(_python_prior.attr("compute_cdf")),
        _informative(extract<bool>(_python_prior.attr("informative")))
    {
        if (! PyCallable_Check(_evaluate.ptr()))
        {
            throw InternalError("ExternalLogPrior encountered a factory that does not yield a callable 'evaluate()' attribute");
        }

        if (! PyCallable_Check(_sample.ptr()))
        {
            throw InternalError("ExternalLogPrior encountered a factory that does not yield a callable 'sample()' attribute");
        }

        if (! PyCallable_Check(_compute_cdf.ptr()))
        {
            throw InternalError("ExternalLogPrior encountered a factory that does not yield a callable 'compute_cdf()' attribute");
        }

        // resolve the names of the varied parameters into the base class' container
        object names = _python_prior.attr("varied_parameters");
        for (stl_input_iterator<std::string> n(names), end; n != end; ++n)
        {
            _varied_parameters.push_back(_parameters[*n]);
        }
    }

    ExternalLogPrior::~ExternalLogPrior() = default;

    LogPriorPtr
    ExternalLogPrior::make(const Parameters & parameters, object factory)
    {
        return LogPriorPtr(new ExternalLogPrior(parameters, factory));
    }

    std::string
    ExternalLogPrior::as_string() const
    {
        return "ExternalLogPrior";
    }

    LogPriorPtr
    ExternalLogPrior::clone(const Parameters & parameters) const
    {
        return LogPriorPtr(new ExternalLogPrior(parameters, this->_factory));
    }

    double
    ExternalLogPrior::operator() () const
    {
        return extract<double>(_evaluate());
    }

    void
    ExternalLogPrior::sample()
    {
        _sample();
    }

    void
    ExternalLogPrior::compute_cdf()
    {
        _compute_cdf();
    }

    bool
    ExternalLogPrior::informative() const
    {
        return _informative;
    }
} // namespace eos
