/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2015-2023 Danny van Dyk
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

#include <config.h>

#include <eos/statistics/log-posterior.hh>
#include <eos/utils/density-impl.hh>
#include <eos/utils/log.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <gsl/gsl_cdf.h>

#include <algorithm>

namespace eos
{
    struct RangeError :
        public Exception
    {
        RangeError(const std::string & message) throw () :
            Exception("Range Error: " + message)
        {
        }
    };

    LogPosterior::LogPosterior(const LogLikelihood & log_likelihood) :
        _log_likelihood(log_likelihood),
        _parameters(log_likelihood.parameters()),
        _informative_priors(0)
    {
    }

    LogPosterior::~LogPosterior()
    {
    }

    bool
    LogPosterior::add(const LogPriorPtr & prior, bool /*nuisance*/)
    {
        // clone has correct Parameters object selected
        LogPriorPtr prior_clone = prior->clone(_parameters);
        _informative_priors += 1 ? prior->informative() : 0;

        // extract parameters, and record their names to check for duplicates
        std::set<QualifiedName> prior_parameter_names;
        for (auto p = prior_clone->begin(), p_end = prior_clone->end() ; p != p_end ; ++p)
        {
            prior_parameter_names.insert(p->name());
        }

        // check if param exists already
        std::set<QualifiedName> intersection;
        std::set_intersection(_parameter_names.begin(), _parameter_names.end(),
                              prior_parameter_names.begin(), prior_parameter_names.end(),
                              std::inserter(intersection, intersection.begin()));

        if (intersection.size() > 0)
        {
            return false;
        }

        // if not, add to prior container and register parameter objects
        _priors.push_back(prior_clone);
        for (auto p = prior_clone->begin(), p_end = prior_clone->end() ; p != p_end ; ++p)
        {
            _varied_parameters.push_back(*p);
        }

        return true;
    }

    LogPosteriorPtr
    LogPosterior::clone() const
    {
        // clone log_likelihood
        LogLikelihood llh = _log_likelihood.clone();
        LogPosterior * result = new LogPosterior(llh);

        // add parameters via prior clones
        for (const auto & _prior : _priors)
        {
            result->add(_prior->clone(result->parameters()));
        }

        return LogPosteriorPtr(result);
    }

    double
    LogPosterior::evaluate() const
    {
        return log_posterior();
    }

    Parameters
    LogPosterior::parameters() const
    {
        return _parameters;
    }

    LogLikelihood
    LogPosterior::log_likelihood() const
    {
        return _log_likelihood;
    }

    double
    LogPosterior::log_posterior() const
    {
        return log_prior() + _log_likelihood();
    }

    double
    LogPosterior::log_prior() const
    {
        if (_priors.empty())
            throw InternalError("LogPosterior::log_prior(): prior is undefined");

        double result = 0.0;

        // all prior components are assumed independent,
        // thus the logs can be simply added up
        for (const auto & _prior : _priors)
        {
            result += (*_prior)();
        }

        return result;
    }

    template <>
    struct WrappedForwardIteratorTraits<LogPosterior::PriorIteratorTag>
    {
        using UnderlyingIterator = std::vector<LogPriorPtr>::const_iterator;
    };
    template class WrappedForwardIterator<LogPosterior::PriorIteratorTag, const LogPriorPtr>;

    LogPosterior::PriorIterator
    LogPosterior::begin_priors() const
    {
        return _priors.begin();
    }

    LogPosterior::PriorIterator
    LogPosterior::end_priors() const
    {
        return _priors.end();
    }

    unsigned
    LogPosterior::informative_priors() const
    {
        return _informative_priors;
    }

    Parameter
    LogPosterior::operator[] (const unsigned & index) const
    {
        return _varied_parameters[index];
    }

    const std::vector<Parameter> &
    LogPosterior::varied_parameters() const
    {
        return _varied_parameters;
    }
}
