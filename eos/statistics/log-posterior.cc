/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <gsl/gsl_cdf.h>

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
    LogPosterior::add(const LogPriorPtr & prior, bool nuisance)
    {
        // clone has correct Parameters object selected
        LogPriorPtr prior_clone = prior->clone(_parameters);
        _informative_priors += 1 ? prior->informative() : 0;

        // check if param exists already
        // read out parameters from prior
        for (auto d = prior->begin(), d_end = prior->end() ; d != d_end ; ++d)
        {
            auto result = _parameter_names.insert(d->parameter->name());
            if (! result.second)
                return false;

            d->nuisance = nuisance;
            _parameter_descriptions.push_back(*d);
        }

        // then add to prior container
        _priors.push_back(prior_clone);

        return true;
    }

    DensityPtr
    LogPosterior::clone() const
    {
        return DensityPtr(private_clone());
    }

    LogPosteriorPtr
    LogPosterior::old_clone() const
    {
        return LogPosteriorPtr(private_clone());
    }

    LogPosterior *
    LogPosterior::private_clone() const
    {
        // clone log_likelihood
        LogLikelihood llh = _log_likelihood.clone();
        LogPosterior * result = new LogPosterior(llh);

        // add parameters via prior clones
        for (const auto & _prior : _priors)
        {
            result->add(_prior->clone(result->parameters()));
        }

        // copy proper range for subspace sampling
        auto j = result->_parameter_descriptions.begin();
        for (auto i = _parameter_descriptions.begin(), i_end = _parameter_descriptions.end() ; i != i_end ; ++i, ++j)
        {
            j->min = i->min;
            j->max = i->max;
            j->nuisance = i->nuisance;
        }

        return result;
    }

    double
    LogPosterior::evaluate() const
    {
        return log_posterior();
    }

    Density::Iterator
    LogPosterior::begin() const
    {
        return Density::Iterator(_parameter_descriptions.cbegin());
    }

    Density::Iterator
    LogPosterior::end() const
    {
        return Density::Iterator(_parameter_descriptions.cend());
    }

    Parameters
    LogPosterior::parameters() const
    {
        return _parameters;
    }

    unsigned
    LogPosterior::index(const std::string & name) const
    {
        unsigned result = 0;

        for (auto d = _parameter_descriptions.cbegin(), d_end = _parameter_descriptions.cend() ; d != d_end ; ++d, ++result)
        {
            if (name == d->parameter->name())
                return result;
        }

        throw InternalError("Implementation<Analysis>::definition: no such parameter '" + name + "'");
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

    LogPriorPtr
    LogPosterior::log_prior(const std::string & name) const
    {
        LogPriorPtr prior;

        // loop over all descriptions of the prior pointers
        for (const auto & _prior : _priors)
        {
            for (auto i = _prior->begin(), i_end = _prior->end() ; i != i_end ; ++i)
            {
                if (i->parameter->name() == name)
                    prior = _prior;
            }
        }

        return prior;
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

    bool
    LogPosterior::nuisance(const std::string& par_name) const
    {
        unsigned index = this->index(par_name);

        if (index >= _parameter_descriptions.size())
        {
            return false;
        }
        else
        {
            return _parameter_descriptions[index].nuisance;
        }
    }

    unsigned
    LogPosterior::informative_priors() const
    {
        return _informative_priors;
    }

    MutablePtr
    LogPosterior::operator[] (const unsigned & index) const
    {
        return _parameter_descriptions[index].parameter;
    }

    const std::vector<ParameterDescription> &
    LogPosterior::parameter_descriptions() const
    {
        return _parameter_descriptions;
    }
}
