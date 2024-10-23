/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2023 Danny van Dyk
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

#include "python/_eos/external-log-likelihood-block.hh"

#include "eos/statistics/test-statistic-impl.hh"

using boost::python::extract;
using boost::python::object;

namespace eos
{
    ExternalLogLikelihoodBlock::ExternalLogLikelihoodBlock(const ObservableCache & cache, object factory) :
        _cache(cache),
        _factory(factory),
        _python_llh_block(_factory(_cache)),
        _evaluate(_python_llh_block.attr("evaluate")),
        _number_of_observations(extract<unsigned>(_python_llh_block.attr("number_of_observations")))
    {
        if (! PyCallable_Check(_evaluate.ptr()))
        {
            throw InternalError("ExternalLogLikelihoodBlock encountered a factory that does not yield a callable 'evaluate()' attribute");
        }
    }

    ExternalLogLikelihoodBlock::~ExternalLogLikelihoodBlock() = default;

    LogLikelihoodBlockPtr
    ExternalLogLikelihoodBlock::make(const ObservableCache & cache, object factory)
    {
        return LogLikelihoodBlockPtr(new ExternalLogLikelihoodBlock(cache, factory));
    }

    std::string
    ExternalLogLikelihoodBlock::as_string() const
    {
        return "ExternalLikelihoodBlock";
    }

    LogLikelihoodBlockPtr
    ExternalLogLikelihoodBlock::clone(ObservableCache cache) const
    {
        return LogLikelihoodBlockPtr(new ExternalLogLikelihoodBlock(cache, this->_factory));
    }

    double
    ExternalLogLikelihoodBlock::evaluate() const
    {
        return extract<double>(_evaluate());
    }

    unsigned int
    ExternalLogLikelihoodBlock::number_of_observations() const
    {
        return _number_of_observations;
    }

    double
    ExternalLogLikelihoodBlock::sample(gsl_rng *) const
    {
        throw InternalError("Not implemented");
        return 0.0;
    }

    double
    ExternalLogLikelihoodBlock::significance() const
    {
        throw InternalError("Not implemented");
        return 0.0;
    }

    TestStatistic
    ExternalLogLikelihoodBlock::primary_test_statistic() const
    {
        return test_statistics::Empty();
    }
} // namespace eos
