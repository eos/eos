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

#include <src/utils/analysis.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

namespace eos
{
   template<>
   struct Implementation<Analysis>
   {
        LogLikelihoodPtr log_likelihood;

        Parameters parameters;

        // prior in N dimensions can decouple
        // at most into N 1D priors
        std::vector<LogPriorPtr> priors;

        // Parameter, minimum, maximum, nuisance
        std::vector<ParameterDescription> parameter_descriptions;

        // names of all parameters. prevent using a parameter twice
        std::set<std::string> parameter_names;

        Implementation(const LogLikelihoodPtr & log_likelihood) :
            log_likelihood(log_likelihood),
            parameters(log_likelihood->parameters())
        {
        }

        bool add_parameter(const LogPriorPtr & prior, bool nuisance)
        {
            // clone has correct Parameters object selected
            LogPriorPtr prior_clone = prior->clone(parameters);

            // check if param exists already
            // read out parameters from prior
            for (auto d = prior->begin(), d_end = prior->end() ; d != d_end ; ++d)
            {
                auto result = parameter_names.insert(d->parameter.name());
                if (! result.second)
                    return false;

                d->nuisance = nuisance;
                parameter_descriptions.push_back(*d);
            }

            // then add to prior container
            priors.push_back(prior_clone);

            return true;
        }

        AnalysisPtr clone() const
        {
            // clone log_likelihood
            LogLikelihoodPtr llh = LogLikelihoodPtr(new LogLikelihood(log_likelihood->clone()));
            AnalysisPtr result = std::make_shared<Analysis>(llh);

            // add parameters via prior clones
            for (auto i = priors.cbegin(), i_end = priors.cend(); i != i_end; ++i)
            {
                result->add((*i)->clone(result->_imp->parameters));
            }

            return result;
        }

        /*
         * Find index of definition of parameter
         * @param name
         * @return index if found, _parameter_descriptions.size() if not found
         */
        unsigned index(const std::string & name) const
        {
            unsigned result = 0;

            for (auto d = parameter_descriptions.cbegin(), d_end = parameter_descriptions.cend() ; d != d_end ; ++d, ++result)
            {
                if (name == d->parameter.name())
                    return result;
            }

            throw InternalError("Implementation<Analysis>::definition: no such parameter '" + name + "'");
        }

        bool nuisance(const std::string & name) const
        {
            unsigned index = this->index(name);

            if (index >= parameter_descriptions.size())
            {
                return false;
            }
            else
            {
                return parameter_descriptions[index].nuisance;
            }
        }

        double log_prior()
        {
            if (priors.empty())
                throw InternalError("Analysis::log_prior(): prior is undefined");

            double result = 0.0;

            // all prior components are assumed independent,
            // thus the logs can be simply added up
            for (auto p = priors.cbegin(), p_end = priors.cend() ; p != p_end; ++p)
            {
                result += (**p)();
            }

            return result;
        }
    };

    Analysis::Analysis(const LogLikelihoodPtr & log_likelihood) :
        PrivateImplementationPattern<Analysis>(new Implementation<Analysis>(log_likelihood))
    {
    }

    Analysis::~Analysis()
    {
    }

    bool
    Analysis::add(const LogPriorPtr & prior, bool nuisance)
    {
        return _imp->add_parameter(prior, nuisance);
    }

    AnalysisPtr
    Analysis::clone() const
    {
        return _imp->clone();
    }

    Parameters
    Analysis::parameters() const
    {
        return _imp->parameters;
    }

    LogLikelihoodPtr &
    Analysis::log_likelihood()
    {
        return _imp->log_likelihood;
    }

    double
    Analysis::log_posterior()
    {
        return _imp->log_prior() + (*_imp->log_likelihood)();
    }

    double
    Analysis::log_prior()
    {
        return _imp->log_prior();
    }

    bool
    Analysis::nuisance(const std::string& par_name) const
    {
        return _imp->nuisance(par_name);
    }

    Parameter
    Analysis::operator[] (const unsigned & index)
    {
        return _imp->parameter_descriptions[index].parameter;
    }

    const std::vector<ParameterDescription> &
    Analysis::parameter_descriptions() const
    {
        return _imp->parameter_descriptions;
    }
}
