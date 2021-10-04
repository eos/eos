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

#ifndef EOS_GUARD_SRC_STATISTICS_LOG_POSTERIOR_HH
#define EOS_GUARD_SRC_STATISTICS_LOG_POSTERIOR_HH 1

#include <config.h>

#include <eos/statistics/log-likelihood.hh>
#include <eos/statistics/log-posterior-fwd.hh>
#include <eos/statistics/log-prior.hh>
#include <eos/utils/density.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/verify.hh>
#include <eos/utils/wrapped_forward_iterator.hh>

#include <set>
#include <vector>

namespace eos
{
    class LogPosterior :
        public Density
    {
        public:
            friend struct Implementation<LogPosterior>;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * Extracts parameters, observables from LogLikelihood.
             * The default prior (flat) is assumed for all parameters.
             *
             * @param log_likelihood  The LogLikelihood functor which shall be analysed.
             *
             * @note LogPosterior assumes ownership of log_likelihood
             */
            LogPosterior(const LogLikelihood & log_likelihood);

            /// Destructor.
            virtual ~LogPosterior();

            /// Clone this LogPosterior
            // todo remove
            LogPosteriorPtr old_clone() const;

            virtual DensityPtr clone() const;

            virtual double evaluate() const;

            virtual Iterator begin() const;
            virtual Iterator end() const;
            ///@}

            ///@name Accessors
            ///@{
            // todo remove
            /// Retrieve a set of all parameters, including ranges
            const std::vector<ParameterDescription> & parameter_descriptions() const;

            // todo remove as functionality available from begin(), end()
            /*!
             * Retrieve a parameter by index.
             *
             * @param index The index of the parameter.
             */
            MutablePtr operator[] (const unsigned & index) const;

            // todo remove
            /// Retrieve our associated Parameters object
            Parameters parameters() const;

            /*!
             * Add one or more parameters and associated prior density
             *
             * @param prior    The logarithmic prior density.
             * @param nuisance False for a parameter of interest
             */
            bool add(const LogPriorPtr & prior, bool nuisance = false);

            /// Retrieve the overall Log(likelihood).
            LogLikelihood log_likelihood() const;

            /// Retrieve the overall Log(prior).
            double log_prior() const;

            /*!
             * Find the prior for a given parameter
             */
            LogPriorPtr log_prior(const std::string & name) const;

            /// Retrieve the overall Log(posterior)
            /// Incorporate normalization constant, the evidence here in getter if available.
            double log_posterior() const;

            /*!
             * Add forward iterator and corresponding helper functions
             */
            struct PriorIteratorTag;
            using PriorIterator = WrappedForwardIterator<PriorIteratorTag, const LogPriorPtr>;

            PriorIterator begin_priors() const;
            PriorIterator end_priors() const;

            /*!
             * Check if a given parameter is a nuisance parameter for this LogPosterior.
             *
             * @param name The name of the parameter we are interested in.
             */
            bool nuisance(const std::string & name) const;

            /*!
             * Returns the number of informative priors.
             */
            unsigned informative_priors() const;
            ///@}

        private:
            /*!
             * Find index of definition of parameter
             * @param name
             * @return index if found, _parameter_descriptions.size() if not found
             */
            unsigned index(const std::string & name) const;

            LogPosterior *
            private_clone() const;

            LogLikelihood _log_likelihood;

            Parameters _parameters;

            /// prior in N dimensions can decouple
            /// at most into N 1D priors
            std::vector<LogPriorPtr> _priors;

            unsigned _informative_priors;

            /// Parameter, minimum, maximum, nuisance
            std::vector<ParameterDescription> _parameter_descriptions;

            /// names of all parameters. prevent using a parameter twice
            std::set<std::string> _parameter_names;
    };

    extern template class WrappedForwardIterator<LogPosterior::PriorIteratorTag, const LogPriorPtr>;
}

#endif
