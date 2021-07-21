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

#include <set>
#include <vector>

namespace eos
{
    class LogPosterior :
        public Density
    {
        public:
            friend struct Implementation<LogPosterior>;

            struct OptimizationOptions;
            struct Output;

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

     // todo move optimization into separate class
     struct LogPosterior::OptimizationOptions
     {
             /// Keep the value of nuisance parameters with a flat prior fixed at the current value during optimization,
             /// to avoid flat directions that cause Migrad to fail.
             bool fix_flat_nuisance;

             /// Fraction of parameter range, in [0,1].
             /// Useful only for simplex method
             VerifiedRange<double> initial_step_size;

             /// If algorithm doesn't converge before, quit
             /// after maximum_iterations.
             unsigned maximum_iterations;

             /*!
              * If non-zero, perform MCMC iterations first,
              * before Minuit2 is invoked from the last point of the chain.
              *
              * @note This is only useful when call from MarkovChainSampler,
              *       further control of the chain is taken from the MarkovChainSampler::Config object.
              */
             bool mcmc_pre_run;

             /*!
              *  Once the algorithm has shrunk the probe
              *  simplex below this size, convergence is declared.
              *
              *  For minuit, it is just their tolerance parameter
              */
             VerifiedRange<double> tolerance;

             /*!
              * When comparing two modes found by minuit to decide whether they correspond
              * to the same mode, the splitting_tolerance decides how far
              * in relative units their distance may be.
              */
             VerifiedRange<double> splitting_tolerance;

             /// 0 - low, 1 - medium, 2 - high precision
             VerifiedRange<unsigned> strategy_level;

             static OptimizationOptions Defaults();

         private:
             OptimizationOptions();
     };

     struct LogPosterior::Output
     {
         using  DescriptionType = hdf5::Composite<hdf5::Scalar<const char *>, hdf5::Scalar<double>, hdf5::Scalar<double>,
                                                  hdf5::Scalar<int>, hdf5::Scalar<const char *>>;
         static DescriptionType description_type();
         static std::tuple<const char *, double, double, int, const char *> description_record();
     };

     /*!
      * Compute an initial guess of the proposal covariance matrix.
      *
      * The variance of each parameter is taken from the prior distribution
      * and scaled if desired for higher efficiency. Zero correlation is assumed a priori.
      *
      * @param log_posterior The log posterior supplying the prior information.
      * @param scale_reduction Value by which sqrt(variance) of parameters is divided.
      * @param scale_nuisance Decide whether only scan parameters or all parameters are scaled.
      * @return The covariance matrix in row major format.
      */
     std::vector<double> proposal_covariance(const LogPosterior & log_posterior,
                                             double scale_reduction=1,
                                             bool scale_nuisance=true);
}

#endif
