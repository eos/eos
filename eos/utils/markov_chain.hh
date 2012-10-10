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

#ifndef EOS_GUARD_SRC_UTILS_MARKOV_CHAIN_HH
#define EOS_GUARD_SRC_UTILS_MARKOV_CHAIN_HH 1

#include <eos/utils/analysis-fwd.hh>
#include <eos/utils/hdf5-fwd.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/stringify.hh>

#include <vector>

#include <gsl/gsl_rng.h>

namespace eos
{
    class MarkovChain :
        public PrivateImplementationPattern<MarkovChain>
    {
        public:
            class History;
            class State;
            class ProposalFunction;
            class HyperParameter;
            class Stats;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor
             *
             * @param analysis The analysis for which we run.
             * @param seed     The initial seed for the RNG.
             */
            MarkovChain(const Analysis & analysis, unsigned long seed, const std::shared_ptr<MarkovChain::ProposalFunction> & proposal_function);

            /// Destructor.
            ~MarkovChain();
            ///@}

            /// Remove existing history of this chain.
            void clear();

            /// Retrieve information regarding the current state.
            const State & current_state() const;

            // todo document
            void dump_description(hdf5::File & file, const std::string & data_set) const;

            /*!
             * Dump a part of the most recent history in the HDF5 file
             * under the given group name.
             *
             * @param file
             * @param data_set_name All output is stored below this directory.
             * @param last_iterations Dump only this many iterations.
             */
            void dump_history(hdf5::File & file, const std::string & data_set_name, const unsigned & last_iterations) const;

            void dump_proposal(hdf5::File & file, const std::string & data_set_name) const;

            MarkovChain::HyperParameter & hyper_parameter(bool current = true) const;

            /// Retrieve the number of iterations used in the last run
            const unsigned & iterations_last_run() const;

            /// Retrieve the chain's detailed history.
            const History & history() const;

            /// Set analysis.
            void set_analysis(const Analysis & analysis);

            /*!
             * Set whether the chain stores samples in runs to come.
             *
             * @param keep_samples If true, the history of the runs to come will be kept.
             * @param keep_observables_and_proposals For debugging purposes, store the proposed points
             *                                       and the observable values.
             */
            void keep_history(bool keep_samples, bool keep_observables_and_proposals = false);

            /*!
             * Compute estimate of normalized posterior density of this chain according to method
             * described in Chib, Jeliazkov (2001), Eq. (9)
             * @param <numerator, denominator> of [CJ2001], Eq. (9)
             * @param point where density is sought
             * @param posterior_evaluations the number of additional posterior evaluations
             *                              needed to estimate the normalized posterior density,
             *                              see denominator of [CJ2001], Eq. (9)
             */
            void normalized_density(std::tuple<double, double> & result, const std::vector<double> & point, const unsigned & posterior_evaluations) const;

            /// Retrieve the descriptions of all parameters that are explored by this chain.
            const std::vector<ParameterDescription> & parameter_descriptions() const;

            /// Check whether the most recently proposed move was accepted.
            bool proposal_accepted() const;

            /// Access to the proposal function.
            std::shared_ptr<MarkovChain::ProposalFunction> proposal_function() const;

            /// Set the proposal function.
            void proposal_function(const std::shared_ptr<MarkovChain::ProposalFunction> &);

            /// Retrieve information regarding the most recently proposed state.
            const State & proposed_state() const;

            /*!
             * Clear all statistics and counter.
             *
             * Do not change current position nor scale.
             * @param hard  If true, rease all statistics (use to clear prerun data before main run)
             * */
            void reset(bool hard = false);

            /*!
             * Read part of the output of a chain's prerun from hdf5 file.
             *
             * @param file
             * @param data_base_name The directory in the file under which the data is parsed.
             * @param history The parsed chain's history is returned by reference, assumed to be empty initially.
             * @param proposal The local proposal density is recreated and returned by reference.
             * @param proposal_type The type of the local proposal function used.
             * @param stats Only the mode of the posterior is restored and returned by reference.
             */
            static void read_data(hdf5::File & file, const std::string & data_base_name,
                                  MarkovChain::History & history,
                                  std::shared_ptr<MarkovChain::ProposalFunction> & proposal,
                                  std::string & proposal_type,
                                  MarkovChain::Stats & stats);

            /*!
             * Read the description part of chain's prerun from hdf5 file.
             *
             * @param file
             * @param data_base_name The directory in the file under which the data is parsed.
             * @param descriptions All parameter ranges etc. Beware, the association to the underlying Parameters object is independent.
             * @param priors The string representation of a prior distribution.
             * @param constraints The string representation of an individual constraint.
             * @param hash The EOS version used to create the file.
             */
            static void read_descriptions(hdf5::File & file, const std::string & data_base_name,
                                          std::vector<ParameterDescription>& descriptions,
                                          std::vector<std::string> & priors,
                                          std::vector<std::string> & constraints,
                                          std::string & hash);

            /*!
             * Perform a number of iterations.
             *
             * @param iterations The number of iterations that shall be performed.
             */
            void run(const unsigned & iterations);

            /*! Set the stats at mode to a point found outside of the chain.
             * Triggers writing to the HDF5 file as another row in the stats section.
             *
             * @param point
             * @param posterior
             */
            void set_mode(hdf5::File & file, const std::string & data_base_name,
                          const std::vector<double> & point, const double & posterior);

            /*!
             * Set the chain to continue its walk from the given point.
             *
             * @param point The point in parameter space.
             */
            void set_point(const std::vector<double> & point, const MarkovChain::HyperParameter & hyper_parameter);

            /// Retrieve statistical data that summarizes the evolution of the chain up to the current point.
            const Stats & statistics() const;
    };

    struct MarkovChain::HyperParameter
    {
            /// Indicating the current component.
            unsigned component;
    };

    /*!
     * Summarize info at current position
     * in parameter space
     */
    struct MarkovChain::State
    {
        typedef std::vector<State>::const_iterator Iterator;

        /// position in parameter space
        std::vector<double> point;

        /// log likelihood at the point
        double log_likelihood;

        /// log prior at the point
        double log_prior;

        /// log posterior at the points
        double log_posterior;

        /// Contains possibly multidimensional hyperparameter information.
        MarkovChain::HyperParameter hyper_parameter;

        State()
        {
        }

        State(const State & other)
        {
            log_likelihood = other.log_likelihood;
            log_prior = other.log_prior;
            log_posterior = other.log_posterior;
            point.resize(other.point.size());
            std::copy(other.point.cbegin(), other.point.cend(), point.begin());
            hyper_parameter = other.hyper_parameter;
        }
    };

    /*!
     * Holds statistical information of a run of a MarkovChain
     */
    struct MarkovChain::Stats
    {
        /*!
         * The total number of iterations per parameter that were collected before the
         * current run, i.e. since the last reset of statistics.
         */
        unsigned iterations_total;

        /*!
         * The number of accepted proposals.
         *
         * @note accepted and rejected only add up to the number of iterations in the current run. The total
         * number of samples used to calculated the variance maybe different from the sum of former two.
         */
        unsigned iterations_accepted;

        unsigned iterations_invalid;

        /*!
         * The number of iterations in which the proposed move of a parameter has been rejected.
         * Reset each time run() is called
         */
        unsigned iterations_rejected;

        /// Maximum value of the posterior
        double mode_of_posterior;

        /// Parameter values at the maximum of the posterior
        std::vector<double> parameters_at_mode;

        /// Sample mean of parameter values
        std::vector<double> mean_of_parameters;

        /// Sample mean of posterior
        double mean_of_posterior;

        /// Sample variance of parameter values
        std::vector<double> variance_of_parameters;

        /// Sample variance of [log] posterior values
        double variance_of_posterior;
    };

    typedef std::shared_ptr<MarkovChain::History> HistoryPtr;

    /*!
     * Holds the entire history of a run of a MarkovChain
     */
    struct MarkovChain::History
    {
        private:
            static bool cmp(const MarkovChain::State & a, const MarkovChain::State & b);

        public:

            /// flag: if false => don't store numbers
            bool keep;

            /// All states.
            std::vector<MarkovChain::State> states;

            /*!
             * Return state with highest posterior probability in selected range
             */
            const MarkovChain::State & local_mode(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end) const;

            /*!
             * Compute mean and variance of the states' parameters between begin and end
             * using Welford's method for all parameters.
             * Results are stored in the vectors.
             *
             * For details, check out http://www.johndcook.com/standard_deviation.html .
             *
             */
            // todo: make static
            void mean_and_variance(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                                   std::vector<double> & mean, std::vector<double> & variance) const;

            /*!
             * Compute the mean vector and the sample covariance in the given range of the chain.
             *
             * @param begin
             * @param end
             * @param mean
             * @param variance
             */
            void mean_and_covariance(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                                   std::vector<double> & mean, std::vector<double> & variance) const;
    };

    typedef std::shared_ptr<MarkovChain::ProposalFunction> ProposalFunctionPtr;

    /*!
     * Interface to proposal functions for a MarkovChain.
     */
    struct MarkovChain::ProposalFunction
    {
        /// (virtual) Destructor.
        virtual ~ProposalFunction() = 0;

        /*!
         *  Adapt the proposal function to the chain's current state and history.
         *  @note adapt() always uses the full history passed as an argument. If only
         *  a subset of an existing history is to be used, the caller is responsible
         *  for removing the unneeded parts from history before calling adapt().
         * @param history
         */
        virtual void adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                           const double & efficiency, const double & efficiency_min, const double & efficiency_max) = 0;

        /// Create independent copy.
        virtual ProposalFunctionPtr clone() const = 0;

        /// Store state in the file under the given base name
        virtual void dump_state(hdf5::File & file, const std::string & data_set_base_name) const = 0;

        /// Evaluate the density to propose x given y.
        virtual double evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const = 0;

        /// Obtain from the density a proposal x given y
        virtual void propose(MarkovChain::State & x, const MarkovChain::State & y, gsl_rng * rng) const = 0;
    };

    std::ostream & operator<< (std::ostream & lhs, const MarkovChain::State & rhs);
}

#endif
