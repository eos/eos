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
#include <eos/utils/parameters.hh>
#include <eos/utils/scan_file.hh>
#include <eos/utils/stringify.hh>

#include <vector>

namespace eos
{
    class MarkovChain :
        public PrivateImplementationPattern<MarkovChain>
    {
        public:
            class History;
            class InfoAtPoint;
            class Stats;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor
             *
             * @param analysis The analysis for which we run.
             * @param seed     The initial seed for the RNG.
             */
            MarkovChain(const std::shared_ptr<Analysis> & analysis, unsigned long seed = 0);

            /// Destructor.
            ~MarkovChain();
            ///@}

            /// Remove existing history of this chain.
            void clear();

            /// Dump the history of this chain into a scan file's data set
            void dump_history(ScanFile::DataSet & data_set);

            /// Retrieve information regarding the current point.
            const InfoAtPoint & info_at_current() const;

            /// Retrive information regarding the most recently proposed point.
            const InfoAtPoint & info_at_proposal() const;

            /// Retrieve the number of iterations used in the last run
            const unsigned & iterations_last_run() const;

            /*!
             * Set the chain to continue its walk from the given point.
             *
             * @param point The point in parameter space.
             */
            void set_point(const std::vector<double> & point);

            /*!
             * Retrieve the scale used for proposals of a given parameter.
             *
             * @param index The index of the parameter.
             */
            double get_scale(const unsigned & index) const;

            /*!
             * Set the scale used for proposals of all parameters.
             *
             * @param scale The scale that shall be used for all parameter proposals.
             */
            void set_scale(const double & scale);

            /*!
             * Set the scale used for proposals of the given parameter.
             *
             * @param index The index of the parameter.
             * @param scale The scale that shall be used to this parameter's proposals.
             */
            void set_scale(const unsigned & index, const double & scale);

            /*!
             * Set whether the chain stores samples in runs to come.
             *
             * @param keep If true, the history of the runs to come will be kept.
             */
            void keep_history(bool keep);

            /// Retrieve the descriptions of all parameters that are explored by this chain.
            const std::vector<ParameterDescription> & parameter_descriptions() const;

            /// Check whether the most recently proposed move was accepted.
            bool proposal_accepted() const;

            /*!
             * Clear all statistics and counter.
             *
             * Do not change current position nor scale.
             * @param hard  If true, rease all statistics (use to clear prerun data before main run)
             * */
            void reset(bool hard = false);

            /*!
             * Perform a number of iterations.
             *
             * @param iterations The number of iterations that shall be performed.
             */
            void run(const unsigned & iterations);

            /// Retrieve statistical data that summarizes the evolution of the chain up to the current point.
            const Stats & statistics() const;
    };

    /*!
     * Holds the entire history of a run of a MarkovChain
     */
    struct MarkovChain::History
    {
        /// flag: if false => don't store numbers
        bool keep;

        /// sequence of points visited by Markov chain
        std::vector<std::vector<double>> points;

        /// at each point in the chain
        std::vector<double> log_likelihood;

        /// at each point in the chain
        std::vector<double> log_posterior;

        /// at each point in the chain
        std::vector<double> log_prior;
    };

    /*!
     * Summarize info at current position
     * in parameter space
     */
    struct MarkovChain::InfoAtPoint
    {
        ///position in parameter space
        std::vector<double> point;

        ///log likelihood at the point
        double log_likelihood;

        ///log prior at the point
        double log_prior;

        ///log posterior at the points
        double log_posterior;
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
        double iterations_total;

        /*!
         * The number of accepted proposals.
         *
         * @note accepted and rejected only add up to the number of iterations in the current run. The total
         * number of samples used to calculated the variance maybe different from the sum of former two.
         *
         * TODO right now we calculate efficiency in each run, and variance for entire chain. Why not calculate efficiency for entire sequence as well?
         */
        std::vector<unsigned> iterations_accepted;

        /*!
         * The number of iterations in which the proposed move of a parameter has been rejected.
         * Reset each time run() is called
         */
        std::vector<unsigned> iterations_rejected;

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

    std::ostream & operator<< (std::ostream & lhs, const MarkovChain::InfoAtPoint & rhs);
}

#endif
