/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2013 Frederik Beaujean
 * Copyright (c) 2011, 2013 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_STATISTICS_MC_SAMPLER_HH
#define EOS_GUARD_SRC_STATISTICS_MC_SAMPLER_HH 1

#include <eos/statistics/markov-chain.hh>
#include <eos/statistics/proposal-functions.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/verify.hh>

namespace eos
{
    class MarkovChainSampler :
        public PrivateImplementationPattern<MarkovChainSampler>
    {
        public:
            class Config;
            struct PreRunInfo;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param density  The density to sample from.
             * @param config   The configuration of the samples.
             */
            MarkovChainSampler(const DensityPtr & density, const MarkovChainSampler::Config & config);

            /// Destructor.
            ~MarkovChainSampler();

            /*!
             * Read the history of Markov chains stored in the input files under `base`.
             */
            static std::vector<HistoryPtr> read_chains(const std::vector<std::shared_ptr<hdf5::File>> & input_files,
                                                       std::string base = "/prerun");
            ///@}

            ///@name Sampling
            ///@{

            /// Retrieve information about the prerun performance
            PreRunInfo pre_run_info();

            /// Start the Markov chain sampling.
            void run();

            /// Retrieve the configuration from which this sampler was constructed.
            const MarkovChainSampler::Config & config();
            ///@}
    };

    /*!
     * Stores all configuration options for a MarkovChainSampler.
     */
    class MarkovChainSampler::Config
    {
        private:
            /// Constructor.
            Config();

        public:
            ///@name Basic Function
            ///@{
            /*!
             * Named constructor
             *
             * MarkovChainSampler settings with reasonably chosen default values.
             */
            static Config Default();

            /*!
             * Named constructor
             *
             * MarkovChainSampler settings with values optimized for quick chain
             * convergence and evaluation.
             *
             * @note The convergence is not very reliable. Use with care!
             * If in doubt, use MarkovChainSampler::Config::Default
             */
            static Config Quick();
            ///@}

            ///@name Basic options
            ///@{
            /// Number of Markov chains
            VerifiedRange<unsigned> number_of_chains;

            /*!
             * The seed that is used to initialize the random number generator.
             * Independent runs with identical seeds will produce identical results.
             */
            unsigned long seed;

            /*!
             * If true, use as many threads as there are cores available.
             * If false, use only one thread.
             */
            bool parallelize;
            ///@}

            ///@name Convergence options
            ///@{

            /// #accepted / #trials should be between 0.15 and 0.50
            VerifiedRange<double> min_efficiency;
            VerifiedRange<double> max_efficiency;

            /// rvalues close to 1 indicates multichain convergence
            VerifiedRange<double> rvalue_criterion_param;
            VerifiedRange<double> rvalue_criterion_posterior;

            /*!
             * Use Gelman/Rubin definition if true, and relaxed BAT definition
             * if false.
             */
            bool use_strict_rvalue_definition;

            /*!
             * often the chains have mixed in parameter space,
             * but not in posterior space. Thus, this R-value
             * impedes convergence declaration.
             */
            bool use_posterior_rvalue;

            /// Rescale multivariate proposal functions' covariance depending
            /// on the dimensionality of the parameter space.
            bool scale_automatic;

            ///@}

            ///@name Prerun options
            ///@{

            bool need_prerun;
            unsigned prerun_iterations_update;
            unsigned prerun_iterations_min;
            unsigned prerun_iterations_max;

            /// Which local proposal function is chosen
            std::string proposal;

            /// Initial covariance matrix for multivariate proposal
            std::vector<double> proposal_initial_covariance;

            /*!
             *  The number of degrees of freedom for a local proposal function
             *  of type MultivariateStudentT
             *  @note a value of one corresponds to a (multivariate) Cauchy function
             */
            VerifiedRange<double> student_t_degrees_of_freedom;

            /// Whether to store prerun samples.
            bool store_prerun;
            ///@}

            ///@name Main run options
            ///@{

            /// Proposal function is adapted during the first iterations.
            /// Default: 0, i.e. don't adapt
            unsigned adapt_iterations;

            /// Number of chunks of sampling.
            unsigned chunks;

            /// Number of iterations per chunk.
            unsigned chunk_size;

            /// Turn off main run, so only prerun is performed
            bool need_main_run;

            /// When computing R-values, one can skip the first (skip_initial)% of the iterations.
            /// This provides more robust results, as it discards an initial burn in where the chain
            /// may be very far off from likely regions
            VerifiedRange<double> skip_initial;

            /// Whether to store collected samples.
            bool store;
            ///@}

            ///@name Output options
            ///@{
            /*!
             * The HDF5 output file to store the markov chains.
             */
            std::string output_file;
            ///@}
    };

    std::ostream & operator<<(std::ostream&, const MarkovChainSampler::Config & config);

    /*!
     * Holds convergence information of the prerun.
     */
    struct MarkovChainSampler::PreRunInfo
    {
        /// Convergence status after performing the prerun.
        bool converged;

        /// The number of iterations that were performed in the prerun.
        unsigned iterations;

        /// The number of iterations after which convergence was declared.
        unsigned iterations_at_convergence;

        /// R-value of posterior
        double rvalue_posterior;

        /// R-values of individual parameters
        std::vector<double> rvalue_parameters;
    };
}

#endif
