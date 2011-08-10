/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2011 Danny van Dyk
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

#ifndef EOS_GUARD_SRC_UTILS_MC_SAMPLER_HH
#define EOS_GUARD_SRC_UTILS_MC_SAMPLER_HH 1

#include <src/utils/analysis.hh>
#include <src/utils/markov_chain.hh>
#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern.hh>
#include <src/utils/verify.hh>

namespace eos
{
    class MarkovChainSampler :
        public PrivateImplementationPattern<MarkovChainSampler>
    {
        public:
            class Config;
            class PreRunInfo;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param analysis The analysis for which we draw samples.
             * @param config   The configuration of the samples.
             */
            MarkovChainSampler(const Analysis & analysis, const MarkovChainSampler::Config & config);

            /// Destructor.
            ~MarkovChainSampler();

            /*!
             * Copy the settings such as proposal density scale and current points
             * from the output of a (successful) prerun
             * @param scan_file the HDF5 output of the prerun
             */
            void resume(const std::shared_ptr<ScanFile> & scan_file);

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
    struct MarkovChainSampler::Config
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
             * Independant runs with identical seeds will produce 
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
             *  but not in posterior space. Thus, this R-value
             *  impedes convergence declaration.
             */
            bool use_posterior_rvalue;

            /// Initial scale of the Cauchy variable used for proposal point generation.
            double scale_initial;

            /// Minimal scale of the Cauchy variable used for proposal point generation
            double scale_min;

            /// Maximal scale of the Cauchy variable used for proposal point generation
            double scale_max;
            ///@}

            ///@name Prerun options
            ///@{
            bool need_prerun;
            unsigned prerun_iterations_update;
            unsigned prerun_iterations_min;
            unsigned prerun_iterations_max;

            /// Whether to store prerun samples.
            bool store_prerun;
            ///@}

            ///@name Main run options
            ///@{
            /// Number of chunks of sampling.
            unsigned chunks;

            /// Number of iterations per chunk.
            unsigned chunk_size;

            /// Whether to store collected samples.
            bool store;
            ///@}

            ///@name Output options
            ///@{
            /*!
             * The HDF5 output file to store the markov chains.
             */
            std::shared_ptr<ScanFile> output_file;
            ///@}
    };

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
