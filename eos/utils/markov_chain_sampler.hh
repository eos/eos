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

#ifndef EOS_GUARD_SRC_UTILS_MC_SAMPLER_HH
#define EOS_GUARD_SRC_UTILS_MC_SAMPLER_HH 1

#include <eos/utils/analysis.hh>
#include <eos/utils/markov_chain.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/proposal_functions.hh>
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
             * @param analysis The analysis for which we draw samples.
             * @param config   The configuration of the samples.
             */
            MarkovChainSampler(const Analysis & analysis, const MarkovChainSampler::Config & config);

            /// Destructor.
            ~MarkovChainSampler();

            /*!
             * Take the output from several independent preruns as input and store a
             * global local proposal function to the output file in directory '/global local'
             * In addition, the meta information about constraints, parameters and their priors
             * is stored in '/descriptions'.
             *
             * @param output_file_name The name of the HDF5 file containing the global local proposal function.
             *                         if empty, do not create the proposal nor any output.
             * @param input_files HDF5 files containing the history, descriptions, stats, and local propopal densities
             *                    as created by EOS.
             * @param config The configuration options influencing how GlobalLocal is created.
             * @param analysis Compare the order of parameters in files vs the one given in analysis.
             * @return The full histories of all chains found in the input files
             */
            static std::vector<HistoryPtr> build_global_local(const std::string & output_file_name,
                                    const std::vector<std::shared_ptr<hdf5::File>> input_files,
                                    const proposal_functions::GlobalLocal::Config & config,
                                    AnalysisPtr analysis = AnalysisPtr());

            /*!
             * Copy the settings such as proposal density
             * from the output of a (successful) prerun to prepare
             * for calling run(), where the main run is started
             * immediately.
             */
            void resume(const hdf5::File &);

            ///@}

            ///@name Sampling
            ///@{

            /*!
             * Use minuit to find local minima.
             *
             * Take as many starting points as there are chains
             *
             * @param options Pass options on to Minuit.
             */
            void massive_mode_finding(const Analysis::OptimizationOptions & options);

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
             * but not in posterior space. Thus, this R-value
             * impedes convergence declaration.
             */
            bool use_posterior_rvalue;

            /// Rescale multivariate proposal functions' covariance depending
            /// on the dimensionality of the parameter space.
            bool scale_automatic;

            /*!
             * Decide whether only scan parameters or all parameters' initial
             * variance should be rescaled in order to increase the efficiency.
             */
            bool scale_nuisance;

            /// Value by which sqrt(covariance) of scan parameters is reduced for initial guesses of the proposal functions.
            double scale_reduction;
            ///@}

            ///@name Prerun options
            ///@{

            /// Find the local maxima after the prerun, using the point
            /// with the highest posterior in each chain as a starting point
            bool find_modes;

            bool need_prerun;
            unsigned prerun_iterations_update;
            unsigned prerun_iterations_min;
            unsigned prerun_iterations_max;

            /// Which local proposal function is chosen
            std::string proposal;

            /// One proposal for the whole space or a decomposition into multivariate part
            /// and 1D, uncorrelated part.
            /// Map a fixed proposal in 1D to a parameter
            std::vector<std::string> block_proposal_parameters;

            /*!
             *  The number of degrees of freedom for a local proposal function
             *  of type MultivariateStudentT
             *  @note a value of one corresponds to a (multivariate) Cauchy function
             */
            VerifiedRange<double> student_t_degrees_of_freedom;

            /// Whether to store prerun samples.
            bool store_prerun;

            /// In case there is only one partition, this is ignored and number_of_chains is used.
            /// Else, this many chains are used per partition.
            unsigned prerun_chains_per_partition;
            std::vector<std::vector<std::tuple<std::string, double, double>>> partitions;
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

            /// Summarize all options pertaining to the GlobalLocal proposal function.
            std::shared_ptr<proposal_functions::GlobalLocal::Config> global_local_config;

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

            bool store_observables_and_proposals;
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
