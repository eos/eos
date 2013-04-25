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

#ifndef EOS_GUARD_SRC_UTILS_PMC_SAMPLER_HH
#define EOS_GUARD_SRC_UTILS_PMC_SAMPLER_HH 1

#include <eos/utils/analysis.hh>
#include <eos/utils/hdf5-fwd.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class PopulationMonteCarloSampler :
        public PrivateImplementationPattern<PopulationMonteCarloSampler>
    {
        public:
            class Config;
            class Status;
            class Output;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param analysis The analysis for which we draw samples.
             * @param config   The configuration of the samples.
             */
            PopulationMonteCarloSampler(const Analysis & analysis, const PopulationMonteCarloSampler::Config & config);

            /*!
             * Initialize the PMC from a HDF5 file.
             *
             * Several kinds of input files are (automatically) recognized:
             *
             * a) GlobalLocal proposal function, stored in file.
             * b) MCMC prerun, used for hierarchical clustering.
             * c) Serialized status of a PMC proposal function from a previous PMC run.
             *
             * All options regarding construction of the PMC in config are ignored.
             */
            PopulationMonteCarloSampler(const Analysis & analysis, const hdf5::File & file,
                                        const PopulationMonteCarloSampler::Config & config, const bool & update = false);

            /// Destructor.
            ~PopulationMonteCarloSampler();

            ///@}

            ///@name Sampling
            ///@{

            /*!
             * Calculate the posterior and importance weights for the range of
             * parameter samples given by min and max indices from the sample file.
             */
            void calculate_weights(const std::string & sample_file, const unsigned & min_index, const unsigned & max_index) const;

            /// Retrieve the configuration from which this sampler was constructed.
            const PopulationMonteCarloSampler::Config & config() const;

            /*!
             * Draw parameter samples from the proposal density and store them,
             * along with the full status of all components.
             *
             * @param output_file Name of the HDF5 output file
             */
            void draw_samples();

            /*!
             * Read in a slice of samples from a previous PMC dump.
             *
             * @param sample_file Name of HDF5 file containing the samples.
             * @param base Directory name within HDF5 file.
             * @param min First element to parse.
             * @param max Index of one-past last element to parse.
             * @param samples Upon return, contains all samples stored. When passing in, there must be at least one sample in there to convey the parameter dimension
             */
            static void read_samples(const std::string & sample_file, const std::string & base,
                                     const unsigned & min, const unsigned & max,
                                     std::vector<std::vector<double>> & samples);

            /// Start the Markov chain sampling.
            void run();

            /// Retrieve the current status.
            const PopulationMonteCarloSampler::Status & status() const;

            /*!
             * Set the current status.
             *
             * @param new_status The new status.
             * @param check_convergence If true, compute convergence based on the new status and return true if converged. This option is for testing only!
             * @return True
             */
            bool status(const PopulationMonteCarloSampler::Status & new_status, bool check_convergence = false);

            ///@}
    };

    /*!
     * Stores all configuration options for a PopulationMonteCarloSampler.
     */
    struct PopulationMonteCarloSampler::Config
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
             * PopulationMonteCarloSampler settings with reasonably chosen default values.
             */
            static Config Default();

            /*!
             * Named constructor
             *
             * PopulationMonteCarloSampler settings with values optimized for quick chain
             * convergence and evaluation.
             *
             * @note The convergence is not very reliable. Use with care!
             * If in doubt, use PopulationMonteCarloSampler::Config::Default
             */
            static Config Quick();
            ///@}

            ///@name Basic options
            ///@{
            /*!
             * The seed that is used to initialize the random number generator.
             * Independent runs with identical seeds will produce
             */
            unsigned long seed;

            /*!
             * If true, use as many threads as there are cores available.
             * If false, use only one thread.
             */
            bool parallelize;


            /*!
             * How many workers to use. If parallelize is true,
             * this is the number of threads used.
             * The default value is 0, implying that the proper
             * number is drawn from the number of available cores on the machine.
             */
            unsigned number_of_workers;

            ///@}

            ///@name Proposal density options
            ///@{

            /*!
             * Perform a decomposition into scan and nuisance parameters
             * and use prior variance for initial proposal of nuisance.
             * The (co)variance in scan direction is unaffected.
             * In addition, set the component's mean to the prior's mean.
             */
            bool block_decomposition;

            /*!
             * Each component has an associated weight,
             * which determines how many samples it may
             * contribute in the next sampling phase.
             */
            std::vector<double> component_weights;

            /// n entries.
            std::vector<std::vector<double>> component_means;

            //TODO better explanation
            /// Shift component by this number of sigmas.
            double component_offset;

            /// n^2 entries.
            std::vector<std::vector<double>> component_variances;

            //todo can we remove this with super_clusters?
            /*!
             * When initializing from global local, draw this
             * many components for each cluster.
             */
            VerifiedRange<unsigned> components_per_cluster;

            /*!
             * Degrees of freedom of a multivariate t-distribution.
             * Special value of -1 corresponds to the Gaussian distribution.
             */
            VerifiedRange<int> degrees_of_freedom;

            /// Used during filtering of components, if too many samples fall outside of box, remove components
            VerifiedRange<double> minimum_overlap;

            /*!
             * Maximum relative distance between two local modes
             * that allows them to be treated as one mode.
             * A rescaling of coordinates to the unit hypercube is
             * performed in order to avoid problems when parameter
             * ranges differ by several orders of magnitude.
             *
             * @example if |x - y|/|y| < eps, they are just one mode
             */
            VerifiedRange<double> mode_distance;

            /// If true, use same local proposal as supplied by GlobalLocal input.
            /// Else use degrees_of_freedom defined above.
            bool override_global_local_proposal;

            /*!
             * Initialize the components randomly at the beginning.
             */
            bool random_start;

            /*!
             * Use only history points from a particular cluster
             * defined by its index. Default value -1: take from all clusters
             */
            int single_cluster;

            /// Skip this percentage from beginning of a chain's history.
            VerifiedRange<double> skip_initial;

            /*!
             * Number of starting positions. Start optimization from each
             * one to identify local modes and covariances.
             */
            unsigned starting_points;

            /*!
             * For random start, take the std. deviation as parameter range divided
             * by std_dev_reduction.
             */
            VerifiedRange<double> std_dev_reduction;

            ///@}

            ///@name Hierarchical clustering options
            ///@{

            /*!
             * If true, use the Gelman-Rubin R-value to combine
             * chains into groups. N initial guess clusters
             * are drawn from each group, no matter how many chains
             * are in it.
             *
             * Default: 1, all chains belong to same group.
             */
            VerifiedRange<double> group_by_r_value;

            /*!
             * Ignore groups (identified by index) from clustering.
             * Default: no groups are ignored.
             * @note Use with great care. Important parts of the posterior may be ignored.
             */
            std::vector<unsigned> ignore_groups;

            /*!
             * When creating patch from a Markov chain history,
             * place the patches around the point with highest posterior
             * in the selected range of point.
             * If false, center on the mean of all points in range.
             */
            bool patch_around_local_mode;

            /*!
             * When grouping with R-value, consider only
             * the scan parameters.
             */
            bool r_value_no_nuisance;

            /*!
             * With hierarchical clustering, use this many
             * samples to form a patch from a single Markov chain.
             */
            unsigned sliding_window;

            /// Store the components created from chain patches.
            bool store_input_components;

            /// Store initial guess for clustering from long patches
            bool store_hc_initial;

            /*!
             * Perform the hierarchical_clustering algorithm on a set of input
             * prerun markov chains to form an initial proposal density.
             * Create this many components for PMC.
             *
             * If group_by_r_value is active, choose this many components per
             * group of chains.
             */
            unsigned super_clusters;

            ///@}

            ///@name Pre run options
            ///@{

            /// Change number of samples according to number of live components.
            bool adjust_sample_size;

            /// The number of updates to proposal functions (mixture components).
            unsigned chunks;

            /// Number of iterations per chunk and component.
            unsigned chunk_size;

            /*!
             *  If > 0, find the samples with highest weight and ignore them during a PMC update.
             *  This can be useful to remove outliers which dominate over the bulk of samples
             *  and lead to many components dying out, as well as to dramatically
             *  shrinking components.
             */
            unsigned crop_highest_weights;

            /// If false, no adaptions to proposal densities are made
            bool need_prerun;

            /// store samples of prerun steps (statistics are always stored)
            bool store_prerun;

            ///@}

            ///@name Convergence options
            ///@{

            /// Declare convergence if both perplexity and effective sample size large enough.
            VerifiedRange<double> convergence_eff_sample_size;
            /// Declare convergence if both perplexity and effective sample size large enough.
            VerifiedRange<double> convergence_perplexity;

            /// Can declare convergence w/o considering the effective sample size,
            /// just by looking at perplexity.
            bool ignore_eff_sample_size;

            /// Declare converge if both perplexity and effective sample size large enough and do not rise anymore.
            VerifiedRange<double> minimum_eff_sample_size;
            /// Declare converge if both perplexity and effective sample size large enough and do not rise anymore.
            VerifiedRange<double> minimum_perplexity;
            /// Consider the last steps to decide if convergence statistics do not rise (significantly) anymore.
            VerifiedRange<unsigned> minimum_steps;
            /// Relative variance of last steps should be small to ascertain convergence.
            VerifiedRange<double> maximum_relative_std_deviation;

            ///@}

            ///@name Main run options
            ///@{

            /// Number of iterations used in final step after adaption is finished.
            long final_chunk_size;

            /// Whether to store collected samples.
            bool store;
            ///@}

            ///@name Output options
            ///@{

            /*!
             * The HDF5 output file to store the markov chains.
             */
            std::string output_file;

            /*!
             *  Determines how often PMC prints outs "Done x%"
             *  during sampling of one chunk.
             *  @example print_steps = 20, then pmc prints
             *  Done 20%
             *  Done 40%
             *  Done 60%
             *  Done 80%
             */
            VerifiedRange<unsigned> print_steps;
            ///@}
    };

    std::ostream & operator<<(std::ostream &, const PopulationMonteCarloSampler::Config &);

    /*!
     * Access to the HDF5 types used in output.
     */
    struct PopulationMonteCarloSampler::Output
    {
         typedef hdf5::Composite<hdf5::Scalar<double>, hdf5::Array<1, double>, hdf5::Array<1, double>> ComponentType;
         typedef hdf5::Scalar<short> IgnoreType;
         typedef hdf5::Array<1, double> SampleType;
         typedef hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>, hdf5::Scalar<double>> StatisticsType;
         typedef hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>> WeightType;

         static ComponentType component_type(const unsigned & dimension);
         static IgnoreType ignore_type();
         static SampleType sample_type(const unsigned & dimension);
         static StatisticsType statistics_type();
         static WeightType weight_type();

         static std::tuple<double, std::vector<double>, std::vector<double>> component_record(const unsigned & dimension);
         static short ignore_record();
         static std::vector<double> sample_record(const unsigned & dimension);
         static std::tuple<double, double, double> statistics_record();
         static std::tuple<double, double> weight_record();
    };

    /*!
     * Holds convergence information of the population sampling.
     */
    struct PopulationMonteCarloSampler::Status
    {
        Status();

        //todo fill and use in [pre]run()
        /// The actual number of samples drawn from mixture proposal density
        /// before the latter is updated.
        unsigned chunk_size;

        /// Whether sampling has converged.
        bool converged;

         //todo fill
        /// The number of iterations after which convergence was declared.
        unsigned iterations_at_convergence;

        /// Evidence.
        double evidence;

        /// Effective sample size.
        double eff_sample_size;

        /// Perplexity.
        double perplexity;
    };
}

#endif
