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

#ifndef EOS_GUARD_EOS_UTILS_PROPROSAL_FUNCTIONS_HH
#define EOS_GUARD_EOS_UTILS_PROPROSAL_FUNCTIONS_HH 1

#include <eos/utils/cluster.hh>
#include <eos/utils/hdf5-fwd.hh>
#include <eos/utils/log_prior-fwd.hh>
#include <eos/utils/markov_chain.hh>
#include <eos/utils/verify.hh>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace eos
{
    namespace proposal_functions
    {
        /*!
         * For a discrete probability distribution defined by its cumulative,
         * draw a random state given by its index in the cumulative.
         */
        unsigned random_index(const std::vector<double> & cumulative, gsl_rng * rng);

        /*!
         * Find the indices [j_min, j_max[ such that they cover range of size
         * within [0, K] around j
         */
        void sliding_window(const unsigned & K, const unsigned & size, const unsigned & j, unsigned & j_min, unsigned & j_max);

        /*!
         * This data type descriptor is needed to identify the proposal function type
         */
        typedef hdf5::Composite<hdf5::Scalar<const char *>, hdf5::Scalar<unsigned>> MetaType;
        MetaType meta_type();

        typedef std::tuple<const char *, unsigned> MetaRecord;
        MetaRecord meta_record();

        /*!
         * UnknownProposalError is thrown when Proposal::make encounters an unknown proposal density name.
         */
        struct UnknownProposalError :
        public Exception
        {
                ///@name Basic Functions
                ///@{
                /*!
                 * Constructor.
                 *
                 * @param name The offending constraint name.
                 */
                UnknownProposalError(const std::string & name);
                ///@}
        };

        /// Abstract factory to read in a proposal density from file
        struct Factory
        {
            static ProposalFunctionPtr make(const hdf5::File & file, const std::string & data_set_base_name,
                                            const std::string & proposal_name, const unsigned & dimension);
        };

        /*!
         * Store long jump vectors and retrieve them efficiently.
         */
        class AdjacencyMatrix
        {
            private:
                std::vector<std::vector<double>> _jump_vectors;
                std::vector<MarkovChain::State> _states;

                unsigned _number_of_clusters;

                /*!
                 * Index of the jump vector related to a jump between i and j,
                 * where i < j is assumed.
                 *
                 * The memory is laid out as in an upper triagonal matrix.
                 * The first element is (0,1) at position 0,
                 * the last is (n-2, n-1) at position n-1.
                 */
                unsigned _index(const unsigned & i, const unsigned & j) const
                {
                    return ((2u * _number_of_clusters - i - 1) * i) / 2 + j - i -1;
                }

            public:
                /// Empty and useless.
                AdjacencyMatrix();

                /// Reserve storage for fixed number of clusters
                void reserve(const unsigned & number_of_clusters);

                /*!
                 *  Add a state to the list, and compute the vector difference
                 *  between the new and all existing ones.
                 *
                 *  @note Fails with an exception if there are enough states already.
                 */

                void add(const MarkovChain::State &);

                /*!
                 * Ignore any differences in dimensions other than those given.
                 *
                 * @example index_list = { 2 }; then only jumps occur only in the 2nd dimension.
                 */
                void indices(const std::vector<unsigned> index_list);

                /*!
                 * Undirected jump excluding the sign.
                 *
                 * @note The caller needs to decide if the sign of the components is to be flipped or not.
                 *
                 * @param h_x Index of the component into which a jump is desired.
                 * @param h_y Index of the component out out which a jump is desired.
                 * @return The vector of the difference in coordinates.
                 */
                const std::vector<double> & jump(const unsigned & h_x, const unsigned & h_y) const;

                const unsigned & number_of_clusters() const
                {
                    return _number_of_clusters;
                }

                /*!
                 * Retrieve the fixed state used for cluster i.
                 */
                const MarkovChain::State & state(const unsigned & i) const;
        };

        /*!
         * After a suitable amount of prerun iterations, combine the information of all
         * chains into one, common proposal function.
         */
        class GlobalLocal :
            public MarkovChain::ProposalFunction
        {
            public:
                /// Types of data sets for reading/writing the state of the proposal density from/to disk
                typedef hdf5::Composite<hdf5::Scalar<unsigned>, hdf5::Array<1, double>> ComponentType;
                typedef hdf5::Composite<hdf5::Array<1, double>, hdf5::Scalar<double>, hdf5::Scalar<double>, hdf5::Scalar<double>> HistoryType;
                typedef hdf5::Composite<hdf5::Array<1, double>, hdf5::Scalar<double>> JumpType;
                typedef hdf5::Composite<hdf5::Array<1, double>, hdf5::Array<1, double>> LocalCovarianceType;

                static ComponentType component_type(const unsigned & dimension);
                static HistoryType history_type(const unsigned & dimension);
                static JumpType jump_type(const unsigned & dimension);
                static LocalCovarianceType local_covariance_type(const unsigned & dimension);

                struct Config
                {
                    private:
                        Config();

                    public:

                        /*!
                         * Choices on how to choose relative weight of
                         * a history point
                         */
                        enum HistoryPointWeighting
                        {
                            posterior,             /// use posterior as weight
                            log_posterior,         /// use log(posterior) as weight
                            equal                  /// all points have same probability
                        };

                        /// Named constructor
                        static Config Default();

                        ///@name Construction options
                        ///@{

                        /// When chains have a common R-value less than this,
                        /// they are considered to overlap sufficiently.
                        VerifiedRange<double> clustering_maximum_r_value;

                        /// Use strict definition of R or relaxed version.
                        bool clustering_strict_r_value;

                        /// Use same probability to draw from each cluster, at least initially.
                        bool equal_weight_components;

                        /*!
                         * Join chains symmetrically such that each chain's history
                         * (excluding the part from skip_initial) has equal contribution
                         * to the covariance matrix
                         */
                        bool join_chains_symmetrically;

                        /// The number of points selected from a single chain's history,
                        /// useful for initializing population Monte Carlo
                        unsigned history_points;

                        /*!
                         * If non-zero, compute covariance around each
                         * history point locally, with the number of samples given.
                         */
                        unsigned history_points_local_covariance_size;

                        /// Choose the points from the ones with highest posterior
                        /// or on a fixed grid (e.g. every 500th point)
                        bool history_points_ordered;

                        /*!
                         * Use log(posterior) to assign relative weight to history points
                         * within each cluster. Else use posterior.
                         */
                        HistoryPointWeighting history_point_weighting;

                        /// If a cluster's weight relative to the "heaviest" cluster
                        /// falls below this threshold, the cluster is discarded.
                        VerifiedRange<double> minimum_relative_cluster_weight;

                        /// Compare chains by their R-values and merge them together.
                        /// Else just put all chains from one partition together.
                        bool perform_clustering;

                        /*!
                         *  In order to increase efficiency, the local proposal
                         *  covariance scale can be divided.
                         *  Default value: 1.0, i.e. don't rescale.
                         *  @note Only meaningful if join_chains_symmetrically is true
                         */
                        VerifiedRange<double> rescale_local_covariance;

                        /// Skip this percentage from beginning of a chain's history.
                        VerifiedRange<double> skip_initial;

                        ///@}

                        ///@name Runtime behavior options
                        ///@{

                        /*!
                         *  Combinining the component probabilities of this and the previous step,
                         *  limit the effect of samples in the distant past by forming a
                         *  weighted average of last component probability, \p_n and the frequency
                         *  encountered in the last chunk, \f_n. In the (n+1)-th step,
                         *
                         *  \p_{n+1} = (1 - 1/n^{cooling_power}) \p_n +  1/n^{cooling_power} * f_n
                         *
                         *  cf. arxiv:0903.0837, Eq. (23)
                         */
                        double cooling_power;

                        /// How often is a local jump proposed instead of a global one.
                        VerifiedRange<double> local_jump_probability;

                        /// If non-empty, all dimensions are masked and changes occur
                        /// only in the dimensions specified by the indices.
                        std::vector<unsigned> long_jump_indices;
                        ///@}
                };

            private:
                unsigned _adaptations;
                Config _config;

                /// Cumulative
                std::vector<double> _component_cumulative;
                std::vector<double> _component_probabilities;

                /// History points selected from the prerun of a chain
                std::vector<std::vector<double>> _history_points_cumulatives;
                std::vector<std::vector<std::vector<double>>> _history_points_local_covariance;
                std::vector<std::vector<double>> _history_points_probabilities;
                std::vector<std::vector<MarkovChain::State>> _history_states;

                /// Stores the jump vectors.
                AdjacencyMatrix _jump_vectors;

                /// the modes of each cluster found during the prerun.
                std::vector<MarkovChain::State> _modes;

                std::vector<std::shared_ptr<MarkovChain::ProposalFunction>> _prop;

                void _sanity_check() const;

                /*!
                 * Choose a number of points from the history of the chains in each
                 * cluster. The points are unique with probability 1 and represent
                 * the points with highest posterior probability.
                 *
                 */
                void _select_history_points(const std::vector<Cluster> & clusters);

                /*!
                 * Find the modes from the points and determine jumps as
                 * the vector differences between the modes.
                 */
                void _select_jump_vectors(const std::vector<Cluster> & clusters);

                typedef std::pair<unsigned, double> index_pair;

                //todo replace by lambda
                bool
                static _cmp(const index_pair & a, const index_pair & b )
                {
                    return a.second > b.second;
                }

                /*!
                 * Copy constructor for cloning only.
                 */
                GlobalLocal(const GlobalLocal &);

            public:

                /*!
                 * Take the important parts of MarkovChains, and combine them into
                 * a global local proposal density. It is assumed that in each vector,
                 * the parts at a given index belong to the same chain.
                 * If their histories only contain the last chunk, variability
                 * of the algorithm is limited.
                 *
                 * @par chains Use pointers to the histories, such that unnecessary copying is avoided.
                 * @par proposal The proposal density of each chain.
                 * @par stats The statistics of each chain.
                 * @par local_jump_probability In order to allow chains to explore more distant regions,
                 *                             a value close to one is advised.
                 * @par history_points the number of points selected from a single chain's history
                 *
                 */
                // todo change to use shared_ptr
                GlobalLocal(const std::vector<HistoryPtr> & chains,
                            const std::vector<ProposalFunctionPtr> & proposals,
                            const std::vector<MarkovChain::Stats> & stats,
                            const Config & config, const unsigned & prerun_chains_per_cluster);

                GlobalLocal(const std::vector<double> & component_probabilities, const unsigned & adaptations,
                            const std::vector<MarkovChain::State> & jump_states, const std::vector<MarkovChain::State> & modes,
                            const std::vector<std::vector<MarkovChain::State>> & history_states,
                            const std::vector<std::vector<double>> & history_point_probabilities,
                            const std::vector<std::vector<std::vector<double>>> & history_local_covariances,
                            const std::vector<ProposalFunctionPtr> & proposals);

                virtual ~GlobalLocal();

                /// No adaptation currently foreseen.
                virtual void adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end, const double & efficiency, const double & efficiency_min, const double & efficiency_max);

                /*!
                 * The state with the highest posterior of all components
                 * with parameter values and correct log posterior.
                 */
                virtual MarkovChain::State mode() const;

                const std::vector<double> & component_probabilities() const
                {
                    return _component_probabilities;
                }

                const std::vector<std::vector<MarkovChain::State>> & history_states() const
                {
                    return _history_states;
                }

                const std::vector<std::vector<std::vector<double>>> & local_covariances() const
                {
                    return _history_points_local_covariance;
                }

                /// Reset config options
                virtual void config(const Config & config);

                virtual ProposalFunctionPtr clone() const;

                virtual void dump_state(hdf5::File & file, const std::string & data_set_base_name) const;

                virtual double evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const;

                virtual void propose(MarkovChain::State & x, const MarkovChain::State & y, gsl_rng * rng) const;
        };

        class Multivariate :
            public MarkovChain::ProposalFunction
        {
            public:
                typedef hdf5::Array<1, double> CovarianceType;
                static CovarianceType covariance_type(const unsigned & dimension);

            protected:
                gsl_vector * _tmp_left;
                gsl_vector * _tmp_right;
                gsl_matrix * _tmp_sample_covariance_current;

                gsl_matrix * _covariance;
                gsl_matrix * _covariance_inverse;
                gsl_matrix * _covariance_chol;

                const unsigned _dimension;

                std::vector<unsigned> _index_list;

                void _compute_cholesky_and_inverse();
                virtual void _compute_norm() = 0;
                void _copy(const Multivariate &);
                void _dump_covariance(hdf5::File & file, const std::string & data_set_base_name, const std::string & proposal_type_name) const;

            public:
                /// The dimension of the space for which samples are proposed
                unsigned dimension() const;

                /// Record how often an adaption to data has been performed.
                unsigned adaptations;

                const gsl_matrix * covariance() const;

                /// Rescale the sample _covariance to form the proposal _covariance.
                double covariance_scale;

                /// Scale enforced to exceed a minimum value.
                static const double covariance_scale_min;
                /// Scale enforced to lie below a maximum value.
                static const double covariance_scale_max;
                /// During an adaptation, the scale is multiplied/divided by this factor
                /// if the efficiency is too high/low .
                static const double covariance_scale_update_factor;

                /*!
                 *  Combinining the proposal _covariance of this and the previous step,
                 *  limit the effect of samples in the distant past by forming a
                 *  weighted average. In the n-th step,
                 *
                 *  \Sigma_n = (1 - 1/n^{cooling_power}) \Sigma_n +  1/n^{cooling_power} * S_n
                 *
                 *  cf. arxiv:0903.0837, Eq. (23)
                 */
                double cooling_power;

                double norm;

                /*!
                 * Construct with given dimension and an initial proposal covariance matrix.
                 * If automatic scaling is desired, the given matrix is interpreted as
                 * an estimate of the target covariance, and it is rescaled by a factor of
                 * 2.38^2 / dimension. Else is the matrix is used as is
                 * for proposing points.
                 *
                 * @note The automatic scaling can be shown to lead to the asymptotically (!)
                 * perfect value of the efficiency of 0.238 for a Gaussian proposal and
                 * a Gaussian target density with unlimited parameter ranges.
                 * These requirements are usually not fullfilled.
                 */
                Multivariate(const unsigned & dimension, const std::vector<double> & covariance, const bool & automatic_scaling/* = false*/);

                virtual ~Multivariate();

                /*!
                 * @note Expect history to have only the most recent part, to which proposal
                 * function has not adapted yet. Accordingly, the efficiency
                 * of that last part is required.
                 * Do not push the entire history of a chain in here repeatedly!
                 */
                virtual void adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                                   const double & efficiency, const double & efficiency_min, const double & efficiency_max);

                /*!
                 * Rescale the covariance scale factor
                 * @param rescale_factor
                 */
                virtual void rescale(const double & rescale_factor);

                /*!
                 * Reset the internal covariance to zero, and recompute it from the histories given.
                 * Apply a new scale factor for proposals.
                 * @param histories
                 * @param scale
                 */
                virtual void reset(const std::vector<HistoryPtr> & histories, const double & scale, const double & skip_initial);

                /*!
                 * Set an index list, such that the proposal will consider only a
                 * subspace for it covariance calculation.
                 */
                virtual void set_indices(const std::vector<unsigned> & index_list);
        };

        typedef std::shared_ptr<Multivariate> MultivariateProposalPtr;

        class MultivariateGaussian :
            public Multivariate
        {
            public:
                typedef hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>, hdf5::Scalar<unsigned>> ScalarsType;
                static ScalarsType scalars_type();

            protected:
                virtual void _compute_norm();

            public:
                MultivariateGaussian(const unsigned & dimension, const std::vector<double> & covariance, const bool & automatic_scaling = true);

                virtual ~MultivariateGaussian();

                virtual ProposalFunctionPtr clone() const;

                virtual void dump_state(hdf5::File & file, const std::string & data_set_base_name) const;

                virtual double evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const;

                virtual void propose(MarkovChain::State & x, const MarkovChain::State & y, gsl_rng * rng) const;
        };

        class MultivariateStudentT :
            public Multivariate
        {
            public:
                typedef hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>, hdf5::Scalar<unsigned>, hdf5::Scalar<double>> ScalarsType;
                static ScalarsType scalars_type();

            protected:

                virtual void _compute_norm();

            public:
                /// Degree of freedom.
                /// Required to be >= 1.
                VerifiedRange<double> dof;

                MultivariateStudentT(const unsigned & dimension, const std::vector<double> & covariance, const double & degree_of_freedom, const bool & automatic_scaling = true);

                virtual ~MultivariateStudentT();

                virtual ProposalFunctionPtr clone() const;

                virtual void dump_state(hdf5::File & file, const std::string & data_set_base_name) const;

                virtual double evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const;

                virtual void propose(MarkovChain::State & x, const MarkovChain::State & y, gsl_rng * rng) const;
        };

        struct MultivariateAccess
        {
            static MultivariateProposalPtr access(const ProposalFunctionPtr & p);
        };

        class BlockDecomposition :
            public MarkovChain::ProposalFunction
        {
            public:
                typedef hdf5::Composite<hdf5::Scalar<const char *>> PriorsType;
                static PriorsType priors_type();
                friend class MultivariateAccess;
            private:
                // use these vectors to represent the full length and ordering
                // as in Analysis
                std::vector<std::shared_ptr<double> > _tmp_vector_x;
                std::vector<std::shared_ptr<double> > _tmp_vector_y;

                mutable MarkovChain::State _tmp_state_x;
                mutable MarkovChain::State _tmp_state_y;

                // the actual proposal functions
                std::vector<MultivariateProposalPtr> _mv;
                std::vector<LogPriorPtr> _priors;

                Parameters p;

                // keep the references right
                std::vector<std::vector<std::shared_ptr<double>>> _mv_values_x;
                std::vector<std::vector<std::shared_ptr<double>>> _mv_values_y;
                std::vector<std::vector<std::shared_ptr<double>>> _priors_values_x;
                std::vector<std::vector<std::shared_ptr<double>>> _priors_values_y;

                /// copy values pointed at into result
                static void _copy_values(const std::vector<std::shared_ptr<double>> & source_ptr_vector, std::vector<double> & result);

                /*!
                 * Copy values from source vector to where the pointers in destination point
                 * @param source
                 * @param destination_ptr_vector
                 */
                static void _copy_values(const std::vector<double> & source, const std::vector<std::shared_ptr<double>> & destination_ptr_vector);
            public:
                BlockDecomposition();
                virtual ~BlockDecomposition();

                virtual void adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                                   const double & efficiency, const double & efficiency_min, const double & efficiency_max);

                void add(const MultivariateProposalPtr &);
                void add(const LogPriorPtr &);

                virtual ProposalFunctionPtr clone() const;

                virtual void dump_state(hdf5::File & file, const std::string & data_set_base_name) const;

                virtual double evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const;

                virtual void propose(MarkovChain::State & x, const MarkovChain::State & y, gsl_rng * rng) const;
        };
    }
}

#endif
