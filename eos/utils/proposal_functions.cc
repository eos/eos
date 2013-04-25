/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
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

#include <eos/utils/proposal_functions.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/histogram.hh>
#include <eos/utils/log.hh>
#include <eos/utils/log_prior.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/rvalue.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <numeric>

namespace eos
{
    namespace proposal_functions
    {
        /*
         *  Assuming  the cumulative of a discrete probability distribution, return an index i
         *  with frequency p[i]
         *
         *  Modeled after the very readable reference
         *  http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary
         *
         *  @par cumulative sorted in ascending order, e.g. cum[0] = 0.2, cum[1] = 0.3. Last value must be = 1.0
         *  @note only works with vectors of length <= max(int)
         *
         */
        unsigned
        random_index(const std::vector<double> & cumulative, gsl_rng * rng)
        {
            // generate random number, uniform on [0,1]
            double u = gsl_ran_flat(rng, 0.0, 1.0);

            // binary search for it
            int high(cumulative.size());
            int low = -1, probe;
            while (high - low > 1)
            {
                // fix overflow bug
                probe = ((unsigned int)low + (unsigned int)high) >> 1;
                if (cumulative[probe] < u)
                    low = probe;
                else
                    high = probe;
            }

            return high;
        }

        void sliding_window(const unsigned & K, const unsigned & size, const unsigned & j, unsigned & j_min, unsigned & j_max)
        {
            if (size >= K)
                throw InternalError("prop::sliding_window: size >= K (" + stringify(size) + " vs " + stringify(K) + ")");
            if (j >= K)
                throw InternalError("prop::sliding_window: j >= K (" + stringify(j) + " vs " + stringify(K) + ")");

            // initial guess
            j_min = j - (size / 2 - 1);
            j_max = j + (size / 2 + 1);

            // overflow
            if (j_min > K)
            {
                j_min = 0;
                j_max = size;
            }
            if (j_max > K)
            {
                j_max = K;
                j_min = K - size;
            }
            if (j_max - j_min != size)
                throw InternalError("prop::sliding_window: Couldn't adjust sizes");
        }

        std::string
        print_matrix(gsl_matrix * const m)
        {
            std::string result("\n");

            for (unsigned i = 0 ; i < m->size1 ; ++i)
            {
                result += '[';
                for (unsigned j = 0 ; j < m->size2 ; ++j)
                {
                    result += stringify(gsl_matrix_get(m, i, j));
                    if (j != m->size2 - 1)
                    {
                        result += ", ";
                    }
                }
                result += "]\n";
            }
            return result;
        }

        UnknownProposalError::UnknownProposalError(const std::string & name) :
            Exception("Proposal '" + name + "' is unknown")
        {
        }

        MetaType
        meta_type()
        {
            return MetaType {"meta", hdf5::Scalar<const char *>("proposal type"), hdf5::Scalar<unsigned>("number of dimensions")};
        }

        MetaRecord
        meta_record()
        {
            return std::make_tuple("prop", 0u);
        }

        /*!
         * Peek inside a proposal function and retrieve the part that acts like a multivariate
         * proposal density
         *
         * todo Find a way around using dynamic casts
         */
        MultivariateProposalPtr
        MultivariateAccess::access(const ProposalFunctionPtr & p)
        {
            MultivariateProposalPtr mv;

            {
                mv = std::dynamic_pointer_cast<Multivariate>(p);
                if (mv.get())
                    return mv;
            }

            {
                auto bd = std::dynamic_pointer_cast<BlockDecomposition>(p);
                if (bd.get())
                    return bd->_mv.front();
            }

            // default case
            throw InternalError("MultivariateAccess: couldn't find type of proposal pointer");
        }

        struct GlobalLocalFactory
        {
            static ProposalFunctionPtr make(hdf5::File & file, const std::string & data_set_base_name, const unsigned & dimension)
            {
                /* read in local proposal functions */

                std::vector<ProposalFunctionPtr> proposals;
                H5E_BEGIN_TRY
                {
                    try
                    {
                        unsigned i = 0;
                        while(true)
                        {
                            std::string sub_directory = data_set_base_name + "/local proposals/" + stringify(i);
                            // need to find out the right type
                            auto meta_data_set = file.open_data_set(sub_directory + "/meta", meta_type());
                            auto meta_record = proposal_functions::meta_record();
                            meta_data_set >> meta_record;

                            if (std::get<1>(meta_record) != dimension)
                                throw InternalError("Factory::make_global_local: current dimension(" + stringify(dimension) +
                                    ") doesn't match that in proposal (" + stringify(std::get<1>(meta_record)) + ").");

                            // use the factory again with the name of the local proposal type
                            proposals.push_back(Factory::make(file, sub_directory, std::get<0>(meta_record), dimension));

                            // update
                            ++i;
                        }
                    }
                    catch (HDF5Error &)
                    {
                    }
                }
                H5E_END_TRY;

                //need to know number of components: get it from number of local proposals
                unsigned n_components = proposals.size();

                /* read components probabilities */
                auto comp_data_set = file.open_data_set(data_set_base_name + "/components", GlobalLocal::component_type(n_components));

                auto record_comp =  std::make_tuple(0u, std::vector<double>(n_components));
                comp_data_set.end();
                comp_data_set >> record_comp;
                unsigned adaptations = std::get<0>(record_comp);
                std::vector<double> component_probabilities(std::get<1>(record_comp));

                /* read the points for long jumps */

                GlobalLocal::JumpType jump_type = GlobalLocal::jump_type(dimension);
                auto jump_data_set = file.open_data_set(data_set_base_name + "/jump", jump_type);
                auto jump_record = std::make_tuple(std::vector<double>(dimension), 1.0);

                std::vector<MarkovChain::State> jump_states;

                for (unsigned i = 0 ; i < n_components ; ++i)
                {
                    jump_data_set >> jump_record;
                    MarkovChain::State state;
                    state.point = std::get<0>(jump_record);
                    state.log_posterior = std::get<1>(jump_record);
                    state.hyper_parameter.component = i;

                    jump_states.push_back(state);
                }

                /* read the local cluster modes */

                auto mode_data_set = file.open_data_set(data_set_base_name + "/modes", jump_type);

                std::vector<MarkovChain::State> mode_states;

                for (unsigned i = 0 ; i < n_components ; ++i)
                {
                    mode_data_set >> jump_record;
                    MarkovChain::State state;
                    state.point = std::get<0>(jump_record);
                    state.log_posterior = std::get<1>(jump_record);
                    state.hyper_parameter.component = i;

                    mode_states.push_back(state);
                }

                /* read the intermediate points */

                GlobalLocal::HistoryType history_type
                {
                    "history point",
                    hdf5::Array<1, double>("point", { dimension}),
                    hdf5::Scalar<double>("log posterior"),
                    hdf5::Scalar<double>("probability"),
                    hdf5::Scalar<double>("cumulative"),
                };
                auto record_history = std::make_tuple(std::vector<double>(dimension), 1.0, 2.0, 3.0);

                std::vector<std::vector<MarkovChain::State>> all_clusters_history_states;
                std::vector<std::vector<double>> all_clusters_point_probabilities;

                for (unsigned i = 0 ; i < n_components ; ++i)
                {
                    std::vector<MarkovChain::State> component_history;
                    std::vector<double> componenent_point_probabilities;

                    MarkovChain::State state;

                    auto history_data_set = file.open_data_set(data_set_base_name + "/history/" + stringify(i), history_type);
                    for (unsigned p = 0 ; p < history_data_set.records() ; ++p)
                    {
                        history_data_set >> record_history;
                        state.point = std::get<0>(record_history);
                        state.log_posterior = std::get<1>(record_history);
                        component_history.push_back(state);
                        componenent_point_probabilities.push_back(std::get<2>(record_history));
                    }
                    all_clusters_history_states.push_back(component_history);
                    all_clusters_point_probabilities.push_back(componenent_point_probabilities);
                }

                /* parse local covariances , they are optional*/

                std::vector<std::vector<std::vector<double>>> local_covariances;
                H5E_BEGIN_TRY
                {
                    try
                    {
                        auto record_local_cov = std::make_tuple(std::vector<double>(dimension), std::vector<double>(dimension * dimension));
                        for (unsigned c = 0 ; c < n_components ; ++c)
                        {
                            auto local_cov_data_set = file.open_data_set(data_set_base_name + "/local covariances/" + stringify(c), GlobalLocal::local_covariance_type(dimension));
                            std::vector<std::vector<double>> component_covariances;
                            for (unsigned i = 0 ; i < local_cov_data_set.records() ; ++i)
                            {
                                local_cov_data_set >> record_local_cov;
                                component_covariances.push_back(std::get<1>(record_local_cov));
                            }
                            local_covariances.push_back(component_covariances);
                        }
                    }
                    catch (HDF5Error &)
                    {
                    }
                }
                H5E_END_TRY;


                // collected all relevant info from the file
                return std::make_shared<GlobalLocal>(component_probabilities, adaptations,
                                                     jump_states, mode_states,
                                                     all_clusters_history_states, all_clusters_point_probabilities,
                                                     local_covariances,
                                                     proposals);
            }
        };

        struct MultivariateGaussianFactory
        {
            static ProposalFunctionPtr make(hdf5::File & file, const std::string & data_set_base_name, const unsigned & dimension)
            {
                // read in covariance
                std::vector<double> covariance (dimension * dimension);
                Multivariate::CovarianceType covariance_type
                { hdf5::Array<1, double>("covariance matrix", { dimension * dimension }) };
                auto cov_data_set = file.open_data_set(data_set_base_name + "/covariance", covariance_type);

                //jump to last record
                cov_data_set.end();
                cov_data_set >> covariance;

                // read in Scalars
                auto scalars_data_set = file.open_data_set(data_set_base_name + "/scalars", MultivariateGaussian::scalars_type());

                auto scalars = std::make_tuple(0.0, 0.0, 0u);
                scalars_data_set.end();
                scalars_data_set >> scalars;

                // create the object and set its properties, but don't rescale again
                MultivariateGaussian * p = new MultivariateGaussian(dimension, covariance, false);

                p->covariance_scale = std::get<0>(scalars);
                p->cooling_power = std::get<1>(scalars);
                p->adaptations = std::get<2>(scalars);

                return ProposalFunctionPtr(p);
            }
        };

        struct MultivariateStudentTFactory
        {
            static ProposalFunctionPtr make(hdf5::File & file, const std::string & data_set_base_name, const unsigned & dimension)
            {
                // read in covariance
                /*
                 * Check that results for evaluate(), propose() agree with and without dumping
                 */
                std::vector<double> covariance (dimension * dimension);
                Multivariate::CovarianceType covariance_type
                { hdf5::Array<1, double>("covariance matrix", { dimension * dimension }) };
                auto cov_data_set = file.open_data_set(data_set_base_name + "/covariance", covariance_type);

                //jump to last record
                for (unsigned i = 0 ; i < cov_data_set.records() ; ++i)
                {
                    cov_data_set >> covariance;
                }

                // read in Scalars
                auto scalars_data_set = file.open_data_set(data_set_base_name + "/scalars", MultivariateStudentT::scalars_type());

                auto scalars = std::make_tuple(0.0, 0.0, 0u, 0.0);
                scalars_data_set >> scalars;

                // create the object and set its properties, but don't rescale covariance again
                MultivariateStudentT * p = new MultivariateStudentT(dimension, covariance, std::get<3>(scalars), false);

                p->covariance_scale = std::get<0>(scalars);
                p->cooling_power = std::get<1>(scalars);
                p->adaptations = std::get<2>(scalars);

                return ProposalFunctionPtr(p);
            }
        };

        /*!
         * Make an instance of BlockDecomposition from information in a HDF5 file.
         */
        struct BlockDecompositionFactory
        {
                static ProposalFunctionPtr make(hdf5::File & file, const std::string & data_set_base_name, const unsigned & )
                {
                    BlockDecomposition * bd = new BlockDecomposition();

                    // read in multivariates
                    //todo support only one multivariate for now, because number_of_objects is not reliable/understood
                    for (unsigned i = 0 ; i < 1 /*file.number_of_objects(data_set_base_name + "/multivariates") - 1*/; ++i)
                    {
                        auto meta_mv_data_set = file.open_data_set(data_set_base_name + "/multivariates/" + stringify(i) + "/meta", meta_type());
                        auto meta_mv_record = meta_record();
                        meta_mv_data_set >> meta_mv_record;
                        ProposalFunctionPtr mv = Factory::make(file, data_set_base_name + "/multivariates/" + stringify(i),
                                                               std::get<0>(meta_mv_record), std::get<1>(meta_mv_record));
                        bd->add(std::static_pointer_cast<Multivariate>(mv));
                    }

                    // read in priors
                    {
                        auto data_set = file.open_data_set(data_set_base_name + "/priors", BlockDecomposition::priors_type());
                        auto record = std::make_tuple("serialized prior");
                        Parameters p = Parameters::Defaults();
                        for (unsigned i = 0 ; i < data_set.records() ; ++i)
                        {
                            data_set >> record;
                            LogPriorPtr prior = LogPrior::Make(p, std::get<0>(record));
                            bd->add(prior);
                        }
                    }

                    return ProposalFunctionPtr(bd);
                }
        };

        typedef std::function< ProposalFunctionPtr (hdf5::File &, const std::string & data_set_base_name, const unsigned & dimension)> ProposalFactory;

        template <typename Factory_> ProposalFactory
        make_factory(const Factory_ &)
        {
            return std::bind(&Factory_::make, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }

        ProposalFunctionPtr
        Factory::make(const hdf5::File & file, const std::string & data_set_base_name,
                      const std::string & proposal_name, const unsigned & dimension)
        {
            static const std::map<std::string, ProposalFactory> factories
            {
                { "BlockDecomposition",  make_factory(BlockDecompositionFactory())},
                { "GlobalLocal",  make_factory(GlobalLocalFactory())},
                { "MultivariateGaussian", make_factory(MultivariateGaussianFactory()) },
                { "MultivariateStudentT", make_factory(MultivariateStudentTFactory())}
            };

            auto f = factories.find(proposal_name);
            if (f == factories.end())
                throw UnknownProposalError(proposal_name);

            return f->second(const_cast<hdf5::File &>(file), data_set_base_name, dimension);
        }

        AdjacencyMatrix::AdjacencyMatrix() :
            _number_of_clusters(0)
        {
        }

        void
        AdjacencyMatrix::reserve(const unsigned & number_of_clusters)
        {
            _jump_vectors.resize(number_of_clusters * (number_of_clusters - 1) / 2);
            _number_of_clusters = number_of_clusters;
        }

        void
        AdjacencyMatrix::add(const MarkovChain::State & state)
        {
            if (_states.size() == _number_of_clusters)
                throw InternalError("AdjacencyMatrix::add: cannot add another state");

            // store jumps row wise, i.e. first (n-1) jumps belong to cluster 0,
            // next (n-2) belong to cluster 1 ...

            std::vector<double> difference(state.point.size());

            // loop over all other points and compute vector difference
            // fill a column until diagonal in adjacency matrix
            for (auto i = _states.cbegin() ; i != _states.cend() ; ++i)
            {
                std::transform(state.point.cbegin(), state.point.cend(), i->point.cbegin(),
                        difference.begin(), std::minus<double>());

                // second index j is the one that point will have once it is added
                _jump_vectors.at(_index(std::distance(_states.cbegin(), i), _states.size())) = difference;
            }

            _states.push_back(state);
        }

        void
        AdjacencyMatrix::indices(const std::vector<unsigned> index_list)
        {
            // loop over parameters
            for (unsigned i = 0 ; i < _jump_vectors.front().size() ; ++i)
            {
                auto res = std::find(index_list.begin(), index_list.end(), i);

                if (res != index_list.end())
                    continue;

                // loop over all jump vectors and eliminate changes in pars not mentioned in index_list
                for (auto jump = _jump_vectors.begin() ; jump != _jump_vectors.end() ; ++jump)
                {
                    (*jump)[i] = 0.0;
                }
            }
        }

        const std::vector<double> &
        AdjacencyMatrix::jump(const unsigned & h_x, const unsigned & h_y) const
        {
            if (h_x == h_y)
                throw InternalError("AdjacencyMatrix::jump: jumps within one component (" + stringify(h_x) + ") are not implemented yet");
            if (h_x < h_y)
                return _jump_vectors[_index(h_x, h_y)];
            else
                return _jump_vectors[_index(h_y, h_x)];
        }

        const MarkovChain::State &
        AdjacencyMatrix::state(const unsigned & i) const
        {
            return _states.at(i);
        }

        void
        GlobalLocal::_sanity_check() const
        {
            unsigned n = _component_cumulative.size();
            if (n != _component_probabilities.size())
            {
                throw InternalError("GlobalLocal.ctor: n != comp_prob ( n = " + stringify(n) + ")");
            }
            if (n != _jump_vectors.number_of_clusters())
            {
                throw InternalError("GlobalLocal.ctor: n != jump_vectors (n = " + stringify(n) + " vs "
                                    + stringify(_jump_vectors.number_of_clusters()) + " )");
            }
            if (n != _prop.size())
            {
                throw InternalError("GlobalLocal.ctor: n != prop ( n = " + stringify(n) + ")");
            }
        }

        void
        GlobalLocal::_select_history_points(const std::vector<Cluster> & clusters)
        {
            unsigned k = 0;
            for (auto cl = clusters.cbegin(), cl_end = clusters.cend() ; cl != cl_end ; ++cl, ++k)
            {
                Log::instance()->message("GlobalLocal.select", ll_debug)
                    << "Selecting " << _config.history_points << " history points from cluster " << k << ":";

                MarkovChain::History history;

                // just count how many there are
                unsigned n_chains_in_cluster = 0;

                // loop over each chain and extract history
                for (auto c = cl->begin(), c_end = cl->end() ; c != c_end ; ++c, ++n_chains_in_cluster)
                {
                    // copy log(posterior) and associate with index
                    if ((**c).states.empty())
                    {
                        throw InternalError("proposal_functions::GlobalLocal: cannot select points from empty history");
                    }

                    std::copy((**c).states.cbegin() + unsigned(_config.skip_initial * (**c).states.size()),
                              (**c).states.cend(), std::back_inserter(history.states));
                }

                // associate points with index
                std::vector<index_pair> posterior_indices(history.states.size());
                for (unsigned i = 0 ; i < posterior_indices.size() ; ++i)
                {
                    posterior_indices[i] = index_pair(i, history.states[i].log_posterior);
                }

                // sort according to posterior in descending order
                std::sort(posterior_indices.begin(), posterior_indices.end(), _cmp);

                std::vector<MarkovChain::State> states(_config.history_points);
                std::vector<double> probabilities(_config.history_points, 0.0);

                // the position of each history state that was actually selected
                // in the list of all states
                std::vector<unsigned> state_indices(_config.history_points);

                // initialize to NaN, so first comparison has to fail
                double previous_posterior = __builtin_nan("");

                if (_config.history_points_ordered)
                {
                    // select points in variable order until a different one is found
                    unsigned i = -1;
                    for (unsigned j = 0 ; j < _config.history_points ; ++j)
                    {
                        while (true)
                        {
                            state_indices[j] = posterior_indices[++i].first;
                            const MarkovChain::State & s = history.states[state_indices[j]];

                            // look for states that differ, since point can appear multiple times in a row
                            if (s.log_posterior != previous_posterior)
                            {
                                previous_posterior = s.log_posterior;
                                states[j] = s;
                                break;
                            }
                        }
                    }
                }
                else
                {
                    // choose points only from upper quantile
                    unsigned increment = history.states.size() / _config.history_points / 2;
                    for (unsigned j = 0 ; j < _config.history_points ; ++j)
                    {
                        state_indices[j] = posterior_indices[j * increment].first;
                        states[j] = history.states[state_indices[j]];
                    }
                }
                // find minimum of points to be chosen, but subtract a few percent s.t. last point doesn't have weight zero
                const double min_posterior = 0.95 * states.back().log_posterior;
                const double max_posterior = states.front().log_posterior;

                for (unsigned j = 0 ; j < _config.history_points ; ++j)
                {
                    // assign relative weights to points
                    switch (_config.history_point_weighting) {
                        case Config::HistoryPointWeighting::equal:
                            probabilities[j] = 1;
                            break;

                        case Config::HistoryPointWeighting::log_posterior:
                            // rescale by minimum of history points(!) so all weights are positive
                            probabilities[j] = states[j].log_posterior - min_posterior;
                            break;

                        case Config::HistoryPointWeighting::posterior:
                        default:
                            // Rescale by entire maximum, s.t. all relative probabilities in [0,1]
                            // to avoid overflows with exp(...)
                            probabilities[j] = std::exp(states[j].log_posterior - max_posterior);
                            break;
                    }
                }

                std::vector<std::vector<double>> component_local_covariances;
                if (_config.history_points_local_covariance_size > 0)
                {
                    const unsigned & n_dim = states.front().point.size();

                    // find local covariance for each history point
                    for (unsigned j = 0 ; j < _config.history_points ; ++j)
                    {
                        // determine actual chain where point is from
                        // assume each chain yielded same amount of samples, K.
                        // Then a sample from chain 1 is in index [K, 2K[
                        const unsigned single_chain_length = history.states.size() / n_chains_in_cluster;
                        auto single_history = cl->begin();
                        for (unsigned n = 0 ; n < state_indices[j] / single_chain_length; ++n)
                            single_history++;

                        // find suitable environment around it:
                        // ideally half of the samples up, half of them down,
                        // but if point at the very end, have to take all samples before the actual point
                        unsigned min_index, max_index;
                        const unsigned chain_index = state_indices[j] / single_chain_length;
                        proposal_functions::sliding_window(single_chain_length, _config.history_points_local_covariance_size,
                                                           state_indices[j] - chain_index * single_chain_length,
                                                           min_index, max_index);

                        // skip first part and shift indices
                        const unsigned offset = _config.skip_initial * (**single_history).states.size();
                        const auto begin_state = (**single_history).states.begin() + offset + min_index;
                        const auto end_state   = (**single_history).states.begin() + offset + max_index;

                        // compute mean and variance along each dimension
                        std::vector<double> means(n_dim);
                        std::vector<double> variances(n_dim);

                        (**single_history).mean_and_variance(begin_state, end_state, means, variances);

                        // compute local_covariance (only off-diagonal elements)
                        std::vector<double> local_covariance(n_dim * n_dim, 0.0);
                        for (auto s = begin_state ;  s != end_state ; ++s)
                        {
                            for (unsigned dim1 = 0 ; dim1 < n_dim ; ++dim1)
                            {
                                // off-diagonal elements
                                for (unsigned dim2 = dim1 + 1 ; dim2 < n_dim ; ++dim2)
                                {
                                    const double summand = (s->point[dim1] - means[dim1]) * (s->point[dim2] - means[dim2]);
                                    local_covariance[dim1 + dim2 * n_dim] += summand;
                                    local_covariance[dim2 + dim1 * n_dim] += summand;
                                }
                            }
                        }

                        // rescale for unbiased estimate 1 / (N - 1)
                        std::transform(local_covariance.begin(), local_covariance.end(), local_covariance.begin(),
                                       std::bind(std::multiplies<double>(), 1.0 / (_config.history_points_local_covariance_size - 1),
                                                 std::placeholders::_1));

                        // diagonal elements are normalized already
                        for (unsigned dim1 = 0 ; dim1 < n_dim ; ++dim1)
                        {
                            local_covariance[dim1 + dim1 * n_dim] = variances[dim1];
                        }
                           component_local_covariances.push_back(local_covariance);
                    }
                }

                // rescale cumulative such that last value is one
                std::vector<double> cumulative;
                std::partial_sum(probabilities.cbegin(), probabilities.cend(),
                                 std::back_inserter(cumulative));
                std::transform(probabilities.begin(), probabilities.end(), probabilities.begin(),
                               std::bind(std::multiplies<double>(), 1.0 / cumulative.back(), std::placeholders::_1));
                std::transform(cumulative.begin(), cumulative.end(), cumulative.begin(),
                               std::bind(std::multiplies<double>(), 1.0 / cumulative.back(), std::placeholders::_1));

                Log::instance()->message("GlobalLocal.select", ll_debug)
                    << "First 5 point probabilities: "
                    << stringify(probabilities.begin(), probabilities.begin() + 5, 5);

                _history_states.push_back(states);
                _history_points_cumulatives.push_back(cumulative);
                _history_points_probabilities.push_back(probabilities);
                if (! component_local_covariances.empty())
                    _history_points_local_covariance.push_back(component_local_covariances);
            }
        }

        void
        GlobalLocal::_select_jump_vectors(const std::vector<Cluster> & clusters)
        {
            if (clusters.size() != _modes.size())
                throw InternalError("prop::gl::select: cluster and modes don't match: "
                                    + stringify(clusters.size() + " vs " + stringify(_modes.size())));

            // create empty adjacency matrix
            _jump_vectors.reserve(clusters.size());

            // compute mean of cluster and use it for translation vectors
            for (auto cl = clusters.cbegin() ; cl != clusters.cend() ; ++cl)
            {
                MarkovChain::State s;
                s.point = cl->mean();
                _jump_vectors.add(s);
            }

            // apply masking if needed
            if ( ! _config.long_jump_indices.empty())
            {
                _jump_vectors.indices(_config.long_jump_indices);
            }

            for (unsigned i = 0 ; i < clusters.size() ; ++i)
            {
                for (unsigned j = i + 1 ; j < clusters.size() ; ++j)
                {
                    Log::instance()->message("GlobalLocal.select", ll_debug)
                        << "vec between " << i << " and " << j << " is " <<
                        stringify(_jump_vectors.jump(i,j).cbegin(), _jump_vectors.jump(i,j).cbegin()
                                     + std::min(3u, unsigned(_jump_vectors.jump(i,j).size())));
                }
            }
        }

        GlobalLocal::GlobalLocal(const GlobalLocal & other) :
            _adaptations(other._adaptations),
            _config(other._config),
            _component_cumulative(other._component_cumulative),
            _component_probabilities(other._component_probabilities),
            _history_points_cumulatives(other._history_points_cumulatives),
            _history_points_local_covariance(other._history_points_local_covariance),
            _history_points_probabilities(other._history_points_probabilities),
            _history_states(other._history_states),
            _jump_vectors(other._jump_vectors)
        {
            for ( auto i = other._prop.begin(), i_end = other._prop.end() ; i != i_end ; ++i)
            {
                _prop.push_back((*i)->clone());
            }
        }

        GlobalLocal::ComponentType
        GlobalLocal::component_type(const unsigned & dimension)
        {
            return ComponentType{"components", hdf5::Scalar<unsigned>("adaptations"),  hdf5::Array<1, double>("probability", { dimension }),};
        }
        GlobalLocal::HistoryType
        GlobalLocal::history_type(const unsigned & dimension)
        {
            return
            HistoryType
            {
                "history point",
                hdf5::Array<1, double>("point", { dimension }),
                hdf5::Scalar<double>("log posterior"),
                hdf5::Scalar<double>("probability"),
                hdf5::Scalar<double>("cumulative"),
            };
        }

        GlobalLocal::JumpType
        GlobalLocal::jump_type(const unsigned & dimension)
        {
            return
            JumpType
            {
                "jump",
                hdf5::Array<1, double>("point", { dimension }),
                hdf5::Scalar<double>("log posterior"),
            };
        }

        GlobalLocal::LocalCovarianceType
        GlobalLocal::local_covariance_type(const unsigned & dimension)
        {
            return
            LocalCovarianceType
            {
                "local covariance",
                hdf5::Array<1, double>("mean", { dimension }),
                hdf5::Array<1, double>("covariance", { dimension * dimension })
            };
        }

        GlobalLocal::GlobalLocal(const std::vector<double> & component_probabilities,
                                 const unsigned & adaptations,
                                 const std::vector<MarkovChain::State> & jump_states,
                                 const std::vector<MarkovChain::State> & modes,
                                 const std::vector<std::vector<MarkovChain::State>> & history_states,
                                 const std::vector<std::vector<double>> & history_point_probabilities,
                                 const std::vector<std::vector<std::vector<double>>> & local_covariances,
                                 const std::vector<ProposalFunctionPtr> & proposals) :
            _adaptations(adaptations),
            _config(Config::Default()),
            _component_probabilities(component_probabilities),
            _history_points_local_covariance(local_covariances),
            _history_points_probabilities(history_point_probabilities),
            _history_states(history_states),
            _modes(modes),
            _prop(proposals)
        {
            // partial sum algorithm does exactly what we want: computes discrete cumulative
            std::partial_sum(_component_probabilities.cbegin(), _component_probabilities.cend(),
                             std::back_inserter(_component_cumulative));

            // rescale so probabilities sum up to exactly 1
            std::transform(_component_cumulative.begin(), _component_cumulative.end(), _component_cumulative.begin(),
                           std::bind(std::multiplies<double>(), 1.0 / _component_cumulative.back(), std::placeholders::_1));

            // copy jump vectors
            _jump_vectors.reserve(jump_states.size());
            for (auto s = jump_states.cbegin() ; s != jump_states.cend() ; ++s)
            {
                _jump_vectors.add(*s);
            }

            _sanity_check();
        }

        GlobalLocal::GlobalLocal(const std::vector<HistoryPtr> & chains,
                                 const std::vector<ProposalFunctionPtr> & proposals,
                                 const std::vector<MarkovChain::Stats> & stats,
                                 const Config & config,
                                 const unsigned & prerun_chains_per_partition) :
            _adaptations(0),
            _config(config)
        {
            Cluster::RValueFunction r = _config.clustering_strict_r_value ? &RValue::gelman_rubin : &RValue::approximation;

            unsigned chain_index = 0;
            std::list<HistoryPtr> available_chains(chains.begin(), chains.end());
            std::vector<Cluster> clusters{ Cluster(r, _config.clustering_maximum_r_value,
                                           available_chains.front(), chain_index, _config.skip_initial) };
            available_chains.pop_front();

            if (_config.perform_clustering)
            {
                Log::instance()->message("GlobalLocal.ctor", ll_informational)
                    << "Merging chains by comparing their R-values (max. allowed value: "
                    << _config.clustering_maximum_r_value << ", initial skip: " << _config.skip_initial << ")"
                    << ", " << (_config.clustering_strict_r_value ? "strict" : "relaxed") << " R-value definition";

                /* cluster chains according to their R-values */
                while (! available_chains.empty() && ++chain_index)
                {
                    Log::instance()->message("GL::ctor", ll_debug)
                        << "chain = " << chain_index;

                    // try to add a single chain to an existing cluster
                    bool added = false;

                    for (auto c = clusters.begin(), c_end = clusters.end() ; c != c_end ; ++c)
                    {
                        if (! c->overlaps(available_chains.front()))
                            continue;

                        c->add(available_chains.front(), chain_index);
                        added = true;

                        Log::instance()->message("GL::ctor", ll_debug)
                            << "Added chain " << chain_index << " to cluster " << std::distance(clusters.begin(), c);

                        break;
                    }

                    if (! added)
                    {
                        clusters.push_back(Cluster(r, _config.clustering_maximum_r_value,
                                            available_chains.front(), chain_index, _config.skip_initial));
                        Log::instance()->message("GL::ctor", ll_debug)
                            << "Created new cluster for chain " << chain_index;

                    }
                    available_chains.pop_front();
                }
            }
            // put each chain of a partition into the same cluster
            // assuming that chains are ordered according to partitions
            else
            {
                Log::instance()->message("GlobalLocal.ctor", ll_informational)
                    << "Merging all chains of a partition together, regardless of whether they fit together";
                // counter #chains in cluster
                unsigned counter = 0;
                while (! available_chains.empty() && ++chain_index && ++counter)
                {
                    if (counter < prerun_chains_per_partition)
                    {
                        // belongs to same cluster as previous chain
                        clusters.back().add(available_chains.front(), chain_index);
                    }
                    else
                    {
                        // we need a new cluster
                        clusters.push_back(Cluster(&RValue::approximation, _config.clustering_maximum_r_value,
                                            available_chains.front(), chain_index, _config.skip_initial));
                        counter = 0;
                    }
                    available_chains.pop_front();
                }
            }

            Log::instance()->message("global_local.ctor", ll_informational)
                << "Found " << clusters.size() << " clusters";

            /* combine proposal_functions from individual chains in cluster */

            for (auto cl = clusters.cbegin(), cl_end = clusters.cend() ; cl != cl_end ; ++cl)
            {
                Log::instance()->message("GlobalLocal.select", ll_debug)
                    << "Forming proposal for cluster " << std::distance(clusters.cbegin(), cl)
                    << " (" << std::distance(cl->begin(), cl->end()) << " chains) :";

                auto prop = proposals[*(cl->begin_indices())]->clone();

                if (config.join_chains_symmetrically)
                {
                    MultivariateProposalPtr mv = MultivariateAccess::access(prop);
                    if ( ! mv)
                        throw InternalError("GlobalLocal: Local proposal not of type Multivariate");

                    // find the covariance scale factor as an average of all chains in cluster
                    double mean_scale = 0.0;
                    std::vector<HistoryPtr> cluster_histories;
                    auto i = cl->begin_indices();
                    double k = 1;
                    for (auto c = cl->begin()++, c_end = cl->end() ; c != c_end ; ++c, ++i, ++k)
                    {
                        cluster_histories.push_back(*c);

                        MultivariateProposalPtr single_prop = MultivariateAccess::access(proposals[*i]);

                        double previous_mean_scale = mean_scale;
                        mean_scale += (single_prop->covariance_scale - previous_mean_scale) / k;

                        // check if scale factors are within two updates of each other, no warning in first step
                        const double ratio = (k > 1) ? single_prop->covariance_scale / previous_mean_scale : 1.0;
                        if (ratio < 1.0 / power_of<2>(Multivariate::covariance_scale_update_factor) ||
                            ratio > power_of<2>(Multivariate::covariance_scale_update_factor) )
                        {
                            Log::instance()->message("GlobalLocal.ctor", ll_warning)
                                << "Covariance scale factors vary significantly";
                        }

                        mean_scale /= _config.rescale_local_covariance;
                    }

                    mv->reset(cluster_histories, mean_scale, _config.skip_initial);
                }
                else
                {
                    MarkovChain::History history;

                    // copy only history of 2nd, 3rd... chain, to adapt proposal of first
                    for (auto c = cl->begin()++, c_end = cl->end() ; c != c_end ; ++c)
                    {
                        // copy log(posterior) and associate with index
                        if ((**c).states.empty())
                        {
                            throw InternalError("proposal_functions::GlobalLocal: cannot select points from empty history");
                        }

                        // Maybe it is better to skip the first, say 20% to remove outliers
                        MarkovChain::State::Iterator it_begin = (**c).states.cbegin() + unsigned(_config.skip_initial * (**c).states.size());
                        MarkovChain::State::Iterator it_end = (**c).states.cend();
                        std::copy(it_begin, it_end, std::back_inserter(history.states));
                    }

                    // hack: use efficiency which leaves scale factor unchanged
                    prop->adapt(history.states.cbegin(), history.states.cend(), 0.238, 0.2, 0.3);
                }

                // the local proposal for this cluster is ready
                _prop.push_back(prop);
            }

            /* determine component weights */

            // store maxima of each cluster
            std::vector<double> max_posterior(clusters.size(), -std::numeric_limits<double>::max());
            _modes.resize(clusters.size());

            {
                // cluster index
                unsigned i = 0;

                auto mp = max_posterior.begin();
                auto mode = _modes.begin();

                for (auto cl = clusters.cbegin(), cl_end = clusters.cend() ; cl != cl_end ; ++cl, ++i, ++mp, ++mode)
                {
                    std::vector<unsigned> chain_indices;
                    // find mode of chains within cluster
                    for (auto index = cl->begin_indices() ; index != cl->end_indices() ; ++index)
                    {
                        chain_indices.push_back(*index);

                        if (stats[*index].mode_of_posterior > *mp)
                        {
                            *mp = stats[*index].mode_of_posterior;
                            MarkovChain::State s;
                            s.log_posterior = stats[*index].mode_of_posterior;
                            s.point = stats[*index].parameters_at_mode;
                            s.hyper_parameter.component = i;
                            *mode = s;
                        }
                    }

                    Log::instance()->message("GlobalLocal.ctor", ll_debug)
                                << "Max posterior for cluster " << i << " = " << *mp << " at "
                                << stringify(mode->point.cbegin(), mode->point.cbegin() + std::min(3u, unsigned(mode->point.size()))) << '\n'
                                << "Chain indices are: " << stringify_container(chain_indices);
                }
            }

            if (_config.equal_weight_components)
            {
                // each cluster gets same weight
                _component_probabilities = std::vector<double>(clusters.size(), 1.0 / clusters.size());
                // compute cumulative
                _component_cumulative.clear();
                std::partial_sum(_component_probabilities.begin(), _component_probabilities.end(), std::back_inserter(_component_cumulative));
            }
            else
            {

                // go back from log scale and find relative weights by max of posterior
                const double global_maximum = *std::max_element(max_posterior.cbegin(), max_posterior.cend());
                std::vector<unsigned> negligible_cluster_indices;
                unsigned i = 0;
                for (auto c = clusters.begin(), c_end = clusters.end() ; c != c_end ; ++c, ++i)
                {
                    /*
                     * todo
                     * more accurate: Multiply with volume, or even compute the integral
                     * but for both we would need to know the parameter ranges (in chain.analysis)
                     */
                    const double relative_component_probability = std::max(std::exp(max_posterior[i] - global_maximum),
                                                                           double(_config.minimum_relative_cluster_weight));
                    if ( relative_component_probability < _config.minimum_relative_cluster_weight)
                    {
                        negligible_cluster_indices.push_back(i);
                        continue;
                    }
                    _component_probabilities.push_back(relative_component_probability);
                    _component_cumulative.push_back( _component_probabilities.back());
                    if (_component_cumulative.size() > 1)
                    {
                        _component_cumulative.back() += *(_component_cumulative.end() - 2);
                    }
                }

                // erase elements from behind, such that indices remain valid
                for (auto i = negligible_cluster_indices.rbegin(), i_end = negligible_cluster_indices.rend() ; i != i_end ; ++i)
                {
                    Log::instance()->message("GlobalLocal.ctor", ll_informational)
                        << "Removing cluster " << *i << " with relative weight " <<  *(_component_probabilities.begin() + *i);

                    clusters.erase(clusters.begin()+ *i);
                    _prop.erase(_prop.begin() + *i);
                    max_posterior.erase(max_posterior.begin() + *i);
                    _modes.erase(_modes.begin() + *i);
                }
            }

            // rescale so probabilities sum up to 1
            std::transform(_component_probabilities.begin(), _component_probabilities.end(), _component_probabilities.begin(),
                           std::bind(std::multiplies<double>(), 1.0 / _component_cumulative.back(), std::placeholders::_1));
            std::transform(_component_cumulative.begin(), _component_cumulative.end(), _component_cumulative.begin(),
                           std::bind(std::multiplies<double>(), 1.0 / _component_cumulative.back(), std::placeholders::_1));

            auto c = _component_cumulative.begin();
            for (auto p = _component_probabilities.begin(), p_end = _component_probabilities.end() ; p != p_end ; ++p, ++c)
            {
                Log::instance()->message("GlobalLocal.ctor", ll_debug)
                    << "comp prob = " << *p
                    << ", cum = " << *c;
            }

            _select_history_points(clusters);
            _select_jump_vectors(clusters);

            _sanity_check();

#if 0
            /* Find the weight of each cluster by integrating the posterior in the cluster volume */

            std::vector<MarkovChain> new_chains;
//            std::vector<std::tuple<double, double>> results;
            // total chain index
            unsigned i = 0;
            //todo replace with config option/argument
            static const unsigned number_of_iterations = 500;
            // extract proposal functions from  chains
            for (auto cl = clusters.cbegin(), cl_end = clusters.cend() ; cl != cl_end ; ++cl)
            {
                double numerator = 0.0;
                double denominator = 0.0;

                std::vector<double> theta_star = cl->begin()->statistics().parameters_at_mode;

                double chains_in_cluster = std::distance(cl->begin(), cl->end());

                for (auto ch = cl->begin(), ch_end = cl->end() ; ch != ch_end ; ++ch)
                {
                    new_chains.push_back(*ch);
                    std::tuple<double, double> result;
                    //todo parallelize
                    ch->normalized_density(result, theta_star, number_of_iterations);
                    numerator += std::get<0>(result);
                    denominator += std::get<1>(result);
                }

                //todo remove?
                // average over number of chains
                numerator /= chains_in_cluster;
                denominator /= chains_in_cluster;

                // divide by number of summands
                numerator /= double(cl->begin()->history().points.size());
                denominator /= double(number_of_iterations);

                double weight = std::exp(cl->begin()->statistics().mode_of_posterior) / numerator * denominator;

                Log::instance()->message("global_local.ctor", ll_debug)
                    << "Estimated weight of cluster " << std::distance(clusters.cbegin(), cl)
                    << " (" << unsigned(chains_in_cluster) << " chains) = " weight;

                // build cumulative with relative(!) weights.
                // Each chain in a cluster is assigned the same weight, its share of the weight
                for (auto ch = cl->begin(), ch_end = cl->end() ; ch != ch_end ; ++ch, ++i)
                {
                    _component_probabilities[i] = weight / chains_in_cluster;
                    _component_cumulative[i] = weight / chains_in_cluster;
                    if (i > 0)
                        _component_cumulative[i] += _component_cumulative[i - 1];

                }
            }
#endif
        }

        GlobalLocal::~GlobalLocal()
        {
        }

        void
        GlobalLocal::adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                           const double & /*efficiency*/, const double & /*efficiency_min*/, const double & /*efficiency_max*/)
        {
            //todo Decide whether we want to adapt or not if efficiency is OK
#if 0
            if(efficiency >= efficiency_min && efficiency <= efficiency_max)
                return;
#endif
            /*
             *  adjust only the component probabilities, but not the local proposal functions
             */
            const double n_bins = _component_probabilities.size();
            // - 0.5 just to have right bins, could use anything in ]0,1[
            auto hist = Histogram<1>::WithEqualBinning(0, n_bins - 0.5, n_bins);
            // compute component frequencies in history
            for (auto s = begin ; s != end ; ++s)
            {
                hist.insert(s->hyper_parameter.component);
            }

            // the first adaptation counts as 1, not 0
            ++_adaptations;

            //  \Sigma_n = (1 - 1/n^{cooling_power}) \Sigma_{n-1} +  1/n^{cooling_power} * S_n
            // add + 1 so the initial guess is not completely ignored
            double weight = 1.0 / std::pow(_adaptations + 1, _config.cooling_power);

            auto prob = _component_probabilities.begin();
            for (auto bin = hist.begin(), last_bin = hist.end() ; bin != last_bin ; ++bin, ++prob)
            {
                *prob = (1.0 - weight) * (*prob) + weight * bin->value / hist.entries();
            }
            // update cumulative as well
            _component_cumulative.clear();
            std::partial_sum(_component_probabilities.begin(), _component_probabilities.end(),
                             std::back_inserter(_component_cumulative));

            Log::instance()->message("GlobalLocal::adapt", ll_debug)
                << "New component probabilities: " << stringify_container(_component_probabilities);
        }

        ProposalFunctionPtr
        GlobalLocal::clone() const
        {
            return ProposalFunctionPtr(new GlobalLocal(*this));
        }

        void
        GlobalLocal::config(const Config & config)
        {
            _config = config;

            if ( ! _config.long_jump_indices.empty())
                _jump_vectors.indices(_config.long_jump_indices);
        }

        /*
         * Create the following structure in the base directory
         * ./components
         * ./history/0 ./history/1 ...
         * ./local proposals/0/... ./local proposals/1/... ...
         * ./meta
         */
        void
        GlobalLocal::dump_state(hdf5::File & file, const std::string & data_set_base_name) const
        {
            H5E_BEGIN_TRY
            {
                // store history and meta only once
                try
                {
                    // dimensionality of parameter space
                    const unsigned dimension = _jump_vectors.state(0).point.size();

                    // specify proposal type and dimension, so opening the rest is easier
                    auto meta_data_set = file.create_data_set(data_set_base_name + "/meta", meta_type());
                    auto meta_record = std::make_tuple("GlobalLocal", dimension);
                    meta_data_set << meta_record;

                    // add points for jump vectors
                    {
                        auto data_set_jump = file.create_data_set(data_set_base_name + "/jump", jump_type(dimension));

                        for (unsigned i = 0 ; i < _jump_vectors.number_of_clusters() ; ++i)
                        {
                            auto s = _jump_vectors.state(i);
                            data_set_jump << std::make_tuple(s.point, s.log_posterior);
                        }
                    }

                    // add modes
                    {
                        auto data_set_mode = file.create_data_set(data_set_base_name + "/modes", jump_type(dimension));

                        for (auto m = _modes.cbegin() ; m != _modes.cend() ; ++m)
                        {
                            data_set_mode << std::make_tuple(m->point, m->log_posterior);
                        }
                    }

                    /*
                     * Add one data set for each component's history points.
                     * NB: If this has not been initialized from the full chains' history,
                     * there will be nothing to dump.
                     */

//                    Log::instance()->message("Gl::dump", ll_debug)
//                                    << "prob.size = " << _history_points_probabilities.size()
//                                    << ", cum.size = " << _history_points_cumulatives.size()
//                                    << ", states.size = " << _history_states.size();

                    const unsigned & n = _history_points_cumulatives.size();
                    if (n > 0 && n == _history_points_probabilities.size() &&
                        n == _history_states.size() && _history_states.front().size() > 0)
                    {
                        Log::instance()->message("Gl::dump", ll_debug)
                            << "Dumping history points";

                        // loop over components
                        for (unsigned c = 0 ; c < _component_probabilities.size() ; ++c)
                        {
                            auto data_set = file.create_data_set(data_set_base_name + "/history/" + stringify(c), history_type(dimension));

                            // loop over points
                            auto prob = _history_points_probabilities[c].cbegin();
                            auto cum = _history_points_cumulatives[c].cbegin();
                            for (auto h = _history_states[c].cbegin(), h_end = _history_states[c].cend() ; h != h_end ; ++h, ++prob, ++cum)
                            {
                                data_set << std::make_tuple(h->point, h->log_posterior, *prob, *cum);
                            }
                        }
                    }

                    // one subgroup for each component's local proposal
                    // note that this will not fail if data sets exist already, so it shouldn't be the first in the try clause
                    unsigned i = 0;
                    for (auto p = _prop.begin(), p_end = _prop.end() ; p != p_end ; ++p, ++i)
                    {
                        (**p).dump_state(file, data_set_base_name + "/local proposals/" + stringify(i));
                    }

                    /* local covariances */

                    if (_config.history_points_local_covariance_size > 0 && ! _history_points_local_covariance.empty())
                    {
                        // loop over components
                        for (unsigned c = 0 ; c < _component_probabilities.size() ; ++c)
                        {
                            auto local_covariance_data_set = file.create_data_set(data_set_base_name + "/local covariances/" + stringify(c), local_covariance_type(dimension));

                            for (unsigned i = 0 ; i < _config.history_points ; ++i)
                            {
                                local_covariance_data_set << std::make_tuple(_history_states[c][i].point, _history_points_local_covariance[c][i]);
                            }
                        }
                    }
                }
                catch (HDF5Error &)
                {
                }
            }
            H5E_END_TRY;

            // comp prob as a data set
            auto comp_data_set = file.create_or_open_data_set(data_set_base_name + "/components", component_type(_component_probabilities.size()));
            comp_data_set << std::make_tuple(_adaptations, _component_probabilities);
        }

        bool nearly_equal(const double & a, const double & b ) {
          return std::fabs(a - b) < 1e-15;
        }

        double
        GlobalLocal::evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const
        {
            const unsigned & h_x = x.hyper_parameter.component;
            const unsigned & h_y = y.hyper_parameter.component;

            const double local = std::exp(_prop.at(h_y)->evaluate(x, y));

            // save some time and ignore non-local part
            if (_config.local_jump_probability == 1.0 || h_x == h_y) //todo doesn't work with overlapping clusters
            {
                return std::log(local);
            }

            double non_local = 0.0;

            std::vector<double> test_point(y.point.size());

            // Check only from y component to come to x
            auto jump = _jump_vectors.jump(h_y, h_x);

            if (h_y < h_x)
                std::transform(y.point.cbegin(), y.point.cend(), jump.cbegin(), test_point.begin(), std::plus<double>());
            else
                std::transform(y.point.cbegin(), y.point.cend(), jump.cbegin(), test_point.begin(), std::minus<double>());

            // Even if only a few dimensions are altered, this is correct
            auto diff = std::mismatch(x.point.cbegin(), x.point.cend(), test_point.cbegin(), &nearly_equal);

            // found a contribution?
            if (diff.second == test_point.cend())
            {
                non_local += _component_probabilities[h_x];
            }
            else
            {
                if (h_x != h_y)
                {
                    Log::instance()->message("prop::GL::evaluate", ll_debug)
                            << "Found mismatch at position" << std::distance(test_point.cbegin(), diff.second)
                            << " with values " << stringify(*(diff.first), 17)
                            << " and " << stringify(*(diff.second), 17);
                }
            }

            //todo check if ignoring local is correct! just a hot fix
            const double result =  (1 - _config.local_jump_probability) * non_local;

//            Log::instance()->message("prop::GL::evaluate", ll_debug)
//                << "Local = " << local << ", non-local = " << non_local;

            return std::log(result);
#if 0
            /* loop over all other clusters and check if they can contribute
             * perhaps only checking on y's component would be enough, but we want to be on the safe side
             */

            // smaller indices first
            for (unsigned i = 0 ; i < x.hyper_parameter.component ; ++i)
            {
                // add jump vector to y and see if we have a mismatch with x
                auto jump = _jump_vectors.jump(i, x.hyper_parameter.component);
                std::transform(y.point.cbegin(), y.point.cend(), jump.cbegin(), test_point.begin(), std::plus<double>());

                // Even if only a few dimensions are altered, this is correct
                auto result = std::mismatch(x.point.cbegin(), x.point.cend(), test_point.cbegin());

                // found a contribution?
                if (result.second == test_point.cend())
                {
                    non_local += _component_probabilities[i];
                }
            }

            // now larger indices
            for (unsigned i = x.hyper_parameter.component + 1 ; i < _prop.size() ; ++i)
            {
                // add jump vector to y and see if we have a mismatch with x
                auto jump = _jump_vectors.jump(x.hyper_parameter.component, i); // order of indices doesn't matter
                std::transform(y.point.cbegin(), y.point.cend(), jump.cbegin(), test_point.begin(), std::minus<double>());

                // Even if only a few dimensions are altered, this is correct
                auto result = std::mismatch(x.point.cbegin(), x.point.cend(), test_point.cbegin());

                // found a contribution?
                if (result.second == test_point.cend())
                {
                    non_local += _component_probabilities[i];
                }
            }
#endif
        }

    MarkovChain::State
    GlobalLocal::mode() const
    {
        if (_modes.empty())
            throw InternalError("prop::GlobalLocal::mode: modes uninitialized");

        /*
         * Implementation assumes that the modes have been used as fixpoints
         * for jump calculations;
         */

        try
        {
            // it is important to set the correct hyper parameter, else it could be out of range
            MarkovChain::State best_state = _modes.front();
            best_state.hyper_parameter.component = 0;

            // loop over clusters
            for (auto m = _modes.cbegin() + 1 ; m != _modes.cend() ; ++m)
            {
                if(m->log_posterior > best_state.log_posterior)
                {
                    best_state = *m;
                    best_state.hyper_parameter.component = std::distance(_modes.cbegin(), m);
                }
            }
            return best_state;
        }
        catch (InternalError & e)
        {
            throw InternalError("prop::GlobalLocal::mode: Couldn't check modes successfully");
        }
    }

    void
        GlobalLocal::propose(MarkovChain::State & proposal, const MarkovChain::State & current, gsl_rng * rng) const
        {
            double u = gsl_ran_flat(rng, 0.0, 1.0);

            // choose between a local and a non-local jump
            if (u < _config.local_jump_probability)
            {
                _prop.at(current.hyper_parameter.component)->propose(proposal, current, rng);
                proposal.hyper_parameter.component = current.hyper_parameter.component;

                return;
            }

            // choose (different!) component non-locally
            unsigned & h_prop = proposal.hyper_parameter.component;
            const unsigned & h_curr = current.hyper_parameter.component;
            do
            {
                h_prop = proposal_functions::random_index(_component_cumulative, rng);
            }
            while (h_prop == h_curr);

            // now add the jump, mind ordering of indices.
            //todo how to assign either std::plus or std::minus to a variable, what type is needed?
            if (h_prop > h_curr)
            {
                std::transform(current.point.cbegin(), current.point.cend(),
                        _jump_vectors.jump(h_curr, h_prop).cbegin(), proposal.point.begin(), std::plus<double>());
            }
            else
            {
                std::transform(current.point.cbegin(), current.point.cend(),
                        _jump_vectors.jump(h_curr, h_prop).cbegin(), proposal.point.begin(), std::minus<double>());
            }

//            Log::instance()->message("prop::GL::propose", ll_debug)
//                << "Proposing non-local jump from comp. " << current.hyper_parameter.component
//                << " to comp. " << proposal.hyper_parameter.component
//                << " with jump " << stringify_container(_jump_vectors.jump(h_curr, h_prop));
        }

        GlobalLocal::Config::Config() :
            clustering_maximum_r_value(1, std::numeric_limits<double>::max(), 1.1),
            clustering_strict_r_value(false),
            equal_weight_components(false),
            join_chains_symmetrically(true),
            history_points(10),
            history_points_local_covariance_size(0),
            history_points_ordered(true),
            history_point_weighting(Config::HistoryPointWeighting::posterior),
            minimum_relative_cluster_weight(0, 1, 1e-3),
            perform_clustering(false),
            rescale_local_covariance(0, std::numeric_limits<double>::max(), 1.0),
            skip_initial(0, 1, 0.1),
            cooling_power(0.5),
            local_jump_probability(0, 1, 0.5)
        {
        }

        GlobalLocal::Config
        GlobalLocal::Config::Default()
        {
            return GlobalLocal::Config();
        }

        // this functions expects the full covariance matrix in covariance prior to invocation
        void
        Multivariate::_compute_cholesky_and_inverse()
        {
            // copy _covariance matrix to _covariance_chol
            gsl_matrix_memcpy(_covariance_chol, _covariance);

            // calculate cholesky decomposition, needed for sampling and one step for inversion
            gsl_error_handler_t * default_gsl_error_handler = gsl_set_error_handler_off();
            if (GSL_EDOM == gsl_linalg_cholesky_decomp(_covariance_chol))
            {
                Log::instance()->message("prop::Multivariate.cholesky", ll_warning)
                    << "Covariance matrix is not positive definite!"
                    << "Proceed by setting off-diagonal elements to zero.";

                // _covariance_chol is potentially changed. Copy again
                gsl_matrix_memcpy(_covariance_chol, _covariance);

                // remove the off-diagonal elements of _covariance_chol
                for (unsigned i = 0 ; i < _dimension ; ++i)
                {
                    for (unsigned j = i + 1 ; j < _dimension ; ++j)
                    {
                        gsl_matrix_set(_covariance_chol, i, j, 0.0);
                        gsl_matrix_set(_covariance_chol, j, i, 0.0);
                    }
                }

                if (GSL_EDOM == gsl_linalg_cholesky_decomp(_covariance_chol))
                {
                    throw InternalError(
                         "prop::Multivariate: GSL couldn't find Cholesky decomposition of " + print_matrix(_covariance)
                        + "Apparently no moves were accepted, so try to increase number of iterations between updates "
                        + "or decrease initial proposal covariance. Proceed by taking square root of covariance manually");
                }
            }
            gsl_set_error_handler(default_gsl_error_handler);

            // copy cholesky decomposition to _covariance_inverse
            gsl_matrix_memcpy(_covariance_inverse, _covariance_chol);

            // calculate the inverse of _covariance
            gsl_linalg_cholesky_invert(_covariance_inverse);

            // remove the upper triangular part of _covariance_chol
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                for (unsigned j = i + 1 ; j < _dimension ; ++j)
                {
                    gsl_matrix_set(_covariance_chol, i, j, 0.0);
                }
            }

            // compute the normalization constant on log scale
            // log_det = 0.5 * ln(det(V))
            double log_det = 0.0;
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                log_det += std::log(gsl_matrix_get(_covariance_chol, i, i));
            }
            // -k/2 * log(2 Pi) - 1/2 log(abs(det(V)))
            // = -k2 * log(2 Pi) - log(det(L))
            norm = -0.5 * _dimension * std::log(2.0 * M_PI) - log_det;
        }

        void
        Multivariate::_copy(const Multivariate & other)
        {
            if (other._dimension != this->_dimension)
                throw InternalError("prop::Multivariate._copy: dimensions do not match ("
                    + stringify(this->_dimension) + " vs " + stringify(other._dimension) + ").");
            adaptations = other.adaptations;
            cooling_power = other.cooling_power;
            covariance_scale = other.covariance_scale;
            gsl_matrix_memcpy(_tmp_sample_covariance_current, other._tmp_sample_covariance_current);
            gsl_matrix_memcpy(_covariance, other._covariance);
            _index_list = other._index_list;
            _compute_cholesky_and_inverse();
        }

        Multivariate::CovarianceType
        Multivariate::covariance_type(const unsigned & dimension)
        {
            return CovarianceType {hdf5::Array<1, double>("covariance matrix", { dimension * dimension })};
        }

        Multivariate::Multivariate(const unsigned & dimension, const std::vector<double> & covariance, const bool & automatic_scaling) :
            _tmp_left(gsl_vector_alloc(dimension)),
            _tmp_right(gsl_vector_alloc(dimension)),
            _tmp_sample_covariance_current(gsl_matrix_alloc(dimension, dimension)),
            _covariance(gsl_matrix_alloc(dimension, dimension)),
            _covariance_inverse(gsl_matrix_alloc(dimension, dimension)),
            _covariance_chol(gsl_matrix_alloc(dimension, dimension)),
            _dimension(dimension),
            _index_list(dimension),
            adaptations(0),
            covariance_scale(2.38 * 2.38 / double(dimension)),
            cooling_power(0.5)
        {
            if (covariance.size() != dimension * dimension)
                throw InternalError("proposal_functions::MultivariateGaussian: covariance and dimension do not match");

            std::copy(covariance.cbegin(), covariance.cend(), this->_covariance->data);
            std::copy(covariance.cbegin(), covariance.cend(), this->_tmp_sample_covariance_current->data);

            // Why not use scale here? Then we would have to interpret arg covariance
            // as an estimate of the sample covariance.
            // The benefit would be that the scale is used in the very first proposals, thus
            // an update in the second step actually is meaningful
            if (automatic_scaling)
            {
                gsl_matrix_scale(_covariance, covariance_scale);
            }

            // basic checking: diagonals > 0
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                if (gsl_matrix_get(_covariance, i, i) <= 0)
                    throw InternalError("proposal_functions::MultivariateGaussian: diagonal covariance elements must be positive"
                            + print_matrix(_covariance));
            }

            _compute_cholesky_and_inverse();

            // indices from 0 to dimension - 1
            std::iota(_index_list.begin(), _index_list.end(), 0);
        }

        const double Multivariate::covariance_scale_min = 1e-4;
        const double Multivariate::covariance_scale_max = 100;
        const double Multivariate::covariance_scale_update_factor = 1.5;

        Multivariate::~Multivariate()
        {
            gsl_matrix_free(_covariance);
            gsl_matrix_free(_covariance_inverse);
            gsl_matrix_free(_covariance_chol);
            gsl_vector_free(_tmp_left);
            gsl_vector_free(_tmp_right);
            gsl_matrix_free(_tmp_sample_covariance_current);
        }

        void
        Multivariate::adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                            const double & efficiency, const double & efficiency_min, const double & efficiency_max)
        {
//            if (history.points.size() < 2)
            unsigned number_of_history_states = std::distance(begin, end);
            if (number_of_history_states < 2)
                throw InternalError("Multivariate: cannot estimate sample covariance for less than two points");


//             Log::instance()->message("prop::Multivariate", ll_debug)
//                 << "Previous sample covariance matrix :"
//                 << proposal_functions::print_matrix(_tmp_sample_covariance_current);

            // the first adaptation counts as 1, not 0
            ++adaptations;
            Log::instance()->message("prop::Multivariate", ll_debug)
                << "Adaptations: " << adaptations;

            // copy previous estimate. Avoid the zero matrix even in first adaptation
            gsl_matrix * tmp_sample_covariance_previous = gsl_matrix_alloc(_dimension, _dimension);
            gsl_matrix_memcpy(tmp_sample_covariance_previous, _tmp_sample_covariance_current);

            // set to zero and compute sample _covariance
            gsl_matrix_set_zero(_tmp_sample_covariance_current);

            // calculate mean in here, so we don't have to rely on the fact that mean from stats
            // is the mean of the last chunk, and not of all previous chunks
            std::vector<double> mean(_dimension, 0.0);
            for (auto s = begin ; s != end ; ++s)
            {
                for (auto i = _index_list.begin(), i_end = _index_list.end() ; i != i_end ; ++i)
                {
                    mean[*i] += s->point[*i];
                }
            }
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                mean[i] *= 1.0 / number_of_history_states;
            }

//            Log::instance()->message("prop::Multivariate", ll_debug)
//                << "mean = " << stringify(mean.begin(), mean.end());

            // _covariance calculation
            for (auto s = begin ; s != end ; ++s)
            {
                for (auto i = _index_list.begin(), i_end = _index_list.end() ; i != i_end ; ++i)
                {
                    // diagonal elements
                    _tmp_sample_covariance_current->data[*i + _dimension * *i] += power_of<2>(s->point[*i] - mean[*i]);

                    // off-diagonal elements
                    for (unsigned j = *i + 1 ; j < _dimension ; ++j)
                    {
                        double summand = (s->point[*i] - mean[*i]) * (s->point[j] -  mean[j]);
                        _tmp_sample_covariance_current->data[*i + _dimension * j] += summand;
                        _tmp_sample_covariance_current->data[j + _dimension * *i] += summand;
                    }
                }
            }

            // unbiased estimate. Enlarge scale.
            gsl_matrix_scale(_tmp_sample_covariance_current, 1.0 / (number_of_history_states - 1.0));

            //  \Sigma_n = (1 - 1/n^{cooling_power}) \Sigma_{n-1} +  1/n^{cooling_power} * S_n
            double weight = 1.0 / std::pow(adaptations + 1, cooling_power);
            gsl_matrix_scale(tmp_sample_covariance_previous, 1.0 - weight);
            gsl_matrix_scale(_tmp_sample_covariance_current, weight);
            gsl_matrix_add(_tmp_sample_covariance_current, tmp_sample_covariance_previous);

//             Log::instance()->message("prop::Multivariate", ll_debug)
//                 << "New sample covariance matrix :"
//                 << proposal_functions::print_matrix(_tmp_sample_covariance_current);

            gsl_matrix_free(tmp_sample_covariance_previous);

            double covariance_scale_old = covariance_scale;

            if (efficiency > efficiency_max)
            {
                if (covariance_scale < covariance_scale_max)
                {
                    covariance_scale *= covariance_scale_update_factor;
                }
            }
            else if (efficiency < efficiency_min)
            {
                if (covariance_scale > covariance_scale_min)
                {
                    covariance_scale /= covariance_scale_update_factor;
                }
            }

            if (covariance_scale > covariance_scale_max)
            {
                Log::instance()->message("prop::Multivariate.adapt", ll_warning)
                    << "Covariance scaling parameter (" << covariance_scale << ") exceeds sensible maximum of " << covariance_scale_max;
            }
            if (covariance_scale < covariance_scale_min)
            {
                Log::instance()->message("prop::Multivariate.adapt", ll_warning)
                    << "Covariance scaling parameter (" << covariance_scale << ") below sensible minimum of " << covariance_scale_min;
            }
            if (covariance_scale != covariance_scale_old)
            {
                Log::instance()->message("prop::Multivariate.adapt", ll_informational)
                   << "Change scale from " << covariance_scale_old << " to " << covariance_scale;
            }

            // proposal_covariance = 2.38^2 / dimension * sample_covariance
            gsl_matrix_memcpy(_covariance, _tmp_sample_covariance_current);
            gsl_matrix_scale(_covariance, covariance_scale);

            // recompute cholesky decomposition and inverse
            _compute_cholesky_and_inverse();

            // polymorphism!
            _compute_norm();
        }

        const gsl_matrix *
        Multivariate::covariance() const
        {
            return _covariance;
        }

        unsigned
        Multivariate::dimension() const
        {
            return _dimension;
        }

        void
        Multivariate::_dump_covariance(hdf5::File & file, const std::string & data_set_base_name,
                                            const std::string & proposal_type_name) const
        {
            /* store _covariance matrix */
            //todo for many parameters and small initial data set capacity, get seg fault in destructor. Why?
            {
                auto data_set = file.create_or_open_data_set(data_set_base_name + "/covariance", covariance_type(_dimension));
                std::vector<double> record(_dimension * _dimension);
                std::copy(_covariance->data, _covariance->data + _dimension * _dimension, record.begin());
                data_set << record;
            }
            H5E_BEGIN_TRY
            {
                // the meta set has only one line, so do nothing if it exists already
                try
                {
                    auto meta_data_set = file.create_data_set(data_set_base_name + "/meta", proposal_functions::meta_type());
                    auto meta_record = std::make_tuple(proposal_type_name.c_str(), _dimension);
                    meta_data_set << meta_record;
                }
                catch (HDF5Error &)
                {
                }
            }
            H5E_END_TRY;
        }

        void
        Multivariate::reset(const std::vector<HistoryPtr> & histories,
                                 const double & scale, const double & skip_initial)
        {
            // set to zero and compute sample _covariance
            gsl_matrix_set_zero(_tmp_sample_covariance_current);

            // calculate mean in here, but skip initial points of each history
            std::vector<double> mean(_dimension, 0.0);
            std::vector<unsigned> lengths;
            for (auto h = histories.cbegin(), h_end = histories.cend() ; h != h_end ; ++h)
            {
                const unsigned number_of_skipped_elements =  skip_initial * (**h).states.size();
                auto s = (**h).states.cbegin() + number_of_skipped_elements;
                lengths.push_back((**h).states.size() - number_of_skipped_elements);

                for (auto s_end = (**h).states.cend() ; s != s_end ; ++s)
                {
                    for (auto i = _index_list.begin(), i_end = _index_list.end() ; i != i_end ; ++i)
                    {
                        mean[*i] += s->point[*i];
                    }
                }
            }

            // rescale the mean
            const unsigned total_length = std::accumulate(lengths.cbegin(), lengths.cend(), 0);
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                mean[i] *= 1.0 / total_length;
            }

            Log::instance()->message("prop::Multivariate::reset", ll_debug)
                << "mean = " << stringify(mean.begin(), mean.end());

            // _covariance calculation
            auto l = lengths.cbegin();
            for (auto h = histories.cbegin(), h_end = histories.cend() ; h != h_end ; ++h, ++l)
            {
                // count from the back
                auto s_end = (**h).states.cend();
                auto s = s_end - (*l);
                for ( ; s != s_end ; ++s)
                {
                    for (auto i = _index_list.begin(), i_end = _index_list.end() ; i != i_end ; ++i)
                    {
                        // diagonal elements
                        _tmp_sample_covariance_current->data[*i + _dimension * *i] += power_of<2>(s->point[*i] - mean[*i]);

                        // off-diagonal elements
                        for (unsigned j = *i + 1 ; j < _dimension ; ++j)
                        {
                            double summand = (s->point[*i] - mean[*i]) * (s->point[j] -  mean[j]);
                            _tmp_sample_covariance_current->data[*i + _dimension * j] += summand;
                            _tmp_sample_covariance_current->data[j + _dimension * *i] += summand;
                        }
                    }
                }
            }

            // unbiased estimate. Enlarge scale.
            gsl_matrix_scale(_tmp_sample_covariance_current, 1.0 / (total_length - 1.0));

            // Ignore the usual update formula
            adaptations = 0;

            if (scale > covariance_scale_max)
            {
                Log::instance()->message("prop::Multivariate.reset", ll_warning)
                    << "Hit maximum of covariance scaling parameter!";
            }
            if (scale < covariance_scale_min)
            {
                Log::instance()->message("prop::Multivariate.reset", ll_warning)
                    << "Hit minimum of covariance scaling parameter!";
            }

            // proposal_covariance = scale * sample_covariance
            gsl_matrix_memcpy(_covariance, _tmp_sample_covariance_current);
            gsl_matrix_scale(_covariance, scale);

//            Log::instance()->message("prop::Multivariate::reset", ll_debug)
//                 << "New  covariance matrix :"
//                 << proposal_functions::print_matrix(_covariance);

            // recompute cholesky decomposition and inverse
            _compute_cholesky_and_inverse();

            // polymorphism!
            _compute_norm();
        }

        void
        Multivariate::rescale(const double & rescale_factor)
        {
            // first divide, then multiply
            gsl_matrix_scale(_covariance, 1 / covariance_scale);
            covariance_scale *= rescale_factor;
            gsl_matrix_scale(_covariance, covariance_scale);

            _compute_cholesky_and_inverse();
            _compute_norm();
        }

        void
        Multivariate::set_indices(const std::vector<unsigned> & index_list)
        {
            if (index_list.size() != _dimension)
                throw InternalError("Multivariate::set_indices: dimension mismatch between dimension (" + stringify(_dimension)
                    + ") and index_list.size (" + stringify(index_list.size()) + ")");

            _index_list = index_list;
        }

        void
        MultivariateGaussian::_compute_norm()
        {
            // compute the normalization constant on log scale
            // log_det = 0.5 * ln(det(V))
            double log_det = 0.0;
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                log_det += std::log(gsl_matrix_get(_covariance_chol, i, i));
            }
            // -k/2 * log(2 Pi) - 1/2 log(abs(det(V)))
            // = -k2 * log(2 Pi) - log(det(L))
            norm = -0.5 * _dimension * std::log(2.0 * M_PI) - log_det;
        }

        MultivariateGaussian::MultivariateGaussian(const unsigned & dimension,
                                                   const std::vector<double> & covariance,
                                                   const bool & automatic_scaling) :
            Multivariate(dimension, covariance, automatic_scaling)
        {
            _compute_norm();
        }

        MultivariateGaussian::ScalarsType
        MultivariateGaussian::scalars_type()
        {
            return
            MultivariateGaussian::ScalarsType
            {
                "single numbers",
                hdf5::Scalar<double>("covariance scale"),
                hdf5::Scalar<double>("cooling factor"),
                hdf5::Scalar<unsigned>("adaptations"),
            };
        }

        MultivariateGaussian::~MultivariateGaussian()
        {
        }

        ProposalFunctionPtr
        MultivariateGaussian::clone() const
        {
            std::vector<double> cov(_covariance->data, _covariance->data + _dimension * _dimension);
            MultivariateGaussian * mvg = new MultivariateGaussian(_dimension, cov);
            mvg->_copy(*this);
            mvg->_compute_norm();

            return ProposalFunctionPtr(mvg);
        }

        void
        MultivariateGaussian::dump_state(hdf5::File & file, const std::string & data_set_base_name) const
        {
            _dump_covariance(file, data_set_base_name, "MultivariateGaussian");

            auto data_set = file.create_or_open_data_set(data_set_base_name + "/scalars", scalars_type());
            auto record = std::make_tuple(covariance_scale, cooling_power, adaptations);
            data_set << record;
        }

        double
        MultivariateGaussian::evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const
        {
            double chi_squared = 0.0;

            std::copy(x.point.cbegin(), x.point.cend(), _tmp_left->data);
            std::copy(y.point.cbegin(), y.point.cend(), _tmp_right->data);
            gsl_vector_sub(_tmp_left, _tmp_right);

            gsl_blas_dgemv(CblasNoTrans, 1.0, _covariance_inverse, _tmp_left, 0.0, _tmp_right);
            gsl_blas_ddot(_tmp_left, _tmp_right, &chi_squared);

//            Log::instance()->message("MultivariateGaussian::evaluate", ll_debug)
//                << "chi^2 = " << chi_squared
//                << ", norm = " << norm;

            return norm - chi_squared / 2.0;
        }

        void
        MultivariateGaussian::propose(MarkovChain::State & proposal, const MarkovChain::State & current, gsl_rng * rng) const
        {
            // generate standard normals
            std::generate(_tmp_left->data, _tmp_left->data + _dimension, std::bind(gsl_ran_ugaussian, rng));

            // transform
            gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, _covariance_chol, _tmp_left);

            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                proposal.point[i] = current.point[i] + gsl_vector_get(_tmp_left, i);
            }
        }

        void
        MultivariateStudentT::_compute_norm()
        {
            // compute the normalization constant on log scale
            // log_det = 0.5 * ln(det(V))
            double log_det = 0.0;
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                log_det += std::log(gsl_matrix_get(_covariance_chol, i, i));
            }
            // see http://en.wikipedia.org/wiki/Multivariate_Student_distribution
            norm = gsl_sf_lngamma(0.5 * (dof + _dimension)) - gsl_sf_lngamma(0.5 * dof)
                   - 0.5 * _dimension * std::log(dof * M_PI) - log_det;
        }

        MultivariateStudentT::MultivariateStudentT(const unsigned int & dimension,
                                                   const std::vector<double> & covariance,
                                                   const double & degree_of_freedom,
                                                   const bool & automatic_scaling) :
            Multivariate(dimension, covariance, automatic_scaling),
            dof(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max(), degree_of_freedom)
        {
            _compute_norm();
        }

        MultivariateStudentT::ScalarsType
        MultivariateStudentT::scalars_type()
        {
            return MultivariateStudentT::ScalarsType
                   {
                        "single numbers",
                        hdf5::Scalar<double>("covariance scale"),
                        hdf5::Scalar<double>("cooling factor"),
                        hdf5::Scalar<unsigned>("adaptations"),
                        hdf5::Scalar<double>("degrees of freedom"),
                   };
        }

        MultivariateStudentT::~MultivariateStudentT()
        {
        }

        ProposalFunctionPtr
        MultivariateStudentT::clone() const
        {
            std::vector<double> cov(_covariance->data, _covariance->data + _dimension * _dimension);
            MultivariateStudentT * mvt = new MultivariateStudentT(_dimension, cov, dof);

            mvt->_copy(*this);
            mvt->_compute_norm();

            return ProposalFunctionPtr(mvt);
        }

        void
        MultivariateStudentT::dump_state(hdf5::File & file, const std::string & data_set_base_name) const
        {
            _dump_covariance(file, data_set_base_name, "MultivariateStudentT");

            auto data_set = file.create_or_open_data_set(data_set_base_name + "/scalars", scalars_type());
            auto record = std::make_tuple(covariance_scale, cooling_power, adaptations, double(dof));
            data_set << record;
        }

        double
        MultivariateStudentT::evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const
        {
            double chi_squared = 0.0;

            // center around zero
            std::copy(x.point.cbegin(), x.point.cend(), _tmp_left->data);
            std::copy(y.point.cbegin(), y.point.cend(), _tmp_right->data);
            gsl_vector_sub(_tmp_left, _tmp_right);

            // \chi^2 from bilinear form
            gsl_blas_dgemv(CblasNoTrans, 1.0, _covariance_inverse, _tmp_left, 0.0, _tmp_right);
            gsl_blas_ddot(_tmp_left, _tmp_right, &chi_squared);

            return norm - 0.5 * (dof + _dimension) * std::log(1.0 + chi_squared / dof);
        }

        void
        MultivariateStudentT::propose(MarkovChain::State & proposal, const MarkovChain::State & current, gsl_rng * rng) const
        {
            // generate standard normals
            std::generate(_tmp_left->data, _tmp_left->data + _dimension, std::bind(gsl_ran_ugaussian, rng));

            // transform to N(0, Sigma)
            gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, _covariance_chol, _tmp_left);

            // correct for degrees of freedom
            gsl_vector_scale(_tmp_left, std::sqrt(dof / gsl_ran_chisq(rng, dof)));

            // add mean
            for (unsigned i = 0 ; i < _dimension ; ++i)
            {
                proposal.point[i] = current.point[i] + gsl_vector_get(_tmp_left, i);
            }
        }

        BlockDecomposition::PriorsType
        BlockDecomposition::priors_type()
        {
            return PriorsType{ "prior", hdf5::Scalar<const char *>("prior description") };
        }

        void
        BlockDecomposition::_copy_values(const std::vector<std::shared_ptr<double>> & ptr_vector, std::vector<double> & result)
        {
            result.resize(ptr_vector.size());
            auto res = result.begin();
            for (auto ptr = ptr_vector.cbegin(), ptr_end = ptr_vector.cend() ; ptr != ptr_end ; ++ptr, ++res)
            {
                *res = **ptr;
            }
        }

        void
        BlockDecomposition::_copy_values(const std::vector<double> & source, const std::vector<std::shared_ptr<double>> & destination_ptr_vector)
        {
            if (source.size() != destination_ptr_vector.size())
                throw InternalError("BlockDecomposition::copy_values: size mismatch");

            auto dest_ptr = destination_ptr_vector.begin();
            for (auto i = source.cbegin(), i_end = source.cend() ; i != i_end ; ++i, ++dest_ptr)
            {
                **dest_ptr = *i;
            }
        }

        BlockDecomposition::BlockDecomposition() :
            p(Parameters::Defaults())
        {
        }

        BlockDecomposition::~BlockDecomposition()
        {
        }

        void
        BlockDecomposition::adapt(const MarkovChain::State::Iterator & begin, const MarkovChain::State::Iterator & end,
                           const double & efficiency, const double & efficiency_min, const double & efficiency_max)
        {
            // adapt only the multivariate components
            for (auto mv = _mv.cbegin(), mv_end = _mv.cend() ; mv != mv_end ; ++mv)
            {
                (**mv).adapt(begin, end, efficiency, efficiency_min, efficiency_max);
            }
        }

        void
        BlockDecomposition::add(const LogPriorPtr & prior)
        {
            _priors.push_back(prior->clone(p));

            // ready for multidimensional priors
            std::vector<std::shared_ptr<double>> pointers_to_parameter_values;
            for (auto d = prior->begin(), d_end = prior->end() ; d != d_end ; ++d)
            {
                auto x = std::make_shared<double>(1);
                _tmp_vector_x.push_back(x);
                pointers_to_parameter_values.push_back(x);
            }
            _priors_values_x.push_back(pointers_to_parameter_values);

            pointers_to_parameter_values.clear();

            for (auto d = prior->begin(), d_end = prior->end() ; d != d_end ; ++d)
            {
                auto y = std::make_shared<double>(1);
                _tmp_vector_y.push_back(y);
                pointers_to_parameter_values.push_back(y);
            }
            _priors_values_y.push_back(pointers_to_parameter_values);
        }

        void
        BlockDecomposition::add(const MultivariateProposalPtr & mv)
        {
            // check that only one MultivariateProposalPtr is added
            if ( ! _mv.empty())
                throw InternalError("BlockDecomposition::add: At the moment, only one multivariate proposal is supported");

            // toggle between the different base class pointers
            ProposalFunctionPtr mv_clone = mv->clone();
            _mv.push_back(std::static_pointer_cast<Multivariate>(mv_clone));

            // loop over dimensions
            unsigned dimension = _mv.back()->dimension();

            std::vector<std::shared_ptr<double>> pointers_to_parameter_values;
            std::vector<unsigned> index_list;
            for (unsigned i = 0 ; i < dimension ; ++i)
            {
                // get size before element added to vector
                index_list.push_back(_tmp_vector_x.size());
                auto x = std::make_shared<double>(1);
                _tmp_vector_x.push_back(x);
                pointers_to_parameter_values.push_back(x);
            }
            _mv_values_x.push_back(pointers_to_parameter_values);
            _mv.back()->set_indices(index_list);
            _tmp_state_x.point.resize(dimension);

            index_list.clear();
            pointers_to_parameter_values.clear();
            for (unsigned i = 0 ; i < dimension ; ++i)
            {
                index_list.push_back(_tmp_vector_y.size());
                auto y = std::make_shared<double>(1);
                _tmp_vector_y.push_back(y);
                pointers_to_parameter_values.push_back(y);
            }
            _mv_values_y.push_back(pointers_to_parameter_values);
            _tmp_state_y.point.resize(dimension);
        }

        ProposalFunctionPtr
        BlockDecomposition::clone() const
        {
            BlockDecomposition * bd = new BlockDecomposition();

            // add multivariates
            for (auto mv = _mv.cbegin(), mv_end = _mv.cend() ; mv != mv_end ; ++mv)
            {
                bd->add(*mv);
            }

            // add priors
            for (auto prior = _priors.cbegin(), prior_end = _priors.cend() ; prior != prior_end ; ++prior)
            {
                bd->add(*prior);
            }

            return ProposalFunctionPtr(bd);
        }

        void
        BlockDecomposition::dump_state(hdf5::File & file, const std::string & data_set_base_name) const
        {
            try
            {
                H5E_BEGIN_TRY
                {
                    auto meta_data_set = file.create_data_set(data_set_base_name + "/meta", proposal_functions::meta_type());
                    auto meta_record = std::make_tuple("BlockDecomposition", _tmp_vector_x.size());
                    meta_data_set << meta_record;

                    // one data set for all priors as they are serialized
                    auto record = std::make_tuple("prior");
                    auto prior_data_set = file.create_or_open_data_set(data_set_base_name + "/priors", priors_type());
                    for (auto prior = _priors.begin(), prior_end = _priors.end() ; prior != prior_end ; ++prior)
                    {
                        std::string serialization((**prior).as_string());
                        std::get<0>(record) = serialization.c_str();
                        prior_data_set << record;
                    }
                }
                H5E_END_TRY;
            }
            catch (HDF5Error &)
            {
            }

            // one subgroup for each multivariate component
            {
                unsigned i = 0;
                for (auto mv = _mv.cbegin(), mv_end = _mv.cend() ; mv != mv_end ; ++mv, ++i)
                {
                    (**mv).dump_state(file, data_set_base_name + "/multivariates/" + stringify(i));
                }
            }
        }

        double
        BlockDecomposition::evaluate(const MarkovChain::State & x, const MarkovChain::State & y) const
        {
            //assign
            for (unsigned i = 0 ; i < y.point.size() ; ++i)
            {
                *_tmp_vector_x.at(i) = x.point[i];
                *_tmp_vector_y.at(i) = y.point[i];
            }
            double result = 0.0;

            // multivariate contribution
            auto mv_values_x = _mv_values_x.cbegin();
            auto mv_values_y = _mv_values_y.cbegin();
            for (auto mv = _mv.cbegin(), mv_end = _mv.cend() ; mv != mv_end ; ++mv, ++mv_values_x, ++mv_values_y)
            {
                _copy_values(*mv_values_x, _tmp_state_x.point);

                _copy_values(*mv_values_y, _tmp_state_y.point);

                result += (**mv).evaluate(_tmp_state_x, _tmp_state_y);
            }

            // prior contribution is non-local, i.e. independent of y
            auto prior_values = _priors_values_x.cbegin();
            for (auto prior = _priors.begin(), prior_end = _priors.end() ; prior != prior_end ; ++prior, ++prior_values)
            {
                // loop over parameters in the prior (usually just a single one)
                // assign to Parameters object
                auto par = prior_values->cbegin();
                for (auto d = (**prior).begin(), d_end = (**prior).end() ; d != d_end ; ++d, ++par)
                {
                    d->parameter = **par;
                }

                result += (**prior)();
            }

            return result;
        }

        void
        BlockDecomposition::propose(MarkovChain::State & x, const MarkovChain::State & y, gsl_rng * rng) const
        {
            //assign
            for (unsigned i = 0 ; i < y.point.size() ; ++i)
            {
                *_tmp_vector_y.at(i) = y.point[i];
            }

            // multivariate part
            auto mv_values_x = _mv_values_x.begin();
            auto mv_values_y = _mv_values_y.cbegin();
            for (auto mv = _mv.cbegin(), mv_end = _mv.cend() ; mv != mv_end ; ++mv, ++mv_values_x, ++mv_values_y)
            {
                _copy_values(*mv_values_y, _tmp_state_y.point);

                // get proposal into _tmp_state_x
                (**mv).propose(_tmp_state_x, _tmp_state_y, rng);

                _copy_values(_tmp_state_x.point, *mv_values_x);
            }

            // prior part: works only with 1D priors
            auto prior_values = _priors_values_x.cbegin();
            for (auto prior = _priors.begin(), prior_end = _priors.end() ; prior != prior_end ; ++prior, ++prior_values)
            {
                *(prior_values->front()) = (**prior).sample(rng);
            }

            // copy final result now that all values have been updated
            _copy_values(_tmp_vector_x, x.point);
        }
    }
}
