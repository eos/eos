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

#include <eos/statistics/proposal-functions.hh>
#include <eos/statistics/histogram.hh>
#include <eos/statistics/log-prior.hh>
#include <eos/statistics/rvalue.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>

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
                    d->parameter->set(**par);
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
