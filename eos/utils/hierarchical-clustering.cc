/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Frederik Beaujean
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

#include <eos/utils/hierarchical-clustering.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <cmath>
#include <limits>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace eos
{
    // TODO: replace by lambda function
    bool is_dead(const HierarchicalClustering::Component & c)
    {
        return c.weight() == 0.0;
    }

    template <> struct Implementation<HierarchicalClustering>
    {
        HierarchicalClustering::Config config;

        // The smallest components
        HierarchicalClustering::MixtureDensity components;

        // The clusters created from components
        HierarchicalClustering::MixtureDensity clusters;

        // keep track of \pi(i), i=1..k
        std::vector<unsigned> mapping;

        // which components map to a given cluster, \pi^{-1}(j)
        std::vector<std::vector<unsigned>> inverse_mapping;

        // matrix, the (i,j) element gives KL(f_i || g_j)
        std::vector<double> divergences;

        Implementation(const HierarchicalClustering::Config & config) :
            config(config)
        {
        }

        ~Implementation()
        {
        }

        void add(const HierarchicalClustering::Component & component)
        {
            components.push_back(component);
        }

        /*
         * Look for dead clusters (weight=0) and remove them.
         * Resize storage.
         */
        void cleanup()
        {
            unsigned active_clusters = clusters.size();
            if (config.kill_clusters)
            {
                auto new_end = std::remove_if(clusters.begin(), clusters.end(), &is_dead);
                active_clusters = std::distance(clusters.begin(), new_end);
                clusters.erase(new_end, clusters.end());
            }
            inverse_mapping.resize(active_clusters);
            divergences.resize(active_clusters * components.size());
        }

        void compute_KL()
        {
            for (unsigned i = 0 ; i < components.size() ; ++i)
            {
                for (unsigned j = 0 ; j < clusters.size() ; ++j)
                {
                    divergences[i * clusters.size() + j] = kullback_leibler_divergence(components[i], clusters[j]);
                }
            }
        }

        /*
         * Compute the distance function d(f,g,\pi), Eq. (3)
         */
        double distance()
        {
            double result = 0.0;
            for (unsigned i = 0 ; i < components.size() ; ++i)
            {
                result += components[i].weight() * divergences[i * clusters.size() + mapping[i]];
            }

            return result;
        }

        /*
         *  KL(c1 || c2) mind the ordering!
         *  Use same notation as in Goldberger, Roweis, ch. 2.
         *
         *  @note: KL(1 || 2) >= 0, and KL(1 || 1) = 0
         */
        static double kullback_leibler_divergence(const HierarchicalClustering::Component & c1, const HierarchicalClustering::Component & c2)
        {
            // first contribution: ratio of determinants
            double d = std::log(c2.determinant() / c1.determinant());
            const unsigned dim = c1.mean()->size;
            gsl_matrix * matrix_res = gsl_matrix_alloc(dim, dim);

            // second contribution: trace of product
            gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, c2.inverse_covariance(), c1.covariance(), 0.0, matrix_res);
            for (unsigned i = 0 ; i < dim ; ++i)
            {
                d += gsl_matrix_get(matrix_res, i, i);
            }
            gsl_matrix_free(matrix_res);

            // third contribution: \chi^2
            gsl_vector * vector_res_left = gsl_vector_alloc(dim);
            gsl_vector * vector_res_right = gsl_vector_alloc(dim);
            gsl_vector_memcpy(vector_res_left, c1.mean());
            gsl_vector_sub(vector_res_left, c2.mean());
            gsl_vector_memcpy(vector_res_right, vector_res_left);
            double chi_squared = 0;

            // a = sigma_2^{-1} * (mu_1 - mu_2)
            gsl_blas_dsymv(CblasUpper, 1.0, c2.inverse_covariance(), vector_res_left, 0.0, vector_res_right);

            // chi^2 = (mu_1 - mu_2) * a
            gsl_blas_ddot(vector_res_left, vector_res_right, &chi_squared);
            gsl_vector_free(vector_res_left);
            gsl_vector_free(vector_res_right);

            d += chi_squared;

            // fourth contribution: dimension
            d -= dim;

            d *= 0.5;

            return d;
        }

        // Eq. (7) and below
        void refit()
        {
            // need temporary vector/matrix for addition
            gsl_vector * mu_diff = gsl_vector_alloc(components[0].mean()->size);
            gsl_matrix * sigma = gsl_matrix_alloc(components[0].mean()->size, components[0].mean()->size);

            for (unsigned j = 0 ; j < clusters.size() ; ++j)
            {
                // initialize values
                clusters[j].weight() = 0;
                gsl_vector_set_all(clusters[j].mean(), 0);
                gsl_matrix_set_all(clusters[j].covariance(), 0);
                gsl_matrix_set_all(sigma, 0);

                // compute total weight and mean
                for (auto i = inverse_mapping[j].cbegin() ; i != inverse_mapping[j].cend() ; ++i)
                {
                    // update weight
                    clusters[j].weight() += components[*i].weight();
                    // update mean
                    gsl_blas_daxpy(components[*i].weight(), components[*i].mean(), clusters[j].mean());
                }
                // rescale by total weight
                gsl_vector_scale(clusters[j].mean(), 1.0 / clusters[j].weight());

                // update covariance
                for (auto i = inverse_mapping[j].cbegin() ; i != inverse_mapping[j].cend() ; ++i)
                {
                    // mu_diff = mu'_j - mu_i
                    gsl_vector_memcpy(mu_diff, clusters[j].mean());
                    gsl_vector_sub(mu_diff, components[*i].mean());
                    // sigma = (mu'_j - mu_i) (mu'_j - mu_i)^T
                    gsl_blas_dsyr(CblasUpper, 1.0, mu_diff, sigma);
                    // sigma += sigma_i
                    gsl_matrix_add(sigma, components[*i].covariance());
                    // multiply with alpha_i
                    gsl_matrix_scale(sigma, components[*i].weight());
                    // sigma_j += alpha_i * (sigma_i + (mu'_j - mu_i) (mu'_j - mu_i)^T
                    gsl_matrix_add(clusters[j].covariance(), sigma);
                }
                // 1 / beta_j
                gsl_matrix_scale(clusters[j].covariance(), 1.0 / clusters[j].weight());
            }
            gsl_matrix_free(sigma);
            gsl_vector_free(mu_diff);

            // reset the components, s.t. inverse and determinant are recomputed
        }

        // Eq. (4)
        void regroup()
        {
            // clean up old maps
            for (unsigned j = 0 ; j < clusters.size() ; ++j)
            {
                inverse_mapping[j].clear();
            }

            // find smallest divergence between component i and any cluster j of the cluster mixture density
            for (unsigned i = 0 ; i < components.size() ; ++i)
            {
                auto first = divergences.cbegin() + i * clusters.size();
                auto last =  divergences.cbegin() + (i + 1) * clusters.size();
                unsigned j = std::distance(first, std::min_element(first, last));
                mapping[i] = j;
                inverse_mapping[j].push_back(i);
            }
#if 0
            std::string output = "New mapping:\n";

            unsigned j = 0;
            for (auto i = inverse_mapping.cbegin() ; i != inverse_mapping.cend() ; ++i, ++j)
            {
                output += stringify(j) + ": ";
                output += print_container(*i);
                output += '\n';
            }
            Log::instance()->message("HierarchicalClustering::regroup", ll_debug)
                << output;
#endif
        }

        void run()
        {
            // sanity checks and setup
            {
                if (this->clusters.empty())
                    throw InternalError("HierarchicalClustering::run: initial guess required");

                if (components.empty())
                    throw InternalError("HierarchicalClustering::run: no components specified");

                if (components.size() <= clusters.size())
                    throw InternalError("HierarchicalClustering::run: cannot reduce #components");

                mapping.resize(components.size());

                if (config.equal_weights)
                {
                    for (unsigned i = 0 ; i < components.size() ; ++i)
                    {
                        components[i].weight() = 1.0 / components.size();
                    }
                }
            }

            double old_distance = std::numeric_limits<double>::max();
            double new_distance = std::numeric_limits<double>::max();
            bool converged = false;
            unsigned step = 0;

            while ((! converged) && (step < config.maximum_steps))
            {
                cleanup();

                // regroup
                compute_KL();
                regroup();

                // refit
                refit();

                // check convergence
                new_distance = distance();

                Log::instance()->message("HierarchicalClustering::run", ll_debug)
                    << "Current distance in step " << step << ": " << stringify(new_distance, 16);

                if (new_distance == old_distance)
                {
                    converged = true;
                    Log::instance()->message("HierarchicalClustering::run", ll_informational)
                        << "Found exact local minimum after " << step << " steps ";
                }

                const double rel_distance_change = (old_distance - new_distance) / old_distance;

                if (new_distance < 0)
                    throw InternalError("HierarchicalClustering::run: found negative distance");

                if ((new_distance - old_distance) / old_distance > 1e-13)
                    throw InternalError("HierarchicalClustering::run: distance increased");

                if (rel_distance_change < config.precision && not converged && step > 0)
                {
                    converged = true;
                    Log::instance()->message("HierarchicalClustering::run", ll_informational)
                        << "Close enough to local minimum after " << step << " steps ";
                }

                old_distance = new_distance;
                step++;
            }
        }
    };

    HierarchicalClustering::HierarchicalClustering(const HierarchicalClustering::Config & config) :
        PrivateImplementationPattern<HierarchicalClustering>(new Implementation<HierarchicalClustering>(config))
    {
    }

    HierarchicalClustering::~HierarchicalClustering()
    {
    }

    void
    HierarchicalClustering::add(const Component & component)
    {
        _imp->add(component);
    }

    void
    HierarchicalClustering::initial_guess(const MixtureDensity & density)
    {
        _imp->clusters = density;
        _imp->inverse_mapping.resize(density.size());

        double total_weight = 0;
        for (auto comp = density.cbegin() ; comp != density.cend() ; ++comp)
        {
            total_weight += comp->weight();
        }

        if (std::abs(total_weight - 1.0) > 1e-8)
            throw InternalError("HierarchicalClustering::initial_guess: Weights are not normalized. Got total weight of " + stringify(total_weight, 16));
    }

    void
    HierarchicalClustering::run()
    {
        _imp->run();
    }

    HierarchicalClustering::ComponentIterator
    HierarchicalClustering::begin_components() const
    {
        return ComponentIterator(_imp->components.begin());
    }

    HierarchicalClustering::ComponentIterator
    HierarchicalClustering::end_components() const
    {
        return ComponentIterator(_imp->components.end());
    }

    HierarchicalClustering::ClusterIterator
    HierarchicalClustering::begin_clusters() const
    {
        return ClusterIterator(_imp->clusters.begin());
    }

    HierarchicalClustering::ClusterIterator
    HierarchicalClustering::end_clusters() const
    {
        return ClusterIterator(_imp->clusters.end());
    }

    HierarchicalClustering::MapIterator
    HierarchicalClustering::begin_map() const
    {
        return MapIterator(_imp->mapping.begin());
    }

    HierarchicalClustering::MapIterator
    HierarchicalClustering::end_map() const
    {
        return MapIterator(_imp->mapping.end());
    }

    template <> struct Implementation<HierarchicalClustering::Component>
    {
        unsigned dimension;
        gsl_matrix * covariance;
        gsl_matrix * inverse_covariance;
        double determinant;

        gsl_vector * mean;

        double weight;

        Implementation(const std::vector<double> & mean, const std::vector<double> & covariance, const double & weight) :
            dimension(mean.size()),
            covariance(gsl_matrix_alloc(dimension, dimension)),
            inverse_covariance(gsl_matrix_alloc(dimension, dimension)),
            mean(gsl_vector_alloc(dimension)),
            weight(weight)
        {
            if (covariance.size() != dimension * dimension)
                throw InternalError("HierarchicalClustering::Component: covariance and dimension do not match");

            std::copy(mean.cbegin(), mean.cend(), this->mean->data);
            std::copy(covariance.cbegin(), covariance.cend(), this->covariance->data);
            gsl_matrix * covariance_chol = gsl_matrix_alloc(dimension, dimension);

            // copy covariance matrix to covariance_chol
            gsl_matrix_memcpy(covariance_chol, this->covariance);

            // calculate cholesky decomposition, needed for sampling and one step for inversion
            gsl_error_handler_t * default_gsl_error_handler = gsl_set_error_handler_off();
            if (GSL_EDOM == gsl_linalg_cholesky_decomp(covariance_chol))
            {
                Log::instance()->message("HierarchicalClustering::Component", ll_warning)
                    << "Covariance matrix is not positive definite!"
                    << "Proceed by setting off-diagonal elements to zero.";

                // covariance_chol is potentially changed. Copy again
                gsl_matrix_memcpy(covariance_chol, this->covariance);

                // remove the off-diagonal elements of covariance_chol
                for (unsigned i = 0 ; i < dimension ; ++i)
                {
                    for (unsigned j = i + 1 ; j < dimension ; ++j)
                    {
                        gsl_matrix_set(covariance_chol, i, j, 0.0);
                        gsl_matrix_set(covariance_chol, j, i, 0.0);
                    }
                }

                if (GSL_EDOM == gsl_linalg_cholesky_decomp(covariance_chol))
                {
                    throw InternalError(
                         "HierarchicalClustering::Component: GSL couldn't find Cholesky decomposition of " + stringify(this->covariance->data, dimension, 4)
                        + "Apparently no moves were accepted, so try to increase number of iterations between updates "
                        + "or decrease initial proposal covariance. Proceed by taking square root of covariance manually");
                }
            }
            gsl_set_error_handler(default_gsl_error_handler);

            // copy cholesky decomposition to inverse_covariance
            gsl_matrix_memcpy(inverse_covariance, covariance_chol);

            // calculate the inverse of covariance
            gsl_linalg_cholesky_invert(inverse_covariance);

            // det(Sigma) = det(L)^2, and det(L) = Prod(diagonal)
            determinant = 1;
            for (unsigned i = 0 ; i < dimension ; ++i)
            {
                determinant *= gsl_matrix_get(covariance_chol, i, i);
            }
            determinant = power_of<2>(determinant);

            gsl_matrix_free(covariance_chol);
        }

        ~Implementation()
        {
            gsl_matrix_free(covariance);
            gsl_matrix_free(inverse_covariance);
            gsl_vector_free(mean);
        }
    };

    HierarchicalClustering::Component::Component(const std::vector<double> & mean, const std::vector<double> & covariance, const double & weight) :
        PrivateImplementationPattern<HierarchicalClustering::Component>(new Implementation<HierarchicalClustering::Component>(mean, covariance, weight))
    {
    }

    HierarchicalClustering::Component::~Component()
    {
    }

    HierarchicalClustering::Component::Component(const gsl_vector * mean, const gsl_matrix * covariance, const double & weight) :
        PrivateImplementationPattern<HierarchicalClustering::Component>(new Implementation<HierarchicalClustering::Component>(
            std::vector<double>(mean->data, mean->data + mean->size),
            std::vector<double>(covariance->data, covariance->data + covariance->size1 * covariance->size2),
            weight))
    {
    }

    gsl_matrix *
    HierarchicalClustering::Component::covariance() const
    {
        return _imp->covariance;
    }

    const gsl_matrix *
    HierarchicalClustering::Component::inverse_covariance() const
    {
        return _imp->inverse_covariance;
    }

    const double &
    HierarchicalClustering::Component::determinant() const
    {
        return _imp->determinant;
    }

    gsl_vector *
    HierarchicalClustering::Component::mean() const
    {
        return _imp->mean;
    }

    double &
    HierarchicalClustering::Component::weight() const
    {
        return _imp->weight;
    }

    std::ostream & operator<< (std::ostream & lhs, const HierarchicalClustering::Component & rhs)
    {
        lhs << "weight = " << rhs.weight();
        lhs << ", mean = " << stringify(rhs.mean()->data, rhs.mean()->data + rhs.mean()->size, 4);
        lhs << ", covariance = " << stringify(rhs.covariance()->data, rhs.mean()->size, 4);

        return lhs;
    }

    HierarchicalClustering::Config::Config() :
        equal_weights(true),
        kill_clusters(true),
        maximum_steps(std::numeric_limits<unsigned>::max()),
        precision(1e-4)
    {
    }

    HierarchicalClustering::Config
    HierarchicalClustering::Config::Default()
    {
        return Config();
    }
}
