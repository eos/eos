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

#ifndef EOS_GUARD_SRC_STATISTICS_HIERARCHICAL_CLUSTERING_HH
#define EOS_GUARD_SRC_STATISTICS_HIERARCHICAL_CLUSTERING_HH 1

#include <eos/utils/wrapped_forward_iterator-impl.hh>
#include <eos/utils/private_implementation_pattern.hh>

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace eos
{
    /*!
     * Implement the hierarchical clustering explained in [GR2004].
     *
     * A mixture density is determined as a reduced representation of input components.
     */
    class HierarchicalClustering :
        public PrivateImplementationPattern<HierarchicalClustering>
    {
        public:
            class Config;
            class Component;
            typedef std::vector<HierarchicalClustering::Component> MixtureDensity;

            ///@name Basic Functions
            ///@{
            /*!
             * Constructor.
             *
             * @param config   The configuration of the samples.
             */
            HierarchicalClustering(const HierarchicalClustering::Config & config);

            /// Destructor.
            ~HierarchicalClustering();
            ///@}

            /*!
             * Add an input component of unit weight.
             */
            void add(const Component & component);

            /*!
             * Add an initial guess for the clusters to be determined.
             * @param density
             */
            void initial_guess(const MixtureDensity & density);

            /*!
             * Perform the clustering.
             */
            void run();

            /// Loop over input components
            struct InputIteratorTag;
            typedef WrappedForwardIterator<InputIteratorTag, Component> InputIterator;
            InputIterator begin_input() const;
            InputIterator end_input() const;

            /// Loop over output components (determined during clustering).
            struct OutputIteratorTag;
            typedef WrappedForwardIterator<OutputIteratorTag, Component> OutputIterator;
            OutputIterator begin_output() const;
            OutputIterator end_output() const;

            /// To which output is an input component mapped?
            struct MapIteratorTag;
            typedef WrappedForwardIterator<MapIteratorTag, unsigned> MapIterator;
            MapIterator begin_map() const;
            MapIterator end_map() const;
    };

    extern template class WrappedForwardIterator<HierarchicalClustering::InputIteratorTag, HierarchicalClustering::Component>;
    extern template class WrappedForwardIterator<HierarchicalClustering::OutputIteratorTag, HierarchicalClustering::Component>;
    extern template class WrappedForwardIterator<HierarchicalClustering::MapIteratorTag, unsigned>;

    /*!
     * Describe a component of a Gaussian mixture density, characterized
     * by mean, covariance and weight.
     */
    class HierarchicalClustering::Component :
        public PrivateImplementationPattern<HierarchicalClustering::Component>
    {
        public:
            Component(const std::vector<double> & mean, const std::vector<double> & covariance, const double & weight);
            Component(const gsl_vector * mean, const gsl_matrix * covariance, const double & weight);
            ~Component();

            gsl_matrix * covariance() const;
            const gsl_matrix * inverse_covariance() const;
            const double & determinant() const;
            gsl_vector * mean() const;
            double & weight() const;
    };

    std::ostream & operator<< (std::ostream & lhs, const HierarchicalClustering::Component & rhs);

    /*!
     * Stores all configuration options for a HierarchicalClustering.
     */
    class HierarchicalClustering::Config
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
             * HierarchicalClustering settings with reasonably chosen default values.
             */
            static Config Default();

            /*!
             * Named constructor
             *
             * HierarchicalClustering settings with values optimized for quick chain
             * convergence and evaluation.
             *
             * @note The convergence is not very reliable. Use with care!
             * If in doubt, use HierarchicalClustering::Config::Default
             */
            static Config Quick();
            ///@}

            /// Set component weights equal before the start of the clustering.
            bool equal_weights;

            /// If a component has zero weight, it is removed;
            bool kill_components;

            /// Perform a maximum number of update steps.
            unsigned maximum_steps;

            /// If relative change of distance between current and last step falls below precision,
            /// declare convergence.
            double precision;
    };
 }

#endif
