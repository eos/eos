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

#include <test/test.hh>
#include <eos/utils/hierarchical-clustering.hh>

#include <algorithm>
#include <numeric>

#include <gsl/gsl_randist.h>

using namespace test;
using namespace eos;

class HierarchicalClusteringTest :
    public TestCase
{
    public:
        HierarchicalClusteringTest() :
            TestCase("hierarchical_clustering_test")
        {
        }

        virtual void run() const
        {
            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);

            // components along the unit circle in 2D
            // draw from clusters and try to find them again
            {
                HierarchicalClustering::MixtureDensity components;
                HierarchicalClustering::MixtureDensity clusters;

                static const unsigned n_clusters = 5;
                static const unsigned n_components = 20 * n_clusters;
                static const double radius = 5;

                std::vector<double> covariance
                {
                    1.0, 0.0,
                    0.0, 1.0
                };

                std::vector<double> mean(2, 0.0);

                gsl_rng_set(rng, n_clusters);

                // initialize clusters/components
                for (unsigned j = 0 ; j < n_clusters ; ++j)
                {
                    // draw an position
                    const double angle = double(j) / n_clusters * 2 * M_PI;
                    mean[0] = radius * std::cos(angle);
                    mean[1] = radius * std::sin(angle);

                    clusters.push_back(HierarchicalClustering::Component(mean, covariance, 1.0 / (n_clusters + 1)));

                    // draw components from clusters
                    for (unsigned i = 0 ; i < n_components / n_clusters ; ++i)
                    {
                        mean[0] = gsl_vector_get(clusters.back().mean(), 0) + gsl_ran_gaussian(rng, std::sqrt(covariance[0]) / 1.0);
                        mean[1] = gsl_vector_get(clusters.back().mean(), 1) + gsl_ran_gaussian(rng, std::sqrt(covariance[3]) / 1.0);

                        components.push_back(HierarchicalClustering::Component(mean, covariance, 1.0 / n_components));
                    }
                }
                // initialize clustering
                HierarchicalClustering::Config config = HierarchicalClustering::Config::Default();
                config.kill_clusters = true;
                HierarchicalClustering hc(config);
                for (unsigned i = 0 ; i < n_components ; ++i)
                {
                    hc.add(components[i]);
                }

                // transform clusters to make it harder
                for (unsigned j = 0 ; j < n_clusters ; ++j)
                {
                    gsl_vector_scale(clusters[j].mean(), 2.0);
                }

                // insert extra cluster in the middle to see it is killed
                mean[0] = radius * 2;
                mean[1] = radius * 2;
                clusters.insert(clusters.begin() + 2, HierarchicalClustering::Component(mean, covariance, 1.0 / 6));
                hc.initial_guess(clusters);
                hc.run();

                /* check the result */
                for (auto cl = hc.begin_clusters(); cl != hc.end_clusters() ; ++cl)
                {
                    TEST_CHECK_RELATIVE_ERROR(cl->weight(), 0.2, 1e-15);
                }

                for (auto map = hc.begin_map() ; map != hc.end_map() ; ++map)
                {
                    // integer division!
                    const unsigned cluster = std::distance(hc.begin_map(), map) / ( n_components / n_clusters);
                    // accounted for dead cluster automatically

                    // check that each component is associated with the cluster it was drawn from
                    TEST_CHECK_EQUAL(cluster, *map);
                }
            }
            gsl_rng_free(rng);
        }
} hierarchical_clustering_test;
