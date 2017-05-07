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

#ifndef EOS_GUARD_SRC_STATISTICS_PRIOR_SAMPLER_HH
#define EOS_GUARD_SRC_STATISTICS_PRIOR_SAMPLER_HH 1

#include <eos/statistics/log-prior-fwd.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/hdf5-fwd.hh>

#include <vector>

namespace eos
{
    /*!
     * Perform simple uncertainty propagation by defining parameters to be varied,
     * and the observables whose variation one is interested in.
     * Parameter values are sampled directly from 1D priors.
     * All observable values are stored to disk.
     */
    class PriorSampler :
        public PrivateImplementationPattern<PriorSampler>
    {
        public:
            struct Config;

            typedef hdf5::Array<1, double> ObservablesType;
            typedef std::vector<std::vector<double>> SamplesList;

            static ObservablesType observables_type(const unsigned & dimension);

            /*!
             * Constructor.
             * @param observables list of observables, can be extended through add().
             * @param config configuration options.
             */
            PriorSampler(const ObservableSet & observables, const Config & config);

            /*!
             * Destructor.
             */
            ~PriorSampler();

            /*!
             * Add a parameter with its prior to sample from.
             * @return false if prior was not added successfully.
             */
            bool add(const LogPriorPtr & prior);

            /*!
             * Add an observable for which the the uncertainty is to be evaluated.
             * @return false if identical to an existing observable.
             */
            bool add(const ObservablePtr & observable);

            /*!
             * Start the sampling process and store the calculated
             * observable values to disk.
             */
            void run();

            /*!
             * Calculate observables at the given parameter samples.
             * The order and meaning of parameters in a sample are
             * specified by the parameter descriptions.
             * @note No new samples are drawn from the priors.
             */
            void run(const SamplesList & samples, const std::vector<ParameterDescription> & );
    };

    /*!
     * Store configuration options
     */
    struct PriorSampler::Config
    {
        private:
            /// Constructor.
            Config();

        public:

            /// Constructor with the default settings.
            static Config Default();

            /// Total number of samples
            unsigned n_samples;

            /// Number of worker threads
            unsigned n_workers;

            /// The file where the observables are stored.
            std::shared_ptr<hdf5::File> output_file;

            /*!
             * If true, use as many threads as there are cores available.
             * If false, use only one thread.
             */
            bool parallelize;

            /// Seed for the random number generator.
            unsigned seed;

            /*!
             * If true, the parameter values are stored together with the observable values.
             * If false, only the observable values are stored.
             */
            bool store_parameters;
    };
}

#endif

