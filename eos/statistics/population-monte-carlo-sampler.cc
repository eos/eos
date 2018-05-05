/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 - 2013 Frederik Beaujean
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
#include <eos/statistics/population-monte-carlo-sampler.hh>

#include <eos/statistics/chain-group.hh>
#include <eos/statistics/hierarchical-clustering.hh>
#include <eos/statistics/log-posterior.hh>
#include <eos/statistics/markov-chain-sampler.hh>
#include <eos/statistics/proposal-functions.hh>
#include <eos/statistics/rvalue.hh>
#include <eos/statistics/welford.hh>
#include <eos/utils/density.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/thread_pool.hh>

#include <Minuit2/FCNGradientBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

extern "C" {
#include <pmclib/pmc.h>
}

#include <algorithm>
#include <math.h>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>

#include <gsl/gsl_randist.h>

using namespace eos::proposal_functions;

namespace eos
{
    // interface routines
    namespace pmc
    {
        /*
         * PMCError is parent to all exceptions thrown by the PMC library.
         */
        class PMC_Error :
            public Exception
        {
            public:
                /*!
                 * Constructor.
                 *
                 * @param message The error message.
                 */
                PMC_Error(const std::string & message) :
                    Exception(message)
                {
                }
        };

        /*
         * Wrapper around C-style error handler
         */
        class ErrorHandler
        {
            public:
                ErrorHandler()
                {
                    error_handler = initError();
                }

                ~ErrorHandler()
                {
                    endError(&error_handler);
                }

                operator error** ()
                {
                    return &error_handler;
                }

            private:
                error * error_handler;
        };

        // replace pmclib C style error handling with a true exception
        void check_error(error ** errorp) throw (PMC_Error)
        {
            if (!_isError(*errorp))
                return;

            error* err = (*errorp)->next;

            // traverse sequence of errors from method to method
            // until the issuing method is found
            while(err->next != NULL)
            {
                err = err->next;
            }
            std::string message =   "pmc_sampler::check_error: Found an error in pmc library.\n"
                                    "Error code is " + stringify(err->errValue) + "\n"
                                    "Error text is '" + err->errText + "'\n"
                                    "Error occurred in " + err->errWhere;
            // clean up before exception
            endError(errorp);
            throw PMC_Error(message);
        }

        /*
         * C-style interface to the unnormalized posterior. No bounds checking is done here,
         * as PMC discards those points during its sampling, and this function
         * is assumed to be called only on valid points.
         */
        double logpdf(void * data, const double * par_point, error ** /*error_handler*/)
        {
            // retrieve data
            Density * density = static_cast<Density *>(data);

            // set parameters
            unsigned j = 0;
            for (auto i = density->begin() ; i != density->end() ; ++i, ++j)
            {
                i->parameter->set(par_point[j]);
            }
            const double value = density->evaluate();
            if ( ! std::isfinite(value))
                throw InternalError("PMC::posterior: not finite " + stringify(value)
                                    + " at " + stringify(par_point, par_point + std::distance(density->begin(), density->end())));
            return value;
        }

        typedef std::pair<unsigned, double> index_pair;

        /*
         * Find the minimal partition of N into K parts, such that the smallest
         * and largest part differ by at most one.
         */
        void minimal_partition(const unsigned & N, const unsigned & K, std::vector<unsigned> & partition)
        {
            partition.resize(K);

            const unsigned remainder = N % K;
            const unsigned minimum = N / K;
            std::fill(partition.begin(), partition.begin() + remainder, minimum + 1);
            std::fill(partition.begin() + remainder, partition.end(), minimum);
        }

        hdf5::DataSet<PopulationMonteCarloSampler::Output::ComponentType>
        open_components(hdf5::File & f, const unsigned & n_dim, const bool & update)
        {
            std::string component_directory = "";
            if (update)
            {
                component_directory = "/data/initial/components";
                auto component_data_set = f.open_data_set(component_directory, PopulationMonteCarloSampler::Output::component_type(n_dim));
                return component_data_set;
            }

            try
            {
                component_directory = "/data/components";
                auto component_data_set = f.open_data_set(component_directory, PopulationMonteCarloSampler::Output::component_type(n_dim));
                return component_data_set;
            }
            catch (HDF5Error &)
            {
                component_directory = "/data/final/components";
                auto component_data_set = f.open_data_set(component_directory, PopulationMonteCarloSampler::Output::component_type(n_dim));
                return component_data_set;
            }
        }

        // Worker allows simple thread parallelization of massive posterior evaluation
       struct Worker
       {
           DensityPtr density;

           // store the posterior values
           std::vector<double> density_values;

           // points at which posterior is evaluated
           std::vector<double> parameter_samples;

           std::shared_ptr<ROOT::Minuit2::FunctionMinimum> minimum;

           Worker(const DensityPtr & density) :
               density(density->clone())
           {
           }

           void clear()
           {
               parameter_samples.resize(0);
               density_values.resize(0);
           }

           // call from main thread before actual work is done
           void setup(double * samples, unsigned n_samples, unsigned n_dim)
           {
               // copy the samples
               parameter_samples = std::vector<double>();
               std::copy(samples, samples + n_samples * n_dim, std::back_inserter(parameter_samples));

               density_values.resize(n_samples);
           }

           // compute log(posterior) at many sample points
           void work()
           {
               pmc::ErrorHandler err;

               if (parameter_samples.empty() || density_values.empty())
                   return;

               unsigned n_dim = parameter_samples.size() / density_values.size();
               unsigned i = 0;
               for (auto p = density_values.begin(), p_end = density_values.end(); p != p_end ; ++i, ++p)
               {
                    *p = pmc::logpdf(density.get(), &parameter_samples[i * n_dim], err);
               }
           }
       };
    }

    template<>
    struct Implementation<PopulationMonteCarloSampler>
    {
        typedef std::vector<unsigned> IndexList;

        // store reference, but don't own log-posterior
        DensityPtr density;

        // our configuration options
        PopulationMonteCarloSampler::Config config;

        // Keep track of the status
        PopulationMonteCarloSampler::Status status;

        // the pmc object
        pmc_simu * pmc;

        // random number generator
        gsl_rng * rng;

        // Workers do the hard part: calculating the posterior.
        std::vector<std::shared_ptr<pmc::Worker>> workers;

        // Posterior of the last sample
        std::vector<double> posterior_values;

        Implementation(const DensityPtr & density, const hdf5::File & file,
                       const PopulationMonteCarloSampler::Config & config, const bool & update) :
            density(density),
            config(config),
            status(),
            pmc(NULL)
        {
            // setup Mersenne-Twister RN generator using custom seed
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, config.seed);

            setup_output();

            // setup PMC library
            initialize_pmc(file, update);

            // initialization of the workers
            const unsigned number_of_workers = config.number_of_workers == 0 ?
                                               ThreadPool::instance()->number_of_threads() :
                                               config.number_of_workers;
            for (unsigned i = 0; i < number_of_workers ; ++i)
                workers.push_back(std::make_shared<pmc::Worker>(density));

            auto f = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);
            density->dump_descriptions(f, "/descriptions");
            dump_proposal("initial");
        }

        ~Implementation<PopulationMonteCarloSampler>()
        {
            // free RN generator
            gsl_rng_free(rng);

            // free pmc object
            pmc_simu_free(&pmc);
        }

        void calculate_weights(const std::string & sample_file,
                                 const unsigned & min_index,
                                 const unsigned & max_index)
        {
            pmc::ErrorHandler err;

            /* parse samples */

            const unsigned n_samples = max_index - min_index;
            pmc_simu_realloc(pmc, n_samples, err);
            pmc::check_error(err);

            auto f = hdf5::File::Open(sample_file);
            auto samples = f.open_data_set("/data/samples",
                PopulationMonteCarloSampler::Output::sample_type(pmc->ndim));
            auto sample_record = PopulationMonteCarloSampler::Output::sample_record(pmc->ndim);
            samples.set_index(min_index);

            for (unsigned i = 0 ; i < n_samples; i++)
            {
                samples >> sample_record;

                std::copy(sample_record.cbegin(), sample_record.cbegin() + pmc->ndim, &pmc->X[i * pmc->ndim]);

                // index of generating component
                pmc->indices[i] = sample_record.at(pmc->ndim);

                // ignore posterior value and weight of record
            }

            /* do the hard computational work */

            calculate_weights();

            // do not normalize, as the total sum of all is not known from this subsample
//            normalize_importance_weight(pmc, err);
//            pmc::check_error(err);

            /* dump samples */

            // open the file whenever writing is desired, so it is in a readable state during lengthy posterior calculations
            auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            auto weights = file.create_data_set("/data/weights",
                PopulationMonteCarloSampler::Output::weight_type());
            auto weights_record = PopulationMonteCarloSampler::Output::weight_record();

            auto ignores = file.create_data_set("/data/broken", PopulationMonteCarloSampler::Output::ignore_type());
            for (unsigned i = 0 ; i < n_samples ; i++)
            {
                // posterior
                std::get<0>(weights_record) = posterior_values[i];

                // weight
                std::get<1>(weights_record) = pmc->weights[i];

                weights << weights_record;

                // don't ignore this value
                ignores << short(0);
            }
        }

        /*
         * Assuming that samples from proposal densities have been created and stored in pmc->X,
         * calculate the posterior values at those samples, optionally in parallel or sequentially.
         */
        void calculate_weights()
        {
            pmc::ErrorHandler err;

            // how many samples per average worker
            unsigned average_samples_per_worker = pmc->nsamples / ThreadPool::instance()->number_of_threads();
            unsigned remainder = pmc->nsamples % ThreadPool::instance()->number_of_threads();

            const unsigned n_dim = std::distance(density->begin(), density->end());

            // tickets for parallel computations
            std::vector<Ticket> tickets;

            Log::instance()->message("PMC_sampler.status", ll_debug)
                << "Workers started";

            // setup workers
            for (unsigned i = 0 ; i < ThreadPool::instance()->number_of_threads() ; ++i)
            {
                unsigned samples_per_worker = average_samples_per_worker;

                // last worker gets the remainder
                if (i == ThreadPool::instance()->number_of_threads() - 1)
                    samples_per_worker += remainder;

                workers[i]->setup(&pmc->X[i * average_samples_per_worker * n_dim], samples_per_worker, n_dim);

                if (config.parallelize)
                    // make sure to pass the pointer, instead of a reference from *w, to bind.
                    // Else temporary copies are created.
                    tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&pmc::Worker::work, workers[i].get())));
                else
                    workers[i]->work();
            }

            // why not clear after weights are assigned?
            posterior_values.clear();

            // wait for job completion
            for (auto t = tickets.begin(), t_end = tickets.end() ; t != t_end ; ++t)
                t->wait();

            // copy results and free memory
            for (auto w = workers.begin(), w_end = workers.end() ; w != w_end ; ++w)
            {
                std::move((**w).density_values.begin(), (**w).density_values.end(), std::back_inserter(posterior_values));
                (**w).clear();
            }

            Log::instance()->message("PMC_sampler.status", ll_debug)
                << "Workers finished";

            double max_rho = 0;
            double max_weight = 0;

            // read sample by sample
            for (unsigned i = 0 ; i < unsigned(pmc->nsamples) ; ++i)
            {
                double * x = &(pmc->X[i * n_dim]);
                pmc->flg[i] = 0;

                /* Compute log(density) according to proposal */

                // rho?, the density of the proposal
                const double rloc = distribution_lkl(pmc->proposal, x, err);
                pmc::check_error(err);

                if ( (i == 0) || rloc > max_rho)
                  max_rho = rloc;
                pmc->log_rho[i] = rloc;

                /* Compute log(weight) = log(posterior) - log(proposal) */
                const double weight = posterior_values[i] - rloc;

                // no support for reduced parameters!

                if ( (i == 0) || weight > max_weight)
                  max_weight = weight;
                pmc->weights[i] = weight;

                if (!finite(rloc))
                    throw InternalError("PMC::calculate_weights: proposal density not finite " + stringify(rloc)
                            + " at " + stringify(x, x + n_dim));
                if (!finite(posterior_values[i]))
                    throw InternalError("PMC::calculate_weights: posterior density not finite " + stringify(posterior_values[i])
                            + " at " + stringify(x, x + n_dim));

                // Everything fine, set flag to one and continue
                pmc->flg[i]  = 1;

                //remember that we only computed log
                pmc->maxW = max_weight;
                pmc->maxR = max_rho;

                pmc->isLog = 1; //todo save in HDF5 as well
            }

            endError(err);
        }

        /*
         * Several criteria are possible
         * 1. 1-perplexity < epsilon && 1 - ess/N < epsilon: very sure about convergence
         * 2. perplexity > 0.5 and perplexity changed only little over last two iterations:
         *    issue warning, perhaps with more components/different DoF perplexity could be raised
         *
         * @param data_set Look only at the values of a particular data set, else scan all previous ones
         *                 in '/data/STEP/statistics'
         */
        bool check_convergence(const std::string & file_name, const std::string & data_set_name = "")
        {
            Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                << "perplexity = " << stringify(status.perplexity, 4)
                << ", effective sample size = " + stringify(status.eff_sample_size, 4);
            if (status.perplexity > config.convergence_perplexity &&
                (config.ignore_eff_sample_size or status.eff_sample_size > config.convergence_eff_sample_size))
            {
                Log::instance()->message("PMC_sampler.check_convergence", ll_debug)
                    << "perplexity (" << stringify(status.perplexity, 4) << ")"
                    << (config.ignore_eff_sample_size ? "" :
                       " and effective sample size " + stringify(status.eff_sample_size, 4))
                   << " large enough";
                return true;
            }

            // read out past perplexity from HDF5
            std::vector<PopulationMonteCarloSampler::Status> past_status;
            auto file = hdf5::File::Open(file_name, H5F_ACC_RDONLY);
            auto statistics_record = std::make_tuple(status.perplexity, status.eff_sample_size, status.evidence);

            H5E_BEGIN_TRY
            {
                try
                {
                    unsigned step = 0;
                    while(true)
                    {
                        std::string sub_directory = data_set_name.empty() ? "/data/" + stringify(step) + "/statistics" : data_set_name;
                        // need to find out the right type
                        auto statistics_data_set = file.open_data_set(sub_directory, PopulationMonteCarloSampler::Output::statistics_type());
                        statistics_data_set.set_index(data_set_name.empty() ? 0 : step);
                        statistics_data_set >> statistics_record;

                        // assign values only to interesting member of Status
                        past_status.push_back(PopulationMonteCarloSampler::Status());
                        past_status.back().perplexity = std::get<0>(statistics_record);
                        past_status.back().eff_sample_size = std::get<1>(statistics_record);

                        // update
                        ++step;

                        // read in all records with massive parallelization
                        if (! data_set_name.empty() && step == statistics_data_set.records())
                            break;
                    }
                }
                catch (HDF5Error &)
                {
                }
            }
            H5E_END_TRY;

            // add current results
            if (! data_set_name.empty())
                past_status.push_back(status);

            // variance for at least two results
            if (past_status.size() < config.minimum_steps)
            {
                Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                    << "Found " << past_status.size()
                    << " status, but need at least " << config.minimum_steps
                    << " to define convergence based on previous steps";

                return false;
            }

            // compute mean and variance of perplexity of last n steps
            Welford welford_perplexity;
            Welford welford_eff_sample_size;

            for (auto s = past_status.rbegin() ; s != past_status.rbegin() + config.minimum_steps ; ++s)
            {
                welford_eff_sample_size.add(s->eff_sample_size);
                welford_perplexity.add(s->perplexity);
            }

            const double relative_std_deviation_perplexity = welford_perplexity.std_deviation() / welford_perplexity.mean();
            const double relative_std_deviation_ess = welford_eff_sample_size.std_deviation() / welford_eff_sample_size.mean();

            // require a minimum value
            if (welford_perplexity.mean() < config.minimum_perplexity)
            {
                Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                    << "perplexity mean too small: " << stringify(welford_perplexity.mean(), 4)
                    << " < " << config.minimum_perplexity;
                return false;
            }

            // require that std_deviation be small enough
            if (relative_std_deviation_perplexity > config.maximum_relative_std_deviation)
            {
                Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                    << "perplexity relative std_deviation too large: " << stringify(relative_std_deviation_perplexity, 4)
                    << " > " << config.maximum_relative_std_deviation;
                return false;
            }

            // same for eff. sample size
            if (not config.ignore_eff_sample_size)
            {
                if (welford_eff_sample_size.mean() < config.minimum_eff_sample_size)
                {
                    Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                        << "ESS mean too small: " << stringify(welford_eff_sample_size.mean(), 4)
                        << " < " << config.minimum_eff_sample_size;
                    return false;
                }
                if (relative_std_deviation_ess > config.maximum_relative_std_deviation)
                {
                    Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                        << "ESS relative std_deviation too large: " << stringify(relative_std_deviation_ess, 4)
                        << " > " << config.maximum_relative_std_deviation;
                    return false;
                }
            }

            Log::instance()->message("PMC_sampler.check_convergence", ll_informational)
                << "Mean and relative std. deviation of perplexity (" << stringify(welford_perplexity.mean(), 4)
                << ", " << stringify(relative_std_deviation_perplexity, 4) << ")"
                << (config.ignore_eff_sample_size ? "" :
                   " and of ESS (" + stringify(welford_eff_sample_size.mean(), 4)
                   + ", " + stringify(relative_std_deviation_ess, 4) + ")")
                << " are OK";

            return true;
        }

        void crop_weights()
        {
            if ( ! config.crop_highest_weights)
                return;

            Log::instance()->message("PMC_sampler.update", ll_informational)
                << "Cropping " << config.crop_highest_weights
                << " highest weights";

            // associate points with index
            std::vector<pmc::index_pair> weight_indices(pmc->nsamples);
            for (unsigned i = 0 ; i < weight_indices.size() ; ++i)
            {
                weight_indices[i] = pmc::index_pair(i, pmc->weights[i]);
            }

            // sort according to posterior in descending order
            std::sort(weight_indices.begin(), weight_indices.end(),
                      [] (const pmc::index_pair & a, const pmc::index_pair & b)
            {
                return a.second > b.second;
            } );

            // set OK flag to false for highest weights
            for (unsigned j = 0 ; j < config.crop_highest_weights ; ++j)
            {
                pmc->flg[weight_indices[j].first] = 0;
            }
        }

        void dump_proposal(const std::string & group)
        {
            // open the file whenever writing is desired, so it is in a readable state during lengthy posterior calculations
            auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            /* dump components */

            // proposal density
            mix_mvdens * prop = static_cast<mix_mvdens *>(pmc->proposal->data);

            auto components = file.create_data_set("/data/" + group + "/components",
                PopulationMonteCarloSampler::Output::component_type(pmc->ndim));
            auto component_record = PopulationMonteCarloSampler::Output::component_record(pmc->ndim);
            auto dof = components.create_attribute("dof", hdf5::Scalar<int>("dof"));
            dof = config.degrees_of_freedom;

            // save whether std contains the actual covariance matrix or the GSL cholesky decomposition
            auto chol = components.create_attribute("chol", hdf5::Scalar<int>("chol"));
            chol = prop->comp[0]->chol;

            unsigned dead_components = 0;

            for (unsigned i = 0 ; i < prop->ncomp ; ++i)
            {
                std::get<0>(component_record) = prop->wght[i];
                std::copy(prop->comp[i]->mean, prop->comp[i]->mean + pmc->ndim, std::get<1>(component_record).begin());
                std::copy(prop->comp[i]->std, prop->comp[i]->std + pmc->ndim * pmc->ndim, std::get<2>(component_record).begin());

                components << component_record;

                if (prop->wght[i] == 0)
                {
                    ++dead_components;
                    continue;
                }
            }

            Log::instance()->message("PMC_sampler.dump", ll_informational)
                << dead_components << " out of " << prop->ncomp << " components died out.";

            if (group == "initial")
                return;
        }

        /*!
         * Dump status to HDF5.
         *
         * @param group
         * @param samples If false, only summary statistics are stored
         */
        void dump(const std::string & group, bool store_samples = true)
        {
            const unsigned dim = std::distance(density->begin(), density->end());

            // open the file whenever writing is desired, so it is in a readable state during lengthy posterior calculations
            auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            /* dump statistics information */

            auto statistics = file.create_data_set("/data/" + group + "/statistics",
                PopulationMonteCarloSampler::Output::statistics_type());
            auto statistics_record = std::make_tuple(status.perplexity, status.eff_sample_size, status.evidence);
            statistics << statistics_record;

            /* dump samples */
            if (not store_samples)
                return;

            auto samples = file.create_data_set("/data/" + group + "/samples",
                PopulationMonteCarloSampler::Output::sample_type(dim));
            auto sample_record = PopulationMonteCarloSampler::Output::sample_record(dim);
            for (int i = 0 ; i < pmc->nsamples ; i++)
            {
                sample_record.clear();

                for (int j = 0 ; j < pmc->ndim ; ++j)
                {
                    sample_record.push_back(pmc->X[i * pmc->ndim + j]);
                }

                // index of generating component
                sample_record.push_back(pmc->indices[i]);

                // posterior
                sample_record.push_back(posterior_values[i]);

                // log weight
                double logw = pmc->weights[i];
                if (!pmc->isLog)
                    logw = std::log(logw);

                /* Back to unnormalized weights (log) */
                logw += pmc->logSum;
                sample_record.push_back(logw);

                samples << sample_record;
            }
        }

        void initialize_pmc(const hdf5::File & file, const bool & update)
        {
            Log::instance()->message("PMC_sampler::initialize", ll_informational)
                << "Reading from file " << file.name();

            // we will only use the file for reading.
            hdf5::File & f = const_cast<hdf5::File &>(file);

            // number of dimensions of parameter cube
            const size_t n_dim = std::distance(density->begin(), density->end());

            // exception handler look-alike
            pmc::ErrorHandler err;

            // parameter cube: copy from log-posterior
            parabox * par_box = init_parabox(n_dim, err);
            int i = 0;
            for (auto & d : *density)
            {
                add_slab(par_box, i, d.min, d.max, err);
                ++i;
            }

            /* setup importance sampling */

            // create posterior distribution
            distribution * target = init_simple_distribution(n_dim, (void *) &density, &pmc::logpdf, NULL, err);
            pmc::check_error(err);

            /* create proposal density from file */

            // use functions in pmc_rc.c:
            // o rcinit_mix_mvdens(...)
            // o init_importance_from_rc(...)
            mix_mvdens * mmv = NULL;

            /* determine type of input file */

            unsigned number_of_live_components = 0;

            if (config.target_ncomponents > 0)
            {
                // initialize mmv by clustering and filtering
                number_of_live_components = hierarchical_clustering(f, mmv);
            }
            // read in from previous PMC output
            else
            {
                auto component_data_set = pmc::open_components(f, n_dim, update);

                mmv = mix_mvdens_alloc(component_data_set.records(), n_dim, err);

                auto record = PopulationMonteCarloSampler::Output::component_record(n_dim);

                for (unsigned i = 0 ; i < mmv->ncomp ; ++i)
                {
                    component_data_set >> record;
                    mmv->wght[i] = std::get<0>(record);
                    if (mmv->wght[i] > 0)
                    {
                        ++number_of_live_components;
                    }
                    std::copy(std::get<1>(record).cbegin(), std::get<1>(record).cend(), mmv->comp[i]->mean);
                    std::copy(std::get<2>(record).cbegin(), std::get<2>(record).cend(), mmv->comp[i]->std);

                    mmv->comp[i]->band_limit = n_dim;
                    mmv->comp[i]->df = component_data_set.open_attribute("dof", hdf5::Scalar<int>("dof")).value();
                    mmv->comp[i]->chol = component_data_set.open_attribute("chol", hdf5::Scalar<int>("chol")).value();

                    // compute determinant of cholesky matrix if needed
                    if (mmv->comp[i]->chol)
                    {
                        mmv->comp[i]->detL = determinant(mmv->comp[i]->std, n_dim);
                    }
                }
            }

            /* final part */

            if (mmv->ndim != n_dim)
                throw InternalError("PMC::ctor: mismatch of parameter dimensions of log-posterior vs proposal ("\
                                    + stringify(n_dim)
                                    + " vs " + stringify(mmv->ndim) + ")");

            // create proposal distribution object
            distribution * proposal = mix_mvdens_distribution(mmv->ndim, (void*) mmv, err);

            // number of samples per chunk. Fixed size for each live component.
            unsigned n_samples = config.samples_per_component * number_of_live_components;

            if (update)
            {
                auto samples_data_set = f.open_data_set("/data/samples", PopulationMonteCarloSampler::Output::sample_type(n_dim));
                n_samples = samples_data_set.records();
            }

            // setup with 0 (?) deduced parameters for now
            pmc = pmc_simu_init_plus_ded(n_samples, target->ndim, target->n_ded, err);

            // setter methods, include safety checks of variable
            pmc_simu_init_target(pmc, target, par_box, err);
            pmc_simu_init_proposal(pmc, proposal, config.print_steps, err);
            pmc_simu_init_pmc(pmc, NULL, NULL, update_prop_rb_void, err);

            if (update)
            {
                this->update(f, n_samples);
            }
        }

        std::vector<ChainGroup> group_chains(const std::vector<HistoryPtr> & chains)
        {
            std::list<HistoryPtr> available_chains(chains.begin(), chains.end());

            ChainGroup::RValueFunction r = &RValue::approximation;

            // create indices of parameters whose R-value is checked
            std::vector<unsigned> parameter_indices;
            {
                unsigned i = 0;
                for (auto & d : *density)
                {
                    if ( d.nuisance && config.r_value_no_nuisance)
                        continue;

                    parameter_indices.push_back(i);
                    ++i;
                }
            }

            // add first chain as first group
            std::vector<ChainGroup> groups{ ChainGroup(r, config.group_by_r_value,
                    available_chains.front(), 0, config.skip_initial) };
            groups.back().parameter_indices(parameter_indices);

            available_chains.pop_front();

            /* add chain to a group if R value with existing chains in group small enough */
            unsigned chain_index = 0;
            while (! available_chains.empty() && ++chain_index)
            {
                // try to add a single chain to an existing group
                bool added = false;

                for (auto c = groups.begin(), c_end = groups.end() ; c != c_end ; ++c)
                {
                    if (! c->overlaps(available_chains.front()))
                        continue;

                    c->add(available_chains.front(), chain_index);
                    added = true;

                    Log::instance()->message("PMC.hierarchical_clustering", ll_debug)
                                    << "Added chain " << chain_index << " to group " << std::distance(groups.begin(), c);

                    break;
                }

                if (! added)
                {
                    groups.push_back(ChainGroup(r, config.group_by_r_value,
                            available_chains.front(), chain_index, config.skip_initial));
                    groups.back().parameter_indices(parameter_indices);

                    Log::instance()->message("PMC.hierarchical_clustering", ll_debug)
                                    << "Created new group for chain " << chain_index;
                }
                available_chains.pop_front();
            }

            /* copy the groups */

            // store lengths of groups for output only
            std::vector<unsigned> group_sizes;
            for (auto g : groups)
            {
                group_sizes.push_back(std::distance(g.begin(), g.end()));
            }

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                            << "Found " << group_sizes.size() << " groups of chains with " << stringify_container(group_sizes)
                            << " members";

            if ( ! config.ignore_groups.empty())
            {
                // sort to iterate in reverse (remove group with largest index first to preserve meaning of index)
                std::vector<unsigned> ignore_groups = config.ignore_groups;
                std::sort(ignore_groups.begin(), ignore_groups.end());

                // remove all duplicate entries
                std::unique(ignore_groups.begin(), ignore_groups.end());

                // remove group
                for (auto i = ignore_groups.crbegin() ; i != ignore_groups.crend() ; ++i)
                {
                    if (*i >= groups.size())
                    {
                        Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_error)
                            << "Skipping invalid ignore group: " << *i;
                        continue;
                    }
                    group_sizes.erase(group_sizes.begin() + *i);

                    Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_warning)
                                    << "Removing group " << *i;
                }

                Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                                << "Using " << group_sizes.size() << " groups of chains with " << stringify_container(group_sizes)
                                << " members";
            }
            return groups;
        }

        /*
         * Finding the initial component guess:
         * 1. Respect chains: choose one component only from patches within a chain
         * 2. Fix desired #components and #samples for a patch
         */
        unsigned hierarchical_clustering(const hdf5::File & file, mix_mvdens * & mmv)
        {
            /* parse chain histories */

            std::vector<HistoryPtr> chains;
            {
                std::vector<std::shared_ptr<hdf5::File>> input_files;
                input_files.push_back(std::make_shared<hdf5::File>(hdf5::File::Open(file.name(), H5F_ACC_RDONLY)));
                chains = MarkovChainSampler::read_chains(input_files);
            }

            // number of parameters
            const unsigned & ndim = chains.front()->states.front().point.size();

            if (ndim != std::distance(density->begin(), density->end()))
            {
                Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_warning)
                    << "The log-posterior in MCMC prerun had dim " << ndim
                    << ", but now the log-posterior has dim " << std::distance(density->begin(), density->end());
            }

            HierarchicalClustering::Config conf = HierarchicalClustering::Config::Default();
            conf.equal_weights = true;
            HierarchicalClustering hc(conf);

            /* group chains according to R-value */

            std::vector<ChainGroup> chain_groups = group_chains(chains);

            /* create initial guess for components by drawing local patches uniformly w/o replacement or large windows */

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Creating initial guess for the " << config.target_ncomponents
                << " target components to be formed from large windows"
                << (config.group_by_r_value > 1 ? " for each of the " + stringify(chain_groups.size()) + " chain groups found" : "");

            HierarchicalClustering::MixtureDensity initial_components;

            // weight of a single component (weights sum up to one)
            const unsigned n_components_total = config.target_ncomponents * chain_groups.size();
            const double weight = 1.0 / n_components_total;

            {
                for (auto g = chain_groups.cbegin() ; g != chain_groups.cend() ; ++g)
                {
                    // how many components should each each in group contribute
                    std::vector<unsigned> components_per_chain;
                    pmc::minimal_partition(config.target_ncomponents, std::distance(g->begin(), g->end()), components_per_chain);

                    auto n_components = components_per_chain.cbegin();
                    for (auto c = g->begin(); c != g->end(); ++c, ++n_components)
                    {
                        // throw away component
                        if (! *n_components)
                            continue;

                        auto first_state = (**c).states.cbegin() + config.skip_initial * (**c).states.size();
                        const int window = std::distance(first_state, (**c).states.cend()) / (*n_components);
                        if (window < 0)
                        {
                            throw InternalError("PMC::hierarchical_clustering: number of components too large for history size and skip initial: " + stringify(window)
                                    + " vs " + stringify(std::distance(first_state, (**c).states.cend())) + " and " + stringify(config.skip_initial));
                        }
                        auto last_state = first_state + window;

                        bool done = false;
                        while (! done)
                        {
                            // add remainder to last component
                            if (std::distance(last_state, (**c).states.cend()) < window)
                            {
                                last_state = (**c).states.cend();
                                done = true;
                            }

                            std::vector<double> mean, covariance;
                            (**c).mean_and_covariance(first_state, last_state, mean, covariance);
                            std::vector<double> center = mean;

                            initial_components.push_back(HierarchicalClustering::Component(center, covariance, weight));

                            // update range for next iteration
                            first_state += window;
                            last_state += window;
                        }
                    }
                }
            }

            hc.initial_guess(initial_components);

            /* create patches from each chain */

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Creating patches of length " << config.patch_length;

            HierarchicalClustering::MixtureDensity local_patches;

            for (auto g = chain_groups.cbegin() ; g != chain_groups.cend() ; ++g)
            {
                for (auto c = g->begin() ; c != g->end() ; ++c)
                {
                    //                local_patches_index_lists.push_back(IndexList());

                    auto first_state = (**c).states.cbegin() + config.skip_initial * (**c).states.size();
                    auto last_state = first_state + config.patch_length;

                    if (std::distance(last_state, (**c).states.cend()) < 0)
                    {
                        throw InternalError("PMC::hierarchical_clustering: sliding window too large for history size and skip initial: " + stringify(config.patch_length)
                                + " vs " + stringify(std::distance(first_state, (**c).states.cend())) + " and " + stringify(config.skip_initial));
                    }

                    bool done = false;
                    while (! done)
                    {
                        // add remainder to last patch
                        if (unsigned(std::distance(last_state, (**c).states.cend())) < config.patch_length)
                        {
                            last_state = (**c).states.cend();
                            done = true;
                        }

                        std::vector<double> mean, covariance;
                        (**c).mean_and_covariance(first_state, last_state, mean, covariance);
                        std::vector<double> center = mean;

                        try
                        {
                            //later on during clustering use equal weights for each component
                            HierarchicalClustering::Component patch(center, covariance, 1.0);
                            local_patches.push_back(patch);
                            hc.add(patch);
                        }
                        catch (InternalError &)
                        {
                            Log::instance()->message("PMC_sampler.hierarchical_clustering.add_patch", ll_debug)
                                    << "Skipping component, probably sliding window too small, and covariance not defined";
                        }

                        // update range for next iteration
                        first_state += config.patch_length;
                        last_state += config.patch_length;
                    }
                }
            }
            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Formed " << local_patches.size() << " input components centered around patch means";

            if (config.store_input_components)
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

                auto components = file.create_data_set("/hc/input-components",
                    PopulationMonteCarloSampler::Output::component_type(ndim));
                auto component_record = PopulationMonteCarloSampler::Output::component_record(ndim);

                for (auto comp = hc.begin_input() ; comp != hc.end_input() ; ++comp)
                {
                    std::get<0>(component_record) = comp->weight();
                    std::copy(comp->mean()->data, comp->mean()->data + ndim,
                              std::get<1>(component_record).begin());
                    std::copy(comp->covariance()->data, comp->covariance()->data + ndim * ndim,
                              std::get<2>(component_record).begin());

                    components << component_record;
                }
            }

            if (config.store_hc_initial)
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

                auto components = file.create_data_set("/hc/initial-guess",
                    PopulationMonteCarloSampler::Output::component_type(ndim));
                auto component_record = PopulationMonteCarloSampler::Output::component_record(ndim);

                for (auto comp = hc.begin_output() ; comp != hc.end_output() ; ++comp)
                {
                    std::get<0>(component_record) = comp->weight();
                    std::copy(comp->mean()->data, comp->mean()->data + ndim,
                              std::get<1>(component_record).begin());
                    std::copy(comp->covariance()->data, comp->covariance()->data + ndim * ndim,
                              std::get<2>(component_record).begin());

                    components << component_record;
                }
            }

            /* create components from patches  */

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Start hierarchical clustering ";

            hc.run();

            /* initialize pmc */

            // count active components
            unsigned active_components = 0;
            for (auto cl = hc.begin_output() ; cl != hc.end_output() ; ++cl, ++active_components);

            // did components die?
            if (active_components != n_components_total)
            {
                Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_warning)
                    << "Using only " << active_components << " components to start PMC. "
                    << n_components_total - active_components << " died out during the hierarchical clustering.";
            }

            pmc::ErrorHandler err;

            mmv = mix_mvdens_alloc(active_components, ndim, err);

            unsigned i_cl = 0;
            for (auto cl = hc.begin_output() ; cl != hc.end_output() ; ++cl, ++i_cl)
            {
                // equal weight to every component
                mmv->wght[i_cl] = 1.0 / active_components;
                mvdens * mv = mmv->comp[i_cl];
                std::copy(cl->mean()->data, cl->mean()->data + ndim, mv->mean);
                std::copy(cl->covariance()->data, cl->covariance()->data + ndim * ndim, mv->std);
                mv->band_limit = ndim;
                mv->df = config.degrees_of_freedom;
                mv->chol = 0;
            }

            return active_components;
        }

        void pre_run()
        {
            Log::instance()->message("PMC_sampler.status", ll_informational)
                << "Starting the prerun";

            pmc::ErrorHandler err;

            // prerun to adapt proposal densities
            for (unsigned i = 0; i < config.max_updates ; ++i)
            {
                dump_proposal(stringify(i));
                //                pmc_simu_pmc_step(pmc_simu *pmc, gsl_rng *r, error **err)
                {
                    // hack my own replacement for
                    // pmc_simu_importance(pmc,r,err);
                    {
                        // create parameter samples from the components
                        Log::instance()->message("PMC_sampler.status", ll_debug)
                            << "Drawing samples";
                        pmc->proposal->simulate(pmc, pmc->proposal->data, rng, pmc->pb, err);
                        pmc::check_error(err);

                        Log::instance()->message("PMC_sampler.status", ll_debug)
                            << "Calculating " << pmc->nsamples << " samples";
                        calculate_weights();
                    }
                    // remove highest weights if desired, needs to come before weight normalization
                    crop_weights();

                    normalize_importance_weight(pmc, err);
                    pmc::check_error(err);

                    Log::instance()->message("PMC_sampler.status", ll_informational)
                        << "Updating the proposal function";
                    pmc->pmc_update(pmc->proposal->data, pmc, err);
                    pmc::check_error(err);
                }

                // both perplexity and ess in [0, 1]
                status.perplexity = perplexity_and_ess(pmc, MC_NORM, &status.eff_sample_size, err);
                status.eff_sample_size /= pmc->nsamples;
                status.evidence = evidence(pmc, NULL, err);
                Log::instance()->message("PMC_sampler.status", ll_informational)
                    << "Status after step " << i + 1 << " of " << config.max_updates
                    << " with " << pmc->nsamples << " samples:";
                Log::instance()->message("PMC_sampler.status", ll_informational)
                    << "perplexity: " << status.perplexity
                    << ", eff. sample size: " << status.eff_sample_size
                    << ", evidence: " << status.evidence;

                // store statistics, but last samples only when requested
                dump(stringify(i), config.store_prerun);

                //check number of live components and adjust sample size accordingly
                unsigned live_components = 0;
                mix_mvdens * mmv = static_cast<mix_mvdens *>(pmc->proposal->data);
                for (unsigned k = 0 ; k < mmv->ncomp ; ++k)
                {
                    if (mmv->wght[k] > 0)
                        ++live_components;
                }

                if (config.adjust_sample_size)
                {
                    pmc_simu_realloc(pmc, config.samples_per_component * live_components, err);
                    pmc::check_error(err);
                }

                if ((status.converged = check_convergence(config.output_file)))
                {
                    Log::instance()->message("PMC_sampler.status", ll_informational)
                        << "Convergence achieved after " << i + 1 << " steps.";
                    status.iterations_at_convergence = i;
                    break;
                }
            }

            if (! status.converged)
            {
                Log::instance()->message("PMC_sampler.status", ll_warning)
                    << "Pre-run did NOT converge!";
            }
        }

        static void read_samples(const std::string & sample_file, const std::string & base,
                                 const unsigned & min, const unsigned & max,
                                 std::vector<std::vector<double>> & samples)
        {
            const unsigned n_dim = samples.front().size();
            samples.clear();

            auto file = hdf5::File::Open(sample_file, H5F_ACC_RDONLY);

            auto data_set = file.open_data_set(base + "/samples", PopulationMonteCarloSampler::Output::sample_type(n_dim));
            auto record = PopulationMonteCarloSampler::Output::sample_record(n_dim);
            data_set.set_index(min);
            for (unsigned i = min ; i < max; ++i)
            {
                data_set >> record;
                samples.push_back(std::vector<double>(record.begin(), record.end() - 3));
            }
        }

        void run()
        {
            pmc::ErrorHandler err;

            if (config.need_prerun)
                pre_run();

            if (config.final_samples == 0)
                return;

            // reserve memory for final chunk
            pmc_simu_realloc(pmc, config.final_samples, err);

            // sample from proposal density
            pmc->proposal->simulate(pmc, pmc->proposal->data, rng, pmc->pb, err);
            pmc::check_error(err);

            // calculate posterior of the samples
            calculate_weights();

            normalize_importance_weight(pmc, err);

            // both perplexity and ess in [0, 1]
            status.perplexity = perplexity_and_ess(pmc, MC_NORM, &status.eff_sample_size, err);
            status.eff_sample_size /= double(pmc->nsamples);
            status.evidence = evidence(pmc, NULL, err);
            Log::instance()->message("PMC_sampler.status", ll_informational)
                << "Status after final step with " << pmc->nsamples << " samples:";
            Log::instance()->message("PMC_sampler.status", ll_informational)
                << "perplexity: " << status.perplexity
                << ", eff. sample size: " << status.eff_sample_size
                << ", evidence: " << status.evidence;

            // store results
            dump_proposal("final");
            if (config.store)
                dump("final");
        }

        void setup_output() const
        {
            if (config.output_file.empty())
            {
                Log::instance()->message("PMC_sampler.setup_output", ll_warning)
                    << "No output file specified, results of sampling will not be stored!";
            }

            //  overwrite existing file
            hdf5::File::Create(config.output_file);
        }

        void draw_samples()
        {
            pmc::ErrorHandler err;

            auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            /* dump components */

            // proposal density
            mix_mvdens * prop = static_cast<mix_mvdens *>(pmc->proposal->data);

            auto components = file.create_data_set("/data/components",
                PopulationMonteCarloSampler::Output::component_type(pmc->ndim));
            auto component_record = PopulationMonteCarloSampler::Output::component_record(pmc->ndim);
            auto dof = components.create_or_open_attribute("dof", hdf5::Scalar<int>("dof"));
            dof = config.degrees_of_freedom;

            // save whether std contains the actual covariance matrix or the GSL cholesky decomposition
            auto chol = components.create_or_open_attribute("chol", hdf5::Scalar<int>("chol"));
            chol = prop->comp[0]->chol;

            unsigned live_components = 0;
            for (unsigned i = 0 ; i < prop->ncomp ; ++i)
            {
                std::get<0>(component_record) = prop->wght[i];
                if (prop->wght[i] > 0.0)
                    ++live_components;

                std::copy(prop->comp[i]->mean, prop->comp[i]->mean + pmc->ndim, std::get<1>(component_record).begin());
                std::copy(prop->comp[i]->std, prop->comp[i]->std + pmc->ndim * pmc->ndim, std::get<2>(component_record).begin());

                components << component_record;
            }

            // need to draw a different number of samples for final run
            if (status.converged)
            {
                pmc_simu_realloc(pmc, config.final_samples, err);
                pmc::check_error(err);
            }
            // reduce overall size if components died out
            else if (config.adjust_sample_size)
            {
                pmc_simu_realloc(pmc, live_components * config.samples_per_component, err);
                pmc::check_error(err);
            }
            // todo is this ever called?
            // or if sample size has been changed externally to differ from the samples before and update
            else if (unsigned(pmc->nsamples) != prop->ncomp * config.samples_per_component)
            {
                Log::instance()->message("PMC_sampler.draw_samples", ll_debug)
                    << "I'm in a surprising place";
                pmc_simu_realloc(pmc, prop->ncomp * config.samples_per_component, err);
                pmc::check_error(err);
            }

            // draw the samples
            pmc->proposal->simulate(pmc, pmc->proposal->data, rng, pmc->pb, err);
            pmc::check_error(err);

            /* dump samples */

            auto samples = file.create_data_set("/data/samples",
                PopulationMonteCarloSampler::Output::sample_type(pmc->ndim));
            auto sample_record = PopulationMonteCarloSampler::Output::sample_record(pmc->ndim);
            for (int i = 0 ; i < pmc->nsamples ; i++)
            {
                sample_record.clear();

                for (int j = 0 ; j < pmc->ndim ; ++j)
                {
                    sample_record.push_back(pmc->X[i * pmc->ndim + j]);
                }

                // index of generating component
                sample_record.push_back(pmc->indices[i]);

                // posterior not known yet
                sample_record.push_back(0);

                // weight not known yet
                sample_record.push_back(0);

                samples << sample_record;
            }

            endError(err);
        }

        void update(hdf5::File & f, const unsigned & n_samples)
        {
            pmc::ErrorHandler err;

            /* parse samples */

            auto samples_data_set = f.open_data_set("/data/samples", PopulationMonteCarloSampler::Output::sample_type(pmc->ndim));
            auto sample_record = PopulationMonteCarloSampler::Output::sample_record(pmc->ndim);
            for (unsigned i = 0 ; i < samples_data_set.records() ; ++i)
            {
                samples_data_set >> sample_record;
                std::copy(sample_record.cbegin(), sample_record.cbegin() + pmc->ndim, &pmc->X[i * pmc->ndim]);
                pmc->indices[i] = sample_record.at(pmc->ndim);
            }

            /* parse weights */

            auto weights_data_set = f.open_data_set("/data/weights", PopulationMonteCarloSampler::Output::weight_type());
            auto weight_record = PopulationMonteCarloSampler::Output::weight_record();

            auto ignores_data_set = f.open_data_set("/data/broken", PopulationMonteCarloSampler::Output::ignore_type());
            auto ignore_record = PopulationMonteCarloSampler::Output::ignore_record();

            if (n_samples != weights_data_set.records())
                throw InternalError("PMC::initialize: mismatch between size of /data/samples and /data/weights ("
                        + stringify(n_samples) + " vs " + stringify(weights_data_set.records()) + ")");
            if (n_samples != ignores_data_set.records())
                throw InternalError("PMC::initialize: mismatch between size of /data/samples and /data/broken ("
                        + stringify(n_samples) + " vs " + stringify(ignores_data_set.records()) + ")");

            for (unsigned i = 0 ; i < n_samples ; ++i)
            {
                weights_data_set >> weight_record;
                ignores_data_set >> ignore_record;

                double * x = &(pmc->X[i * pmc->ndim]);

                if (ignore_record)
                {
                    pmc->flg[i] = 0;
                    continue;
                }

                const double rloc = distribution_lkl(pmc->proposal, x, err);
                pmc::check_error(err);

                pmc->log_rho[i] = rloc;
                pmc->weights[i] = std::get<1>(weight_record);

                // check if rloc matches the difference of posterior and weight
#if 0
                const double & posterior = std::get<0>(weight_record);
                if (rloc != posterior - pmc->weights[i])
                {
                    if (std::abs((rloc - (posterior - pmc->weights[i])) / rloc) > 1e-9)
                    {
                        throw InternalError("PMC::initialize:rloc and info from file don't match for point "
                                + stringify(i) + " with ignore status " + stringify(ignore_record) + ": "
                                + stringify(rloc, 17) + " vs " + stringify(posterior - pmc->weights[i], 17));
                    }
                }
#endif
                // only values with flag == 1 are used in normalization and update
                pmc->flg[i] = 1;
            }

            // remove highest weights if desired, needs to come before weight normalization
            crop_weights();

            pmc->isLog = 1;
            pmc->maxW = * std::max_element(&pmc->weights[0], &pmc->weights[0] + n_samples);
            pmc->maxR = * std::max_element(&pmc->log_rho[0], &pmc->log_rho[0] + n_samples);

            /* update proposal */

            normalize_importance_weight(pmc, err);
            pmc::check_error(err);

            // perform the Rao-Blackwell update, including Cholesky decomposition
            pmc->pmc_update(pmc->proposal->data, pmc, err);
            pmc::check_error(err);

            status.perplexity = perplexity_and_ess(pmc, MC_NORM, &status.eff_sample_size, err);
            pmc::check_error(err);
            status.eff_sample_size /= n_samples;
            status.evidence = evidence(pmc, NULL, err);
            pmc::check_error(err);

            /* dump statistics */

            std::string subdirectory = "/data/statistics";

            auto output_file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            H5E_BEGIN_TRY
            {
                try
                {
                    f.copy(subdirectory, output_file);
                }
                // in first update no previous stats available
                catch (HDF5Error &)
                {
                    output_file.create_data_set(subdirectory, PopulationMonteCarloSampler::Output::statistics_type());
                }
            }
            H5E_END_TRY;

            // append new record to the end
            auto stats_data_set = output_file.open_data_set(subdirectory, PopulationMonteCarloSampler::Output::statistics_type());
            auto stats_record = std::make_tuple(status.perplexity, status.eff_sample_size, status.evidence);
            stats_data_set << stats_record;

            // check convergence based on all available previous values
            status.converged = check_convergence(f.name(), subdirectory);
            auto conv_attr = stats_data_set.create_or_open_attribute("converged", hdf5::Scalar<int>("converged"));
            conv_attr = status.converged;
        }
    };

    PopulationMonteCarloSampler::PopulationMonteCarloSampler(const DensityPtr & density, const hdf5::File & file,
                                                             const PopulationMonteCarloSampler::Config & config,
                                                             const bool & update) :
        PrivateImplementationPattern<PopulationMonteCarloSampler>(new Implementation<PopulationMonteCarloSampler>(density, file, config, update))
    {
    }

    PopulationMonteCarloSampler::~PopulationMonteCarloSampler()
    {
    }

    void
    PopulationMonteCarloSampler::run()
    {
        _imp->run();
    }

    const PopulationMonteCarloSampler::Config &
    PopulationMonteCarloSampler::config() const
    {
        return _imp->config;
    }

    void
    PopulationMonteCarloSampler::calculate_weights(const std::string & sample_file,
                                                     const unsigned & min_index,
                                                     const unsigned & max_index) const
    {
        _imp->calculate_weights(sample_file, min_index, max_index);
    }

    void
    PopulationMonteCarloSampler::draw_samples()
    {
        _imp->draw_samples();
    }

    void
    PopulationMonteCarloSampler::read_samples(const std::string & sample_file, const std::string & base,
                                              const unsigned & min, const unsigned & max,
                                              std::vector<std::vector<double>> & samples)
    {
        Implementation<PopulationMonteCarloSampler>::read_samples(sample_file, base, min, max, samples);
    }

    const PopulationMonteCarloSampler::Status &
    PopulationMonteCarloSampler::status() const
    {
        return _imp->status;
    }

    bool
    PopulationMonteCarloSampler::status(const PopulationMonteCarloSampler::Status & new_status, bool check_convergence)
    {
        _imp->status = new_status;
        if (check_convergence)
            return _imp->check_convergence(_imp->config.output_file);
        else
            return true;
    }

    /* PopulationMonteCarloSampler::Config */

     PopulationMonteCarloSampler::Config::Config() :
         seed(0),
         parallelize(true),
         number_of_workers(0),
         degrees_of_freedom(-1, std::numeric_limits<int>::max(), -1),
         group_by_r_value(1, std::numeric_limits<double>::max(), 1),
         patch_length(1000),
         r_value_no_nuisance(true),
         skip_initial(0, 1, 0.2),
         store_hc_initial(false),
         store_input_components(false),
         target_ncomponents(0),
         adjust_sample_size(false),
         max_updates(10),
         samples_per_component(10000),
         crop_highest_weights(0),
         need_prerun(true),
         store_prerun(true),
         convergence_eff_sample_size(0, 1, 0.80),
         convergence_perplexity(0, 1, 0.90),
         ignore_eff_sample_size(true),
         minimum_eff_sample_size(0, 1, 0.1),
         minimum_perplexity(0, 1, 0.1),
         minimum_steps(2, std::numeric_limits<unsigned>::max(), 2),
         maximum_relative_std_deviation(0, 1, 0.10),
         final_samples(20000),
         store(true),
         print_steps(0, 100, 5)
     {
     }

     PopulationMonteCarloSampler::Config
     PopulationMonteCarloSampler::Config::Default()
     {
         return PopulationMonteCarloSampler::Config();
     }

     PopulationMonteCarloSampler::Config
     PopulationMonteCarloSampler::Config::Quick()
     {
         PopulationMonteCarloSampler::Config config;

         config.samples_per_component = 3000;

         return config;
     }

     std::ostream & operator<<(std::ostream & stream, const PopulationMonteCarloSampler::Config & c)
     {
         stream << std::boolalpha
                << "Clustering options: " << std::endl
                << "critical R value = " << c.group_by_r_value
                << ", ignore groups = " << stringify_container(c.ignore_groups)
                << ", R value no nuisance = " << c.r_value_no_nuisance << std::endl
                << "sliding window = " << c.patch_length
                << ", number of components = " << c.target_ncomponents << std::endl
                << "Prerun options:" << std::endl
                << "chunk size = " << c.samples_per_component
                << ", max #updates = " << c.max_updates
                << ", adjust sample size = " << c.adjust_sample_size << std::endl
                << "degrees of freedom = " << c.degrees_of_freedom << std::endl
                << "Convergence options:" << std::endl
                << "ignore ESS = " << c.ignore_eff_sample_size
                << ", allowed std. dev = " << c.maximum_relative_std_deviation << std::endl
                << "Main run options: " << std::endl
                << "chunk size = " << c.final_samples;
         return stream;
     }

     PopulationMonteCarloSampler::Output::ComponentType
     PopulationMonteCarloSampler::Output::component_type(const unsigned & dimension)
     {
         return
         ComponentType
         {
             "component",
             hdf5::Scalar<double>("weight"),
             hdf5::Array<1, double>("mean", { dimension }),
             hdf5::Array<1, double>("covariance", { dimension * dimension })
         };
     }

     PopulationMonteCarloSampler::Output::IgnoreType
     PopulationMonteCarloSampler::Output::ignore_type()
     {
         return IgnoreType(hdf5::Scalar<short>("ignore"));
     }

     PopulationMonteCarloSampler::Output::SampleType
     PopulationMonteCarloSampler::Output::sample_type(const unsigned & dimension)
     {
         return SampleType("sample", { dimension + 3 });
     }

     PopulationMonteCarloSampler::Output::StatisticsType
     PopulationMonteCarloSampler::Output::statistics_type()
     {
         return
         StatisticsType
         {
             "statistics",
             hdf5::Scalar<double>("perplexity"),
             hdf5::Scalar<double>("effective-sample-size"),
             hdf5::Scalar<double>("evidence"),
         };
     }

     PopulationMonteCarloSampler::Output::WeightType
     PopulationMonteCarloSampler::Output::weight_type()
     {
         return
         WeightType
         {
             "weight",
             hdf5::Scalar<double>("posterior"),
             hdf5::Scalar<double>("weight"),
         };
     }

     std::tuple<double, std::vector<double>, std::vector<double>>
     PopulationMonteCarloSampler::Output::component_record(const unsigned & dimension)
     {
         return std::make_tuple(0.0, std::vector<double>(dimension), std::vector<double>(dimension * dimension));
     }

     short
     PopulationMonteCarloSampler::Output::ignore_record()
     {
         return short(0);
     }

     std::vector<double>
     PopulationMonteCarloSampler::Output::sample_record(const unsigned & dimension)
     {
         return std::vector<double>(dimension + 3);
     }

     std::tuple<double, double, double>
     PopulationMonteCarloSampler::Output::statistics_record()
     {
         return std::make_tuple(0.0, 1.0, 2.0);
     }

     std::tuple<double, double>
     PopulationMonteCarloSampler::Output::weight_record()
     {
         return std::make_tuple(0.0, 1.0);
     }

     /* PopulationMonteCarloSampler::Status */

     PopulationMonteCarloSampler::Status::Status() :
         chunk_size(1000),
         converged(false),
         iterations_at_convergence(std::numeric_limits<unsigned>::max()),
         evidence(0),
         eff_sample_size(0),
         perplexity(0)
     {
     }
}
