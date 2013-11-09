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
#include <eos/utils/population_monte_carlo_sampler.hh>
#include <eos/utils/cluster.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/hierarchical-clustering.hh>
#include <eos/utils/log.hh>
#include <eos/utils/markov_chain_sampler.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/proposal_functions.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/rvalue.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/thread_pool.hh>
#include <eos/utils/welford.hh>

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
            Analysis * ana = (Analysis *) data;

            // set parameters
            unsigned j = 0;
            for (auto i = ana->parameter_descriptions().begin(), i_end = ana->parameter_descriptions().end() ; i != i_end ; ++i, ++j)
            {
                // direct assignment doesn't work. Why?
                // population_monte_carlo_sampler.cc:40: error: no match for ‘operator=’ in ‘i.__gnu_cxx::__normal_iterator<_Iterator, _Container>::operator-> [with _Iterator = const eos::ParameterDescription*, _Container = std::vector<eos::ParameterDescription, std::allocator<eos::ParameterDescription> >]()->eos::ParameterDescription::parameter = 2.29999999999999982236431605997495353221893310547e+0’
                //                    i->parameter = par_point[j];
                Parameter p = i->parameter;
                p = par_point[j];
            }
            const double post = ana->log_posterior();
            if ( ! std::isfinite(post))
                throw InternalError("PMC::posterior: not finite " + stringify(post)
                                    + " at " + stringify(par_point, par_point + ana->parameter_descriptions().size()) );
            return post;
        }

      // todo remove
        double* mvdens_ran_extreme(double* dest, mvdens * g, gsl_rng * r, error** err, double offset) {
            size_t i;
            double *res;
            gsl_vector * res_view;
            gsl_vector_view res_view_container;
            double val,corr,u;

            mvdens_cholesky_decomp(g,err);
            pmc::check_error(err);

            MALLOC_IF_NEEDED(res,dest,sizeof(double)*g->ndim,err);
            pmc::check_error(err);

            gsl_set_error_handler_off();
            /* Generate a N(O,I) distributed vector               */
            for (i = 0; i < g->ndim; i++) {
                val= -offset + 2.0 * offset * gsl_rng_uniform(r);
                res[i]=val;
            }

            /* Make it N(O,Sigma)                         */
            res_view_container=gsl_vector_view_array(res,g->ndim);
            res_view=&(res_view_container.vector);
            gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, g->std_view, res_view);

            corr=1;
            if(g->df!=-1) {
                // make it student-t
                u = gsl_ran_chisq(r, g->df);
                corr=sqrt(g->df/u);
            }

            // add mean and correct if student-t
            for(i=0;i<g->ndim;i++) {
                res[i]*=corr;
                res[i]+=g->mean[i];
            }
            return res;
        }

        typedef std::pair<unsigned, double> index_pair;
        //todo replace by lambda
        bool
        cmp(const index_pair & a, const index_pair & b )
        {
            return a.second > b.second;
        }

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
           std::shared_ptr<Analysis> analysis;

           // store the posterior values
           std::vector<double> posterior_values;

           // points at which posterior is evaluated
           std::vector<double> parameter_samples;

           std::shared_ptr<ROOT::Minuit2::FunctionMinimum> minimum;

           Worker(const Analysis & analysis) :
               analysis(analysis.clone())
           {
           }

           void clear()
           {
               parameter_samples.resize(0);
               posterior_values.resize(0);
           }

           // call from main thread before actual work is done
           void setup(double * samples, unsigned n_samples, unsigned n_dim)
           {
               // copy the samples
               parameter_samples = std::vector<double>();
               std::copy(samples, samples + n_samples * n_dim, std::back_inserter(parameter_samples));

               posterior_values.resize(n_samples);
           }

           // compute log(posterior) at many sample points
           void work()
           {
               pmc::ErrorHandler err;

               if (parameter_samples.empty() || posterior_values.empty())
                   return;

               unsigned n_dim = parameter_samples.size() / posterior_values.size();
               unsigned i = 0;
               for (auto p = posterior_values.begin(), p_end = posterior_values.end(); p != p_end ; ++i, ++p)
               {
                    *p = pmc::logpdf(analysis.get(), &parameter_samples[i * n_dim], err);
               }
           }

           std::vector<double> mode() const
           {
               std::vector<double> result;

               if (! (minimum && minimum->IsValid()) )
               {
                   return result;
               }

               for (unsigned i = 0 ; i < analysis->parameter_descriptions().size() ; ++i)
               {
                   result.push_back(minimum->UserParameters().Value(i));
               }

               return result;
           }

           void optimize(std::vector<double> initial_point)
           {
               minimum.reset(new ROOT::Minuit2::FunctionMinimum(analysis->optimize_minuit(initial_point, Analysis::OptimizationOptions::Defaults())));
           }
       };
    }

    template<>
    struct Implementation<PopulationMonteCarloSampler>
    {
        typedef std::vector<HistoryPtr> ChainGroup;
        typedef std::vector<unsigned> IndexList;

        // store reference, but don't own analysis
        Analysis analysis;

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

        Implementation(const Analysis & analysis, const PopulationMonteCarloSampler::Config & config) :
            analysis(analysis),
            config(config),
            status(),
            pmc(NULL)
        {
            // setup Mersenne-Twister RN generator using custom seed
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, config.seed);

            if (config.component_weights.empty())
                throw InternalError("PMC_sampler.ctor: No weights for components specified");

            setup_output();

            // setup PMC library
            initialize_pmc();

            // initialization of the workers
            const unsigned number_of_workers = config.number_of_workers == 0 ?
                                               ThreadPool::instance()->number_of_threads() :
                                               config.number_of_workers;
            for (unsigned i = 0; i < number_of_workers; ++i)
                workers.push_back(std::make_shared<pmc::Worker>(analysis));

            auto f = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);
            this->analysis.dump_descriptions(f, "descriptions");
            dump("initial");
        }

        Implementation(const Analysis & analysis, const hdf5::File & file,
                       const PopulationMonteCarloSampler::Config & config, const bool & update) :
            analysis(analysis),
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
                workers.push_back(std::make_shared<pmc::Worker>(analysis));

            auto f = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);
            this->analysis.dump_descriptions(f, "descriptions");
            dump("initial");
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

            unsigned n_dim = analysis.parameter_descriptions().size();

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
                std::move((**w).posterior_values.begin(), (**w).posterior_values.end(), std::back_inserter(posterior_values));
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
            std::sort(weight_indices.begin(), weight_indices.end(), &pmc::cmp);

            // set OK flag to false for highest weights
            for (unsigned j = 0 ; j < config.crop_highest_weights ; ++j)
            {
                pmc->flg[weight_indices[j].first] = 0;
            }
        }

        /*!
         * Dump status to HDF5.
         *
         * @param group
         * @param samples If false, only summary statistics are stored
         */
        void dump(const std::string & group, bool store_samples = true)
        {
            const unsigned & dim = analysis.parameter_descriptions().size();

            // open the file whenever writing is desired, so it is in a readable state during lengthy posterior calculations
            auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            /* dump components */

            // proposal density
            mix_mvdens* prop = static_cast<mix_mvdens *>(pmc->proposal->data);

            auto components = file.create_data_set("/data/" + group + "/components",
                PopulationMonteCarloSampler::Output::component_type(dim));
            auto component_record = PopulationMonteCarloSampler::Output::component_record(dim);
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

        /*
         * Filter out components which don't overlap with the integration domain.
         *
         * Draw samples from a single component. If too few inside the domain,
         * set the component's weight to zero
         * @return The number of live components.
         * @note mind the flag mix_mvdens::init_cwght (cumulative weight?) when calling this method
         */
        unsigned filter_components(mix_mvdens * & mmv, const parabox * pb)
        {
            pmc::ErrorHandler err;

            std::vector<double> point(mmv->ndim);

            unsigned n_live_components = 0;
            for (unsigned i = 0 ; i < mmv->ncomp ; ++i)
            {
                unsigned n_in = 0;
                for (unsigned j = 0 ; j < config.chunk_size ; ++j)
                {
                    mvdens_ran(&point[0], mmv->comp[i], rng, err);
                    if (isinBox(pb, &point[0], err))
                        ++n_in;
                }

                // check if enough were in
                const double ratio = 1.0 * n_in / config.chunk_size;
                if (ratio < config.minimum_overlap)
                {
                    mmv->wght[i] = 0.0;
                }
                else
                {
                    ++n_live_components;
                }
            }

            Log::instance()->message("PMC::filter_components", ll_informational)
                << mmv->ncomp - n_live_components << " components were removed "
                << "by minimum overlap of " << config.minimum_overlap;

            if (n_live_components == mmv->ncomp)
            {
                return n_live_components;
            }

            if (n_live_components == 0)
                throw InternalError("PMC::filter_components: removed all components. Check parameter ranges!");

            // return new proposal with components removed
            mix_mvdens * mmv_old = mmv;
            mmv = mix_mvdens_alloc(n_live_components, mmv_old->ndim, err);

            unsigned i_new = 0;
            for (unsigned i_old = 0 ; i_old < mmv_old->ncomp ; ++i_old)
            {
                if (mmv_old->wght[i_old] == 0.0)
                {
                    // leads to  munmap_chunk(): invalid pointer
                    //                    auto mv = &mmv_old->comp[i_old];
                    //                    mvdens_free(mv);
                    mvdens_empty(mmv_old->comp[i_old]);
                    continue;
                }

                // equal weight to every cluster
                mmv->wght[i_new] = 1.0 / n_live_components;
                mmv->comp[i_new] = mmv_old->comp[i_old];

                ++i_new;
            }

            return n_live_components;
        }

        void initialize_pmc(const hdf5::File & file, const bool & update)
        {
            Log::instance()->message("PMC_sampler::initialize", ll_informational)
                << "Reading from file " << file.name();

            // we will only use the file for reading.
            hdf5::File & f = const_cast<hdf5::File &>(file);

            // number of dimensions of parameter cube
            int n_dim = int(analysis.parameter_descriptions().size());

            // exception handler look-alike
            pmc::ErrorHandler err;

            // parameter cube: copy from analysis
            parabox * par_box = init_parabox(n_dim, err);
            int i = 0;
            for (auto d = analysis.parameter_descriptions().begin(), d_end = analysis.parameter_descriptions().end() ; d != d_end ; ++d, ++i)
            {
                add_slab(par_box, i, d->min, d->max, err);
            }

            /* setup importance sampling */

            // create posterior distribution
            distribution * target = init_simple_distribution(n_dim, (void *) &analysis, &pmc::logpdf, NULL, err);
            pmc::check_error(err);

            /* create proposal density from file */

            // use functions in pmc_rc.c:
            // o rcinit_mix_mvdens(...)
            // o init_importance_from_rc(...)
            mix_mvdens * mmv = NULL;

            /* determine type of input file */

            unsigned number_of_live_components = 0;

            // read from a global local proposal function (actually for MCMC)
            if (f.group_exists("/global local"))
            {
                std::string directory_base_name = "/global local";
                // read in proposal
                ProposalFunctionPtr prop;
                try
                {
                    prop = Factory::make(f, directory_base_name, "GlobalLocal", analysis.parameter_descriptions().size());
                }
                catch (HDF5Error & e)
                {
                    Log::instance()->message("population_monte_carlo_sampler.initialize", ll_error)
                        << "Errors in reading from the HDF5 file can be due to a mismatch in the "
                        << "analysis definition. Check that the same number of parameters is defined now "
                        << "and when building the GlobalLocal proposal function";
                    throw e;
                }
                auto gl = dynamic_cast<GlobalLocal *>(prop.get());
                if ( ! gl)
                    throw InternalError("population_monte_carlo_sampler::initialize: couldn't read GlobalLocal from disk");

                // get components from global local. Initially, all components are alive
                const unsigned n_clusters = config.single_cluster > -1 ? 1 : gl->component_probabilities().size();
                number_of_live_components = n_clusters * config.components_per_cluster;
                mmv = mix_mvdens_alloc(number_of_live_components, unsigned(n_dim), err);

                // copy weights: each component of a cluster gets same weight
                // NB: weights add up to one
                if (config.single_cluster > -1)
                {
                    Log::instance()->message("population_monte_carlo_sampler.initialize", ll_debug)
		                << "Using single component: " << config.single_cluster;
                    std::vector<double> weights(config.components_per_cluster, 1.0 / config.components_per_cluster);
                    std::copy(weights.cbegin(), weights.cend(), mmv->wght);
                }
                else
                {
                	// loop over gl clusters and
                    for (auto w = gl->component_probabilities().cbegin() ; w != gl->component_probabilities().cend() ; ++w)
                    {
                        std::vector<double> weights(config.components_per_cluster, *w / config.components_per_cluster);
                        std::copy(weights.cbegin(), weights.cend(),
                                  mmv->wght + std::distance(gl->component_probabilities().cbegin(), w) * config.components_per_cluster);
                    }
                }

                if (std::abs(std::accumulate(mmv->wght, mmv->wght + mmv->ncomp, 0.0) - 1.0) > 1e-13)
                    throw InternalError("Could adjust component weights to sum up to one, the weights are " + stringify(mmv->wght, mmv->wght + mmv->ncomp, 4));

                /* initialize positions and covariances of individual components */

                if (config.components_per_cluster > gl->history_states().front().size())
                    throw InternalError("PMC_sampler::initialize: mismatch between desired number of components"
                                        " per cluster (" + stringify(config.components_per_cluster) +
                                        ") and available history points (" + stringify(gl->history_states().size()));

                // loop over clusters ( or just the single cluster)
                unsigned cl_i = config.single_cluster > -1 ? config.single_cluster : 0;
                auto cl_begin = gl->history_states().cbegin() + cl_i;
                auto cl = cl_begin;
                auto cl_end = config.single_cluster > -1 ? cl + 1 : gl->history_states().cend();
                for ( ; cl !=  cl_end; ++cl, ++cl_i)
                {
                    /* find target covariance within cluster */

                    std::string sub_directory = directory_base_name + "/local proposals/" + stringify(cl_i);
                    // need to find out the right type
                    auto meta_data_set = f.open_data_set(sub_directory + "/meta", meta_type());
                    auto meta_record = proposal_functions::meta_record();
                    meta_data_set >> meta_record;

                    if (std::get<1>(meta_record) != mmv->ndim)
                    {
                        throw InternalError("PMC_sampler::initialize: current dimension(" + stringify(mmv->ndim) +
                                            ") doesn't match that in proposal (" + stringify(std::get<1>(meta_record)) + ").");
                    }

                    // default case: MultivariateGaussian
                    const bool proposal_type_student = (std::get<0>(meta_record) == std::string("MultivariateStudentT"));

                    // use the factory again with the name of the local proposal type
                    MultivariateProposalPtr local_prop = MultivariateAccess::access(
                            Factory::make(f, sub_directory, std::get<0>(meta_record), mmv->ndim));

                    if (! local_prop)
                        throw InternalError("PMC_sampler::initialize: Couldn't cast the Multivariate proposal, found type '" + std::string(std::get<0>(meta_record)) + "'");

                    // default degree of freedom value for Gaussian
                    int dof = -1;

                    if (proposal_type_student)
                    {
                        try
                        {
                            dof = dynamic_cast<MultivariateStudentT &>(*local_prop.get()).dof;
                        }
                        catch (std::bad_cast &)
                        {
                            throw InternalError("PMC_sampler::initialize: Couldn't cast the student T proposal");
                        }
                    }

                    // choose points randomly
                    std::vector<unsigned> history_point_indices(cl->size());
                    std::iota(history_point_indices.begin(), history_point_indices.end(), 0);
                    if (config.random_start)
                        std::random_shuffle(history_point_indices.begin(), history_point_indices.end());

                    /* Extract the target covariance:
                     * a) local covariance for each history point
                     * b) one covariance for each point in the cluster
                     */

                    // need to apply inverse of scale factor to obtain the estimate of the target variance
                    local_prop->rescale(1.0 / local_prop->covariance_scale);
                    const gsl_matrix * covariance = local_prop->covariance();

                    // loop over states within cluster: each gets different position, but same covariance

                    for(unsigned s_j = 0; s_j < config.components_per_cluster ; ++s_j)
                    {
                        mvdens* mv = mmv->comp[std::distance(cl_begin, cl) * config.components_per_cluster + s_j];
                        auto point = cl->at(history_point_indices[s_j]).point;
                        std::copy(point.cbegin(), point.cend(), mv->mean);
                        if ( ! gl->local_covariances().empty())
                        {
                            auto local_cov = gl->local_covariances().at(cl_i).at(history_point_indices[s_j]);
                            std::copy(local_cov.cbegin(), local_cov.cend(), mv->std);
                        }
                        else
                        {
                            std::copy(covariance->data, covariance->data + mmv->ndim * mmv->ndim, mv->std);
                            if (covariance->size1 != mmv->ndim || covariance->size2 != mmv->ndim)
                                throw InternalError("Cannot initialize components which are not in the multivariate block");
                        }
                        // cholesky not yet computed
                        mv->chol = 0;
                        mv->band_limit = n_dim;
                        mv->df = config.override_global_local_proposal ? config.degrees_of_freedom : dof;
                        // we can ignore mv->detL, it is computed when chol. decomp. is performed
                    }
                }
#if 0
                std::string component_means;
                for (unsigned i = 0 ; i < mmv->ncomp ; ++i)
                {
                    component_means += stringify(mmv->comp[i]->mean, mmv->comp[i]->mean + std::min(4u, unsigned(mmv->ndim))) + '\n';
                }
#endif
                Log::instance()->message("PMC_sampler::initialize", ll_informational)
                    << "Initialized " << mmv->ncomp << " components"
                    << " with proposal dof " << mmv->comp[0]->df;
            }
            else if (config.super_clusters > 0)
            {
                // initialize mmv by clustering and filtering
                number_of_live_components = hierarchical_clustering(f, mmv, par_box);
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

            if (mmv->ndim != analysis.parameter_descriptions().size())
                throw InternalError("PMC::ctor: mismatch of parameter dimensions of analysis vs proposal ("\
                                    + stringify(analysis.parameter_descriptions().size())
                                    + " vs " + stringify(mmv->ndim) + ")");

            if (config.block_decomposition)
            {
                unsigned par_i = 0;
                for (auto par = analysis.parameter_descriptions().begin(), d_end = analysis.parameter_descriptions().end() ; par != d_end ; ++par, ++par_i)
                {
                    if (! par->nuisance)
                        continue;

                    auto prior = analysis.log_prior(par->parameter.name());

                    // todo dirty hack
                    std::string s = prior->as_string();
                    if (s.find("flat") != std::string::npos)
                        continue;

                    // set row/column with this parameter to zero in all components
                    // only set diagonal element of variance and change mean
                    for (unsigned c = 0 ; c < mmv->ncomp ; ++c)
                    {
                        gsl_vector_set(mmv->comp[c]->mean_view, par_i, prior->mean());

                        auto row = gsl_matrix_row(mmv->comp[c]->std_view, par_i);
                        gsl_vector_set_all(&row.vector, 0.0);
                        auto column = gsl_matrix_column(mmv->comp[c]->std_view, par_i);
                        gsl_vector_set_all(&column.vector, 0.0);
                        gsl_matrix_set(mmv->comp[c]->std_view, par_i, par_i, prior->variance());
                    }
                }
            }

            // create proposal distribution object
            distribution * proposal = mix_mvdens_distribution(mmv->ndim, (void*) mmv, err);

            // number of samples per chunk. Fixed size for each live component.
            unsigned n_samples = config.chunk_size * number_of_live_components;

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

        // Following the lines of pmclib/pmc_rc.c
        void initialize_pmc()
        {
            Log::instance()->message("PMC_sampler::initialize", ll_informational)
                << " Using random points / Minuit for initialization";

            // number of dimensions of parameter cube
            int n_dim = int(analysis.parameter_descriptions().size());

            // setup importance sampling
            pmc::ErrorHandler err;

            // create posterior distribution
            distribution * target = init_simple_distribution(n_dim, (void *) &analysis, &pmc::logpdf, NULL, err);
            pmc::check_error(err);

            /*  create proposal density */

            // use functions in pmc_rc.c:
            // o rcinit_mix_mvdens(...)
            // o init_importance_from_rc(...)
            mix_mvdens * mmv = NULL;

            // randomize under first component
            if (config.random_start)
            {
                mmv = mix_mvdens_alloc(config.component_weights.size(), unsigned(n_dim), err);
                mvdens* mv = mmv->comp[0];

                // copy weights
                std::copy(config.component_weights.cbegin(), config.component_weights.cend(), mmv->wght);

                // is first mean and variance specified?
                if ( ! (config.component_means.size() == 1 && config.component_variances.size() == 1) )
                {
                    // mean = center or parameter cube
                    std::vector<double> mean(n_dim, 0.0);
                    for (int i = 0; i < n_dim ; ++i)
                    {
                        mean[i] = (analysis.parameter_descriptions()[i].max + analysis.parameter_descriptions()[i].min) / 2.0;
                    }
                    config.component_means.clear();
                    config.component_means.push_back(mean);

                    // variance = (parameter range / 2)^2 (diagonal matrix)
                    std::vector<double> var(n_dim, 0.0);
                    for (int i = 0; i < n_dim ; ++i)
                    {
                        var[i] = (analysis.parameter_descriptions()[i].max - analysis.parameter_descriptions()[i].min) /
                                  config.std_dev_reduction;
                        var[i] *= var[i];
                    }
                    config.component_variances.clear();
                    config.component_variances.push_back(var);
                }

                if ( ! (config.component_means.front().size() == unsigned(n_dim) &&
                        config.component_variances.front().size() == unsigned(n_dim)) )
                    throw InternalError("PMC_sampler.cc: Need to specify mean and variance for exactly one component in all " + stringify(n_dim) + " dimensions.");

                // copy means if they are all specified
                std::copy(config.component_means.front().cbegin(), config.component_means.front().cend(), mv->mean);

                /* create covariance matrix */

                // case 1: only diagonal specified. Set off-diagonal = zero
                if (config.component_variances.front().size() == unsigned(n_dim) )
                {
                    for (int i = 0; i < n_dim ; ++i)
                    {
                        for (int j = i ; j < n_dim ; ++j)
                        {
                            if ( i == j)
                                mv->std[i*n_dim + j] = config.component_variances.front()[i];
                            else
                                mv->std[i*n_dim + j] = 0.0;
                        }
                    }
                }
                // case 2: full matrix specified
                else if (config.component_variances.front().size() == unsigned(n_dim * n_dim) )
                {
                    std::copy(config.component_variances.front().cbegin(), config.component_variances.front().cend(), mv->std);
                }
                // case 3: misspecified
                else
                {
                    throw InternalError("PMC_sampler.cc: Covariance matrix doesn't have right dimensions.");
                }

                // chol. and determinant
                /*
                 *  If set to 0, PMC will use GSL to compute the cholesky decomposition
                 *  as soon as samples are drawn.
                 *  and the result is stored in the same memory block, thus the
                 *  interpretation of mmv->comp[0]->std is changed
                 */
                mv->chol = 0;

                mv->band_limit = n_dim;
                mv->df = config.degrees_of_freedom;

                // specify additional proposal density components
                for (unsigned i = 1 ; i < config.component_weights.size() ; ++i)
                {
                    mvdens_ran(mmv->comp[i]->mean, mmv->comp[0], rng, err);
                    mmv->comp[i]->chol = mmv->comp[0]->chol;
                    mmv->comp[i]->band_limit = mmv->comp[0]->band_limit;
                    mmv->comp[i]->df = mmv->comp[0]->df;
                    mmv->comp[i]->detL = mmv->comp[0]->detL;
                    std::copy(mmv->comp[0]->std, mmv->comp[0]->std + n_dim * n_dim, mmv->comp[i]->std);
                }
            }
            else
            {
                /*
                 * 1. User defines number of starting positions, either explicitly by giving the points.
                 *    Number can be different from #components.
                 * or implicitly by the number. Latter case: draw points from the priors.
                 * 2. Use MINUIT to find local modes and covariance matrices.
                 * 3. Estimate weight of one mode by max. posterior * sqrt(determinant of covariance)
                 * 4. Find unique modes: |mode_i - mode_j| < eps => i and j are same modes
                 * 5. Initialize components: at least two per mode,
                 *    for nuisance parameters, leave position at mode.
                 *    for pars. of interest, draw a sample from multivariate gaussian, and ignore the nuisance values.
                 *    This way, components can distribute along a banana shaped posterior.
                 * 6. Perhaps one could estimate the non-Gaussianity (DoF < inf) from comparing the
                 *    Hessian with std errors. If coincide closely, DoF = inf is good, else chose value close to 0
                 */

                std::vector<double> starting_point(n_dim, 0.0);

                std::vector<Ticket> tickets;
                std::vector<std::shared_ptr<pmc::Worker>> optimizers;

                // draw starting point from priors
                for (unsigned i = 0 ; i < config.starting_points ; ++i)
                {
                    auto s = starting_point.begin();
                    for (auto d = analysis.parameter_descriptions().begin(), d_end = analysis.parameter_descriptions().end() ; d != d_end ; ++d, ++s)
                    {
                        if (d->nuisance)
                        {
                            *s = d->parameter();
                        }
                        else
                        {
                            //todo not ready for multidimensional priors
                            LogPriorPtr prior = analysis.log_prior(d->parameter.name());
                            *s = prior->sample(rng);
                        }
                    }

                    // parallelize computations
                    auto worker = std::make_shared<pmc::Worker>(analysis);
                    optimizers.push_back(worker);

                    if (config.parallelize)
                    {
                        tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&pmc::Worker::optimize, worker.get(), starting_point)));
                    }
                    else
                    {
                        worker->optimize(starting_point);
                    }
                }

                // wait for job completion
                for (auto t = tickets.begin(), t_end = tickets.end() ; t != t_end ; ++t)
                    t->wait();

                /* find unique modes */
                std::vector<std::vector<double>> unique_modes;
                std::vector<unsigned> worker_indices;

                // add first mode
                auto o = optimizers.begin();

                for (auto o_end = optimizers.end() ; o != o_end ; ++o)
                {
                    auto mode = (**o).mode();
                    if( mode.empty() )
                    {
                        Log::instance()->message("PMC_sampler.initialize", ll_warning)
                            << "worker couldn't find mode using minuit";
                        continue;
                    }

                    // first valid mode
                    if (unique_modes.empty())
                    {
                        unique_modes.push_back(mode);
                        worker_indices.push_back(0);
                        continue;
                    }

                    // compare with every predecessor
                    // by calculating the Euclidean norm between this mode and
                    // previously found unique modes, rescaled by allowed range
                    double difference, length1, length2;
                    difference = length1 = length2 = 0;
                    auto m = unique_modes.begin(), m_end = unique_modes.end();
                    for ( ; m != m_end ; ++m, difference = 0)
                    {
                        auto j = mode.begin();
                        auto d = analysis.parameter_descriptions().begin();
                        for ( auto i = m->begin(), i_end = m->end() ; i != i_end ; ++i, ++j, ++d)
                        {
                            double range = d->max - d->min;
                            difference += power_of<2>((*i - *j)  / range);
                            length1 += power_of<2>((*j - d->min) / range);
                            length2 += power_of<2>((*i - d->min) / range);
                        }

                        // rescale to account for dimensionality
                        // maximum distance of two points in n-dim unit hypercube is sqrt(n)
                        difference /= mode.size();
                        length1 /= mode.size();
                        length2 /= mode.size();

                        Log::instance()->message("PMC_sampler.initialize", ll_debug)
                            << "Length1=" << length1 << ", length2=" << length2
                            << ", diff=" << difference;
                        // check that one vector is not of length zero
                        if (std::sqrt(difference) < config.mode_distance)
                        {
                            break;
                        }
                    }
                    // made it to the end?
                    if ( m == m_end)
                    {
                        // so we can use the posterior info there after the loop
                        // for ranking the modes
                        unique_modes.push_back(mode);
                        worker_indices.push_back(std::distance(optimizers.begin(), o));
                    }

//                    std::string output("(");
//                    for (auto i = mode.begin(), i_end = mode.end() ; i != i_end ; ++i)
//                        output += " " + stringify(*i);
//                    output += " )";
                    Log::instance()->message("PMC_sampler.initialize", ll_debug)
                                << "Found mode: " << stringify(mode.begin(), mode.end(), 4);
                }

                // print out modes
                Log::instance()->message("PMC_sampler.initialize", ll_informational)
                    << "Identified " << unique_modes.size() << " unique mode(s) of posterior.";
                for (auto m = unique_modes.begin(), m_end = unique_modes.end() ; m != m_end ; ++m)
                {
                    Log::instance()->message("PMC_sampler.initialize", ll_debug)
                        << stringify(m->begin(), m->end());
                }

                // estimate weight of each mode by posterior
//                std::vector<double> weight_estimates

                // issue warning if there appear to be not enough components

                    // free to change number of components
                unsigned m = config.components_per_cluster;
                mmv = mix_mvdens_alloc(config.components_per_cluster * unique_modes.size(), n_dim, err);
//                }

                // a temp density, needed to initialize the actual components
                mvdens * guide = mvdens_alloc(n_dim, err);

                // create component from minuit result
                for (unsigned i = 0 ; i < mmv->ncomp ; ++i)
                {
                    // equal weight for each component
                    mmv->wght[i] = 1.0 / double(mmv->ncomp);

                    // set mean to the unique mode
                    std::copy(unique_modes[i / m].begin(), unique_modes[i / m].end(), guide->mean);

                    // rescale volume/determinant
                    double scale = std::pow(config.std_dev_reduction, 1.0 / double(n_dim));

                    ROOT::Minuit2::FunctionMinimum min = *optimizers[worker_indices[i / m]]->minimum;
                    Log::instance()->message("pmc-sampler.minimum", ll_debug)
                        << "comp #" << i << ": min = " << min.UserCovariance();

                    // copy covariance matrix
                    for (int j = 0 ; j < n_dim ; ++j)
                    {
                        for (int k = 0 ; k < n_dim ; ++k)
                        {
                            auto w = optimizers[worker_indices[i / m]];
                            mmv->comp[i]->std[j * n_dim + k ] = w->minimum->UserCovariance()(j, k);
                            guide->std[j * n_dim + k ] = scale * w->minimum->UserCovariance()(j, k);
                        }
                    }

                    // sample a mean from the distribution
                    /* but not from first comp. of this mode, since overwriting itself lead to problem:
                     * 'Cholesky decomposition failed' when drawing samples
                     */

                    // sample mean of component around the mode
                    // this way components have a better chance
                    pmc::mvdens_ran_extreme(mmv->comp[i]->mean, guide, rng, err, config.component_offset);
                    pmc::check_error(err);

                    Log::instance()->message("PMC_sampler.initialize", ll_debug)
                        << "comp initialized to " << stringify(mmv->comp[i]->mean, mmv->comp[i]->mean + n_dim);

                    // copy remaining parameters
                    mmv->comp[i]->chol = 0;
                    mmv->comp[i]->band_limit = n_dim;
                    mmv->comp[i]->df = config.degrees_of_freedom;
                    mmv->comp[i]->detL = mmv->comp[0]->detL;
                }

                mvdens_free(&guide);
            }

            // create proposal distribution object
            distribution * proposal = mix_mvdens_distribution(mmv->ndim, (void*) mmv,err);

            // number of samples per chunk. Fixed size for each component.
            long n_samples = config.chunk_size * mmv->ncomp;

            // setup with 0 (?) deduced parameters for now
            pmc = pmc_simu_init_plus_ded(n_samples, target->ndim, target->n_ded, err);

            // parameter cube: copy from analysis
            parabox * par_box = init_parabox(n_dim, err);
            int i = 0;
            for (auto d = analysis.parameter_descriptions().begin(), d_end = analysis.parameter_descriptions().end() ; d != d_end ; ++d, ++i)
            {
                add_slab(par_box, i, d->min, d->max, err);
            }

            // setter methods, include safety checks of variable
            pmc_simu_init_target(pmc, target, par_box, err);
            pmc_simu_init_proposal(pmc, proposal, config.print_steps, err);
            pmc_simu_init_pmc(pmc, NULL, NULL, update_prop_rb_void, err);
            pmc::check_error(err);
        }


        void group_chains(std::vector<ChainGroup> & chain_groups)
        {
            std::list<HistoryPtr> available_chains(chain_groups.front().cbegin(), chain_groups.front().cend());

            // for historical reasons, the R-value grouping was called 'clustering'
            Cluster::RValueFunction r = &RValue::approximation;

            // create indices of parameters whose R-value is checked
            std::vector<unsigned> parameter_indices;
            {
                unsigned i = 0;
                for (auto d = analysis.parameter_descriptions().cbegin() ; d != analysis.parameter_descriptions().cend() ; ++d, ++i)
                {
                    if ( d->nuisance && config.r_value_no_nuisance)
                        continue;

                    parameter_indices.push_back(i);
                }
            }

            // add first chain as first group
            std::vector<Cluster> groups{ Cluster(r, config.group_by_r_value,
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
                    groups.push_back(Cluster(r, config.group_by_r_value,
                            available_chains.front(), chain_index, config.skip_initial));
                    groups.back().parameter_indices(parameter_indices);

                    Log::instance()->message("PMC.hierarchical_clustering", ll_debug)
                                    << "Created new group for chain " << chain_index;
                }
                available_chains.pop_front();
            }

            /* copy the groups */

            chain_groups.clear();

            // store lengths of groups for output onlye
            std::vector<unsigned> sizes_groups;
            for (auto g = groups.begin() ; g != groups.end() ; ++g)
            {
                chain_groups.push_back(ChainGroup());
                unsigned group_size = 0;
                for (auto c = g->begin() ; c != g->end() ; ++c, ++group_size)
                {
                    chain_groups.back().push_back(*c);
                }
                sizes_groups.push_back(group_size);
            }

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                            << "Found " << sizes_groups.size() << " groups of chains with " << stringify_container(sizes_groups)
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
                    if (*i >= chain_groups.size())
                    {
                        Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_error)
                            << "Skipping invalid ignore group: " << *i;
                        continue;
                    }
                    chain_groups.erase(chain_groups.begin() + *i);
                    sizes_groups.erase(sizes_groups.begin() + *i);

                    Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_warning)
                                    << "Removing group " << *i;
                }

                Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                                << "Using " << sizes_groups.size() << " groups of chains with " << stringify_container(sizes_groups)
                                << " members";
            }
        }

        /*
         * Finding the initial cluster guess:
         * 1. Respect chains: choose one cluster only from patches within a chain
         * 2. Fix desired #clusters and #samples for a patch
         */
        unsigned hierarchical_clustering(const hdf5::File & file, mix_mvdens * & mmv, const parabox * pb)
        {
            /* parse chain histories */

            ChainGroup chains;
            {
                std::vector<std::shared_ptr<hdf5::File>> input_files;
                input_files.push_back(std::make_shared<hdf5::File>(hdf5::File::Open(file.name(), H5F_ACC_RDONLY)));
                auto gl_config = proposal_functions::GlobalLocal::Config::Default();
                gl_config.skip_initial = config.skip_initial;
                chains = MarkovChainSampler::build_global_local("", input_files, gl_config, analysis.clone());
            }

            if (chains.front()->states.front().point.size() != analysis.parameter_descriptions().size())
            {
                Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_warning)
                    << "The analysis in MCMC prerun had dim " << chains.front()->states.front().point.size()
                    << ", but now the analysis has dim " << analysis.parameter_descriptions().size();
            }

            HierarchicalClustering::Config conf = HierarchicalClustering::Config::Default();
            conf.equal_weights = true;
            HierarchicalClustering hc(conf);

            /* group chains according to R-value */

            std::vector<ChainGroup> chain_groups{ chains };

            if (config.group_by_r_value > 1)
            {
                group_chains(chain_groups);
            }

            /* create initial guess for super clusters by drawing local patches uniformly w/o replacement or large windows */

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Creating initial guess for the " << config.super_clusters
                << " clusters to be formed from large windows"
                << (config.group_by_r_value > 1 ? " for each of the " + stringify(chain_groups.size()) + " chain groups found" : "");

            HierarchicalClustering::MixtureDensity initial_clusters;

            // weight of a single cluster (weights sum up to one)
            const unsigned n_clusters_total = config.super_clusters * chain_groups.size();
            const double weight = 1.0 / n_clusters_total;

            {
                for (auto g = chain_groups.cbegin() ; g != chain_groups.cend() ; ++g)
                {
                    // how many clusters should each each in group contribute
                    std::vector<unsigned> super_clusters_per_chain;
                    pmc::minimal_partition(config.super_clusters, g->size(), super_clusters_per_chain);

                    auto n_clusters = super_clusters_per_chain.cbegin();
                    for (auto c = g->cbegin() ; c != g->cend() ; ++c, ++n_clusters)
                    {
                        //todo don't throw away clusters. Think about it!
                        if (! *n_clusters)
                            continue;

                        auto first_state = (**c).states.cbegin() + config.skip_initial * (**c).states.size();
                        const int window = std::distance(first_state, (**c).states.cend()) / (*n_clusters);
                        if (window < 0)
                        {
                            throw InternalError("PMC::hierarchical_clustering: number of super clusters too large for history size and skip initial: " + stringify(window)
                                    + " vs " + stringify(std::distance(first_state, (**c).states.cend())) + " and " + stringify(config.skip_initial));
                        }
                        auto last_state = first_state + window;

                        bool done = false;
                        while (! done)
                        {
                            // add remainder to last super_cluster
                            if (std::distance(last_state, (**c).states.cend()) < window)
                            {
                                last_state = (**c).states.cend();
                                done = true;
                            }

                            std::vector<double> mean, covariance;
                            (**c).mean_and_covariance(first_state, last_state, mean, covariance);
                            std::vector<double> center = mean;

                            if (config.patch_around_local_mode)
                            {
                                center = (**c).local_mode(first_state, last_state).point;
                            }

                            HierarchicalClustering::Component super_cluster(center, covariance, weight);

                            initial_clusters.push_back(super_cluster);

                            // update range for next iteration
                            first_state += window;
                            last_state += window;
                        }
                    }
                }
            }

            hc.initial_guess(initial_clusters);

            /* create patches from each chain */

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Creating patches of length " << config.sliding_window;

            HierarchicalClustering::MixtureDensity local_patches;

            for (auto g = chain_groups.cbegin() ; g != chain_groups.cend() ; ++g)
            {
                for (auto c = g->cbegin() ; c != g->cend() ; ++c)
                {
                    //                local_patches_index_lists.push_back(IndexList());

                    auto first_state = (**c).states.cbegin() + config.skip_initial * (**c).states.size();
                    auto last_state = first_state + config.sliding_window;

                    if (std::distance(last_state, (**c).states.cend()) < 0)
                    {
                        throw InternalError("PMC::hierarchical_clustering: sliding window too large for history size and skip initial: " + stringify(config.sliding_window)
                                + " vs " + stringify(std::distance(first_state, (**c).states.cend())) + " and " + stringify(config.skip_initial));
                    }

                    bool done = false;
                    while (! done)
                    {
                        // add remainder to last patch
                        if (unsigned(std::distance(last_state, (**c).states.cend())) < config.sliding_window)
                        {
                            last_state = (**c).states.cend();
                            done = true;
                        }

                        std::vector<double> mean, covariance;
                        (**c).mean_and_covariance(first_state, last_state, mean, covariance);
                        std::vector<double> center = mean;

                        if (config.patch_around_local_mode)
                        {
                            center = (**c).local_mode(first_state, last_state).point;
                        }
                        try
                        {
                            //later on during clustering use equal weights for each patch
                            HierarchicalClustering::Component patch(center, covariance, 1.0);
                            local_patches.push_back(patch);
                            hc.add(patch);
                            //                        local_patches_index_lists.back().push_back(index);
                        }
                        catch (InternalError &)
                        {
                            Log::instance()->message("PMC_sampler.hierarchical_clustering.add_patch", ll_debug)
                                    << "Skipping component, probably sliding window too small, and covariance not defined";
                        }

                        // update range for next iteration
                        first_state += config.sliding_window;
                        last_state += config.sliding_window;
                    }
                }
            }
            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Formed " << local_patches.size() << " local patches centered around"
                << (config.patch_around_local_mode ? " local modes" : " patch means");

            // number of parameters
            const unsigned & ndim = chains.front()->states.front().point.size();

            if (config.store_input_components)
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

                auto components = file.create_data_set("/hc/input-components",
                    PopulationMonteCarloSampler::Output::component_type(ndim));
                auto component_record = PopulationMonteCarloSampler::Output::component_record(ndim);

                for (auto comp = hc.begin_components() ; comp != hc.end_components() ; ++comp)
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

                for (auto comp = hc.begin_clusters() ; comp != hc.end_clusters() ; ++comp)
                {
                    std::get<0>(component_record) = comp->weight();
                    std::copy(comp->mean()->data, comp->mean()->data + ndim,
                              std::get<1>(component_record).begin());
                    std::copy(comp->covariance()->data, comp->covariance()->data + ndim * ndim,
                              std::get<2>(component_record).begin());

                    components << component_record;
                }
            }

            /* create clusters from patches  */

            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_informational)
                << "Start hierarchical clustering ";

            hc.run();

            /* initialize pmc */

            // count active clusters
            unsigned active_clusters = 0;
            for (auto cl = hc.begin_clusters() ; cl != hc.end_clusters() ; ++cl, ++active_clusters);

            // did clusters die?
            if (active_clusters != n_clusters_total)
            {
                Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_warning)
                    << "Using only " << active_clusters << " components to start PMC. "
                    << n_clusters_total - active_clusters << " died out during the hierarchical clustering.";
            }

            pmc::ErrorHandler err;

            mmv = mix_mvdens_alloc(active_clusters, ndim, err);

            unsigned i_cl = 0;
            for (auto cl = hc.begin_clusters() ; cl != hc.end_clusters() ; ++cl, ++i_cl)
            {
                // equal weight to every cluster
                mmv->wght[i_cl] = 1.0 / active_clusters;
                mvdens * mv = mmv->comp[i_cl];
                std::copy(cl->mean()->data, cl->mean()->data + ndim, mv->mean);
                std::copy(cl->covariance()->data, cl->covariance()->data + ndim * ndim, mv->std);
                mv->band_limit = ndim;
                mv->df = config.degrees_of_freedom;
                mv->chol = 0;
            }

            // check if components are in valid region
            Log::instance()->message("PMC_sampler.hierarchical_clustering", ll_debug)
                << "Filtering components that don't match";
            active_clusters = filter_components(mmv, pb);

            return active_clusters;
        }

        void pre_run()
        {
            Log::instance()->message("PMC_sampler.status", ll_informational)
                << "Starting the prerun";

            pmc::ErrorHandler err;

            // prerun to adapt proposal densities
            for (unsigned i = 0; i < config.chunks ; ++i)
            {
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
                    << "Status after step " << i + 1 << " of " << config.chunks
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
                    pmc_simu_realloc(pmc, config.chunk_size * live_components, err);
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

            if (config.final_chunk_size == 0)
                return;

            // reserve memory for final chunk
            pmc_simu_realloc(pmc, config.final_chunk_size, err);

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
                pmc_simu_realloc(pmc, config.final_chunk_size, err);
                pmc::check_error(err);
            }
            // reduce overall size if components died out
            else if (config.adjust_sample_size)
            {
                pmc_simu_realloc(pmc, live_components * config.chunk_size, err);
                pmc::check_error(err);
            }
            // todo is this ever called?
            // or if sample size has been changed externally to differ from the samples before and update
            else if (unsigned(pmc->nsamples) != prop->ncomp * config.chunk_size)
            {
            	Log::instance()->message("PMC_sampler.draw_samples", ll_debug)
                    << "I'm in a surprising place";
                pmc_simu_realloc(pmc, prop->ncomp * config.chunk_size, err);
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

    PopulationMonteCarloSampler::PopulationMonteCarloSampler(const Analysis & analysis, const PopulationMonteCarloSampler::Config & config) :
        PrivateImplementationPattern<PopulationMonteCarloSampler>(new Implementation<PopulationMonteCarloSampler>(analysis, config))
    {
    }

    PopulationMonteCarloSampler::PopulationMonteCarloSampler(const Analysis & analysis, const hdf5::File & file,
                                                             const PopulationMonteCarloSampler::Config & config,
                                                             const bool & update) :
        PrivateImplementationPattern<PopulationMonteCarloSampler>(new Implementation<PopulationMonteCarloSampler>(analysis, file, config, update))
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
         block_decomposition(false),
         component_weights(10, 1/10.0),
         component_offset(3),
         components_per_cluster(1, std::numeric_limits<unsigned>::max(), 4),
         degrees_of_freedom(-1, std::numeric_limits<int>::max(), -1),
         minimum_overlap(0, 1, 0.0),
         mode_distance(std::numeric_limits<double>::epsilon(), 1, 1e-2),
         override_global_local_proposal(false),
         random_start(true),
         single_cluster(-1),
         skip_initial(0, 1, 0.1),
         starting_points(15),
         std_dev_reduction(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max(), 1),
         group_by_r_value(1, std::numeric_limits<double>::max(), 1),
         patch_around_local_mode(false),
         r_value_no_nuisance(true),
         sliding_window(1000),
         store_input_components(false),
         store_hc_initial(false),
         super_clusters(0),
         adjust_sample_size(false),
         chunks(10),
         chunk_size(10000),
         crop_highest_weights(0),
         need_prerun(true),
         store_prerun(true),
         convergence_eff_sample_size(0, 1, 0.92),
         convergence_perplexity(0, 1, 0.92),
         ignore_eff_sample_size(false),
         minimum_eff_sample_size(0, 1, 0.1),
         minimum_perplexity(0, 1, 0.1),
         minimum_steps(2, std::numeric_limits<unsigned>::max(), 3),
         maximum_relative_std_deviation(0, 1, 0.01),
         final_chunk_size(20000),
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

         config.chunk_size = 3000;

         return config;
     }

     std::ostream & operator<<(std::ostream & stream, const PopulationMonteCarloSampler::Config & c)
     {
         stream << std::boolalpha
                << "Clustering options: " << std::endl
                << "critical R value = " << c.group_by_r_value
                << ", ignore groups = " << stringify_container(c.ignore_groups)
                << ", R value no nuisance = " << c.r_value_no_nuisance << std::endl
                << "sliding window = " << c.sliding_window
                << ", number of clusters = " << c.super_clusters << std::endl
                << "Prerun options:" << std::endl
                << "chunk size = " << c.chunk_size
                << ", max #updates = " << c.chunks
                << ", adjust sample size = " << c.adjust_sample_size << std::endl
                << "degrees of freedom = " << c.degrees_of_freedom << std::endl
                << "Convergence options:" << std::endl
                << "ignore ESS = " << c.ignore_eff_sample_size
                << ", allowed std. dev = " << c.maximum_relative_std_deviation << std::endl
                << "Main run options: " << std::endl
                << "chunk size = " << c.final_chunk_size;
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
