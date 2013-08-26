/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012 Frederik Beaujean
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

#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/log_prior.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/prior_sampler.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/thread_pool.hh>

#include <gsl/gsl_randist.h>

namespace eos
{

    template<>
    struct Implementation<PriorSampler>
    {
        struct Worker
    {
            ObservableSet observables;

            // The priors for all parameters to be varied.
            std::vector<LogPriorPtr> priors;

            // Parameter, minimum, maximum, nuisance, discrete
            std::vector<ParameterDescription> parameter_descriptions;

            // Random number generator seed
            unsigned seed;

            // clear cache every N iterations
            const unsigned clear_cache;

            // the sampling output, one vector for each iteration
            std::vector<std::vector<double>> observable_values;
            std::vector<std::vector<double>> parameter_values;

            hdf5::Array<1, double> observable_type;
            hdf5::Array<1, double> parameter_type;

            bool store_parameters;

            Worker(const ObservableSet & observables,
                   const std::vector<LogPriorPtr> & priors,
                   const std::vector<ParameterDescription> & parameter_descriptions,
                   unsigned seed, bool store_parameters) :
                       seed(seed),
                       clear_cache(500),
                       observable_type("observables", { observables.size() }),
                       parameter_type("parameters", { parameter_descriptions.size() }),
                       store_parameters(store_parameters)
            {
                // need to clone, so parameters that are fixed by hand have correct value
                // cast away const-ness
                //            ObservableSet * ptr = const_cast<ObservableSet *>(&observables);
                Parameters p = const_cast<ObservableSet *>(&observables)->parameters().clone();

                // clone observables
                for (auto i = observables.begin(), i_end = observables.end() ; i != i_end ; ++i)
                {
                    this->observables.add((*i)->clone(p));
                }

                // clone priors
                for (auto i = priors.begin(), i_end = priors.end() ; i != i_end ; ++i)
                {
                    this->priors.push_back((*i)->clone(p));
                }

                // clone descriptions
                for (auto i = parameter_descriptions.begin(), i_end = parameter_descriptions.end() ; i != i_end ; ++i)
                {
                    this->parameter_descriptions.push_back(
                        ParameterDescription{ p[i->parameter.name()], i->min, i->max, i->nuisance, i->discrete });
                }
            }

            void dump_history(const std::shared_ptr<hdf5::File> & file)
            {
                // write observables
                auto observable_data_set = file->create_or_open_data_set("/data/observables", observable_type);
                for (auto o = observable_values.cbegin(), o_end = observable_values.cend() ; o != o_end ; ++o)
                {
                    observable_data_set << *o;
                }

                // write parameters
                if ( ! store_parameters)
                    return;

                auto parameter_data_set = file->create_or_open_data_set("/data/parameters", parameter_type);
                for (auto p = parameter_values.cbegin(), p_end = parameter_values.cend() ; p != p_end ; ++p)
                {
                    parameter_data_set << *p;
                }
            }

            void run(unsigned iterations)
            {
                // setup random number generator
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, seed);

                for (unsigned i = 0 ; i < iterations ; ++i)
                {
                    std::vector<double> parameter_sample;
                    std::vector<double> observable_sample;

                    // draw a sample
                    unsigned index = 0;
                    for (auto prior = priors.begin(), i_end = priors.end() ; prior != i_end ; ++prior, ++index)
                    {
                        double p = (**prior).sample(rng);
                        parameter_descriptions[index].parameter = p;
                        if (store_parameters)
                        {
                            parameter_sample.push_back(p);
                        }
                    }
                    if (store_parameters)
                    {
                        parameter_values.push_back(parameter_sample);
                    }

                    // evaluate observables
                    for (auto obs = observables.begin(), i_end = observables.end(); obs != i_end; ++obs)
                    {
                        observable_sample.push_back((**obs).evaluate());
                    }
                    observable_values.push_back(observable_sample);

                    // clear cache, but never with less than a 100 iterations
                    if (i > 0 && i % clear_cache == 0)
                    {
                        Log::instance()->message("prior_sampler.run_clean_up", ll_informational)
                        << "Clearing the memoise cache.";
                        MemoisationControl::instance()->clear();
                    }
                }
                // free RN generator
                gsl_rng_free(rng);
            }

            /*
             * Calculate observables for the given parameter samples
             */
            void run_samples(const std::vector<std::vector<double>> & samples)
            {
                std::vector<double> observable_sample;

                const unsigned & n_dim = samples.front().size();

                // loop over samples
                unsigned sample_index = 0;
                for (auto s = samples.begin() ; s != samples.end() ; ++s, ++sample_index)
                {
#if 1
                    if (s->size() != n_dim)
                        throw InternalError("prior_sampler.run: Found mismatch between parameter dimensions for sample "
                            + stringify(sample_index) + ": " + stringify(n_dim) + " vs " + stringify(s->size()));
#endif
                    // read and update parameter values
                    unsigned par_index = 0;
                    for (auto p = s->begin() ; p != s->end() ; ++p, ++par_index)
                    {
                        parameter_descriptions[par_index].parameter = *p;
                    }

                    // evaluate observables
                    observable_sample.clear();
                    for (auto obs = observables.begin(), i_end = observables.end(); obs != i_end; ++obs)
                    {
                        observable_sample.push_back((**obs).evaluate());
                    }
                    observable_values.push_back(observable_sample);

                    // clear cache every fixed number of iterations, but not at beginning
                    if (sample_index % std::max(1u, clear_cache) == 0)
                    {
                        Log::instance()->message("prior_sampler.run_clean_up", ll_informational)
                                                << "Clearing the memoise cache.";
                        MemoisationControl::instance()->clear();
                    }
                }
            }
        };

        // Holds the configuration options.
        PriorSampler::Config config;

        // Keep all distinct observables.
        ObservableSet observables;

        // s_min, s_max
        hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>> observable_description_type;

        // Names of all parameters. Prevents using a parameter twice
        std::set<std::string> parameter_names;

        // The priors for all parameters to be varied.
        std::vector<LogPriorPtr> priors;

        // Parameter, minimum, maximum, nuisance, discrete
        std::vector<ParameterDescription> parameter_descriptions;

        // min, max, nuisance
        hdf5::Composite<hdf5::Scalar<double>, hdf5::Scalar<double>> parameter_descriptions_type;

        // tickets for parallel computations
        std::vector<Ticket> tickets;

        Implementation(const ObservableSet & observables, const PriorSampler::Config & config) :
            config(config),
            observables(observables),
            observable_description_type
            {
                "kinematics",
                hdf5::Scalar<double>("s_min"),
                hdf5::Scalar<double>("s_max"),
            },
            parameter_descriptions_type
            {
                "parameter description",
                hdf5::Scalar<double>("min"),
                hdf5::Scalar<double>("max"),
            }
        {
            if ( ! config.output_file)
                throw InternalError("PriorSampler(): Missing valid output file");
        }

        bool add(const LogPriorPtr & prior)
        {
            // clone has correct Parameters object selected
            LogPriorPtr prior_clone = prior->clone(observables.parameters());

            // check if param exists already
            // read out parameters from prior
            unsigned counter = 0;
            for (auto d = prior->begin(), d_end = prior->end() ; d != d_end ; ++d, ++counter)
            {
                if (counter > 0)
                    throw InternalError("PriorSampler::add(): eos is not ready for multidimensional priors.");
                auto result = parameter_names.insert(d->parameter.name());
                if ( ! result.second)
                    return false;

                parameter_descriptions.push_back(*d);
            }

            // then add to prior container
            priors.push_back(prior_clone);

            return true;
        }

        bool add(const ObservablePtr & observable)
        {
            return observables.add(observable).second;
        }

        void run()
        {
            // setup the scan file
            setup_output();

            Log::instance()->message("prior_sampler.run", ll_informational)
                << "Starting to generate " << config.n_samples << " samples.";

            // start with empty ticket queue
            tickets.clear();

             // create one Worker per chunk
            std::vector<std::shared_ptr<Worker>> workers;
            const unsigned average_samples_per_worker = config.n_samples / config.n_workers;
            const unsigned remainder = config.n_samples % config.n_workers;

            for (unsigned chunk = 0 ; chunk < config.n_workers; ++chunk)
            {
                workers.push_back(std::make_shared<Worker>(observables, priors, parameter_descriptions,
                                                           config.seed + chunk, config.store_parameters));

                unsigned samples_per_worker = average_samples_per_worker;

                // last worker gets the remainder
                if (chunk == config.n_workers - 1)
                    samples_per_worker += remainder;

                if (config.parallelize)
                {
                    // make sure to pass the pointer, instead of a reference from *w, to bind.
                    // Else copies are created, which lead to a double freeing upon calling the destructor the 2nd time.
                    tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&Worker::run, workers.back().get(), samples_per_worker)));
                }
                else
                {
                    workers.back()->run(samples_per_worker);
                }
            }

            // wait for job completion
            for (auto t = tickets.begin(), t_end = tickets.end() ; t != t_end ; ++t)
            {
                t->wait();
            }

            // retrieve data and delete workers
            for (auto w = workers.begin(), w_end = workers.end() ; w != w_end ; ++w)
            {
                (**w).dump_history(config.output_file);
            }

            // all tickets finished
            tickets.clear();

            Log::instance()->message("prior_sampler.run", ll_informational)
                << "Sampling completed.";
        }

        void
        run(const std::vector<std::vector<double>> & samples, const std::vector<ParameterDescription> & defs)
        {
            // no need to store parameters
            config.store_parameters = false;

            // setup the scan file
            setup_output();

            // start with empty ticket queue
            tickets.clear();

            Worker worker(observables, priors, defs, config.seed, config.store_parameters);
            worker.run_samples(samples);
            worker.dump_history(config.output_file);

            Log::instance()->message("prior_sampler.run", ll_informational)
                << "Calculations completed.";
        }

        void
        setup_output()
        {
            // write parameter descriptions
            unsigned counter = 0;
            for (auto d = parameter_descriptions.cbegin(), d_end = parameter_descriptions.cend() ; d != d_end ; ++d, ++counter)
            {
                auto components = config.output_file->create_data_set("/descriptions/parameters/" + stringify(counter),
                    parameter_descriptions_type);
                auto record = std::make_tuple(d->min, d->max);
                components << record;

                auto attr_name = components.create_attribute("name", hdf5::Scalar<const char *>("name"));
                attr_name = d->parameter.name().c_str();

                auto attr_prior = components.create_attribute("prior", hdf5::Scalar<const char *>("prior"));
                attr_prior = priors[std::distance(parameter_descriptions.cbegin(), d)]->as_string().c_str();
            }

            Log::instance()->message("prior_sampler.run", ll_informational)
                << "Computing the SM prediction for each observable with fixed parameter values "
                << "before drawing random parameter samples.";

            // write observable descriptions
            counter = 0;
            for (auto o = observables.begin(), o_end = observables.end() ; o != o_end ; ++o, ++counter)
            {
                auto components = config.output_file->create_data_set("/descriptions/observables/" + stringify(counter),
                    observable_description_type);
                auto record = std::make_tuple(0.0, 0.0);
                try
                {
                    record = std::make_tuple((**o).kinematics()["s_min"],(**o).kinematics()["s_max"]);
                }
                catch (...)
                {
                }
                components << record;

                auto attr_name = components.create_attribute("name", hdf5::Scalar<const char *>("name"));
                attr_name = (**o).name().c_str();

                auto attr_options = components.create_attribute("options", hdf5::Scalar<const char *>("options"));
                attr_options = (**o).options().as_string().c_str();

                auto attr_kinematics = components.create_attribute("kinematics", hdf5::Scalar<const char *>("kinematics"));
                attr_kinematics = (**o).kinematics().as_string().c_str();

                auto attr_sm_prediction = components.create_attribute("SM prediction", hdf5::Scalar<double>("SM prediction"));
                double value = (**o).evaluate();
                attr_sm_prediction = value;
            }
        }
    };

    PriorSampler::PriorSampler(const ObservableSet & observables,  const Config & config) :
        PrivateImplementationPattern<PriorSampler>(new Implementation<PriorSampler>(observables, config))
    {
    }

    PriorSampler::~PriorSampler()
    {
    }

    bool
    PriorSampler::add(const LogPriorPtr & prior)
    {
        return _imp->add(prior);
    }

    bool
    PriorSampler::add(const ObservablePtr & observable)
    {
        return _imp->add(observable);
    }

    PriorSampler::ObservablesType
    PriorSampler::observables_type(const unsigned & dimension)
    {
        return ObservablesType("observables", { dimension });
    }

    void
    PriorSampler::run()
    {
        _imp->run();
    }

    void
    PriorSampler::run(const std::vector<std::vector<double>> & samples, const std::vector<ParameterDescription> & defs)
    {
        _imp->run(samples, defs);
    }

    PriorSampler::Config::Config() :
        n_samples(100000),
        n_workers(4),
        parallelize(true),
        seed(1234623),
        store_parameters(false)
    {
    }

    PriorSampler::Config
    PriorSampler::Config::Default()
    {
        return PriorSampler::Config();
    }
}
