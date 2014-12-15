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
namespace
{
typedef PriorSampler::SamplesList SamplesList;
}

    template<>
    struct Implementation<PriorSampler>
    {
        typedef std::function<void (void)> Function;

        struct Worker
        {
            ObservableSet observables;

            // The priors for all parameters to be varied.
            std::vector<LogPriorPtr> priors;

            // Parameter, minimum, maximum, nuisance, discrete
            std::vector<ParameterDescription> parameter_descriptions;

            // Random number generator seed
            unsigned seed;

            // the sampling output, one vector for each iteration
            SamplesList observable_samples;
            SamplesList parameter_samples;

            hdf5::Array<1, double> observable_type;
            hdf5::Array<1, double> parameter_type;

            Worker(const ObservableSet & observables,
                   const std::vector<LogPriorPtr> & priors,
                   const std::vector<ParameterDescription> & parameter_descriptions,
                   unsigned seed) :
                       seed(seed),
                       observable_type("observables", { observables.size() }),
                       parameter_type("parameters", { parameter_descriptions.size() })
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

            void dump_history(const std::shared_ptr<hdf5::File> & file, const bool & store_parameters)
            {
                // write observables
                auto observable_data_set = file->create_or_open_data_set("/data/observables", observable_type);
                for (auto o = observable_samples.cbegin(), o_end = observable_samples.cend() ; o != o_end ; ++o)
                {
                    observable_data_set << *o;
                }

                // write parameters
                if ( ! store_parameters)
                    return;

                auto parameter_data_set = file->create_or_open_data_set("/data/parameters", parameter_type);
                for (auto p = parameter_samples.cbegin(), p_end = parameter_samples.cend() ; p != p_end ; ++p)
                {
                    parameter_data_set << *p;
                }
            }

            /*!
             * Compute observables for every sample in range
             */
            void compute_observables(      SamplesList::const_iterator & first,
                    const SamplesList::const_iterator & last)
            {
                Log::instance()->message("prior_sampler.run", ll_informational)
                            << "Computing " << observables.size() << " observables for "
                            << std::distance(first, last) << " parameter samples";

                // setup random number generator
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, seed);

                // loop over samples
                for (; first != last; ++first)
                {
                    // read and update parameter values, one at a time
                    auto def = parameter_descriptions.begin();
                    for (auto p = first->cbegin() ; p != first->cend() ; ++p, ++def)
                    {
                        def->parameter = *p;
                    }

                    auto p = priors.cbegin();
                    std::advance(p, std::distance(first->cbegin(), first->cend()));
                    for (auto p_end = priors.cend() ; p != p_end ; ++p, ++def)
                    {
                        def->parameter = (*p)->sample(rng);
                    }

                    // calculate all observables
                    std::vector<double> observable_sample;
                    for (auto & o : observables)
                        observable_sample.push_back(o->evaluate());
                    observable_samples.push_back(observable_sample);
                }

                // free RN generator
                gsl_rng_free(rng);
            }

            /*!
             * Draw random vector from the priors.
             *
             * @param iterations The number of samples.
             */
            void draw_samples(unsigned iterations)
            {
                Log::instance()->message("prior_sampler.run", ll_informational)
                            << "Drawing " << iterations << " parameter samples";

                // setup random number generator
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, seed);

                for (unsigned i = 0 ; i < iterations ; ++i)
                {
                    std::vector<double> parameter_sample;

                    // draw a sample
                    unsigned index = 0;
                    for (auto prior = priors.begin(), i_end = priors.end() ; prior != i_end ; ++prior, ++index)
                    {
                        double p = (**prior).sample(rng);
                        parameter_descriptions[index].parameter = p;
                        parameter_sample.push_back(p);
                    }
                    parameter_samples.push_back(parameter_sample);
                }

                // free RN generator
                gsl_rng_free(rng);
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

        void run(const SamplesList & samples, const std::vector<ParameterDescription> & defs)
        {
            this->parameter_descriptions.insert(this->parameter_descriptions.begin(), defs.begin(), defs.end());

            {
                std::vector<LogPriorPtr> priors;
                for (auto d = defs.cbegin(), d_end = defs.cend() ; d != d_end ; ++d)
                {
                    priors.push_back(LogPrior::Flat(Parameters::Defaults(), d->parameter.name(), ParameterRange{ d->min, d->max }));
                }

                this->priors.insert(this->priors.begin(), priors.begin(), priors.end());
            }

            // setup the scan file
            setup_output();

            // start with empty ticket queue
            tickets.clear();

            const bool draw = samples.empty();

            if (! draw)
            {
                config.n_samples = samples.size();
                config.store_parameters = false;
            }

            // create one Worker per chunk
            std::vector<std::shared_ptr<Worker>> workers;
            const unsigned average_samples_per_worker = config.n_samples / config.n_workers;
            const unsigned remainder = config.n_samples % config.n_workers;

            for (unsigned chunk = 0 ; chunk < config.n_workers; ++chunk)
            {
                workers.push_back(std::make_shared<Worker>(observables, this->priors, this->parameter_descriptions,
                        config.seed + chunk));

                unsigned samples_per_worker = average_samples_per_worker;

                // last worker gets the remainder
                if (chunk == config.n_workers - 1)
                    samples_per_worker += remainder;

                auto first = samples.cbegin() + chunk * average_samples_per_worker;
                auto last  = first + samples_per_worker;
                if (draw)
                {
                    workers.back()->draw_samples(samples_per_worker);
                    first = workers.back()->parameter_samples.cbegin();
                    last  = workers.back()->parameter_samples.cend();
                }

                Function f = std::bind(&Worker::compute_observables, workers.back().get(), first, last);

                if (config.parallelize)
                {
                    // make sure to pass the pointer, instead of a reference from *w, to bind.
                    // Else copies are created, which lead to a double freeing upon calling the destructor the 2nd time.
                    tickets.push_back(ThreadPool::instance()->enqueue(f));
                }
                else
                {
                    f();
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
                (**w).dump_history(config.output_file, config.store_parameters);
            }

            // all tickets finished
            tickets.clear();

            Log::instance()->message("prior_sampler.run", ll_informational)
                        << "Observable computations completed.";
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
                << "Computing the SM prediction for each observable with fixed parameter values ";

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
        _imp->run(SamplesList(), std::vector<ParameterDescription>());
    }

    void
    PriorSampler::run(const SamplesList & samples, const std::vector<ParameterDescription> & defs)
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
