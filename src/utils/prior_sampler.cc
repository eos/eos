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

#include <src/utils/prior_sampler.hh>
#include <src/utils/log.hh>
#include <src/utils/log_prior.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/thread_pool.hh>

#include <gsl/gsl_randist.h>

namespace eos
{
    struct Worker
    {
        ObservableSet observables;

        // The priors for all parameters to be varied.
        std::vector<LogPriorPtr> priors;

        // Parameter, minimum, maximum, nuisance, discrete
        std::vector<ParameterDescription> parameter_descriptions;

        // Random number generator
        gsl_rng * rng;

        // list of records, with observables and, if requested, the parameter values
        std::vector<std::vector<double>> data;

        bool store_parameters;

        Worker(const ObservableSet & observables,
               const std::vector<LogPriorPtr> & priors,
               const std::vector<ParameterDescription> & parameter_descriptions,
               unsigned seed, bool store_parameters) :
                   store_parameters(store_parameters)
        {
            Parameters p = Parameters::Defaults();

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

            // setup random number generator
            rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, seed);
        }

        ~Worker()
        {
            // free RN generator
            gsl_rng_free(rng);
        }

        // store the values of observables and optionally parameters to hdf5 scan file
        void dump_history(ScanFile::DataSet & data_set)
        {
            ScanFile::WriteBuffer buffer(data_set.fields(), data.size());

            for (auto i = data.begin(), i_end = data.end(); i != i_end; ++i)
            {
                buffer << *i;

                // once the buffer is full, write to disk
                if (buffer.capacity() == buffer.size())
                {
                    data_set << buffer;
                    buffer.clear();
                }
            }

            // write the remainder to disk
            data_set << buffer;
            buffer.clear();
        }

        void run(unsigned iterations)
        {
            for (unsigned i = 0 ; i < iterations ; ++i)
            {
                std::vector<double> record;

                // draw a sample
                unsigned index = 0;
                for (auto i = priors.begin(), i_end = priors.end() ; i != i_end ; ++i, ++index)
                {
                    double p = (**i).sample(rng);
                    parameter_descriptions[index].parameter = p;
                    if (store_parameters)
                        record.push_back(p);
                }

                // evaluate observables
                for (auto i = observables.begin(), i_end = observables.end(); i != i_end; ++i)
                {
                    record.push_back((**i).evaluate());
                }
                data.push_back(record);
            }
        }
    };

    template<>
    struct Implementation<PriorSampler>
    {
        // Holds the configuration options.
        PriorSampler::Config config;

        // Keep all distinct observables.
        ObservableSet observables;

        // Names of all parameters. Prevents using a parameter twice
        std::set<std::string> parameter_names;

        // The priors for all parameters to be varied.
        std::vector<LogPriorPtr> priors;

        // Parameter, minimum, maximum, nuisance, discrete
        std::vector<ParameterDescription> parameter_descriptions;

        // tickets for parallel computations
        std::vector<Ticket> tickets;

        Implementation(const ObservableSet & observables, const PriorSampler::Config & config) :
            config(config),
            observables(observables)
        {
            if (!config.output_file)
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
                if (! result.second)
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
            Log::instance()->message("prior_sampler.run", ll_informational)
                << "Starting to generate " << config.chunks * config.chunk_size << " samples.";

            // setup the scan file
            ScanFile::DataSet data_set = setup_output();

            // start with empty ticket queue
            tickets.clear();

            // create one Worker per chunk
            std::vector<std::shared_ptr<Worker>> workers;

            for (unsigned chunk = 0 ; chunk < config.chunks ; ++chunk)
            {
                auto w = std::make_shared<Worker>(observables, priors, parameter_descriptions, config.seed + chunk, config.store_parameters);
                workers.push_back(w);
                if (config.parallelize)
                {
                    // make sure to pass the pointer, instead of a reference from *w, to bind.
                    // Else copies are created, which lead to a double freeing upon calling the destructor the 2nd time.
                    tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&Worker::run, w.get(), config.chunk_size)));
                }
                else
                {
                    w->run(config.chunk_size);
                }
            }

            // wait for job completion
            for (auto t = tickets.begin(), t_end = tickets.end() ; t != t_end ; ++t)
            {
                t->wait();
            }

            // retrieve data
            for (auto i = workers.begin(), i_end = workers.end() ; i != i_end ; ++i)
            {
                (**i).dump_history(data_set);
            }

            // all tickets finished
            tickets.clear();

            Log::instance()->message("prior_sampler.run", ll_informational)
                << "Sampling completed.";
        }

        ScanFile::DataSet
        setup_output()
        {
            // TODO adjust for multidimensional priors?
            ScanFile::DataSet data_set = ScanFile::DataSet(config.output_file->add("data",
                observables.size() + (config.store_parameters ? priors.size() : 0)));

            // set the attributes
            auto f = data_set.begin_fields();

            if (config.store_parameters)
            {
                for (auto i = priors.begin(), i_end = priors.end(); i != i_end; ++i, ++f)
                {
                    // iterator to pointer to iterator. Puh!
                    auto description = *((**i).begin());

                    f->name(description.parameter.name());

                    f->set("min", description.min);
                    f->set("max", description.max);
                    f->set("nuisance", description.nuisance);
                    f->set("discrete", description.discrete);
                }
            }

            // store observable names
            for (auto i = observables.begin(), i_end = observables.end(); i != i_end; ++i, ++f)
            {
                f->name((**i).name());
            }

            if (data_set.end_fields() != f)
            {
                throw InternalError("PriorSampler::setup_output: insufficient number of fields "
                        + stringify(std::distance(data_set.begin_fields(), data_set.end_fields())) + " in data set, for "
                        + stringify(priors.size()) + " parameters and " + stringify(observables.size()) + " observables.");
            }

            return data_set;
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

    void
    PriorSampler::run()
    {
        _imp->run();
    }

    PriorSampler::Config::Config() :
        chunks(4),
        chunk_size(10000),
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
