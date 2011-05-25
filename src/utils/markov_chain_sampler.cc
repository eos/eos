/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2011 Danny van Dyk
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

#include <src/utils/log.hh>
#include <src/utils/markov_chain_sampler.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/power_of.hh>
#include <src/utils/rvalue.hh>
#include <src/utils/stringify.hh>
#include <src/utils/thread_pool.hh>

#include <limits>
#include <fstream>
#include <sys/stat.h>

namespace eos
{
    template<>
    struct Implementation<MarkovChainSampler>
    {
        // store reference, but don't own analysis
        const Analysis & analysis;

        // our configuration options
        MarkovChainSampler::Config config;

        // pointer to the HDF5 output file.
        std::shared_ptr<ScanFile> output_file;

        // tickets for parallel computations
        std::vector<Ticket> tickets;

        // number of scan parameters
        unsigned number_of_parameters;

        // independent chains
        std::vector<MarkovChain> chains;
        std::vector<ScanFile::DataSet> data_sets;

        // collected information regarding the prerun
        MarkovChainSampler::PreRunInfo pre_run_info;

        std::function<double (const std::vector<double> &, const std::vector<double> &, const unsigned &)> compute_rvalue;

        Implementation(const Analysis & analysis, const MarkovChainSampler::Config & config) :
            analysis(analysis),
            config(config),
            compute_rvalue(config.use_strict_rvalue_definition ? &RValue::gelman_rubin : &RValue::approximation)
        {
            initialize();
        }

        /*
         * Checks efficiencies, adjusts if needed
         * return true if all efficiencies in ranges defined by MarkovChainConfig::min_efficiency, MarkovChainConfig::max_efficiency
         */
        bool adjust_scales()
        {
            bool all_efficiencies_ok = true;

            // store parameter efficiencies for printing at the end
            std::vector<double> efficiencies(number_of_parameters, 0.0);

            // loop over chains
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                // loop over parameters
                for (unsigned p = 0; p < number_of_parameters; ++p)
                {
                    // read out efficiencies
                    double accepted(c->statistics().iterations_accepted[p]);
                    double rejected(c->statistics().iterations_rejected[p]);
                    efficiencies[p] = accepted / (accepted + rejected);

                    Log::instance()->message("markov_chain_sampler.efficiencies", ll_debug)
                        << "Current efficiency for \"" << c->parameter_descriptions()[p].parameter.name()
                        << "\" in chain " << std::distance(chains.begin(), c) << ": " << efficiencies[p];

                    //do not adjust scales if it has no effect on proposal function
                    if (c->parameter_descriptions()[p].discrete)
                        continue;

                    double new_scale = c->get_scale(p);

                    // up/down of scale

                    // decrease scale to raise efficiency
                    if (efficiencies[p] < config.min_efficiency)
                    {
                        new_scale /= 2.0;
                    }
                    if (efficiencies[p] > config.max_efficiency)
                    {
                        new_scale *= 2.0;
                    }

                    // impose min/max scale values
                    // spit out warning if scale exceed allowed ranges
                    if ((new_scale < config.scale_min) || (new_scale > config.scale_max))
                    {
                        Log::instance()->message("markov_chain_sampler.rescale_beyond_limits", ll_warning)
                            << "Attempting to update scale factor of parameter '" << c->parameter_descriptions()[p].parameter.name()
                            << "' in chain " << std::distance(chains.begin(), c) << " beyond its limits [" << config.scale_min << ","
                            << config.scale_max << "]. Scale kept at current value. May be the chosen parameter range is wrong?";

                        new_scale = c->get_scale(p);
                    }

                    // did value change?
                    if (new_scale != c->get_scale(p))
                    {
                        Log::instance()->message("markov_chain_sampler.update_scale", ll_informational)
                            << "Updating scale of parameter '" << c->parameter_descriptions()[p].parameter.name()
                            << "' in chain " << std::distance(chains.begin(), c) << " to " << new_scale;

                        // update scale
                        c->set_scale(p, new_scale);

                        all_efficiencies_ok = false;
                    }
                }
            }

            return all_efficiencies_ok;
        }

        bool check_convergence()
        {
            bool efficiencies_ok = adjust_scales();

            bool rvalues_ok = true;

            // no R-value for single chain
            if (config.number_of_chains > 1)
            {
                rvalues_ok = check_rvalues();
            }

            if (efficiencies_ok && rvalues_ok)
                return true;
            else
                return false;
        }

        bool check_rvalues()
        {
            bool all_rvalues_small = true;

            // loop over all parameters to check and get R-values
            for (unsigned i = 0 ; i < number_of_parameters ; ++i)
            {
                std::vector<double> chain_means, chain_variances;

                // read out statistics
                for (auto c = chains.begin(), c_end = chains.end() ;  c != c_end ; ++c)
                {
                    chain_means.push_back(c->statistics().mean_of_parameters[i]);
                    chain_variances.push_back(c->statistics().variance_of_parameters[i]);
                }

                pre_run_info.rvalue_parameters[i] = compute_rvalue(chain_means, chain_variances, pre_run_info.iterations);

                if (pre_run_info.rvalue_parameters[i] > config.rvalue_criterion_param || isnan(pre_run_info.rvalue_parameters[i]))
                {
                    all_rvalues_small = false;

                    Log::instance()->message("markov_chain_sampler.parameter_rvalue_too_large", ll_informational)
                        << "R-value of parameter '" << chains.front().parameter_descriptions()[i].parameter.name() << "' is too large: "
                        << pre_run_info.rvalue_parameters[i] << " > " << config.rvalue_criterion_param;
                }
            }

            // posterior converged?
            if (config.use_posterior_rvalue)
            {
                // read out statistics
                std::vector<double> chain_means, chain_variances;
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end; ++c)
                {
                    chain_means.push_back(c->statistics().mean_of_posterior);
                    chain_variances.push_back(c->statistics().variance_of_posterior);
                }

                pre_run_info.rvalue_posterior = compute_rvalue(chain_means, chain_variances, pre_run_info.iterations);
                if ((pre_run_info.rvalue_posterior > config.rvalue_criterion_posterior) || isnan(pre_run_info.rvalue_posterior))
                {
                    all_rvalues_small = false;

                    Log::instance()->message("markov_chain_sampler.posterior_rvalue_too_large", ll_informational)
                        << "R-value of posterior is too large: "
                        << pre_run_info.rvalue_posterior << " > " << config.rvalue_criterion_posterior;
                }
            }

            if (all_rvalues_small)
            {
                Log::instance()->message("markov_chain_sampler.convergence", ll_informational)
                    << "Convergence achieved";
            }

            return all_rvalues_small;
        }

        void dump_hdf5()
        {
            Log::instance()->message("markov_chain_sampler.dump_hdf5", ll_debug)
                << "Dumping all chains to HDF5";

            auto c = chains.begin();
            for (auto d = data_sets.begin(), d_end = data_sets.end() ; d != d_end ; ++c, ++d)
            {
                c->dump_history(*d);
                c->clear();
            }
        }

        // common method to call from multiple constructors
        void initialize()
        {
            // create independent chains -> different seeds
            for (unsigned i = 0 ; i < config.number_of_chains ; ++i)
            {
                AnalysisPtr analysis_clone(analysis.clone());

                unsigned long seed = config.seed + i;

                MarkovChain chain(analysis_clone, seed);
                chain.set_scale(config.scale_initial);
                chains.push_back(chain);
            }

            number_of_parameters = analysis.parameter_descriptions().size();

            // setup prerun info
            pre_run_info = MarkovChainSampler::PreRunInfo
            {
                false, 0,0, std::numeric_limits<double>::max(), std::vector<double>(analysis.parameter_descriptions().size(),std::numeric_limits<double>::max())
            };

            if (config.output_file)
            {
                for (unsigned i = 0 ; i < config.number_of_chains ; ++i)
                {
                    // TODO:
                    //  * more than posterior?
                    ScanFile::DataSet data_set = ScanFile::DataSet(config.output_file->add("chain #" + stringify(i), number_of_parameters + 1));

                    auto f = data_set.begin_fields();
                    auto pp = analysis.parameter_descriptions();
                    for (auto p = pp.begin(), p_end = pp.end() ; p != p_end ; ++p, ++f)
                    {
                        f->name = p->parameter.name();
                        f->min = p->min;
                        f->max = p->max;
                        f->nuisance = p->nuisance;
                        f->discrete = p->discrete;
                    }

                    if (data_set.end_fields() == f)
                    {
                        throw InternalError("MarkovChainSampler::initialize: insufficient number of fields "
                                + stringify(std::distance(data_set.begin_fields(), data_set.end_fields())) + " in data set, for "
                                + stringify(number_of_parameters) + " parameters");
                    }

                    f->name = std::string("posterior");
                    data_sets.push_back(data_set);
                }
            }
        }

        /*
         * Collect samples from posterior and check for convergence.
         */
        void pre_run()
        {
            Log::instance()->message("markov_chain_sampler.prerun_start", ll_informational)
                << "Commencing the pre-run";

            // require:chains are initialized
            pre_run_info.converged = false;
            pre_run_info.iterations = 0;

            // set up chains
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                // save history
                c->keep_history(true);
            }

            // keep going till maxIter or  break when convergence estimated

            while (pre_run_info.iterations < config.prerun_iterations_min || (!pre_run_info.converged && pre_run_info.iterations
                            < config.prerun_iterations_max))
            {
                // start with empty ticket queue
                tickets.clear();

                // loop over chains
                // run each chain for N iterations
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                {
                    if (config.parallelize)
                    {
                        tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&MarkovChain::run, *c, config.prerun_iterations_update)));
                    }
                    else
                    {
                        c->run(config.prerun_iterations_update);
                    }
                }

                // wait for job completion
                for (auto t = tickets.begin(), t_end = tickets.end(); t != t_end; ++t)
                {
                    t->wait();
                }

                // all tickets finished
                tickets.clear();

                pre_run_info.iterations += config.prerun_iterations_update;
                Log::instance()->message("markov_chain_sampler.prerun_progress", ll_informational)
                    << "Pre-run has completed " << pre_run_info.iterations << " iterations";

                pre_run_info.converged = check_convergence();
            }

            if (pre_run_info.converged)
            {
                Log::instance()->message("markov_chain_sampler.prerun_converged", ll_informational)
                    << "Pre-run has converged after " << pre_run_info.iterations << " iterations";

                if (config.number_of_chains < 2)
                {
                    Log::instance()->message("markov_chain_sampler.single_chain", ll_warning)
                        << "R-values are undefined for a single chain, so only efficiencies were adjusted";
                }

                pre_run_info.iterations_at_convergence = pre_run_info.iterations;
            }
            else
            {
                Log::instance()->message("markov_chain_sampler.no_convergence", ll_error)
                    << "Pre-run did NOT converge!";
            }

            // clear the history after the prerun
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                c->clear();
            }
        }

        /*
         * Collect samples from posterior for further analysis.
         * No convergence checking here.
         */
        void main_run()
        {
            Log::instance()->message("markov_chain_sampler.mainrun_start", ll_informational)
                << "Commencing the main-run";

            // set up chains
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end; ++c)
            {
                // save history?
                c->keep_history(config.store);
            }

            for (unsigned chunk = 0 ; chunk < config.chunks ; ++chunk)
            {

                // start with empty ticket queue
                tickets.clear();

                // loop over chains
                // run each chain for N iterations
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                {
                    if (config.parallelize)
                    {
                        tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&MarkovChain::run, *c, config.chunk_size)));
                    }
                    else
                    {
                        c->run(config.chunk_size);
                    }
                }

                // wait for job completion
                for (auto t = tickets.begin(), t_end = tickets.end(); t != t_end; ++t)
                {
                    t->wait();
                }

                // all tickets finished
                tickets.clear();

                Log::instance()->message("markov_chain_sampler.mainrun_progress", ll_informational)
                    << "Main-run has completed " << (chunk + 1) * config.chunk_size << " iterations";

                if (config.store)
                {
                    dump_hdf5();
                }
            }
            Log::instance()->message("markov_chain_sampler.mainrun_end", ll_informational)
                << "Finished the main-run";
        }

        void run()
        {
            if (config.need_prerun)
            {
                pre_run();

                //  clear up
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                {
                    c->clear();
                }
            }

            main_run();
        }
    };

    MarkovChainSampler::MarkovChainSampler(const Analysis & analysis, const MarkovChainSampler::Config & config) :
        PrivateImplementationPattern<MarkovChainSampler>(new Implementation<MarkovChainSampler>(analysis, config))
    {
    }

    MarkovChainSampler::~MarkovChainSampler()
    {
    }

    MarkovChainSampler::PreRunInfo
    MarkovChainSampler::pre_run_info()
    {
        return _imp->pre_run_info;
    }

    void
    MarkovChainSampler::run()
    {
        _imp->run();
    }

    const MarkovChainSampler::Config &
    MarkovChainSampler::config()
    {
        return _imp->config;
    }

    /* MarkovChainSampler::Config */

    MarkovChainSampler::Config::Config() :
        number_of_chains(1, std::numeric_limits<unsigned>::max(), 3),
        seed(0),
        parallelize(true),
        min_efficiency(0,1,0.15), // compatible with BAT defaults
        max_efficiency(0,1,0.50), // compatible with BAT defaults
        rvalue_criterion_param(1, 100, 1.1),
        rvalue_criterion_posterior(1, 100, 1.1),
        use_strict_rvalue_definition(true),
        use_posterior_rvalue(true),
        scale_initial(1),
        scale_min(std::pow(2, -15)),
        scale_max(3),
        need_prerun(true),
        prerun_iterations_update(1000),
        prerun_iterations_min(prerun_iterations_update),
        prerun_iterations_max(1e6),
        chunks(100),
        chunk_size(1000),
        store(true)
    {
    }

    MarkovChainSampler::Config
    MarkovChainSampler::Config::Default()
    {
        return MarkovChainSampler::Config();
    }

    MarkovChainSampler::Config
    MarkovChainSampler::Config::Quick()
    {
        MarkovChainSampler::Config config;

        config.number_of_chains = 1;
        config.use_strict_rvalue_definition = false;
        config.use_posterior_rvalue = false;
        config.scale_initial = 0.4;
        config.need_prerun = true;
        config.prerun_iterations_max = 5e4;
        config.prerun_iterations_update = 400;
        config.prerun_iterations_min = config.prerun_iterations_update;
        config.chunks = 10;
        config.chunk_size = 100;

        return config;
    }
}
