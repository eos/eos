/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2011, 2013 Danny van Dyk
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

#include <config.h>

#include <eos/statistics/log-posterior.hh>
#include <eos/statistics/markov-chain-sampler.hh>
#include <eos/statistics/rvalue.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/thread_pool.hh>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#include <algorithm>
#include <map>
#include <limits>
#include <sys/stat.h>

namespace eos
{
    template<>
    struct Implementation<MarkovChainSampler>
    {
        // the target density to sample from
        DensityPtr density;

        // our configuration options
        MarkovChainSampler::Config config;

        // tickets for parallel computations
        std::vector<Ticket> tickets;

        // number of scan parameters
        unsigned number_of_parameters;

        // independent chains
        std::vector<MarkovChain> chains;

        //todo remove? We would need one for each cluster actually
        // collected information regarding the prerun
        MarkovChainSampler::PreRunInfo pre_run_info;

        ChainGroup::RValueFunction compute_rvalue;

        Implementation(const DensityPtr & density, const MarkovChainSampler::Config & config) :
            density(density),
            config(config),
            compute_rvalue(config.use_strict_rvalue_definition ? &RValue::gelman_rubin : &RValue::approximation)
        {
            initialize();
        }

        /*
         * Checks efficiencies, adjusts if needed
         * return true if all efficiencies in ranges defined by MarkovChainConfig::min_efficiency, MarkovChainConfig::max_efficiency
         *
         * @param iterations Consider this many iterations (counting from the end of the history) for efficiency and adaption
         */
        bool adjust_scales(const unsigned & iterations)
        {
            bool efficiencies_ok = true;

            // loop over chains and proposal functions
            for (unsigned c = 0 ; c < chains.size() ; ++c)
            {
                const MarkovChain::Stats & statistics = chains[c].statistics();
                if (chains[c].history().states.empty())
                {
                    throw InternalError("MarkovChainSampler::adjust_scales: cannot adapt from empty history");
                }

                // rely on the fact that counters are reset in each chunk
                double efficiency = 1.0 * statistics.iterations_accepted / (statistics.iterations_accepted + statistics.iterations_rejected);
                if ((efficiency < config.min_efficiency) || (efficiency > config.max_efficiency))
                    efficiencies_ok = false;

                // consider only the last chunk
                MarkovChain::State::Iterator states_begin = chains[c].history().states.end() - iterations;
                MarkovChain::State::Iterator states_end = chains[c].history().states.end();

                chains[c].proposal_function()->adapt(states_begin, states_end, efficiency, config.min_efficiency, config.max_efficiency);

                Log::instance()->message("markov_chain_sampler.efficiencies", ll_debug)
                        << "Current efficiency for chain " << c << ": " << stringify(efficiency, 4);

                Log::instance()->message("markov_chain_sampler.efficiencies", ll_debug)
                        << "invalid/rejected proposals = " << stringify(1.0 * statistics.iterations_invalid / statistics.iterations_rejected, 4);
            }
            if (efficiencies_ok)
                Log::instance()->message("markov_chain_sampler.efficiencies", ll_informational)
                    << "All efficiencies OK";

            return efficiencies_ok;
        }

        bool check_convergence(const unsigned & iterations)
        {
            bool efficiencies_ok = adjust_scales(iterations);

            bool rvalues_ok = true;

            // no R-value for single chain
            if ( chains.size() > 1 )
            {
                rvalues_ok = check_rvalues();
            }

            if (efficiencies_ok && rvalues_ok)
            {
                Log::instance()->message("markov_chain_sampler.convergence", ll_informational)
                    << "Convergence achieved";

                return true;
            }

            return false;
        }

        /*!
         * Check Gelman-Rubin R-value in clusters separately
         * @return
         */
        bool check_rvalues()
        {
            bool all_rvalues_small = true;

            // calculate statistics
            std::vector<std::vector<double>> all_chains_means;
            std::vector<std::vector<double>> all_chains_variances;

            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end; ++c)
            {
                std::vector<double> means, variances;
                MarkovChain::State::Iterator it_begin = c->history().states.cbegin();
                it_begin += unsigned(config.skip_initial * c->history().states.size());
                c->history().mean_and_variance(it_begin, c->history().states.cend(), means, variances);
                all_chains_means.push_back(means);
                all_chains_variances.push_back(variances);
            }

            // loop over all parameters to check and get R-values
            for (unsigned pmtr = 0 ; pmtr < number_of_parameters ; ++pmtr)
            {
                // keep subset of chain statistics in here
                std::vector<double> chain_means, chain_variances;

                // read out statistics
                for (unsigned c = 0 ; c < config.number_of_chains ;  ++c)
                {
                    chain_means.push_back(all_chains_means[c][pmtr]);
                    chain_variances.push_back(all_chains_variances[c][pmtr]);
                }

                double rvalue = compute_rvalue(chain_means, chain_variances, pre_run_info.iterations);
                pre_run_info.rvalue_parameters[pmtr] = rvalue;

                if (rvalue > config.rvalue_criterion_param || std::isnan(rvalue))
                {
                    all_rvalues_small = false;

                    Log::instance()->message("markov_chain_sampler.parameter_rvalue_too_large", ll_informational)
                        << "R-value of parameter '" << chains.front().parameter_descriptions()[pmtr].parameter->name()
                        << "' is too large: "
                        << rvalue << " > " << config.rvalue_criterion_param;
                }
            }

            // posterior converged?
#if 0
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
                if ((pre_run_info.rvalue_posterior > config.rvalue_criterion_posterior) || std::isnan(pre_run_info.rvalue_posterior))
                {
                    all_rvalues_small = false;

                    Log::instance()->message("markov_chain_sampler.posterior_rvalue_too_large", ll_informational)
                        << "R-value of posterior is too large: "
                        << pre_run_info.rvalue_posterior << " > " << config.rvalue_criterion_posterior;
                }
            }
#endif
            if (all_rvalues_small)
            {
                Log::instance()->message("markov_chain_sampler.convergence", ll_informational)
                    << "All R-values OK";
            }

            return all_rvalues_small;
        }

        void check_rvalues_main()
         {
            if (chains.size() < 2)
                return;

            Log::instance()->message("markov_chain_sampler.convergence", ll_informational)
                << "Checking R-values for the last chunk of size " << config.chunk_size;

             bool all_rvalues_small = true;

             // calculate statistics
            std::vector<std::vector<double>> all_chains_means;
            std::vector<std::vector<double>> all_chains_variances;
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end; ++c)
            {
                std::vector<double> means, variances;
                MarkovChain::State::Iterator it_begin = c->history().states.cend() - config.chunk_size;
                c->history().mean_and_variance(it_begin, c->history().states.cend(), means, variances);
                all_chains_means.push_back(means);
                all_chains_variances.push_back(variances);
            }

             // loop over all parameters to check and get R-values
             for (unsigned par = 0 ; par < number_of_parameters ; ++par)
             {
                 // keep subset of chain statistics for single parameter in here
                 std::vector<double> chain_means, chain_variances;

                 // read out statistics
                 for (unsigned c = 0 ; c < chains.size() ;  ++c)
                 {
                     chain_means.push_back(all_chains_means[c][par]);
                     chain_variances.push_back(all_chains_variances[c][par]);
                 }

                 double rvalue = compute_rvalue(chain_means, chain_variances, config.chunk_size);

                 if (rvalue > config.rvalue_criterion_param || std::isnan(rvalue))
                 {
                     all_rvalues_small = false;

                     Log::instance()->message("markov_chain_sampler.main_run", ll_informational)
                         << "R-value of parameter '" << chains.front().parameter_descriptions()[par].parameter->name() << "' is too large: "
                         << rvalue << " > " << config.rvalue_criterion_param;
                 }
             }

             if (all_rvalues_small)
             {
                 Log::instance()->message("markov_chain_sampler.main_run", ll_informational)
                     << "All R-values OK";
             }
         }


        /*
         * Dump MCMC samples and proposal density state to HDF5 file.
         *
         * @param output_base The root directory name within the HDF5 file under which all samples are stored.
         */
        void dump_hdf5(const std::string & output_base, const unsigned & last_iterations)
        {
            hdf5::File file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

            Log::instance()->message("markov_chain_sampler.dump_hdf5", ll_debug)
                << "Dumping all " << chains.size() <<" chains to HDF5 file " << config.output_file;

            unsigned i = 0;
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++i)
            {
                c->dump_history(file, output_base + "/chain #" + stringify(i), last_iterations);
                c->dump_proposal(file, output_base + "/chain #" + stringify(i));
            }
        }

        // common method to call from multiple constructors
        void initialize()
        {
            // the number of scan and nuisance parameters
            number_of_parameters = std::distance(density->begin(), density->end());

            // proposal covariance
            if (config.proposal_initial_covariance.size() != power_of<2>(number_of_parameters))
            {
                Log::instance()->message("markov_chain_sampler.initialize", ll_informational)
                    << "Determining initial proposal covariance assuming flat priors";

                config.proposal_initial_covariance.assign(power_of<2>(number_of_parameters), 0.0);

                unsigned par = 0;
                for (auto & def : *density)
                {
                    config.proposal_initial_covariance[par + number_of_parameters * par] =
                            power_of<2>(def.max - def.min) / 12.0;
                    ++par;
                }
            }

            // create independent chains -> different seeds
            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, config.seed);

            /* setup chains */

            for (unsigned c = 0 ; c < config.number_of_chains ; ++c)
            {
                // todo draw initial point from prior

                //todo use factory for proposal_function
                std::shared_ptr<MarkovChain::ProposalFunction> prop;
                if (config.proposal == "MultivariateGaussian")
                {
                    if (c == 0)
                    {
                        Log::instance()->message("markov_chain_sampler.initialize", ll_informational)
                            << "Using proposal_functions::MultivariateGaussian";
                    }
                    prop.reset(new proposal_functions::MultivariateGaussian(number_of_parameters, config.proposal_initial_covariance,
                                                                            config.scale_automatic));
                }
                if (config.proposal == "MultivariateStudentT")
                {
                    if (c == 0)
                    {
                        Log::instance()->message("markov_chain_sampler.initialize", ll_informational)
                            << "Using proposal_functions::MultivariateStudentT";
                    }
                    prop.reset(new proposal_functions::MultivariateStudentT(number_of_parameters, config.proposal_initial_covariance,
                                                                            config.student_t_degrees_of_freedom, config.scale_automatic));
                }
                // default behavior
                if (! prop)
                {
                    if (c == 0)
                    {
                        Log::instance()->message("markov_chain_sampler.initialize", ll_warning)
                                        << "No proposal function of name '" << config.proposal << "' registered."
                                        << "Falling back to MultivariateGaussian.";
                    }
                    prop.reset(new proposal_functions::MultivariateGaussian(number_of_parameters, config.proposal_initial_covariance,
                                                                            config.scale_automatic));
                }

                MarkovChain chain(density, config.seed + c, prop);
                chains.push_back(chain);
            }
            gsl_rng_free(rng);

            // setup prerun info
            pre_run_info =
            MarkovChainSampler::PreRunInfo
            {
                false, 0,0, std::numeric_limits<double>::max(),
                std::vector<double>(std::distance(density->begin(), density->end()),
                std::numeric_limits<double>::max())
            };
        }

        /*
         * Collect samples from posterior and check for convergence.
         */
        void pre_run()
        {
            Log::instance()->message("markov_chain_sampler.prerun_start", ll_informational)
                << "Commencing the pre-run with " << config.prerun_iterations_min << ", "
                << config.prerun_iterations_max << ", " << config.prerun_iterations_update
                << " (min, max, update) iterations.";

            // write parameter descriptions
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);
                for (unsigned i = 0; i < chains.size(); ++i)
                {
                    density->dump_descriptions(file, "/descriptions/prerun/chain #" + stringify(i));
                }
            }

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

            unsigned number_of_updates = 0;
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
                number_of_updates++;

                // store state before adjusting proposal!
                if (config.store_prerun)
                    dump_hdf5("/prerun", config.prerun_iterations_update);

                // check efficiency in last chunk and overall R-value: typically changes proposal
                pre_run_info.converged = check_convergence(config.prerun_iterations_update);

                Log::instance()->message("markov_chain_sampler.prerun_progress", ll_informational)
                    << "Pre-run has completed " << pre_run_info.iterations << " iterations";
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
                Log::instance()->message("markov_chain_sampler.no_convergence", ll_warning)
                    << "Pre-run did NOT converge!";
            }
        }

        /*
         * Collect samples from posterior for further analysis.
         * No convergence checking here.
         * It is assumed that chains, including their proposal,
         * are set up already.
         */
        void main_run()
        {
            Log::instance()->message("markov_chain_sampler.mainrun_start", ll_informational)
                << "Commencing the main-run";

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
                    dump_hdf5("/main run", config.chunk_size);
                }

                check_rvalues_main();

                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                {
                    double efficiency = 1.0 * c->statistics().iterations_accepted / (c->statistics().iterations_accepted +  c->statistics().iterations_rejected);

                    Log::instance()->message("markov_chain_sampler.mainrun_efficiencies", ll_debug)
                            << "Current efficiency for chain " << std::distance(chains.begin(), c) << ": " << efficiency;

                    Log::instance()->message("markov_chain_sampler.mainrun_invalid", ll_debug)
                            << "invalid/rejected proposals = " << 1.0 * c->statistics().iterations_invalid / c->statistics().iterations_rejected;
                }

                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                {
                    c->clear();
                }

            }
            Log::instance()->message("markov_chain_sampler.mainrun_end", ll_informational)
                << "Finished the main-run";
        }

        void run()
        {
            // overwrite file only if sampling is requested
            setup_output();

            if (config.need_prerun)
            {
                pre_run();
            }

            if (config.need_main_run)
            {
                // set up chains
                setup_main_run();

                main_run();
            }
        }

        /*
         *
         */
        void setup_main_run()
        {
            //  clear up
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                c->clear();

                // save history?
                c->keep_history(config.store);
            }

            // write parameter descriptions
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);
                for (unsigned i = 0; i < chains.size(); ++i)
                {
                    density->dump_descriptions(file, "/descriptions/main run/chain #" + stringify(i));
                }
            }
        }

        void setup_output()
        {
            if (config.output_file.empty())
            {
                Log::instance()->message("markov_chain_sampler.setup_output", ll_warning)
                    << "No output file specified, results of sampling will not be stored!";
            }

            //  overwrite existing file
            hdf5::File::Create(config.output_file);
        }
    };

    MarkovChainSampler::MarkovChainSampler(const DensityPtr & density, const MarkovChainSampler::Config & config) :
        PrivateImplementationPattern<MarkovChainSampler>(new Implementation<MarkovChainSampler>(density, config))
    {
    }

    MarkovChainSampler::~MarkovChainSampler()
    {
    }

    std::vector<HistoryPtr>
    MarkovChainSampler::read_chains(const std::vector<std::shared_ptr<hdf5::File>> & input_files, std::string base)
    {
        std::vector<HistoryPtr> result;

        // loop over files
        bool found_chain = false;
        for (auto f = input_files.begin(), f_end = input_files.end() ; f != f_end ; ++f)
        {
            std::string group_name;
            // loop over chains
            unsigned c = 0;
            while (true)
            {
                group_name = base + "/chain #" + stringify(c);
                if (! (*f)->group_exists(group_name))
                    break;

                found_chain = true;

                HistoryPtr history(new MarkovChain::History());
                ProposalFunctionPtr prop;
                std::string proposal_type;
                MarkovChain::Stats stat;

                // fill the objects
                MarkovChain::read_data(**f, group_name, *history, prop, proposal_type, stat);

                // store them
                result.push_back(history);

                // go to next chain
                ++c;
            }
        }
        if (! found_chain)
        {
            throw InternalError("read_chains: Did not find any usable data in the files given");
        }
        return result;
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
        number_of_chains(1, std::numeric_limits<unsigned>::max(), 4),
        seed(0),
        parallelize(true),
        min_efficiency(0, 1, 0.15), // incompatible with BAT defaults [0.15, 0.5]
        max_efficiency(0, 1, 0.35),
        rvalue_criterion_param(1, 100, 1.1),
        rvalue_criterion_posterior(1, 100, 1.1),
        use_strict_rvalue_definition(true),
        use_posterior_rvalue(false),
        scale_automatic(true),
        need_prerun(true),
        prerun_iterations_update(1000),
        prerun_iterations_min(prerun_iterations_update),
        prerun_iterations_max(1e6),
        proposal("MultivariateGaussian"),
        student_t_degrees_of_freedom(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max(), 1.0),
        store_prerun(true),
        adapt_iterations(0),
        chunks(100),
        chunk_size(1000),
        need_main_run(true),
        skip_initial(0, 1, 0.1),
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
        config.need_prerun = true;
        config.prerun_iterations_max = 1e5;
        config.prerun_iterations_update = 400;
        config.prerun_iterations_min = config.prerun_iterations_update;
        config.chunks = 10;
        config.chunk_size = 100;

        return config;
    }

    std::ostream & operator<<(std::ostream & stream, const MarkovChainSampler::Config & c)
    {
        stream << std::boolalpha
               << "Prerun settings:" << std::endl
               << "nchains = " << c.number_of_chains
               << ", seed = " << c.seed
               << ", parallelize = " << c.parallelize
               << ", prerun min iterations = " << c.prerun_iterations_min << std::endl
               << ", prerun max iterations = " << c.prerun_iterations_max
               << ", prerun update iterations = " << c.prerun_iterations_update
               << ", skip initial = " << c.skip_initial;
        return stream;
    }
}
