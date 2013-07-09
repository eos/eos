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

#include <eos/utils/markov_chain_sampler.hh>
#include <config.h>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/rvalue.hh>
#include <eos/utils/thread_pool.hh>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#include <algorithm>
#include <map>
#include <limits>
#include <sys/stat.h>

namespace eos
{
    // TODO replace by lambda function
    struct SameName
    {
        const std::string name;
        SameName(const std::string & name) :
            name(name)
        {
        }

        bool operator() (const std::tuple<std::string, double, double> & t)
        {
            return name == std::get<0>(t);
        }
    };

    template<>
    struct Implementation<MarkovChainSampler>
    {
        // Worker allows simple thread parallelization of massive mode finding
        struct Worker
        {
            std::shared_ptr<Analysis> analysis;

            std::shared_ptr<ROOT::Minuit2::FunctionMinimum> minimum;

            Worker(const Analysis & analysis) :
                analysis(analysis.clone())
            {
            }

            /*
             * Return empty vector if something went wrong
             */
            std::vector<double> mode() const
            {
                std::vector<double> result;

                if (! minimum)
                {
                    Log::instance()->message("Worker.mode", ll_warning)
                        << "No search conducted yet.";
                    return result;
                }

                if (! minimum->IsValid())
                {
                    return result;
                }

                for (unsigned i = 0 ; i < analysis->parameter_descriptions().size() ; ++i)
                {
                    result.push_back(minimum->UserParameters().Value(i));
                }

                return result;
            }

            void optimize(std::vector<double> initial_point, const Analysis::OptimizationOptions & options)
            {
                Log::instance()->message("Worker.optimize", ll_informational)
                    << "Starting minuit optimization at "
                    << stringify(initial_point.cbegin(), initial_point.cend());

                minimum.reset(new ROOT::Minuit2::FunctionMinimum(analysis->optimize_minuit(initial_point, options)));

                Log::instance()->message("Worker.optimize", ll_informational)
                    << "Finished minuit optimization";
                Log::instance()->message("Worker.optimize", ll_debug)
                    << print_status();
            }

            std::string print_status()
            {
                std::stringstream lhs;

                if ( ! minimum)
                {
                    return lhs.str();
                }

                const ROOT::Minuit2::FunctionMinimum & m = *minimum;

                std::vector<double> parameters(m.UserParameters().Params());
                if (m.IsValid())
                {
                    lhs << "|Success|: found mode after " << m.NFcn() << " calls with log(post) at "
                        << stringify(parameters.begin(), parameters.end())
                        << " = " << -1.0 * m.Fval() << "; ";
                    return lhs.str();
                }

                lhs << "|Failure|, stopped after " << m.NFcn() << " calls with log(post) at "
                    << stringify(parameters.begin(), parameters.end())
                    << " = " << -1.0 * m.Fval() <<", listing the symptoms: ";

                // look for a reason of failure
                if (! m.HasValidParameters())
                {
                    lhs << "invalid parameters; ";
                }
                if (! m.HasValidCovariance())
                {
                    lhs << "invalid covariance; ";
                }
                if (! m.HasAccurateCovar())
                {
                    lhs << "inaccurate covariance; ";
                }
                if (! m.HasPosDefCovar())
                {
                    lhs << "covariance not positive definite; ";
                }
                if (m.HasMadePosDefCovar())
                {
                    lhs << "covariance was made positive definite; ";
                }
                if (m.HesseFailed())
                {
                    lhs << "Hesse failed; ";
                }
                if (! m.HasCovariance())
                {
                    lhs << "has no covariance; ";
                }
                if (m.IsAboveMaxEdm())
                {
                    lhs << "estimated distance to minimum " << m.Edm() << " too large; ";
                }
                if (m.HasReachedCallLimit())
                {
                    lhs << "exceeded function call limit with " << m.NFcn() << " calls; ";
                }
                return lhs.str();
            }
        };

        // store reference, but don't own analysis
        const Analysis & analysis;

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
            if ( (config.partitions.size() >= 2 && config.prerun_chains_per_partition > 1) ||
                 (config.partitions.size() == 1 && chains.size() > 1) )
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
                // keep subset of chain statistics in on partition in here
                std::vector<double> chain_means, chain_variances;

                for (unsigned part = 0 ; part < config.partitions.size() ; ++part)
                {
                    chain_means.clear();
                    chain_variances.clear();

                    // read out statistics
                    for (unsigned c = 0 ; c < config.prerun_chains_per_partition ;  ++c)
                    {
                        chain_means.push_back(all_chains_means[part * config.prerun_chains_per_partition + c][pmtr]);
                        chain_variances.push_back(all_chains_variances[part * config.prerun_chains_per_partition + c][pmtr]);
                    }

                    double rvalue = compute_rvalue(chain_means, chain_variances, pre_run_info.iterations);
                    pre_run_info.rvalue_parameters[pmtr] = rvalue;

                    if (rvalue > config.rvalue_criterion_param || std::isnan(rvalue))
                    {
                        all_rvalues_small = false;

                        Log::instance()->message("markov_chain_sampler.parameter_rvalue_too_large", ll_informational)
                            << "R-value of parameter '" << chains.front().parameter_descriptions()[pmtr].parameter.name()
                            << "' in partition " << part << " is too large: "
                            << rvalue << " > " << config.rvalue_criterion_param;
                    }
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
                         << "R-value of parameter '" << chains.front().parameter_descriptions()[par].parameter.name() << "' is too large: "
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
            number_of_parameters = analysis.parameter_descriptions().size();

            // create independent chains -> different seeds
            std::vector<double> covariance(number_of_parameters * number_of_parameters, 0.0);

            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, config.seed);

            /* setup chains */

            // add one partition that contains the entire parameter cube
            if (config.partitions.empty())
            {
                std::vector<std::tuple<std::string, double, double>> partition;

                partition.push_back(std::make_tuple(
                    analysis.parameter_descriptions().begin()->parameter.name(),
                    analysis.parameter_descriptions().begin()->min,
                    analysis.parameter_descriptions().begin()->max));

                config.partitions.push_back(partition);

                // need to ensure that we have the same number of chains in prerun and main run
                if (config.need_main_run)
                {
                    config.prerun_chains_per_partition = config.number_of_chains;
                }
            }

            // loop over partitions
            for (auto p = config.partitions.cbegin(), p_end = config.partitions.cend() ; p != p_end ; ++p)
            {
                AnalysisPtr ana = analysis.clone();

                for (auto r = p->cbegin(), r_end = p->cend() ; r != r_end ; ++r)
                {
                    ana->restrict(std::get<0>(*r), std::get<1>(*r), std::get<2>(*r));
                }
                for (unsigned c = 0 ; c < config.prerun_chains_per_partition ; ++c)
                {
                    // initialize MVG with means/variance from priors on the diagonal
                    for (unsigned par = 0 ; par < number_of_parameters ; ++par)
                    {
                        std::string name = analysis.parameter_descriptions()[par].parameter.name();
                        LogPriorPtr prior = analysis.log_prior(name);
                        covariance[par + number_of_parameters * par] = prior->variance();

                        // rescale the variance of partitioned parameters by (new_range/old_range)^2.
                        auto r = std::find_if(p->begin(), p->end(), SameName(name));
                        if (p->end() != r)
                        {
                            double min = analysis.parameter_descriptions()[par].min;
                            double max = analysis.parameter_descriptions()[par].max;
                            covariance[par + number_of_parameters * par] *= power_of<2>((std::get<2>(*r) - std::get<1>(*r)) / (max - min));
                        }

                        // rescale variance of scan parameters with a configurable value, in order
                        // to avoid drawing too many samples outside the allowed range.
                        if (! ana->parameter_descriptions()[par].nuisance || config.scale_nuisance)
                        {
                            covariance[par + number_of_parameters * par] /= power_of<2>(config.scale_reduction);
                        }
                    }
                    //todo use factory for proposal_function
                    std::shared_ptr<MarkovChain::ProposalFunction> prop;
                    if (config.proposal == "MultivariateGaussian")
                    {
                        if (p == config.partitions.cbegin() && c == 0)
                        {
                            Log::instance()->message("markov_chain_sampler.initialize", ll_informational)
                                << "Using proposal_functions::MultivariateGaussian";
                        }
                        prop.reset(new proposal_functions::MultivariateGaussian(number_of_parameters, covariance, config.scale_automatic));
                    }
                    if (config.proposal == "MultivariateStudentT")
                    {
                        if (p == config.partitions.cbegin() && c == 0)
                        {
                            Log::instance()->message("markov_chain_sampler.initialize", ll_informational)
                                << "Using proposal_functions::MultivariateStudentT";
                        }
                        prop.reset(new proposal_functions::MultivariateStudentT(number_of_parameters, covariance,
                                config.student_t_degrees_of_freedom, config.scale_automatic));
                    }
                    if (! config.block_proposal_parameters.empty())
                    {
                        if (p == config.partitions.cbegin() && c == 0)
                        {
                            Log::instance()->message("markov_chain_sampler.initialize", ll_informational)
                                        << "Using proposal_functions::BlockDecomposition";
                        }

                        // empty initially
                        proposal_functions::BlockDecomposition* bd = new proposal_functions::BlockDecomposition();

                        // which parameters need special treatment
                        std::vector<bool> parameter_proposal_list;

                        // count pars for multivariate proposal
                        std::vector<unsigned> index_list;

                        // find runs of consecutive parameters
                        for (unsigned par = 0 ; par < number_of_parameters ; ++par)
                        {
                            // is it a special parameter?
                            auto res = std::find(config.block_proposal_parameters.begin(),
                                    config.block_proposal_parameters.end(),
                                    ana->parameter_descriptions()[par].parameter.name());
                            if (res != config.block_proposal_parameters.end())
                            {
                                parameter_proposal_list.push_back(true);
                            }
                            else
                            {
                                parameter_proposal_list.push_back(false);
                                index_list.push_back(par);
                            }
                        }
                        // last points to the conceptual end of new range, but the actual size is untouched
                        auto parameter_proposal_list_copy(parameter_proposal_list);
                        auto last = std::unique(parameter_proposal_list_copy.begin(), parameter_proposal_list_copy.end());
                        if (std::distance(parameter_proposal_list_copy.begin(), last) > 2 || (std::distance(parameter_proposal_list_copy.begin(), last) == 2 && parameter_proposal_list_copy.front()) )
                        {
                            Log::instance()->message("MC_sampler.initialize_decomposition", ll_debug)
                                        << "parameter_proposal_list_copy: " << stringify_container(parameter_proposal_list_copy);
                            throw InternalError("With block decomposition, all parameters with fixed 1D proposal should come after the parameters with a Multivariate proposal");
                        }


                        // flush the accumulated scan parameters
                        if ( ! index_list.empty())
                        {
                            std::vector<double> covariance(power_of<2>(index_list.size()), 0.0);
                            // initialize multivariate with means/variance from priors on the diagonal
                            for (auto multi_par = index_list.cbegin(), scan_par_end = index_list.cend() ;
                                    multi_par != scan_par_end ; ++multi_par)
                            {
                                std::string name = ana->parameter_descriptions()[*multi_par].parameter.name();
                                LogPriorPtr prior = ana->log_prior(name);
                                covariance[*multi_par + index_list.size() * (*multi_par)] = prior->variance();

                                // rescale the variance of partitioned parameters by (new_range/old_range)^2.
                                auto r = std::find_if(p->begin(), p->end(), SameName(name));
                                if (p->end() != r)
                                {
                                    double min = ana->parameter_descriptions()[*multi_par].min;
                                    double max = ana->parameter_descriptions()[*multi_par].max;
                                    covariance[*multi_par + index_list.size() * *multi_par] *=
                                            power_of<2>((std::get<2>(*r) - std::get<1>(*r)) / (max - min));
                                }

                                // rescale variance of scan parameters with a configurable value, in order
                                // to avoid drawing too many samples outside the allowed range.
                                covariance[*multi_par + index_list.size() * *multi_par] /=
                                        power_of<2>(config.scale_reduction);
                            }
                            Log::instance()->message("MC_sampler.initialize_decomposition", ll_debug)
                                                    << "Add scan block with " << index_list.size() << " dimensions";

                            // we can construct the multivariate proposal and add it as one block
                            proposal_functions::MultivariateProposalPtr mv;
                            if (config.proposal == "MultivariateGaussian")
                                mv.reset(new proposal_functions::MultivariateGaussian(index_list.size(), covariance, config.scale_automatic));
                            if (config.proposal == "MultivariateStudentT")
                                mv.reset(new proposal_functions::MultivariateStudentT(index_list.size(), covariance,
                                        config.student_t_degrees_of_freedom, config.scale_automatic));
                            if ( ! mv)
                                throw InternalError("Invalid local proposal function: " + config.proposal);
                            bd->add(mv);
                            index_list.clear();
                        }

                        // add 1D proposals
                        for (unsigned par = 0 ; par < number_of_parameters ; ++par)
                        {
                            if (parameter_proposal_list[par])
                            {
                                std::string name = ana->parameter_descriptions()[par].parameter.name();
                                LogPriorPtr prior = ana->log_prior(name);
                                if (! prior)
                                    throw InternalError("MC_sampler.initialize_decomposition: Undefined prior for '" + name + "'");
                                // is cloned internally, so no back propagation
                                bd->add(prior);
                            }
                        }
                        prop.reset(bd);
                    }
                    // default behavior
                    if (! prop)
                    {
                        if (p == config.partitions.cbegin() && c == 0)
                        {
                            Log::instance()->message("markov_chain_sampler.initialize", ll_warning)
                                            << "No proposal function of name '" << config.proposal << "' registered."
                                            << "Falling back to MultivariateGaussian.";
                        }
                        prop.reset(new proposal_functions::MultivariateGaussian(number_of_parameters, covariance, config.scale_automatic));
                    }

                    MarkovChain chain(*ana, config.seed + config.prerun_chains_per_partition * std::distance(config.partitions.cbegin(), p) + c, prop);
                    chains.push_back(chain);
                }
            }
            gsl_rng_free(rng);

            // setup prerun info
            pre_run_info = MarkovChainSampler::PreRunInfo
                    {
                false, 0,0, std::numeric_limits<double>::max(), std::vector<double>(analysis.parameter_descriptions().size(),std::numeric_limits<double>::max())
                    };
        }

        void massive_mode_finding(const Analysis::OptimizationOptions & options, bool dump = false)
        {
            if (options.mcmc_pre_run)
            {
                setup_output();
                pre_run();
            }

            // create workers
            std::vector<std::shared_ptr<Worker>> workers;

            // let them work
            // tickets for parallel computations
            std::vector<Ticket> tickets;

            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                // parallelize computations
                auto worker = std::make_shared<Worker>(analysis);
                workers.push_back(worker);

                // extract best point seen during prerun
                const std::vector<double> & starting_point = c->statistics().parameters_at_mode;

                if (config.parallelize)
                {
                    tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&Worker::optimize, worker.get(), starting_point, options)));
                }
                else
                {
                    worker->optimize(starting_point, options);
                }
            }

            // wait for job completion
            for (auto t = tickets.begin(), t_end = tickets.end() ; t != t_end ; ++t)
            {
                t->wait();
            }
            tickets.clear();

            if (dump)
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);

                auto w = workers.cbegin();
                unsigned i = 0;
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++w, ++i)
                {
                    // Ignore validity, but use only if it is better than MCMC value
                    if ( -1.0 * (**w).minimum->Fval() > c->statistics().mode_of_posterior)
                    {
                        // switch from minimum to maximum again
                        c->set_mode(file, "/prerun/chain #" + stringify(i), (**w).minimum->UserParameters().Params(), -1.0 * (**w).minimum->Fval());
                    }
                }
            }

            /* find unique modes */
            std::vector<std::vector<double>> unique_modes;
            // sort by posterior value
            std::multimap<double, unsigned> posterior_worker_index;
            unsigned invalid = 0;

            // add first mode
            auto w = workers.begin();

            for (auto w_end = workers.end() ; w != w_end ; ++w)
            {
                auto mode = (**w).mode();
                if( mode.empty() )
                {
                    ++invalid;
                    continue;
                }

                // first valid mode
                if (unique_modes.empty())
                {
                    unique_modes.push_back(mode);
                    posterior_worker_index.insert(std::pair<double, unsigned>(-1.0 * (**w).minimum->Fval(), 0));
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

                    if (std::sqrt(difference) < options.splitting_tolerance)
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
                    posterior_worker_index.insert(std::pair<double, unsigned>(-1.0 * (**w).minimum->Fval(), unique_modes.size() - 1));
                }
            }

            // print out covariances and modes
            for (auto i = posterior_worker_index.cbegin(), i_end = posterior_worker_index.cend() ; i != i_end ;  ++i)
            {
                Log::instance()->message("MC_sampler.mode_finding", ll_debug)
                             << "worker " << i->second << ", unique: "
                             << std::distance(posterior_worker_index.cbegin(), i) << workers[i->second]->minimum->UserState();
            }
            for (auto i = posterior_worker_index.cbegin(), i_end = posterior_worker_index.cend() ; i != i_end ;  ++i)
            {
                const std::vector<double> & m = unique_modes[i->second];
                Log::instance()->message("MC_sampler.mode_finding", ll_informational)
                             << "log(post) at "
                             << stringify(m.cbegin(), m.cend())
                             << " = " << i->first;
            }

            Log::instance()->message("MC_sampler.mode_finding", ll_informational)
                     << "Identified " << unique_modes.size() << " unique mode(s) of posterior,"
                     << " minuit failed " << invalid << " times.";

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
                unsigned i = 0;
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++i)
                {
                    c->dump_description(file, "/descriptions/prerun/chain #" + stringify(i));
                }
            }

            // require:chains are initialized
            pre_run_info.converged = false;
            pre_run_info.iterations = 0;

            // set up chains
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                // save history
                c->keep_history(true, config.store_observables_and_proposals);
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

                if (config.number_of_chains < 2 || (config.global_local_config && config.prerun_chains_per_partition < 2))
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

                    if (config.global_local_config && (chunk + 1) * config.chunk_size < config.adapt_iterations)
                    {
                        c->proposal_function()->adapt(c->history().states.end() - config.chunk_size, c->history().states.end(),
                                                      efficiency, config.min_efficiency, config.max_efficiency);
                    }
                }

                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
                {
                    c->clear();
                }

            }
            Log::instance()->message("markov_chain_sampler.mainrun_end", ll_informational)
                << "Finished the main-run";
        }

        void resume(const hdf5::File & file)
        {
            Log::instance()->message("markov_chain_sampler.resume", ll_informational)
                << "Copying settings from " << file.name();

            // the goal is to start the main run immediately when run() is invoked
            config.need_prerun = false;

            // read in proposal
            ProposalFunctionPtr prop;
            try
            {
                prop = proposal_functions::Factory::make(const_cast<hdf5::File &>(file), "/global local",
                                                             "GlobalLocal", analysis.parameter_descriptions().size());
            }
            catch (HDF5Error & e)
            {
                Log::instance()->message("markov_chain_sampler.setup_global_local", ll_error)
                    << "Errors in reading from the HDF5 file can be due to a mismatch in the "
                    << "analysis definition. Check that the same number of parameters is defined now "
                    << "and when building the GlobalLocal proposal function";
                throw e;
            }
            auto gl = dynamic_cast<proposal_functions::GlobalLocal *>(prop.get());
            if ( ! gl)
                throw InternalError("MarkovChainSampler::resume: couldn't read GlobalLocal from disk");

            // reset options affecting the runtime behavior
            // all the options affecting construction have no impact
            gl->config(*config.global_local_config);

            /* seed chains */

            MarkovChain::State best_state = gl->mode();

            Log::instance()->message("markov_chain_sampler.setup_global_local", ll_debug)
                << "Found global mode at "  << best_state << " in component "
                << best_state.hyper_parameter.component;

            std::vector<MarkovChain> new_chains;
            for (unsigned c = 0 ; c < config.number_of_chains ; ++c)
            {
                MarkovChain chain(analysis, config.seed + c, prop);
                chain.set_point(best_state.point, best_state.hyper_parameter);
                new_chains.push_back(chain);
            }

            chains = new_chains;

            Log::instance()->message("markov_chain_sampler.resume", ll_debug)
                << "chains: " << chains.size();

            /* write parameter descriptions and create output file */

            //todo Fix: don't just try once, keep increasing index until file not found
            if (hdf5::File::Exists(config.output_file))
            {
                std::string old_file_name = config.output_file;
                auto dot_pos = old_file_name.find_last_of('.');
                config.output_file.insert(dot_pos, "_1");
                Log::instance()->message("markov_chain_sampler.resume", ll_warning)
                    << "File " << old_file_name << " already exists. Store data in new file "
                    << config.output_file;
            }

            auto file_out = hdf5::File::Create(config.output_file);
            unsigned i = 0;
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++i)
            {
                c->dump_description(file_out, "/descriptions/main run/chain #" + stringify(i));
            }

            Log::instance()->message("markov_chain_sampler.resume", ll_informational)
                << "Checking parameters, priors and constraints from input file vs the current analysis";

            std::vector<std::string> constraints;
            std::vector<ParameterDescription> descriptions;
            std::string hash;
            std::vector<std::string> priors;

            MarkovChain::read_descriptions(const_cast<hdf5::File &>(file), "/descriptions", descriptions, priors, constraints, hash);

            // compare all priors (includes name checking)
            {
                auto j = priors.cbegin();
                std::string prior;
                for (auto d = descriptions.cbegin(), d_end = descriptions.cend() ; d != d_end ; ++d, ++j)
                {
                    prior = analysis.log_prior(d->parameter.name())->as_string().c_str();
                    if (*j != prior)
                        throw InternalError("MarkovChainSampler::resume: mismatch of priors between " + *j + " and " + prior);
                }
            }
            // compare all constraints
            {
                auto j = constraints.cbegin();
                std::string constraint;
                auto l = const_cast<Analysis &>(analysis).log_likelihood();
                for (auto c = l.begin(), c_end = l.end() ; c != c_end ; ++c, ++j)
                {
                    constraint = c->name();
                    if ( *j != constraint)
                    {
                        throw InternalError("MarkovChainSampler::resume: constraint mismatch:"
                                            + *j + " vs " + constraint);
                    }
                    Log::instance()->message("MarkovChainSampler::resume", ll_debug)
                        << "Comparing constraint " << constraint;
                }
            }
            // compare hash
            {
                if (hash != EOS_GITHEAD)
                    Log::instance()->message("MarkovChainSampler::resume", ll_warning)
                        << "EOS version mismatch detected: " << hash << " vs " << EOS_GITHEAD;
            }

            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                // save history?
                c->keep_history(true, config.store_observables_and_proposals);
            }

            // we are ready to go
            main_run();
        }

        void run()
        {
            // overwrite file only if sampling is requested
            setup_output();
            if (config.need_prerun)
            {
                pre_run();
                if (config.find_modes)
                {
                    bool dump = true;
                    Analysis::OptimizationOptions options = Analysis::OptimizationOptions::Defaults();
                    options.fix_flat_nuisance = true;
                    // we have already done the prerun
                    options.mcmc_pre_run = false;
                    options.maximum_iterations = 4000;
                    massive_mode_finding(options, dump);
                }
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
            if (config.global_local_config)
            {
                // one common proposal function for all chains
                std::vector<HistoryPtr> histories;
                std::vector<ProposalFunctionPtr> proposals;
                std::vector<MarkovChain::Stats> stats;

                // read data from file, so resume will produce the same results
                auto file = hdf5::File::Open(config.output_file);
                unsigned i = 0;
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++i)
                {
                    std::string group_name = "/prerun/chain #" + stringify(i);
                    HistoryPtr history(new MarkovChain::History());
                    ProposalFunctionPtr prop;
                    MarkovChain::Stats stat;
                    std::string proposal_type;
                    MarkovChain::read_data(file, group_name, *history, prop, proposal_type, stat);
                    histories.push_back(history);
                    proposals.push_back(prop);
                    stats.push_back(stat);
                }
                Log::instance()->message("MCsampler::setup_global_local", ll_debug)
                    << "Using skip_initial = " << config.global_local_config->skip_initial;

                std::shared_ptr<proposal_functions::GlobalLocal> gl(
                    new proposal_functions::GlobalLocal(histories, proposals, stats, *config.global_local_config, config.prerun_chains_per_partition));

                Log::instance()->message("MCsampler::setup_global_local", ll_debug)
                    << "first chain has " << histories.front()->states.size() << " elements"
                    << ", and its first element is " <<
                    stringify(histories.front()->states.front().point.begin(), histories.front()->states.front().point.end())
                    << ", the max posterior is " << stats.front().mode_of_posterior << " at parameters "
                    << stringify(stats.front().parameters_at_mode.begin(), stats.front().parameters_at_mode.end());

                /*
                 * Seed new chains with the global mode
                 * of all chains found in the prerun
                 */

                MarkovChain::State state_at_mode = gl->mode();
                Log::instance()->message("markov_chain_sampler.setup_global_local", ll_debug)
                    << "Found global mode at "  << state_at_mode << " in component "
                    << state_at_mode.hyper_parameter.component;

                std::vector<MarkovChain> new_chains;

                for (unsigned c = 0 ; c < config.number_of_chains ; ++c)
                {
                    MarkovChain chain(analysis, config.seed + c, gl);

                     chain.set_point(state_at_mode.point, state_at_mode.hyper_parameter);

                     new_chains.push_back(chain);
                }

                chains = new_chains;
            }

            //  clear up
            for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c)
            {
                c->clear();

                // save history?
                c->keep_history(config.store, config.store_observables_and_proposals);
            }

            // write parameter descriptions
            {
                auto file = hdf5::File::Open(config.output_file, H5F_ACC_RDWR);
                unsigned i = 0;
                for (auto c = chains.begin(), c_end = chains.end() ; c != c_end ; ++c, ++i)
                {
                    c->dump_description(file, "/descriptions/main run/chain #" + stringify(i));
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

    MarkovChainSampler::MarkovChainSampler(const Analysis & analysis, const MarkovChainSampler::Config & config) :
        PrivateImplementationPattern<MarkovChainSampler>(new Implementation<MarkovChainSampler>(analysis, config))
    {
    }

    MarkovChainSampler::~MarkovChainSampler()
    {
    }

    std::vector<HistoryPtr>
    MarkovChainSampler::build_global_local(const std::string & output_file_name,
                                                const std::vector<std::shared_ptr<hdf5::File>> input_files,
                                                const proposal_functions::GlobalLocal::Config & config,
                                                AnalysisPtr analysis)
    {
        Log::instance()->message("MarkovChainSampler::build_global_local", ll_informational)
            << "Building the global local proposal function from " << input_files.size()
            << " input files, storing the result in " << output_file_name;

        /* read in chains */

        std::vector<HistoryPtr> histories_shared;
        std::vector<ProposalFunctionPtr> proposals;
        std::vector<std::string> proposal_types;
        std::vector<MarkovChain::Stats> stats;
        std::vector<std::vector<ParameterDescription>> descriptions;
        std::vector<std::vector<std::string>> priors;
        std::vector<std::vector<std::string>> constraints;
        std::vector<std::string> hashes;

        // loop over files
        bool found_chain = false;
        for (auto f = input_files.begin(), f_end = input_files.end() ; f != f_end ; ++f)
        {
            std::string group_name;
            // loop over chains
            unsigned c = 0;
            while (true)
            {
                group_name = "/prerun/chain #" + stringify(c);
                if ( ! (**f).group_exists(group_name))
                    break;

                found_chain = true;

                HistoryPtr history(new MarkovChain::History());
                ProposalFunctionPtr prop;
                std::string proposal_type;
                MarkovChain::Stats stat;
                std::vector<ParameterDescription> descr;
                std::vector<std::string> prior;
                std::vector<std::string> constraint;
                std::string hash;

                // fill the objects
                MarkovChain::read_data(**f, group_name, *history, prop, proposal_type, stat);
                MarkovChain::read_descriptions(**f, "/descriptions/" + group_name, descr, prior, constraint, hash);

                // store them
                histories_shared.push_back(history);
                proposals.push_back(prop);
                proposal_types.push_back(proposal_type);
                stats.push_back(stat);
                descriptions.push_back(descr);
                priors.push_back(prior);
                constraints.push_back(constraint);
                hashes.push_back(hash);

                ++c;
            }
        }
        if ( ! found_chain)
        {
            throw InternalError("build_global_local: Did not find any usable data in the files given");
        }

        /* find partitions: look for consecutive descriptions which match in all parameters */
        std::vector<unsigned> partition_lengths;
        {
            // first chain is first partition. Index marks beginning of chains in same partition
            std::vector<unsigned> partition_indices(1, 0);

            // start to compare first with second
            auto p = priors.cbegin() + 1;
            auto c = constraints.cbegin() + 1;
            for (auto d = descriptions.cbegin() + 1, d_end = descriptions.cend() ; d != d_end ; ++d, ++p, ++c)
            {
                bool found_new_partition = false;

                // compare parameter descriptions
                {
                    // always compare with last partition found
                    auto j = descriptions[partition_indices.back()].cbegin();
                    for (auto i = d->cbegin(), i_end = d->cend() ; i != i_end ; ++i, ++j)
                    {
                        if (i->min != j->min || i->max != j->max)// ! (*i == *j))
                        {
                            Log::instance()->message("MarkovChainSampler::build_global_local", ll_debug)
                                << "Partitions differ in " << i->parameter.name() << ", " << j->parameter.name()
                                << "min = (" << i->min << ", " << j->min << ") "
                                << "max = (" << i->max << ", " << j->max << ") "
                                << "nus = (" << i->nuisance << ", " << j->nuisance << ") "
                                << "dis  = (" << i->discrete << ", " << j->discrete << ")";
                            found_new_partition = true;
                            break;
                        }
                        if (i->parameter.name() != j->parameter.name() || i->nuisance != j->nuisance)
                            throw InternalError("MarkovChainSampler::build_global_local: parameter mismatch:"
                                                + i->parameter.name() + " vs " + j->parameter.name());
                    }
                }
                // compare all priors
                {
                    auto j = priors[partition_indices.back()].cbegin();
                    for (auto i = p->cbegin(), i_end = p->cend() ; i != i_end ; ++i, ++j)
                    {
                        if ( *i != *j)
                        {
                            throw InternalError("MarkovChainSampler::build_global_local: prior mismatch:"
                                                + *i + " vs " + *j);
                        }
                    }
                }
                // compare all constraints
                {
                    auto j = constraints[partition_indices.back()].cbegin();
                    for (auto i = c->cbegin(), i_end = c->cend() ; i != i_end ; ++i, ++j)
                    {
                        if ( *i != *j)
                        {
                            throw InternalError("MarkovChainSampler::build_global_local: constraint mismatch:"
                                                + *i + " vs " + *j);
                        }
                    }
                }
                if (found_new_partition)
                {
                    unsigned index = std::distance(descriptions.cbegin(), d);
                    partition_lengths.push_back(index - partition_indices.back());
                    partition_indices.push_back(index);
                }
            }
            Log::instance()->message("MarkovChainSampler::build_global_local", ll_informational)
                << "The parameter descriptions, priors and constraints of the chains seem to match";

            // now compare first partition against analysis
            if (analysis)
            {
                if (descriptions.front().size() > analysis->parameter_descriptions().size())
                    throw InternalError("MarkovChainSampler::build_global_local: More parameters in file ("
                                        + stringify(descriptions.front().size())
                                        + ") than in analysis (" + stringify(analysis->parameter_descriptions().size())
                                        + ")");

                unsigned i = 0;
                auto d_part = descriptions.front().cbegin();
                for (auto d = analysis->parameter_descriptions().cbegin() ; d != analysis->parameter_descriptions().cend() ; ++d, ++d_part, ++i)
                {
                    if (d_part == descriptions.front().cend())
                        throw InternalError("MarkovChainSampler::build_global_local: parameter mismatch"
                                            + std::string(" at position " + stringify(i))
                                            + ": in analysis: " + d->parameter.name()
                                            + " but no more parameters in file");

                    if (d->parameter.name() != d_part->parameter.name())
                    {
                        throw InternalError("MarkovChainSampler::build_global_local: parameter mismatch"
                                            + std::string(" at position " + stringify(i))
                                            + ": in analysis: " + d->parameter.name()
                                            + " vs in file: " + d_part->parameter.name());
                    }

                    if (d->min != d_part->min)
                    {
                        Log::instance()->message("MarkovChainSampler::build_global_local", ll_warning)
                            << "Mismatch of minimum of '" << d->parameter.name() << "': "
                            << d->min << " vs " << d_part->min;
                    }
                    if (d->max != d_part->max)
                    {
                        Log::instance()->message("MarkovChainSampler::build_global_local", ll_warning)
                            << "Mismatch of maximum of '" << d->parameter.name() << "': "
                            << d->max << " vs " << d_part->max;
                    }
                }
            }

            // add length of last partition
            partition_lengths.push_back(descriptions.size() - partition_indices.back());

            Log::instance()->message("MarkovChainSampler::build_global_local", ll_debug)
                << "Found " << partition_indices.size() << " partitions: "
                << stringify(partition_indices.cbegin(), partition_indices.cend())
                << " with sizes: " << stringify(partition_lengths.cbegin(), partition_lengths.cend());

            // check that all partitions had same number of chains
            // in that case, unique will reduce the vector to effective length one
            auto last = std::unique(partition_lengths.begin(), partition_lengths.end());
            if (last != partition_lengths.begin() + 1)
                Log::instance()->message("MarkovChainSampler::build_global_local", ll_warning)
                << "Numbers of partitions in each file didn't match: "
                << stringify(partition_indices.cbegin(), partition_indices.cend());
        }

        // check proposal types
        {
            auto last = std::unique(proposal_types.begin(), proposal_types.end());
            if (last != proposal_types.begin() + 1)
            {
                Log::instance()->message("MarkovChainSampler::build_global_local", ll_warning)
                    << "Local proposal do not match" << stringify_container(proposal_types);
            }
        }

        /* check EOS versions */
        {
            auto last = std::unique(hashes.begin(), hashes.end());
            if (last != hashes.begin() + 1)
                Log::instance()->message("MarkovChainSampler::build_global_local", ll_warning)
                    << "Hashes to not match: " << stringify(hashes.begin(), last);
        }

        Log::instance()->message("MCsampler::build_global_local", ll_debug)
            << "Using skip_initial = " << config.skip_initial;

        if (! output_file_name.empty())
        {
            // create output file
            hdf5::File file = hdf5::File::Create(output_file_name);

            // create global local proposal
            proposal_functions::GlobalLocal gl(histories_shared, proposals, stats, config, partition_lengths.front());

            // store proposal to disk
            gl.dump_state(file, "/global local");

            // assuming that all descriptions are equivalent, copy the one from the first chain in the first file
            // to the output, so when resuming, the analysis can be checked again
            input_files.front()->copy("/descriptions/prerun/chain #0", file, "/descriptions");
        }
        return histories_shared;
    }

    void
    MarkovChainSampler::massive_mode_finding(const Analysis::OptimizationOptions & options)
    {
        _imp->massive_mode_finding(options);
    }

    MarkovChainSampler::PreRunInfo
    MarkovChainSampler::pre_run_info()
    {
        return _imp->pre_run_info;
    }

    void
    MarkovChainSampler::resume(const hdf5::File & file)
    {
        _imp->resume(file);
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
        min_efficiency(0, 1, 0.15), // incompatible with BAT defaults [0.15, 0.5]
        max_efficiency(0, 1, 0.35),
        rvalue_criterion_param(1, 100, 1.1),
        rvalue_criterion_posterior(1, 100, 1.1),
        use_strict_rvalue_definition(true),
        use_posterior_rvalue(false),
        scale_automatic(true),
        scale_nuisance(true),
        scale_reduction(1),
        find_modes(false),
        need_prerun(true),
        prerun_iterations_update(1000),
        prerun_iterations_min(prerun_iterations_update),
        prerun_iterations_max(1e6),
        proposal("MultivariateGaussian"),
        student_t_degrees_of_freedom(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max(), 1.0),
        store_prerun(false),
        prerun_chains_per_partition(2),
        adapt_iterations(0),
        chunks(100),
        chunk_size(1000),
        need_main_run(true),
        skip_initial(0, 1, 0.1),
        store(true),
        store_observables_and_proposals(false)
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
               << "nchains = " << c.prerun_chains_per_partition
               << ", seed = " << c.seed
               << ", parallelize = " << c.parallelize
               << ", prerun min iterations = " << c.prerun_iterations_min << std::endl
               << ", prerun max iterations = " << c.prerun_iterations_max
               << ", prerun update iterations = " << c.prerun_iterations_update
               << ", skip initial = " << c.skip_initial;
        return stream;
    }
}
