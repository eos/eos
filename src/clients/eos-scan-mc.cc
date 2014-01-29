/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013 Danny van Dyk
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

#include <eos/constraint.hh>
#include <eos/observable.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/markov_chain_sampler.hh>
#include <eos/utils/population_monte_carlo_sampler.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#include <time.h>

#include <config.h>

using namespace eos;

class DoUsage
{
    private:
        std::string _what;

    public:
        DoUsage(const std::string & what) :
            _what(what)
        {
        }

        const std::string & what() const
        {
            return _what;
        }
};

struct ObservableInput
{
        ObservablePtr observable;

        Kinematics kinematics;

        double min, central, max;
};

struct ParameterData
{
        Parameter parameter;

        double min;

        double max;

        std::string prior;
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        Options global_options;

        LogLikelihood likelihood;

        Analysis analysis;

        MarkovChainSampler::Config config;
#if EOS_ENABLE_PMC
        PopulationMonteCarloSampler::Config config_pmc;
#endif

        std::string global_local_output;
        proposal_functions::GlobalLocal::Config config_gl;
        bool use_global_local;
        std::vector<std::shared_ptr<hdf5::File>> prerun_inputs;

        std::vector<ParameterData> scan_parameters;

        std::vector<ParameterData> nuisance_parameters;

        std::vector<ObservableInput> inputs;

        std::vector<Constraint> constraints;

        std::string creator;

        std::shared_ptr<unsigned> partition_index;

        std::string resume_file;
        std::string pmc_initialization_file;
        std::string pmc_sample_file;
        bool pmc_calculate_posterior;
        unsigned pmc_calculate_posterior_min, pmc_calculate_posterior_max;
        bool pmc_draw_samples;
        bool pmc_final;
        bool pmc_update;

        // use MINUIT
        bool massive_mode_finding;
        unsigned massive_maximum_iterations;

        bool optimize;
        std::vector<double> starting_point;

        bool goodness_of_fit;
        std::vector<double> best_fit_point;

        bool use_pmc;

        CommandLine() :
            parameters(Parameters::Defaults()),
            likelihood(parameters),
            analysis(likelihood),
            config(MarkovChainSampler::Config::Quick()),
#if EOS_ENABLE_PMC
            config_pmc(PopulationMonteCarloSampler::Config::Default()),
#endif
            config_gl(proposal_functions::GlobalLocal::Config::Default()),
            use_global_local(false),
            pmc_calculate_posterior(false),
            pmc_calculate_posterior_min(0),
            pmc_calculate_posterior_max(0),
            pmc_draw_samples(false),
            pmc_final(false),
            massive_mode_finding(false),
            massive_maximum_iterations(2000),
            optimize(false),
            use_pmc(false)
        {
            config.number_of_chains = 4;
            config.need_prerun = true;
            config.chunk_size = 1000;
            config.parallelize = true;
            config.use_strict_rvalue_definition = true;
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_log_level(ll_informational);
            Log::instance()->set_program_name("eos-scan-mc");

            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            creator = std::string(argv[0]);
            for (int i = 1 ; i < argc ; ++i)
            {
                creator += ' ' + std::string(argv[i]);
            }

            for (char ** a(argv + 1), **a_end(argv + argc); a != a_end; ++a)
            {
                std::string argument(*a);

                /*
                 * format: N_SIGMAS in [0, 10]
                 * a) --scan PAR N_SIGMAS --prior ...
                 * b) --scan PAR MIN MAX  --prior ...
                 * c) --scan PAR HARD_MIN HARD_MAX N_SIGMAS --prior ...
                 */
                if (("--scan" == argument) || ("--nuisance" == argument))
                {
                    std::string name = std::string(*(++a));

                    double min = -std::numeric_limits<double>::max();
                    double max =  std::numeric_limits<double>::max();

                    // first word has to be a number
                    double number = destringify<double>(*(++a));

                    std::string keyword = std::string(*(++a));

                    VerifiedRange<double> n_sigmas(0, 10, 0);

                    // case a)
                    if ("--prior" == keyword)
                    {
                        n_sigmas = VerifiedRange<double>(0, 10, number);
                        if (n_sigmas == 0)
                            throw DoUsage("number of sigmas: number expected");
                    }
                    else
                    {
                        // case b), c)
                        min = number;
                        max = destringify<double>(keyword);

                        keyword = std::string(*(++a));

                        // watch for case c)
                        if ("--prior" != keyword)
                        {
                            n_sigmas = VerifiedRange<double>(0, 10,  destringify<double>(keyword));
                            if (n_sigmas == 0)
                                throw DoUsage("number of sigmas: number expected");
                            keyword = std::string(*(++a));
                        }
                    }

                    if ("--prior" != keyword)
                        throw DoUsage("Missing correct prior specification for '" + name + "'!");

                    std::string prior_type = std::string(*(++a));

                    LogPriorPtr prior;

                    ParameterRange range{ min, max };

                    if (prior_type == "gaussian" || prior_type == "log-gamma")
                    {
                        double lower = destringify<double> (*(++a));
                        double central = destringify<double> (*(++a));
                        double upper = destringify<double> (*(++a));

                        // adjust range, but always stay within hard bound supplied by the user
                        if (n_sigmas > 0)
                        {
                            range.min = std::max(range.min, central - n_sigmas * (central - lower));
                            range.max = std::min(range.max, central + n_sigmas * (upper - central));
                        }
                        if (prior_type == "gaussian")
                        {
                            prior = LogPrior::Gauss(parameters, name, range, lower, central, upper);
                        }
                        else
                        {
                            prior = LogPrior::LogGamma(parameters, name, range, lower, central, upper);
                        }
                    }
                    else if (prior_type == "flat")
                    {
                        if (n_sigmas > 0)
                            throw DoUsage("Can't specify number of sigmas for flat prior");
                        prior = LogPrior::Flat(parameters, name, range);
                    }
                    else
                    {
                        throw DoUsage("Unknown prior distribution: " + prior_type);
                    }

                    bool nuisance = ("--nuisance" == argument) ? true : false;

                    if (nuisance)
                    {
                        nuisance_parameters.push_back(ParameterData{ parameters[name], range.min, range.max, prior_type });
                    }
                    else
                    {
                        scan_parameters.push_back(ParameterData{ parameters[name], range.min, range.max, prior_type });
                    }

                    // check for error in setting the prior and adding the parameter
                    if (! analysis.add(prior, nuisance))
                        throw DoUsage("Error in assigning " + prior_type + " prior distribution to '" + name +
                                      "'. Perhaps '" + name + "' appears twice in the list of parameters?");

                    continue;
                }

                if ("--build-global-local" == argument)
                {
                    global_local_output = (*(++a));
                    config.global_local_config.reset(new proposal_functions::GlobalLocal::Config(
                                                         proposal_functions::GlobalLocal::Config::Default()));

                    continue;
                }

                if ("--chains" == argument)
                {
                    config.number_of_chains = destringify<unsigned>(*(++a));
                    continue;
                }

                if ("--chunk-size" == argument)
                {
                    config.chunk_size = destringify<unsigned>(*(++a));
#if EOS_ENABLE_PMC
                    config_pmc.chunk_size = config.chunk_size;
#endif

                    continue;
                }

                if ("--chunks" == argument)
                {
                    config.chunks = destringify<unsigned>(*(++a));
#if EOS_ENABLE_PMC
                    config_pmc.chunks = config.chunks;
#endif

                    continue;
                }

                if ("--constraint" == argument)
                {
                    std::string constraint_name(*(++a));

                    Constraint c(Constraint::make(constraint_name, global_options));
                    likelihood.add(c);
                    constraints.push_back(c);

                    continue;
                }

                if ("--debug" == argument)
                {
                    Log::instance()->set_log_level(ll_debug);

                    continue;
                }

                if ("--discrete" == argument)
                {
                    std::string name = std::string(*(++a));

                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                        throw DoUsage("Put set of discrete values in braces {}");

                    std::set<double> values;
                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        double value = destringify<double> (word);
                        values.insert(value);
                    }
                    while (true);

                    LogPriorPtr prior = LogPrior::Discrete(parameters, name, values);

                    // check for error in setting the prior and adding the parameter
                    if (! analysis.add(prior, true))
                        throw DoUsage("Unknown error in assigning discrete prior distribution to " + name);

                    continue;
                }

                if ("--fix" == argument)
                {
                    std::string par_name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    analysis.parameters()[par_name]=value;

                    continue;
                }

                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    kinematics->declare(name);
                    kinematics->set(name, value);

                    continue;
                }

                if ("--global-option" == argument)
                {
                    std::string name(*(++a));
                    std::string value(*(++a));

                    if (! constraints.empty())
                    {
                        Log::instance()->message("eos-scan-mc", ll_warning)
                            << "Global option (" << name << " = " << value <<") only applies to observables/constraints defined from now on, "
                            << "but doesn't affect the " << constraints.size() << " previously defined constraints.";
                    }

                    global_options.set(name, value);

                    continue;
                }

                if ("--global-local" == argument)
                {
                    use_global_local = destringify<unsigned>(*(++a));
                    if (use_global_local)
                        config.store_prerun = true;

                    continue;
                }

                if ("--global-local-adapt-iterations" == argument)
                {
                    config.adapt_iterations = destringify<unsigned> (*(++a));

                    continue;
                }

                if ("--global-local-covariance-window" == argument)
                {
                    config_gl.history_points_local_covariance_size = destringify<unsigned> (*(++a));
#if EOS_ENABLE_PMC
                    config_pmc.sliding_window = config_gl.history_points_local_covariance_size;
#endif

                    continue;
                }

                if ("--global-local-equal-weights" == argument)
                {
                    config_gl.equal_weight_components = destringify<unsigned> (*(++a));

                    continue;
                }

                if ("--global-local-history-points" == argument)
                {
                    config_gl.history_points = destringify<unsigned> (*(++a));

                    continue;
                }

                if ("--global-local-history-points-ordered" == argument)
                {
                    config_gl.history_points_ordered = destringify<unsigned> (*(++a));

                    continue;
                }

                if ("--global-local-input" == argument)
                {
                    // parse all filename between { and }, separated by whitespace
                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                    {
                        --a;
                        throw DoUsage("Need to specify prerun inputs like --global-local-input { file1 file2 } ");
                        continue;
                    }

                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        prerun_inputs.push_back(std::make_shared<hdf5::File>(hdf5::File::Open(word, H5F_ACC_RDONLY)));
                    }
                    while (true);

                    continue;
                }

                if ("--global-local-jump-indices" == argument)
                {
                    // parse all filename between { and }, separated by whitespace
                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                    {
                        --a;
                        throw DoUsage("Need to specify indices like --global-local-jump-indices { index0 index1 } ");
                        continue;
                    }

                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        config_gl.long_jump_indices.push_back(destringify<unsigned>(word));
                    }
                    while (true);

                    continue;
                }

                if ("--global-local-jump-probability" == argument)
                {
                    config_gl.local_jump_probability = destringify<double> (*(++a));

                    continue;
                }

                if ("--global-local-minimum-weight" == argument)
                {
                    config_gl.minimum_relative_cluster_weight = destringify<double> (*(++a));

                    continue;
                }

                if ("--global-local-perform-clustering" == argument)
                {
                    config_gl.perform_clustering = true;

                    config_gl.clustering_maximum_r_value = destringify<double> (*(++a));

                    continue;
                }

                if ("--global-local-points" == argument)
                {
                    config_gl.history_points = destringify<unsigned> (*(++a));

                    continue;
                }

                if ("--global-local-point-weighting" == argument)
                {
                    typedef proposal_functions::GlobalLocal::Config::HistoryPointWeighting hpw;

                    std::string type(*(++a));
                    if ("log" == type)
                    {
                        config_gl.history_point_weighting = hpw::log_posterior;
                    }
                    else if ("linear" == type)
                    {
                        config_gl.history_point_weighting = hpw::posterior;
                    }
                    else if("flat" == type)
                    {
                        config_gl.history_point_weighting = hpw::equal;
                    }
                    else
                    {
                        throw DoUsage("global-local-point-weighting: invalid type " + type);
                    }

                    continue;
                }

                if ("--global-local-rescale-local-covariance" == argument)
                {
                    config_gl.rescale_local_covariance = destringify<double> (*(++a));

                    continue;
                }

                if ("--global-local-skip-initial" == argument)
                {
                    config_gl.skip_initial = destringify<double> (*(++a));
#if EOS_ENABLE_PMC
                    config_pmc.skip_initial = config_gl.skip_initial;
#endif

                    continue;
                }

                if ("--global-local-strict-clustering" == argument)
                {
                    config_gl.clustering_strict_r_value = destringify<unsigned> (*(++a));

                    continue;
                }

                if ("--goodness-of-fit" == argument)
                {
                    // best-fit point is optional
                    goodness_of_fit = true;

                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                    {
                        --a;
                        continue;
                    }

                    // parse starting point
                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        double value = destringify<double>(word);
                        best_fit_point.push_back(value);
                    }
                    while (true);

                    continue;
                }

                if ("--massive-mode-finding" == argument)
                {
                    massive_mode_finding = true;
                    massive_maximum_iterations = destringify<unsigned> (*(++a));
                    if (massive_maximum_iterations == 0)
                    {
                        throw DoUsage("Need to specify maximum number of Minuit iterations for massive mode finding");
                    }

                    continue;
                }

                if ("--no-prerun" == argument)
                {
                    config.need_prerun = false;

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    input.observable = Observable::make(observable_name, parameters,
                            *kinematics, global_options);
                    if (!input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double> (*(++a));
                    input.central = destringify<double> (*(++a));
                    input.max = destringify<double> (*(++a));

                    likelihood.add(input.observable, input.min, input.central, input.max);

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--observable-prior" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    input.observable = Observable::make(observable_name, parameters,
                            *kinematics, global_options);
                    if (!input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double> (*(++a));
                    input.central = destringify<double> (*(++a));
                    input.max = destringify<double> (*(++a));

                    // TODO: Hack! This is only used for putting parts of the
                    // prior into the likelihood for correlated prior information.
                    likelihood.add(input.observable, input.min, input.central, input.max, 0);

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--optimize" == argument)
                {
                    optimize = true;

                    // starting point is optional
                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                    {
                        --a;
                        continue;
                    }

                    // parse starting point
                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        double value = destringify<double>(word);
                        starting_point.push_back(value);
                    }
                    while (true);

                    continue;
                }

                if ("--output" == argument)
                {
                    std::string filename(*(++a));
                    config.output_file = filename;
#if EOS_ENABLE_PMC
                    config_pmc.output_file = filename;
#endif

                    continue;
                }

                if ("--parallel" == argument)
                {
                    config.parallelize = destringify<unsigned>(*(++a));
#if EOS_ENABLE_PMC
                    config_pmc.parallelize = config.parallelize;
#endif

                    continue;
                }

                if ("--partition" == argument)
                {
                    std::string key(*(++a));
                    std::vector<std::tuple<std::string, double, double>> partition;
                    while (key.substr(0, 2) != "--")
                    {
                        std::string name = key;
                        double min = destringify<double>(*(++a));
                        double max = destringify<double>(*(++a));
                        partition.push_back(std::make_tuple(name, min, max));

                        key = std::string(*(++a));
                    }
                    config.partitions.push_back(partition);
                    --a;

                    continue;
                }

                if ("--partition-index" == argument)
                {
                    partition_index.reset(new unsigned);
                    *partition_index = destringify<unsigned> (*(++a));

                    config.need_main_run = false;
                    config.store_prerun = true;

                    continue;
                }

                if ("--use-pmc" == argument)
                {
                    use_pmc = true;

                    continue;
                }

#if EOS_ENABLE_PMC
                if ("--pmc-adjust-sample-size" == argument)
                {
                    config_pmc.adjust_sample_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-block-decomposition" == argument)
                {
                    config_pmc.block_decomposition = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-crop-highest-weights" == argument)
                {
                    config_pmc.crop_highest_weights = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-dof" == argument)
                {
                    config_pmc.degrees_of_freedom = destringify<long>(*(++a));
                    config_pmc.override_global_local_proposal = true;

                    continue;
                }

                if ("--pmc-calculate-posterior" == argument)
                {
                    pmc_calculate_posterior = true;

                    // read samples from this file
                    pmc_sample_file = std::string(*(++a));
                    pmc_calculate_posterior_min = destringify<unsigned>(*(++a));
                    pmc_calculate_posterior_max = destringify<unsigned>(*(++a));

                    // read components from this file
                    pmc_initialization_file = pmc_sample_file;

                    continue;
                }

                if ("--pmc-components-per-cluster" == argument)
                {
                    unsigned n_comp = destringify<unsigned>(*(++a));

                    // with minuit initialization, total number of components depends on local modes found
                    config_pmc.components_per_cluster = n_comp;

                    continue;
                }

                if ("--pmc-draw-samples" == argument)
                {
                    // samples are to be stored in ordinary output file via config
                    pmc_draw_samples = true;

                    continue;
                }

                if ("--pmc-final" == argument)
                {
                    pmc_final = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-ignore-group" == argument)
                {
                    config_pmc.ignore_groups.push_back(destringify<unsigned>(*(++a)));

                    continue;
                }

                if ("--pmc-initialize-from-file" == argument)
                {
                    pmc_initialization_file = std::string(*(++a));

                    continue;
                }

                if ("--pmc-hierarchical-clusters" == argument)
                {
                    config_pmc.super_clusters = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-minimum-overlap" == argument)
                {
                    config_pmc.minimum_overlap = destringify<double>(*(++a));

                    continue;
                }

                if ("--pmc-patch-around-local-mode" == argument)
                {
                    config_pmc.patch_around_local_mode = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-group-by-r-value" == argument)
                {
                    config_pmc.group_by_r_value = destringify<double>(*(++a));

                    continue;
                }

                if ("--pmc-r-value-no-nuisance" == argument)
                {
                    config_pmc.r_value_no_nuisance = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-single-cluster" == argument)
                {
                    config_pmc.single_cluster = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-starting-points" == argument)
                {
                    config_pmc.starting_points = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-final-chunksize" == argument)
                {
                    config_pmc.final_chunk_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-ignore-ess" == argument)
                {
                    config_pmc.ignore_eff_sample_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-relative-std-deviation-over-last-steps" == argument)
                {
                    config_pmc.maximum_relative_std_deviation = destringify<double>(*(++a));
                    config_pmc.minimum_steps = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-std-dev-reduction" == argument)
                {
                    config_pmc.std_dev_reduction = destringify<double>(*(++a));

                    continue;
                }

                if ("--pmc-update" == argument)
                {
                    pmc_update = true;
                    pmc_initialization_file = std::string(*(++a));

                    continue;
                }
#endif

                if ("--prerun-chains-per-partition" == argument)
                {
                    config.prerun_chains_per_partition = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-find-modes" == argument)
                {
                    config.find_modes = true;

                    continue;
                }

                if ("--prerun-max" == argument)
                {
                    config.prerun_iterations_max = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-min" == argument)
                {
                    config.prerun_iterations_min = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-only" == argument)
                {
                    config.need_prerun = true;
                    config.store_prerun = true;
                    config.need_main_run = false;

                    continue;
                }

                if ("--prerun-update" == argument)
                {
                    config.prerun_iterations_update = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--print-args" == argument)
                {
                    // print arguments and quit
                   for (int i = 1 ; i < argc ; i++)
                    {
                       std::cout << "'" << argv[i] << "' ";
                    }

                    std::cout << std::endl;
                    abort();

                    continue;
                }

                if ("--prior-as-proposal" == argument)
                {
                    // [parameter_name]
                    std::string name = *(++a);
                    LogPriorPtr proposal = analysis.log_prior(name);
                    if (! proposal)
                        throw DoUsage("Define parameter " + name + " and its prior before --prior-as-proposal");
                    config.block_proposal_parameters.push_back(name);

                    continue;
                }

                if ("--proposal" == argument)
                {
                    config.proposal = *(++a);

                    if (config.proposal == "MultivariateStudentT")
                    {
                        double dof = destringify<double>(*(++a));
                        if (dof <= 0)
                        {
                            throw DoUsage("No (or non-positive) degree of freedom for MultivariateStudentT specified");
                        }
                        config.student_t_degrees_of_freedom = dof;
                    }

                    continue;
                }

                if ("--resume" == argument)
                {
                    resume_file = std::string(*(++a));
                    config.need_prerun = false;

                    continue;
                }

                if ("--seed" == argument)
                {
                    std::string value(*(++a));

                    if ("time" == value)
                    {
                        config.seed = ::time(0);
#if EOS_ENABLE_PMC
                        config_pmc.seed = ::time(0);
#endif
                    }
                    else
                    {
                        config.seed = destringify<unsigned long>(value);
#if EOS_ENABLE_PMC
                        config_pmc.seed = destringify<unsigned long>(value);
#endif
                    }

                    continue;
                }

                if ("--scale-reduction" == argument)
                {
                    config.scale_reduction = destringify<double>(*(++a));

                    continue;
                }

                if ("--store-prerun" == argument)
                {
                    config.store_prerun = true;

                    continue;
                }

                if ("--store-observables-and-proposals" == argument)
                {
                    config.store_observables_and_proposals = true;

                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

int main(int argc, char * argv[])
{
    try
    {
        auto inst = CommandLine::instance();
        inst->parse(argc, argv);

        if (inst->inputs.empty() &&
            inst->constraints.empty() &&
            inst->global_local_output.empty())
            throw DoUsage("No inputs, constraints nor build output specified");

        if (inst->nuisance_parameters.empty() &&
            inst->scan_parameters.empty())
           throw  DoUsage("No scan or nuisance parameters defined");

        std::cout << std::scientific;
        std::cout << "# Scan generated by eos-scan-mc" << std::endl;
        if ( ! inst->scan_parameters.empty())
        {
            std::cout << "# Scan parameters (" << inst->scan_parameters.size() << "):" << std::endl;
            for (auto d = inst->analysis.parameter_descriptions().cbegin(), d_end = inst->analysis.parameter_descriptions().cend() ;
                 d != d_end ; ++d)
            {
                if (d->nuisance)
                    continue;
                std::cout << "#   " << inst->analysis.log_prior(d->parameter.name())->as_string() << std::endl;
            }
        }

        if ( ! inst->nuisance_parameters.empty())
        {
            std::cout << "# Nuisance parameters (" << inst->nuisance_parameters.size() << "):" << std::endl;
            for (auto d = inst->analysis.parameter_descriptions().cbegin(), d_end = inst->analysis.parameter_descriptions().cend() ;
                 d != d_end ; ++d)
            {
                if ( ! d->nuisance)
                    continue;
                std::cout << "#   " << inst->analysis.log_prior(d->parameter.name())->as_string() << std::endl;
            }
        }

        if ( ! inst->inputs.empty())
        {
            std::cout << "# Manual inputs (" << inst->inputs.size() << "):" << std::endl;
            for (auto i = inst->inputs.cbegin(), i_end = inst->inputs.cend() ; i != i_end ; ++i)
            {
                std::cout << "#   " << i->observable->name() << '['
                    << i->kinematics.as_string() << "] = (" << i->min << ", "
                    << i->central << ", " << i->max << ')' << std::endl;
            }
        }

        if ( ! inst->constraints.empty())
        {
            std::cout << "# Constraints (" << inst->constraints.size() << "):" << std::endl;
            for (auto c = inst->constraints.cbegin(), c_end = inst->constraints.cend() ; c != c_end ; ++c)
            {
                std::cout << "#  " << c->name() << ": ";
                for (auto o = c->begin_observables(), o_end = c->end_observables(); o != o_end ; ++o)
                {
                    std::cout << (**o).name() << '['
                        << (**o).kinematics().as_string() << ']'
                        << " with options: " << (**o).options().as_string();
                }
                for (auto b = c->begin_blocks(), b_end = c->end_blocks(); b != b_end ; ++b)
                {
                    std::cout << ", " << (**b).as_string();
                }
                std::cout << std::endl;
            }
        }

        // run optimization. Use starting point if given, else sample a point from the prior.
        // Optionally calculate a p-value at the mode.
        if (inst->optimize)
        {
            Analysis & ana(inst->analysis);
            if (inst->starting_point.empty())
            {
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, ::time(0));
                for (auto i = ana.parameter_descriptions().begin(), i_end = ana.parameter_descriptions().end() ; i != i_end ; ++i)
                {
                    LogPriorPtr prior = ana.log_prior(i->parameter.name());
                    inst->starting_point.push_back(prior->sample(rng));
                }
            }

            if (inst->starting_point.size() != ana.parameter_descriptions().size())
            {
                throw DoUsage("Starting point size of" + stringify(inst->starting_point.size())
                              + " doesn't match with analysis size of " + stringify(ana.parameter_descriptions().size()));
            }

            std::cout << std::endl;
            std::cout << "# Starting optimization at " << stringify_container(inst->starting_point, 4) << std::endl;
            std::cout << std::endl;

            auto options = Analysis::OptimizationOptions::Defaults();
            auto ret = ana.optimize_minuit(inst->starting_point, options);

            Log::instance()->message("eos-scan-mc", ll_informational)
                << "Result from minuit:" << ret << ret.UserCovariance();
            Log::instance()->message("eos-scan-mc", ll_informational)
                << "Best result: log(posterior) at "
                << stringify_container(ret.UserParameters().Params(), 6)
                << " = " << -1.0 * ret.Fval();

            if (inst->goodness_of_fit && inst->best_fit_point.empty())
                ana.goodness_of_fit(ret.UserParameters().Params(), 1e5);

            return EXIT_SUCCESS;
        }

        // goodness-of-fit for user specified parameter point
        if (inst->goodness_of_fit)
        {
            inst->analysis.goodness_of_fit(inst->best_fit_point, 1e5, inst->config.output_file);

            return EXIT_SUCCESS;
        }

#if EOS_ENABLE_PMC
        if (CommandLine::instance()->use_pmc)
        {
            std::shared_ptr<PopulationMonteCarloSampler> pop_sampler;
            if (inst->pmc_initialization_file != "")
                pop_sampler.reset(
                    new PopulationMonteCarloSampler(inst->analysis, hdf5::File::Open(inst->pmc_initialization_file),
                                                    inst->config_pmc, inst->pmc_update));
            else
                pop_sampler.reset(new PopulationMonteCarloSampler(inst->analysis, inst->config_pmc));

            if (inst->pmc_final)
            {
                PopulationMonteCarloSampler::Status status =  pop_sampler->status();
                status.converged = true;
                pop_sampler->status(status);
            }
            if (inst->pmc_draw_samples)
            {
                pop_sampler->draw_samples();
            }
            else if (inst->pmc_calculate_posterior)
            {
                pop_sampler->calculate_weights(inst->pmc_sample_file,
                    inst->pmc_calculate_posterior_min, inst->pmc_calculate_posterior_max);
            }
            else if (inst->pmc_update)
                return EXIT_SUCCESS;
            else
                pop_sampler->run();

            return EXIT_SUCCESS;
        }
#endif

        // remove unwanted partitions and select only one
        if (unsigned * i = inst->partition_index.get())
        {
            MarkovChainSampler::Config & c = inst->config;
            if (c.partitions.empty())
                throw DoUsage("Can't select partition " + stringify(*i) + " from no partitions!");

            auto temp = c.partitions;
            c.partitions.clear();
            c.partitions.push_back(temp.at(*i));
        }

        if (! inst->global_local_output.empty())
        {
            MarkovChainSampler::build_global_local(inst->global_local_output, inst->prerun_inputs, inst->config_gl);
            return EXIT_SUCCESS;
        }

        // create the sampler with correct options
        if (inst->use_global_local)
        {
            inst->config.global_local_config.reset(new proposal_functions::GlobalLocal::Config(inst->config_gl));
        }
        MarkovChainSampler sampler(inst->analysis, inst->config);

        if (inst->massive_mode_finding)
        {
            // try to find just anything
            Analysis::OptimizationOptions options = Analysis::OptimizationOptions::Defaults();
            options.algorithm = "minimize";
            options.maximum_iterations = inst->massive_maximum_iterations;
            options.mcmc_pre_run = inst->config.need_prerun;
            options.strategy_level = 0;
            sampler.massive_mode_finding(options);
            return EXIT_SUCCESS;
        }

        // extract scales and starting points from completed prerun
        if (inst->resume_file != "")
        {
            sampler.resume(hdf5::File::Open(inst->resume_file));
            return EXIT_SUCCESS;
        }

        sampler.run();
    }
    catch (DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan-mc" << std::endl;
        std::cout << "  [ [--kinematics NAME VALUE]* --observable NAME LOWER CENTRAL UPPER]+" << std::endl;
        std::cout << "  [--constraint NAME]+" << std::endl;
        std::cout << "  [ [ [--scan PARAMETER MIN MAX] | [--nuisance PARAMETER MIN MAX] ] --prior [flat | [gaussian LOWER CENTRAL UPPER] ] ]+" << std::endl;
        std::cout << "  [--chains VALUE]" << std::endl;
        std::cout << "  [--chunks VALUE]" << std::endl;
        std::cout << "  [--chunksize VALUE]" << std::endl;
        std::cout << "  [--debug]" << std::endl;
        std::cout << "  [--discrete PARAMETER { VALUE1 VALUE2 ...}]+" << std::endl;
        std::cout << "  [--fix PARAMETER VALUE]+" << std::endl;
        std::cout << "  [--goodness_of_fit [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--no-prerun]" << std::endl;
        std::cout << "  [--optimize [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--output FILENAME]" << std::endl;
        std::cout << "  [--resume FILENAME]" << std::endl;
        std::cout << "  [--scale VALUE]" << std::endl;
        std::cout << "  [--seed LONG_VALUE]" << std::endl;
        std::cout << "  [--store-prerun]" << std::endl;


        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-scan-mc --kinematics s_min 14.18 --kinematics s_max 16.00 \\" << std::endl;
        std::cout << "      --observable \"B->K^*ll::BR@LowRecoil\" 0.5e-7 1.25e-7 2.0e-7 \\" << std::endl;
        std::cout << "      --constraint \"B^0->K^*0gamma::BR@BaBar-2009\" \\" << std::endl;
        std::cout << "      --scan     \"Abs{c9}\"        0.0 15.0     --prior flat\\" << std::endl;
        std::cout << "      --scan     \"Arg{c9}\"        0.0  6.28319 --prior flat\\" << std::endl;
        std::cout << "      --nuisance \"mass::b(MSbar)\" 3.8  5.0     --prior gaussian 4.14 4.27 4.37" << std::endl;

        return EXIT_FAILURE;
    }
    catch (Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
