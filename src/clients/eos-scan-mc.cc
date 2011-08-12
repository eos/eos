/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <src/constraint.hh>
#include <src/observable.hh>
#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/log.hh>
#include <src/utils/stringify.hh>
#include <src/utils/markov_chain_sampler.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <time.h>

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

        AnalysisPtr analysis;

        MarkovChainSampler::Config config;

        std::vector<ParameterData> scan_parameters;

        std::vector<ParameterData> nuisance_parameters;

        std::vector<ObservableInput> inputs;

        std::vector<Constraint> constraints;

        std::string creator;

        std::string resume_file;

        bool optimize;
        std::vector<double> starting_point;

        bool goodness_of_fit;
        std::vector<double> best_fit_point;

        CommandLine() :
            parameters(Parameters::Defaults()),
            likelihood(parameters),
            analysis(new Analysis(likelihood)),
            config(MarkovChainSampler::Config::Quick()),
            optimize(false)
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

                if (("--scan" == argument) || ("--nuisance" == argument))
                {
                    std::string name = std::string(*(++a));
                    double min = destringify<double> (*(++a));
                    double max = destringify<double> (*(++a));
                    std::string keyword = std::string(*(++a));

                    if ("--prior" != keyword)
                        throw DoUsage("Missing correct prior specification for '" + name + "'!");

                    std::string prior_type = std::string(*(++a));

                    LogPriorPtr prior;

                    ParameterRange range{ min, max };

                    if (prior_type == "gaussian")
                    {
                        double lower = destringify<double> (*(++a));
                        double central = destringify<double> (*(++a));
                        double upper = destringify<double> (*(++a));
                        prior = LogPrior::Gauss(parameters, name, range, lower, central, upper);
                    }
                    else if (prior_type == "flat")
                    {
                        prior = LogPrior::Flat(parameters, name, range);
                    }
                    else
                    {
                        throw DoUsage("Unknown prior distribution: " + prior_type);
                    }

                    bool nuisance = ("--nuisance" == argument) ? true : false;

                    if (nuisance)
                    {
                        nuisance_parameters.push_back(ParameterData{ parameters[name], min, max, prior_type });
                    }
                    else
                    {
                        scan_parameters.push_back(ParameterData{ parameters[name], min, max, prior_type });
                    }

                    // check for error in setting the prior and adding the parameter
                    if (! analysis->add(prior, nuisance))
                        throw DoUsage("Unknown error in assigning " + prior_type + " prior distribution to " + name);

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
                    if (! analysis->add(prior, true))
                        throw DoUsage("Unknown error in assigning discrete prior distribution to " + name);

                    continue;
                }

                if ("--fix" == argument)
                {
                    std::string par_name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    analysis->parameters()[par_name]=value;

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
                    if (! constraints.empty())
                        throw DoUsage("All global options must be specified before the first --observable/--constraint");

                    std::string name(*(++a));
                    std::string value(*(++a));

                    global_options.set(name, value);

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

                if ("--constraint" == argument)
                {
                    std::string constraint_name(*(++a));

                    Constraint c(Constraint::make(constraint_name, global_options));
                    likelihood.add(c);
                    constraints.push_back(c);

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

                if ("--no-prerun" == argument)
                {
                    config.need_prerun = false;

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
                    config.output_file.reset(new ScanFile(ScanFile::Create(filename, creator)));

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
                    }
                    else
                    {
                        config.seed = destringify<unsigned long>(value);
                    }

                    continue;
                }

                if ("--chunk-size" == argument)
                {
                    config.chunk_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--chunks" == argument)
                {
                    config.chunks = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--scale" == argument)
                {
                    config.scale_initial = destringify<double>(*(++a));

                    continue;
                }

                if ("--update" == argument)
                {
                    config.prerun_iterations_update = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--debug" == argument)
                {
                    Log::instance()->set_log_level(ll_debug);

                    // report exact call
                    for (int i = 0 ; i < argc ; i++)
                    {
                        std::cout << '"' << argv[i] << '"' << " ";
                    }

                    std::cout << std::endl;

                    continue;
                }

                if ("--chains" == argument)
                {
                    config.number_of_chains = destringify<unsigned>(*(++a));
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
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->inputs.empty() && (CommandLine::instance()->constraints.empty()))
            throw DoUsage("No inputs or constraints specified");

        std::cout << std::scientific;
        std::cout << "# Scan generated by eos-scan-mc" << std::endl;
        std::cout << "# Scan parameters:" << std::endl;
        for (auto p = CommandLine::instance()->scan_parameters.cbegin(), p_end = CommandLine::instance()->scan_parameters.cend() ;
                p != p_end ; ++p)
        {
            std::cout << "#   " << p->parameter.name() << " with " << p->prior << " prior "
                         " on range [" << p->min << "," << p->max << "]" << std::endl;
        }

        std::cout << "# Nuisance parameters:" << std::endl;
        for (auto p = CommandLine::instance()->nuisance_parameters.cbegin(), p_end = CommandLine::instance()->nuisance_parameters.cend() ;
                p != p_end ; ++p)
        {
            std::cout << "#   " << p->parameter.name() << " with " << p->prior << " prior "
                         " on range [" << p->min << "," << p->max << "]" << std::endl;
        }

        std::cout << "# Manual inputs:" << std::endl;
        for (auto i = CommandLine::instance()->inputs.cbegin(), i_end = CommandLine::instance()->inputs.cend() ; i != i_end ; ++i)
        {
            std::cout << "#   " << i->observable->name() << '['
                    << i->kinematics.as_string() << "] = (" << i->min << ", "
                    << i->central << ", " << i->max << ')' << std::endl;
        }

        std::cout << "# Constraints:" << std::endl;
        for (auto c = CommandLine::instance()->constraints.cbegin(), c_end = CommandLine::instance()->constraints.cend() ; c != c_end ; ++c)
        {
            std::cout << "#  " << c->name() << std::endl;
        }

        // run optimization. Use starting point if given, else sample a point from the prior.
        // Optionally calculate a p-value at the mode.
        if (CommandLine::instance()->optimize)
        {
            AnalysisPtr ana = CommandLine::instance()->analysis;
            if (CommandLine::instance()->starting_point.empty())
            {
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, ::time(0));
                for (auto i = ana->parameter_descriptions().begin(), i_end = ana->parameter_descriptions().end() ; i != i_end ; ++i)
                {
                    LogPriorPtr prior = ana->log_prior(i->parameter.name());
                    CommandLine::instance()->starting_point.push_back(prior->sample(rng));
                }
            }
            auto options = Analysis::OptimizationOptions::Defaults();
            auto ret = ana->optimize(CommandLine::instance()->starting_point, options);
            if (CommandLine::instance()->goodness_of_fit && CommandLine::instance()->best_fit_point.empty())
                ana->goodness_of_fit(ret.first, 1e5);

	        return EXIT_SUCCESS;
        }

        // goodness-of-fit for user specified parameter point
        if (CommandLine::instance()->goodness_of_fit)
        {
            CommandLine::instance()->analysis->goodness_of_fit(CommandLine::instance()->best_fit_point, 1e5);

            return EXIT_SUCCESS;
        }

        MarkovChainSampler sampler(*CommandLine::instance()->analysis, CommandLine::instance()->config);

        // extract scales and starting points from completed prerun
        if (CommandLine::instance()->resume_file != "")
        {
            std::shared_ptr<ScanFile> f(new ScanFile(ScanFile::Open(CommandLine::instance()->resume_file)));
            sampler.resume(f);
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
