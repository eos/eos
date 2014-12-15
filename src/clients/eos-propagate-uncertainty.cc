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
#include <config.h>

#include <eos/utils/destringify.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/log_prior.hh>
#include <eos/observable.hh>
#include <eos/utils/prior_sampler.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/verify.hh>

#ifdef EOS_ENABLE_PMC
#  include <eos/utils/population_monte_carlo_sampler.hh>
#endif

#include <iomanip>
#include <iostream>
#include <limits>

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
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        PriorSampler::Config config;

        ObservableSet unique_observables;

        std::vector<ObservableInput> inputs;

        Parameters parameters;

        Options global_options;

        std::string pmc_sample_file;
        unsigned pmc_sample_min, pmc_sample_max;

        std::string pmc_sample_directory;

        std::vector<LogPriorPtr> priors;

        CommandLine() :
            config(PriorSampler::Config::Default()),
            parameters(Parameters::Defaults()),
            pmc_sample_min(0),
            pmc_sample_max(0),
            pmc_sample_directory("/data")
        {
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_log_level(ll_informational);
            Log::instance()->set_program_name("eos-propagate-uncertainty");

            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

                if ("--samples" == argument)
                {
                    config.n_samples = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--workers" == argument)
                {
                    config.n_workers = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--debug" == argument)
                {
                    Log::instance()->set_log_level(ll_debug);

                    continue;
                }

                if ("--fix" == argument)
                {
                    std::string par_name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    parameters[par_name] = value;

                    continue;
                }

                if ("--global-option" == argument)
                {
                    std::string name(*(++a));
                    std::string value(*(++a));

                    global_options.set(name, value);

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

                if ("--observable" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    input.observable = Observable::make(observable_name, parameters,
                            *kinematics, global_options);
                    if (!input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    inputs.push_back(input);
                    unique_observables.add(input.observable);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--output" == argument)
                {
                    std::string filename(*(++a));

                    config.output_file.reset(new hdf5::File(hdf5::File::Create(filename)));

                    continue;
                }

                if ("--parallel" == argument)
                {
                    config.parallelize = destringify<unsigned>(*(++a));

                    continue;
                }

#if EOS_ENABLE_PMC
                if ("--pmc-sample-directory" == argument)
                {
                        pmc_sample_directory = std::string(*(++a));

                        continue;
                }

                if ("--pmc-input" == argument)
                {
                    // read samples from this file
                    pmc_sample_file = std::string(*(++a));
                    pmc_sample_min = destringify<unsigned>(*(++a));
                    pmc_sample_max = destringify<unsigned>(*(++a));

                    continue;
                }
#endif

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

                if ("--store-parameters" == argument)
                {
                    config.store_parameters = destringify<unsigned>(*(++a));

                    continue;
                }

                /*
                 * format: N_SIGMAS in [0, 10]
                 * a) --scan PAR N_SIGMAS --prior ...
                 * b) --scan PAR MIN MAX  --prior ...
                 * c) --scan PAR HARD_MIN HARD_MAX N_SIGMAS --prior ...
                 */
                if ("--vary" == argument )
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

                    priors.push_back(prior);

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

        if (inst->inputs.empty())
            throw DoUsage("No inputs specified");

        if (inst->priors.empty() and inst->pmc_sample_file.empty())
            throw DoUsage("Either specify \n a) parameters to vary\n b) a PMC input file");

        std::cout << "Determining the uncertainty on the following observables:" << std::endl;
        for (auto o = inst->inputs.begin(), o_end = inst->inputs.end() ; o != o_end ; ++o)
        {
            std::cout << o->observable->name() << "[" << o->kinematics.as_string() << "]"
                      << " with options: " << o->observable->options().as_string() << std::endl;
        }

        std::cout << std::endl;

        PriorSampler sampler(inst->unique_observables, inst->config);

#if EOS_ENABLE_PMC
        // read in parameter samples from the file and calculate observables for them
        if ( ! inst->pmc_sample_file.empty() && inst->pmc_sample_min < inst->pmc_sample_max)
        {
            auto f = hdf5::File::Open(inst->pmc_sample_file);
            auto descriptions = Analysis::read_descriptions(f);
            std::vector<std::vector<double>> samples;
            samples.push_back(std::vector<double>(descriptions.size()));
            PopulationMonteCarloSampler::read_samples(inst->pmc_sample_file, inst->pmc_sample_directory, inst->pmc_sample_min, inst->pmc_sample_max, samples);

            if (! inst->priors.empty())
            {
                std::cout << "Varying the following parameters:" << std::endl;
            }

            for (auto i = inst->priors.begin(), i_end = inst->priors.end() ; i != i_end ; ++i)
            {
                sampler.add(*i);

                std::cout << (**i).as_string() << std::endl;
            }

            sampler.run(samples, descriptions);

            return EXIT_SUCCESS;
        }
#endif
        // default: draw from priors
        {
            std::cout << "Varying the following parameters:" << std::endl;

            for (auto i = inst->priors.begin(), i_end = inst->priors.end() ; i != i_end ; ++i)
            {
                sampler.add(*i);

                std::cout << (**i).as_string() << std::endl;
            }

            sampler.run();
        }
    }
    catch (DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-propagate-uncertainty" << std::endl;
        std::cout << "  [ [--kinematics NAME VALUE]* --observable]+" << std::endl;
        std::cout << "  [--vary PARAMETER MIN MAX --prior [flat | [gaussian LOWER CENTRAL UPPER] ] ]+" << std::endl;
        std::cout << "  [--chunks VALUE]" << std::endl;
        std::cout << "  [--chunk-size VALUE]" << std::endl;
        std::cout << "  [--fix PARAMETER VALUE]" << std::endl;
        std::cout << "  [--output FILENAME]" << std::endl;
        std::cout << "  [--parallel [0|1]]" << std::endl;
#if EOS_ENABLE_PMC
        std::cout << "  [--pmc-sample-directory DIRECTORY]" << std::endl;
        std::cout << "  [--pmc-input FILENAME MIN_INDEX MAX_INDEX]" << std::endl;
#endif
        std::cout << "  [--seed LONG_VALUE]" << std::endl;
        std::cout << "  [--store-parameters]" << std::endl;
        std::cout << std::endl;
        std::cout << "Vary (nuisance) parameters in consistent way to estimate the uncertainty" << std::endl;
        std::cout << "on theory prediction of observables. Parameter samples are drawn from" << std::endl;
        std::cout << "prior distributions and the observables are calculated and stored to disk." << std::endl;
        std::cout << "One thread is created for each chunk." << std::endl;
        std::cout << "Optionally, the drawn parameters are stored as well." << std::endl;
#if EOS_ENABLE_PMC
        std::cout << std::endl;
        std::cout << "PMC options:" << std::endl;
        std::cout << "If an input file is specified, a slice of the samples is taken from there, and no new samples are drawn." << std::endl;
        std::cout << "Add a sample directory to extract samples from there within the hdf5 file. Else the default is to look for 'samples' in '/data'" << std::endl;
#endif

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
