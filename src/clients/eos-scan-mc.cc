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

#include <src/observable.hh>
#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
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

        LogLikelihoodPtr likelihood;

        AnalysisPtr analysis;

        MarkovChainSampler::Config config;

        std::vector<ParameterData> scan_parameters;

        std::vector<ParameterData> nuisance_parameters;

        std::vector<ObservableInput> inputs;

        CommandLine() :
            parameters(Parameters::Defaults()),
            likelihood(new LogLikelihood(parameters)),
            analysis(new Analysis(likelihood)),
            config(MarkovChainSampler::Config::Quick())
        {
            config.number_of_chains = 4;
            config.need_prerun = true;
            config.chunk_size = 1000;
            config.parallelize = true;
            config.use_strict_rvalue_definition = true;
        }

        void parse(int argc, char ** argv)
        {
            std::shared_ptr<Kinematics> kinematics(new Kinematics);

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

                if ("--observable" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    input.observable = Observable::make(observable_name, parameters,
                            *kinematics, Options());
                    if (!input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double> (*(++a));
                    input.central = destringify<double> (*(++a));
                    input.max = destringify<double> (*(++a));

                    likelihood->add(input.observable, input.min, input.central, input.max);

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--output" == argument)
                {
                    std::string filename(*(++a));

                    config.output_file.reset(new ScanFile(ScanFile::Create(filename, "eos-scan-mc")));

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
                        config.seed = destringify<unsigned long>(*(++a));
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
                    // only sequential execution
                    config.parallelize = false;

                    // report exact call
                    for(int i = 0 ; i < argc ; i++)
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

        if (CommandLine::instance()->inputs.empty())
            throw DoUsage("No inputs specified");

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

        std::cout << "# Inputs:" << std::endl;
        for (auto i = CommandLine::instance()->inputs.cbegin(), i_end = CommandLine::instance()->inputs.cend() ; i != i_end ; ++i)
        {
            std::cout << "#   " << i->observable->name() << '['
                    << i->kinematics.as_string() << "] = (" << i->min << ", "
                    << i->central << ", " << i->max << ')' << std::endl;
        }

        MarkovChainSampler sampler(*CommandLine::instance()->analysis, CommandLine::instance()->config);
        sampler.run();
    }
    catch (DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan-mc" << std::endl;
        std::cout << "  [ [--kinematics NAME VALUE]* --observable NAME LOWER CENTRAL UPPER]+" << std::endl;
        std::cout << "  [ [ [--scan PARAMETER MIN MAX] | [--nuisance PARAMETER MIN MAX] ] --prior [flat | [gaussian LOWER CENTRAL UPPER] ] ]+" << std::endl;
        std::cout << "  [--fix PARAMETER VALUE]" << std::endl;
        std::cout << "  [--chunksize VALUE]+" << std::endl;
        std::cout << "  [--chunks VALUE]+" << std::endl;
        std::cout << "  [--n_chains VALUE]+" << std::endl;
        std::cout << "  [--seed LONG_VALUE]+" << std::endl;
        std::cout << "  [--scale VALUE]+" << std::endl;
        std::cout << "  [--output FILENAME]+" << std::endl;
        std::cout << "  [--debug]" << std::endl;

        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-scan-mc --kinematics s_min 14.18 --kinematics s_max 16.00 \\" << std::endl;
        std::cout << "      --observable \"B->K^*ll::BR@LowRecoil\" 0.5e-7 1.25e-7 2.0e-7 \\" << std::endl;
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
