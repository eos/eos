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

#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/log.hh>
#include <src/utils/log_prior.hh>
#include <src/observable.hh>
#include <src/utils/prior_sampler.hh>
#include <src/utils/scan_file.hh>
#include <src/utils/stringify.hh>

#include <iomanip>
#include <iostream>

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

        std::vector<LogPriorPtr> priors;

        CommandLine() :
            config(PriorSampler::Config::Default()),
            parameters(Parameters::Defaults())
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

                if ("--fix" == argument)
                {
                    std::string par_name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    parameters[par_name] = value;

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

                    inputs.push_back(input);
                    unique_observables.add(input.observable);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--output" == argument)
                {
                    std::string filename(*(++a));

                    config.output_file.reset(new ScanFile(ScanFile::Create(filename, "eos-propagate-uncertainty")));

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

                if ("--store-parameters" == argument)
                {
                    config.store_parameters = true;

                    continue;
                }

                if ("--vary" == argument )
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
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->inputs.empty())
            throw DoUsage("No inputs specified");

        PriorSampler sampler(CommandLine::instance()->unique_observables, CommandLine::instance()->config);
        for (auto i = CommandLine::instance()->priors.begin(), i_end = CommandLine::instance()->priors.end() ; i != i_end ; ++i)
        {
            sampler.add(*i);
        }

        sampler.run();
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
        std::cout << "  [--seed LONG_VALUE]" << std::endl;
        std::cout << "  [--store-parameters]" << std::endl;
        std::cout << std::endl;
        std::cout << "Vary (nuisance) parameters in consistent way to estimate the uncertainty" << std::endl;
        std::cout << "on theory prediction of observables. Parameter samples are drawn from" << std::endl;
        std::cout << "prior distributions and the observables are calculated and stored to disk." << std::endl;
        std::cout << "One thread is created for each chunk." << std::endl;
        std::cout << "Optionally, the drawn parameters are stored as well." << std::endl;

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
