/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Danny van Dyk
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

#include <eos/signal-pdf.hh>
#include <eos/statistics/markov-chain-sampler.hh>
#include <eos/utils/density.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>

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

struct KinematicsData
{
    KinematicVariable kinematic_variable;

    double min;

    double max;
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        DensityPtr density;

        Parameters parameters;

        Options global_options;

        Kinematics kinematics;

        std::vector<KinematicsData> kinematics_data;

        MarkovChainSampler::Config mcmc_config;

        std::string creator;

        double scale_reduction;

        CommandLine() :
            parameters(Parameters::Defaults()),
            mcmc_config(MarkovChainSampler::Config::Quick()),
            scale_reduction(1)
        {
            mcmc_config.number_of_chains = 4;
            mcmc_config.need_prerun = true;
            mcmc_config.chunk_size = 1000;
            mcmc_config.parallelize = true;
            mcmc_config.use_strict_rvalue_definition = true;
            mcmc_config.rvalue_criterion_param = 1.05;
            mcmc_config.rvalue_criterion_posterior = 1.05;
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_log_level(ll_informational);
            Log::instance()->set_program_name("eos-sample-events-mcmc");

            creator = std::string(argv[0]);
            for (int i = 1 ; i < argc ; ++i)
            {
                creator += ' ' + std::string(argv[i]);
            }

            for (char ** a(argv + 1), **a_end(argv + argc); a != a_end; ++a)
            {
                std::string argument(*a);

                /*
                 * format:
                 *    --kinematics NAME MIN MAX
                 */
                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));

                    double min = destringify<double>(*(++a));
                    double max = destringify<double>(*(++a));

                    kinematics_data.push_back(KinematicsData{ kinematics.declare(name, (max + min) / 2.0), min, max});
                }

                if ("--chains" == argument)
                {
                    mcmc_config.number_of_chains = destringify<unsigned>(*(++a));
                    continue;
                }

                if ("--chunk-size" == argument)
                {
                    mcmc_config.chunk_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--chunks" == argument)
                {
                    mcmc_config.chunks = destringify<unsigned>(*(++a));

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
                    parameters.set(par_name, value);

                    continue;
                }

                if ("--signal-pdf" == argument)
                {
                    std::string signal_pdf_name(*(++a));
                    density = SignalPDF::make(signal_pdf_name, parameters, kinematics, global_options);

                    continue;
                }

                if ("--output" == argument)
                {
                    std::string filename(*(++a));
                    mcmc_config.output_file = filename;

                    continue;
                }

                if ("--parallel" == argument)
                {
                    mcmc_config.parallelize = destringify<unsigned>(*(++a));

                    continue;
                }

                // todo rename here and in scripts
                if ("--prerun-chains-per-partition" == argument)
                {
                    mcmc_config.number_of_chains = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-max" == argument)
                {
                    mcmc_config.prerun_iterations_max = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-min" == argument)
                {
                    mcmc_config.prerun_iterations_min = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-only" == argument)
                {
                    mcmc_config.need_prerun = true;
                    mcmc_config.store_prerun = true;
                    mcmc_config.need_main_run = false;

                    continue;
                }

                if ("--prerun-update" == argument)
                {
                    mcmc_config.prerun_iterations_update = destringify<unsigned>(*(++a));

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

                if ("--proposal" == argument)
                {
                    mcmc_config.proposal = *(++a);

                    if (mcmc_config.proposal == "MultivariateStudentT")
                    {
                        double dof = destringify<double>(*(++a));
                        if (dof <= 0)
                        {
                            throw DoUsage("No (or non-positive) degree of freedom for MultivariateStudentT specified");
                        }
                        mcmc_config.student_t_degrees_of_freedom = dof;
                    }

                    continue;
                }

                if ("--seed" == argument)
                {
                    std::string value(*(++a));

                    if ("time" == value)
                    {
                        mcmc_config.seed = ::time(0);
                    }
                    else
                    {
                        mcmc_config.seed = destringify<unsigned long>(value);
                    }

                    continue;
                }

                if ("--scale-reduction" == argument)
                {
                    scale_reduction = destringify<double>(*(++a));

                    continue;
                }

                if ("--store-prerun" == argument)
                {
                    mcmc_config.store_prerun = true;

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

        std::cout << std::scientific;
        std::cout << "# Samples generated by eos-sample-events-mcmc" << std::endl;

        if (! inst->density.get())
        {
            throw DoUsage("Need to specify a signal PDF to sample from");
        }

        unsigned size = std::distance(inst->density->begin(), inst->density->end());
        inst->mcmc_config.proposal_initial_covariance = std::vector<double>(size * size, 0.0);
        for (unsigned i = 0 ; i < size ; ++i)
        {
            inst->mcmc_config.proposal_initial_covariance[i + size * i] = 0.1;
        }

        MarkovChainSampler sampler(inst->density, inst->mcmc_config);

        sampler.run();
    }
    catch (DoUsage & e)
    {
    	// todo update
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan-mc" << std::endl;
        std::cout << "  [ [--kinematics NAME VALUE]* --observable NAME LOWER CENTRAL UPPER]+" << std::endl;
        std::cout << "  [--constraint NAME]+" << std::endl;
        std::cout << "  [ [ [--scan PARAMETER MIN MAX] | [--nuisance PARAMETER MIN MAX] ] --prior [flat | [gaussian LOWER CENTRAL UPPER] ] ]+" << std::endl;
        std::cout << "  [--chains VALUE]" << std::endl;
        std::cout << "  [--chunks VALUE]" << std::endl;
        std::cout << "  [--chunksize VALUE]" << std::endl;
        std::cout << "  [--debug]" << std::endl;
        std::cout << "  [--fix PARAMETER VALUE]+" << std::endl;
        std::cout << "  [--goodness_of_fit [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--no-prerun]" << std::endl;
        std::cout << "  [--optimize [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--output FILENAME]" << std::endl;
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
