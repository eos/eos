/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015, 2016 Danny van Dyk
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

#include <config.h>

#include <eos/constraint.hh>
#include <eos/observable.hh>
#include <eos/optimize/optimizer-gsl.hh>
#include <eos/statistics/log-posterior.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#include <yaml-cpp/yaml.h>

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

        LogPosterior log_posterior;

        std::vector<std::shared_ptr<hdf5::File>> prerun_inputs;

        std::vector<ParameterData> scan_parameters;

        std::vector<ParameterData> nuisance_parameters;

        std::vector<ParameterData> fix_parameters;

        std::vector<ObservableInput> inputs;

        std::vector<Constraint> constraints;

        std::string creator;

        std::vector<std::vector<double>> starting_points;

        std::shared_ptr<std::fstream> output;

        unsigned max_iterations;

        double target_precision;

        CommandLine() :
            parameters(Parameters::Defaults()),
            likelihood(parameters),
            log_posterior(likelihood),
            max_iterations(500),
            target_precision(1e-8)
        {
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
                    if (! log_posterior.add(prior, nuisance))
                        throw DoUsage("Error in assigning " + prior_type + " prior distribution to '" + name +
                                      "'. Perhaps '" + name + "' appears twice in the list of parameters?");

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

                if ("--fix" == argument)
                {
                    auto parameters = log_posterior.parameters();

                    std::string par_name = std::string(*(++a));
                    double value = destringify<double> (*(++a));

                    Parameter p = parameters[par_name];
                    p = value;
                    fix_parameters.push_back(ParameterData{ p, value, value, "fixed" });

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

                if ("--starting-point" == argument)
                {
                    // parse starting point
                    std::string lbrace(*(++a));
                    std::vector<double> point;
                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        double value = destringify<double>(word);
                        point.push_back(value);
                    }
                    while (true);
                    starting_points.push_back(point);

                    continue;
                }

                if ("--max-iterations" == argument)
                {
                    max_iterations = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--target-precision" == argument)
                {
                    target_precision = destringify<double>(*(++a));

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

                if ("--output" == argument)
                {
                    std::string filename(*(++a));
                    output = std::shared_ptr<std::fstream>(new std::fstream(filename, std::fstream::out | std::fstream::trunc));

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

        if (inst->inputs.empty() && inst->constraints.empty())
            throw DoUsage("Neither inputs nor constraints specified");

        if (inst->nuisance_parameters.empty() &&
            inst->scan_parameters.empty())
           throw  DoUsage("Neither scan nor nuisance parameters defined");

        std::cout << std::scientific;
        std::cout << "# Starting mode finding through eos-find-mode" << std::endl;
        if ( ! inst->scan_parameters.empty())
        {
            std::cout << "# Scan parameters (" << inst->scan_parameters.size() << "):" << std::endl;
            for (auto d = inst->log_posterior.parameter_descriptions().cbegin(), d_end = inst->log_posterior.parameter_descriptions().cend() ;
                 d != d_end ; ++d)
            {
                if (d->nuisance)
                    continue;
                std::cout << "#   " << inst->log_posterior.log_prior(d->parameter->name())->as_string() << std::endl;
            }
        }

        if ( ! inst->nuisance_parameters.empty())
        {
            std::cout << "# Nuisance parameters (" << inst->nuisance_parameters.size() << "):" << std::endl;
            for (auto d = inst->log_posterior.parameter_descriptions().cbegin(), d_end = inst->log_posterior.parameter_descriptions().cend() ;
                 d != d_end ; ++d)
            {
                if ( ! d->nuisance)
                    continue;
                std::cout << "#   " << inst->log_posterior.log_prior(d->parameter->name())->as_string() << std::endl;
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
        if (inst->starting_points.empty())
        {
            gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(rng, ::time(0));

            for (auto n = 0 ; n < 1000000 ; ++n)
            {
                std::vector<double> point;

                for (auto i = inst->log_posterior.parameter_descriptions().begin(), i_end = inst->log_posterior.parameter_descriptions().end() ; i != i_end ; ++i)
                {
                    LogPriorPtr prior = inst->log_posterior.log_prior(i->parameter->name());
                    point.push_back(prior->sample(rng));
                }

                inst->starting_points.push_back(point);
            }

            gsl_rng_free(rng);
        }

        for (const auto & point: inst->starting_points)
        {
            if (point.size() != inst->log_posterior.parameter_descriptions().size())
            {
                throw DoUsage("Starting point size of" + stringify(point.size())
                              + " doesn't match with analysis size of " + stringify(inst->log_posterior.parameter_descriptions().size()));
            }

            std::cout << std::endl;
            std::cout << "# Starting optimization at " << stringify_container(point, 4) << std::endl;
            std::cout << std::endl;

            DensityPtr density(new LogPosterior(inst->log_posterior));
            OptimizerPtr optimizer(new OptimizerGSL(density, inst->max_iterations, inst->target_precision));
            try
            {
                double maximum = optimizer->maximize();

                std::cout << "# Found maximum at: " << std::endl;
                std::cout << "#   (";
                for (auto && p : inst->log_posterior)
                {
                    std::cout << ' ' << p.parameter->evaluate();
                }
                std::cout << " )" << std::endl;
                std::cout << "#   value = " << maximum << std::endl;


                std::cout << "# Primary test statistics: " << std::endl;
                auto log_likelihood = inst->log_posterior.log_likelihood();
                for (auto c = log_likelihood.begin(), c_end = log_likelihood.end() ; c != c_end ; ++c)
                {
                    for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
                    {
                        std::cout << *(*b)->primary_test_statistic() << std::endl;
                    }
                }
            }
            catch (OptimizerError & e)
            {
                std::cout << "# Optimization failed:" << std::endl;
                std::cout << "\t" << e.what() << std::endl;
            }

            if (inst->output.get())
            {
                std::cout << "# writing mode into yaml output file" << std::endl;
                YAML::Emitter out;
                out.SetIndent(4);
                out << YAML::Comment("file generated by eos-find-mode");

                out << YAML::BeginMap;
                for (auto && d : inst->scan_parameters)
                {
                    out << YAML::Key << d.parameter.name();
                    out << YAML::Value;

                    out << YAML::BeginMap;

                    out << YAML::Key << "central";
                    out << YAML::Value << d.parameter.evaluate();

                    out << YAML::Key << "min";
                    out << YAML::Value << d.min;

                    out << YAML::Key << "max";
                    out << YAML::Value << d.max;

                    out << YAML::Key << "prior";
                    out << YAML::Value << d.prior;

                    out << YAML::EndMap;
                }
                for (auto && d : inst->nuisance_parameters)
                {
                    out << YAML::Key << d.parameter.name();
                    out << YAML::Value;

                    out << YAML::BeginMap;

                    out << YAML::Key << "central";
                    out << YAML::Value << d.parameter.evaluate();

                    out << YAML::Key << "min";
                    out << YAML::Value << d.min;

                    out << YAML::Key << "max";
                    out << YAML::Value << d.max;

                    out << YAML::Key << "prior";
                    out << YAML::Value << d.prior;

                    out << YAML::Key << "nuisance";
                    out << YAML::Value << true;

                    out << YAML::EndMap;
                }
                for (auto && d : inst->fix_parameters)
                {
                    out << YAML::Key << d.parameter.name();
                    out << YAML::Value;

                    out << YAML::BeginMap;

                    out << YAML::Key << "central";
                    out << YAML::Value << d.parameter.evaluate();

                    out << YAML::Key << "min";
                    out << YAML::Value << d.min;

                    out << YAML::Key << "max";
                    out << YAML::Value << d.max;

                    out << YAML::Key << "prior";
                    out << YAML::Value << d.prior;

                    out << YAML::EndMap;
                }
                *inst->output << out.c_str() << std::flush;
            }
        }
        inst->output->close();

        return EXIT_SUCCESS;
    }
    catch (DoUsage & e)
    {
        // todo update
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-find-mode" << std::endl;
        std::cout << "  [ [--kinematics NAME VALUE]* --observable NAME LOWER CENTRAL UPPER]+" << std::endl;
        std::cout << "  [--constraint NAME]+" << std::endl;
        std::cout << "  [ [ [--scan PARAMETER MIN MAX] | [--nuisance PARAMETER MIN MAX] ] --prior [flat | [gaussian LOWER CENTRAL UPPER] ] ]+" << std::endl;
        std::cout << "  [--debug]" << std::endl;
        std::cout << "  [--fix PARAMETER VALUE]+" << std::endl;
        std::cout << "  [--starting-point [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--max-iterations VALUE]" << std::endl;
        std::cout << "  [--target-precision VALUE]" << std::endl;

        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-find-mode \\" << std::endl;
        std::cout << "      --kinematics s_min 14.18 --kinematics s_max 16.00 \\" << std::endl;
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
