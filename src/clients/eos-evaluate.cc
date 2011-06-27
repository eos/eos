/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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
#include <src/utils/cartesian-product.hh>
#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/log.hh>
#include <src/utils/power_of.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace eos;

double sign(const double & x)
{
    if (x < 0.0)
        return -1.0;

    return +1.0;
}

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

struct EvaluationInput
{
    ObservablePtr observable;

    CartesianProduct<std::vector<double>> ranges;

    std::shared_ptr<Kinematics> kinematics;

    std::vector<std::string> kinematic_names;

    EvaluationInput() :
        kinematics(new Kinematics)
    {
    }
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        std::vector<std::shared_ptr<EvaluationInput>> evaluation_inputs;

        std::vector<std::tuple<std::string, std::vector<Parameter>>> budgets;

        bool use_budget;

        CommandLine() :
            parameters(Parameters::Defaults()),
            budgets{std::make_tuple(std::string("delta"), std::vector<Parameter>())},
            use_budget(false)

        {
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_program_name("eos-evaluate");

            std::shared_ptr<EvaluationInput> evaluation_input(new EvaluationInput);

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));
                    double value = destringify<double>(*(++a));
                    evaluation_input->kinematics->declare(name);
                    evaluation_input->kinematic_names.push_back(name);

                    // create a vector with a single entry
                    std::vector<double> range(1, value);
                    evaluation_input->ranges.over(range);

                    continue;
                }

                if ("--range" == argument)
                {
                    std::string name = std::string(*(++a));
                    double min_value = destringify<double>(*(++a));
                    double max_value = destringify<double>(*(++a));
                    unsigned points  = destringify<unsigned>(*(++a));
                    evaluation_input->kinematics->declare(name);
                    evaluation_input->kinematic_names.push_back(name);

                    std::vector<double> range;
                    double increment = (max_value - min_value) / points;
                    for (std::size_t i = 0; i <= points; ++i)
                    {
                        range.push_back(min_value + i * increment);
                    }
                    evaluation_input->ranges.over(range);

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string name(*(++a));
                    evaluation_input->observable = Observable::make(name, parameters, *evaluation_input->kinematics, Options());
                    if (! evaluation_input->observable )
                        throw DoUsage("Unknown observable '" + name + "'");

                    evaluation_inputs.push_back(evaluation_input);

                    evaluation_input.reset(new EvaluationInput());

                    continue;
                }

                if ("--budget" == argument)
                {
                    std::string name(*(++a));

                    if (! use_budget)
                    {
                        use_budget = true;
                        std::get<0>(budgets.back()) = name;
                    }
                    else
                    {
                        budgets.push_back(std::make_tuple(name, std::vector<Parameter>()));
                    }

                    continue;
                }

                if ("--vary" == argument)
                {
                    std::string variation_name(*(++a));

                    try
                    {
                        Parameter variation = parameters[variation_name];
                        std::get<1>(budgets.back()).push_back(variation);
                    }
                    catch (UnknownParameterError & e)
                    {
                        throw DoUsage("Unknown parameter '" + variation_name + "'");
                    }

                    continue;
                }

                if ("--parameter" == argument)
                {
                    std::string parameter_name(*(++a));
                    double parameter_value = destringify<double>(*(++a));

                    try
                    {
                        Parameter parameter = parameters[parameter_name];
                        parameter = parameter_value;
                    }
                    catch (UnknownParameterError & e)
                    {
                        throw DoUsage("Unknown parameter '" + parameter_name + "'");
                    }

                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

void evaluate_with_sum_of_squares(const std::shared_ptr<EvaluationInput> evaluation_input)
{
    // print headlines
    std::cout << "# " << evaluation_input->observable->name()
              << ": " << evaluation_input->observable->options().as_string() << std::endl;

    std::cout << "# ";
    for (std::size_t i = 0 ; i < evaluation_input->kinematic_names.size() ; ++i)
    {
        std::cout << evaluation_input->kinematic_names[i] << '\t';
    }
    std::cout << "central";
    for (auto b = CommandLine::instance()->budgets.begin() ; b != CommandLine::instance()->budgets.end(); ++b)
    {
        std::cout << '\t' << std::get<0>(*b) << "_min\t" << std::get<0>(*b) << "_max";
    }
    std::cout << "\tdelta_min\tdelta_max" << std::endl;

    // iterate over all kinematical ranges
    for (auto r = evaluation_input->ranges.begin() ; r != evaluation_input->ranges.end() ; ++r)
    {
        // set the kinematics
        // for every dimension
        for (std::size_t i = 0 ; i < (*r).size() ; ++i)
        {
            evaluation_input->kinematics->set(evaluation_input->kinematic_names[i], (*r)[i]);
            std::cout << (*r)[i] << '\t';
        }

        double central = evaluation_input->observable->evaluate();

        std::cout << central;

        // do the variations
        double delta_max = 0.0, delta_min = 0.0;
        for (auto b = CommandLine::instance()->budgets.begin() ; b != CommandLine::instance()->budgets.end() ; ++b)
        {
            double budget_min = 0.0;
            double budget_max = 0.0;

            std::vector<Parameter> variations = std::get<1>(*b);

            for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
            {
                double old_v = *v;

                // raise value
                *v = v->max();

                double value = evaluation_input->observable->evaluate();

                if (value > central)
                {
                    budget_max += power_of<2>(value - central);
                }
                else if (value < central)
                {
                    budget_min += power_of<2>(value - central);
                }

                // lower value
                *v = v->min();

                value = evaluation_input->observable->evaluate();

                if (value > central)
                {
                    budget_max += power_of<2>(value - central);
                }
                else if (value < central)
                {
                    budget_min += power_of<2>(value - central);
                }

                *v = old_v;
            }

            delta_min += budget_min;
            delta_max += budget_max;

            std::cout << '\t' << std::sqrt(budget_min) << '\t' << std::sqrt(budget_max);
        }

        std::cout
            << '\t' << std::sqrt(delta_min) << '\t' << std::sqrt(delta_max)
            << "   (-" << std::abs(std::sqrt(delta_min) / central) * 100 << "% / +" << std::abs(std::sqrt(delta_max) / central) * 100 << "%)"
            << std::endl;
    }
}


int
main(int argc, char * argv[])
{
    try
    {
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->evaluation_inputs.empty())
            throw DoUsage("No input specified");

        for (auto i = CommandLine::instance()->evaluation_inputs.cbegin(), i_end = CommandLine::instance()->evaluation_inputs.cend() ; i != i_end ; ++i)
        {
            evaluate_with_sum_of_squares(*i);
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-evaluate" << std::endl;
        std::cout << "  --samples NUMBER" << std::endl;
        std::cout << "  [--vary PARAMETER]*" << std::endl;
        std::cout << "  [[--kinematics NAME VALUE]* --observable NAME MIN CENTRAL MAX]+" << std::endl;
    }
    catch(Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
