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
#include <src/utils/destringify.hh>
#include <src/utils/histogram.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/log.hh>
#include <src/utils/one-of.hh>
#include <src/utils/power_of.hh>
#include <src/utils/random_number_engine.hh>
#include <src/utils/stringify.hh>
#include <src/utils/wilson-polynomial.hh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <random>
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

void evaluate_with_sum_of_squares(const ObservablePtr & observable);

void evaluate_with_samples(const ObservablePtr & observable);

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        std::list<ObservablePtr> observables;

        unsigned samples;

        double resolution;

        unsigned bins;

        bool show_histogram;

        std::list<Parameter> variations;

        std::function<void (const ObservablePtr &)> evaluate;

        CommandLine() :
            parameters(Parameters::Defaults()),
            samples(0),
            resolution(0.2),
            bins(40),
            show_histogram(false),
            evaluate(&evaluate_with_sum_of_squares)
        {
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_program_name("eos-evaluate");

            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);

                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));
                    double value = destringify<double>(*(++a));
                    kinematics->declare(name);
                    kinematics->set(name, value);

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string name(*(++a));
                    ObservablePtr observable = Observable::make(name, parameters, *kinematics, Options());
                    if (! observable)
                        throw DoUsage("Unknown observable '" + name + "'");

                    observables.push_back(observable);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--samples" == argument)
                {
                    samples = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--resolution" == argument)
                {
                    resolution = destringify<double>(*(++a));

                    continue;
                }

                if ("--bins" == argument)
                {
                    bins = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--show-histogram" == argument)
                {
                    show_histogram = true;

                    continue;
                }

                if ("--vary" == argument)
                {
                    std::string variation_name(*(++a));

                    try
                    {
                        Parameter variation = parameters[variation_name];
                        variations.push_back(variation);
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

                if ("--mode" == argument)
                {
                    std::string mode(*(++a));

                    if ("samples" == mode)
                    {
                        evaluate = &evaluate_with_samples;
                        continue;
                    }
                    else if ("sum-of-squares" == mode)
                    {
                        evaluate = &evaluate_with_sum_of_squares;
                        continue;
                    }

                    throw DoUsage("Unknown evaluation mode '" + mode + "'");
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

void evaluate_with_sum_of_squares(const ObservablePtr & observable)
{
    double central = observable->evaluate();

    std::cout << "# " << observable->name() << '[' << observable->kinematics().as_string() << "]: " << observable->options().as_string() << std::endl;

    double delta_max = 0.0, delta_min = 0.0;
    for (auto v = CommandLine::instance()->variations.begin(), v_end = CommandLine::instance()->variations.end() ; v != v_end ; ++v)
    {
        double old_v = *v;

        // raise value
        *v = v->max();

        double value = observable->evaluate();

        if (value > central)
        {
            delta_max += power_of<2>(value - central);
        }
        else if (value < central)
        {
            delta_min += power_of<2>(value - central);
        }

        // lower value
        *v = v->min();

        value = observable->evaluate();

        if (value > central)
        {
            delta_max += power_of<2>(value - central);
        }
        else if (value < central)
        {
            delta_min += power_of<2>(value - central);
        }

        *v = old_v;
    }

    delta_min = std::sqrt(delta_min);
    delta_max = std::sqrt(delta_max);

    std::cout
        << "#   " << central
        << " -" << delta_min
        << " +" << delta_max
        << "    (-" << std::abs(delta_min / central) * 100 << "% / +" << std::abs(delta_max / central) * 100 << "%)"
        << std::endl;
}

void evaluate_with_samples(const ObservablePtr & observable)
{
    double central = observable->evaluate();

    std::cout << "# " << observable->name() << '[' << observable->kinematics().as_string() << "]: " << observable->options().as_string() << std::endl;
    std::cout << "#   central = " << central << std::endl;

    double delta = CommandLine::instance()->resolution * std::abs(central);
    Histogram<1> histogram = Histogram<1>::WithEqualBinning(
            central - delta,
            central + delta,
            CommandLine::instance()->bins);
    histogram.insert(Histogram<1>::Bin(-std::numeric_limits<double>::max(), central - delta));
    histogram.insert(Histogram<1>::Bin(central + delta, std::numeric_limits<double>::max()));

    RandomNumberEngine engine;
    for (unsigned s = 0 ; s != CommandLine::instance()->samples ; ++s)
    {
        for (auto v = CommandLine::instance()->variations.begin(), v_end = CommandLine::instance()->variations.end() ; v != v_end ; ++v)
        {
            *v = v->sample(engine);
        }

        double value = observable->evaluate();
        if (isnan(value))
            continue;

        histogram.insert(value);
    }

    if (CommandLine::instance()->show_histogram)
    {
        std::cout << "lower\tupper\tentries" << std::endl;
        for (auto b = histogram.begin(), b_end = histogram.end() ; b != b_end ; ++b)
        {
            std::cout << b->lower << '\t' << b->upper << '\t' << b->value << std::endl;
        }
    }

    Histogram<1> ecdf = estimate_cumulative_distribution(histogram);

    if (CommandLine::instance()->show_histogram)
    {
        std::cout << std::endl << "ECDF" << std::endl;
        std::cout << "lower\tupper\tentries" << std::endl;
        for (auto b = ecdf.begin(), b_end = ecdf.end() ; b != b_end ; ++b)
        {
            std::cout << b->lower << '\t' << b->upper << '\t' << b->value << std::endl;
        }
    }

    std::vector<double> confidence_levels{ 0.683, 0.954, 0.997 };
    for (auto c = confidence_levels.cbegin(), c_end = confidence_levels.cend() ; c != c_end ; ++c)
    {
        double lower = -std::numeric_limits<double>::max(), upper = -std::numeric_limits<double>::max();
        for (auto b = ecdf.begin(), b_end = ecdf.end() ; b != b_end ; ++b)
        {
            if (b->value <= (1.0 - *c) / 2.0)
                lower = b->upper;

            if (b->value <= 1.0 - (1.0 - *c) / 2.0)
                upper = b->upper;
        }

        std::cout
            << "#   " << *c << "% CL: " << lower << " .. " << upper
            << "    (-" << std::abs((central - lower) / central) * 100 << "% / +" << std::abs((upper - central) / central) * 100 << "%)"
            << std::endl;
    }
}

int
main(int argc, char * argv[])
{
    try
    {
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->observables.empty())
            throw DoUsage("No input specified");

        std::function<void (const ObservablePtr &)> evaluate = CommandLine::instance()->evaluate;

        for (auto o = CommandLine::instance()->observables.cbegin(), o_end = CommandLine::instance()->observables.cend() ; o != o_end ; ++o)
        {
            evaluate(*o);
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
