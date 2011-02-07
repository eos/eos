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

#include <src/rare-b-decays/factory.hh>
#include <src/utils/destringify.hh>
#include <src/utils/histogram.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/one-of.hh>
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
        Parameters parameters;

        std::list<ObservableInput> inputs;

        unsigned samples;

        double resolution;

        unsigned bins;

        bool show_histogram;

        std::list<Parameter> variations;

        CommandLine() :
            parameters(Parameters::Defaults()),
            samples(0),
            resolution(0.2),
            bins(40),
            show_histogram(false)
        {
        }

        void parse(int argc, char ** argv)
        {
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
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    kinematics.reset(new Kinematics);

                    ObservableOptions options;
                    std::string::size_type pos;
                    while (std::string::npos != (pos = observable_name.rfind(',')))
                    {
                        std::string::size_type sep(observable_name.find('=', pos + 1));
                        if (std::string::npos == sep)
                            throw DoUsage("Invalid observable option: '" + observable_name.substr(pos + 1) + "'");

                        std::string key(observable_name.substr(pos + 1, sep - pos - 1));
                        std::string value(observable_name.substr(sep + 1));

                        options.set(key, value);
                        observable_name.erase(pos);
                    }

                    input.observable = RareBFactory::make(observable_name, parameters, input.kinematics, options);
                    if (! input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    inputs.push_back(input);

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

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

int
main(int argc, char * argv[])
{
    try
    {
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->inputs.empty())
            throw DoUsage("No input specified");

        for (auto i = CommandLine::instance()->inputs.cbegin(), i_end = CommandLine::instance()->inputs.cend() ; i != i_end ; ++i)
        {
            double central = i->observable->evaluate();

            std::cout << "# " << i->observable->name() << '[' << i->kinematics.as_string() << "]: " << i->observable->options().as_string() << std::endl;
            std::cout << "#   central = " << central << std::endl;

            Histogram<1> histogram = Histogram<1>::WithEqualBinning(
                    (1.0 - CommandLine::instance()->resolution) * central,
                    (1.0 + CommandLine::instance()->resolution) * central,
                    CommandLine::instance()->bins);
            histogram.insert(Histogram<1>::Bin(0.0, (1.0 - CommandLine::instance()->resolution) * central));
            histogram.insert(Histogram<1>::Bin((1.0 + CommandLine::instance()->resolution) * central, std::numeric_limits<double>::max()));

            RandomNumberEngine engine;
            for (unsigned s = 0 ; s != CommandLine::instance()->samples ; ++s)
            {
                for (auto v = CommandLine::instance()->variations.begin(), v_end = CommandLine::instance()->variations.end() ; v != v_end ; ++v)
                {
                    std::normal_distribution<double> distribution(v->central(), v->max() - v->min());
                    *v = distribution(engine);
                }

                double value = i->observable->evaluate();
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
                double lower, upper;
                for (auto b = ecdf.begin(), b_end = ecdf.end() ; b != b_end ; ++b)
                {
                    if (b->value <= (1.0 - *c) / 2.0)
                        lower = b->upper;

                    if (b->value <= 1.0 - (1.0 - *c) / 2.0)
                        upper = b->upper;
                }

                std::cout << "#   " << *c << "% CL: " << lower << " .. " << upper << std::endl;
            }
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
