/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <eos/observable.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <iomanip>
#include <iostream>
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

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        std::vector<ObservablePtr> observables;

        bool scan_format;

        CommandLine() :
            parameters(Parameters::Defaults()),
            scan_format(false)
        {
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_program_name("eos-list-parameters");

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

                if ("--scan-format" == argument)
                {
                    scan_format = true;

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

        std::string::size_type max_name_length = 20;
        std::vector<Parameter> parameters;

        if (CommandLine::instance()->observables.empty())
        {
            parameters.insert(parameters.begin(), CommandLine::instance()->parameters.begin(), CommandLine::instance()->parameters.end());
        }
        else
        {
            std::set<Parameter::Id> ids;

            for (auto o = CommandLine::instance()->observables.begin(), o_end = CommandLine::instance()->observables.end() ; o != o_end ; ++o)
            {
                ids.insert((*o)->begin(), (*o)->end());
            }

            for (auto i = ids.cbegin(), i_end = ids.cend() ; i != i_end ; ++i)
            {
                parameters.push_back(CommandLine::instance()->parameters[*i]);
            }
        }

        for (auto p = parameters.begin(), p_end = parameters.end() ; p != p_end ; ++p)
        {
            max_name_length = std::max(max_name_length, p->name().length());
        }

        if (CommandLine::instance()->scan_format)
        {
            // how big to choose the range around mode
            double number_of_sigmas = 2.0;
            std::string::size_type max_prior_length = std::string("log-gamma").length();

            for (auto p = parameters.begin(), p_end = parameters.end() ; p != p_end ; ++p)
            {
                std::string prior;
                std::string min, max;

                // upper and lower ranges
                double delta_down = (*p)() - p->min();
                double delta_up = p->max() - (*p)();

                if ((delta_down == 0) && (delta_up == 0))
                {
                    prior = "flat";
                    min = "MIN\t";
                    max = "MAX\t";
                }
                else
                {
                    // for asymmetries larger than 5%, use log-gamma
                    prior =  std::fabs(delta_up / delta_down - 1.0) < 0.05 ? "gaussian" : "log-gamma";
                    std::stringstream ss;
                    ss
                    << std::setprecision(4)
                    << double(*p) - number_of_sigmas * delta_down;
                    min = ss.str();
                    ss.str("");
                    ss
                    << std::setprecision(4)
                    << double(*p) + number_of_sigmas * delta_up;
                    max = ss.str();
                }

                std::cout
                << "    --scan" << '\t'
                << std::setw(max_name_length) << std::setiosflags(std::ios::left)
                << "\"" + p->name() + "\"" << '\t'
                << min << '\t' << max << '\t'
                << "--prior" << '\t'
                << std::setw(max_prior_length) << std::setiosflags(std::ios::left)
                << prior;

                if (prior != "flat")
                {
                    std::cout
                    << std::setw(7) << std::setprecision(4) << std::setiosflags(std::ios::left | std::ios::showpos)
                    << '\t' << p->min() << '\t' << (*p)() << '\t' << p->max();
                }

                std::cout << " \\" << std::endl;
            }
        }
        else
        {
            for (auto p = parameters.begin(), p_end = parameters.end() ; p != p_end ; ++p)
            {
                std::cout
                << std::setw(max_name_length) << std::setiosflags(std::ios::right)
                << p->name() << '\t'
                << std::setw(7) << std::scientific << std::setprecision(4) << std::setiosflags(std::ios::left | std::ios::showpos)
                << p->min() << '\t' << (*p)() << '\t' << p->max()
                << std::endl;
            }
        }
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-list-parameters" << std::endl;
        std::cout << "  [[--kinematics NAME VALUE]* --observable NAME]*" << std::endl;
        std::cout << "  [--scan-format]" << std::endl;

        std::cout << "Print the parameter dependencies of a given observable (inclusive" << std::endl;
        std::cout << "its mandatory kinematics). If the scan-format option is given," << std::endl;
        std::cout << "the output is formatted in such a way that it can be used" << std::endl;
        std::cout << "as input to a call to eos-scan-mc." << std::endl;
        std::cout << "If no observable is specified, all parameters are listed." << std::endl;
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
