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

#include <algorithm>
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

        const std::string &
        what() const
        {
            return _what;
        }
};

class CommandLine : public InstantiationPolicy<CommandLine, Singleton>
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

        void
        parse(int argc, char ** argv)
        {
            Log::instance()->set_program_name("eos-list-parameters");

            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            for (char **a(argv + 1), **a_end(argv + argc); a != a_end; ++a)
            {
                std::string argument(*a);

                if ("--kinematics" == argument)
                {
                    std::string name  = std::string(*(++a));
                    double      value = destringify<double>(*(++a));
                    kinematics->declare(name);
                    kinematics->set(name, value);

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string   name(*(++a));
                    ObservablePtr observable = Observable::make(name, parameters, *kinematics, Options());
                    if (! observable)
                    {
                        throw DoUsage("Unknown observable '" + name + "'");
                    }

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

        std::string::size_type  max_name_length = 20;
        std::set<Parameter::Id> ids;

        if (CommandLine::instance()->observables.empty())
        {
            for (const auto & p : CommandLine::instance()->parameters)
            {
                ids.insert(p.id());
            }
        }
        else
        {
            for (auto & observable : CommandLine::instance()->observables)
            {
                ids.insert(observable->begin(), observable->end());
            }
        }

        for (const auto & p : CommandLine::instance()->parameters)
        {
            if (ids.end() == ids.find(p.id()))
            {
                continue;
            }

            max_name_length = std::max(max_name_length, p.name().length());
        }

        if (CommandLine::instance()->scan_format)
        {
            // how big to choose the range around mode
            double                 number_of_sigmas = 2.0;
            std::string::size_type max_prior_length = std::string("log-gamma").length();

            for (const auto & p : CommandLine::instance()->parameters)
            {
                std::string prior;
                std::string min, max;

                // upper and lower ranges
                double delta_down = p() - p.min();
                double delta_up   = p.max() - p();

                if ((delta_down == 0) && (delta_up == 0))
                {
                    prior = "flat";
                    min   = "MIN\t";
                    max   = "MAX\t";
                }
                else
                {
                    // for asymmetries larger than 5%, use log-gamma
                    prior = std::fabs(delta_up / delta_down - 1.0) < 0.05 ? "gaussian" : "log-gamma";
                    std::stringstream ss;
                    ss << std::setprecision(4) << double(p) - number_of_sigmas * delta_down;
                    min = ss.str();
                    ss.str("");
                    ss << std::setprecision(4) << double(p) + number_of_sigmas * delta_up;
                    max = ss.str();
                }

                std::cout << "    --scan" << '\t' << std::setw(max_name_length) << std::setiosflags(std::ios::left) << "\"" + p.name() + "\"" << '\t' << min << '\t' << max << '\t'
                          << "--prior" << '\t' << std::setw(max_prior_length) << std::setiosflags(std::ios::left) << prior;

                if (prior != "flat")
                {
                    std::cout << std::setw(7) << std::setprecision(4) << std::setiosflags(std::ios::left | std::ios::showpos) << '\t' << p.min() << '\t' << p() << '\t' << p.max();
                }

                std::cout << " \\" << std::endl;
            }
        }
        else
        {
            for (auto s = CommandLine::instance()->parameters.begin_sections(), s_end = CommandLine::instance()->parameters.end_sections(); s != s_end; ++s)
            {
                std::string section_title(s->name());

                std::cout << std::string(section_title.length(), '=') << '\n' << section_title << '\n' << std::string(section_title.length(), '=') << '\n' << '\n';

                for (const auto & g : *s)
                {
                    std::string group_title(g.name());

                    std::cout << group_title << '\n' << std::string(group_title.length(), '-') << '\n' << '\n';


                    std::vector<Parameter> group_parameters(g.begin(), g.end());
                    // nasty hack to sort
                    // TODO: remove entirely once all parameter names are QualifiedName compatible
                    std::sort(group_parameters.begin(),
                              group_parameters.end(),
                              [](const Parameter & x, const Parameter & y) -> bool
                    {
                        bool x_is_qualified_name = false;
                        try
                        {
                            QualifiedName qnx(x.name());
                            x_is_qualified_name = true;
                            QualifiedName qny(y.name());

                            if (qnx.prefix_part() < qny.prefix_part())
                            {
                                return true;
                            }
                            else if (qny.prefix_part() < qnx.prefix_part())
                            {
                                return false;
                            }

                            if (qnx.suffix_part() < qny.suffix_part())
                            {
                                return true;
                            }
                            else if (qny.suffix_part() < qnx.suffix_part())
                            {
                                return false;
                            }

                            return qnx.name_part() < qny.name_part();
                        }
                        catch (QualifiedNameSyntaxError & e)
                        {
                            if (x_is_qualified_name)
                            {
                                return false;
                            }

                            return x.name() < y.name();
                        }
                    });

                    for (const auto & p : group_parameters)
                    {
                        if (ids.end() == ids.find(p.id()))
                        {
                            continue;
                        }

                        std::cout << std::setw(max_name_length) << std::setiosflags(std::ios::right) << p.name() << '\t' << std::setw(7) << std::scientific << std::setprecision(4)
                                  << std::setiosflags(std::ios::left | std::ios::showpos) << p.min() << '\t' << p() << '\t' << p.max() << std::endl;
                    }

                    std::cout << std::endl;
                }
            }
        }
    }
    catch (DoUsage & e)
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
